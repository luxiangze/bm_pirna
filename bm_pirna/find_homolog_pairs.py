#!/usr/bin/env python3
"""Find homologous protein pairs with DIAMOND.

The command accepts two protein FASTA files, runs reciprocal protein searches,
maps protein accessions back to genes, and writes a deterministic tabular
summary.  FASTA headers, the protein download manifest, or explicit mapping
files can provide the protein-to-gene mapping.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import UTC, datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Iterable, Literal

from Bio import SeqIO
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

DEFAULT_MIN_IDENTITY = 30.0
DEFAULT_QUERY_COVER = 50.0
DEFAULT_SUBJECT_COVER = 50.0
DEFAULT_EVALUE = 1e-5
DEFAULT_SENSITIVITY = "very-sensitive"
DEFAULT_THREADS = 4
DEFAULT_MANIFEST_NAME = "protein_download_manifest.tsv"

DIAMOND_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qcovhsp",
    "scovhsp",
    "evalue",
    "bitscore",
)
SENSITIVITY_LEVELS = {
    "default",
    "fast",
    "mid-sensitive",
    "sensitive",
    "more-sensitive",
    "very-sensitive",
    "ultra-sensitive",
}
OUTPUT_COLUMNS = (
    "gene1",
    "gene2",
    "homology",
    "gene1_label",
    "gene2_label",
    "protein1",
    "protein2",
    "alignment_length",
    "query_coverage_pct",
    "subject_coverage_pct",
    "evalue",
    "bitscore",
    "pair_method",
)

NCBI_GENE_RE = re.compile(r"\[GeneID=(\d+)\]", re.IGNORECASE)
FLYBASE_PARENT_RE = re.compile(r"\bparent=((?:FBgn)\d+)", re.IGNORECASE)
FLYBASE_GENE_RE = re.compile(r"^FBgn\d+$", re.IGNORECASE)
NUMERIC_GENE_RE = re.compile(r"^\d+$")


class HomologyError(RuntimeError):
    """Raised when input validation or a DIAMOND run fails."""


@dataclass(frozen=True)
class FastaRecord:
    """A validated FASTA record."""

    protein_id: str
    sequence_length: int
    description: str


@dataclass(frozen=True)
class ProteinMapping:
    """A protein-to-gene mapping and an optional human-readable label."""

    gene_id: str
    label: str


@dataclass(frozen=True)
class Hit:
    """One DIAMOND tabular alignment result."""

    qseqid: str
    sseqid: str
    pident: float
    length: int
    qlen: int
    slen: int
    qcov: float
    scov: float
    evalue: float
    bitscore: float


@dataclass(frozen=True)
class PairResult:
    """A selected gene pair represented by its best protein alignment."""

    gene1: str
    gene2: str
    hit: Hit
    pair_method: str


@dataclass(frozen=True)
class AnalysisSummary:
    """Small summary returned by :func:`run_analysis`."""

    output: Path
    metadata: Path
    pair_count: int
    forward_hits: int
    reverse_hits: int


def _normalize_field(value: str) -> str:
    """Normalize a mapping or manifest field name."""

    return re.sub(r"[^a-z0-9]+", "", value.removeprefix("\ufeff").strip().lower())


def _normalize_gene_id(value: str) -> str:
    """Canonicalize common NCBI and FlyBase gene identifiers."""

    value = value.strip()
    if NUMERIC_GENE_RE.fullmatch(value):
        return str(int(value))
    if FLYBASE_GENE_RE.fullmatch(value):
        return f"FBgn{value[4:]}"
    return value


def _clean_value(value: str) -> str:
    """Strip whitespace and optional FASTA-style leading punctuation."""

    value = value.removeprefix("\ufeff").strip()
    return value[1:].strip() if value.startswith(">") else value


def _read_fasta(path: Path) -> dict[str, FastaRecord]:
    """Read and validate a protein FASTA file."""

    if not path.exists():
        raise HomologyError(f"FASTA file does not exist: {path}")
    if not path.is_file():
        raise HomologyError(f"FASTA path is not a file: {path}")

    records: dict[str, FastaRecord] = {}
    try:
        with path.open("rt", encoding="utf-8") as handle:
            iterator = SeqIO.parse(handle, "fasta")
            for record in iterator:
                protein_id = record.id.strip()
                if not protein_id:
                    raise HomologyError(f"FASTA contains a record without an ID: {path}")
                if protein_id in records:
                    raise HomologyError(f"Duplicate protein accession {protein_id!r} in {path}")
                sequence = str(record.seq).replace(" ", "").replace("\r", "").replace("\n", "")
                if not sequence:
                    raise HomologyError(f"FASTA record {protein_id!r} has an empty sequence")
                records[protein_id] = FastaRecord(
                    protein_id=protein_id,
                    sequence_length=len(sequence),
                    description=record.description,
                )
    except (OSError, ValueError) as exc:
        raise HomologyError(f"Could not parse FASTA {path}: {exc}") from exc
    if not records:
        raise HomologyError(f"FASTA contains no records: {path}")
    return records


def _header_mapping(record: FastaRecord) -> ProteinMapping | None:
    """Extract a gene mapping from common NCBI and FlyBase FASTA headers."""

    description = record.description
    match = NCBI_GENE_RE.search(description)
    if match:
        gene_id = match.group(1)
        prefix = description.split("[", 1)[0].strip()
        label_tokens = prefix.split(None, 1)
        label = label_tokens[1].strip() if len(label_tokens) == 2 else gene_id
        return ProteinMapping(gene_id=gene_id, label=label or gene_id)

    match = FLYBASE_PARENT_RE.search(description)
    if match:
        gene_id = _normalize_gene_id(match.group(1))
        name_match = re.search(r"(?:^|;)\s*name=([^;]+)", description, re.IGNORECASE)
        label = name_match.group(1).strip() if name_match else gene_id
        return ProteinMapping(gene_id=gene_id, label=label or gene_id)

    if FLYBASE_GENE_RE.fullmatch(record.protein_id):
        gene_id = _normalize_gene_id(record.protein_id)
        return ProteinMapping(gene_id=gene_id, label=gene_id)
    if NUMERIC_GENE_RE.fullmatch(record.protein_id):
        gene_id = _normalize_gene_id(record.protein_id)
        return ProteinMapping(gene_id=gene_id, label=gene_id)
    return None


def _merge_mapping(
    mapping: dict[str, ProteinMapping],
    protein_id: str,
    gene_id: str,
    label: str,
) -> None:
    """Add one mapping and reject conflicting gene assignments."""

    gene_id = _normalize_gene_id(_clean_value(gene_id))
    label = _clean_value(label) or gene_id
    if not gene_id:
        return
    previous = mapping.get(protein_id)
    if previous and previous.gene_id != gene_id:
        raise HomologyError(
            f"Protein {protein_id!r} maps to multiple genes: {previous.gene_id!r} and {gene_id!r}"
        )
    if previous:
        labels = sorted({part for part in previous.label.split(";") if part} | {label})
        mapping[protein_id] = ProteinMapping(gene_id=gene_id, label=";".join(labels))
    else:
        mapping[protein_id] = ProteinMapping(gene_id=gene_id, label=label)


def _mapping_from_headers(records: dict[str, FastaRecord]) -> dict[str, ProteinMapping]:
    """Build mappings from FASTA descriptions."""

    return {
        protein_id: header_mapping
        for protein_id, record in records.items()
        if (header_mapping := _header_mapping(record)) is not None
    }


def _detect_delimiter(lines: list[str], path: Path) -> str:
    """Detect CSV or TSV delimiters for mapping files."""

    if path.suffix.lower() == ".csv":
        return ","
    has_comma = any("," in line for line in lines[:20])
    has_tab = any("\t" in line for line in lines[:20])
    return "," if has_comma and not has_tab else "\t"


def _read_mapping_file(path: Path) -> dict[str, ProteinMapping]:
    """Read a two-column protein-to-gene mapping file."""

    if not path.exists() or not path.is_file():
        raise HomologyError(f"Mapping file does not exist: {path}")
    text = path.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    if not lines:
        raise HomologyError(f"Mapping file is empty: {path}")
    delimiter = _detect_delimiter(lines, path)
    rows = list(csv.reader(lines, delimiter=delimiter))
    rows = [row for row in rows if any(cell.strip() for cell in row)]
    if not rows:
        raise HomologyError(f"Mapping file is empty: {path}")

    aliases = {
        "protein": {"protein", "proteinid", "proteinaccession", "accession", "seqid"},
        "gene": {"gene", "geneid", "resolvedgeneid", "geneidentifier"},
        "label": {"label", "genelabel", "genename", "genesymbol", "rawid"},
    }
    normalized_first = [_normalize_field(cell) for cell in rows[0]]
    has_header = any(
        normalized in aliases[category] for normalized in normalized_first for category in aliases
    )
    if has_header:
        columns = {normalized: index for index, normalized in enumerate(normalized_first)}
        protein_index = next(
            (columns[name] for name in aliases["protein"] if name in columns), None
        )
        gene_index = next((columns[name] for name in aliases["gene"] if name in columns), None)
        label_index = next((columns[name] for name in aliases["label"] if name in columns), None)
        if protein_index is None or gene_index is None:
            raise HomologyError(f"Mapping file {path} must contain protein_id and gene_id columns")
        data_rows = rows[1:]
    else:
        protein_index, gene_index, label_index = 0, 1, None
        data_rows = rows

    mapping: dict[str, ProteinMapping] = {}
    for row in data_rows:
        if max(protein_index, gene_index) >= len(row):
            continue
        protein_id = _clean_value(row[protein_index])
        gene_id = _clean_value(row[gene_index])
        label = (
            _clean_value(row[label_index])
            if label_index is not None and label_index < len(row)
            else gene_id
        )
        if protein_id and gene_id:
            _merge_mapping(mapping, protein_id, gene_id, label)
    if not mapping:
        raise HomologyError(f"Mapping file contains no usable mappings: {path}")
    return mapping


def _source_stem(path_value: str | Path) -> str:
    """Return the input stem shared by downloader outputs and source tables."""

    name = Path(path_value).name
    for suffix in (".proteins.fasta", ".proteins.faa", ".proteins.fa", ".fasta", ".faa", ".fa"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _manifest_mapping(path: Path, fasta_path: Path) -> dict[str, ProteinMapping]:
    """Read protein mappings from a downloader manifest."""

    if not path.exists() or not path.is_file():
        raise HomologyError(f"Manifest does not exist: {path}")
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = reader.fieldnames or []
            fields = {_normalize_field(field): field for field in fieldnames}
            accession_field = fields.get("proteinaccession")
            gene_field = fields.get("resolvedgeneid")
            input_field = fields.get("inputfile")
            raw_field = fields.get("rawid")
            status_field = fields.get("status")
            if not accession_field or not gene_field:
                raise HomologyError(
                    f"Manifest {path} must contain protein_accession and resolved_gene_id"
                )
            rows = list(reader)
    except OSError as exc:
        raise HomologyError(f"Could not read manifest {path}: {exc}") from exc

    expected_stem = _source_stem(fasta_path)
    matching_rows = [
        row
        for row in rows
        if input_field and row.get(input_field) and _source_stem(row[input_field]) == expected_stem
    ]
    if matching_rows:
        rows = matching_rows

    mapping: dict[str, ProteinMapping] = {}
    for row in rows:
        if status_field and row.get(status_field, "").strip().lower() not in {"", "downloaded"}:
            continue
        protein_id = _clean_value(row.get(accession_field, ""))
        gene_id = _clean_value(row.get(gene_field, ""))
        label = _clean_value(row.get(raw_field, "")) if raw_field else gene_id
        if protein_id and gene_id:
            _merge_mapping(mapping, protein_id, gene_id, label)
    if not mapping:
        raise HomologyError(f"Manifest has no mappings for FASTA {fasta_path}")
    return mapping


def _default_manifest(fasta_path: Path) -> Path | None:
    """Find the downloader manifest beside a FASTA file, if present."""

    candidates = (
        fasta_path.parent / DEFAULT_MANIFEST_NAME,
        fasta_path.parent.parent / DEFAULT_MANIFEST_NAME,
    )
    return next((candidate for candidate in candidates if candidate.exists()), None)


def _resolve_mapping(
    fasta_path: Path,
    records: dict[str, FastaRecord],
    *,
    manifest: Path | None,
    explicit_map: Path | None,
    allow_unmapped: bool,
) -> tuple[dict[str, ProteinMapping], str, list[str]]:
    """Resolve all FASTA proteins to genes and report the mapping source."""

    header_mapping = _mapping_from_headers(records)
    if explicit_map:
        mapping = _read_mapping_file(explicit_map)
        source = f"mapping:{explicit_map}"
    else:
        manifest_path = manifest or _default_manifest(fasta_path)
        if manifest_path:
            mapping = _manifest_mapping(manifest_path, fasta_path)
            source = f"manifest:{manifest_path}"
            for protein_id, value in header_mapping.items():
                mapping.setdefault(protein_id, value)
        else:
            mapping = header_mapping
            source = "fasta_headers"

    unmapped = sorted(set(records) - set(mapping))
    if unmapped and allow_unmapped:
        for protein_id in unmapped:
            mapping[protein_id] = ProteinMapping(gene_id=protein_id, label=protein_id)
    if unmapped and not allow_unmapped:
        preview = ", ".join(unmapped[:5])
        suffix = " ..." if len(unmapped) > 5 else ""
        raise HomologyError(
            f"Could not map {len(unmapped)} protein(s) in {fasta_path} to genes: "
            f"{preview}{suffix}. Provide --manifest or an explicit --map file."
        )
    return mapping, source, unmapped


def _run_command(command: list[str]) -> subprocess.CompletedProcess[str]:
    """Run a subprocess without invoking a shell."""

    try:
        return subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            env=os.environ.copy(),
        )
    except FileNotFoundError as exc:
        raise HomologyError(
            "The 'diamond' executable was not found; install diamond with pixi"
        ) from exc
    except OSError as exc:
        raise HomologyError(f"Could not start DIAMOND: {exc}") from exc


def _check_command(result: subprocess.CompletedProcess[str], command: list[str]) -> None:
    """Raise a concise error for a failed external command."""

    if result.returncode:
        detail = (result.stderr or result.stdout or "DIAMOND returned a non-zero status").strip()
        raise HomologyError(f"DIAMOND command failed ({result.returncode}): {detail}")


def _diamond_version(executable: str) -> str:
    """Return the installed DIAMOND version string."""

    result = _run_command([executable, "version"])
    _check_command(result, [executable, "version"])
    text = (result.stdout or result.stderr).strip()
    return text.splitlines()[0] if text else "unknown"


def _run_diamond_search(
    executable: str,
    query: Path,
    target: Path,
    database_prefix: Path,
    output: Path,
    *,
    evalue: float,
    threads: int,
    sensitivity: str,
) -> None:
    """Build a DIAMOND database and run one protein search."""

    makedb_command = [
        executable,
        "makedb",
        "--in",
        str(target),
        "--db",
        str(database_prefix),
    ]
    result = _run_command(makedb_command)
    _check_command(result, makedb_command)

    search_command = [
        executable,
        "blastp",
        "--query",
        str(query),
        "--db",
        str(database_prefix),
        "--out",
        str(output),
        "--outfmt",
        "6",
        *DIAMOND_FIELDS,
        "--evalue",
        f"{evalue:.12g}",
        "--max-target-seqs",
        "0",
        "--max-hsps",
        "1",
        "--threads",
        str(threads),
    ]
    if sensitivity != "default":
        search_command.append(f"--{sensitivity}")
    result = _run_command(search_command)
    _check_command(result, search_command)
    if not output.exists():
        raise HomologyError(f"DIAMOND completed without creating output: {output}")


def _float(value: str, field: str, line_number: int) -> float:
    """Parse a finite DIAMOND numeric field."""

    try:
        result = float(value)
    except ValueError as exc:
        raise HomologyError(f"Invalid {field} on DIAMOND line {line_number}: {value!r}") from exc
    if not math.isfinite(result):
        raise HomologyError(f"Non-finite {field} on DIAMOND line {line_number}: {value!r}")
    return result


def _parse_hits(path: Path) -> list[Hit]:
    """Parse the fixed DIAMOND tabular output fields."""

    hits: list[Hit] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise HomologyError(f"Could not read DIAMOND output {path}: {exc}") from exc
    for line_number, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < len(DIAMOND_FIELDS):
            raise HomologyError(
                f"Expected {len(DIAMOND_FIELDS)} DIAMOND fields on line {line_number}, "
                f"got {len(fields)}"
            )
        qseqid, sseqid = fields[0], fields[1]
        try:
            length = int(fields[3])
            qlen = int(fields[4])
            slen = int(fields[5])
        except ValueError as exc:
            raise HomologyError(f"Invalid alignment length on DIAMOND line {line_number}") from exc
        if min(length, qlen, slen) <= 0:
            raise HomologyError(f"Non-positive sequence length on DIAMOND line {line_number}")
        qcov = _float(fields[6], "qcovhsp", line_number)
        scov = _float(fields[7], "scovhsp", line_number)
        if qcov < 0:
            qcov = 100.0 * length / qlen
        if scov < 0:
            scov = 100.0 * length / slen
        hits.append(
            Hit(
                qseqid=qseqid,
                sseqid=sseqid,
                pident=_float(fields[2], "pident", line_number),
                length=length,
                qlen=qlen,
                slen=slen,
                qcov=qcov,
                scov=scov,
                evalue=_float(fields[8], "evalue", line_number),
                bitscore=_float(fields[9], "bitscore", line_number),
            )
        )
    return hits


def _hit_key(hit: Hit) -> tuple[Any, ...]:
    """Return the deterministic best-hit sort key."""

    return (
        -hit.bitscore,
        hit.evalue,
        -hit.pident,
        -hit.qcov,
        -hit.scov,
        -hit.length,
        hit.qseqid,
        hit.sseqid,
    )


def _filter_hits(
    hits: Iterable[Hit],
    *,
    min_identity: float,
    min_query_cover: float,
    min_subject_cover: float,
    evalue: float,
) -> list[Hit]:
    """Apply identity, coverage, and E-value thresholds."""

    return [
        hit
        for hit in hits
        if hit.pident >= min_identity
        and hit.qcov >= min_query_cover
        and hit.scov >= min_subject_cover
        and hit.evalue <= evalue
    ]


def _gene_pair_hits(
    hits: Iterable[Hit],
    query_mapping: dict[str, ProteinMapping],
    subject_mapping: dict[str, ProteinMapping],
) -> dict[tuple[str, str], Hit]:
    """Collapse protein alignments to the best alignment for each gene pair."""

    grouped: dict[tuple[str, str], Hit] = {}
    for hit in hits:
        try:
            gene1 = query_mapping[hit.qseqid].gene_id
            gene2 = subject_mapping[hit.sseqid].gene_id
        except KeyError as exc:
            raise HomologyError(
                f"DIAMOND output contains an unknown protein ID: {exc.args[0]}"
            ) from exc
        key = (gene1, gene2)
        if key not in grouped or _hit_key(hit) < _hit_key(grouped[key]):
            grouped[key] = hit
    return grouped


def _best_by_gene(pair_hits: dict[tuple[str, str], Hit]) -> dict[str, tuple[str, Hit]]:
    """Select one best target gene for each query gene."""

    best: dict[str, tuple[str, Hit]] = {}
    for (query_gene, subject_gene), hit in pair_hits.items():
        previous = best.get(query_gene)
        candidate = (subject_gene, hit)
        if previous is None or _hit_key(hit) < _hit_key(previous[1]):
            best[query_gene] = candidate
    return best


def select_pairs(
    forward_hits: Iterable[Hit],
    reverse_hits: Iterable[Hit],
    mapping1: dict[str, ProteinMapping],
    mapping2: dict[str, ProteinMapping],
    *,
    mode: Literal["rbh", "best", "all"],
    min_identity: float,
    min_query_cover: float,
    min_subject_cover: float,
    evalue: float,
) -> tuple[list[PairResult], int, int]:
    """Filter alignments and select all, one-way-best, or reciprocal-best pairs."""

    forward_filtered = _filter_hits(
        forward_hits,
        min_identity=min_identity,
        min_query_cover=min_query_cover,
        min_subject_cover=min_subject_cover,
        evalue=evalue,
    )
    reverse_filtered = _filter_hits(
        reverse_hits,
        min_identity=min_identity,
        min_query_cover=min_query_cover,
        min_subject_cover=min_subject_cover,
        evalue=evalue,
    )
    forward_pairs = _gene_pair_hits(forward_filtered, mapping1, mapping2)
    reverse_pairs = _gene_pair_hits(reverse_filtered, mapping2, mapping1)

    if mode == "all":
        selected = [
            PairResult(gene1, gene2, hit, "all") for (gene1, gene2), hit in forward_pairs.items()
        ]
    else:
        forward_best = _best_by_gene(forward_pairs)
        if mode == "best":
            selected = [
                PairResult(gene1, gene2, hit, "best")
                for gene1, (gene2, hit) in forward_best.items()
            ]
        else:
            reverse_best = _best_by_gene(reverse_pairs)
            selected = []
            for gene1, (gene2, hit) in forward_best.items():
                reciprocal = reverse_best.get(gene2)
                if reciprocal and reciprocal[0] == gene1:
                    selected.append(PairResult(gene1, gene2, hit, "rbh"))

    selected.sort(key=lambda pair: (pair.gene1, pair.gene2, pair.hit.qseqid, pair.hit.sseqid))
    return selected, len(forward_filtered), len(reverse_filtered)


def _default_output(fasta1: Path, fasta2: Path) -> Path:
    """Construct the default output path."""

    directory = fasta1.parent / "homolog_pairs"
    name = f"{_source_stem(fasta1)}__vs__{_source_stem(fasta2)}.homolog_pairs.tsv"
    return directory / name


def _sha256(path: Path) -> str:
    """Calculate a file checksum for metadata."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _format_float(value: float) -> str:
    """Format a numeric output field without unnecessary trailing zeroes."""

    return f"{value:.8g}"


def _write_outputs(
    output: Path,
    metadata_path: Path,
    pairs: list[PairResult],
    mapping1: dict[str, ProteinMapping],
    mapping2: dict[str, ProteinMapping],
    metadata: dict[str, Any],
    *,
    force: bool,
) -> None:
    """Write the pair table and JSON metadata."""

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for pair in pairs:
            hit = pair.hit
            writer.writerow(
                {
                    "gene1": pair.gene1,
                    "gene2": pair.gene2,
                    "homology": _format_float(hit.pident),
                    "gene1_label": mapping1[hit.qseqid].label,
                    "gene2_label": mapping2[hit.sseqid].label,
                    "protein1": hit.qseqid,
                    "protein2": hit.sseqid,
                    "alignment_length": hit.length,
                    "query_coverage_pct": _format_float(hit.qcov),
                    "subject_coverage_pct": _format_float(hit.scov),
                    "evalue": _format_float(hit.evalue),
                    "bitscore": _format_float(hit.bitscore),
                    "pair_method": pair.pair_method,
                }
            )
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def run_analysis(
    fasta1: Path,
    fasta2: Path,
    *,
    output: Path | None = None,
    manifest: Path | None = None,
    map1: Path | None = None,
    map2: Path | None = None,
    mode: Literal["rbh", "best", "all"] = "rbh",
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_query_cover: float = DEFAULT_QUERY_COVER,
    min_subject_cover: float = DEFAULT_SUBJECT_COVER,
    evalue: float = DEFAULT_EVALUE,
    threads: int = DEFAULT_THREADS,
    sensitivity: str = DEFAULT_SENSITIVITY,
    force: bool = False,
    keep_intermediate: bool = False,
    allow_unmapped: bool = False,
    diamond: str | None = None,
) -> AnalysisSummary:
    """Run the complete two-way DIAMOND analysis."""

    fasta1 = fasta1.resolve()
    fasta2 = fasta2.resolve()
    if fasta1 == fasta2:
        raise HomologyError("FASTA1 and FASTA2 must be different files")
    if mode not in {"rbh", "best", "all"}:
        raise HomologyError(f"Unsupported pair mode: {mode}")
    if not 0 <= min_identity <= 100:
        raise HomologyError("--min-identity must be between 0 and 100")
    if not 0 <= min_query_cover <= 100 or not 0 <= min_subject_cover <= 100:
        raise HomologyError("Coverage thresholds must be between 0 and 100")
    if evalue <= 0 or not math.isfinite(evalue):
        raise HomologyError("--evalue must be a positive finite number")
    if threads < 1:
        raise HomologyError("--threads must be at least 1")
    if sensitivity not in SENSITIVITY_LEVELS:
        raise HomologyError(
            f"Unsupported sensitivity {sensitivity!r}; choose from {sorted(SENSITIVITY_LEVELS)}"
        )

    output = (output or _default_output(fasta1, fasta2)).resolve()
    metadata_path = output.with_name(f"{output.stem}.metadata.json")
    if not force and (output.exists() or metadata_path.exists()):
        raise HomologyError(f"Output already exists: {output} (use --force to overwrite)")

    records1 = _read_fasta(fasta1)
    records2 = _read_fasta(fasta2)
    mapping1, mapping_source1, unmapped1 = _resolve_mapping(
        fasta1,
        records1,
        manifest=manifest,
        explicit_map=map1,
        allow_unmapped=allow_unmapped,
    )
    mapping2, mapping_source2, unmapped2 = _resolve_mapping(
        fasta2,
        records2,
        manifest=manifest,
        explicit_map=map2,
        allow_unmapped=allow_unmapped,
    )

    executable = diamond or shutil.which("diamond") or "diamond"
    version = _diamond_version(executable)
    temp_path: Path | None = None
    if keep_intermediate:
        output.parent.mkdir(parents=True, exist_ok=True)
        temp_path = Path(tempfile.mkdtemp(prefix=f"{output.stem}.diamond.", dir=output.parent))
        cleanup = False
    else:
        temp_path = Path(tempfile.mkdtemp(prefix="diamond_homology_"))
        cleanup = True

    try:
        forward_raw = temp_path / "forward.tsv"
        reverse_raw = temp_path / "reverse.tsv"
        _run_diamond_search(
            executable,
            fasta1,
            fasta2,
            temp_path / "database2",
            forward_raw,
            evalue=evalue,
            threads=threads,
            sensitivity=sensitivity,
        )
        _run_diamond_search(
            executable,
            fasta2,
            fasta1,
            temp_path / "database1",
            reverse_raw,
            evalue=evalue,
            threads=threads,
            sensitivity=sensitivity,
        )
        forward_raw_hits = _parse_hits(forward_raw)
        reverse_raw_hits = _parse_hits(reverse_raw)
        pairs, forward_count, reverse_count = select_pairs(
            forward_raw_hits,
            reverse_raw_hits,
            mapping1,
            mapping2,
            mode=mode,
            min_identity=min_identity,
            min_query_cover=min_query_cover,
            min_subject_cover=min_subject_cover,
            evalue=evalue,
        )
        metadata = {
            "created_at_utc": datetime.now(UTC).isoformat(),
            "command": sys.argv,
            "diamond": {"executable": executable, "version": version},
            "inputs": [
                {
                    "path": str(fasta1),
                    "sha256": _sha256(fasta1),
                    "protein_count": len(records1),
                    "mapping_source": mapping_source1,
                    "unmapped_proteins": unmapped1,
                },
                {
                    "path": str(fasta2),
                    "sha256": _sha256(fasta2),
                    "protein_count": len(records2),
                    "mapping_source": mapping_source2,
                    "unmapped_proteins": unmapped2,
                },
            ],
            "parameters": {
                "mode": mode,
                "min_identity": min_identity,
                "min_query_cover": min_query_cover,
                "min_subject_cover": min_subject_cover,
                "evalue": evalue,
                "threads": threads,
                "sensitivity": sensitivity,
                "allow_unmapped": allow_unmapped,
            },
            "counts": {
                "forward_raw_hits": len(forward_raw_hits),
                "reverse_raw_hits": len(reverse_raw_hits),
                "forward_filtered_hits": forward_count,
                "reverse_filtered_hits": reverse_count,
                "selected_pairs": len(pairs),
            },
            "intermediate_directory": str(temp_path) if keep_intermediate else None,
        }
        _write_outputs(
            output,
            metadata_path,
            pairs,
            mapping1,
            mapping2,
            metadata,
            force=force,
        )
    finally:
        if cleanup and temp_path is not None:
            shutil.rmtree(temp_path, ignore_errors=True)

    return AnalysisSummary(
        output=output,
        metadata=metadata_path,
        pair_count=len(pairs),
        forward_hits=forward_count,
        reverse_hits=reverse_count,
    )


@app.command()
def main(
    fasta1: Path = typer.Argument(..., help="First protein FASTA (gene1 side)."),
    fasta2: Path = typer.Argument(..., help="Second protein FASTA (gene2 side)."),
    output: Path | None = typer.Option(None, "--output", "-o", help="Output TSV path."),
    manifest: Path | None = typer.Option(
        None,
        "--manifest",
        help="Shared protein_download_manifest.tsv; auto-discovered beside FASTA files.",
    ),
    map1: Path | None = typer.Option(None, "--map1", help="Protein-to-gene map for FASTA1."),
    map2: Path | None = typer.Option(None, "--map2", help="Protein-to-gene map for FASTA2."),
    mode: str = typer.Option("rbh", "--mode", help="Pair mode: rbh, best, or all."),
    min_identity: float = typer.Option(
        DEFAULT_MIN_IDENTITY, "--min-identity", help="Minimum amino-acid identity percentage."
    ),
    min_query_cover: float = typer.Option(
        DEFAULT_QUERY_COVER,
        "--query-cover",
        "--min-query-cover",
        help="Minimum query coverage percentage.",
    ),
    min_subject_cover: float = typer.Option(
        DEFAULT_SUBJECT_COVER,
        "--subject-cover",
        "--min-subject-cover",
        help="Minimum subject coverage percentage.",
    ),
    evalue: float = typer.Option(DEFAULT_EVALUE, "--evalue", help="Maximum E-value."),
    threads: int = typer.Option(DEFAULT_THREADS, "--threads", "-p", min=1),
    sensitivity: str = typer.Option(DEFAULT_SENSITIVITY, "--sensitivity"),
    force: bool = typer.Option(False, "--force", help="Overwrite existing output files."),
    keep_intermediate: bool = typer.Option(
        False, "--keep-intermediate", help="Keep DIAMOND databases and raw tabular results."
    ),
    allow_unmapped: bool = typer.Option(
        False,
        "--allow-unmapped",
        help="Use protein accessions as gene IDs when no mapping is available.",
    ),
) -> None:
    """Run a two-way DIAMOND protein homology search."""

    try:
        summary = run_analysis(
            fasta1,
            fasta2,
            output=output,
            manifest=manifest,
            map1=map1,
            map2=map2,
            mode=mode,  # type: ignore[arg-type]
            min_identity=min_identity,
            min_query_cover=min_query_cover,
            min_subject_cover=min_subject_cover,
            evalue=evalue,
            threads=threads,
            sensitivity=sensitivity,
            force=force,
            keep_intermediate=keep_intermediate,
            allow_unmapped=allow_unmapped,
        )
    except HomologyError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(
        f"Wrote {summary.pair_count} homologous gene pair(s) to {summary.output}\n"
        f"Filtered hits: forward={summary.forward_hits}, reverse={summary.reverse_hits}\n"
        f"Metadata: {summary.metadata}"
    )


if __name__ == "__main__":
    app()
