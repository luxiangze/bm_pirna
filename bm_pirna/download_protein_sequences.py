#!/usr/bin/env python3
"""Download protein sequences from NCBI Datasets and FlyBase.

The command accepts one or more ID files or directories.  It detects common
NCBI and FlyBase identifiers from table columns or single-column files,
downloads the corresponding proteins, and writes FASTA plus an audit manifest.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from datetime import UTC, datetime
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import threading
import time
from typing import Any, Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen
import zipfile

from Bio import SeqIO
from loguru import logger
import typer

app = typer.Typer(add_completion=False)

DEFAULT_OUTPUT_DIR = Path("data/external/homology/protein_sequences")
DEFAULT_GENE_MAP = Path("data/external/gene_info_table.tsv")
FLYBASE_API_BASE = "https://api.flybase.org/api/v1.0"
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DATASET_TOOL = "bm_pirna_download_proteins"
INPUT_SUFFIXES = {".txt", ".tsv", ".csv"}
MISSING_VALUES = {"", "na", "n/a", "nan", "none", "null", ".", "-"}

FLYBASE_GENE_RE = re.compile(r"^FBgn\d+$", re.IGNORECASE)
FLYBASE_PROTEIN_RE = re.compile(r"^FBpp\d+$", re.IGNORECASE)
FLYBASE_TRANSCRIPT_RE = re.compile(r"^FBtr\d+$", re.IGNORECASE)
NCBI_PROTEIN_RE = re.compile(r"^(?:AP|NP|XP|YP|WP|ZP)_\d+(?:\.\d+)?$", re.IGNORECASE)
LOC_RE = re.compile(r"^LOC(\d+)$", re.IGNORECASE)
GENE_SYMBOL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$")
GENE_ID_RE = re.compile(r"\[GeneID=(\d+)\]", re.IGNORECASE)

HEADER_ALIASES = {
    "id": "id",
    "gene": "gene",
    "geneid": "gene_id",
    "gene_id": "gene_id",
    "geneidentifier": "gene_id",
    "genename": "gene_name",
    "gene_name": "gene_name",
    "genesymbol": "gene_symbol",
    "gene_symbol": "gene_symbol",
    "symbol": "gene_symbol",
    "proteinid": "protein_id",
    "protein_id": "protein_id",
    "proteinaccession": "protein_accession",
    "protein_accession": "protein_accession",
    "accession": "accession",
    "accessionversion": "accession",
    "locus": "gene_symbol",
    "locustag": "gene_symbol",
    "locus_tag": "gene_symbol",
}

COLUMN_PRIORITY = (
    "protein_id",
    "protein_accession",
    "accession",
    "gene_id",
    "gene",
    "gene_name",
    "gene_symbol",
    "id",
)


class DownloadError(RuntimeError):
    """Raised when a remote download cannot be completed."""


@dataclass
class ProteinSequence:
    """A protein sequence and its source metadata."""

    accession: str
    sequence: str
    description: str = ""
    source: str = ""
    gene_id: str = ""
    data_version: str = ""


@dataclass
class InputEntry:
    """One ID occurrence in an input file."""

    input_file: str
    row_number: int
    id_column: str
    raw_id: str
    normalized_id: str
    detected_type: str
    source: str
    request_kind: str = ""
    request_id: str = ""
    resolved_gene_id: str = ""
    status: str = "pending"
    message: str = ""
    proteins: list[ProteinSequence] = field(default_factory=list)


@dataclass
class ParsedInput:
    """Parsed IDs from one input file."""

    path: Path
    id_column: str
    entries: list[InputEntry]


@dataclass
class Resolution:
    """Result of resolving one canonical request."""

    status: str
    proteins: list[ProteinSequence] = field(default_factory=list)
    resolved_gene_id: str = ""
    message: str = ""
    data_version: str = ""


def normalize_column_name(value: str) -> str:
    """Normalize a column name for alias matching."""

    value = value.removeprefix("\ufeff").strip().lower()
    return re.sub(r"[^a-z0-9]+", "", value)


def normalize_value(value: str) -> str:
    """Remove CSV quoting and whitespace without changing accession case."""

    value = value.removeprefix("\ufeff").strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        value = value[1:-1].strip()
    if value.startswith(">"):
        value = value[1:].strip()
    return value.split()[0] if value.split() else ""


def normalize_gene_id(value: str) -> str:
    """Normalize a numeric NCBI GeneID while preserving the original in reports."""

    try:
        return str(int(value))
    except ValueError:
        return value


def _is_missing(value: str) -> bool:
    return value.strip().lower() in MISSING_VALUES


def _column_semantics(column_name: str) -> str:
    """Return the normalized semantic category of a selected column."""

    return HEADER_ALIASES.get(normalize_column_name(column_name), "")


def classify_id(
    value: str,
    *,
    column_name: str = "",
    id_type: str = "auto",
    source: str = "auto",
) -> tuple[str, str, str]:
    """Classify an ID as ``(detected_type, source, normalized_id)``.

    Numeric values default to NCBI GeneIDs.  Legacy numeric protein GIs are
    deliberately not guessed because they are ambiguous without explicit
    ``--id-type protein`` context.
    """

    normalized = normalize_value(value)
    if _is_missing(normalized):
        return "empty", "", normalized

    semantics = _column_semantics(column_name)
    explicit_type = id_type.lower()

    if FLYBASE_GENE_RE.fullmatch(normalized):
        detected_type, detected_source = "flybase_gene_id", "flybase"
    elif FLYBASE_PROTEIN_RE.fullmatch(normalized):
        detected_type, detected_source = "flybase_protein_id", "flybase"
    elif FLYBASE_TRANSCRIPT_RE.fullmatch(normalized):
        detected_type, detected_source = "unsupported_flybase_transcript", "flybase"
    elif NCBI_PROTEIN_RE.fullmatch(normalized):
        detected_type, detected_source = "ncbi_protein_accession", "ncbi"
    elif normalized.isdigit():
        detected_type, detected_source = "ncbi_gene_id", "ncbi"
        normalized = normalize_gene_id(normalized)
    elif (
        LOC_RE.fullmatch(normalized)
        or semantics
        in {
            "gene",
            "gene_id",
            "gene_name",
            "gene_symbol",
            "id",
        }
        or (column_name == "column_1" and GENE_SYMBOL_RE.fullmatch(normalized))
    ):
        detected_type, detected_source = "ncbi_gene_symbol", "ncbi"
    elif explicit_type == "protein" and semantics in {
        "protein_id",
        "protein_accession",
        "accession",
    }:
        detected_type, detected_source = "ncbi_protein_accession", "ncbi"
    else:
        detected_type, detected_source = "unrecognized", ""

    if explicit_type == "gene" and detected_type in {
        "ncbi_protein_accession",
        "flybase_protein_id",
    }:
        detected_type, detected_source = "unrecognized", ""
    elif explicit_type == "protein" and detected_type in {
        "ncbi_gene_id",
        "ncbi_gene_symbol",
        "flybase_gene_id",
    }:
        detected_type, detected_source = "unsupported_gene_protein", detected_source

    if source != "auto" and detected_source and detected_source != source.lower():
        detected_type = "source_mismatch"

    return detected_type, detected_source, normalized


def _detect_delimiter(lines: list[str], suffix: str) -> str:
    """Detect a delimiter using the extension and a small content sample."""

    if suffix.lower() == ".csv":
        return ","
    sample = [line for line in lines[:30] if line.strip()]
    if any("\t" in line for line in sample):
        return "\t"
    if any("," in line for line in sample):
        return ","
    return "\t"


def _structured_id_score(value: str) -> int:
    """Score IDs with unambiguous structural patterns for column inference."""

    value = normalize_value(value)
    if (
        FLYBASE_GENE_RE.fullmatch(value)
        or FLYBASE_PROTEIN_RE.fullmatch(value)
        or NCBI_PROTEIN_RE.fullmatch(value)
        or value.isdigit()
        or LOC_RE.fullmatch(value)
    ):
        return 1
    return 0


def _looks_like_header(rows: list[list[str]]) -> bool:
    """Detect a header without mistaking a first ID for one."""

    if not rows:
        return False
    first = rows[0]
    normalized = [normalize_column_name(cell) for cell in first]
    if any(cell in HEADER_ALIASES for cell in normalized):
        return True
    if len(first) == 1:
        return normalize_column_name(first[0]) in HEADER_ALIASES
    if len(first) > 1 and all(_structured_id_score(cell) == 0 for cell in first):
        return any(any(_structured_id_score(cell) > 0 for cell in row) for row in rows[1:])
    return False


def _select_id_column(
    rows: list[list[str]],
    *,
    has_header: bool,
    id_column: str | None,
) -> tuple[int, str, list[str]]:
    """Select an ID column and return ``(index, name, data_rows)``."""

    if not rows:
        raise ValueError("Input file contains no rows")

    header = rows[0] if has_header else []
    data_rows = rows[1:] if has_header else rows
    width = max(len(row) for row in rows)
    if id_column:
        requested = id_column.strip()
        try:
            index = int(requested) - 1
        except ValueError:
            normalized_requested = normalize_column_name(requested)
            matches = [
                i
                for i, name in enumerate(header)
                if normalize_column_name(name) == normalized_requested
            ]
            if not matches:
                raise ValueError(f"ID column {id_column!r} was not found")
            index = matches[0]
        if index < 0 or index >= width:
            raise ValueError(f"ID column index {id_column!r} is outside the table")
        name = header[index] if has_header and index < len(header) else f"column_{index + 1}"
        return index, name, data_rows

    if has_header:
        semantic_by_index = {
            i: HEADER_ALIASES.get(normalize_column_name(name), "") for i, name in enumerate(header)
        }
        for semantic in COLUMN_PRIORITY:
            matches = [i for i, value in semantic_by_index.items() if value == semantic]
            if matches:
                return matches[0], header[matches[0]], data_rows

    if width == 1:
        return 0, header[0] if has_header else "column_1", data_rows

    scores = [
        sum(_structured_id_score(row[i]) for row in data_rows if i < len(row))
        for i in range(width)
    ]
    best = max(scores, default=0)
    best_indices = [i for i, score in enumerate(scores) if score == best and score > 0]
    if len(best_indices) == 1:
        index = best_indices[0]
        name = header[index] if has_header and index < len(header) else f"column_{index + 1}"
        return index, name, data_rows
    raise ValueError(
        "Could not infer a unique ID column; use --id-column to select one explicitly"
    )


def parse_input_file(
    path: Path,
    *,
    id_column: str | None = None,
    id_type: str = "auto",
    source: str = "auto",
) -> ParsedInput:
    """Parse one ID file, including table headers and single-column lists."""

    text = path.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    delimiter = _detect_delimiter(lines, path.suffix)
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    rows = [row for row in rows if any(cell.strip() for cell in row)]
    if not rows:
        raise ValueError(f"Input file is empty: {path}")
    has_header = _looks_like_header(rows)
    index, selected_column, data_rows = _select_id_column(
        rows, has_header=has_header, id_column=id_column
    )
    entries: list[InputEntry] = []
    first_data_row = 2 if has_header else 1
    for offset, row in enumerate(data_rows):
        raw = row[index] if index < len(row) else ""
        detected_type, detected_source, normalized = classify_id(
            raw,
            column_name=selected_column,
            id_type=id_type,
            source=source,
        )
        entry = InputEntry(
            input_file=str(path),
            row_number=first_data_row + offset,
            id_column=selected_column,
            raw_id=raw.strip(),
            normalized_id=normalized,
            detected_type=detected_type,
            source=detected_source,
        )
        if detected_type == "empty":
            entry.status = "empty"
            entry.message = "Empty or missing ID"
        elif detected_type in {
            "unrecognized",
            "source_mismatch",
            "unsupported_flybase_transcript",
            "unsupported_numeric_protein",
            "unsupported_gene_protein",
        }:
            entry.status = "invalid"
            entry.message = f"Unsupported or ambiguous identifier: {normalized or raw}"
        elif detected_type == "ncbi_gene_id":
            entry.request_kind = "ncbi_gene"
            entry.request_id = normalized
        elif detected_type == "ncbi_protein_accession":
            entry.request_kind = "ncbi_protein"
            entry.request_id = normalized
        elif detected_type == "ncbi_gene_symbol":
            entry.request_kind = "ncbi_symbol"
            entry.request_id = normalized
        elif detected_type == "flybase_gene_id":
            entry.request_kind = "flybase_gene"
            entry.request_id = f"FBgn{normalized[4:]}"
        elif detected_type == "flybase_protein_id":
            entry.request_kind = "flybase_protein"
            entry.request_id = f"FBpp{normalized[4:]}"
        entries.append(entry)
    return ParsedInput(path=path, id_column=selected_column, entries=entries)


def collect_input_files(inputs: Iterable[Path]) -> list[Path]:
    """Expand files/directories into a deterministic list of supported files."""

    paths: list[Path] = []
    for input_path in inputs:
        if input_path.is_file():
            paths.append(input_path)
        elif input_path.is_dir():
            paths.extend(
                path
                for path in sorted(input_path.iterdir())
                if path.is_file() and path.suffix.lower() in INPUT_SUFFIXES
            )
        else:
            raise ValueError(f"Input path does not exist: {input_path}")
    unique = list(dict.fromkeys(path.resolve() for path in paths))
    if not unique:
        raise ValueError("No .txt, .tsv, or .csv input files were found")
    return unique


def load_gene_map(path: Path | None) -> dict[str, set[str]]:
    """Load a name-to-GeneID map from a tabular annotation file."""

    if path is None:
        return {}
    if not path.exists():
        raise ValueError(f"Gene map does not exist: {path}")
    text = path.read_text(encoding="utf-8-sig")
    delimiter = _detect_delimiter(text.splitlines(), path.suffix)
    rows = list(csv.DictReader(io.StringIO(text), delimiter=delimiter))
    if not rows:
        return {}
    fields = {normalize_column_name(name): name for name in rows[0]}
    id_field = next(
        (fields[key] for key in ("geneid", "gene_id", "entrezid", "entrez_id") if key in fields),
        None,
    )
    name_fields = [
        fields[key]
        for key in ("genename", "gene_name", "genesymbol", "gene_symbol", "symbol")
        if key in fields
    ]
    if not id_field or not name_fields:
        raise ValueError(
            f"Gene map must contain a GeneID column and a gene name/symbol column: {path}"
        )
    mapping: dict[str, set[str]] = {}
    for row in rows:
        gene_id = normalize_gene_id((row.get(id_field) or "").strip())
        if not gene_id.isdigit():
            continue
        for name_field in name_fields:
            name = normalize_value(row.get(name_field) or "")
            if name and not _is_missing(name):
                mapping.setdefault(name.casefold(), set()).add(gene_id)
    return mapping


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
        raise DownloadError(
            "The 'datasets' executable was not found; install ncbi-datasets-cli with pixi"
        ) from exc


def datasets_version() -> str:
    """Return the installed NCBI Datasets CLI version."""

    executable = shutil.which("datasets") or "datasets"
    result = _run_command([executable, "--version"])
    if result.returncode != 0:
        raise DownloadError(result.stderr.strip() or "Could not determine datasets version")
    return (result.stdout.strip() or result.stderr.strip()).splitlines()[0]


def _parse_summary_lines(stdout: str) -> list[dict[str, Any]]:
    """Parse JSON-lines output from ``datasets summary``."""

    records: list[dict[str, Any]] = []
    for line in stdout.splitlines():
        line = line.strip()
        if not line or not line.startswith("{"):
            continue
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            records.append(value)
    return records


def resolve_symbols_with_datasets(
    symbols: Iterable[str],
    *,
    taxid: str,
) -> dict[str, tuple[str | None, str]]:
    """Resolve symbols to unique NCBI GeneIDs using the requested taxon."""

    executable = shutil.which("datasets") or "datasets"
    resolved: dict[str, tuple[str | None, str]] = {}
    for symbol in dict.fromkeys(symbols):
        result = _run_command(
            [
                executable,
                "summary",
                "gene",
                "symbol",
                symbol,
                "--taxon",
                str(taxid),
                "--as-json-lines",
                "--limit",
                "all",
            ]
        )
        if result.returncode != 0:
            resolved[symbol] = (None, result.stderr.strip() or "NCBI symbol lookup failed")
            continue
        candidates: set[str] = set()
        for record in _parse_summary_lines(result.stdout):
            candidate_id = str(record.get("gene_id") or record.get("geneId") or "")
            if not candidate_id.isdigit():
                continue
            names = {str(record.get("symbol") or "").casefold()}
            synonyms = record.get("synonyms") or []
            if isinstance(synonyms, str):
                synonyms = [synonyms]
            names.update(str(name).casefold() for name in synonyms)
            if symbol.casefold() in names:
                candidates.add(normalize_gene_id(candidate_id))
        if len(candidates) == 1:
            resolved[symbol] = (next(iter(candidates)), "")
        elif not candidates:
            resolved[symbol] = (None, f"No unique NCBI GeneID found for symbol {symbol!r}")
        else:
            resolved[symbol] = (
                None,
                f"Ambiguous NCBI symbol {symbol!r}: {', '.join(sorted(candidates))}",
            )
    return resolved


def _zip_member(zip_file: zipfile.ZipFile, suffix: str) -> str | None:
    """Find an exact known data-package member by suffix."""

    for name in zip_file.namelist():
        normalized = name.replace("\\", "/")
        parts = normalized.split("/")
        if all(part not in {"", ".", ".."} for part in parts) and (
            normalized == suffix or parts[-2:] == ["data", suffix]
        ):
            return name
    return None


def _parse_ncbi_package(
    zip_path: Path,
) -> tuple[dict[str, list[ProteinSequence]], set[str], dict[str, str]]:
    """Parse proteins and gene metadata from an NCBI Datasets zip package."""

    proteins_by_gene: dict[str, list[ProteinSequence]] = {}
    proteins_by_accession: dict[str, list[ProteinSequence]] = {}
    gene_ids: set[str] = set()
    data_versions: dict[str, str] = {}
    with zipfile.ZipFile(zip_path) as archive:
        report_name = _zip_member(archive, "data_report.jsonl")
        if report_name:
            with archive.open(report_name) as raw:
                for line in io.TextIOWrapper(raw, encoding="utf-8"):
                    if not line.strip():
                        continue
                    try:
                        report = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    gene_id = str(report.get("geneId") or report.get("gene_id") or "")
                    if not gene_id.isdigit():
                        continue
                    gene_id = normalize_gene_id(gene_id)
                    gene_ids.add(gene_id)
                    annotations = report.get("annotations") or []
                    if annotations and isinstance(annotations[0], dict):
                        version = str(
                            annotations[0].get("annotationName")
                            or annotations[0].get("annotation_name")
                            or ""
                        )
                        if version:
                            data_versions[gene_id] = version

        protein_name = _zip_member(archive, "protein.faa")
        if protein_name:
            with archive.open(protein_name) as raw:
                for record in SeqIO.parse(io.TextIOWrapper(raw, encoding="utf-8"), "fasta"):
                    accession = record.id.split()[0]
                    gene_match = GENE_ID_RE.search(record.description)
                    gene_id = normalize_gene_id(gene_match.group(1)) if gene_match else ""
                    protein = ProteinSequence(
                        accession=accession,
                        sequence=str(record.seq).upper(),
                        description=record.description,
                        source="NCBI",
                        gene_id=gene_id,
                        data_version=data_versions.get(gene_id, ""),
                    )
                    proteins_by_accession.setdefault(accession.casefold(), []).append(protein)
                    if gene_id:
                        proteins_by_gene.setdefault(gene_id, []).append(protein)

    return proteins_by_gene, gene_ids, data_versions


def _run_ncbi_batches(
    request_ids: list[str],
    *,
    mode: str,
    batch_size: int,
    retries: int,
) -> tuple[dict[str, list[ProteinSequence]], set[str], dict[str, str], dict[str, str]]:
    """Download NCBI packages in batches and return proteins and errors."""

    all_gene_proteins: dict[str, list[ProteinSequence]] = {}
    all_accession_proteins: dict[str, list[ProteinSequence]] = {}
    all_gene_ids: set[str] = set()
    versions: dict[str, str] = {}
    errors: dict[str, str] = {}
    executable = shutil.which("datasets") or "datasets"
    unique_ids = list(dict.fromkeys(request_ids))
    with tempfile.TemporaryDirectory(prefix="bm_pirna_ncbi_") as temp_dir:
        temp_root = Path(temp_dir)
        for start in range(0, len(unique_ids), batch_size):
            chunk = unique_ids[start : start + batch_size]
            id_file = temp_root / f"ids_{start}.txt"
            id_file.write_text("\n".join(chunk) + "\n", encoding="utf-8")
            zip_path = temp_root / f"package_{start}.zip"
            command = [
                executable,
                "download",
                "gene",
                "gene-id" if mode == "gene" else "accession",
                "--inputfile",
                str(id_file),
                "--include",
                "protein",
                "--filename",
                str(zip_path),
                "--no-progressbar",
            ]
            if mode == "accession":
                command.extend(["--fasta-filter-file", str(id_file)])

            completed: subprocess.CompletedProcess[str] | None = None
            for attempt in range(max(1, retries)):
                zip_path.unlink(missing_ok=True)
                completed = _run_command(command)
                if completed.returncode == 0 and zip_path.exists():
                    break
                if attempt + 1 < max(1, retries):
                    time.sleep(min(8.0, 2**attempt))
            if completed is None or completed.returncode != 0 or not zip_path.exists():
                message = (
                    completed.stderr.strip()
                    if completed is not None and completed.stderr.strip()
                    else "NCBI Datasets download failed"
                )
                if len(chunk) > 1:
                    nested = _run_ncbi_batches(
                        chunk,
                        mode=mode,
                        batch_size=max(1, len(chunk) // 2),
                        retries=retries,
                    )
                    gene_map, gene_ids, nested_versions, nested_errors = nested
                    if mode == "accession":
                        all_accession_proteins.update(gene_map)
                    else:
                        all_gene_proteins.update(gene_map)
                    all_gene_ids.update(gene_ids)
                    versions.update(nested_versions)
                    errors.update(nested_errors)
                else:
                    errors[chunk[0]] = message
                continue

            try:
                gene_map, gene_ids, package_versions = _parse_ncbi_package(zip_path)
            except (OSError, ValueError, zipfile.BadZipFile) as exc:
                message = f"Could not parse NCBI package: {exc}"
                for request_id in chunk:
                    errors[request_id] = message
                continue
            all_gene_proteins.update(gene_map)
            all_gene_ids.update(gene_ids)
            versions.update(package_versions)
            if mode == "accession":
                with zipfile.ZipFile(zip_path) as archive:
                    protein_name = _zip_member(archive, "protein.faa")
                    if protein_name:
                        with archive.open(protein_name) as raw:
                            for record in SeqIO.parse(
                                io.TextIOWrapper(raw, encoding="utf-8"), "fasta"
                            ):
                                accession = record.id.split()[0]
                                gene_match = GENE_ID_RE.search(record.description)
                                gene_id = (
                                    normalize_gene_id(gene_match.group(1)) if gene_match else ""
                                )
                                protein = ProteinSequence(
                                    accession=accession,
                                    sequence=str(record.seq).upper(),
                                    description=record.description,
                                    source="NCBI",
                                    gene_id=gene_id,
                                    data_version=versions.get(gene_id, ""),
                                )
                                all_accession_proteins.setdefault(accession.casefold(), []).append(
                                    protein
                                )
    if mode == "accession":
        return all_accession_proteins, all_gene_ids, versions, errors
    return all_gene_proteins, all_gene_ids, versions, errors


def _fetch_ncbi_efetch_batch(
    accessions: list[str],
    *,
    timeout: float,
    retries: int,
    limiter: RateLimiter,
) -> tuple[dict[str, list[ProteinSequence]], str]:
    """Fetch a batch of protein FASTA records from NCBI EFetch."""

    params = {
        "db": "protein",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
        "tool": DATASET_TOOL,
    }
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    request = Request(
        f"{NCBI_EFETCH_URL}?{urlencode(params)}",
        headers={
            "Accept": "text/plain",
            "User-Agent": f"{DATASET_TOOL}/1.0",
        },
    )
    last_error = ""
    for attempt in range(max(1, retries)):
        limiter.wait()
        try:
            with urlopen(request, timeout=timeout) as response:
                text = response.read().decode("utf-8")
            if not text.strip():
                raise DownloadError("NCBI EFetch returned an empty response")
            proteins: dict[str, list[ProteinSequence]] = {}
            for record in SeqIO.parse(io.StringIO(text), "fasta"):
                accession = record.id.split()[0]
                gene_match = GENE_ID_RE.search(record.description)
                gene_id = normalize_gene_id(gene_match.group(1)) if gene_match else ""
                protein = ProteinSequence(
                    accession=accession,
                    sequence=str(record.seq).upper(),
                    description=record.description,
                    source="NCBI",
                    gene_id=gene_id,
                    data_version="NCBI EFetch",
                )
                proteins.setdefault(accession.casefold(), []).append(protein)
            if not proteins:
                raise DownloadError("NCBI EFetch returned no FASTA records")
            return proteins, ""
        except HTTPError as exc:
            last_error = f"HTTP {exc.code}"
            if exc.code not in {408, 425, 429} and not 500 <= exc.code <= 599:
                break
        except (URLError, TimeoutError, OSError, ValueError, DownloadError) as exc:
            last_error = str(exc)
        if attempt + 1 < max(1, retries):
            time.sleep(min(8.0, 2**attempt))
    return {}, last_error or "NCBI EFetch request failed"


def _fetch_ncbi_proteins_efetch(
    accessions: Iterable[str],
    *,
    timeout: float,
    retries: int,
    batch_size: int = 100,
    limiter: RateLimiter | None = None,
) -> tuple[dict[str, list[ProteinSequence]], dict[str, str]]:
    """Fetch direct protein accessions that Datasets no longer resolves."""

    limiter = limiter or RateLimiter()
    unique = list(dict.fromkeys(accessions))
    proteins: dict[str, list[ProteinSequence]] = {}
    errors: dict[str, str] = {}

    def fetch_chunk(chunk: list[str]) -> None:
        fetched, error = _fetch_ncbi_efetch_batch(
            chunk,
            timeout=timeout,
            retries=retries,
            limiter=limiter,
        )
        if not error:
            for key, values in fetched.items():
                proteins.setdefault(key, []).extend(values)
            return
        if len(chunk) > 1:
            midpoint = max(1, len(chunk) // 2)
            fetch_chunk(chunk[:midpoint])
            fetch_chunk(chunk[midpoint:])
            return
        errors[chunk[0]] = error

    for start in range(0, len(unique), max(1, batch_size)):
        fetch_chunk(unique[start : start + max(1, batch_size)])
    return proteins, errors


def _choose_proteins(proteins: list[ProteinSequence], policy: str) -> list[ProteinSequence]:
    """Select all or the deterministic longest protein(s)."""

    unique: dict[str, ProteinSequence] = {}
    for protein in proteins:
        unique.setdefault(protein.accession.casefold(), protein)
    ordered = sorted(unique.values(), key=lambda item: item.accession)
    if policy == "all":
        return ordered
    if not ordered:
        return []
    return [max(ordered, key=lambda item: len(item.sequence))]


class RateLimiter:
    """Simple monotonic rate limiter for FlyBase requests."""

    def __init__(self, minimum_interval: float = 0.34) -> None:
        self.minimum_interval = minimum_interval
        self.last_request = 0.0
        self._lock = threading.Lock()

    def wait(self) -> None:
        with self._lock:
            now = time.monotonic()
            delay = self.minimum_interval - (now - self.last_request)
            if delay > 0:
                time.sleep(delay)
            self.last_request = time.monotonic()


def _request_flybase(
    flybase_id: str,
    *,
    subtype: str | None,
    timeout: float,
    retries: int,
    limiter: RateLimiter,
) -> dict[str, Any]:
    """Fetch one FlyBase sequence response with bounded retries."""

    path = f"/sequence/id/{quote(flybase_id, safe='')}"
    if subtype:
        path += f"/{quote(subtype, safe='')}"
    url = f"{FLYBASE_API_BASE}{path}"
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": f"{DATASET_TOOL}/1.0",
        },
    )
    last_error = ""
    for attempt in range(max(1, retries)):
        limiter.wait()
        try:
            with urlopen(request, timeout=timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
            if not isinstance(payload, dict):
                raise DownloadError("FlyBase returned a non-object JSON response")
            return payload
        except HTTPError as exc:
            last_error = f"HTTP {exc.code}"
            if exc.code == 404:
                return {"not_found": [flybase_id], "data_version": ""}
            if exc.code not in {408, 425, 429} and not 500 <= exc.code <= 599:
                break
        except (URLError, TimeoutError, OSError, json.JSONDecodeError, DownloadError) as exc:
            last_error = str(exc)
        if attempt + 1 < max(1, retries):
            time.sleep(min(8.0, 2**attempt))
    raise DownloadError(f"FlyBase request failed for {flybase_id}: {last_error}")


def _proteins_from_flybase_payload(
    payload: dict[str, Any],
    *,
    source_gene_id: str = "",
) -> tuple[list[ProteinSequence], str, str]:
    """Extract protein records, data version, and a status message."""

    resultset = payload.get("resultset")
    if not isinstance(resultset, dict):
        resultset = payload
    data_version = str(resultset.get("data_version") or payload.get("data_version") or "")
    result = resultset.get("result") or []
    proteins: list[ProteinSequence] = []
    if isinstance(result, list):
        for item in result:
            if not isinstance(item, dict):
                continue
            accession = str(item.get("id") or "").strip()
            sequence = "".join(str(item.get("sequence") or "").split()).upper()
            if not accession or not sequence:
                continue
            proteins.append(
                ProteinSequence(
                    accession=accession,
                    sequence=sequence,
                    description=str(item.get("description") or "").strip(),
                    source="FlyBase",
                    gene_id=source_gene_id,
                    data_version=data_version,
                )
            )
    if proteins:
        return proteins, data_version, ""
    not_found = resultset.get("not_found") or payload.get("not_found") or []
    if not_found:
        return [], data_version, f"FlyBase did not find: {', '.join(map(str, not_found))}"
    return [], data_version, "No protein sequence was returned"


def fetch_flybase_resolutions(
    entries: Iterable[InputEntry],
    *,
    timeout: float,
    retries: int,
    limiter: RateLimiter | None = None,
) -> dict[tuple[str, str], Resolution]:
    """Fetch unique FlyBase gene/protein requests."""

    limiter = limiter or RateLimiter()
    resolutions: dict[tuple[str, str], Resolution] = {}
    requests = {
        (entry.request_kind, entry.request_id)
        for entry in entries
        if entry.status == "pending" and entry.source == "flybase"
    }

    def fetch_one(key: tuple[str, str]) -> tuple[tuple[str, str], Resolution]:
        kind, request_id = key
        try:
            payload = _request_flybase(
                request_id,
                subtype="FBpp" if kind == "flybase_gene" else None,
                timeout=timeout,
                retries=retries,
                limiter=limiter,
            )
            proteins, version, message = _proteins_from_flybase_payload(
                payload,
                source_gene_id=request_id if kind == "flybase_gene" else "",
            )
            status = "downloaded" if proteins else "no_protein"
            resolution = Resolution(
                status=status,
                proteins=proteins,
                resolved_gene_id=request_id if kind == "flybase_gene" else "",
                message=message,
                data_version=version,
            )
        except DownloadError as exc:
            resolution = Resolution(
                status="download_error",
                message=str(exc),
            )
        return key, resolution

    for key in sorted(requests):
        key, resolution = fetch_one(key)
        resolutions[key] = resolution
    return resolutions


def _resolve_symbol_entries(
    entries: Iterable[InputEntry],
    *,
    gene_map: dict[str, set[str]],
    taxid: str | None,
) -> None:
    """Convert gene symbols to canonical NCBI GeneID requests."""

    unresolved_symbols: list[str] = []
    for entry in entries:
        if entry.status != "pending" or entry.request_kind != "ncbi_symbol":
            continue
        mapped = gene_map.get(entry.request_id.casefold(), set())
        if len(mapped) == 1:
            entry.request_kind = "ncbi_gene"
            entry.request_id = next(iter(mapped))
            entry.resolved_gene_id = entry.request_id
        elif len(mapped) > 1:
            entry.status = "ambiguous"
            entry.message = f"Gene map contains multiple GeneIDs: {', '.join(sorted(mapped))}"
        else:
            unresolved_symbols.append(entry.request_id)

    if not unresolved_symbols:
        return
    if not taxid:
        for symbol in unresolved_symbols:
            for entry in entries:
                if entry.status == "pending" and entry.request_id == symbol:
                    entry.status = "ambiguous"
                    entry.message = "A taxid is required to resolve a gene symbol"
        return
    resolutions = resolve_symbols_with_datasets(unresolved_symbols, taxid=taxid)
    for symbol in unresolved_symbols:
        gene_id, message = resolutions.get(symbol, (None, "Symbol lookup returned no result"))
        for entry in entries:
            if entry.status != "pending" or entry.request_id != symbol:
                continue
            if gene_id:
                entry.request_kind = "ncbi_gene"
                entry.request_id = gene_id
                entry.resolved_gene_id = gene_id
            else:
                entry.status = "ambiguous"
                entry.message = message


def _apply_resolutions(
    entries: Iterable[InputEntry],
    resolutions: dict[tuple[str, str], Resolution],
    *,
    isoform: str,
) -> None:
    """Attach downloaded proteins and statuses to input occurrences."""

    for entry in entries:
        if entry.status != "pending":
            continue
        key = (entry.request_kind, entry.request_id)
        resolution = resolutions.get(key)
        if resolution is None:
            entry.status = "download_error"
            entry.message = "No resolution was produced for this ID"
            continue
        entry.status = resolution.status
        entry.message = resolution.message
        entry.resolved_gene_id = entry.resolved_gene_id or resolution.resolved_gene_id
        if resolution.data_version:
            for protein in resolution.proteins:
                protein.data_version = protein.data_version or resolution.data_version
        entry.proteins = _choose_proteins(resolution.proteins, isoform)
        if entry.status == "downloaded" and not entry.proteins:
            entry.status = "no_protein"
            entry.message = "No protein remained after deduplication"


def _manifest_rows(entries: Iterable[InputEntry]) -> list[dict[str, str]]:
    """Expand one-to-many protein mappings into manifest rows."""

    rows: list[dict[str, str]] = []
    fieldnames = [
        "input_file",
        "row_number",
        "id_column",
        "raw_id",
        "normalized_id",
        "detected_type",
        "source",
        "resolved_gene_id",
        "protein_accession",
        "sequence_length",
        "status",
        "message",
        "data_version",
    ]
    for entry in entries:
        proteins = entry.proteins or [None]
        for protein in proteins:
            row = {name: "" for name in fieldnames}
            row.update(
                {
                    "input_file": entry.input_file,
                    "row_number": str(entry.row_number),
                    "id_column": entry.id_column,
                    "raw_id": entry.raw_id,
                    "normalized_id": entry.normalized_id,
                    "detected_type": entry.detected_type,
                    "source": entry.source,
                    "resolved_gene_id": entry.resolved_gene_id
                    or (protein.gene_id if protein else ""),
                    "status": entry.status,
                    "message": entry.message,
                }
            )
            if protein:
                row["protein_accession"] = protein.accession
                row["sequence_length"] = str(len(protein.sequence))
                row["data_version"] = protein.data_version
            rows.append(row)
    return rows


def _write_fasta(path: Path, proteins: Iterable[ProteinSequence]) -> None:
    """Write deterministic FASTA records with wrapped sequences."""

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for protein in proteins:
            description = protein.description.strip()
            if description and not description.startswith(protein.accession):
                description = f"{protein.accession} {description}"
            header = description or protein.accession
            handle.write(f">{header}\n")
            sequence = protein.sequence.upper()
            for start in range(0, len(sequence), 80):
                handle.write(f"{sequence[start : start + 80]}\n")


def _write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    """Write the audit manifest as UTF-8 TSV."""

    fieldnames = [
        "input_file",
        "row_number",
        "id_column",
        "raw_id",
        "normalized_id",
        "detected_type",
        "source",
        "resolved_gene_id",
        "protein_accession",
        "sequence_length",
        "status",
        "message",
        "data_version",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _unique_proteins(entries: Iterable[InputEntry]) -> list[ProteinSequence]:
    """Collect unique proteins in first-seen order."""

    proteins: dict[str, ProteinSequence] = {}
    for entry in entries:
        for protein in entry.proteins:
            proteins.setdefault(protein.accession.casefold(), protein)
    return list(proteins.values())


def _output_paths(output_dir: Path, parsed_inputs: list[ParsedInput]) -> list[Path]:
    """Return all files that the command may overwrite."""

    return [
        *(output_dir / f"{parsed.path.stem}.proteins.faa" for parsed in parsed_inputs),
        output_dir / "all.proteins.faa",
        output_dir / "protein_download_manifest.tsv",
        output_dir / "download_metadata.json",
    ]


def run_download(
    inputs: Iterable[Path],
    *,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    id_column: str | None = None,
    id_type: str = "auto",
    source: str = "auto",
    taxid: str | None = None,
    gene_map_path: Path | None = None,
    isoform: str = "longest",
    batch_size: int = 200,
    retries: int = 3,
    timeout: float = 30.0,
    force: bool = False,
    allow_partial: bool = False,
) -> dict[str, Any]:
    """Run parsing, remote retrieval, and output writing."""

    if id_type not in {"auto", "gene", "protein"}:
        raise ValueError("id_type must be auto, gene, or protein")
    if source not in {"auto", "ncbi", "flybase"}:
        raise ValueError("source must be auto, ncbi, or flybase")
    if isoform not in {"longest", "all"}:
        raise ValueError("isoform must be longest or all")
    if batch_size < 1 or retries < 1 or timeout <= 0:
        raise ValueError("batch_size, retries, and timeout must be positive")

    input_paths = collect_input_files(inputs)
    parsed_inputs = [
        parse_input_file(
            path,
            id_column=id_column,
            id_type=id_type,
            source=source,
        )
        for path in input_paths
    ]
    stems = [parsed.path.stem for parsed in parsed_inputs]
    if len(stems) != len(set(stems)):
        raise ValueError("Input files have duplicate stems; use distinct filenames")
    all_entries = [entry for parsed in parsed_inputs for entry in parsed.entries]
    if gene_map_path is None and DEFAULT_GENE_MAP.exists():
        gene_map_path = DEFAULT_GENE_MAP
    gene_map = load_gene_map(gene_map_path)
    _resolve_symbol_entries(all_entries, gene_map=gene_map, taxid=taxid)

    resolutions: dict[tuple[str, str], Resolution] = {}
    ncbi_gene_ids = sorted(
        {
            (entry.request_id)
            for entry in all_entries
            if entry.status == "pending" and entry.request_kind == "ncbi_gene"
        }
    )
    ncbi_protein_ids = sorted(
        {
            entry.request_id
            for entry in all_entries
            if entry.status == "pending" and entry.request_kind == "ncbi_protein"
        }
    )
    if ncbi_gene_ids and source != "flybase":
        try:
            gene_proteins, returned_gene_ids, versions, errors = _run_ncbi_batches(
                ncbi_gene_ids,
                mode="gene",
                batch_size=batch_size,
                retries=retries,
            )
            for gene_id in ncbi_gene_ids:
                if gene_id in errors:
                    resolutions[("ncbi_gene", gene_id)] = Resolution(
                        status="download_error", message=errors[gene_id]
                    )
                elif gene_id in gene_proteins:
                    resolutions[("ncbi_gene", gene_id)] = Resolution(
                        status="downloaded",
                        proteins=gene_proteins[gene_id],
                        resolved_gene_id=gene_id,
                        data_version=versions.get(gene_id, ""),
                    )
                elif gene_id in returned_gene_ids:
                    resolutions[("ncbi_gene", gene_id)] = Resolution(
                        status="no_protein",
                        resolved_gene_id=gene_id,
                        message="NCBI Gene record has no protein.faa member",
                        data_version=versions.get(gene_id, ""),
                    )
                else:
                    resolutions[("ncbi_gene", gene_id)] = Resolution(
                        status="not_found", message="NCBI GeneID was not returned"
                    )
        except DownloadError as exc:
            for gene_id in ncbi_gene_ids:
                resolutions[("ncbi_gene", gene_id)] = Resolution(
                    status="download_error", message=str(exc)
                )
    if ncbi_protein_ids and source != "flybase":
        try:
            protein_map, _, versions, errors = _run_ncbi_batches(
                ncbi_protein_ids,
                mode="accession",
                batch_size=batch_size,
                retries=retries,
            )

            def accession_candidates(accession: str) -> list[ProteinSequence]:
                candidates = protein_map.get(accession.casefold(), [])
                if candidates:
                    return candidates
                base = accession.split(".", 1)[0].casefold()
                return [
                    protein
                    for key, values in protein_map.items()
                    if key.split(".", 1)[0] == base
                    for protein in values
                ]

            fallback_ids = [
                accession for accession in ncbi_protein_ids if not accession_candidates(accession)
            ]
            if fallback_ids:
                fallback_map, fallback_errors = _fetch_ncbi_proteins_efetch(
                    fallback_ids,
                    timeout=timeout,
                    retries=retries,
                    limiter=RateLimiter(),
                )
                for key, values in fallback_map.items():
                    protein_map.setdefault(key, []).extend(values)
                for accession, message in fallback_errors.items():
                    errors.setdefault(accession, message)

            for accession in ncbi_protein_ids:
                candidates = accession_candidates(accession)
                if candidates:
                    resolutions[("ncbi_protein", accession)] = Resolution(
                        status="downloaded",
                        proteins=_choose_proteins(candidates, "all"),
                        resolved_gene_id=candidates[0].gene_id,
                    )
                elif accession in errors:
                    resolutions[("ncbi_protein", accession)] = Resolution(
                        status="download_error", message=errors[accession]
                    )
                else:
                    resolutions[("ncbi_protein", accession)] = Resolution(
                        status="not_found", message="NCBI protein accession was not returned"
                    )
        except DownloadError as exc:
            for accession in ncbi_protein_ids:
                resolutions[("ncbi_protein", accession)] = Resolution(
                    status="download_error", message=str(exc)
                )

    flybase_entries = [
        entry for entry in all_entries if entry.status == "pending" and entry.source == "flybase"
    ]
    if flybase_entries and source != "ncbi":
        resolutions.update(
            fetch_flybase_resolutions(
                flybase_entries,
                timeout=timeout,
                retries=retries,
            )
        )
    _apply_resolutions(all_entries, resolutions, isoform=isoform)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = _output_paths(output_dir, parsed_inputs)
    if not force:
        existing = [path for path in output_paths if path.exists()]
        if existing:
            raise ValueError(
                "Output exists; use --force to overwrite: "
                + ", ".join(str(path) for path in existing)
            )
    for parsed in parsed_inputs:
        _write_fasta(
            output_dir / f"{parsed.path.stem}.proteins.faa", _unique_proteins(parsed.entries)
        )
    _write_fasta(output_dir / "all.proteins.faa", _unique_proteins(all_entries))
    manifest_rows = _manifest_rows(all_entries)
    _write_manifest(output_dir / "protein_download_manifest.tsv", manifest_rows)
    fatal_statuses = {"empty", "invalid", "ambiguous", "download_error", "not_found"}
    fatal_count = sum(entry.status in fatal_statuses for entry in all_entries)
    flybase_data_versions = sorted(
        {
            protein.data_version
            for entry in all_entries
            for protein in entry.proteins
            if protein.source == "FlyBase" and protein.data_version
        }
    )
    metadata = {
        "created_at": datetime.now(UTC).isoformat(),
        "isoform_policy": isoform,
        "taxid": taxid,
        "gene_map": str(gene_map_path) if gene_map_path else None,
        "inputs": [str(path) for path in input_paths],
        "counts": {
            "input_rows": len(all_entries),
            "manifest_rows": len(manifest_rows),
            "downloaded_rows": sum(entry.status == "downloaded" for entry in all_entries),
            "no_protein_rows": sum(entry.status == "no_protein" for entry in all_entries),
            "fatal_rows": fatal_count,
            "unique_proteins": len(_unique_proteins(all_entries)),
        },
        "flybase_api": FLYBASE_API_BASE,
        "flybase_data_versions": flybase_data_versions,
        "ncbi_datasets_version": datasets_version() if ncbi_gene_ids or ncbi_protein_ids else "",
    }
    (output_dir / "download_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    logger.info(
        "Processed {} rows, wrote {} unique proteins, {} no-protein rows, {} fatal rows",
        len(all_entries),
        metadata["counts"]["unique_proteins"],
        metadata["counts"]["no_protein_rows"],
        fatal_count,
    )
    if fatal_count and not allow_partial:
        raise DownloadError(
            f"{fatal_count} input rows could not be resolved; see protein_download_manifest.tsv"
        )
    return metadata


@app.command()
def main(
    inputs: list[Path] = typer.Argument(..., help="Input ID files or directories"),
    output_dir: Path = typer.Option(
        DEFAULT_OUTPUT_DIR,
        "--output-dir",
        "-o",
        help="Directory for FASTA and manifest outputs",
    ),
    id_column: str | None = typer.Option(
        None,
        "--id-column",
        help="ID column name or one-based column index; infer automatically by default",
    ),
    id_type: str = typer.Option(
        "auto",
        "--id-type",
        help="ID interpretation: auto, gene, or protein",
    ),
    source: str = typer.Option(
        "auto",
        "--source",
        help="Database routing: auto, ncbi, or flybase",
    ),
    taxid: str | None = typer.Option(
        None,
        "--taxid",
        help="NCBI taxonomy ID used for gene-symbol resolution",
    ),
    gene_map: Path | None = typer.Option(
        None,
        "--gene-map",
        help="Optional table containing gene_name/symbol and gene_id columns",
    ),
    isoform: str = typer.Option(
        "longest",
        "--isoform",
        help="Protein selection per gene: longest or all",
    ),
    batch_size: int = typer.Option(
        200,
        "--batch-size",
        help="Number of NCBI IDs per data package request",
    ),
    retries: int = typer.Option(
        3,
        "--retries",
        help="Retries for transient NCBI/FlyBase failures",
    ),
    timeout: float = typer.Option(
        30.0,
        "--timeout",
        help="FlyBase request timeout in seconds",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Overwrite existing output files",
    ),
    allow_partial: bool = typer.Option(
        False,
        "--allow-partial",
        help="Exit successfully even when IDs remain unresolved",
    ),
) -> None:
    """Download proteins from mixed NCBI and FlyBase ID files."""

    try:
        run_download(
            inputs,
            output_dir=output_dir,
            id_column=id_column,
            id_type=id_type,
            source=source,
            taxid=taxid,
            gene_map_path=gene_map,
            isoform=isoform,
            batch_size=batch_size,
            retries=retries,
            timeout=timeout,
            force=force,
            allow_partial=allow_partial,
        )
    except (DownloadError, ValueError) as exc:
        logger.error(str(exc))
        raise typer.Exit(1) from exc


if __name__ == "__main__":
    app()
