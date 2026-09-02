"""Extract Bombyx mori x Drosophila melanogaster ortholog pairs from OrthoDB v12.

Performs a single-pass scan of the 22 GB ``odb12v2_OG2genes.tab`` to find
ortholog groups (OGs) containing both *Bombyx mori* (NCBI taxid 7091) and
*Drosophila melanogaster* (taxid 7227), then maps ODB internal gene IDs to
NCBI GeneIDs (Bombyx) and FlyBase FBgn IDs (Drosophila) via
``odb12v2_genes.tab``.  A FlyBase annotation-ID correspondence table provides
a fallback CG -> FBgn mapping when the Ensembl/FlyBase column is empty.

Outputs (``data/processed/homology/orthodb_orthologs/``):
  - ``bm_dm_ortholog_pairs_unique.tsv``     Unique (bm_gene_id, dm_fbgn) pairs
  - ``bm_ortholog_genes_unique.tsv``        Unique Bombyx genes with orthologs
  - ``dm_ortholog_fbgns_unique.tsv``        Unique Drosophila FBgn with orthologs
  - ``subset_cross_species_intersections.tsv``     Cross-species subset intersections
  - ``subset_cross_species_intersection_counts.tsv``   Per-pair summary counts
  - ``ortholog_extraction_summary.json``    Run metadata and counts

Usage::

    pixi run python -m bm_pirna.extract_orthodb_orthologs [--force]
"""

from collections import defaultdict
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import time

from loguru import logger
import typer

app = typer.Typer(add_completion=False, pretty_exceptions_enable=False)

# --- project paths (resolved relative to the package, never hard-coded) ---
PROJ_ROOT = Path(__file__).resolve().parents[1]
OG2GENES_PATH = PROJ_ROOT / "data/external/odb12v2_OG2genes.tab"
GENES_TAB_PATH = PROJ_ROOT / "data/external/odb12v2_genes.tab"
FBGN_MAP_PATH = PROJ_ROOT / "data/external/dm/fbgn_annotation_ID_fb_2026_02.tsv"
SUBSET_DIR = PROJ_ROOT / "data/external/homology/differentially_alternatively_spliced_genes"
OUTPUT_DIR = PROJ_ROOT / "data/processed/homology/orthodb_orthologs"

# --- taxonomy constants ---
BM_TAXID = "7091"  # Bombyx mori
DM_TAXID = "7227"  # Drosophila melanogaster

# --- gene subset registry: (name, file, id_column, species, id_format) ---
_SUBSETS: list[tuple[str, Path, str | None, str, str]] = [
    (
        "bm_AS_genes_1_control_vs_sugp1",
        SUBSET_DIR / "bm_AS_genes_1_control_vs_sugp1.txt",
        "GeneID",
        "bm",
        "geneid",
    ),
    (
        "bm_AS_genes_2_control_vs_sugp1",
        SUBSET_DIR / "bm_AS_genes_2_control_vs_sugp1.txt",
        "gene_id",
        "bm",
        "geneid",
    ),
    (
        "bm_AS_genes_2_p200_vs_sugp1",
        SUBSET_DIR / "bm_AS_genes_2_p200_vs_sugp1.txt",
        "gene_id",
        "bm",
        "geneid",
    ),
    (
        "dm_AS_genes_control_vs_sugp1",
        SUBSET_DIR / "dm_AS_genes_control_vs_sugp1.txt",
        "GeneID",
        "dm",
        "fbgn",
    ),
    (
        "IP-MS_sugp1",
        SUBSET_DIR / "IP-MS_sugp1_gene_ids.txt",
        "GeneID",
        "bm",
        "geneid",
    ),
]

_BASE_COLUMNS = [
    "og_id",
    "og_level_taxid",
    "pair_type",
    "bm_gene_id",
    "bm_odb_gene_id",
    "dm_fbgn",
    "dm_gene_symbol",
    "dm_odb_gene_id",
]

_UNIQUE_PAIR_COLUMNS = [
    "bm_gene_id",
    "dm_fbgn",
    "dm_gene_symbol",
    "og_ids",
    "og_level_taxids",
    "pair_types",
    "n_ogs",
]

_UNIQUE_BM_COLUMNS = [
    "bm_gene_id",
    "n_ortholog_fbgns",
    "ortholog_dm_fbgns",
    "n_ogs",
    "og_ids",
]

_UNIQUE_DM_COLUMNS = [
    "dm_fbgn",
    "dm_gene_symbol",
    "n_ortholog_bm_genes",
    "ortholog_bm_gene_ids",
    "n_ogs",
    "og_ids",
]


# ---------------------------------------------------------------------------
# Large-file scanning (awk / grep for speed on 22 GB files)
# ---------------------------------------------------------------------------


def _run_awk(script: str, file_path: Path) -> str:
    """Run awk on a large file with LC_ALL=C, return captured stdout."""
    env = {**os.environ, "LC_ALL": "C"}
    result = subprocess.run(
        ["awk", "-F", "\\t", script, str(file_path)],
        capture_output=True,
        text=True,
        env=env,
        check=True,
    )
    return result.stdout


def scan_og2genes(path: Path) -> dict[str, tuple[list[str], list[str]]]:
    """Single-pass scan of OG2genes.tab.

    Returns ``{og_id: (bm_odb_ids, dm_odb_ids)}`` for OGs that contain both
    Bombyx mori (7091) and Drosophila melanogaster (7227).
    """
    logger.info(f"Scanning {path.name} (~22 GB, expect ~5 min)...")
    t0 = time.time()
    script = (
        '{ og=$1; split($2,a,"_"); tax=a[1]; '
        'if (tax==7091) bm[og]=bm[og] (bm[og]?",":"")$2; '
        'else if (tax==7227) dm[og]=dm[og] (dm[og]?",":"")$2; } '
        'END { for (o in bm) if (o in dm) print o "\\t" bm[o] "\\t" dm[o]; }'
    )
    stdout = _run_awk(script, path)
    og_pairs: dict[str, tuple[list[str], list[str]]] = {}
    for line in stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        og = parts[0]
        bm_ids = parts[1].split(",") if len(parts) > 1 and parts[1] else []
        dm_ids = parts[2].split(",") if len(parts) > 2 and parts[2] else []
        og_pairs[og] = (bm_ids, dm_ids)
    logger.info(
        f"  Found {len(og_pairs):,} OGs with both species "
        f"in {time.time() - t0:.0f}s"
    )
    return og_pairs


def scan_genes_tab(path: Path) -> tuple[dict[str, dict], dict[str, dict]]:
    """Stream rows starting with ``7091_`` or ``7227_`` from genes.tab.

    Returns ``(bm_map, dm_map)`` where each maps ``odb_gene_id`` to a dict
    with the resolved public gene ID, symbol and description.
    """
    logger.info(
        f"Scanning {path.name} for {BM_TAXID}/{DM_TAXID} rows "
        f"(~21 GB, expect ~5 min)..."
    )
    t0 = time.time()
    env = {**os.environ, "LC_ALL": "C"}
    proc = subprocess.Popen(
        ["grep", "-a", "-E", f"^({BM_TAXID}|{DM_TAXID})_", str(path)],
        stdout=subprocess.PIPE,
        text=True,
        env=env,
    )
    bm_map: dict[str, dict] = {}
    dm_map: dict[str, dict] = {}
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 8:
            continue
        odb_gene_id = parts[0]
        synonyms = parts[3]  # col 4: Bombyx GeneID / Drosophila symbol or CG
        ensembl_col = parts[5]  # col 6: Drosophila FBgn / Bombyx LOC name
        description = parts[7]  # col 8
        tax = odb_gene_id.split("_", 1)[0]
        if tax == BM_TAXID:
            # Bombyx: col 4 is usually the pure-number NCBI GeneID
            gene_id = synonyms if synonyms.isdigit() else ""
            # Fallback: extract from LOC / GeneID_ prefix in col 6
            if not gene_id:
                if ensembl_col.startswith("LOC") and ensembl_col[3:].isdigit():
                    gene_id = ensembl_col[3:]
                elif ensembl_col.startswith("GeneID_") and ensembl_col[7:].isdigit():
                    gene_id = ensembl_col[7:]
            bm_map[odb_gene_id] = {
                "gene_id": gene_id,
                "symbol": ensembl_col,
                "description": description,
            }
        elif tax == DM_TAXID:
            # Drosophila: col 6 is usually FBgn, col 4 is symbol or CG number
            fbgn = ensembl_col if ensembl_col.startswith("FBgn") else ""
            dm_map[odb_gene_id] = {
                "fbgn": fbgn,
                "symbol": synonyms,
                "description": description,
            }
    proc.wait()
    if proc.returncode not in (0, 1):  # grep returns 1 when no matches
        raise RuntimeError(f"grep failed with code {proc.returncode}")
    logger.info(
        f"  Loaded {len(bm_map):,} Bombyx and {len(dm_map):,} Drosophila "
        f"gene mappings in {time.time() - t0:.0f}s"
    )
    return bm_map, dm_map


# ---------------------------------------------------------------------------
# FlyBase CG -> FBgn fallback mapping
# ---------------------------------------------------------------------------


def load_fbgn_map(path: Path) -> dict[str, str]:
    """Load FlyBase annotation-ID correspondence table.

    Returns a dict that maps both annotation IDs (CG####) and gene symbols
    to primary FBgn IDs.  Annotation IDs take priority over symbols when
    both map to different FBgns.
    """
    cg_to_fbgn: dict[str, str] = {}
    symbol_to_fbgn: dict[str, str] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("##") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            symbol, fbgn, annotation_id = parts[0], parts[2], parts[4]
            if not fbgn.startswith("FBgn"):
                continue
            if annotation_id:
                cg_to_fbgn[annotation_id] = fbgn
            if symbol:
                symbol_to_fbgn[symbol] = fbgn
    # CG numbers override symbols
    return {**symbol_to_fbgn, **cg_to_fbgn}


# ---------------------------------------------------------------------------
# Gene subset loading
# ---------------------------------------------------------------------------


def load_gene_subsets() -> dict[str, set[str]]:
    """Load all gene-subset ID sets defined in ``_SUBSETS``."""
    subsets: dict[str, set[str]] = {}
    for name, path, id_column, _species, id_format in _SUBSETS:
        ids = _load_one_subset(path, id_column, id_format)
        subsets[name] = ids
        logger.info(f"  Subset '{name}': {len(ids):,} IDs")
    return subsets


def _load_one_subset(
    path: Path, id_column: str | None, id_format: str
) -> set[str]:
    """Load gene IDs from a subset file or the protein download manifest."""
    ids: set[str] = set()
    if not path.exists():
        logger.warning(f"Subset file not found: {path}")
        return ids
    with open(path) as f:
        header = f.readline().rstrip().split("\t")
        if id_column not in header:
            logger.error(f"Column '{id_column}' not found in {path.name}")
            return ids
        col_idx = header.index(id_column)
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) <= col_idx:
                continue
            raw = parts[col_idx].strip()
            if not raw:
                continue
            if id_format == "geneid":
                # Strip LOC prefix if present (e.g. LOC101742582 -> 101742582)
                gene_id = raw[3:] if raw.startswith("LOC") else raw
                if gene_id.isdigit():
                    ids.add(gene_id)
            elif id_format == "fbgn" and raw.startswith("FBgn"):
                ids.add(raw)
    return ids


# ---------------------------------------------------------------------------
# Join & output
# ---------------------------------------------------------------------------


def build_ortholog_pairs(
    og_pairs: dict[str, tuple[list[str], list[str]]],
    bm_map: dict[str, dict],
    dm_map: dict[str, dict],
    fbgn_map: dict[str, str],
) -> list[dict]:
    """Join OGs with gene mappings to produce flat ortholog-pair rows."""
    rows: list[dict] = []
    for og_id, (bm_odb_ids, dm_odb_ids) in og_pairs.items():
        level_taxid = og_id.rsplit("at", 1)[-1] if "at" in og_id else ""
        n_bm, n_dm = len(bm_odb_ids), len(dm_odb_ids)
        if n_bm == 1 and n_dm == 1:
            pair_type = "1:1"
        elif n_bm == 1:
            pair_type = "1:many"
        elif n_dm == 1:
            pair_type = "many:1"
        else:
            pair_type = "many:many"

        for bm_odb_id in bm_odb_ids:
            bm_entry = bm_map.get(bm_odb_id)
            if not bm_entry or not bm_entry["gene_id"]:
                continue
            for dm_odb_id in dm_odb_ids:
                dm_entry = dm_map.get(dm_odb_id)
                if not dm_entry:
                    continue
                fbgn = dm_entry["fbgn"]
                symbol = dm_entry["symbol"]
                # Fallback: CG / symbol -> FBgn via FlyBase table
                if not fbgn and symbol and symbol in fbgn_map:
                    fbgn = fbgn_map[symbol]
                if not fbgn:
                    continue
                rows.append(
                    {
                        "og_id": og_id,
                        "og_level_taxid": level_taxid,
                        "pair_type": pair_type,
                        "bm_gene_id": bm_entry["gene_id"],
                        "bm_odb_gene_id": bm_odb_id,
                        "dm_fbgn": fbgn,
                        "dm_gene_symbol": symbol,
                        "dm_odb_gene_id": dm_odb_id,
                    }
                )
    return rows


def deduplicate_ortholog_pairs(
    rows: list[dict], subset_names: list[str]
) -> tuple[list[dict], list[dict], list[dict]]:
    """Collapse rows to unique (bm_gene_id, dm_fbgn) pairs and derive gene lists.

    Returns ``(unique_pairs, bm_genes, dm_fbgns)`` where ``unique_pairs`` is
    one row per unique Bombyx-Drosophila gene pair, and the two gene tables
    aggregate all ortholog partners for each unique gene.
    """
    pair_index: dict[tuple[str, str], dict] = defaultdict(
        lambda: {
            "og_ids": set(),
            "og_level_taxids": set(),
            "pair_types": set(),
            "dm_symbols": set(),
            "subsets": set(),
        }
    )
    for row in rows:
        key = (row["bm_gene_id"], row["dm_fbgn"])
        pair = pair_index[key]
        pair["og_ids"].add(row["og_id"])
        pair["og_level_taxids"].add(row["og_level_taxid"])
        pair["pair_types"].add(row["pair_type"])
        if row.get("dm_gene_symbol"):
            pair["dm_symbols"].add(row["dm_gene_symbol"])
        for s in subset_names:
            if row.get(f"in_{s}") == "yes":
                pair["subsets"].add(s)

    unique_pairs: list[dict] = []
    for (bm, dm), v in pair_index.items():
        unique_pairs.append(
            {
                "bm_gene_id": bm,
                "dm_fbgn": dm,
                "dm_gene_symbol": ";".join(sorted(v["dm_symbols"])),
                "og_ids": ";".join(sorted(v["og_ids"])),
                "og_level_taxids": ";".join(
                    sorted(t for t in v["og_level_taxids"] if t)
                ),
                "pair_types": ";".join(sorted(v["pair_types"])),
                "n_ogs": len(v["og_ids"]),
                **{
                    f"in_{s}": "yes" if s in v["subsets"] else ""
                    for s in subset_names
                },
            }
        )

    bm_index: dict[str, dict] = defaultdict(
        lambda: {"dm_fbgns": set(), "og_ids": set()}
    )
    dm_index: dict[str, dict] = defaultdict(
        lambda: {"bm_gene_ids": set(), "og_ids": set(), "symbols": set()}
    )
    for p in unique_pairs:
        bm = p["bm_gene_id"]
        dm = p["dm_fbgn"]
        bm_index[bm]["dm_fbgns"].add(dm)
        bm_index[bm]["og_ids"].update(p["og_ids"].split(";"))
        dm_index[dm]["bm_gene_ids"].add(bm)
        dm_index[dm]["og_ids"].update(p["og_ids"].split(";"))
        if p["dm_gene_symbol"]:
            dm_index[dm]["symbols"].update(p["dm_gene_symbol"].split(";"))

    bm_genes = [
        {
            "bm_gene_id": g,
            "n_ortholog_fbgns": len(v["dm_fbgns"]),
            "ortholog_dm_fbgns": ";".join(sorted(v["dm_fbgns"])),
            "n_ogs": len(v["og_ids"]),
            "og_ids": ";".join(sorted(v["og_ids"])),
        }
        for g, v in bm_index.items()
    ]

    dm_fbgns = [
        {
            "dm_fbgn": g,
            "dm_gene_symbol": ";".join(sorted(v["symbols"])),
            "n_ortholog_bm_genes": len(v["bm_gene_ids"]),
            "ortholog_bm_gene_ids": ";".join(sorted(v["bm_gene_ids"])),
            "n_ogs": len(v["og_ids"]),
            "og_ids": ";".join(sorted(v["og_ids"])),
        }
        for g, v in dm_index.items()
    ]

    return unique_pairs, bm_genes, dm_fbgns


def annotate_subsets(
    rows: list[dict], subsets: dict[str, set[str]]
) -> list[dict]:
    """Add ``in_<subset>`` flags to each row."""
    # Pre-build name -> species mapping for efficiency
    name_to_species = {name: species for name, _p, _c, species, _f in _SUBSETS}
    for row in rows:
        for name, ids in subsets.items():
            species = name_to_species[name]
            if species == "bm":
                row[f"in_{name}"] = "yes" if row["bm_gene_id"] in ids else ""
            else:
                row[f"in_{name}"] = "yes" if row["dm_fbgn"] in ids else ""
    return rows


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    """Write rows to a tab-separated file."""
    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write(
                "\t".join(str(row.get(col, "")) for col in columns) + "\n"
            )
    logger.info(f"  Wrote {len(rows):,} rows -> {path}")


# ---------------------------------------------------------------------------
# Subset pairwise intersections
# ---------------------------------------------------------------------------

_INTERSECTION_COLUMNS = [
    "subset_a",
    "subset_b",
    "intersection_type",
    "bm_gene_id",
    "dm_fbgn",
    "dm_gene_symbol",
    "og_ids",
    "pair_types",
]

_INTERSECTION_COUNT_COLUMNS = [
    "subset_a",
    "subset_b",
    "intersection_type",
    "intersection_rows",
    "unique_bm_genes",
    "unique_dm_fbgns",
    "unique_ogs",
]


def _intersection_type(name_a: str, name_b: str) -> str:
    """Classify a subset pair as cross-species or same-species."""
    species_a = next(
        s for n, _p, _c, s, _f in _SUBSETS if n == name_a
    )
    species_b = next(
        s for n, _p, _c, s, _f in _SUBSETS if n == name_b
    )
    if species_a != species_b:
        return "cross_species"
    return f"same_species_{species_a}"


def compute_subset_intersections(
    unique_pairs: list[dict], subset_names: list[str]
) -> list[dict]:
    """Compute Bombyx vs Drosophila cross-species subset intersections.

    Each output row is a unique ortholog pair (``bm_gene_id``, ``dm_fbgn``)
    that belongs to one Bombyx subset and one Drosophila subset simultaneously.
    Same-species (bm-bm or dm-dm) intersections are skipped.
    """
    # Pre-compute which unique pairs belong to each subset to avoid
    # repeated scans.
    subset_pair_idxs: dict[str, list[int]] = {}
    for name in subset_names:
        flag = f"in_{name}"
        subset_pair_idxs[name] = [
            idx for idx, p in enumerate(unique_pairs) if p.get(flag) == "yes"
        ]

    intersection_rows: list[dict] = []
    species_map = {n: s for n, _p, _c, s, _f in _SUBSETS}
    for i, name_a in enumerate(subset_names):
        species_a = species_map[name_a]
        idxs_a = set(subset_pair_idxs[name_a])
        for name_b in subset_names[i + 1 :]:
            species_b = species_map[name_b]
            if species_a == species_b:
                continue
            for idx in subset_pair_idxs[name_b]:
                if idx in idxs_a:
                    pair = unique_pairs[idx]
                    intersection_rows.append(
                        {
                            "subset_a": name_a,
                            "subset_b": name_b,
                            "intersection_type": "cross_species",
                            "bm_gene_id": pair["bm_gene_id"],
                            "dm_fbgn": pair["dm_fbgn"],
                            "dm_gene_symbol": pair["dm_gene_symbol"],
                            "og_ids": pair["og_ids"],
                            "pair_types": pair["pair_types"],
                        }
                    )
    return intersection_rows


def summarize_intersections(
    intersection_rows: list[dict],
) -> list[dict]:
    """Aggregate per-pair counts for the summary table."""
    summary_rows: list[dict] = []
    # Group rows by (subset_a, subset_b, intersection_type)
    groups: dict[tuple[str, str, str], list[dict]] = {}
    for row in intersection_rows:
        key = (row["subset_a"], row["subset_b"], row["intersection_type"])
        groups.setdefault(key, []).append(row)

    for (name_a, name_b, itype), group in groups.items():
        summary_rows.append(
            {
                "subset_a": name_a,
                "subset_b": name_b,
                "intersection_type": itype,
                "intersection_rows": len(group),
                "unique_bm_genes": len({r["bm_gene_id"] for r in group}),
                "unique_dm_fbgns": len({r["dm_fbgn"] for r in group}),
                "unique_ogs": len({
                    og for r in group for og in r["og_ids"].split(";") if og
                }),
            }
        )
    return summary_rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@app.command()
def main(
    force: bool = typer.Option(False, "--force", help="Overwrite existing outputs."),
    og2genes: Path = typer.Option(OG2GENES_PATH, "--og2genes"),
    genes_tab: Path = typer.Option(GENES_TAB_PATH, "--genes-tab"),
    fbgn_map: Path = typer.Option(FBGN_MAP_PATH, "--fbgn-map"),
    output_dir: Path = typer.Option(OUTPUT_DIR, "--output-dir"),
) -> None:
    """Extract Bombyx mori x Drosophila melanogaster ortholog pairs from OrthoDB v12."""
    t_start = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check inputs
    for p in [og2genes, genes_tab, fbgn_map]:
        if not p.exists():
            logger.error(f"Input file not found: {p}")
            raise typer.Exit(1)

    # Check existing outputs
    all_out = output_dir / "bm_dm_orthologs_all.tsv"
    summary_out = output_dir / "ortholog_extraction_summary.json"
    if all_out.exists() and not force:
        logger.error(f"Output exists: {all_out}. Use --force to overwrite.")
        raise typer.Exit(1)

    # Step 1: scan OG2genes for OGs with both species
    logger.info("Step 1: Scan OG2genes.tab for ortholog groups with both species")
    og_pairs = scan_og2genes(og2genes)

    # Step 2: scan genes.tab for 7091/7227 gene mappings
    logger.info("Step 2: Scan genes.tab for gene ID mappings")
    bm_map, dm_map = scan_genes_tab(genes_tab)

    # Step 3: load FlyBase CG -> FBgn fallback
    logger.info("Step 3: Load FlyBase annotation-ID mapping")
    fbgn_fallback = load_fbgn_map(fbgn_map)
    logger.info(f"  FlyBase mapping: {len(fbgn_fallback):,} entries")

    # Step 4: build ortholog pairs
    logger.info("Step 4: Build ortholog pairs")
    rows = build_ortholog_pairs(og_pairs, bm_map, dm_map, fbgn_fallback)
    logger.info(f"  Total resolved pairs: {len(rows):,}")

    # Step 5: load gene subsets and annotate
    logger.info("Step 5: Load gene subsets and annotate")
    subsets = load_gene_subsets()
    rows = annotate_subsets(rows, subsets)

    # Step 6: deduplicate and write unique outputs
    logger.info("Step 6: Deduplicate ortholog pairs and genes")
    subset_columns = [f"in_{name}" for name, *_ in _SUBSETS]
    subset_names = [name for name, *_ in _SUBSETS]
    unique_pairs, bm_genes, dm_fbgns = deduplicate_ortholog_pairs(
        rows, subset_names
    )
    logger.info(f"  Unique ortholog pairs: {len(unique_pairs):,}")
    logger.info(f"  Unique Bombyx genes with orthologs: {len(bm_genes):,}")
    logger.info(f"  Unique Drosophila FBgn with orthologs: {len(dm_fbgns):,}")

    # Remove old expanded outputs that are superseded by unique files
    old_files = [
        output_dir / "bm_dm_orthologs_all.tsv",
        output_dir / "subset_pairwise_intersections.tsv",
        output_dir / "subset_pairwise_intersection_counts.tsv",
    ]
    for old in old_files:
        if old.exists():
            old.unlink()
            logger.info(f"  Removed old file: {old.name}")

    all_columns = _UNIQUE_PAIR_COLUMNS + subset_columns
    unique_out = output_dir / "bm_dm_ortholog_pairs_unique.tsv"
    write_tsv(unique_out, unique_pairs, all_columns)

    bm_genes_out = output_dir / "bm_ortholog_genes_unique.tsv"
    write_tsv(bm_genes_out, bm_genes, _UNIQUE_BM_COLUMNS)

    dm_fbgns_out = output_dir / "dm_ortholog_fbgns_unique.tsv"
    write_tsv(dm_fbgns_out, dm_fbgns, _UNIQUE_DM_COLUMNS)

    # Step 7: compute cross-species subset intersections on unique pairs
    logger.info("Step 7: Compute cross-species subset intersections")
    intersection_rows = compute_subset_intersections(unique_pairs, subset_names)
    intersection_out = output_dir / "subset_cross_species_intersections.tsv"
    write_tsv(intersection_out, intersection_rows, _INTERSECTION_COLUMNS)

    intersection_counts = summarize_intersections(intersection_rows)
    intersection_counts_out = (
        output_dir / "subset_cross_species_intersection_counts.tsv"
    )
    write_tsv(
        intersection_counts_out, intersection_counts, _INTERSECTION_COUNT_COLUMNS
    )
    logger.info(
        f"  {len(intersection_counts)} cross-species subset pairs, "
        f"{sum(r['intersection_rows'] for r in intersection_counts):,} total rows"
    )

    # Summary JSON
    pair_type_counts: dict[str, int] = {}
    for r in rows:
        pair_type_counts[r["pair_type"]] = (
            pair_type_counts.get(r["pair_type"], 0) + 1
        )

    # Unique-pair pair_type summary for consistency
    unique_pair_type_counts: dict[str, int] = {}
    for p in unique_pairs:
        for pt in p["pair_types"].split(";"):
            if pt:
                unique_pair_type_counts[pt] = (
                    unique_pair_type_counts.get(pt, 0) + 1
                )

    summary = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "input_files": {
            "og2genes": str(og2genes.relative_to(PROJ_ROOT)),
            "genes_tab": str(genes_tab.relative_to(PROJ_ROOT)),
            "fbgn_map": str(fbgn_map.relative_to(PROJ_ROOT)),
        },
        "taxonomy": {
            "bombyx_mori": BM_TAXID,
            "drosophila_melanogaster": DM_TAXID,
        },
        "counts": {
            "ogs_with_both_species": len(og_pairs),
            "bm_genes_mapped": len(bm_map),
            "dm_genes_mapped": len(dm_map),
            "ortholog_pairs_total": len(rows),
            "by_pair_type": pair_type_counts,
            "unique_ortholog_pairs": len(unique_pairs),
            "unique_bm_genes_with_orthologs": len(bm_genes),
            "unique_dm_fbgns_with_orthologs": len(dm_fbgns),
            "by_unique_pair_type": unique_pair_type_counts,
        },
        "subsets": {
            name: {
                "subset_ids": len(subsets[name]),
                "pairs_matched": sum(
                    1 for p in unique_pairs if p.get(f"in_{name}") == "yes"
                ),
            }
            for name, *_ in _SUBSETS
        },
        "cross_species_subset_intersections": {
            f"{r['subset_a']}__x__{r['subset_b']}": {
                "intersection_type": r["intersection_type"],
                "intersection_rows": r["intersection_rows"],
                "unique_bm_genes": r["unique_bm_genes"],
                "unique_dm_fbgns": r["unique_dm_fbgns"],
                "unique_ogs": r["unique_ogs"],
            }
            for r in intersection_counts
        },
        "elapsed_seconds": round(time.time() - t_start, 1),
    }
    with open(summary_out, "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"  Summary -> {summary_out}")
    logger.info(f"Done in {summary['elapsed_seconds']}s")


if __name__ == "__main__":
    app()
