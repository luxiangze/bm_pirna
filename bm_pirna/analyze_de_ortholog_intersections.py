"""Cross-species ortholog intersections for differentially expressed genes.

Reads the unique Bombyx x Drosophila ortholog pairs extracted from OrthoDB v12
and computes cross-species intersections between differentially expressed
gene sets stored under ``data/external/homology/differentially_expressed_genes``.

Outputs (under ``data/processed/homology/orthodb_orthologs/``):
  - ``de_cross_species_intersections.tsv``
  - ``de_cross_species_intersection_counts.tsv``
  - ``de_intersection_summary.json``
"""

from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path

from loguru import logger
import pandas as pd
import typer

PROJ_ROOT = Path(__file__).resolve().parents[1]
ORTHODB_DIR = PROJ_ROOT / "data/processed/homology/orthodb_orthologs"
ORTHODB_PAIRS = ORTHODB_DIR / "bm_dm_ortholog_pairs_unique.tsv"
GENE_INFO = PROJ_ROOT / "data/external/gene_info_table.tsv"
DE_DIR = PROJ_ROOT / "data/external/homology/differentially_expressed_genes"
OUTPUT_DIR = ORTHODB_DIR

app = typer.Typer(add_completion=False, pretty_exceptions_enable=False)

# Ordered list of contrasts to analyze. Each tuple is:
#   (name, species, file_path, id_column, id_style)
# species: "bm" or "dm"; id_style: "loc" (LOC prefix), "symbol", or "fbgn"
_CONTRASTS: list[tuple[str, str, Path, str, str]] = [
    (
        "bm_DE_1st_control_vs_sugp1",
        "bm",
        DE_DIR
        / "RNAdiff_workflow_ncbi_20260311/differential/condition_control_treated.deseq2.results_filtered.tsv",
        "gene_id",
        "loc",
    ),
    (
        "bm_DE_2nd_control_vs_sugp1",
        "bm",
        DE_DIR
        / "bm_differentialabundance_20260821/differential/contrasts/condition_control_treated_bm_differentialabundance_202607025.deseq2.results_filtered.tsv",
        "gene_id",
        "symbol",
    ),
    (
        "bm_DE_2nd_p200_vs_sugp1",
        "bm",
        DE_DIR
        / "bm_differentialabundance_20260821/differential/contrasts/condition_p200_treated_bm_differentialabundance_202607025.deseq2.results_filtered.tsv",
        "gene_id",
        "symbol",
    ),
    (
        "dm_DE_control_vs_sugp1",
        "dm",
        DE_DIR
        / "dm_differentialabundance_20260419/differential/condition_control_treated.deseq2.results_filtered.tsv",
        "gene_id",
        "fbgn",
    ),
]


def load_gene_symbol_map(path: Path) -> dict[str, str]:
    """Build symbol -> NCBI GeneID map for Bombyx mori."""
    df = pd.read_csv(path, sep="\t", usecols=["gene_name", "gene_id"], dtype=str)
    # Keep the first mapping when a symbol maps to multiple gene IDs
    return dict(zip(df["gene_name"], df["gene_id"]))


def load_de_set(path: Path, id_column: str, id_style: str, symbol_map: dict[str, str]) -> set[str]:
    """Load a DE result file and return normalized gene IDs."""
    if not path.exists():
        logger.warning(f"DE file not found: {path}")
        return set()
    df = pd.read_csv(path, sep="\t", dtype=str)
    if id_column not in df.columns:
        logger.warning(f"Column '{id_column}' not found in {path.name}")
        return set()

    raw_ids = df[id_column].dropna().astype(str).str.strip()
    normalized: set[str] = set()
    for raw in raw_ids:
        if not raw or raw == "NA":
            continue
        if id_style == "loc":
            if raw.startswith("LOC") and raw[3:].isdigit():
                normalized.add(raw[3:])
            elif raw.isdigit():
                normalized.add(raw)
        elif id_style == "symbol":
            if raw.isdigit():
                normalized.add(raw)
            elif raw.startswith("LOC") and raw[3:].isdigit():
                normalized.add(raw[3:])
            elif raw in symbol_map:
                normalized.add(symbol_map[raw])
        elif id_style == "fbgn":
            if raw.startswith("FBgn"):
                normalized.add(raw)
    logger.info(f"  Loaded {len(normalized)} IDs from {path.name}")
    return normalized


def compute_cross_species_intersections(
    pairs_df: pd.DataFrame, de_sets: dict[str, set[str]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute all bm vs dm cross-species intersections from unique ortholog pairs."""
    bm_sets = {name: ids for name, ids in de_sets.items() if name.startswith("bm_")}
    dm_sets = {name: ids for name, ids in de_sets.items() if name.startswith("dm_")}

    rows: list[dict] = []
    counts: list[dict] = []

    for bm_name, bm_ids in bm_sets.items():
        bm_pairs = pairs_df[pairs_df["bm_gene_id"].isin(bm_ids)]
        for dm_name, dm_ids in dm_sets.items():
            subset = bm_pairs[bm_pairs["dm_fbgn"].isin(dm_ids)]
            if subset.empty:
                counts.append(
                    {
                        "subset_a": bm_name,
                        "subset_b": dm_name,
                        "intersection_type": "cross_species",
                        "intersection_rows": 0,
                        "unique_bm_genes": 0,
                        "unique_dm_fbgns": 0,
                        "unique_ogs": 0,
                    }
                )
                continue
            for _, row in subset.iterrows():
                rows.append(
                    {
                        "subset_a": bm_name,
                        "subset_b": dm_name,
                        "bm_gene_id": row["bm_gene_id"],
                        "dm_fbgn": row["dm_fbgn"],
                        "dm_gene_symbol": row["dm_gene_symbol"],
                        "og_ids": row["og_ids"],
                        "pair_types": row["pair_types"],
                    }
                )
            counts.append(
                {
                    "subset_a": bm_name,
                    "subset_b": dm_name,
                    "intersection_type": "cross_species",
                    "intersection_rows": len(subset),
                    "unique_bm_genes": subset["bm_gene_id"].nunique(),
                    "unique_dm_fbgns": subset["dm_fbgn"].nunique(),
                    "unique_ogs": len({
                        og for ogs in subset["og_ids"].dropna()
                        for og in str(ogs).split(";") if og
                    }),
                }
            )

    return pd.DataFrame(rows), pd.DataFrame(counts)


def compute_ip_ms_intersections(
    pairs_df: pd.DataFrame, de_sets: dict[str, set[str]], ipms_ids: set[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute dm_DE vs IP-MS cross-species intersections."""
    dm_sets = {name: ids for name, ids in de_sets.items() if name.startswith("dm_")}
    rows: list[dict] = []
    counts: list[dict] = []

    for dm_name, dm_ids in dm_sets.items():
        dm_pairs = pairs_df[pairs_df["dm_fbgn"].isin(dm_ids)]
        subset = dm_pairs[dm_pairs["bm_gene_id"].isin(ipms_ids)]
        if subset.empty:
            counts.append(
                {
                    "subset_a": dm_name,
                    "subset_b": "IP-MS_sugp1",
                    "intersection_type": "cross_species",
                    "intersection_rows": 0,
                    "unique_bm_genes": 0,
                    "unique_dm_fbgns": 0,
                    "unique_ogs": 0,
                }
            )
            continue
        for _, row in subset.iterrows():
            rows.append(
                {
                    "subset_a": dm_name,
                    "subset_b": "IP-MS_sugp1",
                    "bm_gene_id": row["bm_gene_id"],
                    "dm_fbgn": row["dm_fbgn"],
                    "dm_gene_symbol": row["dm_gene_symbol"],
                    "og_ids": row["og_ids"],
                    "pair_types": row["pair_types"],
                }
            )
        counts.append(
            {
                "subset_a": dm_name,
                "subset_b": "IP-MS_sugp1",
                "intersection_type": "cross_species",
                "intersection_rows": len(subset),
                "unique_bm_genes": subset["bm_gene_id"].nunique(),
                "unique_dm_fbgns": subset["dm_fbgn"].nunique(),
                "unique_ogs": len({
                    og for ogs in subset["og_ids"].dropna()
                    for og in str(ogs).split(";") if og
                }),
            }
        )

    return pd.DataFrame(rows), pd.DataFrame(counts)


@app.command()
def main(
    ortholog_pairs: Path = typer.Option(ORTHODB_PAIRS, "--ortholog-pairs"),
    gene_info: Path = typer.Option(GENE_INFO, "--gene-info"),
    ipms_file: Path = typer.Option(
        PROJ_ROOT
        / "data/external/homology/differentially_alternatively_spliced_genes/IP-MS_sugp1_gene_ids.txt",
        "--ipms-file",
    ),
    output_dir: Path = typer.Option(OUTPUT_DIR, "--output-dir"),
) -> None:
    """Analyze cross-species ortholog intersections for DE gene sets."""
    t_start = datetime.now()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not ortholog_pairs.exists():
        logger.error(f"Ortholog pairs file not found: {ortholog_pairs}")
        raise typer.Exit(1)

    logger.info("Loading unique ortholog pairs")
    pairs_df = pd.read_csv(ortholog_pairs, sep="\t", dtype=str)

    logger.info("Loading Bombyx gene symbol map")
    symbol_map = load_gene_symbol_map(gene_info)

    logger.info("Loading DE gene sets")
    de_sets: dict[str, set[str]] = {}
    for name, species, path, col, style in _CONTRASTS:
        ids = load_de_set(path, col, style, symbol_map)
        de_sets[name] = ids
        logger.info(f"  {name}: {len(ids)} {species} IDs")

    # IP-MS set uses the same file as the AS analysis; column is GeneID
    logger.info("Loading IP-MS gene set")
    ipms_df = pd.read_csv(ipms_file, sep="\t", dtype=str)
    ipms_ids: set[str] = set(
        ipms_df["GeneID"].dropna().astype(str).str.strip()
    ) & {gid for gid in ipms_df["GeneID"].dropna().astype(str).str.strip() if gid.isdigit()}
    logger.info(f"  IP-MS_sugp1: {len(ipms_ids)} Bombyx IDs")

    logger.info("Computing bm_DE vs dm_DE cross-species intersections")
    de_rows, de_counts = compute_cross_species_intersections(pairs_df, de_sets)

    logger.info("Computing dm_DE vs IP-MS cross-species intersections")
    ipms_rows, ipms_counts = compute_ip_ms_intersections(pairs_df, de_sets, ipms_ids)

    all_rows = pd.concat([de_rows, ipms_rows], ignore_index=True)
    all_counts = pd.concat([de_counts, ipms_counts], ignore_index=True)

    rows_out = output_dir / "de_cross_species_intersections.tsv"
    counts_out = output_dir / "de_cross_species_intersection_counts.tsv"
    all_rows.to_csv(rows_out, sep="\t", index=False)
    all_counts.to_csv(counts_out, sep="\t", index=False)
    logger.info(f"  Wrote {len(all_rows)} intersection rows -> {rows_out}")
    logger.info(f"  Wrote {len(all_counts)} summary rows -> {counts_out}")

    summary = {
        "timestamp": t_start.isoformat(timespec="seconds"),
        "input_files": {
            "ortholog_pairs": str(ortholog_pairs.relative_to(PROJ_ROOT)),
            "gene_info": str(gene_info.relative_to(PROJ_ROOT)),
        },
        "de_sets": {
            name: len(ids) for name, ids in de_sets.items()
        },
        "ipms_ids": len(ipms_ids),
        "cross_species_intersections": {
            f"{r['subset_a']}__x__{r['subset_b']}": {
                "intersection_rows": int(r["intersection_rows"]),
                "unique_bm_genes": int(r["unique_bm_genes"]),
                "unique_dm_fbgns": int(r["unique_dm_fbgns"]),
                "unique_ogs": int(r["unique_ogs"]),
            }
            for _, r in all_counts.iterrows()
        },
    }
    summary_out = output_dir / "de_intersection_summary.json"
    with open(summary_out, "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"  Summary -> {summary_out}")


if __name__ == "__main__":
    app()
