#!/usr/bin/env python3
"""
Calculate gene-level TPM for all samples using featureCounts + TPM conversion.

Workflow:
  1. Run featureCounts on all sample BAMs simultaneously
  2. Convert raw counts to TPM: TPM = (count/length) / sum(count/length) * 1e6

Usage:
  python calc_tpm.py --bam-dir star --gtf genome/genomic_chr.gtf \
      --output-dir results/tpm --threads 4
"""

import json
from pathlib import Path
import subprocess

from loguru import logger
import pandas as pd
import typer

app = typer.Typer()


def run_featurecounts(
    bams: list[Path],
    gtf: Path,
    out_counts: Path,
    threads: int,
    strand: int,
) -> Path:
    """
    Run featureCounts on all BAMs at once.
    strand: 0=unstranded, 1=stranded, 2=reversely stranded
    """
    if out_counts.exists():
        logger.info(f"Counts file already exists at {out_counts}, skipping")
        return out_counts

    out_counts.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "featureCounts",
        "-T", str(threads),
        "-p",                   # paired-end
        "-B",                   # both mates mapped
        "-s", str(strand),
        "-g", "gene_id",
        "-a", str(gtf),
        "-o", str(out_counts),
        *[str(b) for b in bams],
    ]
    logger.info(f"Running featureCounts on {len(bams)} BAMs ...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"featureCounts failed:\n{result.stderr}")
        raise RuntimeError("featureCounts failed")
    logger.info(f"Counts written to {out_counts}")
    return out_counts


def counts_to_tpm(counts_file: Path, sample_names: list[str]) -> pd.DataFrame:
    """
    Parse featureCounts output and convert to TPM.

    featureCounts output columns:
      Geneid, Chr, Start, End, Strand, Length, <bam1>, <bam2>, ...
    """
    df = pd.read_csv(counts_file, sep="\t", comment="#")
    df = df.set_index("Geneid")

    lengths = df["Length"]
    # BAM columns start after the annotation columns (Chr, Start, End, Strand, Length)
    count_cols = df.columns[5:]
    counts = df[count_cols].copy()

    # Rename columns to sample names (strip path)
    counts.columns = sample_names

    # TPM = (count / length) / sum(count / length) * 1e6
    rpk = counts.div(lengths, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    tpm.index.name = "gene_id"
    return tpm


@app.command()
def main(
    bam_dir:    Path = typer.Option(...,               "--bam-dir",    help="Directory containing per-sample STAR BAM subdirectories"),
    gtf:        Path = typer.Option(...,               "--gtf",        help="Genome annotation GTF file"),
    output_dir: Path = typer.Option("results/tpm",     "--output-dir", help="Output directory for TPM results"),
    config:     Path = typer.Option("config.json",     "--config",     help="Workflow config.json (for sample list)"),
    threads:    int  = typer.Option(4,                 "--threads",    help="Number of threads for featureCounts"),
    strand:     int  = typer.Option(2,                 "--strand",     help="Strandness: 0=unstranded, 1=stranded, 2=reversely stranded"),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Discover samples ──────────────────────────────────────────────────────
    if config.exists():
        with open(config) as f:
            cfg = json.load(f)
        samples = [s for grp in cfg.get("samples", {}).values() for s in grp]
        logger.info(f"Loaded {len(samples)} samples from {config}: {samples}")
    else:
        samples = [d.name for d in sorted(bam_dir.iterdir()) if d.is_dir()]
        logger.info(f"Auto-discovered {len(samples)} samples from {bam_dir}: {samples}")

    bams = []
    for sample in samples:
        bam = bam_dir / sample / "Aligned.sortedByCoord.out.bam"
        if not bam.exists():
            logger.warning(f"BAM not found for {sample}: {bam}, skipping")
        else:
            bams.append(bam)

    if not bams:
        logger.error("No BAM files found, exiting")
        raise typer.Exit(1)

    # Keep only samples that have a BAM
    samples = [b.parent.name for b in bams]

    # ── Run featureCounts ─────────────────────────────────────────────────────
    out_counts = output_dir / "counts.txt"
    run_featurecounts(bams, gtf, out_counts, threads, strand)

    # ── Convert to TPM ────────────────────────────────────────────────────────
    tpm = counts_to_tpm(out_counts, samples)

    out_tsv = output_dir / "tpm_matrix.tsv"
    tpm.to_csv(out_tsv, sep="\t", float_format="%.4f")
    logger.success(
        f"TPM matrix saved to {out_tsv}  "
        f"({tpm.shape[0]} genes × {tpm.shape[1]} samples)"
    )


if __name__ == "__main__":
    app()
