#!/usr/bin/env python3
"""
Analyze piRNA cluster expression changes from featureCounts output.

Workflow:
  1. Parse featureCounts count matrix (with cluster length)
  2. RPK normalization (reads / cluster_length_kb) per sample
  3. Compute mean log2FC (group2 / group1, pseudocount=0.5)
  4. Compute per-REP log2FC and assess REP1 vs REP2 reproducibility
     via Pearson and Spearman correlation
  5. Output results TSV + log2FC scatter + REP correlation plot PDF

Usage:
  python analyze_clusters.py --counts results/cluster_counts.txt \
      --config config.json --group1 control --group2 sugp1 \
      --output-dir results/cluster_analysis
"""

import json
from pathlib import Path
import re

from loguru import logger
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
import typer

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
mpl.rcParams['axes.unicode_minus'] = False

app = typer.Typer()


def load_config(config_path: Path) -> dict:
    with open(config_path) as f:
        return json.load(f)


def parse_counts(counts_path: Path) -> tuple[pd.DataFrame, pd.Series]:
    """Parse featureCounts output; return (count_df, length_series)."""
    df = pd.read_csv(counts_path, sep='\t', comment='#')
    df = df.set_index('Geneid')
    lengths = df['Length'].astype(float)
    count_cols = [c for c in df.columns if c not in ('Chr', 'Start', 'End', 'Strand', 'Length')]
    return df[count_cols], lengths


def bam_col_to_sample(col: str) -> str:
    """Extract sample name from BAM path, e.g. star/CONTROL_REP1/Aligned... -> CONTROL_REP1."""
    m = re.search(r'star/([^/]+)/', col)
    return m.group(1) if m else col


def rpk_normalize(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    """RPK = counts / (cluster_length / 1000). Length-normalizes within each sample."""
    length_kb = lengths / 1000.0
    return counts.divide(length_kb, axis=0)


def log2fc_per_sample(
    rpk: pd.DataFrame,
    g1_sample: str,
    g2_sample: str,
    pseudocount: float = 0.5,
) -> pd.Series:
    """Compute log2FC for a single g2/g1 sample pair."""
    return np.log2((rpk[g2_sample] + pseudocount) / (rpk[g1_sample] + pseudocount))


def pearson_r(x: pd.Series, y: pd.Series) -> float:
    mask = x.notna() & y.notna()
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def spearman_r(x: pd.Series, y: pd.Series) -> float:
    mask = x.notna() & y.notna()
    n = mask.sum()
    if n < 3:
        return float('nan')
    rank_x = x[mask].rank()
    rank_y = y[mask].rank()
    return pearson_r(rank_x, rank_y)


def plot_rep_correlation(
    rep1_log2fc: pd.Series,
    rep2_log2fc: pd.Series,
    group1: str,
    group2: str,
    output_path: Path,
    log2fc_threshold: float,
    pr: float,
    sr: float,
) -> None:
    """Scatter plot of REP1 log2FC vs REP2 log2FC with correlation stats."""
    x, y = rep1_log2fc.dropna(), rep2_log2fc.dropna()
    common = x.index.intersection(y.index)
    x, y = x[common], y[common]

    up   = (x >= log2fc_threshold) & (y >= log2fc_threshold)
    down = (x <= -log2fc_threshold) & (y <= -log2fc_threshold)
    ns   = ~(up | down)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(x[ns],   y[ns],   c='#aaaaaa', s=18, alpha=0.6, label=f'NS ({ns.sum()})')
    ax.scatter(x[down], y[down], c='#3a86ff', s=28, alpha=0.8, label=f'Consistent down ({down.sum()})')
    ax.scatter(x[up],   y[up],   c='#ff4757', s=28, alpha=0.8, label=f'Consistent up ({up.sum()})')

    lim = max(abs(x).max(), abs(y).max()) * 1.1
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axhline(0, color='black', lw=0.6, alpha=0.4)
    ax.axvline(0, color='black', lw=0.6, alpha=0.4)
    ax.plot([-lim, lim], [-lim, lim], color='black', lw=0.8, ls='--', alpha=0.4)

    ax.set_xlabel(f'REP1  log\u2082FC  ({group2}/{group1})', fontsize=13)
    ax.set_ylabel(f'REP2  log\u2082FC  ({group2}/{group1})', fontsize=13)
    ax.tick_params(labelsize=11)
    ax.legend(fontsize=10)

    stats_text = f'Pearson r = {pr:.3f}\nSpearman \u03c1 = {sr:.3f}\nn = {len(common)}'
    ax.text(0.03, 0.97, stats_text, transform=ax.transAxes, fontsize=11,
            va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_path, format='pdf')
    plt.close()
    logger.info(f"REP correlation plot saved to {output_path}")


def plot_log2fc_distribution(
    result_df: pd.DataFrame,
    group1: str,
    group2: str,
    output_path: Path,
    log2fc_threshold: float,
) -> None:
    """Histogram of mean log2FC values."""
    fc = result_df['mean_log2FC'].dropna()
    up   = (fc >= log2fc_threshold).sum()
    down = (fc <= -log2fc_threshold).sum()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(fc, bins=40, color='#4a90e2', edgecolor='white', alpha=0.8)
    ax.axvline( log2fc_threshold, color='#ff4757', lw=1.2, ls='--', label=f'|log\u2082FC| = {log2fc_threshold}')
    ax.axvline(-log2fc_threshold, color='#ff4757', lw=1.2, ls='--')
    ax.set_xlabel(f'Mean log\u2082FC  ({group2} / {group1})', fontsize=13)
    ax.set_ylabel('Number of clusters', fontsize=13)
    ax.tick_params(labelsize=11)
    ax.legend(fontsize=10)

    note = f'Up: {up}   Down: {down}   Total: {len(fc)}'
    ax.text(0.97, 0.97, note, transform=ax.transAxes, fontsize=10,
            ha='right', va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_path, format='pdf')
    plt.close()
    logger.info(f"log2FC distribution plot saved to {output_path}")


@app.command()
def main(
    counts_path:      Path  = typer.Option(...,   '--counts',    help='featureCounts output file'),
    config_path:      Path  = typer.Option('config.json', '--config', help='Workflow config.json'),
    group1:           str   = typer.Option(...,   '--group1',    help='Reference group name (e.g. control)'),
    group2:           str   = typer.Option(...,   '--group2',    help='Treatment group name (e.g. sugp1)'),
    output_dir:       Path  = typer.Option('results/cluster_analysis', '--output-dir'),
    min_rpk:          float = typer.Option(1.0,   '--min-rpk',   help='Min mean RPK threshold for filtering'),
    log2fc_threshold: float = typer.Option(1.0,   '--log2fc',    help='|log2FC| threshold for plots'),
    pseudocount:      float = typer.Option(0.5,   '--pseudocount', help='Pseudocount added before log2 transform'),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load config ───────────────────────────────────────────────────────────
    logger.info(f"Loading config from {config_path}")
    cfg = load_config(config_path)
    sample_groups: dict[str, list[str]] = cfg['samples']

    for grp in (group1, group2):
        if grp not in sample_groups:
            logger.error(f"Group '{grp}' not found in config. Available: {list(sample_groups)}")
            raise typer.Exit(1)

    g1_samples = sample_groups[group1]
    g2_samples = sample_groups[group2]
    logger.info(f"Group1 ({group1}): {g1_samples}")
    logger.info(f"Group2 ({group2}): {g2_samples}")

    if len(g1_samples) != len(g2_samples):
        logger.warning("Groups have unequal numbers of replicates; REP correlation will use min(n1,n2) pairs")

    # ── Parse counts + lengths ────────────────────────────────────────────────
    logger.info(f"Parsing count matrix from {counts_path}")
    raw, lengths = parse_counts(counts_path)
    raw.columns = [bam_col_to_sample(c) for c in raw.columns]

    missing = [s for s in g1_samples + g2_samples if s not in raw.columns]
    if missing:
        logger.error(f"Samples not found in count matrix: {missing}. Available: {list(raw.columns)}")
        raise typer.Exit(1)

    # ── RPK normalization ─────────────────────────────────────────────────────
    logger.info("Normalizing to RPK (reads / cluster_length_kb)")
    rpk = rpk_normalize(raw, lengths)

    # ── Low-expression filter ─────────────────────────────────────────────────
    all_samples = g1_samples + g2_samples
    mean_rpk = rpk[all_samples].mean(axis=1)
    keep = mean_rpk >= min_rpk
    rpk_filt = rpk.loc[keep]
    logger.info(f"Kept {keep.sum()}/{len(keep)} clusters (mean RPK >= {min_rpk})")

    # ── Mean log2FC ───────────────────────────────────────────────────────────
    logger.info("Computing mean log2FC")
    result_rows = []
    for gene in tqdm(rpk_filt.index, desc="log2FC"):
        g1_vals = rpk_filt.loc[gene, g1_samples].values.astype(float)
        g2_vals = rpk_filt.loc[gene, g2_samples].values.astype(float)
        mean1 = g1_vals.mean()
        mean2 = g2_vals.mean()
        mean_fc = float(np.log2((mean2 + pseudocount) / (mean1 + pseudocount)))
        result_rows.append({
            'cluster':        gene,
            f'mean_rpk_{group1}': mean1,
            f'mean_rpk_{group2}': mean2,
            'mean_log2FC':    mean_fc,
        })

    result_df = pd.DataFrame(result_rows).set_index('cluster')

    # ── Per-REP log2FC for reproducibility ───────────────────────────────────
    n_pairs = min(len(g1_samples), len(g2_samples))
    rep_fc_cols: list[str] = []
    for i in range(n_pairs):
        col = f'log2FC_REP{i+1}'
        result_df[col] = log2fc_per_sample(
            rpk_filt, g1_samples[i], g2_samples[i], pseudocount
        )
        rep_fc_cols.append(col)

    # ── Save full results ─────────────────────────────────────────────────────
    out_tsv = output_dir / f"{group2}_vs_{group1}_clusters.tsv"
    result_df.to_csv(out_tsv, sep='\t', float_format='%.6f')
    logger.info(f"Results saved to {out_tsv}")

    # ── REP correlation ───────────────────────────────────────────────────────
    if n_pairs >= 2:
        rep1_fc = result_df['log2FC_REP1']
        rep2_fc = result_df['log2FC_REP2']
        pr = pearson_r(rep1_fc, rep2_fc)
        sr = spearman_r(rep1_fc, rep2_fc)
        logger.info(f"REP1 vs REP2 log2FC — Pearson r = {pr:.4f}, Spearman rho = {sr:.4f}")

        corr_plot = output_dir / f"{group2}_vs_{group1}_rep_correlation.pdf"
        plot_rep_correlation(rep1_fc, rep2_fc, group1, group2, corr_plot, log2fc_threshold, pr, sr)
    else:
        logger.warning("Only 1 replicate pair found; skipping REP correlation plot")

    # ── log2FC distribution plot ──────────────────────────────────────────────
    dist_plot = output_dir / f"{group2}_vs_{group1}_log2fc_distribution.pdf"
    plot_log2fc_distribution(result_df, group1, group2, dist_plot, log2fc_threshold)

    up   = (result_df['mean_log2FC'] >= log2fc_threshold).sum()
    down = (result_df['mean_log2FC'] <= -log2fc_threshold).sum()
    logger.success(
        f"Done. Clusters up: {up}, down: {down} (|log2FC| >= {log2fc_threshold})"
    )


if __name__ == '__main__':
    app()
