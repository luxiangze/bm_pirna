#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Summarize rMATS analysis results and generate a single scatter plot

This script is used for:
1. Reading rMATS analysis results from a single comparison
2. Combining all splicing event types into one dataset
3. Generating a single scatter plot (volcano plot)
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Configure matplotlib to use standard fonts to avoid encoding issues
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'Helvetica', 'sans-serif']
mpl.rcParams['axes.unicode_minus'] = False


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Summarize rMATS analysis results and generate a scatter plot')
    parser.add_argument('--input-dir', required=True, help='rMATS results directory (single comparison)')
    parser.add_argument('--output-file', required=True, help='Output summary file path')
    parser.add_argument('--plot-file', required=True, help='Output plot file path (PDF)')
    parser.add_argument('--fdr', type=float, default=0.05, help='FDR threshold (default: 0.05)')
    parser.add_argument('--annotation-file', default=None, help='piRNA annotation CSV with is_known and Gene_ID columns')
    return parser.parse_args()


def read_rmats_file(file_path):
    """
    Read rMATS result file

    Args:
        file_path: Path to rMATS result file

    Returns:
        DataFrame: DataFrame containing required fields
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        # Keep only key columns
        if 'FDR' in df.columns:
            key_columns = ['ID', 'GeneID', 'geneSymbol',
                          'PValue', 'FDR', 'IncLevel1', 'IncLevel2',
                          'IncLevelDifference']
            return df[key_columns]
        return None
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None


def collect_results(input_dir):
    """
    Collect rMATS results from a single comparison directory

    Args:
        input_dir: rMATS results directory (single comparison)

    Returns:
        DataFrame: DataFrame containing all events with event type information
    """
    all_data = []
    event_types = ['A3SS', 'A5SS', 'MXE', 'RI', 'SE']

    for event in event_types:
        jc_file = os.path.join(input_dir, f"{event}.MATS.JC.txt")
        if os.path.exists(jc_file):
            df = read_rmats_file(jc_file)
            if df is not None:
                df['EventType'] = event
                all_data.append(df)

    if all_data:
        return pd.concat(all_data, ignore_index=True)
    else:
        return pd.DataFrame()


def parse_inc_level(val):
    """Parse IncLevel value (comma-separated) to mean"""
    if pd.isna(val) or val == 'NA' or val == '':
        return np.nan
    try:
        vals = [float(x) for x in str(val).split(',') if x.strip() and x != 'NA']
        return np.mean(vals) if vals else np.nan
    except (ValueError, TypeError):
        return np.nan


def load_known_pirnas(annotation_file):
    """
    Load known piRNA gene IDs from annotation file

    Args:
        annotation_file: Path to piRNA annotation CSV

    Returns:
        dict: mapping Gene_ID -> sample_id for is_known=True entries
    """
    ann = pd.read_csv(annotation_file)
    ann['is_known'] = ann['is_known'].astype(str).str.lower() == 'true'
    known = ann[ann['is_known']]
    # Filter out rows with NaN Gene_ID
    known = known[known['Gene_ID'].notna()]
    # Convert Gene_ID to string
    known['Gene_ID'] = known['Gene_ID'].astype(str)
    return dict(zip(known['Gene_ID'], known['sample_id']))


def plot_volcano(df, output_file, fdr_threshold=0.05, known_pirnas=None):
    """
    Create a volcano plot for all splicing events

    Args:
        df: DataFrame containing all events
        output_file: Output file path (PDF)
        fdr_threshold: FDR threshold (default: 0.05)
        known_pirnas: dict mapping Gene_ID -> label for known piRNAs
    """
    if df.empty or 'FDR' not in df.columns:
        print("No suitable data for volcano plot")
        return

    data = df.copy()
    # Add -log10(FDR) column, handle FDR=0 case
    data['FDR'] = data['FDR'].replace(0, 1e-300)
    data['-log10(FDR)'] = -np.log10(data['FDR'])

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Define colors: upregulated (red), downregulated (blue), non-significant (gray)
    colors = []
    for _, row in data.iterrows():
        if row['FDR'] < fdr_threshold and row['IncLevelDifference'] > 0:
            colors.append('red')
        elif row['FDR'] < fdr_threshold and row['IncLevelDifference'] < 0:
            colors.append('blue')
        else:
            colors.append('gray')
    data['color'] = colors

    # Draw scatter plot by color groups
    for color, label in [('gray', 'Non-significant'),
                         ('blue', 'Down-regulated'),
                         ('red', 'Up-regulated')]:
        subset = data[data['color'] == color]
        ax.scatter(subset['IncLevelDifference'], subset['-log10(FDR)'],
                  c=color, alpha=0.6, s=20, label=f'{label} ({len(subset)})')

    # Add significance threshold line
    ax.axhline(y=-np.log10(fdr_threshold), color='black', linestyle='--', alpha=0.3)

    # Add labels (no title, larger font)
    ax.set_xlabel('Inclusion Level Difference (ΔΨ)', fontsize=16)
    ax.set_ylabel('-log₁₀(FDR)', fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(loc='upper right', fontsize=12)

    # Count significant events
    sig_up = len(data[(data['FDR'] < fdr_threshold) & (data['IncLevelDifference'] > 0)])
    sig_down = len(data[(data['FDR'] < fdr_threshold) & (data['IncLevelDifference'] < 0)])
    total = len(data)

    # Add statistics text box
    textstr = f'Total: {total}\nUp: {sig_up}\nDown: {sig_down}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=12,
           verticalalignment='top', bbox=props)

    # Label known piRNA points
    if known_pirnas:
        data['GeneID_str'] = data['GeneID'].astype(str)
        # Match GeneID directly with Gene_ID from annotation file
        labeled = data[data['GeneID_str'].isin(known_pirnas)]

        # Add labels with better visibility
        for _, row in labeled.iterrows():
            label = known_pirnas[row['GeneID_str']]
            ax.annotate(label,
                       xy=(row['IncLevelDifference'], row['-log10(FDR)']),
                       xytext=(0, 10), textcoords='offset points',
                       fontsize=9, ha='center', va='bottom',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7, edgecolor='none'),
                       arrowprops=dict(arrowstyle='->', lw=1, color='black', alpha=0.8))

    # Save figure as PDF
    plt.tight_layout()
    plt.savefig(output_file, format='pdf')
    plt.close()
    print(f"Saved plot to {output_file}")


def main():
    """Main function"""
    args = parse_args()

    # Ensure output directories exist
    output_dir = os.path.dirname(args.output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    plot_dir = os.path.dirname(args.plot_file)
    if plot_dir:
        os.makedirs(plot_dir, exist_ok=True)

    print("Reading rMATS results...")
    all_results = collect_results(args.input_dir)

    if all_results.empty:
        print("No rMATS results found or failed to read results")
        return

    print(f"Total events: {len(all_results)}")

    # Save summary results
    print(f"Saving summary to {args.output_file}")
    all_results.to_csv(args.output_file, sep='\t', index=False)

    # Load piRNA annotations if provided
    known_pirnas = None
    if args.annotation_file:
        print(f"Loading piRNA annotations from {args.annotation_file}")
        known_pirnas = load_known_pirnas(args.annotation_file)
        print(f"Found {len(known_pirnas)} known piRNA genes for labeling")

    # Generate volcano plot
    print("Generating volcano plot...")
    plot_volcano(all_results, args.plot_file, args.fdr, known_pirnas)

    print("Done!")


if __name__ == "__main__":
    main()
