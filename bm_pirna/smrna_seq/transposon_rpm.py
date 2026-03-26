"""Calculate transposon RPM from piRNA alignments and fold change between samples."""

from collections import defaultdict
from pathlib import Path
import subprocess

from loguru import logger
import numpy as np
import pandas as pd
import typer

from bm_pirna.config import (
    EXTERNAL_DATA_DIR,
    INTERIM_DATA_DIR,
    PROCESSED_DATA_DIR,
)

app = typer.Typer()


def parse_collapsed_fasta_counts(fasta_file: Path) -> dict[str, int]:
    """Parse collapsed fasta file and return read_id -> count mapping.

    Args:
        fasta_file: Path to collapsed fasta file with format >rank-count

    Returns:
        Dict mapping read_id (rank-count) to read count
    """
    read_counts: dict[str, int] = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                # Header format: >rank-count (e.g., >1-40608)
                parts = header.split("-")
                if len(parts) >= 2:
                    count = int(parts[-1])
                else:
                    count = 1
                read_counts[header] = count
    return read_counts


def run_bowtie_and_parse(
    fasta_file: Path,
    index: Path,
    threads: int = 16,
    mismatches: int = 1,
) -> list[tuple[str, str]]:
    """Run bowtie alignment and parse SAM output directly.

    Args:
        fasta_file: Path to input fasta file
        index: Path to bowtie index (without extension)
        threads: Number of threads
        mismatches: Number of allowed mismatches

    Returns:
        List of (read_id, transposon_id) tuples for all alignments
    """
    cmd = [
        "bowtie",
        "-p",
        str(threads),
        "-v",
        str(mismatches),
        "-a",  # report all alignments
        "--best",
        "--strata",
        "-x",
        str(index),
        "-f",
        str(fasta_file),
        "-S",
    ]

    logger.debug(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    # Parse SAM output
    alignments = []
    for line in result.stdout.split("\n"):
        if line.startswith("@") or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue

        read_id = fields[0]
        flag = int(fields[1])
        transposon_id = fields[2]

        # Skip unmapped reads (flag & 4)
        if flag & 4 or transposon_id == "*":
            continue

        alignments.append((read_id, transposon_id))

    return alignments


def calculate_transposon_counts(
    alignments: list[tuple[str, str]],
    read_counts: dict[str, int],
) -> dict[str, float]:
    """Calculate transposon counts with multi-mapping reads assigned randomly.

    Args:
        alignments: List of (read_id, transposon_id) tuples
        read_counts: Dict mapping read_id to read count

    Returns:
        Dict mapping transposon_id to count
    """
    import random

    # Group alignments by read_id
    read_to_transposons: dict[str, list[str]] = defaultdict(list)
    for read_id, transposon_id in alignments:
        read_to_transposons[read_id].append(transposon_id)

    # Randomly assign each read to one transposon
    transposon_counts: dict[str, float] = defaultdict(float)
    for read_id, transposon_list in read_to_transposons.items():
        count = read_counts.get(read_id, 1)
        # Randomly pick one transposon from all alignments
        chosen_transposon = random.choice(transposon_list)
        transposon_counts[chosen_transposon] += count

    return dict(transposon_counts)


def calculate_rpm(counts: dict[str, float], library_size: int) -> dict[str, float]:
    """Calculate RPM (Reads Per Million) for each transposon.

    Args:
        counts: Dict mapping transposon_id to count
        library_size: Total reads in the library

    Returns:
        Dict mapping transposon_id to RPM value
    """
    rpm_factor = 1_000_000 / library_size
    return {te_id: count * rpm_factor for te_id, count in counts.items()}


def load_library_sizes(read_counts_file: Path) -> dict[str, int]:
    """Load library sizes from filtered_read_counts TSV."""
    df = pd.read_csv(read_counts_file, sep="\t")
    return dict(zip(df["sample"], df["filtered_reads"]))


def load_sample_map(sample_map_file: Path) -> list[tuple[str, str]]:
    """Load sample comparison pairs from sample map CSV."""
    df = pd.read_csv(sample_map_file)
    pairs = []
    for _, row in df.iterrows():
        control = row["Control_sample"].strip()
        treated = row["Treated_sample"].strip()
        if control and treated:
            pairs.append((control, treated))
    return pairs


def group_samples_by_base_name(sample_names: list[str]) -> dict[str, list[str]]:
    """Group sample names by base name (without _repN suffix).

    Args:
        sample_names: List of sample names (e.g., ['WT_rep1', 'WT_rep2', 'KO_rep1'])

    Returns:
        Dict mapping base name to list of sample names
    """
    import re

    groups: dict[str, list[str]] = {}
    for name in sample_names:
        # Match patterns like _rep1, _rep2, _Rep1, etc.
        match = re.match(r"(.+?)_[Rr]ep\d+$", name)
        if match:
            base_name = match.group(1)
        else:
            base_name = name
        if base_name not in groups:
            groups[base_name] = []
        groups[base_name].append(name)
    return groups


def average_rpm_across_replicates(
    sample_rpm: dict[str, dict[str, float]],
    sample_groups: dict[str, list[str]],
) -> dict[str, dict[str, float]]:
    """Calculate average RPM across replicates for each group.

    Args:
        sample_rpm: Dict mapping sample_name to {te_id: rpm}
        sample_groups: Dict mapping base_name to list of sample names

    Returns:
        Dict mapping base_name to averaged {te_id: rpm}
    """
    averaged_rpm: dict[str, dict[str, float]] = {}

    for base_name, sample_names in sample_groups.items():
        # Collect all transposon IDs across replicates
        all_tes: set[str] = set()
        for sample_name in sample_names:
            if sample_name in sample_rpm:
                all_tes.update(sample_rpm[sample_name].keys())

        # Calculate average RPM for each transposon
        avg_rpm: dict[str, float] = {}
        for te_id in all_tes:
            values = [sample_rpm[s].get(te_id, 0) for s in sample_names if s in sample_rpm]
            if values:
                avg_rpm[te_id] = sum(values) / len(values)

        averaged_rpm[base_name] = avg_rpm

    return averaged_rpm


def calculate_fold_change(
    control_rpm: dict[str, float],
    treated_rpm: dict[str, float],
    pseudocount: float = 0.01,
) -> dict[str, float]:
    """Calculate fold change (treated / control) for each transposon."""
    all_tes = set(control_rpm.keys()) | set(treated_rpm.keys())
    fold_changes = {}
    for te_id in all_tes:
        ctrl = control_rpm.get(te_id, 0) + pseudocount
        treat = treated_rpm.get(te_id, 0) + pseudocount
        fold_changes[te_id] = treat / ctrl
    return fold_changes


@app.command()
def main(
    input_dir: Path = typer.Argument(
        ...,
        help="Input directory containing collapsed fasta files",
    ),
    bowtie_index: Path = typer.Option(
        INTERIM_DATA_DIR / "smRNA-seq_20260122/bowtie_index/transposon_index",
        "--index",
        "-i",
        help="Path to bowtie index (without extension)",
    ),
    cutadapt_summary: Path = typer.Option(
        INTERIM_DATA_DIR / "smRNA-seq_20260321/structure_rna_filtered/filtered_read_counts.tsv",
        "--cutadapt-summary",
        "-c",
        help="Path to filtered_read_counts TSV file (filtered_reads column used as normalization factor)",
    ),
    sample_map: Path = typer.Option(
        EXTERNAL_DATA_DIR / "sample_map.csv",
        "--sample-map",
        "-s",
        help="Path to sample map CSV file",
    ),
    output_dir: Path = typer.Option(
        PROCESSED_DATA_DIR / "transposon_rpm_analysis",
        "--output-dir",
        "-o",
        help="Output directory for results",
    ),
    threads: int = typer.Option(
        16,
        "--threads",
        "-t",
        help="Number of threads for bowtie",
    ),
    mismatches: int = typer.Option(
        1,
        "--mismatches",
        "-m",
        help="Number of allowed mismatches",
    ),
    pseudocount: float = typer.Option(
        0.01,
        "--pseudocount",
        "-p",
        help="Pseudocount for fold change calculation",
    ),
    pattern: str = typer.Option(
        "*.fa",
        "--pattern",
        help="Glob pattern to match fasta files",
    ),
) -> None:
    """Calculate transposon RPM from piRNA alignments and fold change.

    Aligns collapsed fasta files to transposon index using bowtie,
    calculates RPM using library sizes from fastp summary,
    and computes fold changes based on sample comparisons.

    For multi-mapping reads, counts are randomly assigned to one alignment.

    Examples:
        python -m bm_pirna.smrna_seq.transposon_rpm collapsed/
        python -m bm_pirna.smrna_seq.transposon_rpm collapsed/ -t 32 -m 1
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load library sizes
    if not cutadapt_summary.exists():
        logger.error(f"Cutadapt summary not found: {cutadapt_summary}")
        raise typer.Exit(1)
    library_sizes = load_library_sizes(cutadapt_summary)
    logger.info(f"Loaded library sizes for {len(library_sizes)} samples")

    # Load sample comparisons
    if not sample_map.exists():
        logger.error(f"Sample map not found: {sample_map}")
        raise typer.Exit(1)
    comparisons = load_sample_map(sample_map)
    logger.info(f"Loaded {len(comparisons)} sample comparisons")

    # Get input files
    if not input_dir.is_dir():
        logger.error(f"Input directory does not exist: {input_dir}")
        raise typer.Exit(1)
    fasta_files = sorted(input_dir.glob(pattern))
    if not fasta_files:
        logger.error(f"No fasta files found in {input_dir} with pattern '{pattern}'")
        raise typer.Exit(1)
    logger.info(f"Found {len(fasta_files)} fasta files")

    # Process all samples
    sample_counts: dict[str, dict[str, float]] = {}
    sample_rpm: dict[str, dict[str, float]] = {}

    for fasta_file in fasta_files:
        sample_name = fasta_file.stem
        logger.info(f"Processing {sample_name}")

        # Parse read counts from fasta
        read_counts = parse_collapsed_fasta_counts(fasta_file)
        total_reads = sum(read_counts.values())
        logger.info(f"  {len(read_counts)} unique reads, {total_reads:,} total reads")

        # Run bowtie and get alignments
        alignments = run_bowtie_and_parse(
            fasta_file, bowtie_index, threads=threads, mismatches=mismatches
        )
        logger.info(f"  {len(alignments)} alignments")

        # Calculate transposon counts (with multi-mapping split)
        te_counts = calculate_transposon_counts(alignments, read_counts)
        sample_counts[sample_name] = te_counts
        logger.info(f"  {len(te_counts)} transposons hit")

        # Get library size
        lib_size = library_sizes.get(sample_name)
        if lib_size is None:
            logger.warning(f"No library size found for {sample_name}, using total reads")
            lib_size = total_reads

        # Calculate RPM
        rpm = calculate_rpm(te_counts, lib_size)
        sample_rpm[sample_name] = rpm

    # Build RPM matrix
    all_tes = set()
    for rpm in sample_rpm.values():
        all_tes.update(rpm.keys())
    all_tes = sorted(all_tes)

    rpm_data = {"transposon_ID": all_tes}
    for sample_name, rpm in sample_rpm.items():
        rpm_data[sample_name] = [rpm.get(te, 0) for te in all_tes]

    rpm_df = pd.DataFrame(rpm_data)
    rpm_file = output_dir / "transposon_rpm_matrix.tsv"
    rpm_df.to_csv(rpm_file, sep="\t", index=False)
    logger.info(f"Saved RPM matrix: {rpm_file}")

    # Also save raw counts matrix
    count_data = {"transposon_ID": all_tes}
    for sample_name, counts in sample_counts.items():
        count_data[sample_name] = [counts.get(te, 0) for te in all_tes]

    count_df = pd.DataFrame(count_data)
    count_file = output_dir / "transposon_count_matrix.tsv"
    count_df.to_csv(count_file, sep="\t", index=False)
    logger.info(f"Saved count matrix: {count_file}")

    # Group samples by base name and calculate average RPM
    sample_groups = group_samples_by_base_name(list(sample_rpm.keys()))
    logger.info(f"Detected {len(sample_groups)} sample groups from {len(sample_rpm)} samples")
    for base_name, members in sample_groups.items():
        if len(members) > 1:
            logger.info(f"  {base_name}: {len(members)} replicates -> {members}")

    averaged_rpm = average_rpm_across_replicates(sample_rpm, sample_groups)

    # Calculate fold changes for each comparison (using averaged RPM)
    fc_results = []
    for control_name, treated_name in comparisons:
        # Try to find matching group (exact match or base name match)
        control_rpm_data = averaged_rpm.get(control_name)
        treated_rpm_data = averaged_rpm.get(treated_name)

        if control_rpm_data is None:
            logger.warning(f"Control sample not found: {control_name}")
            continue
        if treated_rpm_data is None:
            logger.warning(f"Treated sample not found: {treated_name}")
            continue

        fc = calculate_fold_change(
            control_rpm_data,
            treated_rpm_data,
            pseudocount=pseudocount,
        )

        for te_id, fold_change in fc.items():
            fc_results.append(
                {
                    "transposon_ID": te_id,
                    "control": control_name,
                    "treated": treated_name,
                    "control_rpm": control_rpm_data.get(te_id, 0),
                    "treated_rpm": treated_rpm_data.get(te_id, 0),
                    "fold_change": fold_change,
                    "log2_fold_change": np.log2(fold_change) if fold_change > 0 else 0,
                }
            )

    fc_df = pd.DataFrame(fc_results)
    fc_file = output_dir / "transposon_fold_change.tsv"
    fc_df.to_csv(fc_file, sep="\t", index=False)
    logger.info(f"Saved fold change results: {fc_file}")

    # Create summary per comparison
    for control_name, treated_name in comparisons:
        comparison_df = fc_df[
            (fc_df["control"] == control_name) & (fc_df["treated"] == treated_name)
        ].copy()
        if comparison_df.empty:
            continue

        comparison_df = comparison_df.sort_values("fold_change", ascending=False)
        comparison_file = output_dir / f"fc_{treated_name}_vs_{control_name}.tsv"
        comparison_df.to_csv(comparison_file, sep="\t", index=False)
        logger.info(f"Saved: {comparison_file}")

        # Summary stats
        up_regulated = (comparison_df["fold_change"] > 2).sum()
        down_regulated = (comparison_df["fold_change"] < 0.5).sum()
        logger.info(
            f"  {treated_name} vs {control_name}: "
            f"{up_regulated} up (FC>2), {down_regulated} down (FC<0.5)"
        )

    logger.success(f"All results saved to: {output_dir}")


if __name__ == "__main__":
    app()
