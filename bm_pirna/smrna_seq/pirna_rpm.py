"""Calculate piRNA RPM and fold change between samples."""

from pathlib import Path

from Bio import SeqIO
from loguru import logger
import numpy as np
import pandas as pd
import typer

from bm_pirna.config import EXTERNAL_DATA_DIR, PROCESSED_DATA_DIR, REPORTS_DIR

app = typer.Typer()


def parse_collapsed_fasta(fasta_file: Path) -> dict[str, int]:
    """Parse collapsed fasta file and return piRNA_ID -> count mapping.

    Args:
        fasta_file: Path to collapsed fasta file with format >piRNA_ID_count

    Returns:
        Dict mapping piRNA_ID to read count
    """
    pirna_counts: dict[str, int] = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Header format: >piR-bmo-XXX_count (e.g., >piR-bmo-5_29509)
        header = record.id
        # Split from the last underscore to get piRNA_ID and count
        last_underscore = header.rfind("_")
        if last_underscore != -1:
            pirna_id = header[:last_underscore]
            count = int(header[last_underscore + 1 :])
        else:
            pirna_id = header
            count = 1

        # Aggregate counts for same piRNA_ID (in case of duplicates)
        pirna_counts[pirna_id] = pirna_counts.get(pirna_id, 0) + count

    return pirna_counts


def load_library_sizes(cutadapt_summary: Path) -> dict[str, int]:
    """Load library sizes from cutadapt summary TSV.

    Args:
        cutadapt_summary: Path to cutadapt summary TSV file

    Returns:
        Dict mapping sample name to reads_written (QC-passed reads)
    """
    df = pd.read_csv(cutadapt_summary, sep="\t")
    return dict(zip(df["sample"], df["reads_written"]))


def load_sample_map(sample_map_file: Path) -> list[tuple[str, str]]:
    """Load sample comparison pairs from sample map CSV.

    Args:
        sample_map_file: Path to sample map CSV file

    Returns:
        List of (control_sample, treated_sample) tuples
    """
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
        sample_rpm: Dict mapping sample_name to {pirna_id: rpm}
        sample_groups: Dict mapping base_name to list of sample names

    Returns:
        Dict mapping base_name to averaged {pirna_id: rpm}
    """
    averaged_rpm: dict[str, dict[str, float]] = {}

    for base_name, sample_names in sample_groups.items():
        # Collect all piRNA IDs across replicates
        all_pirnas: set[str] = set()
        for sample_name in sample_names:
            if sample_name in sample_rpm:
                all_pirnas.update(sample_rpm[sample_name].keys())

        # Calculate average RPM for each piRNA
        avg_rpm: dict[str, float] = {}
        for pirna_id in all_pirnas:
            values = [sample_rpm[s].get(pirna_id, 0) for s in sample_names if s in sample_rpm]
            if values:
                avg_rpm[pirna_id] = sum(values) / len(values)

        averaged_rpm[base_name] = avg_rpm

    return averaged_rpm


def calculate_rpm(counts: dict[str, int], library_size: int) -> dict[str, float]:
    """Calculate RPM (Reads Per Million) for each piRNA.

    Args:
        counts: Dict mapping piRNA_ID to read count
        library_size: Total reads in the library

    Returns:
        Dict mapping piRNA_ID to RPM value
    """
    rpm_factor = 1_000_000 / library_size
    return {pirna_id: count * rpm_factor for pirna_id, count in counts.items()}


def calculate_fold_change(
    control_rpm: dict[str, float],
    treated_rpm: dict[str, float],
    pseudocount: float = 0.01,
) -> dict[str, float]:
    """Calculate fold change (treated / control) for each piRNA.

    Args:
        control_rpm: RPM values for control sample
        treated_rpm: RPM values for treated sample
        pseudocount: Small value to avoid division by zero

    Returns:
        Dict mapping piRNA_ID to fold change
    """
    all_pirnas = set(control_rpm.keys()) | set(treated_rpm.keys())
    fold_changes = {}
    for pirna_id in all_pirnas:
        ctrl = control_rpm.get(pirna_id, 0) + pseudocount
        treat = treated_rpm.get(pirna_id, 0) + pseudocount
        fold_changes[pirna_id] = treat / ctrl
    return fold_changes


@app.command()
def main(
    input_dir: Path = typer.Argument(
        ...,
        help="Input directory containing collapsed fasta files",
    ),
    cutadapt_summary: Path = typer.Option(
        REPORTS_DIR / "cutadapt/cutadapt_summary.tsv",
        "--cutadapt-summary",
        "-c",
        help="Path to cutadapt summary TSV file",
    ),
    sample_map: Path = typer.Option(
        EXTERNAL_DATA_DIR / "sample_map.csv",
        "--sample-map",
        "-s",
        help="Path to sample map CSV file",
    ),
    output_dir: Path = typer.Option(
        PROCESSED_DATA_DIR / "pirna_rpm_analysis",
        "--output-dir",
        "-o",
        help="Output directory for results",
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
    """Calculate piRNA RPM and fold change between samples.

    Reads collapsed fasta files, calculates RPM using library sizes from
    fastp summary, and computes fold changes based on sample comparisons.

    Examples:
        python -m bm_pirna.smrna_seq.pirna_rpm collapsed/
        python -m bm_pirna.smrna_seq.pirna_rpm collapsed/ -o results/
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

    # Parse all samples
    sample_counts: dict[str, dict[str, int]] = {}
    sample_rpm: dict[str, dict[str, float]] = {}

    for fasta_file in fasta_files:
        sample_name = fasta_file.stem
        logger.info(f"Processing {sample_name}")

        # Parse counts
        counts = parse_collapsed_fasta(fasta_file)
        sample_counts[sample_name] = counts

        # Get library size (match by sample name prefix)
        lib_size = None
        for lib_sample, size in library_sizes.items():
            if sample_name.startswith(lib_sample.split("_")[0]):
                lib_size = size
                break
        if lib_size is None:
            lib_size = library_sizes.get(sample_name)
        if lib_size is None:
            logger.warning(f"No library size found for {sample_name}, using sum of counts")
            lib_size = sum(counts.values())

        # Calculate RPM
        rpm = calculate_rpm(counts, lib_size)
        sample_rpm[sample_name] = rpm
        logger.info(
            f"  {len(counts)} piRNAs, {sum(counts.values()):,} reads, lib_size={lib_size:,}"
        )

    # Build RPM matrix
    all_pirnas = set()
    for rpm in sample_rpm.values():
        all_pirnas.update(rpm.keys())
    all_pirnas = sorted(all_pirnas)

    rpm_data = {"piRNA_ID": all_pirnas}
    for sample_name, rpm in sample_rpm.items():
        rpm_data[sample_name] = [rpm.get(pirna, 0) for pirna in all_pirnas]

    rpm_df = pd.DataFrame(rpm_data)
    rpm_file = output_dir / "rpm_matrix.tsv"
    rpm_df.to_csv(rpm_file, sep="\t", index=False)
    logger.info(f"Saved RPM matrix: {rpm_file}")

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

        for pirna_id, fold_change in fc.items():
            fc_results.append(
                {
                    "piRNA_ID": pirna_id,
                    "control": control_name,
                    "treated": treated_name,
                    "control_rpm": control_rpm_data.get(pirna_id, 0),
                    "treated_rpm": treated_rpm_data.get(pirna_id, 0),
                    "fold_change": fold_change,
                    "log2_fold_change": np.log2(fold_change) if fold_change > 0 else 0,
                }
            )

    fc_df = pd.DataFrame(fc_results)
    fc_file = output_dir / "fold_change.tsv"
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
