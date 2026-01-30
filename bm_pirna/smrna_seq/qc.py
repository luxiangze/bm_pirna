"""Quality control for smRNA-seq data using cutadapt."""

import json
from pathlib import Path
import subprocess

from loguru import logger
import pandas as pd
import typer

from bm_pirna.config import EXTERNAL_DATA_DIR, INTERIM_DATA_DIR, REPORTS_DIR

app = typer.Typer()

DEFAULT_ADAPTER = "TGGAATTCTCGGGTGCCAAGG"


def run_fastqc(
    input_file: Path,
    output_dir: Path,
    threads: int = 16,
) -> subprocess.CompletedProcess:
    """Run FastQC on raw data before adapter trimming."""
    cmd = [
        "fastqc",
        str(input_file),
        "-o",
        str(output_dir),
        "-t",
        str(threads),
        "--quiet",
    ]

    logger.debug(f"Running: {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def run_cutadapt(
    input_file: Path,
    output_file: Path,
    json_report: Path | None = None,
    adapter: str | None = None,
    min_length: int = 18,
    max_length: int = 37,
    discard_untrimmed: bool = False,
    cores: int = 16,
) -> subprocess.CompletedProcess:
    """Run cutadapt for adapter trimming and length filtering."""
    cmd = [
        "cutadapt",
        "-q 20",
        "-m",
        str(min_length),
        "-M",
        str(max_length),
        "-j",
        str(cores),
        "-o",
        str(output_file),
    ]

    if json_report:
        cmd.extend(["--json", str(json_report)])

    if adapter:
        cmd.extend(["-a", adapter])
        if discard_untrimmed:
            cmd.append("--discard-untrimmed")

    cmd.append(str(input_file))

    logger.debug(f"Running: {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def run_seqkit_length_filter(
    input_file: Path,
    output_file: Path,
    min_length: int,
    max_length: int,
    threads: int = 16,
) -> subprocess.CompletedProcess:
    """Run seqkit seq for length filtering."""
    cmd = [
        "seqkit",
        "seq",
        "-m",
        str(min_length),
        "-M",
        str(max_length),
        "-j",
        str(threads),
        "-o",
        str(output_file),
        str(input_file),
    ]

    logger.debug(f"Running: {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def run_multiqc(
    input_dir: Path,
    output_dir: Path,
    title: str = "smRNA-seq QC Report",
) -> subprocess.CompletedProcess:
    """Run MultiQC to aggregate reports."""
    cmd = [
        "multiqc",
        str(input_dir),
        "-o",
        str(output_dir),
        "-n",
        "multiqc_report",
        "--title",
        title,
        "--force",
    ]

    logger.debug(f"Running: {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def parse_cutadapt_json(json_file: Path) -> dict:
    """Parse cutadapt JSON report and extract key statistics.

    Args:
        json_file: Path to cutadapt JSON report

    Returns:
        Dictionary with key statistics
    """
    with open(json_file) as f:
        data = json.load(f)

    read_counts = data.get("read_counts", {})
    basepair_counts = data.get("basepair_counts", {})

    total_reads = read_counts.get("input", 0) or 0
    reads_with_adapters = read_counts.get("read1_with_adapter", 0) or 0
    filtered = read_counts.get("filtered", {}) or {}
    reads_too_short = filtered.get("too_short", 0) or 0
    reads_too_long = filtered.get("too_long", 0) or 0
    reads_written = read_counts.get("output", 0) or 0
    total_bp = basepair_counts.get("input", 0) or 0
    written_bp = basepair_counts.get("output", 0) or 0

    def pct(num: int, denom: int) -> float:
        return (num / denom * 100) if denom > 0 else 0.0

    return {
        "total_reads": total_reads,
        "reads_with_adapters": reads_with_adapters,
        "reads_with_adapters_pct": pct(reads_with_adapters, total_reads),
        "reads_too_short": reads_too_short,
        "reads_too_short_pct": pct(reads_too_short, total_reads),
        "reads_too_long": reads_too_long,
        "reads_too_long_pct": pct(reads_too_long, total_reads),
        "reads_written": reads_written,
        "reads_written_pct": pct(reads_written, total_reads),
        "total_bp": total_bp,
        "written_bp": written_bp,
        "written_bp_pct": pct(written_bp, total_bp),
    }


def write_cutadapt_summary(stats_list: list[tuple[str, dict]], output_file: Path) -> None:
    """Write summary statistics from cutadapt to TSV.

    Args:
        stats_list: List of (sample_name, stats_dict) tuples
        output_file: Output TSV file
    """
    with open(output_file, "w") as f:
        f.write(
            "sample\ttotal_reads\treads_with_adapters\treads_with_adapters_pct\t"
            "reads_too_short\treads_too_short_pct\treads_too_long\treads_too_long_pct\t"
            "reads_written\treads_written_pct\ttotal_bp\twritten_bp\twritten_bp_pct\n"
        )

        for sample_name, stats in sorted(stats_list, key=lambda x: x[0]):
            f.write(
                f"{sample_name}\t{stats['total_reads']}\t{stats['reads_with_adapters']}\t"
                f"{stats['reads_with_adapters_pct']:.1f}\t{stats['reads_too_short']}\t"
                f"{stats['reads_too_short_pct']:.1f}\t{stats['reads_too_long']}\t"
                f"{stats['reads_too_long_pct']:.1f}\t{stats['reads_written']}\t"
                f"{stats['reads_written_pct']:.1f}\t{stats['total_bp']}\t"
                f"{stats['written_bp']}\t{stats['written_bp_pct']:.1f}\n"
            )


def qc_single_file(
    input_file: Path,
    output_file: Path,
    fastqc_dir: Path | None = None,
    cutadapt_dir: Path | None = None,
    adapter: str | None = None,
    min_length: int = 18,
    max_length: int = 37,
    discard_untrimmed: bool = False,
    cores: int = 16,
) -> dict:
    """Run quality control on a single fastq file.

    Process:
    1. fastqc: raw data quality check (optional)
    2. cutadapt: adapter trimming (if adapter provided) and length filtering (18-37bp)

    Returns:
        Dictionary with cutadapt statistics
    """
    sample_name = input_file.stem.split(".")[0]

    # Step 1: FastQC on raw data (before adapter trimming)
    if fastqc_dir:
        fastqc_dir.mkdir(parents=True, exist_ok=True)
        fastqc_html = fastqc_dir / f"{sample_name}_fastqc.html"
        if fastqc_html.exists():
            logger.debug(f"[fastqc] Skipping {sample_name} (report exists)")
        else:
            logger.info(f"[fastqc] Checking raw data quality: {input_file.name}")
            run_fastqc(input_file=input_file, output_dir=fastqc_dir, threads=cores)

    # Step 2: cutadapt for adapter trimming and length filtering
    if adapter:
        logger.info(f"[cutadapt] Trimming adapter and filtering length from {input_file.name}")
    else:
        logger.info(f"[cutadapt] Filtering length from {input_file.name}")

    json_report = cutadapt_dir / f"{sample_name}.cutadapt.json" if cutadapt_dir else None
    if cutadapt_dir:
        cutadapt_dir.mkdir(parents=True, exist_ok=True)

    run_cutadapt(
        input_file=input_file,
        output_file=output_file,
        json_report=json_report,
        adapter=adapter,
        min_length=min_length,
        max_length=max_length,
        discard_untrimmed=discard_untrimmed,
        cores=cores,
    )

    stats = parse_cutadapt_json(json_report) if json_report and json_report.exists() else {}
    logger.success(f"QC complete: {output_file}")

    return stats


@app.command()
def main(
    input_path: Path = typer.Argument(
        ...,
        help="Input fastq file or directory containing fastq files",
    ),
    output_dir: Path = typer.Option(
        INTERIM_DATA_DIR / "qc",
        "--output-dir",
        "-o",
        help="Output directory for QC'd files",
    ),
    adapter_file: Path = typer.Option(
        EXTERNAL_DATA_DIR / "sample_adapter.csv",
        "--adapter-file",
        "-a",
        help="CSV file with sample-specific adapters (columns: sample, adapter, rm_adapter)",
    ),
    min_length: int = typer.Option(
        18,
        "--min-length",
        "-m",
        help="Minimum read length to keep",
    ),
    max_length: int = typer.Option(
        37,
        "--max-length",
        "-M",
        help="Maximum read length to keep",
    ),
    length_filter: str = typer.Option(
        None,
        "--length-filter",
        "-l",
        help="Optional secondary length filter range (e.g., '24-35'). If set, generates additional filtered file.",
    ),
    cores: int = typer.Option(
        16,
        "--cores",
        "-j",
        help="Number of CPU cores to use",
    ),
    pattern: str = typer.Option(
        "*.f*q.gz",
        "--pattern",
        "-p",
        help="Glob pattern to match fastq files (when input is a directory)",
    ),
    report_dir: Path = typer.Option(
        REPORTS_DIR,
        "--report-dir",
        "-r",
        help="Base directory for QC reports (cutadapt/, fastqc/, multiqc/ subdirs)",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Force re-run even if output files exist",
    ),
) -> None:
    """Quality control for smRNA-seq data using FastQC and cutadapt.

    Process:
    1. fastqc: raw data quality check
    2. cutadapt: adapter trimming (if rm_adapter=true) and length filtering (18-37bp)
    3. seqkit (optional): secondary length filtering if --length-filter is specified
    4. multiqc: aggregate all FastQC reports

    Output files:
    - {sample}.fq.gz: QC'd reads (always generated)
    - {min}-{max}_length_filtered/{sample}.fq.gz: secondary filtered reads (only if --length-filter is set)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup report subdirectories
    cutadapt_dir = report_dir / "cutadapt"
    fastqc_dir = report_dir / "fastqc"
    cutadapt_dir.mkdir(parents=True, exist_ok=True)
    fastqc_dir.mkdir(parents=True, exist_ok=True)

    # Parse length filter
    length_filter_tuple: tuple[int, int] | None = None
    filtered_dir: Path | None = None
    if length_filter:
        try:
            parts = length_filter.split("-")
            if len(parts) != 2:
                raise ValueError
            length_filter_tuple = (int(parts[0]), int(parts[1]))
            filtered_dir = output_dir / f"{length_filter}_length_filtered"
            filtered_dir.mkdir(parents=True, exist_ok=True)
            logger.info(
                f"Secondary length filter enabled: {length_filter_tuple[0]}-{length_filter_tuple[1]}bp"
            )
            logger.info(f"Filtered output dir: {filtered_dir}")
        except ValueError:
            logger.error(
                f"Invalid length filter format: '{length_filter}'. Expected format: 'min-max' (e.g., '24-35')"
            )
            raise typer.Exit(1)

    if input_path.is_file():
        fastq_files = [input_path]
    elif input_path.is_dir():
        fastq_files = sorted(input_path.glob(pattern))
        if not fastq_files:
            logger.error(f"No fastq files found in {input_path} with pattern '{pattern}'")
            raise typer.Exit(1)
    else:
        logger.error(f"Input path does not exist: {input_path}")
        raise typer.Exit(1)

    # Load adapter mapping from CSV (columns: sample, adapter, rm_adapter)
    adapter_map: dict[str, tuple[str, bool]] = {}
    if adapter_file.exists():
        adapter_df = pd.read_csv(adapter_file)
        for _, row in adapter_df.iterrows():
            rm_adapter = str(row.get("rm_adapter", "true")).lower() == "true"
            adapter_map[row["sample"]] = (row["adapter"], rm_adapter)
        logger.info(f"Loaded {len(adapter_map)} sample configs from {adapter_file}")
    else:
        logger.warning(f"Adapter file not found: {adapter_file}")

    logger.info(f"Processing {len(fastq_files)} file(s)")

    processed_count = 0
    skipped_count = 0
    stats_list: list[tuple[str, dict]] = []

    for fastq_file in fastq_files:
        sample_name = fastq_file.stem.split(".")[0]
        output_file = output_dir / f"{sample_name}.fq.gz"
        filtered_file = filtered_dir / f"{sample_name}.fq.gz" if filtered_dir else None

        # Get sample-specific config
        adapter_config = adapter_map.get(sample_name)
        if adapter_config:
            adapter, rm_adapter = adapter_config
        else:
            adapter, rm_adapter = None, False
            logger.warning(f"No config for {sample_name}, only length filtering")

        # Check if outputs already exist
        outputs_exist = output_file.exists()
        if length_filter_tuple and filtered_file:
            outputs_exist = outputs_exist and filtered_file.exists()

        if not force and outputs_exist:
            logger.info(f"Skipping {sample_name} (outputs exist)")
            skipped_count += 1
            continue

        if rm_adapter:
            logger.info(f"[{sample_name}] Adapter: {adapter}, discard_untrimmed=True")
        else:
            logger.info(f"[{sample_name}] No adapter trimming, length filter only")

        try:
            stats = qc_single_file(
                input_file=fastq_file,
                output_file=output_file,
                fastqc_dir=fastqc_dir,
                cutadapt_dir=cutadapt_dir,
                adapter=adapter if rm_adapter else None,
                min_length=min_length,
                max_length=max_length,
                discard_untrimmed=rm_adapter,
                cores=cores,
            )
            stats_list.append((sample_name, stats))

            # Secondary length filtering with seqkit
            if length_filter_tuple and filtered_file:
                min_len, max_len = length_filter_tuple
                logger.info(
                    f"[seqkit] Secondary length filtering {sample_name} ({min_len}-{max_len}bp)"
                )
                run_seqkit_length_filter(
                    input_file=output_file,
                    output_file=filtered_file,
                    min_length=min_len,
                    max_length=max_len,
                    threads=cores,
                )
                logger.success(f"Secondary filtered: {filtered_file}")

            processed_count += 1
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to process {fastq_file.name}: {e.stderr}")
            raise typer.Exit(1)

    logger.info(f"Processed: {processed_count}, Skipped: {skipped_count}")

    # Generate cutadapt summary statistics
    if stats_list:
        summary_file = cutadapt_dir / "cutadapt_summary.tsv"
        logger.info("[summary] Generating cutadapt summary statistics")
        write_cutadapt_summary(stats_list, summary_file)
        logger.success(f"Summary statistics: {summary_file}")

    # Aggregate all reports with MultiQC
    multiqc_output_dir = report_dir / "multiqc"
    multiqc_output_dir.mkdir(parents=True, exist_ok=True)
    multiqc_report = multiqc_output_dir / "multiqc_report.html"

    if not force and multiqc_report.exists() and processed_count == 0:
        logger.info("Skipping MultiQC (report exists, no new samples)")
    else:
        logger.info("[multiqc] Aggregating FastQC and cutadapt reports")
        try:
            cmd = [
                "multiqc",
                str(fastqc_dir),
                str(cutadapt_dir),
                "-o",
                str(multiqc_output_dir),
                "-n",
                "multiqc_report",
                "--title",
                "smRNA-seq QC Report",
                "--force",
            ]
            logger.debug(f"Running: {' '.join(cmd)}")
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.success(f"MultiQC report: {multiqc_report}")
        except subprocess.CalledProcessError as e:
            logger.warning(f"MultiQC failed (non-critical): {e.stderr}")

    logger.success(f"All QC complete. Output: {output_dir}")


if __name__ == "__main__":
    app()
