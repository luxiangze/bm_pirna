"""Quality control for smRNA-seq data using cutadapt and fastp."""

import json
from pathlib import Path
import subprocess
import tempfile

from loguru import logger
import typer

from bm_pirna.config import INTERIM_DATA_DIR, REPORTS_DIR

app = typer.Typer()

DEFAULT_ADAPTER = "TGGAATTCTCGGGTGCCAAGG"


def run_cutadapt(
    input_file: Path,
    output_file: Path,
    adapter: str,
    overlap: int = 5,
    discard_untrimmed: bool = True,
    cores: int = 16,
) -> subprocess.CompletedProcess:
    """Run cutadapt for adapter trimming only."""
    cmd = [
        "cutadapt",
        "-a",
        adapter,
        "--overlap",
        str(overlap),
        "-j",
        str(cores),
        "-o",
        str(output_file),
    ]

    if discard_untrimmed:
        cmd.append("--discard-untrimmed")

    cmd.append(str(input_file))

    logger.debug(f"Running: {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def run_fastp(
    input_file: Path,
    output_file: Path,
    html_report: Path,
    json_report: Path,
    quality: int = 20,
    min_length: int = 18,
    max_length: int = 37,
    max_n: int = 0,
    cores: int = 16,
) -> subprocess.CompletedProcess:
    """Run fastp for quality control with default length filtering (18-37bp)."""
    cmd = [
        "fastp",
        "-i",
        str(input_file),
        "-o",
        str(output_file),
        "-q",
        str(quality),
        "-l",
        str(min_length),
        "--length_limit",
        str(max_length),
        "-n",
        str(max_n),
        "-w",
        str(cores),
        "-A",  # disable adapter trimming (already done by cutadapt)
        "--html",
        str(html_report),
        "--json",
        str(json_report),
    ]

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
    """Run MultiQC to aggregate fastp reports."""
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


def parse_fastp_json(json_file: Path) -> dict:
    """Parse fastp JSON report and extract key statistics.

    Args:
        json_file: Path to fastp JSON report

    Returns:
        Dictionary with key statistics
    """
    with open(json_file) as f:
        data = json.load(f)

    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    filtering = data.get("filtering_result", {})
    duplication = data.get("duplication", {})

    return {
        "before_total_reads": before.get("total_reads", 0),
        "before_total_bases": before.get("total_bases", 0),
        "before_q20_rate": before.get("q20_rate", 0) * 100,
        "before_q30_rate": before.get("q30_rate", 0) * 100,
        "before_gc_content": before.get("gc_content", 0) * 100,
        "after_total_reads": after.get("total_reads", 0),
        "after_total_bases": after.get("total_bases", 0),
        "after_q20_rate": after.get("q20_rate", 0) * 100,
        "after_q30_rate": after.get("q30_rate", 0) * 100,
        "after_gc_content": after.get("gc_content", 0) * 100,
        "after_mean_length": after.get("read1_mean_length", 0),
        "passed_filter_reads": filtering.get("passed_filter_reads", 0),
        "low_quality_reads": filtering.get("low_quality_reads", 0),
        "too_many_N_reads": filtering.get("too_many_N_reads", 0),
        "too_short_reads": filtering.get("too_short_reads", 0),
        "too_long_reads": filtering.get("too_long_reads", 0),
        "duplication_rate": duplication.get("rate", 0) * 100,
    }


def write_fastp_summary(json_files: list[Path], output_file: Path) -> None:
    """Write summary statistics from fastp JSON reports to TSV.

    Args:
        json_files: List of fastp JSON report files
        output_file: Output TSV file
    """
    with open(output_file, "w") as f:
        # Header
        f.write(
            "sample\tbefore_total_reads\tbefore_total_bases\tbefore_q20_rate\t"
            "before_q30_rate\tbefore_gc_content\tafter_total_reads\tafter_total_bases\t"
            "after_q20_rate\tafter_q30_rate\tafter_gc_content\tafter_mean_length\t"
            "passed_filter_reads\tlow_quality_reads\ttoo_many_N_reads\t"
            "too_short_reads\ttoo_long_reads\tduplication_rate\n"
        )

        for json_file in sorted(json_files):
            sample_name = json_file.stem
            try:
                stats = parse_fastp_json(json_file)
                f.write(
                    f"{sample_name}\t{stats['before_total_reads']}\t{stats['before_total_bases']}\t"
                    f"{stats['before_q20_rate']:.2f}\t{stats['before_q30_rate']:.2f}\t"
                    f"{stats['before_gc_content']:.2f}\t{stats['after_total_reads']}\t"
                    f"{stats['after_total_bases']}\t{stats['after_q20_rate']:.2f}\t"
                    f"{stats['after_q30_rate']:.2f}\t{stats['after_gc_content']:.2f}\t"
                    f"{stats['after_mean_length']:.2f}\t{stats['passed_filter_reads']}\t"
                    f"{stats['low_quality_reads']}\t{stats['too_many_N_reads']}\t"
                    f"{stats['too_short_reads']}\t{stats['too_long_reads']}\t"
                    f"{stats['duplication_rate']:.2f}\n"
                )
            except Exception as e:
                logger.warning(f"Failed to parse {json_file.name}: {e}")


def qc_single_file(
    input_file: Path,
    output_file: Path,
    report_dir: Path,
    adapter: str = DEFAULT_ADAPTER,
    quality: int = 20,
    overlap: int = 5,
    discard_untrimmed: bool = True,
    cores: int = 16,
    length_filter: tuple[int, int] | None = None,
    filtered_output_file: Path | None = None,
) -> None:
    """Run quality control on a single fastq file.

    Process:
    1. cutadapt: adapter trimming
    2. fastp: quality filtering and QC report generation
    3. seqkit (optional): length filtering if length_filter is specified
    """
    sample_name = input_file.stem.split(".")[0]

    with tempfile.TemporaryDirectory() as tmpdir:
        trimmed_file = Path(tmpdir) / f"{sample_name}.trim.fq.gz"

        # Step 1: Adapter trimming with cutadapt
        logger.info(f"[cutadapt] Trimming adapter from {input_file.name}")
        run_cutadapt(
            input_file=input_file,
            output_file=trimmed_file,
            adapter=adapter,
            overlap=overlap,
            discard_untrimmed=discard_untrimmed,
            cores=cores,
        )

        # Step 2: Quality filtering with fastp
        logger.info(f"[fastp] Quality filtering {sample_name}")
        report_dir.mkdir(parents=True, exist_ok=True)
        html_report = report_dir / f"{sample_name}.html"
        json_report = report_dir / f"{sample_name}.json"

        run_fastp(
            input_file=trimmed_file,
            output_file=output_file,
            html_report=html_report,
            json_report=json_report,
            quality=quality,
            cores=cores,
        )

    logger.success(f"QC complete: {output_file}")

    # Step 3: Optional length filtering with seqkit
    if length_filter and filtered_output_file:
        min_len, max_len = length_filter
        logger.info(f"[seqkit] Length filtering {sample_name} ({min_len}-{max_len}bp)")
        run_seqkit_length_filter(
            input_file=output_file,
            output_file=filtered_output_file,
            min_length=min_len,
            max_length=max_len,
            threads=cores,
        )
        logger.success(f"Length filtered: {filtered_output_file}")


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
    adapter: str = typer.Option(
        DEFAULT_ADAPTER,
        "--adapter",
        "-a",
        help="3' adapter sequence to remove",
    ),
    length_filter: str = typer.Option(
        None,
        "--length-filter",
        "-l",
        help="Optional length filter range (e.g., '24-35'). If set, generates additional filtered file.",
    ),
    quality: int = typer.Option(
        20,
        "--quality",
        "-q",
        help="Quality threshold for fastp filtering",
    ),
    overlap: int = typer.Option(
        5,
        "--overlap",
        help="Minimum overlap length for adapter detection",
    ),
    discard_untrimmed: bool = typer.Option(
        True,
        "--discard-untrimmed/--keep-untrimmed",
        help="Discard reads without adapter",
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
        REPORTS_DIR / "fastp",
        "--report-dir",
        "-r",
        help="Directory for fastp QC reports",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Force re-run even if output files exist",
    ),
) -> None:
    """Quality control for smRNA-seq data using cutadapt and fastp.

    Process:
    1. cutadapt: adapter trimming
    2. fastp: quality filtering (Q20) and QC report generation
    3. seqkit (optional): length filtering if --length-filter is specified
    4. multiqc: aggregate all fastp reports

    Output files:
    - {sample}.trim.fq.gz: QC'd reads (always generated)
    - {min}-{max}_length_filtered/{sample}.fq.gz: length-filtered reads (only if --length-filter is set)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

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
                f"Length filter enabled: {length_filter_tuple[0]}-{length_filter_tuple[1]}bp"
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

    logger.info(f"Processing {len(fastq_files)} file(s) with adapter: {adapter}")

    processed_count = 0
    skipped_count = 0

    for fastq_file in fastq_files:
        sample_name = fastq_file.stem.split(".")[0]
        output_file = output_dir / f"{sample_name}.trim.fq.gz"
        filtered_file = filtered_dir / f"{sample_name}.fq.gz" if filtered_dir else None
        html_report = report_dir / f"{sample_name}.html"
        json_report = report_dir / f"{sample_name}.json"

        # Check if outputs already exist
        outputs_exist = output_file.exists() and html_report.exists() and json_report.exists()
        if length_filter_tuple and filtered_file:
            outputs_exist = outputs_exist and filtered_file.exists()

        if not force and outputs_exist:
            logger.info(f"⏭️ Skipping {sample_name} (outputs exist)")
            skipped_count += 1
            continue

        try:
            qc_single_file(
                input_file=fastq_file,
                output_file=output_file,
                report_dir=report_dir,
                adapter=adapter,
                quality=quality,
                overlap=overlap,
                discard_untrimmed=discard_untrimmed,
                cores=cores,
                length_filter=length_filter_tuple,
                filtered_output_file=filtered_file,
            )
            processed_count += 1
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to process {fastq_file.name}: {e.stderr}")
            raise typer.Exit(1)

    logger.info(f"Processed: {processed_count}, Skipped: {skipped_count}")

    # Step 3: Generate fastp summary statistics
    json_files = sorted(report_dir.glob("*.json"))
    if json_files:
        summary_file = report_dir / "fastp_summary.tsv"
        logger.info("[summary] Generating fastp summary statistics")
        write_fastp_summary(json_files, summary_file)
        logger.success(f"Summary statistics: {summary_file}")

    # Step 4: Aggregate reports with MultiQC
    multiqc_report = report_dir / "multiqc_report.html"
    if not force and multiqc_report.exists() and processed_count == 0:
        logger.info("⏭️ Skipping MultiQC (report exists, no new samples)")
    else:
        logger.info("[multiqc] Aggregating QC reports")
        try:
            run_multiqc(
                input_dir=report_dir,
                output_dir=report_dir,
            )
            logger.success(f"MultiQC report: {multiqc_report}")
        except subprocess.CalledProcessError as e:
            logger.warning(f"MultiQC failed (non-critical): {e.stderr}")

    logger.success(f"All QC complete. Output: {output_dir}")


if __name__ == "__main__":
    app()
