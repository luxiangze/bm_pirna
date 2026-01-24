"""Align reads to genome using bowtie."""

from enum import Enum
from pathlib import Path
import subprocess

from loguru import logger
import typer

from bm_pirna.config import PROCESSED_DATA_DIR

app = typer.Typer()


class OutputFormat(str, Enum):
    """Output format for alignment."""

    SAM = "sam"
    BAM = "bam"


def run_bowtie_align(
    input_file: Path,
    index: Path,
    output_file: Path,
    threads: int = 16,
    mismatches: int = 0,
    report_all: bool = False,
    best: bool = True,
) -> subprocess.CompletedProcess:
    """Run bowtie alignment to genome.

    Args:
        input_file: Input fastq file
        index: Bowtie index prefix (without .ebwt extension)
        output_file: Output SAM file
        threads: Number of threads
        mismatches: Number of mismatches allowed (-v mode)
        report_all: Report all alignments (-a)
        best: Report best alignments (--best)
    """
    cmd = [
        "bowtie",
        "-p",
        str(threads),
        "-v",
        str(mismatches),
        "-q",  # input is fastq
        "-S",  # output SAM format
    ]

    if report_all:
        cmd.append("-a")

    if best:
        cmd.append("--best")

    cmd.extend([str(index), str(input_file), str(output_file)])

    logger.debug(f"Running: {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def sam_to_bam(
    sam_file: Path,
    bam_file: Path,
    threads: int = 16,
    sort: bool = True,
    index: bool = True,
) -> None:
    """Convert SAM to sorted BAM and create index.

    Args:
        sam_file: Input SAM file
        bam_file: Output BAM file
        threads: Number of threads
        sort: Sort the BAM file
        index: Create BAM index (.bai)
    """
    if sort:
        # samtools sort directly from SAM to sorted BAM
        cmd = [
            "samtools",
            "sort",
            "-@",
            str(threads),
            "-o",
            str(bam_file),
            str(sam_file),
        ]
        logger.debug(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    else:
        # samtools view to convert SAM to BAM
        cmd = [
            "samtools",
            "view",
            "-@",
            str(threads),
            "-bS",
            "-o",
            str(bam_file),
            str(sam_file),
        ]
        logger.debug(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)

    if index:
        # Create BAM index
        cmd = ["samtools", "index", "-@", str(threads), str(bam_file)]
        logger.debug(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)


def align_single_file(
    input_file: Path,
    output_file: Path,
    index: Path,
    threads: int = 16,
    mismatches: int = 0,
    report_all: bool = False,
    output_format: OutputFormat = OutputFormat.BAM,
) -> None:
    """Align a single fastq file to genome.

    Args:
        input_file: Input fastq file
        output_file: Output SAM/BAM file
        index: Bowtie index prefix
        threads: Number of threads
        mismatches: Number of mismatches allowed (-v mode)
        report_all: Report all alignments
        output_format: Output format (sam or bam)
    """
    logger.info(f"[bowtie] Aligning {input_file.name} to {index.name}")

    if output_format == OutputFormat.SAM:
        # Direct SAM output
        run_bowtie_align(
            input_file=input_file,
            index=index,
            output_file=output_file,
            threads=threads,
            mismatches=mismatches,
            report_all=report_all,
        )
    else:
        # Output BAM: align to temp SAM, then convert
        sam_file = output_file.with_suffix(".sam")
        run_bowtie_align(
            input_file=input_file,
            index=index,
            output_file=sam_file,
            threads=threads,
            mismatches=mismatches,
            report_all=report_all,
        )

        logger.info("[samtools] Converting to sorted BAM")
        sam_to_bam(
            sam_file=sam_file,
            bam_file=output_file,
            threads=threads,
            sort=True,
            index=True,
        )

        # Remove temp SAM file
        sam_file.unlink()

    logger.success(f"Alignment complete: {output_file}")


@app.command()
def main(
    input_path: Path = typer.Argument(
        ...,
        help="Input fastq file or directory containing fastq files",
    ),
    index: Path = typer.Argument(
        ...,
        help="Bowtie index prefix (e.g., data/external/genome/genome_index)",
    ),
    output_dir: Path = typer.Option(
        PROCESSED_DATA_DIR / "aligned",
        "--output-dir",
        "-o",
        help="Output directory for aligned files",
    ),
    output_format: OutputFormat = typer.Option(
        OutputFormat.BAM,
        "--format",
        "-F",
        help="Output format: sam or bam",
    ),
    threads: int = typer.Option(
        16,
        "--threads",
        "-t",
        help="Number of threads",
    ),
    mismatches: int = typer.Option(
        0,
        "--mismatches",
        "-v",
        help="Number of mismatches allowed in entire read (-v mode)",
    ),
    report_all: bool = typer.Option(
        False,
        "--all",
        "-a",
        help="Report all alignments (default: report best)",
    ),
    pattern: str = typer.Option(
        "*.f*q.gz",
        "--pattern",
        "-p",
        help="Glob pattern to match fastq files (when input is a directory)",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Force re-run even if output files exist",
    ),
) -> None:
    """Align reads to genome using bowtie.

    Aligns filtered reads to a reference genome and outputs SAM or sorted BAM files.

    Examples:
        # Align to genome, output BAM
        python -m bm_pirna.smrna_seq.align data/interim/filtered/ genome_index

        # Align with 1 mismatch, output SAM
        python -m bm_pirna.smrna_seq.align data/interim/filtered/ genome_index -v 1 -F sam
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate index exists
    index_file = Path(f"{index}.1.ebwt")
    if not index_file.exists():
        logger.error(f"Bowtie index not found: {index}")
        logger.info("Build index with: bowtie-build reference.fa index_prefix")
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

    logger.info(f"Aligning {len(fastq_files)} file(s) to {index.name}")
    logger.info(f"Output format: {output_format.value}")

    processed_count = 0
    skipped_count = 0

    for fastq_file in fastq_files:
        sample_name = fastq_file.stem.split(".")[0]
        suffix = ".bam" if output_format == OutputFormat.BAM else ".sam"
        output_file = output_dir / f"{sample_name}{suffix}"

        # Check if output already exists
        if not force and output_file.exists():
            logger.info(f"⏭️ Skipping {sample_name} (output exists)")
            skipped_count += 1
            continue

        try:
            align_single_file(
                input_file=fastq_file,
                output_file=output_file,
                index=index,
                threads=threads,
                mismatches=mismatches,
                report_all=report_all,
                output_format=output_format,
            )
            processed_count += 1
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to process {fastq_file.name}: {e.stderr}")
            raise typer.Exit(1)

    logger.info(f"Processed: {processed_count}, Skipped: {skipped_count}")
    logger.success(f"All alignments complete. Output: {output_dir}")


if __name__ == "__main__":
    app()
