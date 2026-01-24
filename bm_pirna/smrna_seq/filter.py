"""Filter reads using bowtie alignment."""

from enum import Enum
from pathlib import Path
import subprocess
import tempfile

from loguru import logger
import typer

from bm_pirna.config import INTERIM_DATA_DIR

app = typer.Typer()


class FilterMode(str, Enum):
    """Filter mode for bowtie alignment."""

    EXCLUDE = "exclude"  # Remove reads that align (e.g., filter out miRNA/rRNA)
    INCLUDE = "include"  # Keep only reads that align (e.g., keep transposon reads)


def run_bowtie(
    input_file: Path,
    index: Path,
    aligned_file: Path | None = None,
    unaligned_file: Path | None = None,
    threads: int = 16,
    mismatches: int = 0,
    best: bool = True,
) -> subprocess.CompletedProcess:
    """Run bowtie alignment using -v mode.

    Args:
        input_file: Input fastq file
        index: Bowtie index prefix (without .ebwt extension)
        aligned_file: Output file for aligned reads (optional)
        unaligned_file: Output file for unaligned reads (optional)
        threads: Number of threads
        mismatches: Number of mismatches allowed in entire read (-v mode)
        best: Report best alignments (--best)
    """
    cmd = [
        "bowtie",
        "-p",
        str(threads),
        "-v",
        str(mismatches),
        "-q",  # input is fastq
        "-S",  # output SAM format (for stats)
    ]

    if best:
        cmd.append("--best")

    if unaligned_file:
        cmd.extend(["--un", str(unaligned_file)])

    if aligned_file:
        cmd.extend(["--al", str(aligned_file)])

    cmd.extend([str(index), str(input_file)])

    logger.debug(f"Running: {' '.join(cmd)}")
    return subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )


def compress_with_pigz(
    input_file: Path,
    output_file: Path,
) -> subprocess.CompletedProcess:
    """Compress a file using pigz."""
    cmd = [
        "pigz",
        "-c",
        str(input_file),
    ]

    logger.debug(f"Running: {' '.join(cmd)} > {output_file}")
    with open(output_file, "wb") as f:
        return subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE)


def filter_single_file(
    input_file: Path,
    output_file: Path,
    index: Path,
    mode: FilterMode,
    threads: int = 16,
    mismatches: int = 0,
) -> None:
    """Filter a single fastq file using bowtie.

    Args:
        input_file: Input fastq file
        output_file: Output fastq file
        index: Bowtie index prefix
        mode: Filter mode (exclude or include)
        threads: Number of threads
        mismatches: Number of mismatches allowed (-v mode)
    """
    logger.info(f"[bowtie] Filtering {input_file.name} ({mode.value} mode)")

    # Check if output needs gzip compression
    needs_gzip = str(output_file).endswith(".gz")

    if needs_gzip:
        # Use temp file for bowtie output, then compress with pigz
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_output = Path(tmpdir) / "bowtie_output.fq"

            if mode == FilterMode.EXCLUDE:
                run_bowtie(
                    input_file=input_file,
                    index=index,
                    unaligned_file=temp_output,
                    threads=threads,
                    mismatches=mismatches,
                )
            else:
                run_bowtie(
                    input_file=input_file,
                    index=index,
                    aligned_file=temp_output,
                    threads=threads,
                    mismatches=mismatches,
                )

            # Compress with pigz
            logger.info(f"[pigz] Compressing to {output_file.name}")
            compress_with_pigz(
                input_file=temp_output,
                output_file=output_file,
            )
    else:
        # No compression needed
        if mode == FilterMode.EXCLUDE:
            run_bowtie(
                input_file=input_file,
                index=index,
                unaligned_file=output_file,
                threads=threads,
                mismatches=mismatches,
            )
        else:
            run_bowtie(
                input_file=input_file,
                index=index,
                aligned_file=output_file,
                threads=threads,
                mismatches=mismatches,
            )

    logger.success(f"Filter complete: {output_file}")


@app.command()
def main(
    input_path: Path = typer.Argument(
        ...,
        help="Input fastq file or directory containing fastq files",
    ),
    index: Path = typer.Argument(
        ...,
        help="Bowtie index prefix (e.g., data/external/miRNA/miRNA_index)",
    ),
    output_dir: Path = typer.Option(
        INTERIM_DATA_DIR / "filtered",
        "--output-dir",
        "-o",
        help="Output directory for filtered files",
    ),
    mode: FilterMode = typer.Option(
        FilterMode.EXCLUDE,
        "--mode",
        "-m",
        help="Filter mode: exclude (remove aligned) or include (keep aligned)",
    ),
    suffix: str = typer.Option(
        ".filtered.fq.gz",
        "--suffix",
        "-s",
        help="Suffix for output files",
    ),
    threads: int = typer.Option(
        16,
        "--threads",
        "-t",
        help="Number of threads for bowtie",
    ),
    mismatches: int = typer.Option(
        0,
        "--mismatches",
        "-v",
        help="Number of mismatches allowed in entire read (-v mode)",
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
    """Filter reads using bowtie alignment.

    Two filter modes:
    - exclude: Remove reads that align to the index (e.g., filter out miRNA/rRNA)
    - include: Keep only reads that align to the index (e.g., keep transposon reads)

    Examples:
        # Filter out miRNA reads
        python -m bm_pirna.smrna_seq.filter input.fq.gz miRNA_index -m exclude

        # Keep only transposon-mapping reads
        python -m bm_pirna.smrna_seq.filter input.fq.gz TE_index -m include
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

    logger.info(f"Filtering {len(fastq_files)} file(s) using index: {index.name}")
    logger.info(f"Mode: {mode.value}")

    processed_count = 0
    skipped_count = 0

    for fastq_file in fastq_files:
        sample_name = fastq_file.stem.split(".")[0]
        output_file = output_dir / f"{sample_name}{suffix}"

        # Check if output already exists
        if not force and output_file.exists():
            logger.info(f"⏭️ Skipping {sample_name} (output exists)")
            skipped_count += 1
            continue

        try:
            filter_single_file(
                input_file=fastq_file,
                output_file=output_file,
                index=index,
                mode=mode,
                threads=threads,
                mismatches=mismatches,
            )
            processed_count += 1
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to process {fastq_file.name}: {e.stderr}")
            raise typer.Exit(1)

    logger.info(f"Processed: {processed_count}, Skipped: {skipped_count}")
    logger.success(f"All filtering complete. Output: {output_dir}")


if __name__ == "__main__":
    app()
