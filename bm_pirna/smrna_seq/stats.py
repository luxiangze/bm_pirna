"""Generate comprehensive statistics for fastq files for piRNA analysis."""

from collections import Counter
from dataclasses import dataclass, field
import gzip
from pathlib import Path

from Bio import SeqIO
from loguru import logger
import typer

from bm_pirna.config import PROCESSED_DATA_DIR

app = typer.Typer()


@dataclass
class SampleStats:
    """Statistics for a single sample."""

    sample_name: str
    total_reads: int = 0
    total_bases: int = 0
    min_length: int = 0
    max_length: int = 0
    mean_length: float = 0.0
    gc_content: float = 0.0
    # Length distribution
    length_counts: Counter = field(default_factory=Counter)
    # First nucleotide distribution (1U bias for piRNA)
    first_nt_counts: Counter = field(default_factory=Counter)
    # 10th position nucleotide (10A bias for ping-pong)
    tenth_nt_counts: Counter = field(default_factory=Counter)
    # Last nucleotide distribution
    last_nt_counts: Counter = field(default_factory=Counter)
    # Nucleotide composition
    nt_composition: Counter = field(default_factory=Counter)
    # Sequence counts for collapse (sequence -> count)
    seq_counts: Counter = field(default_factory=Counter)


def analyze_fastq(input_file: Path) -> SampleStats:
    """Analyze a single fastq file and return statistics.

    Args:
        input_file: Input fastq file (can be gzipped)

    Returns:
        SampleStats object with all statistics
    """
    sample_name = input_file.stem.split(".")[0]
    stats = SampleStats(sample_name=sample_name)

    lengths = []
    gc_count = 0
    total_bases = 0

    # Open file (handle gzip)
    if str(input_file).endswith(".gz"):
        handle = gzip.open(input_file, "rt")
    else:
        handle = open(input_file)

    try:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq).upper()
            seq_len = len(seq)

            stats.total_reads += 1
            total_bases += seq_len
            lengths.append(seq_len)

            # Length distribution
            stats.length_counts[seq_len] += 1

            # Sequence counts for collapse
            stats.seq_counts[seq] += 1

            # First nucleotide (1U bias)
            if seq:
                stats.first_nt_counts[seq[0]] += 1

            # 10th position nucleotide (10A bias for ping-pong)
            if seq_len >= 10:
                stats.tenth_nt_counts[seq[9]] += 1

            # Last nucleotide
            if seq:
                stats.last_nt_counts[seq[-1]] += 1

            # Nucleotide composition
            for nt in seq:
                stats.nt_composition[nt] += 1
                if nt in ("G", "C"):
                    gc_count += 1
    finally:
        handle.close()

    # Calculate summary statistics
    stats.total_bases = total_bases
    if lengths:
        stats.min_length = min(lengths)
        stats.max_length = max(lengths)
        stats.mean_length = sum(lengths) / len(lengths)
    if total_bases > 0:
        stats.gc_content = gc_count / total_bases * 100

    return stats


def write_summary_stats(stats_list: list[SampleStats], output_file: Path) -> None:
    """Write summary statistics to TSV file."""
    with open(output_file, "w") as f:
        # Header
        f.write(
            "sample\ttotal_reads\ttotal_bases\tmin_length\tmax_length\tmean_length\t"
            "gc_content\tfirst_U_pct\ttenth_A_pct\n"
        )
        for s in stats_list:
            # Calculate 1U and 10A percentages
            first_u_pct = (
                s.first_nt_counts.get("T", 0) / s.total_reads * 100 if s.total_reads > 0 else 0
            )
            tenth_a_pct = (
                s.tenth_nt_counts.get("A", 0) / sum(s.tenth_nt_counts.values()) * 100
                if sum(s.tenth_nt_counts.values()) > 0
                else 0
            )
            f.write(
                f"{s.sample_name}\t{s.total_reads}\t{s.total_bases}\t"
                f"{s.min_length}\t{s.max_length}\t{s.mean_length:.2f}\t"
                f"{s.gc_content:.2f}\t{first_u_pct:.2f}\t{tenth_a_pct:.2f}\n"
            )


def write_length_distribution(stats_list: list[SampleStats], output_file: Path) -> None:
    """Write length distribution to TSV file (long format for ggplot2)."""
    with open(output_file, "w") as f:
        f.write("sample\tlength\tcount\n")
        for s in stats_list:
            for length, count in sorted(s.length_counts.items()):
                f.write(f"{s.sample_name}\t{length}\t{count}\n")


def write_first_nt_distribution(stats_list: list[SampleStats], output_file: Path) -> None:
    """Write first nucleotide distribution to TSV file."""
    with open(output_file, "w") as f:
        f.write("sample\tnucleotide\tcount\tpercentage\n")
        for s in stats_list:
            total = sum(s.first_nt_counts.values())
            for nt in ["A", "T", "G", "C"]:
                count = s.first_nt_counts.get(nt, 0)
                pct = count / total * 100 if total > 0 else 0
                f.write(f"{s.sample_name}\t{nt}\t{count}\t{pct:.4f}\n")


def write_tenth_nt_distribution(stats_list: list[SampleStats], output_file: Path) -> None:
    """Write 10th position nucleotide distribution to TSV file."""
    with open(output_file, "w") as f:
        f.write("sample\tnucleotide\tcount\tpercentage\n")
        for s in stats_list:
            total = sum(s.tenth_nt_counts.values())
            for nt in ["A", "T", "G", "C"]:
                count = s.tenth_nt_counts.get(nt, 0)
                pct = count / total * 100 if total > 0 else 0
                f.write(f"{s.sample_name}\t{nt}\t{count}\t{pct:.4f}\n")


def write_nt_composition(stats_list: list[SampleStats], output_file: Path) -> None:
    """Write overall nucleotide composition to TSV file."""
    with open(output_file, "w") as f:
        f.write("sample\tnucleotide\tcount\tpercentage\n")
        for s in stats_list:
            total = sum(s.nt_composition.values())
            for nt in ["A", "T", "G", "C"]:
                count = s.nt_composition.get(nt, 0)
                pct = count / total * 100 if total > 0 else 0
                f.write(f"{s.sample_name}\t{nt}\t{count}\t{pct:.4f}\n")


def write_collapsed_fasta(stats: SampleStats, output_file: Path) -> None:
    """Write collapsed fasta file with format >rank-count.

    Output format:
    >1-1000
    TGCGATCG
    >2-500
    GCGGATCG

    Where the first number is the rank (1 = most abundant) and
    the second number is the read count.
    """
    # Sort sequences by count (descending)
    sorted_seqs = stats.seq_counts.most_common()

    with open(output_file, "w") as f:
        for rank, (seq, count) in enumerate(sorted_seqs, 1):
            f.write(f">{rank}-{count}\n{seq}\n")


@app.command()
def main(
    input_path: Path = typer.Argument(
        ...,
        help="Input fastq file or directory containing fastq files",
    ),
    output_dir: Path = typer.Option(
        PROCESSED_DATA_DIR / "stats",
        "--output-dir",
        "-o",
        help="Output directory for statistics files",
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
    """Generate comprehensive statistics for fastq files for piRNA analysis.

    Outputs multiple TSV files for R visualization:
    - summary_stats.tsv: Basic statistics per sample
    - length_distribution.tsv: Read length distribution
    - first_nt_distribution.tsv: First nucleotide distribution (1U bias)
    - tenth_nt_distribution.tsv: 10th position nucleotide (10A ping-pong)
    - nt_composition.tsv: Overall nucleotide composition
    - collapsed/: Directory containing collapsed fasta files per sample

    Examples:
        # Stats for all fastq files in a directory
        python -m bm_pirna.smrna_seq.stats data/interim/qc/

        # Specify output directory
        python -m bm_pirna.smrna_seq.stats data/interim/qc/ -o reports/pirna_stats/
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if outputs already exist
    summary_file = output_dir / "summary_stats.tsv"
    if not force and summary_file.exists():
        logger.info(f"⏭️ Skipping (outputs exist): {output_dir}")
        logger.info("Use --force to overwrite")
        return

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

    logger.info(f"Analyzing {len(fastq_files)} file(s)")

    # Analyze all files
    stats_list = []
    for i, fastq_file in enumerate(fastq_files, 1):
        logger.info(f"[{i}/{len(fastq_files)}] Analyzing {fastq_file.name}")
        stats = analyze_fastq(fastq_file)
        stats_list.append(stats)
        logger.debug(
            f"  Reads: {stats.total_reads:,}, "
            f"Length: {stats.min_length}-{stats.max_length}, "
            f"GC: {stats.gc_content:.1f}%"
        )

    # Write output files
    logger.info("Writing output files...")

    write_summary_stats(stats_list, output_dir / "summary_stats.tsv")
    write_length_distribution(stats_list, output_dir / "length_distribution.tsv")
    write_first_nt_distribution(stats_list, output_dir / "first_nt_distribution.tsv")
    write_tenth_nt_distribution(stats_list, output_dir / "tenth_nt_distribution.tsv")
    write_nt_composition(stats_list, output_dir / "nt_composition.tsv")

    # Write collapsed fasta files
    collapsed_dir = output_dir / "collapsed"
    collapsed_dir.mkdir(parents=True, exist_ok=True)
    for stats in stats_list:
        collapsed_file = collapsed_dir / f"{stats.sample_name}.fa"
        write_collapsed_fasta(stats, collapsed_file)

    logger.success(f"All stats saved to: {output_dir}")
    logger.info("Output files:")
    logger.info("  - summary_stats.tsv: Basic statistics")
    logger.info("  - length_distribution.tsv: Read length distribution")
    logger.info("  - first_nt_distribution.tsv: 1U bias analysis")
    logger.info("  - tenth_nt_distribution.tsv: 10A ping-pong analysis")
    logger.info("  - nt_composition.tsv: Nucleotide composition")
    logger.info("  - collapsed/: Collapsed fasta files (>rank-count format)")


if __name__ == "__main__":
    app()
