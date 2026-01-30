"""Match reads against known piRNA database and generate statistics."""

from collections import Counter
from dataclasses import dataclass, field
from multiprocessing import Pool
from pathlib import Path

from Bio import SeqIO
from loguru import logger
import typer

from bm_pirna.config import EXTERNAL_DATA_DIR, PROCESSED_DATA_DIR

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
    # Sequence to piRNA ID mapping (sequence -> piRNA_id)
    seq_to_pirna_id: dict = field(default_factory=dict)


NUCLEOTIDES = "ATGC"


def generate_1mm_variants(seq: str) -> set[str]:
    """Generate all 1-mismatch variants of a sequence.

    For a sequence of length L, generates 3*L variants (each position can be
    replaced by 3 other nucleotides).
    """
    variants = set()
    seq_list = list(seq)
    for i in range(len(seq)):
        orig = seq_list[i]
        for nt in NUCLEOTIDES:
            if nt != orig:
                seq_list[i] = nt
                variants.add("".join(seq_list))
        seq_list[i] = orig
    return variants


def load_pirna_database(
    db_file: Path,
    build_variants: bool = False,
) -> tuple[set[str], dict[str, str], set[str]]:
    """Load piRNA sequences from database file.

    Args:
        db_file: Path to piRNA database file (fasta format)
        build_variants: If True, build 1-mismatch variant lookup set

    Returns:
        Tuple of (set of all sequences, dict of sequence to piRNA ID mapping,
                  set of all 1-mismatch variants mapped to original sequence)
    """
    pirna_seqs = set()
    seq_to_id: dict[str, str] = {}

    for record in SeqIO.parse(db_file, "fasta"):
        seq = str(record.seq).upper()
        pirna_seqs.add(seq)
        seq_to_id[seq] = record.id

    logger.info(f"Loaded {len(pirna_seqs)} unique piRNA sequences from database")

    # Build variant lookup: variant -> original sequence
    variant_to_orig: dict[str, str] = {}
    if build_variants:
        logger.info("Building 1-mismatch variant index (this may take a moment)...")
        for seq in pirna_seqs:
            for variant in generate_1mm_variants(seq):
                # Only store if not already a piRNA itself
                if variant not in pirna_seqs:
                    variant_to_orig[variant] = seq
        logger.info(f"Built {len(variant_to_orig):,} variant mappings")

    return pirna_seqs, seq_to_id, variant_to_orig


def find_match(
    query: str,
    db_seqs: set[str],
    variant_to_orig: dict[str, str],
    allow_mismatch: bool = False,
) -> str | None:
    """Find matching sequence in database using O(1) set lookup.

    Args:
        query: Query sequence
        db_seqs: Set of database sequences
        variant_to_orig: Dict mapping 1-mismatch variants to original sequences
        allow_mismatch: If True, also check 1-mismatch variants

    Returns:
        Matched database sequence if found, None otherwise
    """
    # Exact match - O(1)
    if query in db_seqs:
        return query

    # 1-mismatch match via precomputed variants - O(1)
    if allow_mismatch and query in variant_to_orig:
        return variant_to_orig[query]

    return None


def parse_collapsed_fasta(fasta_file: Path) -> list[tuple[str, int]]:
    """Parse collapsed fasta file and return (sequence, count) pairs.

    Args:
        fasta_file: Path to collapsed fasta file with format >rank-count

    Returns:
        List of (sequence, count) tuples
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Header format: >rank-count (e.g., >1-29509)
        parts = record.id.split("-")
        if len(parts) >= 2:
            count = int(parts[-1])
        else:
            count = 1
        sequences.append((str(record.seq).upper(), count))
    return sequences


# Global variables for multiprocessing (set by main before spawning workers)
_PIRNA_DB: set[str] = set()
_VARIANT_TO_ORIG: dict[str, str] = {}
_SEQ_TO_ID: dict[str, str] = {}
_ALLOW_MISMATCH: bool = False


def _init_worker(
    pirna_db: set[str],
    variant_to_orig: dict[str, str],
    seq_to_id: dict[str, str],
    allow_mismatch: bool,
) -> None:
    """Initialize worker process with shared data."""
    global _PIRNA_DB, _VARIANT_TO_ORIG, _SEQ_TO_ID, _ALLOW_MISMATCH
    _PIRNA_DB = pirna_db
    _VARIANT_TO_ORIG = variant_to_orig
    _SEQ_TO_ID = seq_to_id
    _ALLOW_MISMATCH = allow_mismatch


def _process_file_worker(collapsed_file: Path) -> SampleStats:
    """Worker function for multiprocessing."""
    return filter_and_analyze(
        collapsed_file, _PIRNA_DB, _VARIANT_TO_ORIG, _SEQ_TO_ID, _ALLOW_MISMATCH
    )


def filter_and_analyze(
    collapsed_file: Path,
    pirna_db: set[str],
    variant_to_orig: dict[str, str],
    seq_to_id: dict[str, str],
    allow_mismatch: bool = False,
) -> SampleStats:
    """Filter reads that match piRNA database and generate statistics.

    Args:
        collapsed_file: Path to collapsed fasta file
        pirna_db: Set of piRNA sequences
        variant_to_orig: Dict mapping 1-mismatch variants to original sequences
        seq_to_id: Dict mapping sequences to piRNA IDs
        allow_mismatch: If True, allow 1 mismatch

    Returns:
        SampleStats object with statistics of matched reads
    """
    sample_name = collapsed_file.stem
    stats = SampleStats(sample_name=sample_name)

    lengths = []
    gc_count = 0
    total_bases = 0

    sequences = parse_collapsed_fasta(collapsed_file)
    total_input_reads = sum(count for _, count in sequences)

    for seq, count in sequences:
        matched_seq = find_match(seq, pirna_db, variant_to_orig, allow_mismatch)
        if matched_seq is None:
            continue

        # Store the piRNA ID for this sequence
        pirna_id = seq_to_id.get(matched_seq, "unknown")
        stats.seq_to_pirna_id[seq] = pirna_id

        seq_len = len(seq)

        # Each unique sequence is counted by its read count
        stats.total_reads += count
        total_bases += seq_len * count
        lengths.extend([seq_len] * count)

        # Length distribution
        stats.length_counts[seq_len] += count

        # Sequence counts for collapse
        stats.seq_counts[seq] += count

        # First nucleotide (1U bias)
        if seq:
            stats.first_nt_counts[seq[0]] += count

        # 10th position nucleotide (10A bias for ping-pong)
        if seq_len >= 10:
            stats.tenth_nt_counts[seq[9]] += count

        # Last nucleotide
        if seq:
            stats.last_nt_counts[seq[-1]] += count

        # Nucleotide composition
        for nt in seq:
            stats.nt_composition[nt] += count
            if nt in ("G", "C"):
                gc_count += count

    # Calculate summary statistics
    stats.total_bases = total_bases
    if lengths:
        stats.min_length = min(lengths)
        stats.max_length = max(lengths)
        stats.mean_length = sum(lengths) / len(lengths)
    if total_bases > 0:
        stats.gc_content = gc_count / total_bases * 100

    logger.info(
        f"  Matched: {stats.total_reads:,} / {total_input_reads:,} reads "
        f"({stats.total_reads / total_input_reads * 100:.1f}%)"
    )

    return stats


def write_summary_stats(stats_list: list[SampleStats], output_file: Path) -> None:
    """Write summary statistics to TSV file."""
    with open(output_file, "w") as f:
        f.write(
            "sample\ttotal_reads\ttotal_bases\tmin_length\tmax_length\tmean_length\t"
            "gc_content\tfirst_U_pct\ttenth_A_pct\n"
        )
        for s in stats_list:
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
                f.write(f"{s.sample_name}\t{nt}\t{count}\t{pct:.2f}\n")


def write_tenth_nt_distribution(stats_list: list[SampleStats], output_file: Path) -> None:
    """Write 10th position nucleotide distribution to TSV file."""
    with open(output_file, "w") as f:
        f.write("sample\tnucleotide\tcount\tpercentage\n")
        for s in stats_list:
            total = sum(s.tenth_nt_counts.values())
            for nt in ["A", "T", "G", "C"]:
                count = s.tenth_nt_counts.get(nt, 0)
                pct = count / total * 100 if total > 0 else 0
                f.write(f"{s.sample_name}\t{nt}\t{count}\t{pct:.2f}\n")


def write_nt_composition(stats_list: list[SampleStats], output_file: Path) -> None:
    """Write overall nucleotide composition to TSV file."""
    with open(output_file, "w") as f:
        f.write("sample\tnucleotide\tcount\tpercentage\n")
        for s in stats_list:
            total = sum(s.nt_composition.values())
            for nt in ["A", "T", "G", "C"]:
                count = s.nt_composition.get(nt, 0)
                pct = count / total * 100 if total > 0 else 0
                f.write(f"{s.sample_name}\t{nt}\t{count}\t{pct:.2f}\n")


def write_collapsed_fasta(stats: SampleStats, output_file: Path) -> None:
    """Write collapsed fasta file with format >piRNA_ID_count."""
    sorted_seqs = stats.seq_counts.most_common()
    with open(output_file, "w") as f:
        for seq, count in sorted_seqs:
            pirna_id = stats.seq_to_pirna_id.get(seq, "unknown")
            f.write(f">{pirna_id}_{count}\n{seq}\n")


@app.command()
def main(
    input_path: Path = typer.Argument(
        ...,
        help="Input collapsed fasta file or directory containing .fa files",
    ),
    pirna_db: Path = typer.Option(
        EXTERNAL_DATA_DIR / "bmo_piRNA.v3.0.txt",
        "--db",
        "-d",
        help="Path to piRNA database file (fasta format)",
    ),
    output_dir: Path = typer.Option(
        PROCESSED_DATA_DIR / "pirna_matched_stats",
        "--output-dir",
        "-o",
        help="Output directory for statistics files",
    ),
    max_mismatch: int = typer.Option(
        0,
        "--mismatch",
        "-m",
        help="Maximum allowed mismatches (0 for exact match, 1 for 1 mismatch)",
    ),
    pattern: str = typer.Option(
        "*.fa",
        "--pattern",
        "-p",
        help="Glob pattern to match fasta files (when input is a directory)",
    ),
    threads: int = typer.Option(
        4,
        "--threads",
        "-t",
        help="Number of parallel processes",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Force re-run even if output files exist",
    ),
) -> None:
    """Match collapsed reads against piRNA database and generate statistics.

    Filters reads that match known piRNA sequences (with optional mismatches)
    and outputs statistics in the same format as stats.py.

    Examples:
        # Exact match
        python -m bm_pirna.smrna_seq.pirna_match collapsed/

        # Allow 1 mismatch
        python -m bm_pirna.smrna_seq.pirna_match collapsed/ -m 1

        # Specify custom database
        python -m bm_pirna.smrna_seq.pirna_match collapsed/ -d my_pirna.fa
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if outputs already exist
    summary_file = output_dir / "summary_stats.tsv"
    if not force and summary_file.exists():
        logger.info(f"⏭️ Skipping (outputs exist): {output_dir}")
        logger.info("Use --force to overwrite")
        return

    # Load piRNA database
    if not pirna_db.exists():
        logger.error(f"piRNA database not found: {pirna_db}")
        raise typer.Exit(1)
    allow_mismatch = max_mismatch > 0
    pirna_seqs, seq_to_id, variant_to_orig = load_pirna_database(
        pirna_db, build_variants=allow_mismatch
    )

    # Get input files
    if input_path.is_file():
        fasta_files = [input_path]
    elif input_path.is_dir():
        fasta_files = sorted(input_path.glob(pattern))
        if not fasta_files:
            logger.error(f"No fasta files found in {input_path} with pattern '{pattern}'")
            raise typer.Exit(1)
    else:
        logger.error(f"Input path does not exist: {input_path}")
        raise typer.Exit(1)

    logger.info(f"Processing {len(fasta_files)} file(s) with {threads} threads")
    logger.info(f"Max mismatches allowed: {max_mismatch}")

    # Analyze all files in parallel
    if threads > 1 and len(fasta_files) > 1:
        # Initialize global variables for workers
        global _PIRNA_DB, _VARIANT_TO_ORIG, _SEQ_TO_ID, _ALLOW_MISMATCH
        _PIRNA_DB = pirna_seqs
        _VARIANT_TO_ORIG = variant_to_orig
        _SEQ_TO_ID = seq_to_id
        _ALLOW_MISMATCH = allow_mismatch

        with Pool(
            processes=min(threads, len(fasta_files)),
            initializer=_init_worker,
            initargs=(pirna_seqs, variant_to_orig, seq_to_id, allow_mismatch),
        ) as pool:
            stats_list = pool.map(_process_file_worker, fasta_files)
    else:
        stats_list = []
        for i, fasta_file in enumerate(fasta_files, 1):
            logger.info(f"[{i}/{len(fasta_files)}] Processing {fasta_file.name}")
            stats = filter_and_analyze(
                fasta_file, pirna_seqs, variant_to_orig, seq_to_id, allow_mismatch
            )
            stats_list.append(stats)

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
    logger.info("  - collapsed/: Matched reads in collapsed fasta format")


if __name__ == "__main__":
    app()
