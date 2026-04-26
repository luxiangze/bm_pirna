"""Build a Drosophila TE accession -> class CSV.

Maps each TE record in FlyBase ``dmel-all-transposon-*.fasta`` to a
``class/subclass`` label by joining on family name with the bergmanlab
``D_mel_transposon_sequence_set.fa`` (header format: ``>family#class/sub``).
"""

from collections import Counter
import csv
from pathlib import Path
import re
import urllib.request

from loguru import logger
import typer

app = typer.Typer()

DEFAULT_INPUT = Path("data/external/dm/dmel-all-transposon-r6.67.fasta")
DEFAULT_CLASS_FASTA = Path("data/external/dm/D_mel_transposon_sequence_set.fa")
DEFAULT_CLASS_URL = (
    "https://raw.githubusercontent.com/bergmanlab/drosophila-transposons/"
    "master/current/D_mel_transposon_sequence_set.fa"
)
DEFAULT_OUTPUT = Path("data/external/dm_te_class.csv")

# FlyBase TE header: ">FBti0019256 ... name=invader2{}555; ..."
FBTI_RE = re.compile(r"^>(\S+)")
NAME_RE = re.compile(r"name=([^{};\s]+)")
# bergmanlab header: ">family#class/subclass"
CLASS_RE = re.compile(r"^>([^#\s]+)#(\S+)")


def ensure_class_fasta(path: Path, url: str) -> Path:
    """Download the bergmanlab class fasta if missing."""
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Downloading class fasta from {url}")
    urllib.request.urlretrieve(url, path)
    return path


def load_family_class(path: Path) -> dict[str, str]:
    """Parse ``>family#class/subclass`` headers into a dict."""
    mapping: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.startswith(">"):
            continue
        m = CLASS_RE.match(line)
        if not m:
            logger.warning(f"Unparseable class header: {line}")
            continue
        family, klass = m.group(1), m.group(2)
        mapping[family] = klass
    logger.info(f"Loaded {len(mapping)} family->class entries from {path}")
    return mapping


def iter_te_records(path: Path):
    """Yield (fbti_id, family_name) for each FASTA header."""
    for line in path.read_text().splitlines():
        if not line.startswith(">"):
            continue
        m_id = FBTI_RE.match(line)
        m_name = NAME_RE.search(line)
        if not m_id or not m_name:
            logger.warning(f"Unparseable TE header: {line[:80]}")
            continue
        yield m_id.group(1), m_name.group(1)


@app.command()
def main(
    input_fasta: Path = typer.Option(DEFAULT_INPUT, "--input-fasta"),
    class_fasta: Path = typer.Option(DEFAULT_CLASS_FASTA, "--class-fasta"),
    class_url: str = typer.Option(DEFAULT_CLASS_URL, "--class-url"),
    output: Path = typer.Option(DEFAULT_OUTPUT, "--output", "-o"),
    accession: str = typer.Option(
        "fbti", "--accession",
        help="Which ID to use as first column: 'fbti' or 'family'",
    ),
    unknown_label: str = typer.Option("Unknown", "--unknown-label"),
) -> None:
    if accession not in {"fbti", "family"}:
        raise typer.BadParameter("--accession must be 'fbti' or 'family'")

    class_fasta = ensure_class_fasta(class_fasta, class_url)
    fam2class = load_family_class(class_fasta)

    output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter = Counter()
    n_total = 0

    with output.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["accession", "prediction"])
        for fbti, family in iter_te_records(input_fasta):
            # Some FlyBase family names omit the "-element" suffix used by
            # the bergmanlab set (e.g. "Doc3" vs "Doc3-element").
            klass = (
                fam2class.get(family)
                or fam2class.get(f"{family}-element")
                or unknown_label
            )
            counts[klass] += 1
            n_total += 1
            acc = fbti if accession == "fbti" else family
            writer.writerow([acc, klass])

    n_unknown = counts.get(unknown_label, 0)
    logger.info(f"Wrote {n_total} rows to {output}")
    logger.info(f"Unknown: {n_unknown} ({n_unknown / n_total:.1%})")
    for klass, n in counts.most_common(10):
        logger.info(f"  {klass}: {n}")
    logger.success("Done.")


if __name__ == "__main__":
    app()
