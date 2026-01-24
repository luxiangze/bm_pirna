"""Standardize sample names by creating symlinks."""

import csv
from pathlib import Path

from loguru import logger
import typer

from bm_pirna.config import RAW_DATA_DIR

app = typer.Typer()


@app.command()
def main(
    sample_map: Path = typer.Argument(..., help="CSV file with original_name,sample_id"),
    input_dir: Path = typer.Option(
        RAW_DATA_DIR,
        "--input-dir",
        "-i",
        help="Directory containing original files",
    ),
    output_dir: Path = typer.Option(
        RAW_DATA_DIR / "renamed",
        "--output-dir",
        "-o",
        help="Directory for creating symlinks",
    ),
    suffix: str = typer.Option(
        ".fq.gz",
        "--suffix",
        "-s",
        help="Suffix for renamed files",
    ),
) -> None:
    """Create standardized symlinks based on sample mapping.

    The sample_map CSV should have columns: original_name, sample_id
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(sample_map) as f:
        reader = csv.DictReader(f)
        for row in reader:
            original = input_dir / row["original_name"]
            renamed = output_dir / f"{row['sample_id']}{suffix}"

            if not original.exists():
                logger.warning(f"File not found: {original}")
                continue

            if renamed.exists():
                renamed.unlink()

            renamed.symlink_to(original.resolve())
            logger.info(f"{row['original_name']} -> {row['sample_id']}{suffix}")

    logger.success(f"Symlinks created in {output_dir}")


if __name__ == "__main__":
    app()
