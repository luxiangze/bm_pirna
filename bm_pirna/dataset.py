"""Download raw data from rustfs to data/raw."""

from pathlib import Path
import subprocess

from loguru import logger
import pandas as pd
import typer

from bm_pirna.config import (
    DATA_DIR,
    EXTERNAL_DATA_DIR,
    INTERIM_DATA_DIR,
    PROCESSED_DATA_DIR,
    RAW_DATA_DIR,
)

app = typer.Typer()

# Data directory structure
DATA_SUBDIRS = [
    RAW_DATA_DIR,
    INTERIM_DATA_DIR,
    PROCESSED_DATA_DIR,
    EXTERNAL_DATA_DIR,
]


def ensure_data_dirs() -> None:
    """Ensure all data directories exist."""
    for subdir in DATA_SUBDIRS:
        subdir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Ensured data directories exist under {DATA_DIR}")


def download_file(remote_path: str, local_path: Path) -> bool:
    """Download file from rustfs using mc.

    Args:
        remote_path: Remote path in rustfs (e.g., rustfs/bucket/path/file)
        local_path: Local destination path

    Returns:
        True if download successful, False otherwise
    """
    local_path.parent.mkdir(parents=True, exist_ok=True)

    if local_path.exists():
        logger.info(f"  File already exists: {local_path}")
        return True

    cmd = ["mc", "cp", remote_path, str(local_path)]
    logger.info(f"  Downloading: {remote_path}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"  Download failed: {e.stderr}")
        return False


@app.command()
def main(
    dataset_csv: Path = typer.Argument(
        RAW_DATA_DIR / "dataset.csv",
        help="Path to dataset CSV file",
    ),
    output_dir: Path = typer.Option(
        RAW_DATA_DIR,
        "--output-dir",
        "-o",
        help="Directory to store downloaded raw data",
    ),
    dataset_filter: str = typer.Option(
        None,
        "--dataset",
        "-d",
        help="Filter by dataset name (optional)",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        "-n",
        help="Show what would be done without actually doing it",
    ),
) -> None:
    """Download raw data from rustfs to data/raw.

    Reads dataset.csv to determine which files to download from rustfs,
    and stores them directly in data/raw/{dataset}/.

    Also ensures the data directory structure exists (data/raw, data/interim,
    data/processed, data/external) for git collaboration.

    CSV format:
        remote_path,local_name,dataset
        rustfs/bucket/path/file.fq.gz,sample1.fq.gz,smRNA-seq_20260122

    Examples:
        python -m bm_pirna.dataset data/raw/dataset.csv
        python -m bm_pirna.dataset data/raw/dataset.csv -d smRNA-seq_20260122
        python -m bm_pirna.dataset data/raw/dataset.csv --dry-run
    """
    # Ensure data directories exist
    ensure_data_dirs()

    if not dataset_csv.exists():
        logger.error(f"Dataset CSV not found: {dataset_csv}")
        raise typer.Exit(1)

    # Read dataset CSV
    df = pd.read_csv(dataset_csv)
    required_cols = ["remote_path", "local_name", "dataset"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.error(f"Missing columns in CSV: {missing_cols}")
        raise typer.Exit(1)

    # Filter by dataset if specified
    if dataset_filter:
        df = df[df["dataset"] == dataset_filter]
        if df.empty:
            logger.warning(f"No entries found for dataset: {dataset_filter}")
            return

    logger.info(f"Processing {len(df)} files")
    logger.info(f"Output directory: {output_dir}")

    if dry_run:
        logger.info("DRY RUN - no changes will be made")

    success_count = 0
    fail_count = 0

    for _, row in df.iterrows():
        remote_path = row["remote_path"]
        local_name = row["local_name"]
        dataset = row["dataset"]

        # Output path
        output_path = output_dir / dataset / local_name

        logger.info(f"Processing: {local_name} ({dataset})")

        if dry_run:
            logger.info(f"  Would download: {remote_path} -> {output_path}")
            success_count += 1
            continue

        # Download file
        if download_file(remote_path, output_path):
            success_count += 1
        else:
            fail_count += 1

    logger.info(f"Completed: {success_count} success, {fail_count} failed")
    if fail_count == 0:
        logger.success("All files processed successfully")
    else:
        logger.warning(f"{fail_count} files failed")


if __name__ == "__main__":
    app()
