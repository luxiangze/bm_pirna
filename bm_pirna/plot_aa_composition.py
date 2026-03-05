from collections import Counter
from pathlib import Path

from loguru import logger
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import typer

from bm_pirna.config import FIGURES_DIR, RAW_DATA_DIR

app = typer.Typer()

# 20 standard amino acids
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# Use a colormap with enough distinct colors
CMAP = plt.get_cmap("tab20", len(AMINO_ACIDS))
AA_COLORS = {aa: CMAP(i) for i, aa in enumerate(AMINO_ACIDS)}


def parse_fasta(fasta_path: Path) -> dict[str, dict[str, str]]:
    """Parse FASTA file and group sequences by gene name.

    Returns:
        {gene: {species: sequence}}
    """
    genes: dict[str, dict[str, str]] = {}
    current_gene = None
    current_species = None
    seq_parts: list[str] = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous sequence
                if current_gene and current_species:
                    genes.setdefault(current_gene, {})[current_species] = "".join(
                        seq_parts
                    ).upper()
                # Parse header: >gene_species (gene may contain underscores?)
                # Convention: first segment before first underscore is NOT reliable.
                # Actual format: first token is <gene>_<Species>_<name>...
                # From data: ints4_Bombyx_mori -> gene=ints4, species=Bombyx_mori
                header = line[1:]
                idx = header.index("_")
                current_gene = header[:idx]
                current_species = header[idx + 1 :]
                seq_parts = []
            else:
                seq_parts.append(line)

    # Save last sequence
    if current_gene and current_species:
        genes.setdefault(current_gene, {})[current_species] = "".join(
            seq_parts
        ).upper()

    return genes


def calc_aa_counts(sequence: str) -> dict[str, int]:
    """Count occurrences of each of the 20 standard AAs."""
    counts = Counter(sequence)
    return {aa: counts.get(aa, 0) for aa in AMINO_ACIDS}


def plot_gene(gene: str, species_seqs: dict[str, str], output_dir: Path) -> None:
    """Draw a horizontal stacked bar chart for one gene."""
    species_list = list(species_seqs.keys())
    compositions = [calc_aa_counts(species_seqs[sp]) for sp in species_list]

    fig, ax = plt.subplots(figsize=(14, max(3, len(species_list) * 0.6 + 1.5)))

    lefts = np.zeros(len(species_list))
    for aa in AMINO_ACIDS:
        widths = np.array([comp[aa] for comp in compositions])
        ax.barh(
            species_list,
            widths,
            left=lefts,
            color=AA_COLORS[aa],
            label=aa,
            edgecolor="none",
            height=0.7,
        )
        lefts += widths

    ax.set_xlabel("Amino acid count")
    ax.set_title(f"Amino acid composition — {gene}")
    ax.set_xlim(0, lefts.max() * 1.01)
    ax.tick_params(axis="y", labelsize=9)

    # Legend outside plot
    handles = [mpatches.Patch(color=AA_COLORS[aa], label=aa) for aa in AMINO_ACIDS]
    ax.legend(
        handles=handles,
        title="AA",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=8,
        title_fontsize=9,
        ncol=1,
    )

    plt.tight_layout()
    out_path = output_dir / f"aa_composition_{gene}.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved: {out_path}")


@app.command()
def main(
    input_path: Path = RAW_DATA_DIR / "Intersection_genes.fasta",
    output_dir: Path = FIGURES_DIR,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Parsing {input_path} ...")
    genes = parse_fasta(input_path)
    logger.info(f"Found {len(genes)} genes: {list(genes.keys())}")

    for gene, species_seqs in genes.items():
        logger.info(f"Plotting gene: {gene} ({len(species_seqs)} species)")
        plot_gene(gene, species_seqs, output_dir)

    logger.success("All plots done.")


if __name__ == "__main__":
    app()
