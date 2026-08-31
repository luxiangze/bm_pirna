# R Scripts Usage Guide

This directory contains R scripts for *Bombyx mori* (silkworm) piRNA-related
analysis, plus two *Drosophila melanogaster* enrichment workflows.

Run scripts from the repository root so `here::here()` resolves CCDS
paths (`data/processed/...`, `reports/figures/...`):

```bash
Rscript R_scripts/<script_name>.R
```

Shared path constants are in `R_scripts/R/project_paths.R`. Do not write
to Mac-era folders (`graphs/`, `data/output/`, `data/rnadiff_results/`).

## Dependencies

Install common CRAN packages:

```r
install.packages(c(
  "dplyr", "ggplot2", "tidyr", "scales", "patchwork", "ggrepel",
  "RColorBrewer", "stringr", "here", "tidyverse"
))
```

Install Bioconductor packages:

```r
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "clusterProfiler", "enrichplot", "AnnotationForge", "GO.db",
  "AnnotationDbi", "org.Dm.eg.db"
))
```

The *Bombyx mori* OrgDb used by scripts `03` and `04` is bundled as a local
source package (`R_scripts/org.Bmori.eg.db_0.3.0.tar.gz` or
`data/external/bm_ncbi/org.Bmori.eg.db_*.tar.gz`):

```r
install.packages(
  "R_scripts/org.Bmori.eg.db_0.3.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

To rebuild the *Bombyx mori* OrgDb from the NCBI GAF in
`data/external/bm_ncbi/`, run:

```bash
Rscript R_scripts/00_build_orgdb.R
```

## Script Overview

### Data Preparation

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `00_build_orgdb.R` | `data/external/bm_ncbi/GCF_030269925.1-RS_2024_01_gene_ontology.gaf` | `data/external/bm_ncbi/org.Bmori.eg.db_<version>.tar.gz` and installed `org.Bmori.eg.db` | Build a custom *Bombyx mori* OrgDb from NCBI GAF annotations for GO enrichment |
| `01_sig_splicing_annotate.R` | rMATS summary, TPM matrix, gene annotation | `data/output/sig_splicing_genes_annotated.tsv` | Filter significant splicing genes (`FDR < 0.05`), merge TPM expression and gene annotation |

### Visualization

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `02_pirna_gene_psi_dotplot.R` | known piRNA gene list and rMATS event tables | `graphs/pirna_gene_psi_dotplot.*` | Plot PSI dotplots for known piRNA pathway genes |
| `05_te_rpm_visualization.R` | `data/transposon_rpm/te_top_class_rpm_matrix.tsv`, `data/transposon_rpm/te_sub_class_rpm_matrix.tsv` | `graphs/transposon_rpm/` | Visualize transposon RPM by top-level TE class and sub-class |

### Enrichment Analysis

| Script | Species | Input | Output | Description |
|--------|---------|-------|--------|-------------|
| `03_go_kegg_enrichment.R` | *Bombyx mori* | `data/output/sig_splicing_genes_annotated.tsv` | `data/output/enrichment/`, `graphs/enrichment/` | GO (BP/MF/CC) and KEGG enrichment for significant alternative splicing genes |
| `04_deseq2_go_kegg_enrichment.R` | *Bombyx mori* | `workflow/bm_differentialabundance_20260821/results/tables/differential/contrasts/*_filtered.tsv`, `data/external/gene_info_table.tsv` | `data/processed/sugp1_rnaseq_20260821/go_kegg/{control_vs_sugp1,p200_vs_sugp1}/`, `reports/figures/sugp1_rnaseq_20260821/go_kegg/...` | GO and KEGG enrichment for DESeq2 DEGs (`padj < 0.05`, `|log2FoldChange| > 1`) for SUGP1 vs Control and SUGP1 vs P200. Uses the installed `org.Bmori.eg.db` package, not the tar.gz. |
| `06_dm_deseq2_go_kegg_enrichment.R` | *Drosophila melanogaster* | `data/dm_differential/condition_control_treated.deseq2.results.tsv` | `data/output/dm_enrichment_deseq2/{all,up,down}/`, `graphs/dm_enrichment_deseq2/` | GO and KEGG enrichment for all, up-regulated, and down-regulated DESeq2 DEGs (`padj < 0.05`, `|log2FoldChange| > 1`) |
| `07_dm_splicing_go_kegg_enrichment.R` | *Drosophila melanogaster* | `data/dm_splicing_20260425/summary.JC.tsv` | `data/output/dm_enrichment_splicing/`, `graphs/dm_enrichment_splicing/` | GO and KEGG enrichment for significant rMATS JC splicing genes (`FDR < 0.05`) |

Notes:

- `03_go_kegg_enrichment.R` and `04_deseq2_go_kegg_enrichment.R` require the
  local `org.Bmori.eg.db` package and use KEGG organism code `bmor`.
- `06_dm_deseq2_go_kegg_enrichment.R` and `07_dm_splicing_go_kegg_enrichment.R`
  require `org.Dm.eg.db`, use FlyBase IDs as input, and use KEGG organism code
  `dme`.
- Enrichment scripts save tabular results as TSV files and publication-style
  PDF/PNG plots including bubble, lollipop, dot, combined GO, cnetplot, and
  upset-style visualizations where applicable.

### smRNA-seq Analysis (`smRNA-seq/`)

| Script | Description |
|--------|-------------|
| `mirtrace_length_distribution.R` | miRNA length distribution visualization (RPM normalized, with mean +/- SD ribbon) |
| `piRNA_length_distribution.R` | piRNA length distribution comparison plot (24-35 nt range) |
| `piRNA_nt_distribution.R` | Nucleotide distribution comparison at 1st and 10th positions of piRNA |
| `td-piRNA_scatter_plot.R` | Transposon-derived piRNA scatter plot (log2FC color mapping) - abandoned |
| `transposon_scatter_plot.R` | Transposon expression scatter plot (RPM comparison with significant changes colored) |

## Typical Run Order

For *Bombyx mori* splicing enrichment:

```bash
Rscript R_scripts/01_sig_splicing_annotate.R
Rscript R_scripts/03_go_kegg_enrichment.R
```

For *Bombyx mori* DESeq2 enrichment (SUGP1 vs Control and SUGP1 vs P200):

```bash
Rscript R_scripts/04_deseq2_go_kegg_enrichment.R
```

For *Drosophila melanogaster* enrichment:

```bash
Rscript R_scripts/06_dm_deseq2_go_kegg_enrichment.R
Rscript R_scripts/07_dm_splicing_go_kegg_enrichment.R
```
