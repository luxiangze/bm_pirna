# R Scripts Usage Guide

This directory contains R scripts for *Bombyx mori* (silkworm) piRNA-related data analysis.

## Dependencies

```r
install.packages(c("dplyr", "ggplot2", "tidyr", "scales", "patchwork", "ggrepel",
                   "RColorBrewer", "stringr", "here", "tidyverse"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "enrichplot", "AnnotationForge",
                        "GO.db", "AnnotationDbi", "org.Bmori.eg.db"))
```

## Script Overview

### Data Preparation

| Script | Description |
|--------|-------------|
| `00_build_orgdb.R` | Build custom OrgDb database for *Bombyx mori* from NCBI GAF file, used for GO enrichment analysis |
| `01_sig_splicing_annotate.R` | Filter significant splicing genes from rMATS (FDR<0.05), merge with TPM expression and gene annotation |

### Visualization

| Script | Description |
|--------|-------------|
| `02_pirna_gene_psi_dotplot.R` | Plot PSI (Percent Spliced In) dotplot for known piRNA pathway genes, showing alternative splicing event changes |
| `05_te_rpm_visualization.R` | Transposon RPM matrix visualization, including grouped bar charts, normalized plots, and heatmaps |

### Enrichment Analysis

| Script | Description |
|--------|-------------|
| `03_go_kegg_enrichment.R` | GO (BP/MF/CC) and KEGG enrichment analysis for significant differential splicing genes, generating multiple chart styles |
| `04_deseq2_go_kegg_enrichment.R` | Enrichment analysis for DESeq2 differentially expressed genes (p<0.05, \|LFC\|>1) |

### smRNA-seq Analysis (`smRNA-seq/`)

| Script | Description |
|--------|-------------|
| `mirtrace_length_distribution.R` | miRNA length distribution visualization (RPM normalized, with mean±SD ribbon) |
| `piRNA_length_distribution.R` | piRNA length distribution comparison plot (24-35 nt range) |
| `piRNA_nt_distribution.R` | Nucleotide distribution comparison at 1st and 10th positions of piRNA |
| `td-piRNA_scatter_plot.R` | Transposon-derived piRNA scatter plot (log2FC color mapping) - Abandoned|
| `transposon_scatter_plot.R` | Transposon expression scatter plot (RPM comparison with significant changes colored) |