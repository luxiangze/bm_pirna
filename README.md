# Bombyx mori piRNA Analysis

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

A comprehensive pipeline for analyzing piRNA-related genes in *Bombyx mori* (silkworm), covering smRNA-seq preprocessing, alternative splicing analysis, and downstream visualization.

## Features

- **Quality Control**: Adapter trimming and length filtering using cutadapt
- **Read Filtering**: Bowtie-based alignment to exclude/include specific RNA classes (structural RNA, transposons, tRNA)
- **Statistics**: Read length distribution, nucleotide composition, 1U/10A bias analysis
- ~~**piRNA Matching**: Identify known piRNAs from piRNAdb~~
- **RPM Calculation**: Normalize read counts and calculate fold changes between samples
- **Transposon Analysis**: Align piRNAs to transposon sequences and calculate RPM
- **Alternative Splicing**: Snakemake workflow for differential splicing analysis (STAR + rMATS)
- **R Downstream Analysis**: GO/KEGG enrichment, PSI visualization, TE RPM plots

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/luxiangze/bm_pirna.git
cd bm_pirna

# Install dependencies using pixi
pixi install
```

### smRNA-seq Analysis

#### 1. Quality Control

```bash
# Adapter trimming + length filtering (24-35 nt)
pixi run python -m bm_pirna.smrna_seq.qc /path/to/raw_data \
    -o data/interim/qc \
    -a data/external/sample_adapter.csv \
    -l 24-35 \
    --report-dir reports/qc_reports
```

#### 2. Structural RNA Filtering

```bash
# Exclude structural RNAs (tRNA, rRNA, snRNA, snoRNA) via bowtie
pixi run python -m bm_pirna.smrna_seq.filter data/interim/qc/24-35_length_filtered \
    data/interim/bowtie_index/structure_rna_index \
    -v 1 -m exclude --threads 16 \
    --output-dir data/interim/structure_rna_filtered
```

#### 3. Transposon-Derived piRNA

```bash
# Retain reads aligning to transposons (0 mismatches)
pixi run python -m bm_pirna.smrna_seq.filter data/interim/structure_rna_filtered \
    data/interim/bowtie_index/transposon_index \
    -v 0 -m include --threads 16 \
    --output-dir data/interim/transposon_filtered

# Generate read statistics
pixi run python -m bm_pirna.smrna_seq.stats data/interim/transposon_filtered \
    -o data/processed/transposon_pirna_stats

# Calculate transposon RPM and fold change
pixi run python -m bm_pirna.smrna_seq.transposon_rpm data/processed/transposon_pirna_stats/collapsed \
    -i data/interim/bowtie_index/transposon_index \
    -c data/interim/structure_rna_filtered/filtered_read_counts.tsv \
    -s data/external/sample_map.csv \
    -o data/processed/transposon_rpm
```

#### 4. piRNA Matching

```bash
# Match reads to known piRNAs from piRNAdb
pixi run python -m bm_pirna.smrna_seq.pirna_match data/processed/stats/collapsed \
    -o data/processed/pirna_match
```

#### 5. RPM and Fold Change

```bash
# Calculate piRNA RPM and fold change between samples
pixi run python -m bm_pirna.smrna_seq.pirna_rpm data/processed/pirna_match/collapsed \
    -c reports/qc_reports/cutadapt/cutadapt_summary.tsv \
    -s data/external/sample_map.csv \
    -o data/processed/pirna_rpm
```

### Alternative Splicing Workflow

See [`workflow/rnasplicing/README.md`](workflow/rnasplicing/README.md) for full documentation.

```bash
cd workflow/rnasplicing

# Dry run
snakemake -n

# Execute (4 parallel jobs)
snakemake -j 4 --use-conda
```

### R Downstream Analysis

See [`R_scripts/README.md`](R_scripts/README.md) for full documentation.

```bash
# Example: Build OrgDb → annotate splicing genes → enrichment analysis
Rscript R_scripts/00_build_orgdb.R
Rscript R_scripts/01_sig_splicing_annotate.R
Rscript R_scripts/03_go_kegg_enrichment.R
```

## Input Files

### sample_adapter.csv

```csv
sample,adapter,rm_adapter
WT_rep1,TGGAATTCTCGGGTGCCAAGG,true
WT_rep2,TGGAATTCTCGGGTGCCAAGG,true
SUGP1_rep1,AACTGTAGGCACCATCAAT,false
```

### sample_map.csv

```csv
Control_sample,Treated_sample
WT,SUGP1
```

## Project Organization

```
├── LICENSE                <- MIT License
├── README.md              <- This file
├── pixi.toml              <- Pixi package manager configuration
├── pyproject.toml         <- Python project metadata
│
├── notebooks/             <- Jupyter notebooks for exploratory analysis
│   ├── bm_pirna_smrna-seq.ipynb     <- smRNA-seq analysis walkthrough
│   └── bm_pirna_mirna_stats.ipynb   <- miRNA statistics
│
├── workflow/
│   └── rnasplicing/       <- Snakemake workflow for alternative splicing
│       ├── Snakefile      <- Workflow definition (STAR + rMATS)
│       ├── config.json    <- Sample groups and comparison configuration
│       └── env.yml        <- Conda environment
│
├── R_scripts/             <- R scripts for downstream analysis
│   ├── 00_build_orgdb.R           <- Build Bombyx mori OrgDb
│   ├── 01_sig_splicing_annotate.R <- Annotate significant splicing genes
│   ├── 02_pirna_gene_psi_dotplot.R <- PSI dot plot for piRNA genes
│   ├── 03_go_kegg_enrichment.R    <- GO/KEGG enrichment (splicing genes)
│   ├── 04_deseq2_go_kegg_enrichment.R <- GO/KEGG enrichment (DEGs)
│   ├── 05_te_rpm_visualization.R  <- Transposon RPM visualization
│   └── smRNA-seq/                 <- smRNA-seq R visualization scripts
│
└── bm_pirna/              <- Python source package
    ├── config.py          <- Project paths configuration
    ├── dataset.py         <- Data download utilities (S3/rustfs)
    ├── extract_gene_info.py  <- Extract gene info from GTF
    ├── plot_aa_composition.py <- Amino acid composition plots
    ├── utils/             <- Shared utilities (run_cmd)
    └── smrna_seq/         <- smRNA-seq analysis modules
        ├── qc.py              <- Quality control (cutadapt + fastp)
        ├── filter.py          <- Bowtie-based read filtering
        ├── stats.py           <- Read statistics and collapsed FASTA
        ├── rename.py          <- Sample renaming utilities
        ├── pirna_match.py     <- piRNA identification against piRNAdb
        ├── pirna_rpm.py       <- piRNA RPM and fold change
        └── transposon_rpm.py  <- Transposon-level RPM calculation
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

Yongkang Guo

--------
