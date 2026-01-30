# Bombyx mori piRNA Analysis

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

A comprehensive pipeline for analyzing piRNA-related genes in *Bombyx mori* (silkworm) using smRNA-seq data.

## Features

- **Quality Control**: Adapter trimming and length filtering using cutadapt
- **Statistics**: Read length distribution, nucleotide composition, 1U/10A bias analysis
- **piRNA Matching**: Identify known piRNAs from piRNAdb
- **RPM Calculation**: Normalize read counts and calculate fold changes between samples
- **Transposon Analysis**: Align piRNAs to transposon sequences and calculate RPM

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/your-username/bm_pirna.git
cd bm_pirna

# Install dependencies using pixi
pixi install
```

### Basic Usage

#### 1. Quality Control

```bash
# Run QC on smRNA-seq data (adapter trimming + length filtering 18-37nt)
pixi run python -m bm_pirna.smrna_seq.qc /path/to/raw_data \
    -o data/interim/qc \
    -a data/external/sample_adapter.csv \
    --report-dir reports/qc_reports
```

#### 2. Generate Statistics

```bash
# Generate read statistics and collapsed fasta files
pixi run python -m bm_pirna.smrna_seq.stats data/interim/qc \
    -o data/processed/stats
```

#### 3. piRNA Matching

```bash
# Match reads to known piRNAs from piRNAdb
pixi run python -m bm_pirna.smrna_seq.pirna_match data/processed/stats/collapsed \
    -o data/processed/pirna_match
```

#### 4. Calculate RPM and Fold Change

```bash
# Calculate piRNA RPM and fold change between samples
pixi run python -m bm_pirna.smrna_seq.pirna_rpm data/processed/pirna_match \
    -c reports/qc_reports/cutadapt/cutadapt_summary.tsv \
    -s data/external/sample_map.csv \
    -o data/processed/pirna_rpm
```

#### 5. Transposon Analysis

```bash
# Align piRNAs to transposons and calculate RPM
pixi run python -m bm_pirna.smrna_seq.transposon_rpm data/processed/stats/collapsed \
    -i data/interim/bowtie_index/transposon_index \
    -c reports/qc_reports/cutadapt/cutadapt_summary.tsv \
    -o data/processed/transposon_rpm
```

## Input Files

### sample_adapter.csv

```csv
sample,adapter,rm_adapter
WT_rep1,TGGAATTCTCGGGTGCCAAGG,true
WT_rep2,TGGAATTCTCGGGTGCCAAGG,true
KO_rep1,AACTGTAGGCACCATCAAT,false
```

### sample_map.csv

```csv
Control_sample,Treated_sample
WT,KO
WT,Mutant
```

## Project Organization

```
├── LICENSE            <- MIT License
├── README.md          <- This file
├── pixi.toml          <- Pixi package manager configuration
├── data
│   ├── external       <- Sample metadata and reference files
│   ├── interim        <- Intermediate data (QC output, indexes)
│   ├── processed      <- Final analysis results
│   └── raw            <- Raw smRNA-seq data
│
├── notebooks          <- Jupyter notebooks for analysis
│
├── reports            <- QC reports and figures
│
└── bm_pirna           <- Source code
    ├── config.py      <- Project configuration
    ├── dataset.py     <- Data download utilities
    └── smrna_seq/     <- smRNA-seq analysis modules
        ├── qc.py          <- Quality control (cutadapt)
        ├── stats.py       <- Read statistics
        ├── pirna_match.py <- piRNA identification
        ├── pirna_rpm.py   <- piRNA RPM calculation
        └── transposon_rpm.py <- Transposon analysis
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

Yongkang Guo

--------

