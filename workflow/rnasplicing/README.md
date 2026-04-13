# rnasplicing_20260309

A Snakemake-based workflow for transcriptome-level alternative splicing analysis.

Input FASTQ files are specified in `samplesheet.csv`. The workflow supports retrieval of sequencing data from S3-compatible object storage using `mc`.

## Overview

This workflow performs the following major steps:

1. Read sample metadata from `samplesheet.csv`
2. Download paired-end FASTQ files from object storage
3. Build a STAR genome index
4. Align reads to the reference genome
5. Generate group-specific BAM lists
6. Run rMATS for differential splicing analysis
7. Summarize the results and generate plots

## Project Structure

- `Snakefile`: Main workflow definition
- `config.json`: Workflow configuration
- `samplesheet.csv`: Sample sheet
- `env.yml`: Runtime environment
- `genome/`: Genome file directory
- `scripts/`: Helper scripts

## Input Files

### Sample Sheet

The input sample manifest should be provided as `samplesheet.csv` with the following structure:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,s3://bucket/path/sample1_R1.fq.gz,s3://bucket/path/sample1_R2.fq.gz
CONTROL_REP2,s3://bucket/path/sample2_R1.fq.gz,s3://bucket/path/sample2_R2.fq.gz
```

Column description:

- `sample`: Sample name
- `fastq_1`: S3 path to the R1 FASTQ file
- `fastq_2`: S3 path to the R2 FASTQ file

### Configuration File

Key parameters in `config.json` include:

```json
{
  "samplesheet": "samplesheet.csv",
  "genome_fa": "genome/genomic_chr.fa",
  "genome_gtf": "genome/genomic_chr.gtf",
  "mc_alias": "rustfs"
}
```

Parameter description:

- `samplesheet`: Path to the sample sheet
- `genome_fa`: Reference genome FASTA
- `genome_gtf`: Primary annotation GTF
- `extra_gtf`: *(optional)* Secondary annotation GTF to be merged with `genome_gtf` prior to rMATS analysis (e.g. a piRNA cluster annotation)
- `cluster_gtf`: *(optional)* GTF file of piRNA cluster regions (proTRAC output); when set, enables the `quantify_clusters` branch to produce a per-cluster count matrix across all samples
- `mc_alias`: Preconfigured object storage alias in `mc`
- `samples`: Group definitions
- `comparisons`: Group comparisons

## Requirements and Preparation

### Software Environment

Create the conda environment if it has not been initialized:

```bash
mamba env create -f env.yml
```

Activate the environment:

```bash
mamba activate rnasplicing
```

### Object Storage Configuration

The workflow uses `mc cp` to retrieve input data from object storage.

Verify that the configured alias is available:

```bash
mc alias list
```

The current workflow is configured to use:

```bash
rustfs
```

If another alias should be used, such as `minio`, update `mc_alias` in `config.json` accordingly.

### Reference Files

Ensure that the following reference files are present:

```bash
genome/genomic_chr.fa
genome/genomic_chr.gtf
```

## How to Run

Run the workflow from the project directory:

```bash
snakemake -j 4 --use-conda
```

If the environment has already been activated manually, the workflow can also be started with:

```bash
snakemake -j 4
```

A dry run is recommended prior to the first execution:

```bash
snakemake -n
```

To display the shell commands that would be executed:

```bash
snakemake -n -p
```

To resume interrupted or incomplete jobs:

```bash
snakemake -j 4 --rerun-incomplete
```

## Analytical Workflow

The analytical procedure implemented by this workflow consists of the following steps:

1. Read `samplesheet.csv`
2. Download FASTQ files from object storage using `mc cp`
3. Build the genome index with STAR
4. Align reads with STAR and generate BAM files
5. Generate group files from `config.json`
6. Run alternative splicing analysis with rMATS
7. Summarize results and generate the final report

## Output Files

### Primary Output

- `results/{comparison}/summary.txt`: Summary report for each group comparison
- `results/{comparison}/volcano.pdf`: Volcano plot for each group comparison

### Intermediate and Auxiliary Outputs

- `fastq/{sample}_R1.fastq.gz`
- `fastq/{sample}_R2.fastq.gz`
- `star/{sample}/Aligned.sortedByCoord.out.bam`
- `groups/{group}.txt`
- `groups/comparison_index.txt`
- `rmats/{comparison}/*.MATS.JC.txt`
- `rmats/{comparison}/*.MATS.JCEC.txt`
- `results/{comparison}/`
- `results/cluster_counts.txt` *(produced only when `cluster_gtf` is set)*: raw count matrix (clusters × samples) for downstream differential expression analysis
- `logs/`

Notes:

- Files under `fastq/` are declared as `temp()` and may be removed automatically by Snakemake after successful completion
- The principal workflow deliverables are `results/{comparison}/summary.txt` and `results/{comparison}/volcano.pdf`

## Utility Commands

### Visualize the Workflow DAG

```bash
snakemake --dag | dot -Tpng > dag.png
```

### Build a Specific Target

For example, to generate a single aligned BAM file:

```bash
snakemake -j 4 star/CONTROL_REP1/Aligned.sortedByCoord.out.bam
```

### Resume Incomplete Runs

```bash
snakemake -j 4 --rerun-incomplete
```

## Troubleshooting

### `mc cp` Download Failure

Check the object storage alias configuration:

```bash
mc alias list
mc ls rustfs
```

If the alias is not `rustfs`, update `config.json` accordingly.

### Reference Genome Files Not Found

Check the genome directory:

```bash
ls genome/
```

Make sure `genomic_chr.fa` and `genomic_chr.gtf` are present.

### Sample-Level Alignment Failure

Inspect the STAR log file:

```bash
less logs/star/sample_name.log
```

### rMATS Failure

Inspect the rMATS log file:

```bash
less logs/rmats/comparison_1_control_vs_sugp1.log
```

## Current Experimental Design

The current `config.json` defines:

- `control`: `CONTROL_REP1`, `CONTROL_REP2`
- `sugp1`: `SUGP1_REP1`, `SUGP1_REP2`

Current comparison:

- `control vs sugp1`
