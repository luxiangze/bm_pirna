# Scripts for RNA Splicing Workflow

This directory contains three helper scripts used by the RNA splicing workflow.

## Overview

- **`group_files.py`**
  - Generate sample-group BAM list files and comparison index files from `config.json` and STAR alignment outputs.

- **`summarize_splicing.py`**
  - Read rMATS event tables from one comparison, merge event types, export a summary table, and draw a volcano plot.

- **`analyze_clusters.py`**
  - Analyze piRNA cluster expression from `featureCounts` output, compute RPK and log2 fold change, and generate summary plots.

---

## 1. `group_files.py`

### Purpose

This script scans the STAR BAM output directory, groups BAM files according to the sample definitions in `config.json`, and writes text files required by downstream rMATS steps.

### Main inputs

- **`--bam-dir`**
  - Directory containing STAR alignment results.
  - Expected layout:

```text
star/
  CONTROL_REP1/
    Aligned.sortedByCoord.out.bam
  CONTROL_REP2/
    Aligned.sortedByCoord.out.bam
  SUGP1_REP1/
    Aligned.sortedByCoord.out.bam
```

- **`--config`**
  - JSON workflow configuration file.
  - Required fields used by this script:

```json
{
  "samples": {
    "control": ["CONTROL_REP1", "CONTROL_REP2"],
    "sugp1": ["SUGP1_REP1", "SUGP1_REP2"]
  },
  "comparisons": [
    ["control", "sugp1"]
  ]
}
```

- **`--output-dir`**
  - Directory for generated group and comparison files.

### Main outputs

- **`{output_dir}/{group}.txt`**
  - One file per sample group.
  - Format: a single line of comma-separated BAM file paths.

```text
star/CONTROL_REP1/Aligned.sortedByCoord.out.bam,star/CONTROL_REP2/Aligned.sortedByCoord.out.bam
```

- **`{output_dir}/comparison_{i}_{group1}_vs_{group2}.txt`**
  - One file per comparison.
  - Format: one line with two group names separated by a tab.

```text
control	sugp1
```

- **`{output_dir}/comparison_index.txt`**
  - Index of all comparisons.
  - Format: one comparison per line.

```text
comparison_1_control_vs_sugp1	control	sugp1
```

### Example

```bash
python group_files.py \
  --bam-dir star \
  --config config.json \
  --output-dir groups
```

---

## 2. `summarize_splicing.py`

### Purpose

This script reads rMATS result tables for one comparison, combines all supported event types into a single summary table, and generates a volcano plot based on `IncLevelDifference` and `FDR`.

It also optionally labels known piRNA-related genes using an annotation CSV file with `is_known` and `ENTREZ.ID` columns.

### Supported rMATS event files

The script checks these files under `--input-dir`:

- `A3SS.MATS.JC.txt`
- `A5SS.MATS.JC.txt`
- `MXE.MATS.JC.txt`
- `RI.MATS.JC.txt`
- `SE.MATS.JC.txt`

### Main inputs

- **`--input-dir`**
  - Directory containing rMATS result files for one comparison.

- **`--output-file`**
  - Path to the merged summary TSV file.

- **`--plot-file`**
  - Path to the output volcano plot PDF.

- **`--fdr`**
  - FDR threshold used to color points.
  - Default: `0.05`

- **`--annotation-file`**
  - Optional CSV file for point labeling.
  - Expected columns:

```text
sample_id,silkbase_id,ENTREZ.ID,Gene_Symbol,is_known,Gene_Name
```

### Expected input table format

Each rMATS table must contain at least these columns:

```text
ID
GeneID
geneSymbol
PValue
FDR
IncLevel1
IncLevel2
IncLevelDifference
```

`GeneID` may be either:

- a numeric identifier such as `101741153`
- a `LOC`-prefixed identifier such as `LOC101741153`

When `--annotation-file` is provided, the script extracts the numeric part from `LOC...` identifiers and matches it to `ENTREZ.ID`. Only rows with `is_known=true` in the annotation file are labeled.

### Main outputs

- **`summary.txt` / custom summary path**
  - Tab-separated merged event table.
  - Output columns:

```text
ID	GeneID	geneSymbol	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
```

- **`volcano.pdf` / custom PDF path**
  - Volcano plot with:
    - gray points for non-significant events
    - blue points for negative significant events
    - red points for positive significant events
    - optional labels for known piRNA genes

### Example

```bash
python summarize_splicing.py \
  --input-dir rmats/comparison_1_control_vs_sugp1 \
  --output-file results/comparison_1_control_vs_sugp1/summary.txt \
  --plot-file results/comparison_1_control_vs_sugp1/volcano.pdf \
  --annotation-file ../../data/external/piRNA_annotations_id.csv
```

---

## 3. `analyze_clusters.py`

### Purpose

This script analyzes piRNA cluster expression from a `featureCounts` count matrix. It performs length normalization to RPK, computes group-level mean log2 fold change, evaluates replicate reproducibility, and generates summary plots.

### Main inputs

- **`--counts`**
  - `featureCounts` output table.

- **`--config`**
  - Workflow `config.json` containing sample group definitions.

- **`--group1`**
  - Reference group name.

- **`--group2`**
  - Treatment group name.

- **`--output-dir`**
  - Directory for result files.

- **`--min-rpk`**
  - Minimum mean RPK used for filtering.
  - Default: `1.0`

- **`--log2fc`**
  - Absolute log2FC threshold used in plots.
  - Default: `1.0`

- **`--pseudocount`**
  - Value added before log2 transform.
  - Default: `0.5`

### Expected input count matrix format

The script expects a standard `featureCounts` output with:

- metadata columns:
  - `Geneid`
  - `Chr`
  - `Start`
  - `End`
  - `Strand`
  - `Length`
- one count column per BAM file

Example header:

```text
Geneid	Chr	Start	End	Strand	Length	star/CONTROL_REP1/Aligned.sortedByCoord.out.bam	star/CONTROL_REP2/Aligned.sortedByCoord.out.bam
```

### Main outputs

- **`{group2}_vs_{group1}_clusters.tsv`**
  - Tab-separated result table.
  - Includes:
    - cluster ID
    - mean RPK in each group
    - mean log2FC
    - per-replicate log2FC columns such as `log2FC_REP1`, `log2FC_REP2`

- **`{group2}_vs_{group1}_rep_correlation.pdf`**
  - Scatter plot comparing replicate-level log2FC values.
  - Generated only when at least two replicate pairs are available.

- **`{group2}_vs_{group1}_log2fc_distribution.pdf`**
  - Histogram of mean log2FC values across clusters.

### Example

```bash
python analyze_clusters.py \
  --counts results/cluster_counts.txt \
  --config config.json \
  --group1 control \
  --group2 sugp1 \
  --output-dir results/cluster_analysis
```

---

## Notes

- `group_files.py` and `summarize_splicing.py` use `argparse`.
- `analyze_clusters.py` uses `typer`.
- Output plots are saved as PDF files.
- `summarize_splicing.py` currently uses `JC` rMATS files only.
