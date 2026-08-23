## 1. Language & Formatting

- **Output**: Respond exclusively in Chinese.
- **Comments**: All code comments (including docstrings and inline comments) must be written in English.
- **Style**: Follow the idiomatic style and best practices of the target language.

## 2. Code Quality

- **Simplicity First**: Write concise, readable, and practical code. Avoid over-engineering and premature optimization.
- **Reuse & Complexity**: Keep functions focused and single-purpose; extract reusable logic to reduce duplication and cyclomatic complexity.
- **Design Patterns**: Apply patterns only when they genuinely reduce complexity — never force them into simple logic.
- **Code formatting and linting**: Use ruff

## 3. Refactoring & Modification

- **Least Modification**: When changing existing code, confine changes to the target module and avoid unintended side effects on adjacent modules.

## 4. Tooling

- **Package Management**: Use `pixi` by default; use `renv` for R projects. Override only when explicitly instructed.

## 5. Bioinformatics-Specific

- **Reproducibility**: Always set random seeds where applicable; pin tool and package versions explicitly.
- **Data Assumptions**: Never silently assume input format, genome build, or sample metadata — validate inputs and state assumptions explicitly in comments.
- **Pipeline Awareness**: Structure code to be modular and compatible with workflow managers (e.g., Snakemake, Nextflow) when appropriate.
- **Result Transparency**: Include brief output summaries (e.g., number of genes, samples filtered) at key steps via logging or inline comments.

## 6. Directory Layout (CCDS + Bioinformatics)

This repository follows [Cookiecutter Data Science](https://cookiecutter-data-science.drivendata.org/) with bioinformatics additions. Put new files in the matching directory; do not invent extra top-level data/figure folders.

### Canonical layout

| Path | Role | Mutable? |
|------|------|----------|
| `data/raw/` | Original downloads (FASTQ, vendor dumps). Organized as `{dataset}/{assay}/`. Recorded in `data/raw/dataset.csv`. | No (treat as immutable) |
| `data/interim/` | Intermediate files that are transformed but not analysis-ready (trimmed FASTQ, bowtie hits, QC intermediates). Date-stamp batches: `{experiment}_{YYYYMMDD}/`. | Yes, regenerable |
| `data/processed/` | Analysis-ready tables and derived matrices (counts used by notebooks/R, enrichment TSVs, RPM tables). Same batch naming as interim. | Yes, regenerable |
| `data/external/` | Third-party references that are not raw sequencing: GTF/GFF/GAF, gene info tables, TE libraries, OrgDb sources, sample maps. | Rarely; versioned by source |
| `reports/` | Human-facing QC and pipeline reports (FastQC, MultiQC, cutadapt JSON). | Yes |
| `reports/figures/` | Publication figures produced by R/Python. Mirror the processed batch name: `reports/figures/{experiment}_{YYYYMMDD}/`. | Yes |
| `notebooks/` | Exploratory analysis only. Do not write final tables or figures from notebooks into `data/raw` or repo root. | Yes |
| `bm_pirna/` | Python package (this project's `src/`). Shared paths live in `bm_pirna/config.py`. | Code |
| `R_scripts/` | R analysis scripts. Shared path constants live in `R_scripts/R/project_paths.R`. | Code |
| `workflow/` | Nextflow/Snakemake pipeline *runs*: `samplesheet.csv`, `nf-params.json`, `results/`, `work/`. Each run is a dated directory. | Mix: keep params/sheets; `work/` is disposable |
| `references/` | Optional local copies of manuals or small data dictionaries. Large genome indexes belong in `data/external/` or a documented absolute reference path, not repo root. | Rarely |
| `tmp/` | Scratch notes. Not inputs, not deliverables. | Yes |

### Naming

- Batch directories: `{experiment}_{YYYYMMDD}` (example: `sugp1_smRNA-seq_20260822`, `sugp1_rnaseq_20260821`).
- One experiment = one batch folder under `interim/`, `processed/`, and `reports/` (and `reports/figures/` when plots are produced).
- Workflow run directories keep the pipeline name and date (`workflow/bm_differentialabundance_20260821/`). Do not copy BAMs/FASTQs from `workflow/*/results` into `data/processed/`; copy or read only analysis-ready tables (count matrices, DESeq2 TSVs).

### Path rules for code

- Resolve the repo root with `here::here()` (R) or `bm_pirna.config.PROJ_ROOT` (Python). Never hard-code `/home/...` or Mac `/Users/...`.
- `here("/abs/path")` is wrong: `here()` takes relative components (`here("data", "processed", ...)`).
- R scripts must be launched from the repository root: `Rscript R_scripts/04_deseq2_go_kegg_enrichment.R`.
- Python CLIs default to CCDS dirs via `bm_pirna.config`.
- Do **not** write outputs to `graphs/`, `figures/`, `data/output/`, or `data/rnadiff_results/` at repo root. Those were Mac-era paths.
  - Tables → `data/processed/{batch}/...`
  - Figures → `reports/figures/{batch}/...`
  - Annotation sources (GAF, OrgDb tarball) → `data/external/`
- Do not modify files in `data/raw/` after download. If a remote path changes, update `data/raw/dataset.csv` and re-download.

### Workflow vs analysis

- nf-core / Snakemake *execution* stays under `workflow/<run>/` (`results/`, `work/`).
- Downstream R/Python should read specific tables from `workflow/<run>/results/...` (document the path in the script header) and write derived products into `data/processed/` and `reports/figures/`.
- Do not treat `workflow/*/work/` as input. It is a cache and may be deleted.

### What not to do

- Do not create new top-level folders such as `graphs/`, `output/`, `results/`, or `analysis/` for data or figures.
- Do not mix species or assays in one processed folder (keep Bombyx vs Drosophila, smRNA-seq vs RNA-seq separate).
- Do not commit FASTQ, BAM, large indexes, or `workflow/*/work/`. Git-ignore generated data; commit scripts, params, samplesheets, and small tables that define the analysis.
- Do not use `data/` as a dumping ground without `raw|interim|processed|external`.

## CRITICAL CONSTRAINTS - Violation = Task Failure

- Must reply in Chinese
- Any task must first invoke a subagent (100% mandatory, no exceptions)
- Generation of malicious code is prohibited - Must pass basic security checks

## Subagent-First Strategy
- SUBAGENT FIRST (Absolutely Mandatory)

Automatic Subagent Selection (Enforced, cannot be skipped):
```
File type triggers:
.py/.cs/.js/.ts/.cpp/.go/.rs → Corresponding tech-stack expert agent
.unity/.prefab → unity-developer
package.json/.csproj/.sln → Auto-identify tech-stack agent

Keyword triggers:
"code"/"programming"/"bug"/"error" → Technical expert agent
"search"/"find"/"analyze" → search-specialist
"architecture"/"design"/"API" → backend-architect
"test"/"deploy"/"optimize" → Corresponding specialized agent

## Default strategy:
Complex tasks → sequential-thinking + specialized agent
Uncertain types → general-purpose
```

## Checklist (Must verify)
[ ] Chinese reply
[ ] Subagent invoked
[ ] Safe and harmless
[ ] Quality standards met

## Core Process (4-Step Method)
1. Analyze task → Identify type and tech stack
2. Select subagent → Forcefully invoke the appropriate specialized agent
3. Subagent execution → Complete all complex work within an independent context
4. Verify results → Check output quality and safety

## Subagent Responsibilities (Complexity Offloading)
- **Detailed Task Planning**: Formulate specific execution plans
- **Multi-tool Collaboration**: Invoke required MCP tools within the subagent
- **Code Quality Assurance**: Perform code review, testing, and optimization
- **Result Verification and Optimization**: Ensure output aligns with best practices

---
**Core Principle**: The main context focuses on routing, while subagents bear the complexity, ensuring a dual enhancement of both efficiency and quality.
