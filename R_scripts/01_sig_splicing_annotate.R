#' Annotate significant splicing genes with TPM and GTF metadata
#'
#' Filter rMATS summary events at FDR < 0.05, keep the most significant
#' event per gene, join TPM from the splicing workflow and gene annotation
#' extracted from the NCBI GTF.
#'
#' Run from the repository root:
#'   Rscript R_scripts/01_sig_splicing_annotate.R
#'
#' @author Yongkang Guo

library(dplyr)
library(readr)

FDR_CUTOFF <- 0.05

SPLICING_DIR <- "workflow/rnasplicing_bm_ncbi_20260824"
TPM_FILE <- file.path(SPLICING_DIR, "results/tpm/tpm_matrix.tsv")
GENE_INFO_FILE <- "data/external/gene_info_table.tsv"
OUTPUT_DIR <- file.path("data/processed", basename(SPLICING_DIR))

COMPARISONS <- list(
  control_vs_sugp1 = file.path(
    SPLICING_DIR, "results/comparison_1_control_vs_sugp1/summary.txt"
  ),
  p200_vs_sugp1 = file.path(
    SPLICING_DIR, "results/comparison_2_p200_vs_sugp1/summary.txt"
  )
)

#' Load rMATS summary, keep significant events, one row per gene
#'
#' @param summary_file Path to summarize_splicing.py summary.txt
#' @param comparison Comparison label
#' @return Tibble of significant genes
load_sig_summary <- function(summary_file, comparison) {
  if (!file.exists(summary_file)) {
    stop("rMATS summary not found: ", summary_file)
  }

  summary_tbl <- read_tsv(
    summary_file,
    col_types = cols(
      ID = col_character(),
      GeneID = col_character(),
      geneSymbol = col_character(),
      PValue = col_double(),
      FDR = col_double(),
      IncLevel1 = col_character(),
      IncLevel2 = col_character(),
      IncLevelDifference = col_double(),
      EventType = col_character()
    ),
    quote = ""
  )

  sig <- summary_tbl |>
    filter(!is.na(GeneID), GeneID != "", !is.na(FDR), FDR < FDR_CUTOFF) |>
    arrange(FDR, desc(abs(IncLevelDifference))) |>
    distinct(GeneID, .keep_all = TRUE) |>
    mutate(comparison = comparison)

  message(sprintf(
    "%s: %d events, %d genes with FDR < %s",
    comparison, nrow(summary_tbl), nrow(sig), FDR_CUTOFF
  ))
  sig
}

#' Annotate one comparison with TPM and gene info
#'
#' @param sig Significant splicing table from load_sig_summary()
#' @param tpm TPM matrix (gene_id + sample columns)
#' @param gene_info GTF gene table from extract_gene_info.py
#' @return Annotated tibble
annotate_sig <- function(sig, tpm, gene_info) {
  sample_cols <- setdiff(names(tpm), "gene_id")

  result <- sig |>
    left_join(tpm, by = c("GeneID" = "gene_id")) |>
    left_join(gene_info, by = c("GeneID" = "gene_name")) |>
    transmute(
      comparison,
      gene_name = GeneID,
      gene_id = gene_id,
      gene_symbol = na_if(geneSymbol, ""),
      event_type = EventType,
      pvalue = PValue,
      fdr = FDR,
      inc_level_1 = IncLevel1,
      inc_level_2 = IncLevel2,
      inc_level_difference = IncLevelDifference,
      across(all_of(sample_cols)),
      description,
      gene_biotype,
      chromosome,
      start,
      end,
      strand,
      gene_synonym
    )

  n_tpm <- sum(!is.na(result[[sample_cols[1]]]))
  n_ann <- sum(!is.na(result$gene_id) & result$gene_id != "")
  message(sprintf(
    "  joined TPM for %d/%d genes, annotation for %d/%d genes",
    n_tpm, nrow(result), n_ann, nrow(result)
  ))
  result
}

main <- function() {
  if (!file.exists(TPM_FILE)) {
    stop("TPM matrix not found: ", TPM_FILE)
  }
  if (!file.exists(GENE_INFO_FILE)) {
    stop(
      "Gene annotation not found: ", GENE_INFO_FILE,
      "\nGenerate it from the NCBI GTF with bm_pirna.extract_gene_info"
    )
  }

  tpm <- read_tsv(TPM_FILE, show_col_types = FALSE)
  if (!"gene_id" %in% names(tpm)) {
    stop("TPM matrix must have a gene_id column")
  }

  gene_info <- read_tsv(
    GENE_INFO_FILE,
    col_types = cols(
      gene_name = col_character(),
      gene_id = col_character(),
      gene_synonym = col_character(),
      description = col_character(),
      gbkey = col_character(),
      gene_biotype = col_character(),
      chromosome = col_character(),
      start = col_double(),
      end = col_double(),
      strand = col_character()
    ),
    quote = ""
  )

  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

  annotated <- lapply(names(COMPARISONS), function(name) {
    sig <- load_sig_summary(COMPARISONS[[name]], name)
    annotate_sig(sig, tpm, gene_info)
  })
  names(annotated) <- names(COMPARISONS)

  for (name in names(annotated)) {
    out_file <- file.path(
      OUTPUT_DIR, sprintf("sig_splicing_genes_annotated_%s.tsv", name)
    )
    write_tsv(annotated[[name]], out_file)
    message("Wrote ", out_file, " (", nrow(annotated[[name]]), " genes)")
  }

  combined <- bind_rows(annotated)
  combined_file <- file.path(OUTPUT_DIR, "sig_splicing_genes_annotated.tsv")
  write_tsv(combined, combined_file)
  message("Wrote ", combined_file, " (", nrow(combined), " rows)")
  message("Columns: ", paste(colnames(combined), collapse = ", "))
}

main()
