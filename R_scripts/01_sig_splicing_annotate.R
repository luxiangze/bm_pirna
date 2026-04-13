# Summarize significant splicing genes (FDR < 0.05) with TPM and annotation
library(dplyr)

# 1. Load rMATS summary, filter FDR < 0.05, deduplicate by GeneID
summary_sig <- read.table(
  "data/rmats_comparison_1_control_vs_sugp1/summary.txt",
  header = TRUE, sep = "\t", quote = "", comment.char = ""
) |>
  filter(FDR < 0.05) |>
  distinct(GeneID, .keep_all = TRUE)

cat("Significant genes (FDR < 0.05):", nrow(summary_sig), "\n")

# 2. Load TPM matrix and join
tpm <- read.table(
  "data/tpm/tpm_matrix.tsv",
  header = TRUE, sep = "\t", check.names = FALSE
)

# 3. Load gene annotation
gene_info <- read.table(
  "data/gene_info_table.tsv",
  header = TRUE, sep = "\t", quote = "", comment.char = ""
)

# 4. Merge: summary -> TPM -> annotation, keep informative columns only
result <- summary_sig |>
  left_join(tpm, by = c("GeneID" = "gene_id")) |>
  left_join(gene_info, by = c("GeneID" = "gene_name")) |>
  select(gene_id, gene_name = GeneID, FDR, IncLevelDifference,
         CONTROL_REP1, CONTROL_REP2, SUGP1_REP1, SUGP1_REP2,
         description, gene_biotype, chromosome, start, end, strand)

# 5. Save
dir.create("data/output", showWarnings = FALSE)
write.table(result, "data/output/sig_splicing_genes_annotated.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

cat("Output saved to data/output/sig_splicing_genes_annotated.tsv\n")
cat("Columns:", paste(colnames(result), collapse = ", "), "\n")
