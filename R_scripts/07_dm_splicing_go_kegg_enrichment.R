# GO and KEGG enrichment analysis for Drosophila significant alternative splicing genes
# Filter: FDR < 0.05 (rMATS JC mode), one row per GeneID
# OrgDb: org.Dm.eg.db (input keytype: FLYBASE; KEGG organism: dme)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(GO.db)
library(org.Dm.eg.db)
library(AnnotationDbi)

# ---- 1. Load and filter splicing summary ----
splicing_df <- read.table(
  "data/dm_splicing_20260425/summary.JC.tsv",
  header = TRUE, sep = "\t", quote = "", comment.char = ""
)

sig_df <- splicing_df |>
  filter(!is.na(FDR), FDR < 0.05) |>
  distinct(GeneID, .keep_all = TRUE)

sig_fbgn <- unique(sig_df$GeneID)
cat("Significant splicing genes (FDR<0.05):", length(sig_fbgn), "\n")

# ---- 2. ID mapping helpers ----
# enrichKEGG(organism="dme") expects KEGG-style IDs like "Dmel_CG11348".
fbgn_to_kegg <- function(fbgn) {
  cg <- suppressMessages(
    AnnotationDbi::mapIds(org.Dm.eg.db,
                          keys      = fbgn,
                          keytype   = "FLYBASE",
                          column    = "FLYBASECG",
                          multiVals = "first")
  )
  cg <- cg[!is.na(cg)]
  paste0("Dmel_", unname(cg))
}

.cg_keys <- keys(org.Dm.eg.db, keytype = "FLYBASECG")
.cg_sym  <- suppressMessages(
  AnnotationDbi::mapIds(org.Dm.eg.db,
                        keys      = .cg_keys,
                        keytype   = "FLYBASECG",
                        column    = "SYMBOL",
                        multiVals = "first")
)
KEGG_ID_TO_SYMBOL <- setNames(unname(.cg_sym), paste0("Dmel_", .cg_keys))

kegg_geneID_to_symbol <- function(s) {
  if (is.na(s) || !nzchar(s)) return(s)
  parts <- strsplit(s, "/", fixed = TRUE)[[1]]
  mapped <- ifelse(parts %in% names(KEGG_ID_TO_SYMBOL) &
                     !is.na(KEGG_ID_TO_SYMBOL[parts]),
                   KEGG_ID_TO_SYMBOL[parts], parts)
  paste(mapped, collapse = "/")
}

make_lookup <- function(keys, values) {
  keys <- as.character(keys)
  values <- as.character(values)
  keep <- !is.na(keys) & nzchar(keys) & !duplicated(keys)
  setNames(values[keep], keys[keep])
}

split_gene_ids <- function(s) {
  if (is.na(s) || !nzchar(s)) return(character())
  trimws(strsplit(as.character(s), "/", fixed = TRUE)[[1]])
}

collapse_lookup <- function(tokens, lookup, fallback = tokens) {
  if (length(tokens) == 0) return("")
  values <- unname(lookup[tokens])
  missing <- is.na(values) | !nzchar(values)
  values[missing] <- fallback[missing]
  paste(values, collapse = "/")
}

make_dm_gene_lookup <- function(fbgn) {
  raw <- suppressMessages(
    AnnotationDbi::select(org.Dm.eg.db,
                          keys    = unique(fbgn),
                          keytype = "FLYBASE",
                          columns = c("SYMBOL", "FLYBASECG", "GENENAME"))
  )
  cg_tokens <- ifelse(is.na(raw$FLYBASECG) | !nzchar(raw$FLYBASECG),
                      NA_character_, paste0("Dmel_", raw$FLYBASECG))
  lookup_keys <- c(raw$SYMBOL, cg_tokens)
  list(
    source_ids = make_lookup(lookup_keys, c(raw$FLYBASE, raw$FLYBASE)),
    symbols = make_lookup(lookup_keys, c(raw$SYMBOL, raw$SYMBOL)),
    descriptions = make_lookup(lookup_keys, c(raw$GENENAME, raw$GENENAME))
  )
}

add_gene_columns <- function(df, gene_lookup) {
  if (nrow(df) == 0 || !"geneID" %in% colnames(df)) return(df)
  token_rows <- lapply(df$geneID, split_gene_ids)
  df$source_gene_ids <- vapply(token_rows, collapse_lookup, character(1),
                               lookup = gene_lookup$source_ids)
  df$gene_symbols <- vapply(token_rows, collapse_lookup, character(1),
                            lookup = gene_lookup$symbols)
  df$gene_descriptions <- vapply(token_rows, collapse_lookup, character(1),
                                 lookup = gene_lookup$descriptions)
  df
}

# ---- 3. GO enrichment (FLYBASE keytype) ----
run_go <- function(fbgn, ont) {
  enrichGO(
    gene          = fbgn,
    OrgDb         = org.Dm.eg.db,
    keyType       = "FLYBASE",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
}

go_bp <- run_go(sig_fbgn, "BP")
go_mf <- run_go(sig_fbgn, "MF")
go_cc <- run_go(sig_fbgn, "CC")

# ---- 3b. KEGG enrichment (organism = "dme", needs Dmel_CGxxx IDs) ----
kegg_ids <- fbgn_to_kegg(sig_fbgn)
cat("Mapped to KEGG (Dmel_CG) IDs:", length(kegg_ids), "\n")

kegg_res <- enrichKEGG(
  gene          = kegg_ids,
  organism      = "dme",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  kegg_res@result$geneID <- vapply(kegg_res@result$geneID,
                                   kegg_geneID_to_symbol, character(1))
}
gene_lookup <- make_dm_gene_lookup(sig_fbgn)

# ---- 4. Save TSV ----
dir.create("data/output/dm_enrichment_splicing", showWarnings = FALSE, recursive = TRUE)
dir.create("graphs/dm_enrichment_splicing",      showWarnings = FALSE, recursive = TRUE)

save_tsv <- function(res, name) {
  df <- if (is.null(res)) data.frame() else as.data.frame(res)
  df <- add_gene_columns(df, gene_lookup)
  if (nrow(df) > 0) {
    write.table(df, paste0("data/output/dm_enrichment_splicing/", name, ".tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  cat(name, ":", nrow(df), "terms\n")
  invisible(df)
}

bp_df   <- save_tsv(go_bp,    "GO_BP")
mf_df   <- save_tsv(go_mf,    "GO_MF")
cc_df   <- save_tsv(go_cc,    "GO_CC")
kegg_df <- save_tsv(kegg_res, "KEGG")

# ============================================================
# Publication-quality plots
# ============================================================

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(color = "black"),
      axis.title         = element_text(color = "black"),
      plot.title         = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle      = element_text(hjust = 0.5, size = base_size - 1, color = "grey40"),
      legend.key.size    = unit(0.4, "cm"),
      legend.title       = element_text(size = 9),
      legend.text        = element_text(size = 8),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
    )
}

parse_ratio <- function(x) {
  vapply(x, function(r) {
    v <- as.numeric(strsplit(r, "/")[[1]]); v[1] / v[2]
  }, numeric(1))
}

save_plot <- function(p, name, w = 8, h = 6) {
  if (is.null(p)) { cat("Skipped:", name, "(no data)\n"); return(invisible(NULL)) }
  ggsave(paste0("graphs/dm_enrichment_splicing/", name, ".pdf"), p,
         width = w, height = h, useDingbats = FALSE)
  ggsave(paste0("graphs/dm_enrichment_splicing/", name, ".png"), p,
         width = w, height = h, dpi = 300)
  cat("Saved:", name, "\n")
}

plot_lollipop <- function(df, title, subtitle = NULL, top_n = 15) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df |>
    slice_min(p.adjust, n = top_n) |>
    mutate(
      Description = stringr::str_wrap(Description, 45),
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_n = parse_ratio(GeneRatio)
    )
  ggplot(df, aes(x = GeneRatio_n, y = Description, color = p.adjust)) +
    geom_segment(aes(x = 0, xend = GeneRatio_n, yend = Description),
                 color = "grey75", linewidth = 0.7) +
    geom_point(aes(size = Count)) +
    scale_color_gradientn(
      colors = rev(brewer.pal(9, "RdYlBu")),
      name   = "p.adjust",
      guide  = guide_colorbar(barwidth = 0.6, barheight = 4)
    ) +
    scale_size_continuous(range = c(3, 8), name = "Count") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    labs(title = title, subtitle = subtitle, x = "Gene Ratio", y = NULL) +
    theme_pub()
}

plot_dot <- function(df, title, subtitle = NULL, top_n = 15) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df |>
    slice_min(p.adjust, n = top_n) |>
    mutate(
      Description = stringr::str_wrap(Description, 45),
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_n = parse_ratio(GeneRatio)
    )
  ggplot(df, aes(x = GeneRatio_n, y = Description,
                 color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradientn(
      colors = rev(brewer.pal(9, "RdYlBu")),
      name   = "p.adjust",
      guide  = guide_colorbar(barwidth = 0.6, barheight = 4)
    ) +
    scale_size_continuous(range = c(3, 8), name = "Count") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = title, subtitle = subtitle, x = "Gene Ratio", y = NULL) +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))
}

plot_bubble <- function(df, title, subtitle = NULL, top_n = 15) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df |>
    slice_min(p.adjust, n = top_n) |>
    mutate(
      Description = stringr::str_wrap(Description, 45),
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_n = parse_ratio(GeneRatio),
      neg_log10   = -log10(p.adjust)
    )
  ggplot(df, aes(x = GeneRatio_n, y = Description,
                 size = Count, color = neg_log10)) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(
      colors = c("#3288BD", "#FDAE61", "#D53E4F"),
      name   = expression(-log[10](p.adjust))
    ) +
    scale_size_continuous(range = c(4, 10), name = "Count") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = title, subtitle = subtitle, x = "Gene Ratio", y = NULL) +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))
}

plot_go_combined <- function(bp, mf, cc, subtitle = "FDR < 0.05",
                             top_n = 8) {
  prepare <- function(df, ont_label) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df |>
      slice_min(p.adjust, n = top_n) |>
      mutate(
        ont         = ont_label,
        Description = stringr::str_wrap(Description, 40),
        Description = factor(Description, levels = rev(Description)),
        neg_log10   = -log10(p.adjust)
      )
  }
  combined <- bind_rows(prepare(bp, "BP"), prepare(mf, "MF"), prepare(cc, "CC"))
  if (is.null(combined) || nrow(combined) == 0) return(NULL)

  pal <- c(BP = "#4393C3", MF = "#E6550D", CC = "#74C476")
  combined$ont <- factor(combined$ont, levels = c("BP", "MF", "CC"))

  ggplot(combined, aes(x = neg_log10, y = Description, fill = ont)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = paste0("n=", Count)),
              hjust = -0.1, size = 3, color = "grey30") +
    facet_grid(ont ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = pal) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "GO Enrichment Summary",
         subtitle = subtitle,
         x = expression(-log[10](p.adjust)), y = NULL) +
    theme_pub() +
    theme(
      strip.text         = element_text(face = "bold", size = 10, color = "white"),
      strip.background   = element_rect(fill = "grey40", color = NA),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
    )
}

plot_cnet <- function(res, title, top_n = 5) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  tryCatch(
    cnetplot(res, showCategory = top_n, node_label = "all") +
      labs(title = title) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5)),
    error = function(e) { cat("cnetplot failed:", conditionMessage(e), "\n"); NULL }
  )
}

plot_upset <- function(res, title, top_n = 8) {
  if (is.null(res) || nrow(as.data.frame(res)) < 3) return(NULL)
  tryCatch(
    upsetplot(res, n = top_n) +
      labs(title = title) +
      theme_pub(),
    error = function(e) { cat("upsetplot failed:", conditionMessage(e), "\n"); NULL }
  )
}

# ---- 5. Generate all plots ----
sub_label <- "Significant AS genes (FDR<0.05, JC mode)"

# Bubble (primary)
save_plot(plot_bubble(bp_df,   "GO Biological Process",   sub_label),     "GO_BP_bubble",  9, 6)
save_plot(plot_bubble(mf_df,   "GO Molecular Function",   sub_label),     "GO_MF_bubble",  9, 6)
save_plot(plot_bubble(cc_df,   "GO Cellular Component",   sub_label),     "GO_CC_bubble",  9, 6)
save_plot(plot_bubble(kegg_df, "KEGG Pathway Enrichment", sub_label, 20), "KEGG_bubble",  10, 8)

# Lollipop
save_plot(plot_lollipop(bp_df,   "GO Biological Process",   sub_label),     "GO_BP_lollipop",  9, 5)
save_plot(plot_lollipop(mf_df,   "GO Molecular Function",   sub_label),     "GO_MF_lollipop",  9, 5)
save_plot(plot_lollipop(cc_df,   "GO Cellular Component",   sub_label),     "GO_CC_lollipop",  9, 5)
save_plot(plot_lollipop(kegg_df, "KEGG Pathway Enrichment", sub_label, 20), "KEGG_lollipop", 10, 8)

# Dot
save_plot(plot_dot(bp_df,   "GO Biological Process",   sub_label),     "GO_BP_dotplot",  9, 5)
save_plot(plot_dot(mf_df,   "GO Molecular Function",   sub_label),     "GO_MF_dotplot",  9, 5)
save_plot(plot_dot(cc_df,   "GO Cellular Component",   sub_label),     "GO_CC_dotplot",  9, 5)
save_plot(plot_dot(kegg_df, "KEGG Pathway Enrichment", sub_label, 20), "KEGG_dotplot",  10, 8)

# Combined GO
save_plot(plot_go_combined(bp_df, mf_df, cc_df, sub_label, top_n = 8),
          "GO_combined_barplot", 9, 7)

# Gene-concept network
save_plot(plot_cnet(go_bp,    "GO BP — Gene-Concept Network"),  "GO_BP_cnetplot",   10, 9)
save_plot(plot_cnet(go_mf,    "GO MF — Gene-Concept Network"),  "GO_MF_cnetplot",   10, 9)
save_plot(plot_cnet(kegg_res, "KEGG — Gene-Concept Network"),   "KEGG_cnetplot",    10, 9)

# Upset
save_plot(plot_upset(go_bp,    "GO BP — Term Overlap"), "GO_BP_upset", 10, 6)
save_plot(plot_upset(go_mf,    "GO MF — Term Overlap"), "GO_MF_upset", 10, 6)
save_plot(plot_upset(kegg_res, "KEGG — Term Overlap"),  "KEGG_upset",  10, 6)

cat("\nDone. Results in data/output/dm_enrichment_splicing/ and graphs/dm_enrichment_splicing/\n")
