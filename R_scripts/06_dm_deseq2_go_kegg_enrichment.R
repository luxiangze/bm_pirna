# GO and KEGG enrichment analysis for Drosophila DESeq2 differentially expressed genes
# Three gene sets: all DEGs, up-regulated, down-regulated
# Filter: padj < 0.05, |log2FoldChange| > 1
# OrgDb: org.Dm.eg.db (input keytype: FLYBASE; KEGG organism: dme)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(GO.db)
library(org.Dm.eg.db)
library(AnnotationDbi)

# ---- 1. Load and split DESeq2 results ----
deseq_df <- read.table(
  "data/dm_differential/condition_control_treated.deseq2.results.tsv",
  header = TRUE, sep = "\t"
)

sig_all <- deseq_df |>
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
sig_up   <- sig_all |> filter(log2FoldChange >  1)
sig_down <- sig_all |> filter(log2FoldChange < -1)

cat("All DEGs:",  nrow(sig_all),  "\n")
cat("Up-reg:",    nrow(sig_up),   "\n")
cat("Down-reg:",  nrow(sig_down), "\n")

# ---- 2. ID mapping helpers ----
# enrichKEGG(organism="dme") expects KEGG-style IDs like "Dmel_CG11348",
# i.e. FlyBase CG number prefixed with "Dmel_". We also prebuild a
# CG -> SYMBOL lookup so KEGG geneID columns can be made readable.
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

# Build a Dmel_CG -> SYMBOL lookup (full OrgDb, used for KEGG readability)
.cg_keys <- keys(org.Dm.eg.db, keytype = "FLYBASECG")
.cg_sym  <- suppressMessages(
  AnnotationDbi::mapIds(org.Dm.eg.db,
                        keys      = .cg_keys,
                        keytype   = "FLYBASECG",
                        column    = "SYMBOL",
                        multiVals = "first")
)
KEGG_ID_TO_SYMBOL <- setNames(unname(.cg_sym), paste0("Dmel_", .cg_keys))

# Replace "Dmel_CGxxx/Dmel_CGyyy" in a geneID string with "sym1/sym2"
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

# ---- 3. Enrichment runners ----
run_go <- function(fbgn, ont) {
  if (length(fbgn) == 0) return(NULL)
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

run_kegg <- function(fbgn) {
  kegg_ids <- fbgn_to_kegg(fbgn)
  if (length(kegg_ids) == 0) return(NULL)
  res <- enrichKEGG(
    gene          = kegg_ids,
    organism      = "dme",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(res)
  # Make geneID column human-readable (Dmel_CGxxx -> SYMBOL)
  res@result$geneID <- vapply(res@result$geneID, kegg_geneID_to_symbol,
                              character(1))
  res
}

# ============================================================
# Publication-quality plot helpers (ported from munge/03 / 04)
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

plot_go_combined <- function(bp, mf, cc, subtitle = "padj<0.05, |LFC|>1",
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

# ---- 4. Per-group runner ----
# group_label: "all" / "up" / "down"; subtitle distinguishes them on the plots
run_one_group <- function(sig_df, group_label, subtitle) {
  if (nrow(sig_df) == 0) {
    cat("\n[", group_label, "] No genes — skipped\n")
    return(invisible(NULL))
  }

  out_data <- file.path("data/output/dm_enrichment_deseq2", group_label)
  out_fig  <- "graphs/dm_enrichment_deseq2"
  dir.create(out_data, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig,  showWarnings = FALSE, recursive = TRUE)

  fbgn <- unique(sig_df$gene_id)
  cat("\n[", group_label, "] FBgn input:", length(fbgn), "\n")
  gene_lookup <- make_dm_gene_lookup(fbgn)

  go_bp <- run_go(fbgn, "BP")
  go_mf <- run_go(fbgn, "MF")
  go_cc <- run_go(fbgn, "CC")
  kegg  <- run_kegg(fbgn)

  save_tsv <- function(res, name) {
    df <- if (is.null(res)) data.frame() else as.data.frame(res)
    df <- add_gene_columns(df, gene_lookup)
    if (nrow(df) > 0) {
      write.table(df, file.path(out_data, paste0(name, ".tsv")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
    cat("  ", name, ":", nrow(df), "terms\n")
    invisible(df)
  }
  bp_df   <- save_tsv(go_bp, "GO_BP")
  mf_df   <- save_tsv(go_mf, "GO_MF")
  cc_df   <- save_tsv(go_cc, "GO_CC")
  kegg_df <- save_tsv(kegg,  "KEGG")

  save_plot <- function(p, name, w = 8, h = 6) {
    if (is.null(p)) { cat("  Skipped:", name, "(no data)\n"); return(invisible(NULL)) }
    ggsave(file.path(out_fig, paste0(group_label, "_", name, ".pdf")), p,
           width = w, height = h, useDingbats = FALSE)
    ggsave(file.path(out_fig, paste0(group_label, "_", name, ".png")), p,
           width = w, height = h, dpi = 300)
    cat("  Saved:", group_label, "_", name, "\n", sep = "")
  }

  # Bubble (primary, requested by user)
  save_plot(plot_bubble(bp_df,   "GO Biological Process",   subtitle),     "GO_BP_bubble",  9, 6)
  save_plot(plot_bubble(mf_df,   "GO Molecular Function",   subtitle),     "GO_MF_bubble",  9, 6)
  save_plot(plot_bubble(cc_df,   "GO Cellular Component",   subtitle),     "GO_CC_bubble",  9, 6)
  save_plot(plot_bubble(kegg_df, "KEGG Pathway Enrichment", subtitle, 20), "KEGG_bubble",  10, 8)

  # Lollipop
  save_plot(plot_lollipop(bp_df,   "GO Biological Process",   subtitle),     "GO_BP_lollipop",  9, 5)
  save_plot(plot_lollipop(mf_df,   "GO Molecular Function",   subtitle),     "GO_MF_lollipop",  9, 5)
  save_plot(plot_lollipop(cc_df,   "GO Cellular Component",   subtitle),     "GO_CC_lollipop",  9, 5)
  save_plot(plot_lollipop(kegg_df, "KEGG Pathway Enrichment", subtitle, 20), "KEGG_lollipop", 10, 8)

  # Dot
  save_plot(plot_dot(bp_df,   "GO Biological Process",   subtitle),     "GO_BP_dotplot",  9, 5)
  save_plot(plot_dot(mf_df,   "GO Molecular Function",   subtitle),     "GO_MF_dotplot",  9, 5)
  save_plot(plot_dot(cc_df,   "GO Cellular Component",   subtitle),     "GO_CC_dotplot",  9, 5)
  save_plot(plot_dot(kegg_df, "KEGG Pathway Enrichment", subtitle, 20), "KEGG_dotplot",  10, 8)

  # Combined GO
  save_plot(plot_go_combined(bp_df, mf_df, cc_df, subtitle, top_n = 8),
            "GO_combined_barplot", 9, 7)

  # Gene-concept network
  save_plot(plot_cnet(go_bp, "GO BP — Gene-Concept Network"), "GO_BP_cnetplot", 10, 9)
  save_plot(plot_cnet(go_mf, "GO MF — Gene-Concept Network"), "GO_MF_cnetplot", 10, 9)
  save_plot(plot_cnet(kegg,  "KEGG — Gene-Concept Network"),  "KEGG_cnetplot",  10, 9)

  # Upset
  save_plot(plot_upset(go_bp, "GO BP — Term Overlap"), "GO_BP_upset", 10, 6)
  save_plot(plot_upset(go_mf, "GO MF — Term Overlap"), "GO_MF_upset", 10, 6)
  save_plot(plot_upset(kegg,  "KEGG — Term Overlap"),  "KEGG_upset",  10, 6)
}

# ---- 5. Run for each group ----
run_one_group(sig_all,  "all",  "All DEGs (padj<0.05, |LFC|>1)")
run_one_group(sig_up,   "up",   "Up-regulated DEGs (padj<0.05, LFC>1)")
run_one_group(sig_down, "down", "Down-regulated DEGs (padj<0.05, LFC<-1)")

cat("\nDone. Results in data/output/dm_enrichment_deseq2/{all,up,down}/ and graphs/dm_enrichment_deseq2/\n")
