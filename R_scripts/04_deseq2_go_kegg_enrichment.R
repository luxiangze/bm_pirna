# GO and KEGG enrichment for Bombyx mori DESeq2 DEGs.
# Input:  nf-core/differentialabundance filtered tables (padj < 0.05, |LFC| > 1)
#         plus data/external/gene_info_table.tsv for NCBI GeneID mapping
# Output: data/processed/sugp1_rnaseq_20260821/go_kegg/<contrast>/
#         reports/figures/sugp1_rnaseq_20260821/go_kegg/<contrast>/
# KEGG organism code: bmor
library(here)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(GO.db)
library(org.Bmori.eg.db)
library(AnnotationDbi)

source(here("R_scripts", "R", "project_paths.R"))

# ---- 0. Paths and contrasts -------------------------------------------------
DESEQ_DIR <- file.path(
  WORKFLOW_DIR, "bm_differentialabundance_20260821",
  "results", "tables", "differential", "contrasts"
)
GENE_INFO_FILE <- file.path(EXTERNAL_DATA_DIR, "gene_info_table.tsv")
ANALYSIS_ID <- "sugp1_rnaseq_20260821"

# Filtered DESeq2 tables already use padj < 0.05 and |log2FoldChange| > 1.
CONTRASTS <- list(
  control_vs_sugp1 = list(
    file = "condition_control_treated_bm_differentialabundance_202607025.deseq2.results_filtered.tsv",
    title = "SUGP1 vs Control"
  ),
  p200_vs_sugp1 = list(
    file = "condition_p200_treated_bm_differentialabundance_202607025.deseq2.results_filtered.tsv",
    title = "SUGP1 vs P200"
  )
)

PADJ_CUTOFF <- 0.2
SUB_LABEL <- "DEGs: padj<0.05 & |LFC|>1; enrichment FDR<0.2"

normalize_gene_ids <- function(x) {
  sub("\\.0$", "", as.character(x))
}

make_lookup <- function(keys, values) {
  keys <- normalize_gene_ids(keys)
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

make_go_maps <- function(gene2go_df, ont) {
  sub <- gene2go_df[gene2go_df$ONTOLOGY == ont, ]
  term2gene <- sub[, c("GO", "GID")]
  go_ids <- unique(sub$GO)
  term2name <- suppressMessages(
    AnnotationDbi::select(
      GO.db,
      keys = go_ids,
      keytype = "GOID",
      columns = "TERM"
    )[, c("GOID", "TERM")]
  )
  list(term2gene = term2gene, term2name = term2name)
}

run_go <- function(genes, gene2go_df, ont, padj_cutoff = PADJ_CUTOFF) {
  maps <- make_go_maps(gene2go_df, ont)
  enricher(
    gene = genes,
    TERM2GENE = maps$term2gene,
    TERM2NAME = maps$term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = padj_cutoff,
    qvalueCutoff = padj_cutoff
  )
}

parse_ratio <- function(x) {
  vapply(x, function(r) {
    v <- as.numeric(strsplit(r, "/")[[1]])
    v[1] / v[2]
  }, numeric(1))
}

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "grey40"),
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
    )
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
      name = "p.adjust",
      guide = guide_colorbar(barwidth = 0.6, barheight = 4)
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
      name = "p.adjust",
      guide = guide_colorbar(barwidth = 0.6, barheight = 4)
    ) +
    scale_size_continuous(range = c(3, 8), name = "Count") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = title, subtitle = subtitle, x = "Gene Ratio", y = NULL) +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))
}

plot_go_combined <- function(bp, mf, cc, title_prefix, top_n = 8) {
  prepare <- function(df, ont_label) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df |>
      slice_min(p.adjust, n = top_n) |>
      mutate(
        ont = ont_label,
        Description = stringr::str_wrap(Description, 40),
        Description = factor(Description, levels = rev(Description)),
        neg_log10 = -log10(p.adjust)
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
    labs(
      title = paste(title_prefix, "GO Enrichment Summary"),
      subtitle = SUB_LABEL,
      x = expression(-log[10](p.adjust)),
      y = NULL
    ) +
    theme_pub() +
    theme(
      strip.text = element_text(face = "bold", size = 10, color = "white"),
      strip.background = element_rect(fill = "grey40", color = NA),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
    )
}

plot_bubble <- function(df, title, subtitle = NULL, top_n = 15) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df |>
    slice_min(p.adjust, n = top_n) |>
    mutate(
      Description = stringr::str_wrap(Description, 45),
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_n = parse_ratio(GeneRatio),
      neg_log10 = -log10(p.adjust)
    )
  ggplot(df, aes(x = GeneRatio_n, y = Description,
                  size = Count, color = neg_log10)) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(
      colors = c("#3288BD", "#FDAE61", "#D53E4F"),
      name = expression(-log[10](p.adjust))
    ) +
    scale_size_continuous(range = c(4, 10), name = "Count") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = title, subtitle = subtitle, x = "Gene Ratio", y = NULL) +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))
}

plot_cnet <- function(res, title, top_n = 5) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  tryCatch(
    cnetplot(res, showCategory = top_n, node_label = "all") +
      labs(title = title) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5)),
    error = function(e) {
      cat("cnetplot failed:", conditionMessage(e), "\n")
      NULL
    }
  )
}

plot_upset <- function(res, title, top_n = 8) {
  if (is.null(res) || nrow(as.data.frame(res)) < 3) return(NULL)
  tryCatch(
    upsetplot(res, n = top_n) +
      labs(title = title) +
      theme_pub(),
    error = function(e) {
      cat("upsetplot failed:", conditionMessage(e), "\n")
      NULL
    }
  )
}

# ---- 1. Shared annotation objects ------------------------------------------
gene_info <- read.delim(
  GENE_INFO_FILE,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)
gene_info$gene_id <- normalize_gene_ids(gene_info$gene_id)
gene_info <- gene_info[!duplicated(gene_info$gene_name), ]

all_gids <- normalize_gene_ids(keys(org.Bmori.eg.db, keytype = "GID"))
gene2go <- AnnotationDbi::select(
  org.Bmori.eg.db,
  keys = keys(org.Bmori.eg.db, keytype = "GID"),
  keytype = "GID",
  columns = c("GO", "ONTOLOGY")
)
gene2go$GID <- normalize_gene_ids(gene2go$GID)
gene2go <- gene2go[!is.na(gene2go$GO), ]

run_one_contrast <- function(contrast_id, spec) {
  cat("\n========", spec$title, "(", contrast_id, ") ========\n")
  deseq_path <- file.path(DESEQ_DIR, spec$file)
  if (!file.exists(deseq_path)) {
    stop("Missing DESeq2 table: ", deseq_path)
  }

  sig_df <- read.delim(deseq_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sig_df <- sig_df |>
    left_join(
      gene_info |> dplyr::select(gene_name, ncbi_gene_id = gene_id, description),
      by = c("gene_id" = "gene_name")
    )

  n_deg <- nrow(sig_df)
  n_mapped <- sum(!is.na(sig_df$ncbi_gene_id) & nzchar(sig_df$ncbi_gene_id))
  cat("DEGs:", n_deg, " mapped NCBI GeneIDs:", n_mapped, "\n")

  sig_gids <- unique(sig_df$ncbi_gene_id[!is.na(sig_df$ncbi_gene_id) & nzchar(sig_df$ncbi_gene_id)])
  n_in_orgdb <- sum(sig_gids %in% all_gids)
  cat("GIDs in OrgDb:", n_in_orgdb, "/", length(sig_gids), "\n")

  GID_TO_SYMBOL <- make_lookup(sig_df$ncbi_gene_id, sig_df$gene_id)
  GID_TO_DESCRIPTION <- make_lookup(sig_df$ncbi_gene_id, sig_df$description)

  add_gene_columns <- function(df) {
    if (nrow(df) == 0 || !"geneID" %in% colnames(df)) return(df)
    token_rows <- lapply(df$geneID, function(s) normalize_gene_ids(split_gene_ids(s)))
    df$source_gene_ids <- vapply(token_rows, paste, character(1), collapse = "/")
    df$gene_symbols <- vapply(
      token_rows, collapse_lookup, character(1),
      lookup = GID_TO_SYMBOL
    )
    df$gene_descriptions <- vapply(
      token_rows, collapse_lookup, character(1),
      lookup = GID_TO_DESCRIPTION
    )
    df
  }

  go_bp <- run_go(sig_gids, gene2go, "BP")
  go_mf <- run_go(sig_gids, gene2go, "MF")
  go_cc <- run_go(sig_gids, gene2go, "CC")

  kegg_res <- tryCatch(
    enrichKEGG(
      gene = sig_gids,
      organism = "bmor",
      pAdjustMethod = "BH",
      pvalueCutoff = PADJ_CUTOFF,
      qvalueCutoff = PADJ_CUTOFF
    ),
    error = function(e) {
      cat("KEGG enrichment failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  out_data <- file.path(PROCESSED_DATA_DIR, ANALYSIS_ID, "go_kegg", contrast_id)
  out_fig <- file.path(FIGURES_DIR, ANALYSIS_ID, "go_kegg", contrast_id)
  dir.create(out_data, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

  save_tsv <- function(res, name) {
    df <- if (is.null(res)) data.frame() else as.data.frame(res)
    df <- add_gene_columns(df)
    if (nrow(df) > 0) {
      write.table(
        df, file.path(out_data, paste0(name, ".tsv")),
        sep = "\t", row.names = FALSE, quote = FALSE
      )
    }
    cat(name, ":", nrow(df), "terms\n")
    invisible(df)
  }

  bp_df <- save_tsv(go_bp, "GO_BP")
  mf_df <- save_tsv(go_mf, "GO_MF")
  cc_df <- save_tsv(go_cc, "GO_CC")
  kegg_df <- save_tsv(kegg_res, "KEGG")

  deg_anno <- sig_df |>
    dplyr::select(
      gene_id, ncbi_gene_id, description,
      baseMean, log2FoldChange, lfcSE, pvalue, padj
    )
  write.table(
    deg_anno, file.path(out_data, "DEG_annotated.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  cat("DEG_annotated:", nrow(deg_anno), "genes\n")

  save_plot <- function(p, name, w = 8, h = 6) {
    if (is.null(p)) {
      cat("Skipped:", name, "(no data)\n")
      return(invisible(NULL))
    }
    ggsave(file.path(out_fig, paste0(name, ".pdf")), p,
           width = w, height = h, useDingbats = FALSE)
    ggsave(file.path(out_fig, paste0(name, ".png")), p,
           width = w, height = h, dpi = 300)
    cat("Saved:", name, "\n")
  }

  prefix <- spec$title
  save_plot(plot_lollipop(bp_df, paste(prefix, "GO Biological Process"), SUB_LABEL),
            "GO_BP_lollipop", 9, 5)
  save_plot(plot_lollipop(mf_df, paste(prefix, "GO Molecular Function"), SUB_LABEL),
            "GO_MF_lollipop", 9, 5)
  save_plot(plot_lollipop(cc_df, paste(prefix, "GO Cellular Component"), SUB_LABEL),
            "GO_CC_lollipop", 9, 5)
  save_plot(plot_lollipop(kegg_df, paste(prefix, "KEGG Pathway Enrichment"), SUB_LABEL, 20),
            "KEGG_lollipop", 10, 8)

  save_plot(plot_dot(bp_df, paste(prefix, "GO Biological Process"), SUB_LABEL),
            "GO_BP_dotplot", 9, 5)
  save_plot(plot_dot(mf_df, paste(prefix, "GO Molecular Function"), SUB_LABEL),
            "GO_MF_dotplot", 9, 5)
  save_plot(plot_dot(cc_df, paste(prefix, "GO Cellular Component"), SUB_LABEL),
            "GO_CC_dotplot", 9, 5)
  save_plot(plot_dot(kegg_df, paste(prefix, "KEGG Pathway Enrichment"), SUB_LABEL, 20),
            "KEGG_dotplot", 10, 8)

  save_plot(plot_bubble(mf_df, paste(prefix, "GO Molecular Function"), SUB_LABEL),
            "GO_MF_bubble", 9, 5)
  save_plot(plot_bubble(kegg_df, paste(prefix, "KEGG Pathway Enrichment"), SUB_LABEL, 20),
            "KEGG_bubble", 10, 8)

  save_plot(plot_go_combined(bp_df, mf_df, cc_df, prefix, top_n = 8),
            "GO_combined_barplot", 9, 7)

  save_plot(plot_cnet(go_bp, paste(prefix, "GO BP — Gene-Concept Network")),
            "GO_BP_cnetplot", 10, 9)
  save_plot(plot_cnet(go_mf, paste(prefix, "GO MF — Gene-Concept Network")),
            "GO_MF_cnetplot", 10, 9)
  save_plot(plot_cnet(kegg_res, paste(prefix, "KEGG — Gene-Concept Network")),
            "KEGG_cnetplot", 10, 9)

  save_plot(plot_upset(go_bp, paste(prefix, "GO BP — Term Overlap")),
            "GO_BP_upset", 10, 6)
  save_plot(plot_upset(go_mf, paste(prefix, "GO MF — Term Overlap")),
            "GO_MF_upset", 10, 6)
  save_plot(plot_upset(kegg_res, paste(prefix, "KEGG — Term Overlap")),
            "KEGG_upset", 10, 6)

  cat("Wrote tables:", out_data, "\n")
  cat("Wrote figures:", out_fig, "\n")
}

invisible(lapply(names(CONTRASTS), function(id) {
  run_one_contrast(id, CONTRASTS[[id]])
}))

cat(
  "\nDone. Results in",
  file.path("data/processed", ANALYSIS_ID, "go_kegg"),
  "and",
  file.path("reports/figures", ANALYSIS_ID, "go_kegg"),
  "\n"
)
