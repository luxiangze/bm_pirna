# GO and KEGG enrichment analysis for DESeq2 differentially expressed genes
# Filter: pvalue < 0.05, |log2FoldChange| > 1
# KEGG organism code: bmor (Bombyx mori)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(GO.db)
library(org.Bmori.eg.db)

# ---- 1. Load and filter DESeq2 results ----
deseq_df <- read.table(
  "data/rnadiff_results/tables/differential/condition_control_treated.deseq2.results.tsv",
  header = TRUE, sep = "\t"
)

sig_df <- deseq_df |>
  filter(!is.na(pvalue), pvalue < 0.05, abs(log2FoldChange) > 1)

cat("DEGs (pvalue<0.05, |lfc|>1):", nrow(sig_df), "\n")

# Convert SYMBOL -> GID (required by org.Bmori.eg.db and enrichKEGG)
id_map <- bitr(sig_df$gene_id,
               fromType = "Symbol",
               toType   = "GID",
               OrgDb    = org.Bmori.eg.db)
sig_gids <- unique(id_map$GID)
cat("Mapped GIDs:", length(sig_gids), "\n")

# ---- 2. Build gene2GO mapping from OrgDb ----
all_gids <- keys(org.Bmori.eg.db, keytype = "GID")
gene2go  <- select(org.Bmori.eg.db,
                   keys    = all_gids,
                   keytype = "GID",
                   columns = c("GO", "ONTOLOGY"))
gene2go  <- gene2go[!is.na(gene2go$GO), ]

# Build term2gene + term2name for a given ontology
make_go_maps <- function(gene2go_df, ont) {
  sub       <- gene2go_df[gene2go_df$ONTOLOGY == ont, ]
  term2gene <- sub[, c("GO", "GID")]
  go_ids    <- unique(sub$GO)
  term2name <- suppressMessages(
    AnnotationDbi::select(GO.db,
                          keys    = go_ids,
                          keytype = "GOID",
                          columns = "TERM")[, c("GOID", "TERM")]
  )
  list(term2gene = term2gene, term2name = term2name)
}

# ---- 3. GO enrichment ----
run_go <- function(genes, gene2go_df, ont, padj_cutoff = 0.2) {
  maps <- make_go_maps(gene2go_df, ont)
  enricher(
    gene          = genes,
    TERM2GENE     = maps$term2gene,
    TERM2NAME     = maps$term2name,
    pAdjustMethod = "BH",
    pvalueCutoff  = padj_cutoff,
    qvalueCutoff  = padj_cutoff
  )
}

go_bp <- run_go(sig_gids, gene2go, "BP")
go_mf <- run_go(sig_gids, gene2go, "MF")
go_cc <- run_go(sig_gids, gene2go, "CC")

# ---- 3b. KEGG enrichment ----
# enrichKEGG requires bare NCBI Gene IDs; strip "LOC" prefix from gene symbols
kegg_ids <- unique(sub("^LOC", "", sig_df$gene_id))
kegg_res <- enrichKEGG(
  gene          = kegg_ids,
  organism      = "bmor",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2
)

# ---- 4. Annotated DEG table ----
# Step 1: query GENENAME from OrgDb for genes with a valid Symbol mapping
orgdb_sym  <- keys(org.Bmori.eg.db, keytype = "Symbol")
has_symbol <- sig_df$gene_id %in% orgdb_sym

orgdb_anno <- data.frame(gene_id = character(), description = character(),
                         stringsAsFactors = FALSE)
if (any(has_symbol)) {
  raw <- suppressMessages(
    AnnotationDbi::select(org.Bmori.eg.db,
                          keys    = sig_df$gene_id[has_symbol],
                          keytype = "Symbol",
                          columns = "GENENAME")
  )
  # keep first hit per gene
  orgdb_anno <- raw[!duplicated(raw$Symbol), c("Symbol", "GENENAME")]
  colnames(orgdb_anno) <- c("gene_id", "description")
}

# Step 2: for genes missing from OrgDb, pull 'product' from the annotation file
anno_file <- read.table(
  "data/rnadiff_results/tables/annotation/genomic_chr.anno.tsv",
  header = TRUE, sep = "\t", quote = "", fill = TRUE
)
# deduplicate: keep one row per gene_id
anno_file <- anno_file[!duplicated(anno_file$gene_id), ]

missing_ids  <- sig_df$gene_id[!has_symbol]
fallback_raw <- anno_file[anno_file$gene_id %in% missing_ids,
                           c("gene_id", "product")]
colnames(fallback_raw) <- c("gene_id", "description")

# Step 3: merge and attach to sig_df
desc_df <- rbind(orgdb_anno, fallback_raw)
deg_anno <- sig_df |>
  left_join(desc_df, by = "gene_id") |>
  dplyr::select(gene_id, description, baseMean, log2FoldChange, lfcSE, pvalue, padj)

# ---- 5. Save TSV ----
dir.create("data/output/enrichment_deseq2", showWarnings = FALSE, recursive = TRUE)
dir.create("graphs/enrichment_deseq2",      showWarnings = FALSE, recursive = TRUE)

save_tsv <- function(res, name) {
  df <- as.data.frame(res)
  if (nrow(df) > 0) {
    write.table(df, paste0("data/output/enrichment_deseq2/", name, ".tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  cat(name, ":", nrow(df), "terms\n")
  invisible(df)
}

bp_df   <- save_tsv(go_bp,    "GO_BP")
mf_df   <- save_tsv(go_mf,    "GO_MF")
cc_df   <- save_tsv(go_cc,    "GO_CC")
kegg_df <- save_tsv(kegg_res, "KEGG")

write.table(deg_anno, "data/output/enrichment_deseq2/DEG_annotated.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("DEG_annotated:", nrow(deg_anno), "genes\n")

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
  ggsave(paste0("graphs/enrichment_deseq2/", name, ".pdf"), p,
         width = w, height = h, useDingbats = FALSE)
  ggsave(paste0("graphs/enrichment_deseq2/", name, ".png"), p,
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

plot_go_combined <- function(bp, mf, cc, top_n = 8) {
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
         subtitle = "DEGs: pvalue<0.05 & |LFC|>1, FDR<0.2",
         x = expression(-log[10](p.adjust)), y = NULL) +
    theme_pub() +
    theme(
      strip.text       = element_text(face = "bold", size = 10, color = "white"),
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

# ---- Generate all plots ----
sub_label <- "DEGs: pvalue<0.05 & |LFC|>1, FDR<0.2"

save_plot(plot_lollipop(bp_df,   "GO Biological Process",    sub_label),        "GO_BP_lollipop",  9, 5)
save_plot(plot_lollipop(mf_df,   "GO Molecular Function",    sub_label),        "GO_MF_lollipop",  9, 5)
save_plot(plot_lollipop(cc_df,   "GO Cellular Component",    sub_label),        "GO_CC_lollipop",  9, 5)
save_plot(plot_lollipop(kegg_df, "KEGG Pathway Enrichment",  sub_label, 20),    "KEGG_lollipop",  10, 8)

save_plot(plot_dot(bp_df,   "GO Biological Process",   sub_label),        "GO_BP_dotplot",  9, 5)
save_plot(plot_dot(mf_df,   "GO Molecular Function",   sub_label),        "GO_MF_dotplot",  9, 5)
save_plot(plot_dot(cc_df,   "GO Cellular Component",   sub_label),        "GO_CC_dotplot",  9, 5)
save_plot(plot_dot(kegg_df, "KEGG Pathway Enrichment", sub_label, 20),    "KEGG_dotplot",  10, 8)

save_plot(plot_bubble(mf_df,   "GO Molecular Function",   sub_label),     "GO_MF_bubble",  9, 5)
save_plot(plot_bubble(kegg_df, "KEGG Pathway Enrichment", sub_label, 20), "KEGG_bubble",  10, 8)

save_plot(plot_go_combined(bp_df, mf_df, cc_df, top_n = 8), "GO_combined_barplot", 9, 7)

save_plot(plot_cnet(go_bp,    "GO BP — Gene-Concept Network"),  "GO_BP_cnetplot",  10, 9)
save_plot(plot_cnet(go_mf,    "GO MF — Gene-Concept Network"),  "GO_MF_cnetplot",  10, 9)
save_plot(plot_cnet(kegg_res, "KEGG — Gene-Concept Network"),   "KEGG_cnetplot",   10, 9)

save_plot(plot_upset(go_bp,    "GO BP — Term Overlap"), "GO_BP_upset", 10, 6)
save_plot(plot_upset(go_mf,    "GO MF — Term Overlap"), "GO_MF_upset", 10, 6)
save_plot(plot_upset(kegg_res, "KEGG — Term Overlap"),  "KEGG_upset",  10, 6)

cat("\nDone. Results in data/output/enrichment_deseq2/ and graphs/enrichment_deseq2/\n")
