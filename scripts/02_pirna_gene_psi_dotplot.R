#' PSI dot plot for known piRNA pathway genes
#'
#' Read rMATS JC tables from the current Bombyx mori splicing workflow,
#' keep events whose GeneID matches data/external/known_pirna_gene_id.csv,
#' and plot mean PSI per event for each comparison group.
#'
#' Run from the repository root:
#'   Rscript --vanilla R_scripts/02_pirna_gene_psi_dotplot.R
#'
#' @author Yongkang Guo

library(dplyr)
library(readr)
library(ggplot2)

SPLICING_DIR <- "workflow/rnasplicing_bm_ncbi_20260824"
GENE_MAP_FILE <- "data/external/known_pirna_gene_id.csv"
OUTPUT_DIR <- file.path("data/processed", basename(SPLICING_DIR), "plots")
AS_TYPES <- c("SE", "RI", "A3SS", "A5SS", "MXE")

COMPARISONS <- list(
  control_vs_sugp1 = list(
    rmats_dir = file.path(SPLICING_DIR, "rmats/comparison_1_control_vs_sugp1"),
    group1 = "Control",
    group2 = "SUGP1"
  ),
  p200_vs_sugp1 = list(
    rmats_dir = file.path(SPLICING_DIR, "rmats/comparison_2_p200_vs_sugp1"),
    group1 = "P200",
    group2 = "SUGP1"
  )
)

AS_COLORS <- c(
  "SE"   = "#E41A1C",
  "RI"   = "#377EB8",
  "A3SS" = "#4DAF4A",
  "A5SS" = "#FF7F00",
  "MXE"  = "#984EA3"
)

EVENT_SPACING <- 1.0
GENE_GAP <- 2.2
GROUP_OFFSET <- 0.08

#' Load known piRNA GeneID -> display name
#'
#' Duplicate GeneIDs (e.g. Qin/DDx43) are collapsed with "/".
load_gene_map <- function(path) {
  if (!file.exists(path)) {
    stop("Known piRNA table not found: ", path)
  }
  read_csv(path, show_col_types = FALSE) |>
    filter(!is.na(GeneID), GeneID != "") |>
    group_by(GeneID) |>
    summarise(gene_name = paste(unique(gene_name), collapse = "/"), .groups = "drop")
}

#' Collect rMATS JC events for known piRNA genes
read_known_events <- function(rmats_dir, gene_map) {
  known_ids <- unique(c(gene_map$GeneID, gene_map$gene_name))

  read_as_file <- function(as_type) {
    path <- file.path(rmats_dir, paste0(as_type, ".MATS.JC.txt"))
    if (!file.exists(path)) {
      return(NULL)
    }
    # rMATS JC tables have two ID columns; check.names makes them ID and ID.1
    df <- read.table(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = TRUE
    )
    df$GeneID <- gsub('"', "", as.character(df$GeneID))
    if ("geneSymbol" %in% names(df)) {
      df$geneSymbol <- gsub('"', "", as.character(df$geneSymbol))
    } else {
      df$geneSymbol <- NA_character_
    }
    df |>
      filter(GeneID %in% known_ids | geneSymbol %in% known_ids) |>
      transmute(
        ID,
        GeneID,
        IncLevel1 = as.character(IncLevel1),
        IncLevel2 = as.character(IncLevel2),
        FDR,
        IncLevelDifference,
        as_type = as_type,
        event_id = paste0(as_type, "_", ID)
      )
  }

  events <- bind_rows(lapply(AS_TYPES, read_as_file))
  if (nrow(events) == 0) {
    return(events)
  }

  events |>
    left_join(gene_map, by = "GeneID") |>
    mutate(
      gene_name = if_else(
        is.na(gene_name) | gene_name == "",
        GeneID,
        gene_name
      )
    )
}

#' Expand comma-separated PSI values then average replicates
expand_mean_psi <- function(events, group1, group2) {
  expand_psi <- function(df, psi_col, group_label) {
    bind_rows(lapply(seq_len(nrow(df)), function(i) {
      vals <- suppressWarnings(
        as.numeric(strsplit(as.character(df[[psi_col]][i]), ",", fixed = TRUE)[[1]])
      )
      tibble(
        gene_name = df$gene_name[i],
        GeneID = df$GeneID[i],
        event_id = df$event_id[i],
        as_type = df$as_type[i],
        PSI = vals,
        group = group_label
      )
    }))
  }

  bind_rows(
    expand_psi(events, "IncLevel1", group1),
    expand_psi(events, "IncLevel2", group2)
  ) |>
    filter(!is.na(PSI)) |>
    group_by(gene_name, GeneID, event_id, as_type, group) |>
    summarise(PSI = mean(PSI), .groups = "drop")
}

build_layout <- function(long_df) {
  event_layout <- long_df |>
    distinct(gene_name, event_id, as_type) |>
    arrange(gene_name, as_type, event_id) |>
    group_by(gene_name) |>
    mutate(event_rank = row_number()) |>
    ungroup()

  gene_widths <- event_layout |>
    group_by(gene_name) |>
    summarise(n_ev = n(), .groups = "drop") |>
    mutate(
      block_width = n_ev * EVENT_SPACING,
      block_start = cumsum(lag(block_width + GENE_GAP, default = 0))
    )

  event_layout <- event_layout |>
    left_join(gene_widths, by = "gene_name") |>
    mutate(x_event = block_start + (event_rank - 1) * EVENT_SPACING)

  gene_axis <- gene_widths |>
    left_join(
      event_layout |> group_by(gene_name) |> summarise(x_center = mean(x_event), .groups = "drop"),
      by = "gene_name"
    )

  separator_x <- gene_widths |>
    mutate(
      x_end = block_start + (n_ev - 1) * EVENT_SPACING,
      x_sep = x_end + GENE_GAP / 2
    ) |>
    filter(row_number() < n()) |>
    pull(x_sep)

  list(
    event_layout = event_layout,
    gene_axis = gene_axis,
    separator_x = separator_x
  )
}

plot_psi <- function(long_df, events, group1, group2, title_suffix) {
  layout <- build_layout(long_df)
  group_levels <- c(group1, group2)

  event_annot <- events |>
    select(event_id, FDR, IncLevelDifference) |>
    left_join(layout$event_layout |> select(event_id, x_event), by = "event_id") |>
    mutate(label = sprintf("dPSI=%.2f  FDR=%.2g", IncLevelDifference, FDR))

  plot_df <- long_df |>
    left_join(
      layout$event_layout |> select(gene_name, event_id, x_event),
      by = c("gene_name", "event_id")
    ) |>
    mutate(
      group = factor(group, levels = group_levels),
      x_final = x_event + if_else(group == group1, -GROUP_OFFSET, GROUP_OFFSET)
    )

  shape_vals <- setNames(c(16, 17), group_levels)
  n_sig <- sum(events$FDR < 0.05, na.rm = TRUE)

  ggplot(plot_df, aes(x = x_final, y = PSI, color = as_type, shape = group)) +
    geom_vline(
      xintercept = layout$separator_x,
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.4
    ) +
    geom_line(aes(group = event_id), color = "grey55", linewidth = 0.8, alpha = 0.7) +
    geom_point(size = 5, alpha = 0.9) +
    geom_text(
      data = event_annot,
      aes(x = x_event, y = 1.06, label = label),
      hjust = 0,
      size = 4,
      color = "grey35",
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = AS_COLORS, name = "AS Type") +
    scale_shape_manual(values = shape_vals, name = "Group") +
    scale_x_continuous(
      breaks = layout$gene_axis$x_center,
      labels = layout$gene_axis$gene_name,
      expand = expansion(add = 0.6)
    ) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
    coord_flip(clip = "off") +
    labs(
      title = "PSI of AS events in known piRNA pathway genes",
      subtitle = sprintf(
        "%s (circle) vs %s (triangle); %d / %d events with FDR < 0.05",
        group1, group2, n_sig, nrow(events)
      ),
      x = NULL,
      y = "PSI (IncLevel)"
    ) +
    theme_bw(base_size = 26) +
    theme(
      axis.text.y = element_text(size = 24),
      axis.text.x = element_text(size = 24),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.box = "horizontal"
    )
}

save_plot <- function(p, out_base, n_events) {
  plot_w <- 18
  plot_h <- max(8, n_events * 0.75)
  ggsave(
    paste0(out_base, ".pdf"), p,
    width = plot_w, height = plot_h, limitsize = FALSE
  )
  ggsave(
    paste0(out_base, ".png"), p,
    width = plot_w, height = plot_h, dpi = 150, limitsize = FALSE
  )
  message("Saved ", out_base, ".pdf / .png")
}

main <- function() {
  gene_map <- load_gene_map(GENE_MAP_FILE)
  message("Known piRNA GeneIDs: ", nrow(gene_map))
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

  for (comp_name in names(COMPARISONS)) {
    spec <- COMPARISONS[[comp_name]]
    message("Comparison: ", comp_name)
    events <- read_known_events(spec$rmats_dir, gene_map)
    if (nrow(events) == 0) {
      warning("No events found for ", comp_name, "; skipping")
      next
    }

    long_df <- expand_mean_psi(events, spec$group1, spec$group2)
    p <- plot_psi(long_df, events, spec$group1, spec$group2, comp_name)
    out_base <- file.path(OUTPUT_DIR, paste0("pirna_gene_psi_dotplot_", comp_name))
    save_plot(p, out_base, n_distinct(long_df$event_id))
    message(
      "  events: ", n_distinct(long_df$event_id),
      " | genes: ", paste(sort(unique(long_df$gene_name)), collapse = ", ")
    )
  }
}

main()
