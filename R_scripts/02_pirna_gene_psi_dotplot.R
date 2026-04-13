# PSI dot plot for known piRNA genes — grouped by gene, colored by AS type
# Input:  data/rmats/comparison_1_control_vs_sugp1/*.MATS.JC.txt
#         data/known_pirna_gene_id.csv  (columns: GeneID, gene_name)
# Output: graphs/pirna_gene_psi_dotplot.{pdf,png}

library(dplyr)
library(ggplot2)

# ── 1. Load gene mapping ───────────────────────────────────────────────────────
gene_map <- read.csv("data/known_pirna_gene_id.csv",
                     stringsAsFactors = FALSE) |>
  distinct(GeneID, .keep_all = TRUE)          # drop duplicate GeneID rows

# ── 2. Read & filter rMATS files ──────────────────────────────────────────────
AS_TYPES  <- c("SE", "RI", "A3SS", "A5SS", "MXE")
RMATS_DIR <- "data/rmats/comparison_1_control_vs_sugp1"

read_as_file <- function(as_type) {
  path <- file.path(RMATS_DIR, paste0(as_type, ".MATS.JC.txt"))
  df   <- read.table(path, header = TRUE, sep = "\t",
                     quote = "", comment.char = "", check.names = TRUE)
  # GeneID column has surrounding quotes — strip them
  df$GeneID_clean <- gsub('"', '', df$GeneID)
  df |>
    filter(GeneID_clean %in% gene_map$GeneID) |>
    select(ID, GeneID = GeneID_clean, IncLevel1, IncLevel2,
           FDR, IncLevelDifference) |>
    mutate(as_type  = as_type,
           event_id = paste0(as_type, "_", ID))
}

events <- lapply(AS_TYPES, read_as_file) |>
  do.call(rbind, args = _) |>
  left_join(gene_map, by = "GeneID")

if (nrow(events) == 0) stop("No events found for the specified piRNA genes.")

# ── 3. Expand comma-separated PSI → long format ───────────────────────────────
expand_psi <- function(df, psi_col, group_label) {
  rows <- lapply(seq_len(nrow(df)), function(i) {
    vals <- suppressWarnings(
      as.numeric(strsplit(as.character(df[[psi_col]][i]), ",")[[1]])
    )
    data.frame(
      gene_name = df$gene_name[i],
      GeneID    = df$GeneID[i],
      event_id  = df$event_id[i],
      as_type   = df$as_type[i],
      PSI       = vals,
      group     = group_label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

long_df <- rbind(
  expand_psi(events, "IncLevel1", "Control"),
  expand_psi(events, "IncLevel2", "SUGP1_KD")
) |>
  filter(!is.na(PSI)) |>
  # Average replicates: one mean PSI per (event, group)
  group_by(gene_name, GeneID, event_id, as_type, group) |>
  summarise(PSI = mean(PSI), .groups = "drop")

# ── 4. Compute numeric x positions ────────────────────────────────────────────
# Layout: events within a gene are spread at 0.4 spacing;
#         genes are separated by a gap of 1.2.
EVENT_SPACING <- 1.0
GENE_GAP      <- 2.2
GROUP_OFFSET  <- 0.08   # Control left / KD right within each event

event_layout <- long_df |>
  distinct(gene_name, event_id, as_type) |>
  arrange(gene_name, as_type, event_id) |>
  group_by(gene_name) |>
  mutate(event_rank = row_number()) |>
  ungroup()

gene_widths <- event_layout |>
  group_by(gene_name) |>
  summarise(n_ev = n(), .groups = "drop") |>
  mutate(block_width = n_ev * EVENT_SPACING,
         block_start = cumsum(lag(block_width + GENE_GAP, default = 0)))

event_layout <- event_layout |>
  left_join(gene_widths, by = "gene_name") |>
  mutate(x_event = block_start + (event_rank - 1) * EVENT_SPACING)

gene_axis <- gene_widths |>
  left_join(
    event_layout |> group_by(gene_name) |> summarise(x_center = mean(x_event)),
    by = "gene_name"
  )

# Separator x positions (midpoints between last event of one gene and first of next)
separator_x <- gene_widths |>
  mutate(x_end   = block_start + (n_ev - 1) * EVENT_SPACING,
         x_sep   = x_end + GENE_GAP / 2) |>
  filter(row_number() < n()) |>
  pull(x_sep)

# Annotation: ΔΨ and FDR per event (one row per event, pre-flip coordinates)
event_annot <- events |>
  select(event_id, FDR, IncLevelDifference) |>
  left_join(event_layout |> select(event_id, x_event), by = "event_id") |>
  mutate(label = sprintf("dPSI=%.2f  FDR=%.2f", IncLevelDifference, FDR))

long_df <- long_df |>
  left_join(event_layout |> select(gene_name, event_id, x_event),
            by = c("gene_name", "event_id")) |>
  mutate(x_final = x_event + ifelse(group == "Control", -GROUP_OFFSET, GROUP_OFFSET))

# ── 5. Plot ────────────────────────────────────────────────────────────────────
AS_COLORS <- c(
  "SE"   = "#E41A1C",
  "RI"   = "#377EB8",
  "A3SS" = "#4DAF4A",
  "A5SS" = "#FF7F00",
  "MXE"  = "#984EA3"
)

p <- ggplot(long_df, aes(x = x_final, y = PSI,
                          color = as_type, shape = group)) +
  # Gene separator lines
  geom_vline(xintercept = separator_x,
             linetype = "dashed", color = "grey70", linewidth = 0.4) +
  # Connecting line between Control and SUGP1_KD: flat line = no PSI change
  geom_line(aes(group = event_id), color = "grey55", linewidth = 0.8, alpha = 0.7) +
  # Points (mean PSI per event × group, no jitter needed)
  geom_point(size = 5, alpha = 0.9) +
  # ΔΨ and FDR annotation to the right of each event
  geom_text(data = event_annot,
            aes(x = x_event, y = 1.06, label = label),
            hjust = 0, size = 4, color = "grey35", inherit.aes = FALSE) +
  scale_color_manual(values = AS_COLORS, name = "AS Type") +
  scale_shape_manual(values = c("Control" = 16, "SUGP1_KD" = 17),
                     name = "Group") +
  scale_x_continuous(
    breaks = gene_axis$x_center,
    labels = gene_axis$gene_name,
    expand = expansion(add = 0.6)
  ) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
  coord_flip(clip = "off") +
  labs(
    title    = "PSI of AS Events in Known piRNA Pathway Genes",
    subtitle = "Control (circle) vs SUGP1 KD (triangle); connecting line shows PSI change per event; all FDR > 0.05",
    x        = NULL,
    y        = "PSI (IncLevel)"
  ) +
  theme_bw(base_size = 26) +
  theme(
    axis.text.y        = element_text(size = 24),
    axis.text.x        = element_text(size = 24),
    panel.grid.major.y = element_blank(),
    legend.position    = "top",
    legend.box         = "horizontal"
  )

# ── 6. Save PDF and PNG ───────────────────────────────────────────────────────
dir.create("graphs", showWarnings = FALSE)
n_events <- n_distinct(long_df$event_id)
plot_w   <- 18
plot_h   <- max(8, n_events * 0.75)

out_base <- "graphs/pirna_gene_psi_dotplot"
ggsave(paste0(out_base, ".pdf"), p,
       width = plot_w, height = plot_h, limitsize = FALSE)
ggsave(paste0(out_base, ".png"), p,
       width = plot_w, height = plot_h, dpi = 150, limitsize = FALSE)

cat("Saved:", out_base, ".pdf / .png\n")
cat("Events plotted:", n_events, "| AS types:", n_distinct(long_df$as_type), "\n")
cat("Genes:", paste(unique(long_df$gene_name), collapse = ", "), "\n")
