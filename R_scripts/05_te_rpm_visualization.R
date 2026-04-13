# TE RPM Matrix Visualization
# Figure 1: Grouped bar chart - top-level TE class RPM (SUGP1-KD vs Wild Type)
# Figure 2: Heatmap - TE sub-class RPM across all four samples
# Style: Nature journal inspired, matching project conventions

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ---- Configuration ----
TOP_CLASS_FILE <- here("data/transposon_rpm/te_top_class_rpm_matrix.tsv")
SUB_CLASS_FILE <- here("data/transposon_rpm/te_sub_class_rpm_matrix.tsv")
OUTPUT_DIR     <- here("graphs/transposon_rpm")

# Keep rows where max RPM across all samples is at least this value
RPM_THRESHOLD <- 1000

# Sample column names (fixed order for all plots)
SAMPLE_COLS <- c("SUGP1_rep1", "SUGP1_rep2", "WT_rep1", "WT_rep2")

# Condition color palette: high-contrast, modern
COND_COLORS <- c("SUGP1-KD" = "#E63946", "Wild Type" = "#264653")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: discard rows where all samples fall below RPM_THRESHOLD ----
filter_by_max_rpm <- function(df, threshold = RPM_THRESHOLD) {
  df |> filter(pmax(SUGP1_rep1, SUGP1_rep2, WT_rep1, WT_rep2) >= threshold)
}


# ==============================================================================
# Figure 1: Top-class grouped bar chart
# ==============================================================================

top_raw <- read.delim(TOP_CLASS_FILE, sep = "\t", header = TRUE)
top_df  <- filter_by_max_rpm(top_raw) |> filter(te_class != "Unknown")

cat("te_top_class: retained", nrow(top_df), "/", nrow(top_raw),
    "rows (max RPM >=", RPM_THRESHOLD, ")\n")

# Build long format with condition label
top_long <- top_df |>
  pivot_longer(all_of(SAMPLE_COLS), names_to = "sample", values_to = "rpm") |>
  mutate(
    condition = factor(
      if_else(startsWith(sample, "SUGP1"), "SUGP1-KD", "Wild Type"),
      levels = c("Wild Type", "SUGP1-KD")
    )
  )

# Per-class per-condition mean (bar heights)
top_mean <- top_long |>
  group_by(te_class, condition) |>
  summarise(mean_rpm = mean(rpm), .groups = "drop")

# Order x-axis by Wild Type mean descending
class_order <- top_mean |>
  filter(condition == "Wild Type") |>
  arrange(desc(mean_rpm)) |>
  pull(te_class)

top_long <- mutate(top_long, te_class = factor(te_class, levels = class_order))
top_mean <- mutate(top_mean, te_class = factor(te_class, levels = class_order))

p1 <- ggplot() +
  # Bars: group means
  geom_col(
    data     = top_mean,
    aes(x    = te_class, y = mean_rpm, fill = condition),
    position = position_dodge(width = 0.75),
    width    = 0.65,
    alpha    = 0.92
  ) +
  # Dots: individual replicates overlaid on bars
  geom_point(
    data     = top_long,
    aes(x    = te_class, y = rpm, fill = condition),
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.07),
    size     = 3,
    shape    = 21,
    color    = "white",
    stroke   = 0.6,
    alpha    = 0.95
  ) +
  scale_fill_manual(values = COND_COLORS) +
  scale_y_continuous(
    labels = label_number(scale = 1e-3, suffix = "K"),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(x = NULL, y = "Mean RPM", fill = NULL) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 16),
    axis.text.x        = element_text(size = 15, face = "bold"),
    axis.text.y        = element_text(size = 14),
    axis.title.y       = element_text(size = 16, face = "bold"),
    plot.margin        = margin(12, 10, 8, 10)
  )


# ==============================================================================
# Figure 3: Normalized bar chart (Wild Type = 1 per class)
# ==============================================================================

# Normalize each class by its Wild Type mean so all WT bars reach the same height
wt_ref <- top_mean |>
  filter(condition == "Wild Type") |>
  select(te_class, wt_mean = mean_rpm)

top_mean_norm <- top_mean |>
  left_join(wt_ref, by = "te_class") |>
  mutate(rel_rpm = mean_rpm / wt_mean)

top_long_norm <- top_long |>
  left_join(wt_ref, by = "te_class") |>
  mutate(rel_rpm = rpm / wt_mean)

p3 <- ggplot() +
  # Dashed reference line at Wild Type level
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  # Bars: normalized group means
  geom_col(
    data     = top_mean_norm,
    aes(x    = te_class, y = rel_rpm, fill = condition),
    position = position_dodge(width = 0.75),
    width    = 0.65,
    alpha    = 0.92
  ) +
  # Dots: individual replicates
  geom_point(
    data     = top_long_norm,
    aes(x    = te_class, y = rel_rpm, fill = condition),
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.07),
    size     = 3,
    shape    = 21,
    color    = "white",
    stroke   = 0.6,
    alpha    = 0.95
  ) +
  scale_fill_manual(values = COND_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "Relative RPM (Wild Type = 1)", fill = NULL) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 16),
    axis.text.x        = element_text(size = 15, face = "bold"),
    axis.text.y        = element_text(size = 14),
    axis.title.y       = element_text(size = 16, face = "bold"),
    plot.margin        = margin(12, 10, 8, 10)
  )


# ==============================================================================
# Figure 2: Sub-class heatmap
# ==============================================================================

sub_raw <- read.delim(SUB_CLASS_FILE, sep = "\t", header = TRUE)
sub_df  <- filter_by_max_rpm(sub_raw) |> filter(te_class != "Unknown")

cat("te_sub_class: retained", nrow(sub_df), "/", nrow(sub_raw),
    "rows (max RPM >=", RPM_THRESHOLD, ")\n")

# Geometric midpoint of all RPM values (for adaptive text color on log10 scale)
rpm_vals <- unlist(sub_df[SAMPLE_COLS])
log_mid  <- 10^(mean(log10(rpm_vals[rpm_vals > 0])))

# Row order: ascending WT mean so highest WT class appears at top of heatmap
sub_order <- sub_df |>
  mutate(wt_mean = (WT_rep1 + WT_rep2) / 2) |>
  arrange(wt_mean) |>
  pull(te_class)

sub_long <- sub_df |>
  pivot_longer(all_of(SAMPLE_COLS), names_to = "sample", values_to = "rpm") |>
  mutate(
    te_class  = factor(te_class, levels = sub_order),
    condition = factor(
      if_else(startsWith(sample, "SUGP1"), "SUGP1-KD", "Wild Type"),
      levels = c("Wild Type", "SUGP1-KD")
    ),
    # Label within each condition facet: "rep1" / "rep2"
    rep_label  = factor(gsub(".*_", "", sample), levels = c("rep1", "rep2")),
    # Text color: white on dark tiles, dark on light tiles
    text_color = if_else(rpm > log_mid, "white", "grey15"),
    tile_label = sprintf("%.1fK", rpm / 1000)
  )

p2 <- ggplot(sub_long, aes(x = rep_label, y = te_class, fill = rpm)) +
  geom_tile(color = "white", linewidth = 0.9) +
  # RPM values in K units; color adapts to tile brightness
  geom_text(aes(label = tile_label, color = text_color), size = 3.5) +
  scale_color_identity(guide = "none") +
  # Condition facets serve as the top annotation bar
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  scale_fill_gradient(
    low       = "#FFF5F0",
    high      = "#B91C1C",
    name      = "RPM",
    transform = "log10",
    labels    = label_number(scale = 1e-3, suffix = "K")
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid        = element_blank(),
    strip.background  = element_rect(fill = "#2D2D2D", color = NA),
    strip.text        = element_text(color = "white", face = "bold", size = 13),
    axis.text.x       = element_text(size = 12),
    axis.text.y       = element_text(size = 12),
    legend.key.height = unit(2.5, "lines"),
    legend.title      = element_text(size = 13),
    legend.text       = element_text(size = 12),
    panel.spacing     = unit(0.8, "lines"),
    plot.margin       = margin(10, 15, 10, 10)
  )


# ==============================================================================
# Figure 4: Sub-class normalized bar chart (Wild Type = 1 per class)
# ==============================================================================

# Build per-class per-condition mean from sub_long
sub_mean <- sub_long |>
  group_by(te_class, condition) |>
  summarise(mean_rpm = mean(rpm), .groups = "drop")

wt_ref_sub <- sub_mean |>
  filter(condition == "Wild Type") |>
  select(te_class, wt_mean = mean_rpm)

# Order x-axis by Wild Type mean descending (highest on left)
sub_bar_order <- wt_ref_sub |>
  arrange(desc(wt_mean)) |>
  pull(te_class)

sub_mean_norm <- sub_mean |>
  left_join(wt_ref_sub, by = "te_class") |>
  mutate(
    rel_rpm  = mean_rpm / wt_mean,
    te_class = factor(te_class, levels = sub_bar_order)
  )

sub_long_norm <- sub_long |>
  left_join(wt_ref_sub, by = "te_class") |>
  mutate(
    rel_rpm  = rpm / wt_mean,
    te_class = factor(te_class, levels = sub_bar_order)
  )

p4 <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_col(
    data     = sub_mean_norm,
    aes(x    = te_class, y = rel_rpm, fill = condition),
    position = position_dodge(width = 0.75),
    width    = 0.65,
    alpha    = 0.92
  ) +
  geom_point(
    data     = sub_long_norm,
    aes(x    = te_class, y = rel_rpm, fill = condition),
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.07),
    size     = 2.5,
    shape    = 21,
    color    = "white",
    stroke   = 0.6,
    alpha    = 0.95
  ) +
  scale_fill_manual(values = COND_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "Relative RPM (Wild Type = 1)", fill = NULL) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 14),
    axis.text.x        = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 13),
    axis.title.y       = element_text(size = 14, face = "bold"),
    plot.margin        = margin(12, 10, 8, 10)
  )


# ==============================================================================
# Save figures (PNG 300 dpi + PDF)
# ==============================================================================

ggsave(
  file.path(OUTPUT_DIR, "te_top_class_grouped_bar.png"),
  plot = p1, width = 8, height = 6, units = "in", dpi = 300
)
ggsave(
  file.path(OUTPUT_DIR, "te_top_class_grouped_bar.pdf"),
  plot = p1, width = 8, height = 6, units = "in"
)

ggsave(
  file.path(OUTPUT_DIR, "te_top_class_normalized_bar.png"),
  plot = p3, width = 8, height = 6, units = "in", dpi = 300
)
ggsave(
  file.path(OUTPUT_DIR, "te_top_class_normalized_bar.pdf"),
  plot = p3, width = 8, height = 6, units = "in"
)

ggsave(
  file.path(OUTPUT_DIR, "te_sub_class_normalized_bar.png"),
  plot = p4, width = 12, height = 6, units = "in", dpi = 300
)
ggsave(
  file.path(OUTPUT_DIR, "te_sub_class_normalized_bar.pdf"),
  plot = p4, width = 12, height = 6, units = "in"
)

ggsave(
  file.path(OUTPUT_DIR, "te_sub_class_heatmap.png"),
  plot = p2, width = 7, height = 5, units = "in", dpi = 300
)
ggsave(
  file.path(OUTPUT_DIR, "te_sub_class_heatmap.pdf"),
  plot = p2, width = 7, height = 5, units = "in"
)

cat("\nFigures saved to:", OUTPUT_DIR, "\n")
cat("  te_top_class_grouped_bar.png / .pdf\n")
cat("  te_top_class_normalized_bar.png / .pdf\n")
cat("  te_sub_class_normalized_bar.png / .pdf\n")
cat("  te_sub_class_heatmap.png / .pdf\n")
