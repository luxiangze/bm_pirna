#' Visualize miRTrace Length Distribution (RPM-normalized)
#'
#' Reads four per-sample mirtrace-stats-length.tsv files, normalizes each
#' sample to RPM using the total read count (LENGTH 0-50) as the normalization
#' factor, then plots mean +/- SD ribbon lines for SUGP1-KD vs WT groups.
#'
#' Input assumption: each TSV has 2 columns (LENGTH 0-50, raw count); 51 data rows.
#' Output: PNG (300 dpi) and PDF saved to graphs/mirtrace/.

library(here)
library(tidyverse)
library(scales)

# ---- 1. Configuration -------------------------------------------------------

SAMPLES <- tibble::tribble(
  ~sample_id,    ~group,    ~file_path,
  "SUGP1_rep1",  "SUGP1",   "data/mirtrace/SUGP1_rep1/mirtrace-stats-length.tsv",
  "SUGP1_rep2",  "SUGP1",   "data/mirtrace/SUGP1_rep2/mirtrace-stats-length.tsv",
  "WT_rep1",     "WT",      "data/mirtrace/WT_rep1/mirtrace-stats-length.tsv",
  "WT_rep2",     "WT",      "data/mirtrace/WT_rep2/mirtrace-stats-length.tsv"
)

X_RANGE       <- c(15L, 45L)
COLOR_PALETTE <- c("WT" = "#264653", "SUGP1" = "#E63946")
OUTPUT_DIR    <- here("graphs", "mirtrace")

# ---- 2. I/O -----------------------------------------------------------------

#' Read a single mirtrace-stats-length.tsv and attach sample metadata.
#' The second column name varies per file, so we rename by position.
read_length_tsv <- function(file_path, sample_id) {
  read_tsv(here(file_path), col_types = cols(.default = col_double()),
           show_col_types = FALSE) |>
    set_names(c("length", "count")) |>
    mutate(
      length = as.integer(length),
      count  = as.integer(count),
      sample = sample_id
    )
}

load_all_samples <- function(sample_tbl) {
  pmap(
    sample_tbl |> select(sample_id, file_path),
    \(sample_id, file_path) read_length_tsv(file_path, sample_id)
  ) |>
    bind_rows() |>
    left_join(sample_tbl |> select(sample_id, group),
              by = c("sample" = "sample_id"))
}

# ---- 3. Normalization -------------------------------------------------------

#' Normalize raw counts to RPM.
#' Normalization factor = sum of all reads per sample (LENGTH 0-50).
normalize_to_rpm <- function(raw_df) {
  norm_factors <- raw_df |>
    group_by(sample) |>
    summarise(total_reads = sum(count), .groups = "drop")

  message("Normalization factors (total reads per sample):")
  print(norm_factors)

  raw_df |>
    left_join(norm_factors, by = "sample") |>
    mutate(rpm = count / total_reads * 1e6)
}

# ---- 4. Group statistics ----------------------------------------------------

#' Compute per-group mean and SD of RPM for the specified length range.
compute_group_stats <- function(rpm_df, x_range = X_RANGE) {
  rpm_df |>
    filter(length >= x_range[1], length <= x_range[2]) |>
    mutate(group = factor(group, levels = c("WT", "SUGP1"))) |>
    group_by(group, length) |>
    summarise(
      rpm_mean = mean(rpm, na.rm = TRUE),
      rpm_sd   = sd(rpm,  na.rm = TRUE),
      .groups  = "drop"
    ) |>
    mutate(rpm_sd = replace_na(rpm_sd, 0))
}

# ---- 5. Plotting ------------------------------------------------------------

plot_length_dist <- function(plot_df) {
  ggplot(plot_df, aes(x = length, y = rpm_mean,
                      color = group, fill = group)) +
    # Draw ribbon first so lines and points render on top
    geom_ribbon(
      aes(ymin = rpm_mean - rpm_sd, ymax = rpm_mean + rpm_sd),
      alpha = 0.15, color = NA
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3, shape = 16) +
    scale_color_manual(values = COLOR_PALETTE) +
    scale_fill_manual(values  = COLOR_PALETTE) +
    scale_x_continuous(
      breaks       = seq(X_RANGE[1], X_RANGE[2], by = 5),
      minor_breaks = seq(X_RANGE[1], X_RANGE[2], by = 1),
      limits       = X_RANGE
    ) +
    scale_y_continuous(labels = comma) +
    labs(
      x     = "RNA Length (nt)",
      y     = "Normalized counts (RPM)",
      color = "Condition",
      fill  = "Condition"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text          = element_text(size = 20, color = "#1D3557"),
      axis.title         = element_text(size = 22, color = "#1D3557"),
      legend.position    = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.text        = element_text(size = 20, color = "#1D3557"),
      legend.title       = element_text(size = 22, face = "bold",
                                        color = "#1D3557"),
      panel.grid.major   = element_line(color = "#E8E8E8", linewidth = 0.5),
      panel.grid.minor   = element_line(color = "#F4F4F4", linewidth = 0.3),
      panel.background   = element_rect(fill = "#FFFFFF", color = "#1D3557",
                                        linewidth = 0.3),
      plot.margin        = margin(10, 15, 10, 10)
    ) +
    guides(
      color = guide_legend(override.aes = list(linewidth = 1.2, size = 4)),
      fill  = "none"
    )
}

# ---- 6. Save ----------------------------------------------------------------

save_plot <- function(p, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  base <- file.path(output_dir, "mirtrace_length_distribution_sugp1_vs_wt")

  ggsave(paste0(base, ".png"), plot = p,
         width = 9, height = 5.5, units = "in", dpi = 300)
  message("Saved PNG: ", paste0(base, ".png"))

  ggsave(paste0(base, ".pdf"), plot = p,
         width = 9, height = 5.5, units = "in", device = pdf)
  message("Saved PDF: ", paste0(base, ".pdf"))
}

# ---- 7. Main ----------------------------------------------------------------

main <- function() {
  message("Loading mirtrace length data...")
  raw <- load_all_samples(SAMPLES)

  # Sanity check: 51 rows per sample, LENGTH 0-50
  stopifnot(all(table(raw$sample) == 51L))
  stopifnot(all(range(raw$length) == c(0L, 50L)))

  message("Normalizing to RPM...")
  rpm <- normalize_to_rpm(raw)

  # Verify: total RPM per sample should equal 1e6
  rpm_check <- rpm |> group_by(sample) |> summarise(rpm_total = sum(rpm))
  message("RPM totals per sample (should be ~1,000,000):")
  print(rpm_check)

  message("Computing group mean +/- SD...")
  stats <- compute_group_stats(rpm)

  message("Plotting...")
  p <- plot_length_dist(stats)

  message("Saving outputs to: ", OUTPUT_DIR)
  save_plot(p, OUTPUT_DIR)

  message("Done.")
}

main()
