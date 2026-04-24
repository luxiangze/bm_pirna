#' Plot piRNA Length Distribution (RPM-normalized) - Paired Comparison
#'
#' This script generates paired comparison plots for piRNA length distribution.
#' Each pair from sample_map.csv generates a separate comparison plot.
#' Reads are normalized using RPM (Reads Per Million) based on filtered_reads from filtered_read_counts.tsv.
#'
#' @author Yongkang Guo
#' @date 2025-01-23

# Load required packages
library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Define file paths
DATA_DIR <- here("data")
PIRNA_STATS_DIR <- file.path(DATA_DIR, "piRNA_mirna_merged_stats")
SAMPLE_MAP_FILE <- here("data", "sample_map.csv")
FILTER_SUMMARY_FILE <- file.path(DATA_DIR, "structure_np_miRNA_filtered_read_counts.tsv")
OUTPUT_DIR <- here("graphs", "piRNA_mirna_merged_stats")

# Create output directory if not exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

#' Load data files
#'
#' @return List containing loaded data frames
load_data <- function() {
  # Load length distribution
  length_dist <- read_tsv(
    file.path(PIRNA_STATS_DIR, "length_distribution.tsv"),
    col_types = cols(
      sample = col_character(),
      length = col_integer(),
      count = col_integer()
    )
  )

  # Load filter summary for effective reads (filtered_reads)
  filter_summary <- read_tsv(
    FILTER_SUMMARY_FILE,
    col_types = cols(
      sample = col_character(),
      input_reads = col_integer(),
      filtered_reads = col_integer(),
      retained_pct = col_double()
    )
  )

  # Load sample map
  sample_map <- read_csv(SAMPLE_MAP_FILE, col_types = cols())

  list(
    length_dist = length_dist,
    filter_summary = filter_summary,
    sample_map = sample_map
  )
}

#' Calculate RPM (Reads Per Million) normalization
#'
#' @param count Raw count
#' @param total_reads Total effective reads for normalization
#' @return RPM value
calculate_rpm <- function(count, total_reads) {
  (count / total_reads) * 1e6
}

#' Process length distribution data with RPM normalization
#'
#' @param data Data list from load_data()
#' @return Processed data frame with mean ± SD RPM per group
process_length_dist <- function(data) {
  data$length_dist %>%
    left_join(
      data$filter_summary %>% select(sample, filtered_reads),
      by = "sample"
    ) %>%
    mutate(
      rpm = calculate_rpm(count, filtered_reads),
      # Extract group name by removing _repX suffix
      group_name = str_replace(sample, "_rep\\d+$", "")
    ) %>%
    group_by(group_name, length) %>%
    summarise(
      rpm_mean = mean(rpm, na.rm = TRUE),
      rpm_sd   = sd(rpm, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(rpm_sd = replace_na(rpm_sd, 0))
}

#' Plot piRNA length distribution for a group comparison with modern styling
#'
#' @param length_dist_df Full length distribution data (mean ± SD per group)
#' @param control_group Control group name (e.g., "WT")
#' @param treated_group Treated group name (e.g., "SUGP1")
#' @return ggplot object
plot_length_group <- function(length_dist_df, control_group, treated_group) {
  plot_data <- length_dist_df %>%
    filter(group_name %in% c(control_group, treated_group)) %>%
    mutate(
      group = case_when(
        group_name == control_group ~ "Control",
        group_name == treated_group ~ "Treated"
      ),
      group = factor(group, levels = c("Control", "Treated"))
    )

  color_palette <- c(
    "Control" = "#457B9D",  # Steel Blue
    "Treated" = "#E63946"   # Vivid Red
  )

  ggplot(plot_data, aes(x = length, y = rpm_mean, color = group, fill = group)) +
    geom_ribbon(
      aes(ymin = rpm_mean - rpm_sd, ymax = rpm_mean + rpm_sd),
      alpha = 0.2, color = NA
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5, shape = 16) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    labs(
      x = "piRNA Length (nt)",
      y = "Normalized counts (RPM)"
    ) +
    theme_minimal() +
    theme(
      axis.text        = element_text(size = 20, color = "#1D3557"),
      axis.title       = element_text(size = 22, color = "#1D3557"),
      legend.position  = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.text      = element_text(size = 20, color = "#1D3557"),
      legend.title     = element_text(size = 22, face = "bold", color = "#1D3557"),
      panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.5),
      panel.grid.major.x = element_line(color = "#E8E8E8", linewidth = 0.5),
      panel.background = element_rect(fill = "#FFFFFF", color = "#1D3557", linewidth = 0.3),
      plot.margin      = margin(10, 10, 10, 10)
    ) +
    scale_x_continuous(breaks = 24:35, limits = c(24, 35)) +
    guides(
      color = guide_legend(override.aes = list(linewidth = 1.2, size = 4)),
      fill  = "none"
    )
}

#' Generate all group comparison plots
#'
#' @param data Data list from load_data()
generate_all_plots <- function(data) {
  # Process data
  length_dist_df <- process_length_dist(data)

  # Generate plot for each group pair in sample_map
  for (i in seq_len(nrow(data$sample_map))) {
    control_group <- data$sample_map$Control_sample[i]
    treated_group <- data$sample_map$Treated_sample[i]

    # Base filename (lowercase)
    base_name <- str_to_lower(treated_group)

    message(sprintf("Generating plot for group: %s vs %s", control_group, treated_group))

    # Create plot (all replicates in each group)
    p <- plot_length_group(length_dist_df, control_group, treated_group)

    # Save as PNG (300ppi)
    png_file <- file.path(OUTPUT_DIR, sprintf("pirna_length_distribution_%s_vs_%s.png", base_name, str_to_lower(control_group)))
    ggsave(
      filename = png_file,
      plot = p,
      width = 8,
      height = 5,
      units = "in",
      dpi = 300
    )
    message(sprintf("  Saved PNG: %s", png_file))

    # Save as PDF (vector format)
    pdf_file <- file.path(OUTPUT_DIR, sprintf("pirna_length_distribution_%s_vs_%s.pdf", base_name, str_to_lower(control_group)))
    ggsave(
      filename = pdf_file,
      plot = p,
      width = 8,
      height = 5,
      units = "in",
      device = pdf
    )
    message(sprintf("  Saved PDF: %s", pdf_file))
  }

  message("All plots saved to ", OUTPUT_DIR)
}

# Main execution
main <- function() {
  message("Loading data...")
  data <- load_data()

  message("Processing length distribution...")
  generate_all_plots(data)

  message("Done!")
}

# Run main function
main()
