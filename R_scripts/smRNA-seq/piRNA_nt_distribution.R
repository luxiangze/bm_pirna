#' Plot Combined First and Tenth Nucleotide Distribution of piRNA - Paired Comparison
#'
#' This script generates elegant combined comparison plots for nucleotide distribution at 1st and 10th positions.
#' Each pair from sample_map.csv generates separate modern and aesthetically pleasing comparison plots.
#'
#' @author Claude
#' @date 2025-01-23

# Load required packages
library(here)
library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)

# Define file paths
DATA_DIR <- here("data")
PIRNA_STATS_DIR <- file.path(DATA_DIR, "transposon_dedrived_pirna_stats")
SAMPLE_MAP_FILE <- here("data", "sample_map.csv")
OUTPUT_DIR <- here("graphs", "transposon_dedrived_pirna_stats")

# Create output directory if not exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

#' Load data files
#'
#' @return List containing loaded data frames
load_data <- function() {
  # Load first nucleotide distribution
  first_nt <- read_tsv(
    file.path(PIRNA_STATS_DIR, "first_nt_distribution.tsv"),
    col_types = cols(
      sample = col_character(),
      nucleotide = col_character(),
      count = col_integer(),
      percentage = col_double()
    )
  )

  # Load tenth nucleotide distribution
  tenth_nt <- read_tsv(
    file.path(PIRNA_STATS_DIR, "tenth_nt_distribution.tsv"),
    col_types = cols(
      sample = col_character(),
      nucleotide = col_character(),
      count = col_integer(),
      percentage = col_double()
    )
  )

  # Load sample map
  sample_map <- read_csv(SAMPLE_MAP_FILE, col_types = cols())

  list(
    first_nt = first_nt,
    tenth_nt = tenth_nt,
    sample_map = sample_map
  )
}

#' Process nucleotide distribution data for group comparison (calculate mean across replicates)
#'
#' @param nt_data Nucleotide distribution data frame
#' @param control_group Control group name (e.g., "WT")
#' @param treated_group Treated group name (e.g., "SUGP1")
#' @return Processed data frame with mean values per group
process_nt_group_data <- function(nt_data, control_group, treated_group) {
  nt_data |>
    mutate(group_name = str_replace(sample, "_rep\\d+$", "")) |>
    filter(group_name %in% c(control_group, treated_group)) |>
    # Recode T -> U for display (piRNA are RNA molecules)
    mutate(nucleotide = recode(nucleotide, "T" = "U")) |>
    # Normalize percentages to sum to 100% within each sample
    group_by(sample) |>
    mutate(percentage = percentage / sum(percentage) * 100) |>
    ungroup() |>
    # Calculate mean across replicates per group and nucleotide
    group_by(nucleotide, group_name) |>
    summarise(mean_percentage = mean(percentage), .groups = "drop") |>
    mutate(
      group = case_when(
        group_name == control_group ~ "Control",
        group_name == treated_group ~ "Treated"
      ),
      group = factor(group, levels = c("Control", "Treated"))
    )
}

#' Create modern combined nucleotide distribution plot using mean values
#'
#' @param first_data Processed first nucleotide data (mean per group)
#' @param tenth_data Processed tenth nucleotide data (mean per group)
#' @param treatment_name Treatment name
#' @return Patchwork plot object
create_modern_plot <- function(first_data, tenth_data, treatment_name) {
  color_palette <- c(
    "A" = "#D32F2F",   # Deep Red
    "U" = "#457B9D",   # Steel Blue
    "C" = "#2A9D8F",   # Persian Green
    "G" = "#F4A261"    # Sandy Orange
  )

  bar_theme <- theme_minimal() +
    theme(
      axis.text.x    = element_text(size = 22, color = "#1D3557"),
      axis.text.y    = element_text(size = 22, color = "#1D3557"),
      axis.title.y   = element_text(size = 24, color = "#1D3557", margin = margin(r = 5)),
      axis.title.x   = element_blank(),
      legend.position  = "right",
      legend.text      = element_text(size = 22, color = "#1D3557"),
      legend.title     = element_text(size = 24, face = "bold", color = "#1D3557"),
      legend.key.size  = unit(1.5, "lines"),
      panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.5),
      panel.grid.major.x = element_blank(),
      panel.background   = element_rect(fill = "#FFFFFF", color = "#1D3557", linewidth = 0.3),
      plot.margin        = margin(10, 10, 10, 10)
    )

  p_1st <- ggplot(first_data, aes(x = group, y = mean_percentage, fill = nucleotide)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = color_palette, name = "Nucleotide", drop = FALSE) +
    scale_y_continuous(labels = percent_format(scale = 1), breaks = seq(0, 100, 20)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(y = "Percentage (%)") +
    bar_theme +
    guides(fill = guide_legend(nrow = 4, byrow = FALSE))

  p_10th <- ggplot(tenth_data, aes(x = group, y = mean_percentage, fill = nucleotide)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = color_palette, name = "Nucleotide", drop = FALSE) +
    scale_y_continuous(labels = percent_format(scale = 1), breaks = seq(0, 100, 20)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(y = "Percentage (%)") +
    bar_theme +
    guides(fill = "none")

  (p_1st + labs(tag = "1st Position")) / (p_10th + labs(tag = "10th Position")) +
    plot_annotation(
      theme = theme(plot.tag = element_text(size = 28, face = "bold", color = "#1D3557"))
    )
}

#' Generate all group comparison plots for nucleotide distribution
#'
#' @param data Data list from load_data()
generate_all_plots <- function(data) {
  # Generate plot for each group pair in sample_map
  for (i in seq_len(nrow(data$sample_map))) {
    control_group <- data$sample_map$Control_sample[i]
    treated_group <- data$sample_map$Treated_sample[i]

    # Treatment name (capitalize first letter)
    treatment_name <- paste0(toupper(substr(treated_group, 1, 1)), substr(treated_group, 2, nchar(treated_group)))

    # Base filename (lowercase for consistency)
    base_name <- str_to_lower(treated_group)
    control_name <- str_to_lower(control_group)

    message(sprintf("Generating plots for group: %s vs %s", control_group, treated_group))

    # Process 1st nucleotide data (mean across replicates)
    first_nt_group <- process_nt_group_data(data$first_nt, control_group, treated_group)

    # Process 10th nucleotide data (mean across replicates)
    tenth_nt_group <- process_nt_group_data(data$tenth_nt, control_group, treated_group)

    # Create combined modern plot using patchwork
    p_combined <- create_modern_plot(first_nt_group, tenth_nt_group, treatment_name)

    # Save as PNG (300ppi)
    png_file <- file.path(OUTPUT_DIR, sprintf("pirna_nt_distribution_%s_vs_%s_combined.png", base_name, control_name))
    ggsave(
      filename = png_file,
      plot = p_combined,
      width = 10,
      height = 10,
      units = "in",
      dpi = 300
    )
    message(sprintf("  Saved PNG: %s", png_file))

    # Save as PDF (vector format)
    pdf_file <- file.path(OUTPUT_DIR, sprintf("pirna_nt_distribution_%s_vs_%s_combined.pdf", base_name, control_name))
    ggsave(
      filename = pdf_file,
      plot = p_combined,
      width = 10,
      height = 10,
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

  message("Generating nucleotide distribution paired comparison plots (mean of replicates)...")
  generate_all_plots(data)

  message("Done!")
}

# Run main function
main()
