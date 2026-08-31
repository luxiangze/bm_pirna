#' Plot TE RPM Visualization - Paired Comparison
#'
#' This script generates paired comparison plots for transposon RPM at
#' top-level class and sub-class resolution. Each pair from sample_map.csv
#' produces grouped bar, Control-normalized bar, and heatmap figures.
#'
#' @author Yongkang Guo
#' @date 2026-08-23

# Load required packages
library(here)
library(tidyverse)
library(ggplot2)
library(scales)

# Define file paths
DATA_DIR <- here()
TE_RPM_DIR <- file.path(DATA_DIR, "data/processed/sugp1_smRNA-seq_20260822/transposon_rpm")
SAMPLE_MAP_FILE <- file.path(DATA_DIR, "data/external/sample_map.csv")
OUTPUT_DIR <- file.path(DATA_DIR, "data/processed/sugp1_smRNA-seq_20260822/plots")

# Keep rows where max RPM across compared samples is at least this value
RPM_THRESHOLD <- 1000

#' Build the condition color palette for one comparison
#'
#' A single control keeps the legacy "Control"/"Treated" labels and the
#' blue/red palette. When the sample map lists several controls separated by
#' semicolons, each control receives its own color and the real group names
#' are used as legend/facet labels so every control gets its own bar/facet.
#'
#' @param control_groups Character vector of control group names
#' @param treated_group Treated group name (scalar)
#' @return Named character vector mapping group label -> color
create_cond_colors <- function(control_groups, treated_group) {
  group_names <- c(control_groups, treated_group)
  group_labels <- if (length(control_groups) == 1L) {
    c("Control", "Treated")
  } else {
    group_names
  }

  base_control_colors <- c("#457B9D", "#F4A261", "#2A9D8F", "#6A4C93", "#E9C46A")
  if (length(control_groups) <= length(base_control_colors)) {
    control_colors <- base_control_colors[seq_along(control_groups)]
  } else {
    control_colors <- grDevices::hcl.colors(
      length(control_groups), palette = "Dark 3"
    )
  }

  setNames(c(control_colors, "#E63946"), group_labels)
}

#' Parse group names from a sample map cell
#'
#' Control_sample may contain multiple group names separated by semicolons
#' (e.g. "control;p200"). Treated_sample must contain exactly one group.
#'
#' @param value Scalar sample map value
#' @param column_name Column name used in validation messages
#' @param allow_multiple Whether semicolon-separated values are allowed
#' @return Character vector of validated group names
parse_group_names <- function(value, column_name, allow_multiple = TRUE) {
  if (length(value) != 1L || is.na(value)) {
    stop(column_name, " must contain a non-missing group name.")
  }

  groups <- str_split(as.character(value), fixed(";"))[[1]] %>%
    str_trim()

  if (length(groups) == 0L || any(groups == "")) {
    stop(column_name, " contains an empty group name.")
  }
  if (!allow_multiple && length(groups) != 1L) {
    stop(column_name, " must contain exactly one group name.")
  }

  duplicate_groups <- unique(groups[duplicated(groups)])
  if (length(duplicate_groups) > 0) {
    stop(
      column_name,
      " contains duplicate group name(s): ",
      paste(duplicate_groups, collapse = ", ")
    )
  }

  groups
}

#' Validate that comparison groups exist in the available sample columns
#'
#' @param control_groups Character vector of control group names
#' @param treated_group Treated group name (scalar)
#' @param available_groups Character vector of groups present in the data
validate_comparison_groups <- function(
    control_groups,
    treated_group,
    available_groups
) {
  if (treated_group %in% control_groups) {
    stop("Treated group cannot also be listed as a control group: ", treated_group)
  }

  missing_groups <- setdiff(c(control_groups, treated_group), available_groups)
  if (length(missing_groups) > 0) {
    stop(
      "Sample map group(s) not found in TE RPM sample columns: ",
      paste(missing_groups, collapse = ", ")
    )
  }
}

#' Build a filesystem-safe slug from one or more group names
#'
#' Multiple group names are joined with "_and_" so that a combined control
#' such as "control;p200" produces "control_and_p200" rather than a filename
#' containing a semicolon.
#'
#' @param groups Character vector of group names
#' @return Lowercase slug string safe for use in filenames
make_group_slug <- function(groups) {
  slug_parts <- groups %>%
    str_to_lower() %>%
    str_replace_all("[^[:alnum:]_.-]+", "_") %>%
    str_replace_all("^_+|_+$", "")

  if (any(slug_parts == "")) {
    stop("Cannot create an output filename from an empty group name.")
  }

  paste(slug_parts, collapse = "_and_")
}

# Create output directory if not exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

#' Load data files
#'
#' @return List containing loaded data frames
load_data <- function() {
  top_class <- read_tsv(
    file.path(TE_RPM_DIR, "te_top_class_rpm_matrix.tsv"),
    col_types = cols(te_class = col_character(), .default = col_double())
  )

  sub_class <- read_tsv(
    file.path(TE_RPM_DIR, "te_sub_class_rpm_matrix.tsv"),
    col_types = cols(te_class = col_character(), .default = col_double())
  )

  sample_map <- read_csv(SAMPLE_MAP_FILE, col_types = cols())
  required_map_columns <- c("Control_sample", "Treated_sample")
  missing_map_columns <- setdiff(required_map_columns, names(sample_map))
  if (length(missing_map_columns) > 0) {
    stop(
      "Missing required sample map column(s): ",
      paste(missing_map_columns, collapse = ", ")
    )
  }

  list(
    top_class = top_class,
    sub_class = sub_class,
    sample_map = sample_map
  )
}

#' Sample columns matching one or more group names
#'
#' @param df RPM matrix with te_class plus sample columns
#' @param groups Character vector of group names (e.g. "control", "sugp1")
#' @return Character vector of matching sample column names
get_sample_cols <- function(df, groups) {
  pattern <- paste0("^(", paste(str_escape(groups), collapse = "|"), ")_rep\\d+$")
  names(df)[str_detect(names(df), pattern)]
}

#' Discard rows where all compared samples fall below RPM_THRESHOLD
#'
#' @param df RPM matrix
#' @param sample_cols Sample columns to consider
#' @param threshold Minimum max RPM to retain a row
#' @return Filtered data frame
filter_by_max_rpm <- function(df, sample_cols, threshold = RPM_THRESHOLD) {
  df %>%
    filter(pmax(!!!syms(sample_cols), na.rm = TRUE) >= threshold)
}

#' Convert a RPM matrix subset to long format for one group pair
#'
#' @param df Filtered RPM matrix
#' @param sample_cols Sample columns to pivot
#' @param control_groups Control group name(s); may be a vector when the
#'   sample map lists multiple controls separated by semicolons
#' @param treated_group Treated group name
#' @return Long data frame with group and replicate labels
to_long <- function(df, sample_cols, control_groups, treated_group) {
  group_names <- c(control_groups, treated_group)
  # Single control keeps the legacy "Control"/"Treated" labels; multiple
  # controls fall back to the real group names so each control is plotted
  # as its own bar / facet.
  group_labels <- if (length(control_groups) == 1L) {
    c("Control", "Treated")
  } else {
    group_names
  }

  df %>%
    pivot_longer(all_of(sample_cols), names_to = "sample", values_to = "rpm") %>%
    mutate(
      group_name = str_replace(sample, "_rep\\d+$", ""),
      group = factor(group_name, levels = group_names, labels = group_labels),
      rep_label = factor(
        str_extract(sample, "rep\\d+$"),
        levels = unique(str_extract(sample_cols, "rep\\d+$"))
      )
    )
}

#' Shared bar-chart theme
#'
#' @param base_size Base font size
#' @param x_angle X-axis label rotation
#' @return ggplot2 theme
bar_theme <- function(base_size = 20, x_angle = 0) {
  x_hjust <- if (x_angle > 0) 1 else 0.5
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position    = "top",
      legend.text        = element_text(size = base_size - 4),
      axis.text.x        = element_text(
        size = base_size - 5, face = "bold", angle = x_angle, hjust = x_hjust
      ),
      axis.text.y        = element_text(size = base_size - 6),
      axis.title.y       = element_text(size = base_size - 4, face = "bold"),
      plot.margin        = margin(12, 10, 8, 10)
    )
}

#' Run a two-group Wilcoxon rank-sum / Mann-Whitney U test
#'
#' @param x Numeric values for the first group
#' @param y Numeric values for the second group
#' @return A tibble containing the P value and test method or status
run_two_group_test <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < 2 || length(y) < 2) {
    return(tibble(P_value = NA_real_, Method = "insufficient replicates"))
  }

  test_result <- tryCatch(
    wilcox.test(x, y),
    error = function(e) NULL
  )

  if (is.null(test_result)) {
    return(tibble(P_value = NA_real_, Method = "test failed"))
  }

  tibble(P_value = test_result$p.value, Method = "Wilcoxon rank-sum")
}

#' Print per-class significance results for a group comparison
#'
#' When multiple controls are supplied, each control is tested against the
#' treated group independently (controls are never pooled).
#'
#' @param long_df Long-format RPM data for one pair
#' @param cond_colors Named color vector whose names are the group labels
#'   (control label(s) first, treated label last)
print_significance_per_class <- function(long_df, cond_colors) {
  group_labels <- names(cond_colors)
  control_labels <- head(group_labels, -1)
  treated_label <- tail(group_labels, 1)

  for (ctrl in control_labels) {
    significance_results <- long_df %>%
      group_by(te_class) %>%
      summarise(
        control_values = list(rpm[group == ctrl]),
        treated_values = list(rpm[group == treated_label]),
        .groups = "drop"
      ) %>%
      rowwise() %>%
      mutate(test_result = list(run_two_group_test(control_values, treated_values))) %>%
      unnest_wider(test_result) %>%
      mutate(
        Significance = case_when(
          is.na(P_value) ~ "NA",
          P_value < 0.001 ~ "***",
          P_value < 0.01 ~ "**",
          P_value < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      ) %>%
      select(te_class, P_value, Significance, Method)

    message(sprintf(
      "Significance per class: %s vs %s", ctrl, treated_label
    ))
    print(significance_results)
    message("Significance levels: *** p < 0.001, ** p < 0.01, * p < 0.05, ns not significant, NA not testable.")
  }
}

#' Grouped bar chart of mean RPM with replicate points
#'
#' @param long_df Long-format RPM data
#' @param cond_colors Named color vector (control label(s) first, treated last)
#' @param y_scale_k If TRUE, y-axis labels use K units
#' @param x_angle X-axis label rotation
#' @return ggplot object
plot_grouped_bar <- function(long_df, cond_colors, y_scale_k = TRUE, x_angle = 0) {
  control_label <- names(cond_colors)[1]

  mean_df <- long_df %>%
    group_by(te_class, group) %>%
    summarise(mean_rpm = mean(rpm), .groups = "drop")

  # Order classes by the first control's mean RPM (descending)
  class_order <- mean_df %>%
    filter(group == control_label) %>%
    arrange(desc(mean_rpm)) %>%
    pull(te_class)

  mean_df <- mutate(mean_df, te_class = factor(te_class, levels = class_order))
  long_df <- mutate(long_df, te_class = factor(te_class, levels = class_order))

  p <- ggplot() +
    geom_col(
      data     = mean_df,
      aes(x    = te_class, y = mean_rpm, fill = group),
      position = position_dodge(width = 0.75),
      width    = 0.65,
      alpha    = 0.92
    ) +
    geom_point(
      data     = long_df,
      aes(x    = te_class, y = rpm, fill = group),
      position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.07),
      size     = 3,
      shape    = 21,
      color    = "white",
      stroke   = 0.6,
      alpha    = 0.95,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = cond_colors) +
    labs(x = NULL, y = "Mean RPM", fill = NULL) +
    bar_theme(x_angle = x_angle)

  if (y_scale_k) {
    p <- p + scale_y_continuous(
      labels = label_number(scale = 1e-3, suffix = "K"),
      expand = expansion(mult = c(0, 0.06))
    )
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.06)))
  }

  p
}

#' Normalized bar chart (first Control mean = 1 per class)
#'
#' The first control acts as the reference. When multiple controls are
#' supplied, every other control is also expressed relative to that same
#' reference so all bars share one consistent baseline.
#'
#' @param long_df Long-format RPM data
#' @param cond_colors Named color vector (control label(s) first, treated last)
#' @param x_angle X-axis label rotation
#' @param point_size Replicate point size
#' @param base_size Theme base size
#' @return ggplot object
plot_normalized_bar <- function(long_df, cond_colors, x_angle = 0, point_size = 3, base_size = 20) {
  control_label <- names(cond_colors)[1]

  mean_df <- long_df %>%
    group_by(te_class, group) %>%
    summarise(mean_rpm = mean(rpm), .groups = "drop")

  control_ref <- mean_df %>%
    filter(group == control_label) %>%
    select(te_class, control_mean = mean_rpm)

  class_order <- control_ref %>%
    arrange(desc(control_mean)) %>%
    pull(te_class)

  mean_norm <- mean_df %>%
    left_join(control_ref, by = "te_class") %>%
    mutate(
      rel_rpm  = mean_rpm / control_mean,
      te_class = factor(te_class, levels = class_order)
    )

  long_norm <- long_df %>%
    left_join(control_ref, by = "te_class") %>%
    mutate(
      rel_rpm  = rpm / control_mean,
      te_class = factor(te_class, levels = class_order)
    )

  ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_col(
      data     = mean_norm,
      aes(x    = te_class, y = rel_rpm, fill = group),
      position = position_dodge(width = 0.75),
      width    = 0.65,
      alpha    = 0.92
    ) +
    geom_point(
      data     = long_norm,
      aes(x    = te_class, y = rel_rpm, fill = group),
      position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.07),
      size     = point_size,
      shape    = 21,
      color    = "white",
      stroke   = 0.6,
      alpha    = 0.95,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = cond_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(x = NULL, y = "Relative RPM (Control = 1)", fill = NULL) +
    bar_theme(base_size = base_size, x_angle = x_angle)
}

#' Sub-class heatmap faceted by group (control(s) / treated)
#'
#' @param long_df Long-format RPM data
#' @param cond_colors Named color vector used only to determine the reference
#'   control label (the first name) for row ordering
#' @return ggplot object
plot_heatmap <- function(long_df, cond_colors) {
  control_label <- names(cond_colors)[1]

  rpm_vals <- long_df$rpm
  pos_vals <- rpm_vals[rpm_vals > 0]
  log_mid <- if (length(pos_vals) > 0) {
    10^(mean(log10(pos_vals)))
  } else {
    1
  }

  class_order <- long_df %>%
    filter(group == control_label) %>%
    group_by(te_class) %>%
    summarise(control_mean = mean(rpm), .groups = "drop") %>%
    arrange(control_mean) %>%
    pull(te_class)

  plot_data <- long_df %>%
    mutate(
      te_class   = factor(te_class, levels = class_order),
      text_color = if_else(rpm > log_mid, "white", "grey15"),
      tile_label = sprintf("%.1fK", rpm / 1000)
    )

  ggplot(plot_data, aes(x = rep_label, y = te_class, fill = rpm)) +
    geom_tile(color = "white", linewidth = 0.9) +
    geom_text(aes(label = tile_label, color = text_color), size = 3.5) +
    scale_color_identity(guide = "none") +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
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
      legend.key.height = grid::unit(2.5, "lines"),
      legend.title      = element_text(size = 13),
      legend.text       = element_text(size = 12),
      panel.spacing     = grid::unit(0.8, "lines"),
      plot.margin       = margin(10, 15, 10, 10)
    )
}

#' Save a plot as PNG (300 dpi) and PDF
#'
#' @param plot ggplot object
#' @param filename File stem without extension
#' @param width Width in inches
#' @param height Height in inches
save_plot <- function(plot, filename, width, height) {
  png_file <- file.path(OUTPUT_DIR, paste0(filename, ".png"))
  pdf_file <- file.path(OUTPUT_DIR, paste0(filename, ".pdf"))

  ggsave(png_file, plot = plot, width = width, height = height, units = "in", dpi = 300)
  message(sprintf("  Saved PNG: %s", png_file))

  ggsave(
    pdf_file, plot = plot, width = width, height = height, units = "in", device = pdf
  )
  message(sprintf("  Saved PDF: %s", pdf_file))
}

#' Prepare long-format data for one matrix and one group pair
#'
#' @param matrix_df Top-class or sub-class RPM matrix
#' @param control_groups Control group name(s); may be a vector
#' @param treated_group Treated group name
#' @param label Label used in log messages
#' @return Long data frame, or NULL if no rows / samples remain
prepare_pair_data <- function(matrix_df, control_groups, treated_group, label) {
  sample_cols <- get_sample_cols(matrix_df, c(control_groups, treated_group))
  if (length(sample_cols) == 0) {
    message(sprintf("  No sample columns for %s vs %s in %s, skip",
                    paste(control_groups, collapse = ", "), treated_group, label))
    return(NULL)
  }

  filtered <- matrix_df %>%
    filter(te_class != "Unknown") %>%
    filter_by_max_rpm(sample_cols)

  message(sprintf(
    "  %s: retained %d / %d rows (max RPM >= %s)",
    label, nrow(filtered), nrow(matrix_df), RPM_THRESHOLD
  ))

  if (nrow(filtered) == 0) {
    message(sprintf("  No %s rows pass RPM threshold, skip", label))
    return(NULL)
  }

  to_long(filtered, sample_cols, control_groups, treated_group)
}

#' Generate all plots for one Control / Treated pair
#'
#' @param data Data list from load_data()
#' @param control_groups Control group name(s); may be a vector
#' @param treated_group Treated group name
generate_pair_plots <- function(data, control_groups, treated_group) {
  cond_colors <- create_cond_colors(control_groups, treated_group)
  treated_slug <- make_group_slug(treated_group)
  control_slug <- make_group_slug(control_groups)
  base_name <- sprintf("%s_vs_%s", treated_slug, control_slug)

  control_label <- paste(control_groups, collapse = ", ")
  message(sprintf("Generating plots for group: %s vs %s", control_label, treated_group))

  top_long <- prepare_pair_data(
    data$top_class, control_groups, treated_group, "te_top_class"
  )
  if (!is.null(top_long)) {
    print_significance_per_class(top_long, cond_colors)
    save_plot(
      plot_grouped_bar(top_long, cond_colors),
      sprintf("te_top_class_grouped_bar_%s", base_name),
      width = 8, height = 6
    )
    save_plot(
      plot_normalized_bar(top_long, cond_colors),
      sprintf("te_top_class_normalized_bar_%s", base_name),
      width = 8, height = 6
    )
  }

  sub_long <- prepare_pair_data(
    data$sub_class, control_groups, treated_group, "te_sub_class"
  )
  if (!is.null(sub_long)) {
    print_significance_per_class(sub_long, cond_colors)
    save_plot(
      plot_normalized_bar(sub_long, cond_colors, x_angle = 45, point_size = 2.5, base_size = 16),
      sprintf("te_sub_class_normalized_bar_%s", base_name),
      width = 12, height = 6
    )
    save_plot(
      plot_heatmap(sub_long, cond_colors),
      sprintf("te_sub_class_heatmap_%s", base_name),
      width = 8, height = 5
    )
  }
}

#' Generate all group comparison plots
#'
#' @param data Data list from load_data()
generate_all_plots <- function(data) {
  # Groups available across both RPM matrices (replicate suffix stripped)
  available_groups <- unique(str_replace(
    names(data$top_class)[-1], "_rep\\d+$", ""
  ))

  for (i in seq_len(nrow(data$sample_map))) {
    control_groups <- parse_group_names(
      data$sample_map$Control_sample[i],
      sprintf("Control_sample in row %d", i)
    )
    treated_group <- parse_group_names(
      data$sample_map$Treated_sample[i],
      sprintf("Treated_sample in row %d", i),
      allow_multiple = FALSE
    )[[1]]

    validate_comparison_groups(control_groups, treated_group, available_groups)

    generate_pair_plots(data, control_groups, treated_group)
  }

  message("All plots saved to ", OUTPUT_DIR)
}

# Main execution
main <- function() {
  message("Loading data...")
  data <- load_data()

  message("Processing TE RPM visualization...")
  generate_all_plots(data)

  message("Done!")
}

# Run main function
main()
