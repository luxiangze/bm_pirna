# Transposon Scatter Plot Visualization
# Plots control_rpm vs treated_rpm with log2_fold_change color mapping
# Style: Nature journal inspired

library(here)
library(ggplot2)
library(scales)
library(dplyr)
library(ggrepel)

# Configuration
INPUT_FILE <- here("data/transposon_rpm/fc_SUGP1_vs_WT.tsv")
OUTPUT_DIR <- here("graphs/transposon_rpm")
LABEL_FILE <- NULL  # Set to here("path/to/labels.list") to enable custom labeling

# Thresholds
RPM_THRESHOLD <- 0.1
LOG2FC_THRESHOLD <- 1
TOP_N_LABELS <- 0  # Set to 0 to disable labels

# Create output directory if not exists
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load data
data <- read.delim(INPUT_FILE, sep = "\t", header = TRUE)

# Parse sample names from filename for title generation
# Expected format: fc_Treated_vs_Control.tsv
base_name <- basename(INPUT_FILE)
base_name <- gsub("^fc_", "", base_name)
base_name <- gsub("\\.tsv$", "", base_name)

# Extract treated and control groups from filename
parts <- strsplit(base_name, "_vs_")[[1]]
if (length(parts) == 2) {
  treated_name <- parts[1]
  control_name <- parts[2]
} else {
  # Fallback: use whole filename
  treated_name <- base_name
  control_name <- "Control"
}

# Generate title from parsed names
title_text <- sprintf("Expression level of piRNAs corresponding to transposons: %s vs %s", "SUGP1-KD", control_name)

# Categorize: significant (colored) vs not significant (gray)
data$significant <- !(data$control_rpm < RPM_THRESHOLD |
                       data$treated_rpm < RPM_THRESHOLD |
                       abs(data$log2_fold_change) < LOG2FC_THRESHOLD)

# Replace zeros with small value to avoid log10(-Inf)
# Find minimum positive value in both columns
pos_control <- data[data$control_rpm > 0, "control_rpm"]
pos_treated <- data[data$treated_rpm > 0, "treated_rpm"]

if (length(pos_control) > 0 && length(pos_treated) > 0) {
  min_val <- min(pos_control, pos_treated) / 2
} else {
  min_val <- 0.1  # Fallback pseudo-count
}

data$control_rpm_plot <- ifelse(data$control_rpm > 0, data$control_rpm, min_val)
data$treated_rpm_plot <- ifelse(data$treated_rpm > 0, data$treated_rpm, min_val)

# Calculate reasonable axis limits for placing the n= label
x_max <- max(data$control_rpm_plot, na.rm = TRUE)
y_max <- max(data$treated_rpm_plot, na.rm = TRUE)

# Determine max absolute log2fc for color scale symmetric limits
max_log2fc <- max(abs(data$log2_fold_change), na.rm = TRUE)
max_log2fc <- ceiling(max_log2fc)

# Determine labels to show on plot
labels_df <- NULL

if (!is.null(LABEL_FILE) && file.exists(LABEL_FILE)) {
  # Load custom labels from file (one transposon id per line)
  label_ids <- readLines(LABEL_FILE, warn = FALSE)
  label_ids <- trimws(label_ids)
  label_ids <- label_ids[label_ids != ""]
  labels_df <- data %>% filter(transposon_ID %in% label_ids)
} else {
  # Auto-label top up/down regulated significant transposons
  sig_data <- data %>% filter(significant)
  if (nrow(sig_data) > 0) {
    n_up <- min(TOP_N_LABELS, sum(sig_data$log2_fold_change > 0))
    n_down <- min(TOP_N_LABELS, sum(sig_data$log2_fold_change < 0))

    # Top upregulated (highest log2fc)
    up_top <- sig_data %>%
      filter(log2_fold_change > 0) %>%
      arrange(desc(log2_fold_change)) %>%
      slice_head(n = n_up)

    # Top downregulated (lowest log2fc)
    down_top <- sig_data %>%
      filter(log2_fold_change < 0) %>%
      arrange(log2_fold_change) %>%
      slice_head(n = n_down)

    labels_df <- bind_rows(up_top, down_top)
  }
}

# Custom color gradient: Deep Red-White-Deep Blue
color_gradient <- scale_color_gradient2(
  low = "#00008B",
  mid = "white",
  high = "#8B0000",
  limits = c(-5, 5),
  breaks = seq(-5, 5, by = 2.5),
  name = expression(log[2] * "FC")
)

# Calculate x_range for diagonal lines (y = 2x for log2fc=1, y = 0.5x for log2fc=-1)
# Use log scale compatible range
x_range <- c(min_val, x_max)

# Create the plot
p <- ggplot() +

  # Layer 1: Gray points (low count or |log2fc| < 1)
  geom_point(
    data = subset(data, !significant),
    aes(x = control_rpm_plot, y = treated_rpm_plot),
    color = "gray70",
    size = 3,
    shape = 19,
    stroke = 0
  ) +

  # Layer 2: Colored points (significant changes)
  geom_point(
    data = subset(data, significant),
    aes(x = control_rpm_plot, y = treated_rpm_plot, color = log2_fold_change),
    alpha = 1,
    size = 4,
    shape = 19,
    stroke = 0
  ) +

  # Log scale for both axes with log10 formatted labels (10^0, 10^1, 10^2, etc.)
  scale_x_log10(labels = function(x) sprintf("10^%.0f", log10(x)), expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_log10(labels = function(x) sprintf("10^%.0f", log10(x)), expand = expansion(mult = c(0.02, 0.05))) +

  # Log ticks
  annotation_logticks(sides = "bl", color = "grey40") +

  # Apply color gradient for significant points
  color_gradient +

  # Add diagonal reference lines for log2fc thresholds (log2fc = +-1)
  # y = 2x for log2fc=1, y = 0.5x for log2fc=-1
  annotate("line", x = x_range, y = x_range * 2,
           linetype = "dashed", color = "grey50", linewidth = 0.5) +
  annotate("line", x = x_range, y = x_range * 0.5,
           linetype = "dashed", color = "grey50", linewidth = 0.5) +

  # Add diagonal 1:1 line
  annotate("line", x = x_range, y = x_range,
           linetype = "solid", color = "grey30", linewidth = 0.5)

# Add labels if available
if (!is.null(labels_df) && nrow(labels_df) > 0) {
  p <- p +
    geom_text_repel(
      data = labels_df,
      aes(x = control_rpm_plot, y = treated_rpm_plot, label = transposon_ID),
      size = 6,
      color = "black",
      segment.color = "grey50",
      segment.alpha = 0.6,
      min.segment.length = 0.1,
      max.overlaps = Inf,
      box.padding = 0.3,
      force = 3  # Increase repulsion between labels (default=1)
    )
}

p <- p +

  # Labels
  labs(
    x = "Control RPM (log10)",
    y = "SUGP1-KD RPM (log10)"
  ) +

  theme_bw(base_size = 24) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    plot.title = element_text(size = 28, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 22, hjust = 0.5, color = "grey50"),
    legend.position = "right",
    legend.background = element_rect(color = NA, fill = "white"),
    legend.key.size = unit(1.5, "lines"),
    legend.text = element_text(size = 20),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 24, face = "bold"),
    plot.margin = margin(15, 10, 10, 10)
  ) +

  # n= added in top left corner (log scale compatible position)
  annotate("text", x = min_val, y = y_max * 1.3, label = sprintf("n = %d", nrow(data)),
           hjust = 0, vjust = 0, size = 8, fontface = "bold", color = "black")

# Calculate significant counts for labels
upregulated <- sum(data$significant & data$log2_fold_change > 0)
downregulated <- sum(data$significant & data$log2_fold_change < 0)

# Output filename
output_file <- file.path(OUTPUT_DIR, gsub("\\.tsv$", "_scatter.png", basename(INPUT_FILE)))

# Save PNG (300 dpi)
ggsave(output_file, plot = p, width = 8, height = 7, units = "in", dpi = 300)

# Save PDF
output_pdf <- gsub("\\.png$", ".pdf", output_file)
ggsave(output_pdf, plot = p, width = 8, height = 7, units = "in")

# Summary
cat("=== Summary Statistics ===\n")
cat("Total transposons:", nrow(data), "\n")
cat("Significant (|log2FC| >= 1):", sum(data$significant), "\n")
cat("Not significant (gray):", sum(!data$significant), "\n")
cat("  - Upregulated:", upregulated, "\n")
cat("  - Downregulated:", downregulated, "\n")

cat("\nFigure saved to:\n")
cat("  PNG:", output_file, "\n")
cat("  PDF:", output_pdf, "\n")
