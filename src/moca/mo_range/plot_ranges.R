#!/usr/bin/env Rscript

# plot_ranges.R
# Visualize motif positional distribution analysis from mo_range outputs
# Creates TSS and TTS proximity plots showing regulatory region preferences

library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)
library(cowplot)

# Custom color palette (matching your motif analysis colors)
custom_colors <- c("#AF91AF", "#804a5f", "#478fca", "#3e2c57")

# Set up paths - work from current directory structure
current_dir <- getwd()
cat("Current working directory:", current_dir, "\n")

# Locate output directory for mo_range results
if (basename(current_dir) == "mo_range") {
  # Running from mo_range directory - look in the actual output directory
  input_dir <- "../../../out/moca_results/mo_range"
  output_dir <- "../../../out/moca_results/mo_range"
} else if (file.exists("../out/moca_results/mo_range")) {
  # Running from src directory
  input_dir <- "../out/moca_results/mo_range"
  output_dir <- "../out/moca_results/mo_range"
} else if (file.exists("out/moca_results/mo_range")) {
  # Running from project root
  input_dir <- "out/moca_results/mo_range"
  output_dir <- "out/moca_results/mo_range"
} else {
  stop("Could not locate mo_range output directory. Please run from appropriate directory.")
}

cat("Looking for range data in:", input_dir, "\n")
cat("Outputs will be saved to:", output_dir, "\n")

###########################
# 1. Load range analysis data
###########################

# Look for seqlet pattern data
seqlet_files <- list.files(input_dir, pattern = ".*seqlet_pattern.*\\.csv$", full.names = TRUE)
tss_range_files <- list.files(input_dir, pattern = ".*TSS.*motif_ranges.*\\.csv$", full.names = TRUE)
tts_range_files <- list.files(input_dir, pattern = ".*TTS.*motif_ranges.*\\.csv$", full.names = TRUE)

cat("Found files:\n")
cat("- Seqlet files:", length(seqlet_files), "\n")
cat("- TSS range files:", length(tss_range_files), "\n")
cat("- TTS range files:", length(tts_range_files), "\n")

if (length(seqlet_files) == 0) {
  stop("No seqlet pattern files found. Expected pattern: *seqlet_pattern*.csv")
}

# Load the data
seqlet_data <- read_csv(seqlet_files[1], show_col_types = FALSE)
cat("Loaded seqlet data from:", basename(seqlet_files[1]), "\n")
cat("Data shape:", nrow(seqlet_data), "rows x", ncol(seqlet_data), "columns\n")

# Check required columns
required_cols <- c("start", "pattern")
missing_cols <- required_cols[!required_cols %in% colnames(seqlet_data)]
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Load range summary files if available
tss_ranges <- NULL
tts_ranges <- NULL

if (length(tss_range_files) > 0) {
  tss_ranges <- read_csv(tss_range_files[1], show_col_types = FALSE)
  cat("Loaded TSS range data from:", basename(tss_range_files[1]), "\n")
}

if (length(tts_range_files) > 0) {
  tts_ranges <- read_csv(tts_range_files[1], show_col_types = FALSE)
  cat("Loaded TTS range data from:", basename(tts_range_files[1]), "\n")
}

###########################
# 2. Prepare data for plotting
###########################

cat("\nAnalyzing positional distributions...\n")

# Define regions based on coordinate system
# Assuming 0-1500 is TSS-proximal, 1520+ is TTS-proximal
tss_threshold <- 1500
tts_start <- 1520
total_length <- 3020  # Total sequence length

# Filter data for TSS and TTS regions
tss_seqlets <- seqlet_data %>%
  filter(start <= tss_threshold)

tts_seqlets <- seqlet_data %>%
  filter(start > tts_start)

cat("TSS-proximal seqlets:", nrow(tss_seqlets), "\n")
cat("TTS-proximal seqlets:", nrow(tts_seqlets), "\n")

###########################
# 3. Create TSS-proximal distribution plot
###########################

create_tss_plot <- function(data) {
  # Create bins for histogram
  bin_width <- 100
  bins <- seq(0, tss_threshold + bin_width, by = bin_width)

  # Calculate histogram data
  hist_data <- data %>%
    mutate(
      bin = cut(start, breaks = bins, include.lowest = TRUE, right = FALSE),
      bin_center = as.numeric(gsub("\\[(.*),.*", "\\1", bin)) + bin_width/2
    ) %>%
    filter(!is.na(bin)) %>%
    count(bin_center, name = "count")

  # Create the plot
  p <- ggplot(hist_data, aes(x = bin_center, y = count)) +
    geom_bar(stat = "identity", width = bin_width * 0.8,
             fill = custom_colors[1], alpha = 0.7, color = "black", size = 0.3) +
    geom_vline(xintercept = 750, color = custom_colors[3], linestyle = "dashed",
               size = 1, alpha = 0.8) +
    labs(
      title = "Motif Distribution: TSS-Proximal Region",
      x = "Position (bp)",
      y = "Motif Count"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    ) +
    scale_x_continuous(breaks = seq(0, tss_threshold, 300)) +
    scale_y_continuous(breaks = seq(0, 1500, 500), limits = c(0, 1700)) +
    coord_cartesian(xlim = c(0, tss_threshold))

  # Region annotations removed to avoid overlap with bars

  return(p)
}

###########################
# 4. Create TTS-proximal distribution plot
###########################

create_tts_plot <- function(data) {
  if (nrow(data) == 0) {
    # Create empty plot if no data
    p <- ggplot() +
      annotate("text", x = 750, y = 0.5, label = "No TTS-proximal data available",
               size = 4, hjust = 0.5) +
      labs(
        title = "Motif Distribution: TTS-Proximal Region",
        x = "Position (bp)",
        y = "Motif Count"
      ) +
      theme_minimal() +
      coord_cartesian(xlim = c(0, 1500), ylim = c(0, 1))
    return(p)
  }

  # Transform coordinates for TTS analysis (3020 - position)
  tts_transformed <- data %>%
    mutate(tts_pos = total_length - start) %>%
    filter(tts_pos >= 0 & tts_pos <= 1500)  # Keep reasonable range

  if (nrow(tts_transformed) == 0) {
    # Create empty plot if no data after transformation
    p <- ggplot() +
      annotate("text", x = 750, y = 0.5, label = "No valid TTS data after transformation",
               size = 4, hjust = 0.5) +
      labs(
        title = "Motif Distribution: TTS-Proximal Region",
        x = "Position (bp)",
        y = "Motif Count"
      ) +
      theme_minimal() +
      coord_cartesian(xlim = c(0, 1500), ylim = c(0, 1))
    return(p)
  }

  # Create bins for histogram
  bin_width <- 100
  bins <- seq(0, 1500 + bin_width, by = bin_width)

  # Calculate histogram data
  hist_data <- tts_transformed %>%
    mutate(
      bin = cut(tts_pos, breaks = bins, include.lowest = TRUE, right = FALSE),
      bin_center = as.numeric(gsub("\\[(.*),.*", "\\1", bin)) + bin_width/2
    ) %>%
    filter(!is.na(bin)) %>%
    count(bin_center, name = "count")

  # Create the plot
  p <- ggplot(hist_data, aes(x = bin_center, y = count)) +
    geom_bar(stat = "identity", width = bin_width * 0.8,
             fill = custom_colors[2], alpha = 0.7, color = "black", size = 0.3) +
    geom_vline(xintercept = 750, color = custom_colors[4], linestyle = "dashed",
               size = 1, alpha = 0.8) +
    labs(
      title = "Motif Distribution: TTS-Proximal Region",
      x = "Position (bp)",
      y = "Motif Count"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    ) +
    scale_x_continuous(breaks = seq(0, 1500, 300)) +
    scale_y_continuous(breaks = seq(0, 1500, 500), limits = c(0, 1700)) +
    coord_cartesian(xlim = c(0, 1500))

  # Region annotations removed to avoid overlap with bars

  return(p)
}

###########################
# 5. Generate plots
###########################

cat("\nGenerating plots...\n")

# Create individual plots
p1 <- create_tss_plot(tss_seqlets)
p2 <- create_tts_plot(tts_seqlets)

# Combine plots
combined_plot <- plot_grid(p1, p2, ncol = 2, align = "h", axis = "tb")

# Add overall title
title <- ggdraw() +
  draw_label("Motif Positional Distribution Analysis",
             fontface = 'bold', size = 14)

final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

###########################
# 6. Save outputs
###########################

# Save the combined plot
output_file <- file.path(output_dir, "motif_positional_distribution.png")
ggsave(output_file, final_plot, width = 12, height = 6, dpi = 300)
cat("Combined plot saved to:", output_file, "\n")

# Save individual plots
tss_output <- file.path(output_dir, "tss_proximal_distribution.png")
ggsave(tss_output, p1, width = 6, height = 4, dpi = 300)
cat("TSS plot saved to:", tss_output, "\n")

tts_output <- file.path(output_dir, "tts_proximal_distribution.png")
ggsave(tts_output, p2, width = 6, height = 4, dpi = 300)
cat("TTS plot saved to:", tts_output, "\n")

###########################
# 7. Generate summary statistics
###########################

cat("\n=== Summary Statistics ===\n")
cat("Total seqlets analyzed:", nrow(seqlet_data), "\n")
cat("TSS-proximal seqlets (≤", tss_threshold, "bp):", nrow(tss_seqlets),
    sprintf(" (%.1f%%)", nrow(tss_seqlets)/nrow(seqlet_data)*100), "\n")
cat("TTS-proximal seqlets (>", tts_start, "bp):", nrow(tts_seqlets),
    sprintf(" (%.1f%%)", nrow(tts_seqlets)/nrow(seqlet_data)*100), "\n")

# Pattern distribution
if ("pattern" %in% colnames(seqlet_data)) {
  pattern_summary <- seqlet_data %>%
    count(pattern, sort = TRUE)

  cat("\nTop 5 most frequent patterns:\n")
  print(head(pattern_summary, 5))
}

# Position statistics
pos_stats <- seqlet_data %>%
  summarise(
    min_pos = min(start, na.rm = TRUE),
    max_pos = max(start, na.rm = TRUE),
    mean_pos = round(mean(start, na.rm = TRUE), 1),
    median_pos = round(median(start, na.rm = TRUE), 1)
  )

cat("\nPosition statistics:\n")
cat("Range:", pos_stats$min_pos, "-", pos_stats$max_pos, "bp\n")
cat("Mean position:", pos_stats$mean_pos, "bp\n")
cat("Median position:", pos_stats$median_pos, "bp\n")

###########################
# 8. Save summary data
###########################

# Create summary table
summary_table <- data.frame(
  Region = c("TSS-proximal", "TTS-proximal", "Total"),
  Count = c(nrow(tss_seqlets), nrow(tts_seqlets), nrow(seqlet_data)),
  Percentage = c(
    round(nrow(tss_seqlets)/nrow(seqlet_data)*100, 1),
    round(nrow(tts_seqlets)/nrow(seqlet_data)*100, 1),
    100.0
  )
)

summary_file <- file.path(output_dir, "positional_distribution_summary.csv")
write_csv(summary_table, summary_file)
cat("\nSummary table saved to:", summary_file, "\n")

###########################
# 9. Interpretation notes
###########################

cat("\n=== Interpretation Notes ===\n")
cat("TSS-proximal (left panel): Shows motif distribution around transcription start sites\n")
cat("TTS-proximal (right panel): Shows motif distribution around transcription termination sites\n")
cat("Vertical dashed lines: Approximate TSS/TTS positions\n")
cat("Clustering patterns: Indicate preferred regulatory regions\n\n")

cat("Biological Significance:\n")
cat("- Promoter-proximal clustering suggests transcriptional initiation control\n")
cat("- Downstream patterns may indicate post-transcriptional regulation\n")
cat("- Distribution asymmetry reveals directional regulatory preferences\n\n")

cat("Analysis complete!\n")