library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(readr)

visualize_motif_positions <- function(csv_file, output_dir = "vitis_cre/src/results/moca_blue/mo_range/summary") {
  #' Visualize motif mean positions across TSS and TTS regions
  #' 
  #' @param csv_file 
  #' @param output_dir 
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set up logging
  log_file <- file.path(output_dir, "summary.log")
  log_con <- file(log_file, "w")
  
  # Function to write to both console and log
  log_cat <- function(...) {
    message <- paste(...)
    cat(message)
    cat(message, file = log_con)
  }
  
  # Custom colors
  CUSTOM_COLORS <- c("#AF91AF", "#804a5f", "#478fca", "#3e2c57")
  
  # Read data
  df <- read_csv(csv_file, show_col_types = FALSE)
  
  # Check column names and adjust if needed
  if ("mean_start" %in% names(df) && "mean_end" %in% names(df)) {
    # Using seqlet summary file
    start_col <- "mean_start"
    end_col <- "mean_end"
  } else if ("start" %in% names(df) && "end" %in% names(df)) {
    # Using raw seqlet file
    start_col <- "start"
    end_col <- "end"
  } else {
    stop("Expected columns 'mean_start'/'mean_end' not found in data")
  }
  
  # Print data summary for inspection
  log_cat("=== MOTIF POSITION ANALYSIS ===\n")
  log_cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  log_cat("Input file:", csv_file, "\n")
  log_cat("Output directory:", output_dir, "\n\n")
  
  log_cat("Data Summary:\n")
  log_cat("Position range:", round(min(df[[start_col]])), "-", round(max(df[[end_col]])), "\n")
  log_cat("Metaclusters:", paste(unique(df$metacluster), collapse = ", "), "\n")
  log_cat("Total patterns:", nrow(df), "\n\n")
  
  # Define coordinate system based on input sequence architecture
  # Using values: TSS region (1000+500=1500bp) + Buffer (20bp) + TTS region (500+1000=1500bp)
  TSS_REGION_LENGTH <- 1000 + 500  # 1000bp upstream + 500bp downstream of TSS
  BUFFER_LENGTH <- 20  # 20bp N buffer
  TTS_REGION_LENGTH <- 500 + 1000  # 500bp upstream + 1000bp downstream of TTS
  
  TSS_END <- TSS_REGION_LENGTH  # End of TSS region (position 1500)
  BUFFER_END <- TSS_END + BUFFER_LENGTH  # End of buffer (position 1520) 
  TTS_START <- BUFFER_END  # Start of TTS region (position 1520)
  TOTAL_LENGTH <- TSS_REGION_LENGTH + BUFFER_LENGTH + TTS_REGION_LENGTH  # Total: 3020bp
  
  # Create position bins for visualization
  bin_size <- 200  # 200bp bins
  bins <- seq(0, TOTAL_LENGTH + bin_size, bin_size)
  
  # Calculate mean position and assign bins
  df <- df %>%
    mutate(
      mean_position = (.data[[start_col]] + .data[[end_col]]) / 2,
      motif_length = .data[[end_col]] - .data[[start_col]],
      expression_type = case_when(
        metacluster == "metacluster_0" ~ "High Expression",
        metacluster == "metacluster_1" ~ "Low Expression",
        TRUE ~ "Unknown"
      ),
      pattern_short = gsub("pattern_", "P", pattern),
      region = ifelse(mean_position < TSS_END, "TSS", "TTS")
    )
  
  # Create bin centers vector
  bin_centers <- bins[-length(bins)] + bin_size/2
  
  # Assign bins and bin centers
  df$bin <- cut(df$mean_position, breaks = bins, include.lowest = TRUE)
  df$bin_center <- bin_centers[as.numeric(df$bin)]
  
  # Remove rows with NA bin centers
  df <- df %>% filter(!is.na(bin_center))
  
  # Color mapping
  color_map <- c("High Expression" = CUSTOM_COLORS[1], 
                 "Low Expression" = CUSTOM_COLORS[2])
  
  # Create binned data for histogram
  bin_data <- df %>%
    filter(!is.na(bin_center)) %>%
    group_by(bin_center, expression_type) %>%
    summarise(
      count = n(),
      patterns = paste(pattern_short, collapse = ", "),
      .groups = "drop"
    )
  
  # Plot 1: Histogram of motif positions
  p1 <- ggplot(bin_data, aes(x = bin_center, y = count, fill = expression_type)) +
    geom_col(position = "dodge", alpha = 0.7, width = bin_size * 0.8, 
             color = "white", linewidth = 0.5) +
    geom_text(aes(label = patterns), 
              position = position_dodge(width = bin_size * 0.8),
              angle = 45, hjust = 0, vjust = -0.2, size = 2.5) +
    geom_vline(xintercept = TSS_END, color = "#8B0000", linetype = "dashed", 
               linewidth = 1.2, alpha = 0.9) +
    geom_vline(xintercept = TTS_START, color = "#8B0000", linetype = "dashed", 
               linewidth = 1.2, alpha = 0.9) +
    annotate("text", x = TSS_END/2, y = max(bin_data$count) * 0.9, 
             label = "TSS Region", hjust = 0.5, vjust = 0.5, size = 4, 
             fontface = "bold", color = CUSTOM_COLORS[3]) +
    annotate("rect", xmin = 0, xmax = TSS_END/2 + 200, 
             ymin = max(bin_data$count) * 0.85, ymax = max(bin_data$count) * 0.95,
             fill = CUSTOM_COLORS[3], alpha = 0.2) +
    annotate("text", x = TSS_END + (TOTAL_LENGTH - TSS_END)/2, 
             y = max(bin_data$count) * 0.9, 
             label = "TTS Region", hjust = 0.5, vjust = 0.5, size = 4, 
             fontface = "bold", color = CUSTOM_COLORS[4]) +
    annotate("rect", xmin = TSS_END + (TOTAL_LENGTH - TSS_END)/2 - 200, 
             xmax = TOTAL_LENGTH, 
             ymin = max(bin_data$count) * 0.85, ymax = max(bin_data$count) * 0.95,
             fill = CUSTOM_COLORS[4], alpha = 0.2) +
    scale_fill_manual(values = color_map, name = "Expression Association") +
    labs(
      title = "Distribution of Motif Patterns Across TSS and TTS Regions",
      x = "Genomic Position (bp)",
      y = "Number of Motifs"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  # Plot 2: Individual motif positions as horizontal bars
  df_ordered <- df %>%
    arrange(metacluster, desc(mean_position)) %>%
    mutate(y_pos = row_number())
  
  p2 <- ggplot(df_ordered, aes(y = y_pos)) +
    geom_segment(aes(x = mean_start, xend = mean_end, 
                     color = expression_type), 
                 linewidth = 3, alpha = 0.7) +
    geom_text(aes(x = mean_position, label = pattern_short), 
              size = 2.5, fontface = "bold", color = "white") +
    geom_vline(xintercept = TSS_END, color = "red", linetype = "dashed", 
               linewidth = 1, alpha = 0.8) +
    scale_color_manual(values = color_map, name = "Expression Association") +
    scale_y_continuous(
      breaks = df_ordered$y_pos,
      labels = paste0(df_ordered$pattern_short, " (", 
                     gsub("metacluster_", "MC", df_ordered$metacluster), ")"),
      expand = c(0.01, 0.01)
    ) +
    labs(
      title = "Individual Motif Positions and Lengths",
      x = "Genomic Position (bp)",
      y = "Individual Patterns"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Combine plots using patchwork
  combined_plot <- p1 / p2 + plot_layout(heights = c(2, 1))
  
  # Save plot
  output_file <- file.path(output_dir, "motif_position_distribution.png")
  ggsave(output_file, combined_plot, width = 14, height = 10, dpi = 300, bg = "white")
  log_cat("Plot saved to:", output_file, "\n\n")
  
  # Print summary statistics
  log_cat("Summary Statistics:\n")
  log_cat(paste(rep("=", 50), collapse = ""), "\n")
  
  summary_stats <- df %>%
    group_by(expression_type, metacluster) %>%
    summarise(
      n_patterns = n(),
      min_pos = min(.data[[start_col]]),
      max_pos = max(.data[[end_col]]),
      avg_length = mean(motif_length),
      tss_motifs = sum(region == "TSS"),
      tts_motifs = sum(region == "TTS"),
      .groups = "drop"
    )
  
  for (i in 1:nrow(summary_stats)) {
    row <- summary_stats[i, ]
    log_cat("\n", row$expression_type, " (", row$metacluster, "):\n", sep = "")
    log_cat("  Number of patterns:", row$n_patterns, "\n")
    log_cat("  Mean position range:", round(row$min_pos), "-", round(row$max_pos), "\n")
    log_cat("  Average motif length:", round(row$avg_length, 1), "bp\n")
    log_cat("  TSS region motifs:", row$tss_motifs, "\n")
    log_cat("  TTS region motifs:", row$tts_motifs, "\n")
  }
  
  # Close log file
  close(log_con)
  cat("Analysis complete. Log saved to:", log_file, "\n")
  
  # Return the data for further analysis
  return(invisible(list(
    data = df,
    plot = combined_plot,
    summary = summary_stats,
    output_dir = output_dir,
    files_created = c(output_file, log_file)
  )))
}

# Function to analyze positional preferences
analyze_positional_bias <- function(df, tss_end = 1500, output_dir = "results/moca_blue/mo_range") {
  #' Analyze whether motifs show positional bias
  #' 
  #' @param df Data frame with motif position data
  #' @param tss_end Position marking end of TSS region (default: 1500bp)
  #' @param output_dir Directory to save outputs
  
  # Set up logging
  log_file <- file.path(output_dir, "positional_bias_analysis.log")
  log_con <- file(log_file, "w")
  
  # Function to write to both console and log
  log_cat <- function(...) {
    message <- paste(...)
    cat(message)
    cat(message, file = log_con)
  }
  
  bias_analysis <- df %>%
    mutate(region = ifelse(mean_position < tss_end, "TSS", "TTS")) %>%
    group_by(pattern, expression_type, metacluster) %>%
    summarise(
      total_seqlets = first(seqlet_count),
      mean_pos = first(mean_position),
      region = first(region),
      .groups = "drop"
    ) %>%
    arrange(desc(total_seqlets))
  
  log_cat("=== POSITIONAL BIAS ANALYSIS ===\n")
  log_cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  log_cat("TSS/TTS boundary:", tss_end, "bp\n\n")
  
  # Check if any motifs are region-specific
  region_summary <- bias_analysis %>%
    group_by(region, expression_type) %>%
    summarise(n_motifs = n(), .groups = "drop")
  
  log_cat("Regional distribution of motifs:\n")
  for (i in 1:nrow(region_summary)) {
    row <- region_summary[i, ]
    log_cat("  ", row$region, " region (", row$expression_type, "): ", row$n_motifs, " motifs\n", sep = "")
  }
  
  # Most abundant motifs by region
  log_cat("\nTop 5 motifs by seqlet count in each region:\n")
  top_motifs <- bias_analysis %>%
    group_by(region) %>%
    top_n(5, total_seqlets) %>%
    arrange(region, desc(total_seqlets))
  
  for (region in unique(top_motifs$region)) {
    log_cat("\n", region, " Region:\n", sep = "")
    region_data <- top_motifs[top_motifs$region == region, ]
    for (j in 1:nrow(region_data)) {
      row <- region_data[j, ]
      log_cat("  ", j, ". ", row$pattern, " (", row$expression_type, 
              ", ", row$total_seqlets, " seqlets, pos: ", round(row$mean_pos), ")\n", sep = "")
    }
  }
  
  # Save detailed results
  bias_file <- file.path(output_dir, "positional_bias_detailed.csv")
  write_csv(bias_analysis, bias_file)
  log_cat("\nDetailed results saved to:", bias_file, "\n")
  
  close(log_con)
  
  return(bias_analysis)
}

# Example usage
if (TRUE) {  # Set to TRUE to run
  # Update this path to your CSV file
  csv_file <- "vitis_cre/src/results/moca_blue/mo_range/seqlet_summary_vitis_ssr_20250801.csv"
  
  # Check if file exists
  if (file.exists(csv_file)) {
    # Create visualization (saves to results/moca_blue/mo_range/)
    result <- visualize_motif_positions(csv_file)
    
    # Analyze positional bias (saves to results/moca_blue/mo_range/)
    bias_analysis <- analyze_positional_bias(result$data)
    
  } else {
    cat("File not found:", csv_file, "\n")
    cat("Please update the csv_file path in the script.\n")
  }
}

# To run the script, set the example usage block to TRUE and update the csv_file path