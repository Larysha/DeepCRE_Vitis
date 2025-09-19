#!/usr/bin/env Rscript
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre")
# R script to analyse mo_nom outputs with summary statistics and exploratory visuals
# Run from: vitis_cre/src/moca/mo_nom/
# Outputs go to: ../../out/moca_results/mo_nom/

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(gridExtra)
library(RColorBrewer)
library(cowplot)

# Custom colour palette
custom_colors <- c("#AF91AF", "#804a5f", "#478fca", "#3e2c57")

# Set output directory relative to current working directory
# Check if we're in vitis_cre or mo_nom directory and adjust path accordingly
if (basename(getwd()) == "vitis_cre") {
  output_dir <- "out/moca_results/mo_nom"
} else {
  output_dir <- "../../../out/moca_results/mo_nom"
}

# Function to try loading specific file types
try_load_specific_data <- function(file_type) {
  # Find all CSV files matching the pattern
  all_csv_files <- list.files(output_dir, pattern = ".*vitisssr_.*_summary\\.csv$", full.names = TRUE)

  if (file_type == "cwm_unscaled") {
    target_files <- grep("cwm_unscaled_summary", all_csv_files, value = TRUE)
  } else if (file_type == "cwm_scaled") {
    target_files <- grep("cwm_summary", all_csv_files, value = TRUE)
    target_files <- target_files[!grepl("unscaled", target_files)]  # Exclude unscaled
  } else if (file_type == "pfm") {
    target_files <- grep("pfm_summary", all_csv_files, value = TRUE)
  } else {
    target_files <- c()
  }

  contrib_df <- NULL
  loaded_file <- NULL

  for (file_path in target_files) {
    tryCatch({
      if (file.exists(file_path)) {
        contrib_df <- read_csv(file_path, show_col_types = FALSE)
        loaded_file <- file_path
        cat("Successfully loaded:", basename(file_path), "\n")
        break
      }
    }, error = function(e) {
      cat("Error loading", file_path, ":", e$message, "\n")
    })
  }

  return(list(data = contrib_df, file = loaded_file))
}

# Function to parse motif names
parse_motif_name <- function(motif) {
  parts <- strsplit(as.character(motif), "_")[[1]]

  if (length(parts) >= 6 && parts[1] == "epm") {
    tryCatch({
      pattern <- paste0("pattern_", parts[4])  # pattern number
      metacluster <- parts[5]  # metacluster (0 or 1)
      strand <- ifelse(parts[6] == "F", "fwd", "rev")
      n_seqlets <- if (length(parts) > 6) as.numeric(parts[7]) else NA

      return(data.frame(
        pattern = pattern,
        strand = strand,
        metacluster = metacluster,
        n_seqlets = n_seqlets,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      return(data.frame(
        pattern = NA,
        strand = NA,
        metacluster = NA,
        n_seqlets = NA,
        stringsAsFactors = FALSE
      ))
    })
  } else {
    return(data.frame(
      pattern = NA,
      strand = NA,
      metacluster = NA,
      n_seqlets = NA,
      stringsAsFactors = FALSE
    ))
  }
}

# Function to analyze a single dataset
analyze_single_dataset <- function(contrib_df, loaded_file, analysis_type) {
  if (is.null(contrib_df)) {
    return(FALSE)
  }

  cat("### Motif Contribution Scores -", analysis_type, "(from", basename(loaded_file), ")\n\n")

  # Display basic info
  cat("Available columns:", paste(colnames(contrib_df), collapse = ", "), "\n")
  cat("Data shape:", nrow(contrib_df), "rows x", ncol(contrib_df), "columns\n\n")

  # Show first 10 rows
  cat("#### First 10 rows:\n")
  print(head(contrib_df, 10))
  cat("\n")

  # Parse motif names
  parsed_data <- do.call(rbind, lapply(contrib_df$motif, parse_motif_name))
  contrib_df <- cbind(contrib_df, parsed_data)

  # Remove rows where parsing failed
  contrib_df <- contrib_df[!is.na(contrib_df$pattern) & !is.na(contrib_df$metacluster), ]
  cat("After parsing and filtering:", nrow(contrib_df), "rows\n\n")

  # Convert metacluster to factor for plotting
  contrib_df$metacluster <- as.factor(contrib_df$metacluster)

  cat("### Summary Statistics -", analysis_type, "\n")

  # Identify numeric columns for plotting (only score_sum and seqlet_count)
  plot_cols <- c("score_sum", "seqlet_count")
  available_plot_cols <- plot_cols[plot_cols %in% colnames(contrib_df)]

  # All numeric columns for summary statistics
  numeric_cols <- c("score_sum", "score_max", "score_range", "seqlet_count")
  available_numeric_cols <- numeric_cols[numeric_cols %in% colnames(contrib_df)]

  cat("Using numeric columns:", paste(available_numeric_cols, collapse = ", "), "\n\n")

  if (length(available_numeric_cols) > 0) {
    # Summary statistics
    summary_stats <- contrib_df[available_numeric_cols] %>%
      summarise_all(list(
        mean = ~round(mean(., na.rm = TRUE), 3),
        median = ~round(median(., na.rm = TRUE), 3),
        sd = ~round(sd(., na.rm = TRUE), 3),
        min = ~round(min(., na.rm = TRUE), 3),
        max = ~round(max(., na.rm = TRUE), 3)
      ))

    # Reshape for better display
    summary_long <- summary_stats %>%
      pivot_longer(everything(), names_to = "metric_stat", values_to = "value") %>%
      separate(metric_stat, into = c("metric", "statistic"), sep = "_(?=[^_]*$)") %>%
      pivot_wider(names_from = statistic, values_from = value)

    print(summary_long)
    cat("\n")

    # Create visualizations
    n_metrics <- length(available_plot_cols)
    if (n_metrics > 0) {
      # Determine plot layout
      if (n_metrics <= 2) {
        n_cols <- n_metrics
        n_rows <- 1
      } else if (n_metrics <= 4) {
        n_cols <- 2
        n_rows <- 2
      } else {
        n_cols <- 3
        n_rows <- 2
      }

      # Create plots
      plot_list <- list()

      for (i in seq_along(available_plot_cols)) {
        score_col <- available_plot_cols[i]

        # Aggregate by pattern and metacluster
        pattern_summary <- contrib_df %>%
          group_by(pattern, metacluster) %>%
          summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop")

        p <- ggplot(pattern_summary, aes(x = pattern, y = mean_score, fill = metacluster)) +
          geom_bar(stat = "identity", position = "dodge") +
          scale_fill_manual(values = custom_colors[1:length(levels(contrib_df$metacluster))],
                          labels = paste0("MC", levels(contrib_df$metacluster),
                                        c(" (High Expression)", " (Low Expression)")[1:length(levels(contrib_df$metacluster))])) +
          labs(title = paste("Mean", gsub("_", " ", stringr::str_to_title(score_col)), "by Motif Pattern"),
               y = paste("Mean", gsub("_", " ", stringr::str_to_title(score_col))),
               x = "Motif Pattern",
               fill = "Metacluster") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
                axis.title.x = element_text(size = rel(0.8)),
                axis.title.y = element_text(size = rel(0.8)),
                plot.title = element_text(size = rel(0.8)),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.position = "none")

        plot_list[[i]] <- p
      }

      # Save plots with analysis type in filename
      analysis_suffix <- tolower(gsub(" ", "_", analysis_type))
      output_plot_file <- file.path(output_dir, paste0("motif_summary_plots_", analysis_suffix, ".png"))
      png(output_plot_file, width = 1200 * n_cols / 2, height = 900 * n_rows / 2, res = 150)

      # Create one plot with legend to extract it
      temp_plot <- plot_list[[1]] +
        theme(legend.position = "bottom",
              legend.justification = "left",
              legend.title = element_text(size = rel(0.9)),
              legend.text = element_text(size = rel(0.9))) +
        guides(fill = guide_legend(title = "Metacluster", ncol = 2))

      # Extract legend
      legend <- get_legend(temp_plot)

      # Create final combined plot
      if (length(plot_list) == 1) {
        final_plot <- plot_grid(plot_list[[1]], legend, ncol = 1, rel_heights = c(0.85, 0.15))
      } else {
        plots_combined <- plot_grid(plotlist = plot_list, ncol = n_cols)
        final_plot <- plot_grid(plots_combined, legend, ncol = 1, rel_heights = c(0.85, 0.15))
      }

      print(final_plot)
      dev.off()
      cat("Plots saved to:", output_plot_file, "\n\n")

      # Top motifs analysis
      if ("score_sum" %in% available_numeric_cols) {
        cat("### Top 10 Motifs by Score Sum -", analysis_type, "\n")
        top_motifs_cols <- c("motif", "pattern", "metacluster", "score_sum")
        if ("n_seqlets" %in% colnames(contrib_df)) {
          top_motifs_cols <- c(top_motifs_cols, "n_seqlets")
        }

        top_motifs <- contrib_df %>%
          arrange(desc(score_sum)) %>%
          slice_head(n = 10) %>%
          select(all_of(top_motifs_cols))

        print(top_motifs)
        cat("\n")
      }

      # Metacluster comparison
      cat("### Metacluster Comparison -", analysis_type, "\n")
      metacluster_stats <- contrib_df %>%
        group_by(metacluster) %>%
        summarise(across(all_of(available_numeric_cols),
                        list(mean = ~round(mean(., na.rm = TRUE), 3),
                             sd = ~round(sd(., na.rm = TRUE), 3),
                             count = ~sum(!is.na(.))),
                        .names = "{.col}_{.fn}"),
                 .groups = "drop")

      print(metacluster_stats)
      cat("\n")

      # Save summary data
      analysis_suffix <- tolower(gsub(" ", "_", analysis_type))
      output_summary_file <- file.path(output_dir, paste0("motif_analysis_summary_", analysis_suffix, ".csv"))
      write_csv(contrib_df, output_summary_file)
      cat("Enhanced data with parsed motif information saved to:", output_summary_file, "\n\n")
    }

  } else {
    cat("No numeric columns found for analysis.\n")
  }

  return(TRUE)
}

# Main analysis function
analyze_motif_data <- function() {
  cat("=== Motif Contribution Analysis ===\n\n")

  # Try to load both scaled and unscaled CWM data
  unscaled_result <- try_load_specific_data("cwm_unscaled")
  scaled_result <- try_load_specific_data("cwm_scaled")

  analyses_completed <- 0

  # Analyze unscaled data if available
  if (!is.null(unscaled_result$data)) {
    analyze_single_dataset(unscaled_result$data, unscaled_result$file, "CWM Unscaled")
    analyses_completed <- analyses_completed + 1
  }

  # Analyze scaled data if available
  if (!is.null(scaled_result$data)) {
    analyze_single_dataset(scaled_result$data, scaled_result$file, "CWM Scaled")
    analyses_completed <- analyses_completed + 1
  }

  # If no data was found, show error message
  if (analyses_completed == 0) {
    cat("### Motif Contribution Analysis - Data Not Available\n\n")
    cat("The contribution score analysis requires CSV files generated by get_matrices.R\n")
    cat("Expected files in", output_dir, ":\n")
    cat("- 20250914_vitisssr_cwm_summary.csv\n")
    cat("- 20250914_vitisssr_pfm_summary.csv\n")
    cat("- 20250811_vitisssr_cwm_summary.csv\n")
    cat("- 20250811_vitisssr_pfm_summary.csv\n")
    cat("- 20250811_vitisssr_pwm_summary.csv\n\n")

    # Check for available files
    available_files <- list.files(output_dir, pattern = "*.csv", full.names = TRUE)
    if (length(available_files) > 0) {
      cat("Available CSV files in output directory:\n")
      for (f in available_files) {
        cat("-", basename(f), "\n")
      }
    }
    return(invisible(NULL))
  }

  # Add summary of what was analyzed
  if (analyses_completed > 0) {
    cat("### Analysis Summary\n")
    cat("Completed", analyses_completed, "analysis/analyses\n")
    if (analyses_completed == 2) {
      cat("Compare the scaled vs unscaled results to understand:\n")
      cat("- Scaled: Better for pattern discovery and clustering\n")
      cat("- Unscaled: Shows true biological importance and impact\n")
    }
    cat("\n")
  }

  cat("### Key Observations:\n")
  cat("- Score Sum: Total contribution magnitude across all motif positions\n")
  cat("- Score Max/Min: Peak contribution values (positive and negative)\n")
  cat("- Score Range: Spread of contribution values within each motif\n")
  cat("- Metacluster patterns: MC0 (high expression) vs MC1 (low expression) motifs\n")
  cat("- Supporting evidence: Number of seqlets indicates motif robustness\n")
  cat("- Pattern distribution: Shows how motifs are distributed across discovered patterns\n\n")
}

# Check if script is being run directly or sourced
if (!interactive()) {
  # Set working directory to script location if run from command line
  script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[4]))
  if (!is.na(script_dir) && script_dir != ".") {
    setwd(script_dir)
  }

  # Run analysis
  analyze_motif_data()

  cat("Analysis complete!\n")
  cat("Run from:", getwd(), "\n")
  cat("Outputs in:", normalizePath(file.path(getwd(), output_dir)), "\n")
} else {
  # If sourced interactively, just run the analysis
  analyze_motif_data()
}
