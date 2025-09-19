#!/usr/bin/env Rscript
# Refactored from: motif_predictabilityV1.5.R (original moca_blue project)
#
# This script analyzes motif predictability within GO functional categories,
# focusing on how well different motifs predict expression within specific
# biological functions. Extends the basic performance analysis to include
# functional context-specific predictive power.
#
# ORIGINAL LOGIC: Calculates prediction accuracy metrics for each motif
# within different GO categories, identifies functional categories where
# specific motifs show enhanced or reduced predictive performance.
#
# Adapted for Vitis vinifera project structure with pattern-based file finding
# and outputs to dedicated mo_go results directory.
######################

library(tidyr)
library(dplyr)
library(readr)
library(stringr)

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Find input files by pattern
# Use outputs from the mapping performance analysis if available
perf_pattern <- "../../../out/moca_results/mo_go/*_vitis_go_motif_performance_integrated_data.csv"
perf_files <- Sys.glob(perf_pattern)
default_perf_file <- if (length(perf_files) > 0) sort(perf_files, decreasing = TRUE)[1] else perf_pattern

# Alternatively, rebuild from raw data
motif_pattern <- "../../../out/moca_results/mo_proj/*_vitis_ssr_q10q90_filtered.csv"
motif_files <- Sys.glob(motif_pattern)
default_motif_file <- if (length(motif_files) > 0) sort(motif_files, decreasing = TRUE)[1] else motif_pattern

pred_pattern <- "../../../out/predictions/vitis_*_SSR_deepcre_predict_*.csv"
pred_files <- Sys.glob(pred_pattern)
default_pred_file <- if (length(pred_files) > 0) sort(pred_files, decreasing = TRUE)[1] else pred_pattern

default_go_file <- "../../../vitis_data/gene_ontology/V_vinifera_ont_converted.gmt"
default_expr_file <- "../../../vitis_data/tpm_counts/vitis_drought_leaf_targets.csv"
default_output_dir <- "../../../out/moca_results/mo_go"

# Parse arguments
PERF_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_perf_file
MOTIF_FILE <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_motif_file
PRED_FILE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_pred_file
GO_FILE <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_go_file
EXPR_FILE <- if (length(args) >= 5 && nzchar(args[5])) args[5] else default_expr_file
OUTPUT_DIR <- if (length(args) >= 6 && nzchar(args[6])) args[6] else default_output_dir

# Analysis parameters
DATE_STAMP <- format(Sys.Date(), "%Y%m%d")
PROJECT_NAME <- "vitis_motif_predictability"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

cat("Motif Predictability within GO Categories Analysis:\n")
cat("  Performance file:", PERF_FILE, "\n")
cat("  Motif file:", MOTIF_FILE, "\n")
cat("  Predictions file:", PRED_FILE, "\n")
cat("  GO annotation file:", GO_FILE, "\n")
cat("  Expression file:", EXPR_FILE, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  Date:", DATE_STAMP, "\n\n")

# Source the helper functions from the mapping performance script
source_check_mapping_functions <- function() {
  # Load the functions from our companion script
  mapping_script <- file.path(dirname(rstudioapi::getSourceEditorContext()$path),
                              "mo_check_mapping_performance_vitis.R")

  if (file.exists(mapping_script)) {
    source(mapping_script)
  } else {
    # Fallback: define essential functions locally
    warning("Could not source mapping performance functions, using local definitions")
  }
}

#' Load integrated dataset or build from components
#'
#' @param perf_file Path to pre-computed integrated data
#' @param motif_file Path to motif data
#' @param pred_file Path to predictions
#' @param go_file Path to GO annotations
#' @param expr_file Path to expression data
#' @return Integrated dataset for predictability analysis
load_predictability_data <- function(perf_file, motif_file, pred_file, go_file, expr_file) {

  cat("Loading data for predictability analysis...\n")

  # Try to load pre-computed integrated data first
  if (file.exists(perf_file)) {
    cat("  Using pre-computed integrated dataset...\n")
    integrated_data <- read.csv(perf_file, stringsAsFactors = FALSE)
    cat("  Loaded", nrow(integrated_data), "integrated records\n")
    return(integrated_data)
  }

  # Otherwise, build from scratch (following original logic)
  cat("  Building integrated dataset from components...\n")

  # Load components (requires functions from mapping performance script)
  if (exists("load_go_annotations")) {
    go_data <- load_go_annotations(go_file)
    motif_data <- load_motif_data(motif_file)
    pred_data <- load_predictions(pred_file)
    expr_data <- load_expression_data(expr_file)

    integrated_data <- integrate_datasets(motif_data, pred_data, expr_data, go_data)
    return(integrated_data)
  } else {
    stop("Cannot build integrated dataset - helper functions not available")
  }
}

#' Calculate motif predictability metrics within GO categories
#'
#' Following original motif_predictabilityV1.5.R logic:
#' - For each motif-GO category combination
#' - Calculate prediction accuracy, precision, recall
#' - Identify GO categories where motifs show enhanced performance
#'
#' @param integrated_data Complete dataset with motifs, predictions, and GO
#' @return Predictability metrics by motif and GO category
calculate_go_predictability <- function(integrated_data) {

  cat("Calculating motif predictability within GO categories...\n")

  # Group by motif and GO category
  predictability_results <- integrated_data %>%
    group_by(epm, category_name) %>%
    summarise(
      # Basic counts
      total_genes = n(),

      # Expression class distribution
      high_expr_genes = sum(expr_class == "high"),
      low_expr_genes = sum(expr_class == "low"),

      # Prediction accuracy within this GO category
      correct_predictions = sum(pred_perf == "TRUE"),
      prediction_accuracy = correct_predictions / total_genes,

      # EPM-specific performance (following original logic)
      epm_correct_high = sum(epm_pred_perf == "TRUE_high", na.rm = TRUE),
      epm_correct_low = sum(epm_pred_perf == "TRUE_low", na.rm = TRUE),
      epm_total_correct = epm_correct_high + emp_correct_low,

      # Calculate True Positive Rate for expected class
      tpr_expected = case_when(
        grepl("p0m", epm) ~ epm_correct_high / high_expr_genes,
        grepl("p1m", epm) ~ emp_correct_low / low_expr_genes,
        TRUE ~ NA_real_
      ),

      # Overall precision and recall
      precision = correct_predictions / total_genes,

      # Functional category context
      go_description = first(description),

      .groups = "drop"
    ) %>%
    filter(total_genes >= 5)  # Only analyze categories with sufficient data

  cat("  Analyzed", nrow(predictability_results), "motif-GO combinations\n")
  cat("  Minimum genes per combination:", min(predictability_results$total_genes), "\n")
  cat("  Maximum genes per combination:", max(predictability_results$total_genes), "\n")

  return(predictability_results)
}

#' Identify GO categories with enhanced/reduced motif performance
#'
#' @param predictability_results Results from calculate_go_predictability
#' @param integrated_data Original integrated dataset
#' @return Analysis of GO-specific performance patterns
analyze_go_performance_patterns <- function(predictability_results, integrated_data) {

  cat("Analyzing GO-specific performance patterns...\n")

  # Calculate baseline performance for each motif across all categories
  baseline_performance <- integrated_data %>%
    group_by(epm) %>%
    summarise(
      baseline_accuracy = sum(pred_perf == "TRUE") / n(),
      baseline_tpr_expected = case_when(
        grepl("p0m", epm) ~ sum(epm_pred_perf == "TRUE_high", na.rm = TRUE) / sum(expr_class == "high"),
        grepl("p1m", epm) ~ sum(epm_pred_perf == "TRUE_low", na.rm = TRUE) / sum(expr_class == "low"),
        TRUE ~ NA_real_
      ),
      total_genes_baseline = n(),
      .groups = "drop"
    )

  # Compare GO-specific performance to baseline
  performance_comparison <- predictability_results %>%
    left_join(baseline_performance, by = "epm") %>%
    mutate(
      # Performance relative to baseline
      accuracy_fold_change = prediction_accuracy / baseline_accuracy,
      tpr_fold_change = tpr_expected / baseline_tpr_expected,

      # Classification of performance change
      performance_change = case_when(
        accuracy_fold_change >= 1.2 ~ "enhanced",
        accuracy_fold_change <= 0.8 ~ "reduced",
        TRUE ~ "stable"
      ),

      # Statistical significance (basic threshold)
      is_significant = total_genes >= 10 & abs(log2(accuracy_fold_change)) >= 0.5
    ) %>%
    arrange(desc(accuracy_fold_change))

  # Identify top enhanced and reduced categories
  enhanced_categories <- performance_comparison %>%
    filter(performance_change == "enhanced", is_significant) %>%
    head(20)

  reduced_categories <- performance_comparison %>%
    filter(performance_change == "reduced", is_significant) %>%
    head(20)

  cat("  Enhanced performance combinations:", nrow(enhanced_categories), "\n")
  cat("  Reduced performance combinations:", nrow(reduced_categories), "\n")

  return(list(
    performance_comparison = performance_comparison,
    enhanced_categories = enhanced_categories,
    reduced_categories = reduced_categories,
    baseline_performance = baseline_performance
  ))
}

#' Generate summary statistics by motif type and GO category
#'
#' @param performance_analysis Results from analyze_go_performance_patterns
#' @return Summary tables for reporting
generate_predictability_summaries <- function(performance_analysis) {

  cat("Generating predictability summary statistics...\n")

  # Summary by motif type (p0m vs p1m)
  motif_type_summary <- performance_analysis$performance_comparison %>%
    mutate(motif_type = ifelse(grepl("p0m", epm), "p0m_high_expr", "p1m_low_expr")) %>%
    group_by(motif_type) %>%
    summarise(
      unique_motifs = length(unique(epm)),
      unique_go_categories = length(unique(category_name)),
      total_combinations = n(),
      mean_accuracy = mean(prediction_accuracy, na.rm = TRUE),
      mean_fold_change = mean(accuracy_fold_change, na.rm = TRUE),
      enhanced_combinations = sum(performance_change == "enhanced"),
      reduced_combinations = sum(performance_change == "reduced"),
      significant_combinations = sum(is_significant),
      .groups = "drop"
    )

  # Summary by GO category
  go_category_summary <- performance_analysis$performance_comparison %>%
    group_by(category_name, go_description) %>%
    summarise(
      unique_motifs = length(unique(epm)),
      total_genes = sum(total_genes),
      mean_accuracy = mean(prediction_accuracy, na.rm = TRUE),
      mean_fold_change = mean(accuracy_fold_change, na.rm = TRUE),
      enhanced_motifs = sum(performance_change == "enhanced"),
      reduced_motifs = sum(performance_change == "reduced"),
      .groups = "drop"
    ) %>%
    filter(unique_motifs >= 3) %>%  # Only categories with multiple motifs
    arrange(desc(mean_fold_change))

  # Top performing motif-GO combinations
  top_combinations <- performance_analysis$performance_comparison %>%
    filter(is_significant) %>%
    select(epm, category_name, go_description, total_genes,
           prediction_accuracy, accuracy_fold_change, tpr_expected, performance_change) %>%
    arrange(desc(accuracy_fold_change)) %>%
    head(50)

  return(list(
    motif_type_summary = motif_type_summary,
    go_category_summary = go_category_summary,
    top_combinations = top_combinations
  ))
}

# Main execution
cat("Starting motif predictability within GO categories analysis...\n")

# Load or build integrated dataset
integrated_data <- load_predictability_data(PERF_FILE, MOTIF_FILE, PRED_FILE, GO_FILE, EXPR_FILE)

# Calculate predictability metrics
predictability_results <- calculate_go_predictability(integrated_data)

# Analyze performance patterns
performance_analysis <- analyze_go_performance_patterns(predictability_results, integrated_data)

# Generate summaries
summaries <- generate_predictability_summaries(performance_analysis)

# Write outputs
output_base <- file.path(OUTPUT_DIR, paste0(DATE_STAMP, "_", PROJECT_NAME))

cat("\nWriting output files...\n")

# Main predictability results
write.csv(predictability_results,
          paste0(output_base, "_go_predictability.csv"),
          row.names = FALSE)
cat("  GO predictability:", basename(paste0(output_base, "_go_predictability.csv")), "\n")

# Performance comparison results
write.csv(performance_analysis$performance_comparison,
          paste0(output_base, "_performance_comparison.csv"),
          row.names = FALSE)
cat("  Performance comparison:", basename(paste0(output_base, "_performance_comparison.csv")), "\n")

# Enhanced categories
write.csv(performance_analysis$enhanced_categories,
          paste0(output_base, "_enhanced_categories.csv"),
          row.names = FALSE)
cat("  Enhanced categories:", basename(paste0(output_base, "_enhanced_categories.csv")), "\n")

# Reduced categories
write.csv(performance_analysis$reduced_categories,
          paste0(output_base, "_reduced_categories.csv"),
          row.names = FALSE)
cat("  Reduced categories:", basename(paste0(output_base, "_reduced_categories.csv")), "\n")

# Summary tables
write.csv(summaries$motif_type_summary,
          paste0(output_base, "_motif_type_summary.csv"),
          row.names = FALSE)
cat("  Motif type summary:", basename(paste0(output_base, "_motif_type_summary.csv")), "\n")

write.csv(summaries$go_category_summary,
          paste0(output_base, "_go_category_summary.csv"),
          row.names = FALSE)
cat("  GO category summary:", basename(paste0(output_base, "_go_category_summary.csv")), "\n")

write.csv(summaries$top_combinations,
          paste0(output_base, "_top_combinations.csv"),
          row.names = FALSE)
cat("  Top combinations:", basename(paste0(output_base, "_top_combinations.csv")), "\n")

# Final summary
cat("\n", rep("=", 60), "\n")
cat("MOTIF PREDICTABILITY WITHIN GO CATEGORIES COMPLETED\n")
cat(rep("=", 60), "\n")
cat("Analysis Summary:\n")
cat("  Motif-GO combinations analyzed:", nrow(predictability_results), "\n")
cat("  Enhanced performance combinations:", nrow(performance_analysis$enhanced_categories), "\n")
cat("  Reduced performance combinations:", nrow(performance_analysis$reduced_categories), "\n")
cat("  Unique motifs:", length(unique(predictability_results$epm)), "\n")
cat("  Unique GO categories:", length(unique(predictability_results$category_name)), "\n")

if (nrow(summaries$motif_type_summary) > 0) {
  cat("\nMotif Type Performance:\n")
  print(summaries$motif_type_summary)
}

cat("\nOutput Directory:", OUTPUT_DIR, "\n")
cat("Files Generated: Predictability metrics, performance comparisons, and summary tables\n")

cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")