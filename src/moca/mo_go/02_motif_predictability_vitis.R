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
# OUTPUT STRATEGY:
# The full performance_comparison table (all motif-GO combinations) can be
# extremely large (40GB+) and is NOT written to disk. All biologically relevant
# information is captured in summary outputs:
# - enhanced_categories (top 20 enhanced performance contexts)
# - reduced_categories (top 20 reduced performance contexts)
# - motif_type_summary (aggregated p0m vs p1m statistics)
# - go_category_summary (per-functional-category statistics)
# - top_combinations (top 50 significant motif-GO pairs)
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
# Use outputs from the gene mapper analysis if available (more likely to exist)
perf_pattern <- "../../../out/moca_results/mo_go/*_vitis_gene_motif_mapper_gene_mapping.csv"
perf_files <- Sys.glob(perf_pattern)
default_perf_file <- if (length(perf_files) > 0) sort(perf_files, decreasing = TRUE)[1] else perf_pattern

# Alternatively, rebuild from raw data
motif_pattern <- "../../../out/moca_results/mo_proj/filtering/*_vitis_ssr_q10q90_filtered.csv"
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

# No external helper functions needed - all functionality is self-contained

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
    cat("  Using pre-computed dataset from gene mapper...\n")
    integrated_data <- read.csv(perf_file, stringsAsFactors = FALSE)
    cat("  Loaded", nrow(integrated_data), "records\n")

    # The gene mapping file doesn't have motif data (epm), so we need to add it
    # But we can use this file for genes that have both GO and expression/prediction data

    # Check if we need to add motif data
    if (!"epm" %in% colnames(integrated_data)) {
      cat("  Gene mapping file missing motif data, need to add motifs from motif file...\n")

      # Try to find and load motif data
      if (file.exists(motif_file) && !grepl("\\*", motif_file)) {
        cat("  Loading motif data to merge...\n")

        # Source gene mapper functions to load motif data
        gene_mapper_script <- "01_mo_gene_mapper_vitis.R"
        if (file.exists(gene_mapper_script)) {
          source(gene_mapper_script, local = TRUE)
          motif_data <- load_motif_data(motif_file)

          # Add motif data to the gene mapping (genes can have multiple motifs)
          integrated_data <- integrated_data %>%
            inner_join(motif_data %>% select(gene_id_std, epm),
                      by = "gene_id_std", relationship = "many-to-many")

          cat("  Added motif data, now have", nrow(integrated_data), "records\n")
        } else {
          stop("Cannot find gene mapper script to load motif data")
        }
      } else {
        stop("Cannot find motif file: ", motif_file)
      }
    }

    # Rename columns for compatibility if needed
    if ("primary_category" %in% colnames(integrated_data)) {
      integrated_data$category_name <- integrated_data$primary_category
    }
    if ("primary_description" %in% colnames(integrated_data)) {
      integrated_data$description <- integrated_data$primary_description
    }

    # Add prediction performance columns using existing data
    cat("  Adding prediction performance metrics...\n")
    integrated_data <- integrated_data %>%
      mutate(
        # Overall prediction performance
        pred_perf = case_when(
          (expr_class == "high" & prob_class == "high") |
          (expr_class == "low" & prob_class == "low") ~ "TRUE",
          TRUE ~ "FALSE"
        ),

        # EPM-specific prediction performance
        epm_pred_perf = case_when(
          grepl("p0m", epm) & expr_class == "high" & prob_class == "high" ~ "TRUE_high",
          grepl("p1m", epm) & expr_class == "low" & prob_class == "low" ~ "TRUE_low",
          TRUE ~ "FALSE"
        )
      )

    return(integrated_data)
  }

  # Build from scratch using gene mapper script functions
  cat("  Building integrated dataset from components...\n")

  # Source the gene mapper script which has all the functions we need
  gene_mapper_script <- "01_mo_gene_mapper_vitis.R"
  if (file.exists(gene_mapper_script)) {
    cat("  Sourcing functions from:", gene_mapper_script, "\n")
    source(gene_mapper_script, local = TRUE)
  } else {
    stop("Cannot find gene mapper script: ", gene_mapper_script)
  }

  # Load the raw data using functions from gene mapper
  go_data <- load_go_annotations(go_file)
  motif_data <- load_motif_data(motif_file)
  pred_data <- load_predictions(pred_file)
  expr_data <- load_expression_data(expr_file)

  # Create integrated dataset by joining all data
  cat("  Integrating datasets...\n")

  # Start with motif data (this has the epm column we need)
  integrated_data <- motif_data %>%
    # Add GO annotations
    left_join(go_data %>% select(gene_id_std, category_name, description),
              by = "gene_id_std", relationship = "many-to-many") %>%
    # Add predictions
    left_join(pred_data %>% select(gene_id_std, prob, prob_class),
              by = "gene_id_std") %>%
    # Add expression data
    left_join(expr_data %>% select(gene_id_std, expr_class),
              by = "gene_id_std") %>%
    # Only keep records with all required data
    filter(!is.na(category_name), !is.na(prob), !is.na(expr_class))

  # Add prediction performance columns
  cat("  Adding prediction performance metrics...\n")
  integrated_data <- integrated_data %>%
    mutate(
      # Overall prediction performance
      pred_perf = case_when(
        (expr_class == "high" & prob_class == "high") |
        (expr_class == "low" & prob_class == "low") ~ "TRUE",
        TRUE ~ "FALSE"
      ),

      # EPM-specific prediction performance
      epm_pred_perf = case_when(
        grepl("p0m", epm) & expr_class == "high" & prob_class == "high" ~ "TRUE_high",
        grepl("p1m", epm) & expr_class == "low" & prob_class == "low" ~ "TRUE_low",
        TRUE ~ "FALSE"
      )
    )

  cat("  Integrated dataset created with", nrow(integrated_data), "records\n")
  return(integrated_data)
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
    reframe(
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
      epm_total_correct = epm_correct_high + epm_correct_low,

      # Calculate True Positive Rate for expected class
      tpr_expected = case_when(
        grepl("p0m", epm) ~ epm_correct_high / pmax(high_expr_genes, 1),
        grepl("p1m", epm) ~ epm_correct_low / pmax(low_expr_genes, 1),
        TRUE ~ NA_real_
      ),

      # Overall precision and recall
      precision = correct_predictions / total_genes,

      # Functional category context
      go_description = first(description),

      .groups = "drop"
    ) %>%
    # Filter for manageable analysis size and meaningful results
    filter(total_genes >= 10, total_genes <= 1000) %>%  # Focus on reasonable sample sizes
    # Keep only top combinations by gene count for efficiency
    slice_max(order_by = total_genes, n = 10000, with_ties = FALSE)

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
    reframe(
      baseline_accuracy = sum(pred_perf == "TRUE", na.rm = TRUE) / n(),
      baseline_tpr_expected = case_when(
        grepl("p0m", epm) ~ sum(epm_pred_perf == "TRUE_high", na.rm = TRUE) / pmax(sum(expr_class == "high"), 1),
        grepl("p1m", epm) ~ sum(epm_pred_perf == "TRUE_low", na.rm = TRUE) / pmax(sum(expr_class == "low"), 1),
        TRUE ~ NA_real_
      ),
      total_genes_baseline = n(),
      .groups = "drop"
    )

  # Compare GO-specific performance to baseline (use inner join to control size)
  performance_comparison <- predictability_results %>%
    inner_join(baseline_performance, by = "epm", relationship = "many-to-many") %>%
    mutate(
      # Performance relative to baseline (prevent division by zero)
      accuracy_fold_change = prediction_accuracy / pmax(baseline_accuracy, 0.001),
      tpr_fold_change = tpr_expected / pmax(baseline_tpr_expected, 0.001, na.rm = TRUE),

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

# All functions are self-contained in this script

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

# Performance comparison results (REMOVED - intermediate calculation only)
# The full performance_comparison table contains ~10,000+ rows and can be 40GB+
# All relevant information is captured in the summary outputs below:
# - enhanced_categories (top 20)
# - reduced_categories (top 20)
# - motif_type_summary (aggregated stats)
# - go_category_summary (per-category stats)
# - top_combinations (top 50 significant)
# If full detail is ever needed, can be regenerated by re-running this script
cat("  Performance comparison: Not written (captured in summary outputs)\n")

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