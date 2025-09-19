#!/usr/bin/env Rscript
# Motif Enrichment Analysis for Gene Expression Associations
# Performs statistical tests and chi-squared analysis for motif-feature associations
# Part 1 of 2-script motif analysis pipeline (feeds into motif_performance.R)
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_proj")
######################

library(tidyr)
library(dplyr)
library(readr)
library(stringr)

# EDIT: Added reticulate for Python pickle file loading
library(reticulate)

# Source utility functions
source("../utils.R")

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Default file paths for Vitis analysis - find files by pattern
# Find most recent q10q90 filtered motif file
motif_pattern <- "../../../out/moca_results/mo_proj/*_vitis_ssr_q10q90_filtered.csv"
motif_files <- Sys.glob(motif_pattern)
default_motif_file <- if (length(motif_files) > 0) sort(motif_files, decreasing = TRUE)[1] else motif_pattern

# Find most recent prediction file
pred_pattern <- "../../../out/predictions/vitis_*_SSR_deepcre_predict_*.csv"
pred_files <- Sys.glob(pred_pattern)
default_prediction_file <- if (length(pred_files) > 0) sort(pred_files, decreasing = TRUE)[1] else pred_pattern
default_expression_file <- "../../../vitis_data/tpm_counts/vitis_drought_leaf_targets.csv"
default_output_dir <- "../../../out/moca_results/mo_proj"
# EDIT: Added validation genes file path for proper validation filtering
default_validation_file <- "../../../vitis_data/validation_genes_vitis.pickle"
default_validation_key <- "vitis"

# Parse arguments
MOTIF_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_motif_file
PREDICTION_FILE <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_prediction_file
EXPRESSION_FILE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_expression_file
OUTPUT_DIR <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_output_dir
# EDIT: Added validation file arguments for proper validation gene filtering
VALIDATION_FILE <- if (length(args) >= 5 && nzchar(args[5])) args[5] else default_validation_file
VALIDATION_KEY <- if (length(args) >= 6 && nzchar(args[6])) args[6] else default_validation_key

# Analysis parameters
DATE_STAMP <- format(Sys.Date(), "%Y%m%d")
PROJECT_NAME <- "vitis_motif_enrichment"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

cat("Motif Enrichment Analysis Parameters:\n")
cat("  Motif file:", MOTIF_FILE, "\n")
cat("  Prediction file:", PREDICTION_FILE, "\n") 
cat("  Expression file:", EXPRESSION_FILE, "\n")
cat("  Validation file:", VALIDATION_FILE, "\n")
cat("  Validation key:", VALIDATION_KEY, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  Project:", PROJECT_NAME, "\n")
cat("  Date:", DATE_STAMP, "\n\n")

#' Load validation genes from pickle file
#' 
#' @param validation_file Path to validation genes pickle file
#' @param validation_key Key for validation genes in pickle dictionary
#' @return Vector of validation gene IDs
load_validation_genes <- function(validation_file, validation_key) {
  
  cat("Loading validation genes...\n")
  
  if (!file.exists(validation_file)) {
    stop("Validation genes file not found: ", validation_file)
  }
  
  # EDIT: Use reticulate to load Python pickle file containing validation genes
  py_run_string(paste0("import pickle"))
  py_run_string(paste0("with open('", validation_file, "', 'rb') as f: val_genes = pickle.load(f)"))
  validation_genes <- py$val_genes[[validation_key]]
  
  if (is.null(validation_genes) || length(validation_genes) == 0) {
    stop("No validation genes found for key '", validation_key, "' in ", validation_file)
  }
  
  cat("  Loaded", length(validation_genes), "validation genes\n")
  
  return(validation_genes)
}

#' Load and validate input datasets for enrichment analysis
#' 
#' This function loads motif occurrences, model predictions, and expression data,
#' ensuring proper formatting and compatibility for downstream analysis.
#' 
#' @param motif_file Path to filtered motif occurrences CSV
#' @param prediction_file Path to model prediction results CSV  
#' @param expression_file Path to expression target data CSV
#' @param validation_genes Vector of validation gene IDs for filtering
#' @return List containing loaded and processed datasets
load_enrichment_data <- function(motif_file, prediction_file, expression_file, validation_genes) {
  
  cat("Loading enrichment analysis datasets...\n")
  
  # Load motif occurrences
  if (!file.exists(motif_file)) {
    stop("Motif file not found: ", motif_file)
  }
  motif_data <- read.csv(motif_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(motif_data), "motif occurrences\n")
  
  # EDIT: Load model predictions and filter for validation genes only to prevent circular validation
  if (!file.exists(prediction_file)) {
    stop("Prediction file not found: ", prediction_file, ". Dummy data removed - use proper validation predictions only.")
  }
  
  prediction_data_raw <- read.csv(prediction_file, stringsAsFactors = FALSE)
  # Standardize column names
  colnames(prediction_data_raw) <- c("true_targets", "pred_probs", "gene_id")
  
  cat("  Loaded", nrow(prediction_data_raw), "raw prediction results\n")
  
  # EDIT: Filter predictions to validation genes only to ensure independent validation
  prediction_data <- prediction_data_raw %>%
    filter(gene_id %in% validation_genes)
  
  pred_available <- nrow(prediction_data) > 0
  
  # EDIT: Validation check - warn if prediction file contains training data
  contamination_ratio <- nrow(prediction_data_raw) / length(validation_genes)
  if (contamination_ratio > 2) {
    cat("  WARNING: Prediction file may contain training data!\n")
    cat("    Raw predictions:", nrow(prediction_data_raw), "genes\n")
    cat("    Validation genes:", length(validation_genes), "genes\n")
    cat("    Ratio:", round(contamination_ratio, 2), "(>2 suggests contamination)\n")
  }
  
  cat("  Filtered to", nrow(prediction_data), "validation gene predictions\n")
  cat("  Prediction file state: VALIDATION_FILTERED\n")
  
  # Load expression data
  if (!file.exists(expression_file)) {
    stop("Expression file not found: ", expression_file)
  }
  expression_data_raw <- read.csv(expression_file, stringsAsFactors = FALSE)
  # Standardize column names
  colnames(expression_data_raw) <- c("gene_id", "target")
  
  cat("  Loaded", nrow(expression_data_raw), "raw expression targets\n")
  
  # EDIT: Filter expression data to validation genes only
  expression_data <- expression_data_raw %>%
    filter(gene_id %in% validation_genes)
  
  cat("  Filtered to", nrow(expression_data), "validation gene expression targets\n")
  
  
  # EDIT: Filter motif data to validation genes only
  motif_data <- motif_data %>%
    filter(gene_id %in% validation_genes)
  
  cat("  Filtered motifs to", nrow(motif_data), "validation gene occurrences\n")
  
  # Data validation and summary
  cat("\nValidation-filtered data summary:\n")
  cat("  Validation genes total:", length(validation_genes), "\n")
  cat("  Unique genes in motifs:", length(unique(motif_data$gene_id)), "\n")
  cat("  Unique genes in predictions:", length(unique(prediction_data$gene_id)), "\n")
  cat("  Unique genes in expression:", length(unique(expression_data$gene_id)), "\n")
  
  # Expression class distribution
  expr_dist <- table(expression_data$target)
  cat("  Expression distribution: Low(0):", expr_dist["0"], 
      " High(1):", expr_dist["1"], " Medium(2):", expr_dist["2"], "\n")
  
  # EDIT: Validation integrity check
  gene_overlap <- intersect(intersect(unique(motif_data$gene_id), unique(prediction_data$gene_id)), unique(expression_data$gene_id))
  cat("  Genes with all data types:", length(gene_overlap), "\n")
  
  if (length(gene_overlap) < 100) {
    warning("Very few genes have all data types - check file compatibility")
  }
  
  return(list(
    motifs = motif_data,
    predictions = prediction_data,
    expression = expression_data,
    pred_available = pred_available,
    validation_genes = validation_genes  # EDIT: Include validation genes for downstream checks
  ))
}

#' Integrate datasets and prepare for enrichment analysis
#' 
#' Merges motif occurrences with predictions and expression data,
#' creates performance classifications, and prepares analysis-ready dataset.
#' 
#' @param data_list List of loaded datasets
#' @return Integrated dataset ready for enrichment analysis
integrate_enrichment_data <- function(data_list) {
  
  cat("\nIntegrating datasets for enrichment analysis...\n")
  
  # Debug: Check input data structure
  cat("  Debug: Motif data columns:", colnames(data_list$motifs), "\n")
  cat("  Debug: Sample motif names:", head(unique(data_list$motifs$motif), 5), "\n")
  if ("epm" %in% colnames(data_list$motifs)) {
    cat("  Debug: EPM column already exists in motif data\n")
    cat("  Debug: Sample existing EPM values:", head(unique(data_list$motifs$epm), 5), "\n")
  }
  
  # Merge motif data with predictions
  integrated_data <- data_list$motifs %>%
    left_join(data_list$predictions, by = "gene_id") %>%
    left_join(data_list$expression, by = "gene_id") %>%
    filter(!is.na(true_targets) | !is.na(target))  # Keep genes with expression data
  
  # Use expression data as ground truth, fall back to prediction targets
  integrated_data$expression_class <- ifelse(
    !is.na(integrated_data$target), 
    integrated_data$target,
    integrated_data$true_targets
  )
  
  # EDIT: Filter out medium expression genes following original V1.7 logic - these were excluded from training
  # Comment from original: "!!! CHANGES DATASET SIZE BY HALF !!! # Were not included within the model training.set"
  cat("  Pre-filter dataset size:", nrow(integrated_data), "associations\n")
  
  integrated_data <- integrated_data %>%
    filter(expression_class %in% c(0, 1))
  
  cat("  Post-filter dataset size (medium genes removed):", nrow(integrated_data), "associations\n")
  
  integrated_data <- integrated_data %>%
    mutate(
      # Create binary classifications
      expr_binary = ifelse(expression_class == 0, "low", "high"),
      pred_binary = ifelse(pred_probs <= 0.5, "low", "high"),
      
      # Extract metacluster information from motif names
      metacluster = ifelse(grepl("_1_", motif), "p1m", "p0m"),
            
      # Use full motif names as EPM for individual motif analysis
      epm = motif
    )
  
  # Debug: Check EPM generation
  cat("  Debug: After EPM generation - unique EPMs:", length(unique(integrated_data$epm)), "\n")
  cat("  Debug: Sample generated EPM names:", head(unique(integrated_data$epm), 5), "\n")
  cat("  Debug: Any EPMs still NA or empty:", sum(is.na(integrated_data$epm) | integrated_data$epm == ""), "\n")
  
  # Check if regex is working on motif names
  sample_motifs <- head(unique(integrated_data$motif), 3)
  cat("  Debug: Sample motif -> EPM conversion:\n")
  for (motif_name in sample_motifs) {
    epm_name <- paste0("epm_vitis_ssr_", gsub(".*_([0-9]+)_([01])_.*", "p\\1m0\\2", motif_name))
    cat("    ", motif_name, " -> ", epm_name, "\n")
  }
  
  integrated_data <- integrated_data %>%
    mutate(
      
      # Calculate prediction performance
      pred_correct = ifelse(expr_binary == pred_binary, "TRUE", "FALSE"),
      
      # EPM-specific performance (following original V1.7 logic)
      epm_pred_perf = case_when(
        pred_correct == "FALSE" ~ "NA",
        grepl("p0m", epm) & expr_binary == "high" ~ "TRUE_high",
        grepl("p1m", epm) & expr_binary == "low" ~ "TRUE_low",
        TRUE ~ "NA"
      ),
      
      # Metacluster-specific performance
      metacluster_correct = case_when(
        metacluster == "p0m" & expr_binary == "high" ~ "TRUE_high",
        metacluster == "p1m" & expr_binary == "low" ~ "TRUE_low", 
        TRUE ~ "FALSE"
      )
    )
  
  cat("  Integrated dataset: ", nrow(integrated_data), " motif-gene associations\n")
  cat("  Expression classes: Low=", sum(integrated_data$expr_binary == "low"),
      " High=", sum(integrated_data$expr_binary == "high"), "\n")
  cat("  Metaclusters: p0m=", sum(integrated_data$metacluster == "p0m"),
      " p1m=", sum(integrated_data$metacluster == "p1m"), "\n")
  
  return(integrated_data)
}

#' Calculate comprehensive enrichment metrics for each motif
#' 
#' Computes TPR, chi-square statistics, importance scores, and other
#' statistical indicators for motif-expression associations.
#' 
#' 
#' @param integrated_data Integrated motif-expression dataset
#' @return Data frame with enrichment metrics for each motif
calculate_motif_enrichment <- function(integrated_data) {
  
  cat("\nCalculating motif enrichment metrics...\n")
  
  # Debug: Check data before contingency table creation
  cat("  Debug: Integrated data rows:", nrow(integrated_data), "\n")
  cat("  Debug: Unique EPMs in integrated data:", length(unique(integrated_data$epm)), "\n")
  cat("  Debug: Unique motifs in integrated data:", length(unique(integrated_data$motif)), "\n")
  cat("  Debug: EPM column has NAs:", sum(is.na(integrated_data$epm)), "\n")
  cat("  Debug: Sample EPM values:", head(unique(integrated_data$epm), 10), "\n")
  
  # Check if EPM generation is working correctly
  cat("  Debug: Sample motif names:", head(unique(integrated_data$motif), 5), "\n")
  
  # Create contingency tables
  contingency_expr <- table(integrated_data$epm, integrated_data$expr_binary)
  contingency_pred <- table(integrated_data$epm, integrated_data$pred_binary)
  contingency_correct <- table(integrated_data$epm, integrated_data$pred_correct)
  
  cat("  Debug: Contingency table dimensions - expr:", dim(contingency_expr), "\n")
  cat("  Debug: Contingency table dimensions - pred:", dim(contingency_pred), "\n")
  cat("  Debug: Contingency table dimensions - correct:", dim(contingency_correct), "\n")
  
  # Convert to data frames for analysis
  expr_df <- as.data.frame(contingency_expr) %>%
    pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0) %>%
    rename(epm = Var1, expr_high = high, expr_low = low)
  
  pred_df <- as.data.frame(contingency_pred) %>%
    pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0) %>%
    rename(epm = Var1, pred_high = high, pred_low = low)
  
  correct_df <- as.data.frame(contingency_correct) %>%
    pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0) %>%
    rename(epm = Var1, pred_false = `FALSE`, pred_true = `TRUE`)
  
  # Merge enrichment tables
  enrichment_metrics <- expr_df %>%
    left_join(pred_df, by = "epm") %>%
    left_join(correct_df, by = "epm") %>%
    mutate(
      # Add metacluster information (convert string to numeric for consistency)
      metacluster = ifelse(grepl("_1_", epm), 1, 0),
      
      # Calculate total occurrences
      total_occurrences = expr_high + expr_low,
      
      # Calculate True Positive Rates
      TPR_expr_low = expr_low / (expr_low + expr_high),
      TPR_expr_high = expr_high / (expr_low + expr_high),
      TPR_pred_low = pred_low / (pred_low + pred_high),
      TPR_pred_high = pred_high / (pred_low + pred_high),
      TPR_correct = pred_true / (pred_true + pred_false),
      
      # Metacluster-specific TPR (expected direction)
      TPR_expected = ifelse(metacluster == 0, TPR_expr_high, TPR_expr_low),
      
      # Calculate totals for importance scores
      total_low_expr = sum(integrated_data$expr_binary == "low"),
      total_high_expr = sum(integrated_data$expr_binary == "high"),
      total_low_pred = sum(integrated_data$pred_binary == "low", na.rm = TRUE),
      total_high_pred = sum(integrated_data$pred_binary == "high", na.rm = TRUE)
    ) %>%
    rowwise() %>%
    mutate(
      # Differential importance scores following original V1.7 logic
      # p0m motifs: log2((high_count/total_high) / (low_count/total_low)) - positive when enriched in high expression
      # p1m motifs: log2((low_count/total_low) / (high_count/total_high)) - positive when enriched in low expression
      imp_expr_score = ifelse(
        metacluster == 1,  # p1m motifs (low expression associated)
        log2(((expr_low + 0.5) / total_low_expr) / ((expr_high + 0.5) / total_high_expr)),
        log2(((expr_high + 0.5) / total_high_expr) / ((expr_low + 0.5) / total_low_expr))  # p0m motifs (high expression associated)
      ),
      imp_pred_score = ifelse(
        metacluster == 1,  # p1m motifs
        log2(((pred_low + 0.5) / total_low_pred) / ((pred_high + 0.5) / total_high_pred)),
        log2(((pred_high + 0.5) / total_high_pred) / ((pred_low + 0.5) / total_low_pred))  # p0m motifs
      ),
      # Combined importance score (primary ranking metric)
      importance_score = imp_expr_score
    ) %>%
    ungroup()
  
  # Calculate chi-square p-values
  enrichment_metrics$chi_expr_pval <- apply(
    enrichment_metrics[, c("expr_low", "expr_high")], 1, function(row) {
      if (sum(row) < 5) return(NA)  # Skip if too few observations
      result <- suppressWarnings(chisq.test(row, p = c(0.5, 0.5)))
      result$p.value
    })
  
  enrichment_metrics$chi_pred_pval <- apply(
    enrichment_metrics[, c("pred_low", "pred_high")], 1, function(row) {
      if (sum(row) < 5) return(NA)
      result <- suppressWarnings(chisq.test(row, p = c(0.5, 0.5)))
      result$p.value
    })
  
  enrichment_metrics$chi_correct_pval <- apply(
    enrichment_metrics[, c("pred_false", "pred_true")], 1, function(row) {
      if (sum(row) < 5) return(NA)
      result <- suppressWarnings(chisq.test(row, p = c(0.5, 0.5)))
      result$p.value
    })
  
  # Rank motifs by enrichment
  enrichment_metrics <- enrichment_metrics %>%
    arrange(desc(importance_score)) %>%
    mutate(
      enrichment_rank = row_number(),
      significance_expr = ifelse(chi_expr_pval < 0.05, "significant", "not_significant"),
      significance_pred = ifelse(chi_pred_pval < 0.05, "significant", "not_significant")
    )
  
  cat("  Calculated metrics for", nrow(enrichment_metrics), "motifs\n")
  cat("  Significantly associated with expression:", 
      sum(enrichment_metrics$significance_expr == "significant", na.rm = TRUE), "\n")
  
  return(enrichment_metrics)
}

#' Generate basic enrichment summary statistics by metacluster
#' 
#' Creates summary statistics following the original V1.7 approach for
#' metacluster-based analysis.
#' 
#' @param enrichment_metrics Calculated enrichment metrics
#' @param integrated_data Full integrated dataset
#' @return Summary statistics by metacluster
generate_enrichment_summary <- function(enrichment_metrics, integrated_data) {
  
  cat("\nGenerating enrichment summary statistics...\n")
  
  # Generate comprehensive summary statistics by metacluster (following original V1.7)
  metacluster_summary <- enrichment_metrics %>%
    group_by(metacluster) %>%
    summarise(
      motif_count = n(),
      total_matches = sum(total_occurrences, na.rm = TRUE),
      mean_importance_expr = mean(imp_expr_score, na.rm = TRUE),
      mean_importance_pred = mean(imp_pred_score, na.rm = TRUE),
      mean_tpr_expr_low = mean(TPR_expr_low, na.rm = TRUE),
      mean_tpr_expr_high = mean(TPR_expr_high, na.rm = TRUE),
      mean_tpr_pred_low = mean(TPR_pred_low, na.rm = TRUE),
      mean_tpr_pred_high = mean(TPR_pred_high, na.rm = TRUE),
      mean_tpr_expected = mean(TPR_expected, na.rm = TRUE),
      mean_tpr_correct = mean(TPR_correct, na.rm = TRUE),
      significant_expr_motifs = sum(significance_expr == "significant", na.rm = TRUE),
      significant_pred_motifs = sum(significance_pred == "significant", na.rm = TRUE),
      mean_pval_expr = mean(chi_expr_pval, na.rm = TRUE),
      mean_pval_pred = mean(chi_pred_pval, na.rm = TRUE),
      mean_pval_correct = mean(chi_correct_pval, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      metacluster_label = ifelse(metacluster == 0, "p0m (High-expr)", "p1m (Low-expr)"),
      # Adjust TPR for expected direction (following corrected logic)
      mean_tpr_expected_adj = ifelse(metacluster == 0, 
                                    mean_tpr_expr_high, 
                                    mean_tpr_expr_low),
      mean_tpr_pred_adj = ifelse(metacluster == 0, 
                                mean_tpr_pred_high, 
                                mean_tpr_pred_low)
    )
  
  cat("  Generated summary for", nrow(metacluster_summary), "metaclusters\n")
  
  return(metacluster_summary)
}

# Main execution
cat("Starting motif enrichment analysis pipeline...\n")

# EDIT: Load validation genes first, then filter all datasets to validation genes only
validation_genes <- load_validation_genes(VALIDATION_FILE, VALIDATION_KEY)

# Load all datasets with validation gene filtering
data_list <- load_enrichment_data(MOTIF_FILE, PREDICTION_FILE, EXPRESSION_FILE, validation_genes)

# Integrate data for analysis
integrated_data <- integrate_enrichment_data(data_list)

# Calculate motif enrichment metrics
enrichment_metrics <- calculate_motif_enrichment(integrated_data)

# Generate summary statistics
enrichment_summary <- generate_enrichment_summary(enrichment_metrics, integrated_data)

# Write output files
output_base <- file.path(OUTPUT_DIR, paste0(DATE_STAMP, "_", PROJECT_NAME))

# Enrichment metrics CSV (main output for motif_performance.R)
write.csv(enrichment_metrics, paste0(output_base, "_results.csv"), row.names = FALSE)
cat("\nOutput files generated:\n")
cat("  Enrichment metrics:", basename(paste0(output_base, "_results.csv")), "\n")

# Integrated data CSV (for further analysis)
write.csv(integrated_data, paste0(output_base, "_integrated_data.csv"), row.names = FALSE)
cat("  Integrated data:", basename(paste0(output_base, "_integrated_data.csv")), "\n")

# Summary statistics CSV
write.csv(enrichment_summary, paste0(output_base, "_metacluster_summary.csv"), row.names = FALSE)
cat("  Metacluster summary:", basename(paste0(output_base, "_metacluster_summary.csv")), "\n")

# Analysis state file (for motif_performance.R)
analysis_state <- list(
  pred_available = data_list$pred_available,
  motif_file = MOTIF_FILE,
  prediction_file = PREDICTION_FILE,
  expression_file = EXPRESSION_FILE,
  date_stamp = DATE_STAMP,
  total_motifs = nrow(enrichment_metrics),
  total_associations = nrow(integrated_data)
)

# Save as RDS for easy loading in motif_performance.R
saveRDS(analysis_state, paste0(output_base, "_analysis_state.rds"))
cat("  Analysis state:", basename(paste0(output_base, "_analysis_state.rds")), "\n")

# Final summary
cat("\n", rep("=", 60), "\n", sep = "")
cat("MOTIF ENRICHMENT ANALYSIS COMPLETED\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Enrichment Analysis Summary:\n")
cat("  Total motifs analyzed:", nrow(enrichment_metrics), "\n")
cat("  Significantly associated with expression:", 
    sum(enrichment_metrics$significance_expr == "significant", na.rm = TRUE), "\n")
cat("  Significantly associated with predictions:", 
    sum(enrichment_metrics$significance_pred == "significant", na.rm = TRUE), "\n")
cat("  Top enriched motif:", enrichment_metrics$epm[1], 
    "(importance:", round(enrichment_metrics$imp_expr_score[1], 3), ")\n")

cat("\nMetacluster Summary:\n")
print(enrichment_summary)

cat("\nNext Step: Run motif_performance.R for visualization and performance analysis\n")
cat("  Use output files generated in:", OUTPUT_DIR, "\n")

cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Enrichment analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")