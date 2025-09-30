#!/usr/bin/env Rscript
# Refactored from: mo_check_mapping-performance_V1.7.R (original moca_blue project)
#
# This script performs GO-based functional enrichment analysis using existing
# motif performance results from the mo_proj pipeline. 
#
# ORIGINAL LOGIC: Creates contingency tables between motifs and GO categories,
# and identifies which motifs are enriched in specific biological functions.
# Performance analysis is handled by existing mo_proj scripts.
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

# Find input files by pattern - use existing performance analysis outputs
# Integrated data from motif enrichment analysis (contains motifs + performance)
integrated_pattern <- "../../../out/moca_results/mo_proj/enrichment/*_vitis_motif_enrichment_integrated_data.csv"
integrated_files <- Sys.glob(integrated_pattern)
default_integrated_file <- if (length(integrated_files) > 0) sort(integrated_files, decreasing = TRUE)[1] else integrated_pattern

# Performance metrics from performance analysis
perf_pattern <- "../../../out/moca_results/mo_proj/performance/*_vitis_motif_performance_motif_performance.csv"
perf_files <- Sys.glob(perf_pattern)
default_perf_file <- if (length(perf_files) > 0) sort(perf_files, decreasing = TRUE)[1] else perf_pattern

# GO annotation file (MapMan format)
default_go_file <- "../../../../vitis_cre/vitis_data/gene_ontology/V_vinifera_ont_converted.gmt"

# Output directory
default_output_dir <- "../../../out/moca_results/mo_go"

# Parse arguments
INTEGRATED_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_integrated_file
PERF_FILE <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_perf_file
GO_FILE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_go_file
OUTPUT_DIR <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_output_dir

# Analysis parameters
DATE_STAMP <- format(Sys.Date(), "%Y%m%d")
PROJECT_NAME <- "vitis_go_motif_performance"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

cat("GO-based Functional Enrichment Analysis Parameters:\n")
cat("  Integrated data file:", INTEGRATED_FILE, "\n")
cat("  Performance metrics file:", PERF_FILE, "\n")
cat("  GO annotation file:", GO_FILE, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  Project:", PROJECT_NAME, "\n")
cat("  Date:", DATE_STAMP, "\n\n")

#' Load and parse GMT format GO annotation file
#'
#' Reads MapMan-style GMT file and creates gene-to-function mapping
#' following the original gene ID standardization logic
#'
#' @param gmt_file Path to GMT format GO file (importantly-reformatted to match Vitis IDs)
#' @return Data frame with gene_id and functional annotations
load_go_annotations <- function(gmt_file) {

  cat("Loading GO annotations from GMT file...\n")

  if (!file.exists(gmt_file)) {
    stop("GO annotation file not found: ", gmt_file)
  }

  # Read GMT file
  gmt_lines <- readLines(gmt_file)

  # Parse GMT format: category_name \t description \t gene1 \t gene2 \t ...
  go_data <- data.frame()

  for (line in gmt_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      category_name <- parts[1]
      description <- parts[2]
      genes <- parts[3:length(parts)]

      # Create rows for each gene-category pair
      for (gene in genes) {
        if (nzchar(gene)) {  # Skip empty genes
          go_data <- rbind(go_data, data.frame(
            gene_id = gene,
            category_name = category_name,
            description = description,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  cat("  Loaded", nrow(go_data), "gene-category associations\n")
  cat("  Unique genes:", length(unique(go_data$gene_id)), "\n")
  cat("  Unique categories:", length(unique(go_data$category_name)), "\n")

  return(go_data)
}


#' Load and prepare motif occurrence data
#'
#' @param motif_file Path to filtered motif occurrence data
#' @return Data frame with standardized motif data
load_motif_data <- function(motif_file) {

  cat("Loading motif occurrence data...\n")

  if (!file.exists(motif_file)) {
    stop("Motif file not found: ", motif_file)
  }

  motif_data <- read.csv(motif_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(motif_data), "motif occurrences\n")

  # Standardize gene IDs
  motif_data$gene_id_std <- motif_data$gene_id  # No standardization needed

  # Extract EPM format (following original naming)
  # Assuming epm column exists, otherwise derive from motif column
  if (!"epm" %in% colnames(motif_data)) {
    motif_data$epm <- motif_data$motif  # Use motif column as fallback
  }

  cat("  Unique genes with motifs:", length(unique(motif_data$gene_id_std)), "\n")
  cat("  Unique motifs:", length(unique(motif_data$epm)), "\n")

  return(motif_data)
}

#' Load model predictions
#'
#' @param pred_file Path to model prediction file
#' @return Data frame with predictions and expression classes
load_predictions <- function(pred_file) {

  cat("Loading model predictions...\n")

  if (!file.exists(pred_file)) {
    stop("Prediction file not found: ", pred_file)
  }

  pred_data <- read.csv(pred_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(pred_data), "predictions\n")

  # Standardize gene IDs
  pred_data$gene_id_std <- pred_data$gene_id  # No standardization needed

  # Create probability classes (following original logic)
  pred_data$prob_class <- ifelse(pred_data$prob > 0.5, "high", "low")

  cat("  High probability predictions:", sum(pred_data$prob_class == "high"), "\n")
  cat("  Low probability predictions:", sum(pred_data$prob_class == "low"), "\n")

  return(pred_data)
}

#' Load expression data and create classes
#'
#' @param expr_file Path to expression data file
#' @return Data frame with expression classes
load_expression_data <- function(expr_file) {

  cat("Loading expression data...\n")

  if (!file.exists(expr_file)) {
    stop("Expression file not found: ", expr_file)
  }

  expr_data <- read.csv(expr_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(expr_data), "expression measurements\n")

  # Standardize gene IDs
  expr_data$gene_id_std <- expr_data$gene_id  # No standardization needed

  # Create expression classes (assuming target_class column exists)
  if ("target_class" %in% colnames(expr_data)) {
    expr_data$expr_class <- ifelse(expr_data$target_class == 1, "high", "low")
  } else {
    # Fallback: create classes based on median split of expression values
    if ("tpm" %in% colnames(expr_data)) {
      median_expr <- median(expr_data$tpm, na.rm = TRUE)
      expr_data$expr_class <- ifelse(expr_data$tpm >= median_expr, "high", "low")
    } else {
      stop("Cannot determine expression classes - no target_class or tpm column found")
    }
  }

  cat("  High expression genes:", sum(expr_data$expr_class == "high"), "\n")
  cat("  Low expression genes:", sum(expr_data$expr_class == "low"), "\n")

  return(expr_data)
}

#' Integrate all datasets following original merge logic
#'
#' @param motif_data Motif occurrence data
#' @param pred_data Model predictions
#' @param expr_data Expression data
#' @param go_data GO annotations
#' @return Integrated dataset
integrate_datasets <- function(motif_data, pred_data, expr_data, go_data) {

  cat("Integrating datasets...\n")

  # Standardize GO gene IDs
  go_data$gene_id_std <- go_data$gene_id  # No standardization needed

  # Step 1: Merge motif data with predictions (following original logic)
  merged_step1 <- merge(motif_data, pred_data,
                       by = "gene_id_std",
                       all.x = TRUE, all.y = TRUE)
  merged_step1 <- na.omit(merged_step1)
  cat("  After motif-prediction merge:", nrow(merged_step1), "rows\n")

  # Step 2: Merge with expression data
  merged_step2 <- merge(merged_step1, expr_data,
                       by = "gene_id_std",
                       all.x = TRUE, all.y = TRUE)
  merged_step2 <- na.omit(merged_step2)
  cat("  After expression merge:", nrow(merged_step2), "rows\n")

  # Step 3: Merge with GO annotations
  integrated_data <- merge(merged_step2, go_data,
                          by = "gene_id_std",
                          all = FALSE, ignore.case = TRUE)
  cat("  After GO annotation merge:", nrow(integrated_data), "rows\n")

  # Performance metrics are already calculated in the pre-processed data from mo_proj
  # No need to recalculate pred_perf and epm_pred_perf as they already exist

  cat("  Final integrated dataset:", nrow(integrated_data), "rows\n")
  cat("  Unique genes:", length(unique(integrated_data$gene_id_std)), "\n")
  cat("  Unique motifs:", length(unique(integrated_data$epm)), "\n")
  cat("  Unique GO categories:", length(unique(integrated_data$category_name)), "\n")

  return(integrated_data)
}

#' Create contingency tables and calculate performance metrics
#'
#' Following the original contingency table analysis approach
#'
#' @param integrated_data Complete integrated dataset
#' @return List of analysis results
analyze_motif_go_performance <- function(integrated_data) {

  cat("Analyzing motif-GO performance relationships...\n")

  # Memory check before analysis
  cat("  Input data size:", nrow(integrated_data), "rows,", ncol(integrated_data), "columns\n")

# Note: Data has been pre-filtered for manageable analysis size

  # Create contingency tables (following original script exactly)

  # 1. Motif vs Expression Class
  contingency_expr_class <- table(integrated_data$epm, integrated_data$expr_binary)

  # 2. Motif vs Probability Class
  contingency_prob_class <- table(integrated_data$epm, integrated_data$pred_binary)

  # 3. Motif vs EPM Prediction Performance (key analysis)
  epm_perf_data <- subset(integrated_data, epm_pred_perf != "NA")
  contingency_epm_pred_perf <- table(epm_perf_data$epm, epm_perf_data$epm_pred_perf)

  # 4. Motif vs GO Categories (the main GO analysis)
  contingency_motif_go <- table(integrated_data$epm, integrated_data$category_name)

  cat("  Created contingency tables:\n")
  cat("    Motif vs Expression Class:", nrow(contingency_expr_class), "motifs x", ncol(contingency_expr_class), "classes\n")
  cat("    Motif vs GO Categories:", nrow(contingency_motif_go), "motifs x", ncol(contingency_motif_go), "categories\n")

  # Convert contingency tables to wide format (following original pivot logic)

  # Expression class table
  expr_table <- as.data.frame(contingency_expr_class) %>%
    pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0) %>%
    rename(epm = Var1, expr_class_high = high, expr_class_low = low)

  # Probability class table
  prob_table <- as.data.frame(contingency_prob_class) %>%
    pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0) %>%
    rename(epm = Var1, prob_class_high = high, prob_class_low = low)

  # EPM prediction performance table
  epm_perf_table <- as.data.frame(contingency_epm_pred_perf) %>%
    pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0) %>%
    rename(epm = Var1)

  # Merge performance tables
  performance_summary <- merge(expr_table,
                              merge(prob_table, epm_perf_table, by = "epm", all = TRUE),
                              by = "epm", all = TRUE)

  # Calculate True Positive metrics (following original logic)
  if ("TRUE_high" %in% colnames(performance_summary) && "TRUE_low" %in% colnames(performance_summary)) {
    performance_summary$epm_TPpred <- ifelse(
      is.na(performance_summary$TRUE_high) | performance_summary$TRUE_high == 0,
      performance_summary$TRUE_low,
      performance_summary$TRUE_high
    )
    performance_summary$epm_TNpred <- (performance_summary$expr_class_high + performance_summary$expr_class_low) -
                                      performance_summary$epm_TPpred
  }

  # GO category enrichment analysis
  go_enrichment <- as.data.frame(contingency_motif_go) %>%
    rename(epm = Var1, go_category = Var2, count = Freq) %>%
    filter(count > 0)  # Only keep non-zero associations

  return(list(
    performance_summary = performance_summary,
    go_enrichment = go_enrichment,
    contingency_tables = list(
      expr_class = contingency_expr_class,
      prob_class = contingency_prob_class,
      emp_perf = contingency_epm_pred_perf,
      motif_go = contingency_motif_go
    ),
    integrated_data = integrated_data
  ))
}

#' Load existing performance analysis results
#'
#' @param integrated_file Path to integrated data from motif enrichment analysis
#' @param perf_file Path to performance metrics from performance analysis
#' @return List containing integrated data and performance metrics
load_existing_analysis_results <- function(integrated_file, perf_file) {

  cat("Loading existing performance analysis results...\n")

  # Load integrated data
  if (!file.exists(integrated_file)) {
    stop("Integrated data file not found: ", integrated_file)
  }
  integrated_data <- read.csv(integrated_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(integrated_data), "integrated motif-gene associations\n")

  # Load performance metrics (optional)
  performance_metrics <- NULL
  if (file.exists(perf_file)) {
    performance_metrics <- read.csv(perf_file, stringsAsFactors = FALSE)
    cat("  Loaded", nrow(performance_metrics), "motif performance metrics\n")
  } else {
    cat("  Performance metrics file not found, proceeding without detailed metrics\n")
  }

  return(list(
    integrated_data = integrated_data,
    performance_metrics = performance_metrics
  ))
}

# Main execution
cat("Starting GO-based functional enrichment analysis pipeline...\n")

# Memory optimization: Enable garbage collection and set options
gc()
options(warn = 1)  # Show warnings immediately

# Load existing analysis results
existing_results <- load_existing_analysis_results(INTEGRATED_FILE, PERF_FILE)
integrated_data <- existing_results$integrated_data
performance_metrics <- existing_results$performance_metrics

# Memory optimization: Garbage collection without sampling (full dataset analysis)
gc()

# Load GO annotations
go_annotations <- load_go_annotations(GO_FILE)

# Integrate with GO annotations (add functional categories to existing data)
cat("Integrating existing results with GO annotations...\n")
go_annotations$gene_id_std <- go_annotations$gene_id  # No standardization needed

# Merge integrated data with GO annotations
if ("gene_id" %in% colnames(integrated_data)) {
  integrated_data$gene_id_std <- integrated_data$gene_id  # No standardization needed
} else if ("gene_id_std" %in% colnames(integrated_data)) {
  # Already standardized
} else {
  stop("Cannot find gene ID column in integrated data")
}

# Memory-safe merge with size checking
cat("  Attempting GO integration...\n")
cat("  Integrated data size:", nrow(integrated_data), "rows\n")
cat("  GO annotations size:", nrow(go_annotations), "rows\n")

# Check for reasonable merge size before proceeding
common_genes <- intersect(integrated_data$gene_id_std, go_annotations$gene_id_std)
cat("  Common genes for merge:", length(common_genes), "\n")

if (length(common_genes) == 0) {
  stop("ERROR: No common genes found between integrated data and GO annotations")
}

# Perform merge with explicit relationship specification
go_integrated_data <- merge(integrated_data, go_annotations,
                           by = "gene_id_std", all.x = FALSE,
                           suffixes = c("", "_go"))

# Force garbage collection after merge
gc()

cat("  After GO integration:", nrow(go_integrated_data), "motif-gene-GO associations\n")

# Apply confidence and size-based filtering before analysis (following original approach)
cat("Applying confidence-based filtering to reduce dataset size...\n")

# 1. Probability confidence filtering (from original mo_check_mapping-performance script)
# Original tested prob > 0.8 and prob > 0.01 - use moderate confidence threshold
initial_size <- nrow(go_integrated_data)
go_integrated_data <- go_integrated_data %>%
  filter(abs(pred_probs - 0.5) >= 0.1)  # Only confident predictions (pred_probs <= 0.4 or pred_probs >= 0.6)

cat("  After probability confidence filtering:", nrow(go_integrated_data), "associations (removed",
    initial_size - nrow(go_integrated_data), "uncertain predictions)\n")

# 2. Motif frequency filtering (only analyze motifs with sufficient representation)
motif_frequencies <- go_integrated_data %>%
  group_by(epm) %>%
  summarise(n_genes = n_distinct(gene_id_std), .groups = "drop") %>%
  filter(n_genes >= 10)  # Only motifs appearing in ≥10 genes

go_integrated_data <- go_integrated_data %>%
  filter(epm %in% motif_frequencies$epm)

cat("  After motif frequency filtering:", nrow(go_integrated_data), "associations\n")
cat("  Retained motifs:", nrow(motif_frequencies), "out of",
    length(unique(integrated_data$epm)), "\n")

# Final dataset summary
cat("\nFiltered dataset summary:\n")
cat("  Total associations:", nrow(go_integrated_data), "\n")
cat("  Unique genes:", length(unique(go_integrated_data$gene_id_std)), "\n")
cat("  Unique motifs:", length(unique(go_integrated_data$epm)), "\n")
cat("  Unique GO categories:", length(unique(go_integrated_data$category_name)), "(all retained)\n")
cat("  Data reduction:", sprintf("%.1f%%", 100 * (1 - nrow(go_integrated_data)/initial_size)), "\n")

# Force garbage collection after filtering
gc()

# Perform GO-specific enrichment analysis on filtered data
analysis_results <- analyze_motif_go_performance(go_integrated_data)

# Generate outputs
output_base <- file.path(OUTPUT_DIR, paste0(DATE_STAMP, "_", PROJECT_NAME))

cat("\nWriting output files...\n")

# Main performance summary (equivalent to original CSV output)
write.csv(analysis_results$performance_summary,
          paste0(output_base, "_performance_summary.csv"),
          row.names = FALSE)
cat("  Performance summary:", basename(paste0(output_base, "_performance_summary.csv")), "\n")

# GO enrichment results
write.csv(analysis_results$go_enrichment,
          paste0(output_base, "_go_enrichment.csv"),
          row.names = FALSE)
cat("  GO enrichment:", basename(paste0(output_base, "_go_enrichment.csv")), "\n")

# Complete integrated dataset for further analysis
write.csv(integrated_data,
          paste0(output_base, "_integrated_data.csv"),
          row.names = FALSE)
cat("  Integrated dataset:", basename(paste0(output_base, "_integrated_data.csv")), "\n")

# Summary statistics
cat("\n", rep("=", 60), "\n")
cat("GO-BASED MOTIF PERFORMANCE ANALYSIS COMPLETED\n")
cat(rep("=", 60), "\n")
cat("Analysis Summary:\n")
cat("  Total motif-gene-GO associations:", nrow(integrated_data), "\n")
cat("  Unique motifs analyzed:", length(unique(integrated_data$epm)), "\n")
cat("  Unique GO categories:", length(unique(integrated_data$category_name)), "\n")
cat("  Motifs with performance data:", nrow(analysis_results$performance_summary), "\n")
cat("  Non-zero motif-GO associations:", nrow(analysis_results$go_enrichment), "\n")

cat("\nOutput Directory:", OUTPUT_DIR, "\n")
cat("Files Generated: Performance summary, GO enrichment, and integrated dataset\n")

cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")