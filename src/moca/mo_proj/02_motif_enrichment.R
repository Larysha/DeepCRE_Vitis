#!/usr/bin/env Rscript
# Motif Enrichment Analysis for Gene Expression Associations
# Performs statistical tests and chi-squared analysis for motif-feature associations
# Part 1 of 2-script motif analysis pipeline (feeds into motif_performance.R)
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_proj")
######################
#
# SCRIPT OVERVIEW: This script is the foundation of biologically-informed motif performance analysis
#
# CORE LOGIC:
# 1. Load motifs, model predictions, and expression data (ALL GENES - original approach)
# 2. Integrate datasets to create motif-gene-expression-prediction associations
# 3. Calculate enrichment scores using log2 fold-change relative to background
# 4. Rank motifs by biological relevance (enrichment), not just technical performance
# 5. Provide statistical validation through chi-square tests
#
# KEY BIOLOGICAL HYPOTHESIS:
# - p0m motifs should be enriched in HIGH expression genes
# - p1m motifs should be enriched in LOW expression genes
# - Importance scores measure how well motifs follow these expectations
#
# OUTPUT: Ranked motifs by biological significance + integrated dataset for performance analysis
######################

library(tidyr)
library(dplyr)
library(readr)
library(stringr)

# Source utility functions
source("../utils.R")

# Convert BLAMM format motif names to expected p0m/p1m format
convert_epm_vitis <- function(code) {
  # Convert BLAMM format to seqlet statistics format
  # BLAMM: "epm_vitis_ssr_5_0_R_534" (pattern_metacluster_strand_count)
  # Target: "epm_vitis_ssr_p0m05" (p{metacluster}m{pattern_padded})

  # Split by underscore to get components
  parts <- strsplit(code, "_")

  converted <- sapply(parts, function(p) {
    if (length(p) >= 6) {
      # Extract pattern number (4th element) and metacluster (5th element)
      pattern_num <- as.numeric(p[4])
      metacluster <- as.numeric(p[5])

      # Format as p{metacluster}m{pattern_padded} to match seqlet statistics
      formatted_pattern <- sprintf("p%dm%02d", metacluster, pattern_num)

      # Reconstruct: epm_vitis_ssr_p{metacluster}m{pattern}
      paste("epm_vitis_ssr", formatted_pattern, sep = "_")
    } else {
      # Fallback for unexpected formats
      code
    }
  })

  return(converted)
}

# Enhanced dual naming system using existing clustering data

#' Load motif clustering data to get long format names with sequences
#' @param clustering_file Path to clustering motif summary CSV file
#' @return Data frame with motif name mappings
load_motif_clustering_data <- function(clustering_file) {

  cat("Loading motif clustering data for dual naming system...\n")

  if (!file.exists(clustering_file)) {
    # Try to find the most recent clustering file
    clustering_pattern <- "../../../out/moca_results/mo_clu/*_vitis_ssr_cwm_clustering_motif_summary.csv"
    clustering_files <- Sys.glob(clustering_pattern)

    if (length(clustering_files) > 0) {
      clustering_file <- sort(clustering_files, decreasing = TRUE)[1]
      cat("  Using found clustering file:", clustering_file, "\n")
    } else {
      cat("  Warning: No clustering file found - using basic naming only\n")
      return(NULL)
    }
  }

  # Load clustering data
  clustering_data <- read.csv(clustering_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(clustering_data), "clustered motifs with sequences\n")

  # The clustering file has:
  # name: "epm_vitis_ssr_0_0_F_4163_13.9_SGCCGCAGCGSCSS" (long format)
  # consensus: "SGCCGCAGCGSCSS" (sequence)

  return(clustering_data)
}

#' Create motif name mapping table from clustering data
#' @param clustering_data Clustering data from load_motif_clustering_data
#' @return Data frame with basic_name, long_name, short_name, sequence mappings
create_motif_name_mapping_from_clustering <- function(clustering_data) {

  if (is.null(clustering_data)) {
    return(NULL)
  }

  # Extract components from long format names
  mapping_data <- clustering_data %>%
    mutate(
      long_name = name,
      sequence = consensus,

      # Extract basic BLAMM name components (remove sequence and score from long name)
      # epm_vitis_ssr_0_0_F_4163_13.9_SGCCGCAGCGSCSS -> epm_vitis_ssr_0_0_F_4163
      basic_name = gsub("_[0-9.]+_[ATGCNKRYWSM]+$", "", name),

      # Create short format for statistical analysis (p0m/p1m format with strand)
      short_name = {
        # Parse pattern_metacluster_strand from basic name
        parts_list <- strsplit(basic_name, "_")

        sapply(parts_list, function(parts) {
          if (length(parts) >= 6) {
            pattern_num <- as.numeric(parts[4])
            metacluster <- as.numeric(parts[5])
            strand <- parts[6]

            # Create p0m/p1m format with strand preserved
            formatted_pattern <- sprintf("p%dm%02d%s", metacluster, pattern_num, strand)
            paste("epm_vitis_ssr", formatted_pattern, sep = "_")
          } else {
            # Fallback to basic name
            basic_name[1]  # Use the first element as fallback
          }
        })
      }
    ) %>%
    select(basic_name, long_name, short_name, sequence, icscore, nsites)

  cat("  Created mapping for", nrow(mapping_data), "motifs\n")
  cat("  Sample mappings:\n")
  cat("    Basic: ", mapping_data$basic_name[1], "\n")
  cat("    Long:  ", mapping_data$long_name[1], "\n")
  cat("    Short: ", mapping_data$short_name[1], "\n")

  return(mapping_data)
}

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Default file paths for Vitis analysis - find files by pattern
# Find most recent q10q90 filtered motif file
motif_pattern <- "../../../out/moca_results/mo_proj/filtering/*_vitis_ssr_q10q90_filtered.csv"
motif_files <- Sys.glob(motif_pattern)
default_motif_file <- if (length(motif_files) > 0) sort(motif_files, decreasing = TRUE)[1] else motif_pattern

# Find most recent prediction file
pred_pattern <- "../../../out/predictions/vitis_*_SSR_deepcre_predict_*.csv"
pred_files <- Sys.glob(pred_pattern)
default_prediction_file <- if (length(pred_files) > 0) sort(pred_files, decreasing = TRUE)[1] else pred_pattern
default_expression_file <- "../../../vitis_data/tpm_counts/vitis_drought_leaf_targets.csv"
default_output_dir <- "../../../out/moca_results/mo_proj/enrichment"

# Clustering file for dual naming system (contains sequences already)
clustering_pattern <- "../../../out/moca_results/mo_clu/*_vitis_ssr_cwm_clustering_motif_summary.csv"
clustering_files <- Sys.glob(clustering_pattern)
default_clustering_file <- if (length(clustering_files) > 0) sort(clustering_files, decreasing = TRUE)[1] else clustering_pattern

# Parse arguments (if not default)
MOTIF_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_motif_file
PREDICTION_FILE <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_prediction_file
EXPRESSION_FILE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_expression_file
OUTPUT_DIR <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_output_dir
CLUSTERING_FILE <- if (length(args) >= 5 && nzchar(args[5])) args[5] else default_clustering_file

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
cat("  Clustering file:", CLUSTERING_FILE, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  Project:", PROJECT_NAME, "\n")
cat("  Date:", DATE_STAMP, "\n")
cat("  Analysis scope: ALL GENES (following original moca_blue approach)\n")
cat("  Naming system: DUAL (long format with sequence + short format for stats)\n\n")

#' Load and validate input datasets for enrichment analysis
#'
#' This function loads motif occurrences, model predictions, and expression data,
#' ensuring proper formatting and compatibility for downstream analysis.
#' Following original moca_blue approach: uses ALL available genes.
#'
#' @param motif_file Path to filtered motif occurrences CSV
#' @param prediction_file Path to model prediction results CSV
#' @param expression_file Path to expression target data CSV
#' @param clustering_file Path to motif clustering summary CSV (contains sequences)
#' @return List containing loaded and processed datasets
load_enrichment_data <- function(motif_file, prediction_file, expression_file, clustering_file = NULL) {

  cat("Loading enrichment analysis datasets (ALL GENES approach)...\n")

  # STEP 1: Load motif occurrences (filtered genomic coordinates + gene associations)
  # This contains: gene_id, motif coordinates, EPM names, regulatory regions
  # Each row = one motif occurrence in one gene's regulatory region
  if (!file.exists(motif_file)) {
    stop("Motif file not found: ", motif_file)
  }
  motif_data <- read.csv(motif_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(motif_data), "motif occurrences\n")

  # DUAL NAMING SYSTEM: Use existing clustering data with sequences
  if (!is.null(clustering_file)) {
    cat("  Creating dual naming system from clustering data...\n")

    # Load clustering data and create mapping
    clustering_data <- load_motif_clustering_data(clustering_file)
    motif_mapping <- create_motif_name_mapping_from_clustering(clustering_data)

    if (!is.null(motif_mapping)) {
      # Join motif data with clustering mapping
      motif_data <- motif_data %>%
        left_join(motif_mapping, by = c("motif" = "basic_name")) %>%
        mutate(
          # Fill in missing values for motifs not in clustering
          long_name = ifelse(is.na(long_name), motif, long_name),
          short_name = ifelse(is.na(short_name), motif, short_name),
          sequence = ifelse(is.na(sequence), "NOT_CLUSTERED", sequence),
          basic_name = motif
        )

      cat("  Dual naming applied to", nrow(motif_data), "motif occurrences\n")
      cat("  Sample mappings:\n")
      cat("    Basic: ", motif_data$basic_name[1], "\n")
      cat("    Long:  ", motif_data$long_name[1], "\n")
      cat("    Short: ", motif_data$short_name[1], "\n")
    } else {
      cat("  No clustering mapping created - using basic naming\n")
      motif_data$long_name <- motif_data$motif
      motif_data$short_name <- motif_data$motif
      motif_data$basic_name <- motif_data$motif
      motif_data$sequence <- "NOT_AVAILABLE"
    }

  } else {
    cat("  No clustering file provided - using basic naming only\n")
    motif_data$long_name <- motif_data$motif
    motif_data$short_name <- motif_data$motif
    motif_data$basic_name <- motif_data$motif
    motif_data$sequence <- "NOT_AVAILABLE"
  }

  # STEP 2: Load model predictions (all available predictions)
  # Expects prediction file to always be available
  if (!file.exists(prediction_file)) {
    stop("Prediction file not found: ", prediction_file,
         "\nPlease ensure model predictions are generated before running enrichment analysis.")
  }

  # Load all available predictions (no filtering)
  prediction_data_raw <- read.csv(prediction_file, stringsAsFactors = FALSE)
  # Standardize column names for downstream analysis
  colnames(prediction_data_raw) <- c("true_targets", "pred_probs", "gene_id")
  prediction_data <- prediction_data_raw
  pred_available <- TRUE
  cat("  Loaded", nrow(prediction_data), "prediction results\n")

  # STEP 3: Load expression data (all available genes)
  # This contains actual expression classes: 0=low, 1=high, 2=medium
  if (!file.exists(expression_file)) {
    stop("Expression file not found: ", expression_file)
  }
  expression_data_raw <- read.csv(expression_file, stringsAsFactors = FALSE)
  # Standardize column names for downstream merging
  colnames(expression_data_raw) <- c("gene_id", "target")

  cat("  Loaded", nrow(expression_data_raw), "expression targets\n")

  # Data validation and summary
  cat("\nData summary (ALL GENES approach):\n")
  cat("  Unique genes in motifs:", length(unique(motif_data$gene_id)), "\n")
  cat("  Unique genes in predictions:", length(unique(prediction_data$gene_id)), "\n")
  cat("  Unique genes in expression:", length(unique(expression_data_raw$gene_id)), "\n")

  # Expression class distribution
  expr_dist <- table(expression_data_raw$target)
  cat("  Expression distribution: Low(0):", expr_dist["0"],
      " High(1):", expr_dist["1"], " Medium(2):", expr_dist["2"], "\n")

  # Gene overlap check
  gene_overlap <- intersect(intersect(unique(motif_data$gene_id), unique(prediction_data$gene_id)), unique(expression_data_raw$gene_id))
  cat("  Genes with all data types:", length(gene_overlap), "\n")

  if (length(gene_overlap) < 100) {
    warning("Very few genes have all data types - check file compatibility")
  }

  return(list(
    motifs = motif_data,
    predictions = prediction_data,
    expression = expression_data_raw,
    pred_available = pred_available
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

  # STEP 6: Create the core analysis dataset by merging three data sources
  # This creates the fundamental unit for enrichment analysis:
  # Each row = one motif occurrence + its gene's expression + model's prediction
  integrated_data <- data_list$motifs %>%
    left_join(data_list$predictions, by = "gene_id") %>%  # Add model predictions
    left_join(data_list$expression, by = "gene_id") %>%   # Add actual expression
    filter(!is.na(true_targets) | !is.na(target))        # Keep only genes with expression data

  # STEP 7: Determine ground truth expression classes
  # Use actual expression data if available, otherwise fall back to prediction targets
  integrated_data$expression_class <- ifelse(
    !is.na(integrated_data$target),
    integrated_data$target,                # Preferred: actual expression classes
    integrated_data$true_targets           # Fallback: targets from predictions
  )

  # STEP 8: Filter out medium expression genes (critical for binary analysis)
  # Original model was trained only on high (1) and low (0) expression genes
  # Medium genes (2) were excluded, so we must exclude them here for consistency
  cat("  Pre-filter dataset size:", nrow(integrated_data), "associations\n")

  integrated_data <- integrated_data %>%
    filter(expression_class %in% c(0, 1))  # Keep only binary classes

  cat("  Post-filter dataset size (medium genes removed):", nrow(integrated_data), "associations\n")

  # STEP 9: Create analysis-ready classifications and identifiers
  integrated_data <- integrated_data %>%
    mutate(
      # Create binary expression labels for analysis
      expr_binary = ifelse(expression_class == 0, "low", "high"),
      pred_binary = ifelse(pred_probs <= 0.5, "low", "high"),

      # Extract metacluster from motif names - THIS IS KEY FOR BIOLOGICAL LOGIC
      # For short format: epm_vitis_ssr_p0m05F or epm_vitis_ssr_p1m05F
      # For BLAMM format: epm_vitis_ssr_pattern_metacluster_strand_count
      # metacluster 0 = expect high expression association (p0m)
      # metacluster 1 = expect low expression association (p1m)
      metacluster = ifelse(grepl("p1m|_1_", motif), "p1m", "p0m"),

      # NAMING SYSTEM: Use clustering-derived short_name if available, otherwise keep original
      epm = if ("short_name" %in% colnames(.)) {
        ifelse(is.na(short_name) | short_name == motif,
               motif,  # Keep original BLAMM name as fallback
               short_name)  # Use clustering-derived short name (p0m/p1m format)
      } else {
        motif  # Keep original BLAMM names
      }
    )

  integrated_data <- integrated_data %>%
    mutate(

      # STEP 10: Calculate model prediction accuracy
      # Simple binary: did the model predict the correct expression class?
      pred_correct = ifelse(expr_binary == pred_binary, "TRUE", "FALSE"),

      # STEP 11: Calculate metacluster-specific prediction performance
      # This checks if motifs behave according to their biological hypothesis:
      # Check motif behavior based on metacluster assignment
      # Use metacluster column (derived from original names) rather than epm pattern matching
      epm_pred_perf = case_when(
        pred_correct == "FALSE" ~ "NA",                                # Model wrong overall
        metacluster == "p0m" & expr_binary == "high" ~ "TRUE_high",   # p0m motifs in high expr (expected)
        metacluster == "p1m" & expr_binary == "low" ~ "TRUE_low",     # p1m motifs in low expr (expected)
        TRUE ~ "NA"                                                    # Other combinations
      ),

      # Simplified metacluster correctness for statistical analysis
      metacluster_correct = case_when(
        metacluster == "p0m" & expr_binary == "high" ~ "TRUE_high",   # p0m behaving correctly
        metacluster == "p1m" & expr_binary == "low" ~ "TRUE_low",     # p1m behaving correctly
        TRUE ~ "FALSE"                                                # Motifs not in expected class
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
#' Computes TPR (true positive rate), chi-square statistics, importance scores, and other
#' statistical indicators for motif-expression associations.
#'
#'
#' @param integrated_data Integrated motif-expression dataset
#' @return Data frame with enrichment metrics for each motif
calculate_motif_enrichment <- function(integrated_data) {

  cat("\nCalculating motif enrichment metrics...\n")

  # STEP 12: Create contingency tables - the foundation of enrichment analysis
  # These 2x2 tables count motif occurrences across different classifications:
  # contingency_expr: motif vs. actual expression (high/low)
  # contingency_pred: motif vs. predicted expression (high/low)
  # contingency_correct: motif vs. prediction accuracy (TRUE/FALSE)
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
      # Detect p1m in short format or _1_ in BLAMM format
      metacluster = ifelse(grepl("p1m|_1_", epm), 1, 0),

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
      # STEP 13: Calculate importance scores - THE CORE OF ENRICHMENT ANALYSIS
      # These scores measure log2 fold-change enrichment relative to background frequency
      #
      # BIOLOGICAL LOGIC:
      # - p0m motifs should be enriched in HIGH expression genes
      # - p1m motifs should be enriched in LOW expression genes
      #
      # MATHEMATICAL LOGIC:
      # For p0m: log2((obs_high/total_high) / (obs_low/total_low))
      # - Positive score = enriched in high expression (biologically correct)
      # - Negative score = depleted from high expression (unexpected)
      #
      # For p1m: log2((obs_low/total_low) / (obs_high/total_high))
      # - Positive score = enriched in low expression (biologically correct)
      # - Negative score = depleted from low expression (unexpected)

      imp_expr_score = ifelse(
        metacluster == 1,  # p1m motifs (expect LOW expression enrichment)
        log2(((expr_low + 0.5) / total_low_expr) / ((expr_high + 0.5) / total_high_expr)),
        log2(((expr_high + 0.5) / total_high_expr) / ((expr_low + 0.5) / total_low_expr))  # p0m motifs (expect HIGH expression enrichment)
      ),

      # Same calculation but for model predictions (validation check)
      imp_pred_score = ifelse(
        metacluster == 1,  # p1m motifs
        log2(((pred_low + 0.5) / total_low_pred) / ((pred_high + 0.5) / total_high_pred)),
        log2(((pred_high + 0.5) / total_high_pred) / ((pred_low + 0.5) / total_low_pred))  # p0m motifs
      ),

      # PRIMARY RANKING METRIC: Use actual expression enrichment for motif ranking
      # This ensures motifs are ranked by biological relevance, not model performance
      importance_score = imp_expr_score
    ) %>%
    ungroup()

  # STEP 14: Statistical significance testing using chi-square tests
  # Tests whether motif distributions deviate significantly from random (50:50) expectation
  # testing independence of motif presence and expression/prediction classes,
  # null hypothesis: motif presence is independent of class

  # Test actual expression associations (primary biological validation)
  enrichment_metrics$chi_expr_pval <- apply(
    enrichment_metrics[, c("expr_low", "expr_high")], 1, function(row) {
      if (sum(row) < 5) return(NA)  # Skip motifs with insufficient data
      result <- suppressWarnings(chisq.test(row, p = c(0.5, 0.5)))  # Test against 50:50 null
      result$p.value
    })

  # Test model prediction associations (technical validation)
  enrichment_metrics$chi_pred_pval <- apply(
    enrichment_metrics[, c("pred_low", "pred_high")], 1, function(row) {
      if (sum(row) < 5) return(NA)
      result <- suppressWarnings(chisq.test(row, p = c(0.5, 0.5)))
      result$p.value
    })

  # Test prediction accuracy associations
  enrichment_metrics$chi_correct_pval <- apply(
    enrichment_metrics[, c("pred_false", "pred_true")], 1, function(row) {
      if (sum(row) < 5) return(NA)
      result <- suppressWarnings(chisq.test(row, p = c(0.5, 0.5)))
      result$p.value
    })

  # STEP 15: Final motif ranking by biological enrichment 
  enrichment_metrics <- enrichment_metrics %>%
    arrange(desc(importance_score)) %>%  # Rank motifs by biological enrichment strength
    mutate(
      enrichment_rank = row_number(),
      significance_expr = ifelse(chi_expr_pval < 0.05, "significant", "not_significant"),
      significance_pred = ifelse(chi_pred_pval < 0.05, "significant", "not_significant")
    )

  cat("  Calculated metrics for", nrow(enrichment_metrics), "motifs\n")
  cat("  Significantly associated with expression:",
      sum(enrichment_metrics$significance_expr == "significant", na.rm = TRUE), "\n")
  cat("  Significantly associated with predictions:",
      sum(enrichment_metrics$significance_pred == "significant", na.rm = TRUE), "\n")

  return(enrichment_metrics)
}

# MAIN EXECUTION
cat("Starting motif enrichment analysis pipeline...\n")

# Load all datasets (no validation filtering)
data_list <- load_enrichment_data(MOTIF_FILE, PREDICTION_FILE, EXPRESSION_FILE, CLUSTERING_FILE)

# Integrate data for analysis
integrated_data <- integrate_enrichment_data(data_list)

# Calculate enrichment metrics
enrichment_metrics <- calculate_motif_enrichment(integrated_data)

# Generate summary statistics by metacluster
enrichment_summary <- enrichment_metrics %>%
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

cat("  Generated summary for", nrow(enrichment_summary), "metaclusters\n")

# Create analysis state for downstream scripts
analysis_state <- list(
  date_stamp = DATE_STAMP,
  total_motifs = nrow(enrichment_metrics),
  total_associations = nrow(integrated_data),
  pred_available = data_list$pred_available,
  analysis_scope = "ALL_GENES"
)

# Write output files
output_base <- file.path(OUTPUT_DIR, paste0(DATE_STAMP, "_", PROJECT_NAME))

# Enrichment metrics CSV
write.csv(enrichment_metrics, paste0(output_base, "_results.csv"), row.names = FALSE)
cat("\nOutput files generated:\n")
cat("  Enrichment results:", basename(paste0(output_base, "_results.csv")), "\n")

# Integrated dataset CSV (includes dual naming system)
# Contains: long_name (with sequence), epm (short format for stats), basic_name, motif_sequence
write.csv(integrated_data, paste0(output_base, "_integrated_data.csv"), row.names = FALSE)
cat("  Integrated data:", basename(paste0(output_base, "_integrated_data.csv")), "\n")

# Summary statistics CSV
write.csv(enrichment_summary, paste0(output_base, "_summary.csv"), row.names = FALSE)
cat("  Summary statistics:", basename(paste0(output_base, "_summary.csv")), "\n")

# Analysis state RDS
saveRDS(analysis_state, paste0(output_base, "_analysis_state.rds"))
cat("  Analysis state:", basename(paste0(output_base, "_analysis_state.rds")), "\n")

# Final summary
cat("\n", rep("=", 60), "\n", sep = "")
cat("MOTIF ENRICHMENT ANALYSIS COMPLETED\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Analysis Summary:\n")
cat("  Total motifs analyzed:", nrow(enrichment_metrics), "\n")
cat("  Significantly associated with expression:",
    sum(enrichment_metrics$significance_expr == "significant", na.rm = TRUE), "\n")
cat("  Significantly associated with predictions:",
    sum(enrichment_metrics$significance_pred == "significant", na.rm = TRUE), "\n")
cat("  Top enriched motif:", enrichment_metrics$epm[1],
    "(importance:", round(enrichment_metrics$imp_expr_score[1], 3), ")\n")
cat("\nMetacluster Summary:\n")
print(enrichment_summary)

cat("\nOutput Directory:", OUTPUT_DIR, "\n")
cat("Files Generated: Enrichment results, integrated data, summaries\n")

cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")