#!/usr/bin/env Rscript
# Refactored from: mo_gene_mapper_v0.1-2.R (original moca_blue project)
#
# This script creates comprehensive gene-motif-function mapping tables,
# focusing on cross-species predictions and functional annotation mapping.
# Extends basic motif analysis to include detailed gene-level mapping
# with functional category assignments.
#
# ORIGINAL LOGIC: Maps genes containing specific motifs to their functional
# categories, creates comprehensive tables for cross-reference analysis,
# and enables identification of functionally coherent gene sets that
# share motif patterns.
#
# Adapted for Vitis vinifera project structure with pattern-based file finding
# and outputs to dedicated mo_go results directory.
#setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_go")
######################

library(tidyr)
library(dplyr)
library(readr)
library(stringr)

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Find input files by pattern
# Motif occurrence data (filtered)
motif_pattern <- "../../../out/moca_results/mo_proj/filtering/*_vitis_ssr_q10q90_filtered.csv"
motif_files <- Sys.glob(motif_pattern)
default_motif_file <- if (length(motif_files) > 0) sort(motif_files, decreasing = TRUE)[1] else motif_pattern

# Model predictions (can be cross-species predictions)
pred_pattern <- "../../../out/predictions/vitis_*_SSR_deepcre_predict_*.csv"
pred_files <- Sys.glob(pred_pattern)
default_pred_file <- if (length(pred_files) > 0) sort(pred_files, decreasing = TRUE)[1] else pred_pattern

# GO annotation file
default_go_file <- "../../../vitis_data/gene_ontology/V_vinifera_ont_converted.gmt"

# Expression/probability data
default_expr_file <- "../../../vitis_data/tpm_counts/vitis_drought_leaf_targets.csv"

# Output directory
default_output_dir <- "../../../out/moca_results/mo_go"

# Parse arguments
MOTIF_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_motif_file
PRED_FILE <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_pred_file
GO_FILE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_go_file
EXPR_FILE <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_expr_file
OUTPUT_DIR <- if (length(args) >= 5 && nzchar(args[5])) args[5] else default_output_dir

# Analysis parameters
DATE_STAMP <- format(Sys.Date(), "%Y%m%d")
PROJECT_NAME <- "vitis_gene_motif_mapper"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

cat("Gene-Motif-Function Mapping Analysis:\n")
cat("  Motif file:", MOTIF_FILE, "\n")
cat("  Predictions file:", PRED_FILE, "\n")
cat("  GO annotation file:", GO_FILE, "\n")
cat("  Expression file:", EXPR_FILE, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  Date:", DATE_STAMP, "\n\n")

#' Use gene IDs as-is (no standardization needed - already handled by 00_convert_gene_ids.R)
use_original_gene_ids <- function(gene_ids) {
  return(gene_ids)
}

#' Load GO annotations from GMT file
load_go_annotations <- function(gmt_file) {
  cat("Loading GO annotations...\n")

  if (!file.exists(gmt_file)) {
    stop("GO annotation file not found: ", gmt_file)
  }

  gmt_lines <- readLines(gmt_file)

  # Pre-allocate list for efficiency
  go_list <- list()
  counter <- 1

  for (line in gmt_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      category_name <- parts[1]
      description <- parts[2]
      genes <- parts[3:length(parts)]

      # Filter non-empty genes
      valid_genes <- genes[nzchar(genes)]

      if (length(valid_genes) > 0) {
        go_list[[counter]] <- data.frame(
          gene_id = valid_genes,
          category_name = category_name,
          description = description,
          stringsAsFactors = FALSE
        )
        counter <- counter + 1
      }
    }
  }

  # Combine all at once instead of repeated rbind
  go_data <- do.call(rbind, go_list)

  # Keep gene IDs as-is (already standardized by 00_convert_gene_ids.R)
  go_data$gene_id_std <- use_original_gene_ids(go_data$gene_id)

  cat("  Loaded", nrow(go_data), "gene-category associations\n")
  return(go_data)
}

#' Load motif occurrence data
load_motif_data <- function(motif_file) {
  cat("Loading motif occurrence data...\n")

  if (!file.exists(motif_file)) {
    stop("Motif file not found: ", motif_file)
  }

  motif_data <- read.csv(motif_file, stringsAsFactors = FALSE)
  motif_data$gene_id_std <- use_original_gene_ids(motif_data$gene_id)

  # Extract motif information
  if (!"epm" %in% colnames(motif_data)) {
    motif_data$epm <- motif_data$motif
  }

  cat("  Loaded", nrow(motif_data), "motif occurrences\n")
  cat("  Unique genes with motifs:", length(unique(motif_data$gene_id_std)), "\n")
  return(motif_data)
}

#' Load prediction data
load_predictions <- function(pred_file) {
  cat("Loading model predictions...\n")

  if (!file.exists(pred_file)) {
    stop("Prediction file not found: ", pred_file)
  }

  pred_data <- read.csv(pred_file, stringsAsFactors = FALSE)

  # Handle different column naming conventions
  if ("genes" %in% colnames(pred_data) && !"gene_id" %in% colnames(pred_data)) {
    pred_data$gene_id <- pred_data$genes
  }
  if ("pred_probs" %in% colnames(pred_data) && !"prob" %in% colnames(pred_data)) {
    pred_data$prob <- pred_data$pred_probs
  }

  pred_data$gene_id_std <- use_original_gene_ids(pred_data$gene_id)

  # Create probability classes
  pred_data$prob_class <- ifelse(pred_data$prob > 0.5, "high", "low")

  cat("  Loaded", nrow(pred_data), "predictions\n")
  return(pred_data)
}

#' Load expression/probability data
load_expression_data <- function(expr_file) {
  cat("Loading expression data...\n")

  if (!file.exists(expr_file)) {
    stop("Expression file not found: ", expr_file)
  }

  expr_data <- read.csv(expr_file, stringsAsFactors = FALSE)
  expr_data$gene_id_std <- use_original_gene_ids(expr_data$gene_id)

  # Create expression classes
  if ("target_class" %in% colnames(expr_data)) {
    expr_data$expr_class <- ifelse(expr_data$target_class == 1, "high", "low")
  } else if ("target" %in% colnames(expr_data)) {
    # Handle target column - assuming 1 = high expression, 0 = low expression, 2 = medium (exclude)
    expr_data$expr_class <- ifelse(expr_data$target == 1, "high",
                                   ifelse(expr_data$target == 0, "low", NA))
    # Remove medium expression genes (target == 2) as they weren't used in training
    expr_data <- expr_data[!is.na(expr_data$expr_class), ]
  } else if ("tpm" %in% colnames(expr_data)) {
    median_expr <- median(expr_data$tpm, na.rm = TRUE)
    expr_data$expr_class <- ifelse(expr_data$tpm >= median_expr, "high", "low")
  }

  cat("  Loaded", nrow(expr_data), "expression measurements\n")
  return(expr_data)
}

#' Create comprehensive gene-motif-function mapping
#'
#' Following original mo_gene_mapper logic to create detailed mapping tables
#'
#' @param motif_data Motif occurrence data
#' @param pred_data Model predictions
#' @param expr_data Expression data
#' @param go_data GO annotations
#' @return Comprehensive mapping table
create_gene_motif_mapping <- function(motif_data, pred_data, expr_data, go_data) {

  cat("Creating comprehensive gene-motif-function mapping...\n")

  # Step 1: Create base gene mapping (all genes with any data)
  all_genes <- unique(c(motif_data$gene_id_std, pred_data$gene_id_std,
                       expr_data$gene_id_std, go_data$gene_id_std))

  base_mapping <- data.frame(gene_id_std = all_genes, stringsAsFactors = FALSE)

  # Step 2: Add expression information
  base_mapping <- base_mapping %>%
    left_join(expr_data %>% select(gene_id_std, expr_class), by = "gene_id_std")

  # Step 3: Add prediction information
  base_mapping <- base_mapping %>%
    left_join(pred_data %>% select(gene_id_std, prob, prob_class), by = "gene_id_std")

  # Step 4: Add GO functional annotations (preserve all categories per gene)
  go_summary <- go_data %>%
    group_by(gene_id_std) %>%
    summarise(
      go_category_count = n(),
      go_categories = paste(unique(category_name), collapse = ";"),
      go_descriptions = paste(unique(description), collapse = ";"),
      primary_category = first(category_name),
      primary_description = first(description),
      .groups = "drop"
    )

  base_mapping <- base_mapping %>%
    left_join(go_summary, by = "gene_id_std")

  # Step 5: Add motif information (genes may have multiple motifs)
  # Create a summary of motifs per gene
  gene_motif_summary <- motif_data %>%
    group_by(gene_id_std) %>%
    summarise(
      total_motifs = n(),
      unique_motifs = length(unique(epm)),
      motif_list = paste(unique(epm), collapse = ";"),
      has_p0m_motif = any(grepl("p0m", epm)),
      has_p1m_motif = any(grepl("p1m", epm)),
      .groups = "drop"
    )

  base_mapping <- base_mapping %>%
    left_join(gene_motif_summary, by = "gene_id_std")

  # Add flags for data availability
  base_mapping <- base_mapping %>%
    mutate(
      has_expression = !is.na(expr_class),
      has_prediction = !is.na(prob),
      has_go_annotation = !is.na(primary_category),
      has_motifs = !is.na(total_motifs),

      # Classification of genes
      data_completeness = case_when(
        has_expression & has_prediction & has_go_annotation & has_motifs ~ "complete",
        has_expression & has_prediction & has_motifs ~ "no_go",
        has_expression & has_go_annotation & has_motifs ~ "no_prediction",
        has_motifs & has_go_annotation ~ "motif_go_only",
        TRUE ~ "incomplete"
      )
    )

  cat("  Created mapping for", nrow(base_mapping), "genes\n")
  cat("  Complete data:", sum(base_mapping$data_completeness == "complete"), "genes\n")
  cat("  Genes with motifs:", sum(base_mapping$has_motifs, na.rm = TRUE), "\n")
  cat("  Genes with GO annotations:", sum(base_mapping$has_go_annotation, na.rm = TRUE), "\n")

  return(base_mapping)
}

#' Create motif-function association tables
#'
#' Following original logic to create detailed motif-GO association tables
#'
#' @param gene_mapping Complete gene mapping
#' @param motif_data Original motif data
#' @param go_data GO annotation data
#' @return Motif-function association analysis
analyze_motif_function_associations <- function(gene_mapping, motif_data, go_data) {

  cat("Analyzing motif-function associations...\n")

  # Create detailed motif-GO associations (each motif occurrence)
  # Select available columns dynamically
  available_cols <- c("gene_id_std", "epm", "category_name", "description")
  optional_cols <- c("motif_region", "score")

  # Add optional columns if they exist
  for (col in optional_cols) {
    if (col %in% colnames(motif_data) || col %in% colnames(go_data)) {
      available_cols <- c(available_cols, col)
    }
  }

  motif_go_detailed <- motif_data %>%
    left_join(go_data, by = "gene_id_std") %>%
    filter(!is.na(category_name)) %>%
    select(all_of(available_cols))

  # Summarize associations by motif and GO category
  summary_base <- motif_go_detailed %>%
    group_by(epm, category_name, description) %>%
    summarise(
      gene_count = length(unique(gene_id_std)),
      total_occurrences = n(),
      .groups = "drop"
    )

  # Add mean score if score column exists
  if ("score" %in% colnames(motif_go_detailed)) {
    motif_go_summary <- summary_base %>%
      left_join(
        motif_go_detailed %>%
          group_by(epm, category_name, description) %>%
          summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop"),
        by = c("epm", "category_name", "description")
      )
  } else {
    motif_go_summary <- summary_base
  }

  motif_go_summary <- motif_go_summary %>% arrange(desc(gene_count))

  # Analyze motif enrichment in GO categories
  # Calculate background frequencies
  total_genes_with_go <- length(unique(go_data$gene_id_std))
  go_category_sizes <- go_data %>%
    group_by(category_name) %>%
    summarise(category_size = length(unique(gene_id_std)), .groups = "drop")

  motif_totals <- motif_data %>%
    group_by(epm) %>%
    summarise(total_genes_with_motif = length(unique(gene_id_std)), .groups = "drop")

  # Calculate enrichment statistics
  enrichment_analysis <- motif_go_summary %>%
    left_join(go_category_sizes, by = "category_name") %>%
    left_join(motif_totals, by = "epm") %>%
    mutate(
      # Expected frequency under null hypothesis
      expected_genes = (total_genes_with_motif * category_size) / total_genes_with_go,

      # Enrichment fold change
      enrichment_fold = gene_count / pmax(expected_genes, 1),

      # Simple significance test (could be enhanced with proper statistics)
      is_enriched = gene_count >= 3 & enrichment_fold >= 2,
      is_depleted = gene_count >= 3 & enrichment_fold <= 0.5
    ) %>%
    arrange(desc(enrichment_fold))

  # Functional coherence analysis by motif type
  motif_type_functions <- motif_go_summary %>%
    mutate(motif_type = ifelse(grepl("p0m", epm), "p0m_high_expr", "p1m_low_expr")) %>%
    group_by(motif_type, category_name, description) %>%
    summarise(
      unique_motifs = length(unique(epm)),
      total_genes = sum(gene_count),
      mean_enrichment = mean(gene_count),
      .groups = "drop"
    ) %>%
    filter(unique_motifs >= 2) %>%  # Categories with multiple motifs of same type
    arrange(desc(total_genes))

  cat("  Detailed motif-GO associations:", nrow(motif_go_detailed), "\n")
  cat("  Unique motif-GO combinations:", nrow(motif_go_summary), "\n")
  cat("  Enriched combinations:", sum(enrichment_analysis$is_enriched), "\n")
  cat("  Depleted combinations:", sum(enrichment_analysis$is_depleted), "\n")

  return(list(
    detailed = motif_go_detailed,
    summary = motif_go_summary,
    enrichment = enrichment_analysis,
    motif_type_functions = motif_type_functions
  ))
}

#' Create gene sets for functional analysis
#'
#' Generate gene lists for downstream functional enrichment analysis
#'
#' @param gene_mapping Complete gene mapping
#' @param associations Motif-function associations
#' @return Gene sets organized by different criteria
create_functional_gene_sets <- function(gene_mapping, associations) {

  cat("Creating functional gene sets...\n")

  # Gene sets by motif presence
  motif_gene_sets <- gene_mapping %>%
    filter(has_motifs) %>%
    group_by(motif_list) %>%
    summarise(
      gene_set = paste(unique(gene_id_std), collapse = ","),
      gene_count = length(unique(gene_id_std)),
      .groups = "drop"
    ) %>%
    filter(gene_count >= 5)  # Only sets with sufficient genes

  # Gene sets by motif type and GO category
  motif_go_gene_sets <- associations$detailed %>%
    group_by(epm, category_name) %>%
    summarise(
      gene_set = paste(unique(gene_id_std), collapse = ","),
      gene_count = length(unique(gene_id_std)),
      description = first(description),
      .groups = "drop"
    ) %>%
    filter(gene_count >= 3)

  # High-confidence gene sets (multiple lines of evidence)
  high_confidence_sets <- gene_mapping %>%
    filter(data_completeness == "complete") %>%
    group_by(primary_category, has_p0m_motif, has_p1m_motif) %>%
    summarise(
      gene_set = paste(unique(gene_id_std), collapse = ","),
      gene_count = length(unique(gene_id_std)),
      description = first(primary_description),
      .groups = "drop"
    ) %>%
    filter(gene_count >= 3)

  cat("  Motif-based gene sets:", nrow(motif_gene_sets), "\n")
  cat("  Motif-GO gene sets:", nrow(motif_go_gene_sets), "\n")
  cat("  High-confidence sets:", nrow(high_confidence_sets), "\n")

  return(list(
    motif_sets = motif_gene_sets,
    motif_go_sets = motif_go_gene_sets,
    high_confidence_sets = high_confidence_sets
  ))
}

# Main execution
cat("Starting gene-motif-function mapping analysis...\n")

# Load all datasets
go_annotations <- load_go_annotations(GO_FILE)
motif_data <- load_motif_data(MOTIF_FILE)
pred_data <- load_predictions(PRED_FILE)
expr_data <- load_expression_data(EXPR_FILE)

# Create comprehensive mapping
gene_mapping <- create_gene_motif_mapping(motif_data, pred_data, expr_data, go_annotations)

# Analyze motif-function associations
associations <- analyze_motif_function_associations(gene_mapping, motif_data, go_annotations)

# Create functional gene sets
gene_sets <- create_functional_gene_sets(gene_mapping, associations)

# Write outputs
output_base <- file.path(OUTPUT_DIR, paste0(DATE_STAMP, "_", PROJECT_NAME))

cat("\nWriting output files...\n")

# Complete gene mapping table
write.csv(gene_mapping,
          paste0(output_base, "_gene_mapping.csv"),
          row.names = FALSE)
cat("  Gene mapping:", basename(paste0(output_base, "_gene_mapping.csv")), "\n")

# Motif-GO association tables
write.csv(associations$summary,
          paste0(output_base, "_motif_go_summary.csv"),
          row.names = FALSE)
cat("  Motif-GO summary:", basename(paste0(output_base, "_motif_go_summary.csv")), "\n")

write.csv(associations$enrichment,
          paste0(output_base, "_motif_go_enrichment.csv"),
          row.names = FALSE)
cat("  Motif-GO enrichment:", basename(paste0(output_base, "_motif_go_enrichment.csv")), "\n")

write.csv(associations$motif_type_functions,
          paste0(output_base, "_motif_type_functions.csv"),
          row.names = FALSE)
cat("  Motif type functions:", basename(paste0(output_base, "_motif_type_functions.csv")), "\n")

# Gene sets for functional analysis
write.csv(gene_sets$motif_go_sets,
          paste0(output_base, "_motif_go_gene_sets.csv"),
          row.names = FALSE)
cat("  Motif-GO gene sets:", basename(paste0(output_base, "_motif_go_gene_sets.csv")), "\n")

write.csv(gene_sets$high_confidence_sets,
          paste0(output_base, "_high_confidence_gene_sets.csv"),
          row.names = FALSE)
cat("  High-confidence gene sets:", basename(paste0(output_base, "_high_confidence_gene_sets.csv")), "\n")

# Detailed associations (large file)
write.csv(associations$detailed,
          paste0(output_base, "_detailed_motif_go_associations.csv"),
          row.names = FALSE)
cat("  Detailed associations:", basename(paste0(output_base, "_detailed_motif_go_associations.csv")), "\n")

# Final summary
cat("\n", rep("=", 60), "\n")
cat("GENE-MOTIF-FUNCTION MAPPING COMPLETED\n")
cat(rep("=", 60), "\n")
cat("Analysis Summary:\n")
cat("  Total genes mapped:", nrow(gene_mapping), "\n")
cat("  Genes with complete data:", sum(gene_mapping$data_completeness == "complete"), "\n")
cat("  Motif-GO associations:", nrow(associations$summary), "\n")
cat("  Enriched motif-GO combinations:", sum(associations$enrichment$is_enriched), "\n")
cat("  Functional gene sets created:", nrow(gene_sets$motif_go_sets), "\n")

cat("\nOutput Directory:", OUTPUT_DIR, "\n")
cat("Files Generated: Comprehensive mapping tables and functional gene sets\n")

cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")