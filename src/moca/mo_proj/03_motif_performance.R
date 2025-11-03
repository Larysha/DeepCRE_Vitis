#!/usr/bin/env Rscript
# Motif Performance Analysis and Visualization
# Part 2 of 2-script motif analysis pipeline (uses outputs from motif_enrichment.R)
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_proj")
######################
#
# SCRIPT OVERVIEW: This script transforms enrichment analysis results into performance metrics
#
# CORE INTEGRATION LOGIC:
# 1. Loads enrichment results from motif_enrichment_03.R as "performance metrics"
# 2. Uses enrichment scores (importance_score) as the PRIMARY ranking system
# 3. Validates model predictions on genes containing biologically-ranked motifs
# 4. Creates dual validation: biological relevance + technical performance
# 5. Generates comprehensive visualizations and performance reports
#
# Biological relevance (enrichment) + Technical accuracy (predictions)
# I.e. bases performance metrics on motif-containing genes
######################

library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(caret)
library(pROC)

# Set working directory to script location
script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
if (nchar(script_dir) > 0) {
  setwd(script_dir)
  cat("Working directory set to:", getwd(), "\n")
}

# Source utility functions
source("../utils.R")

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Default file paths for Vitis analysis - find files by pattern
# Find most recent enrichment results files
enrichment_pattern <- "../../../out/moca_results/mo_proj/enrichment/*_vitis_motif_enrichment_results.csv"
enrichment_files <- Sys.glob(enrichment_pattern)
default_enrichment_results <- if (length(enrichment_files) > 0) sort(enrichment_files, decreasing = TRUE)[1] else enrichment_pattern

integrated_pattern <- "../../../out/moca_results/mo_proj/enrichment/*_vitis_motif_enrichment_integrated_data.csv"
integrated_files <- Sys.glob(integrated_pattern)
default_integrated_data <- if (length(integrated_files) > 0) sort(integrated_files, decreasing = TRUE)[1] else integrated_pattern

state_pattern <- "../../../out/moca_results/mo_proj/enrichment/*_vitis_motif_enrichment_analysis_state.rds"
state_files <- Sys.glob(state_pattern)
default_analysis_state <- if (length(state_files) > 0) sort(state_files, decreasing = TRUE)[1] else state_pattern

default_output_dir <- "../../../out/moca_results/mo_proj/performance"

# Parse arguments
ENRICHMENT_RESULTS <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_enrichment_results
INTEGRATED_DATA <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_integrated_data
ANALYSIS_STATE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_analysis_state
OUTPUT_DIR <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_output_dir

# Analysis parameters
DATE_STAMP <- format(Sys.Date(), "%Y%m%d")
PROJECT_NAME <- "vitis_motif_performance"

# Custom color palette
CUSTOM_COLORS <- c("#671436", "#9e757f", "#052b67", "#176272")
# Generate additional colors if needed
EXTENDED_COLORS <- c(CUSTOM_COLORS, "#ad597e", "#66967b", "#808000", "#657b9e", "#8d77ab")

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

cat("Motif Performance Analysis Parameters:\n")
cat("  Enrichment results:", ENRICHMENT_RESULTS, "\n")
cat("  Integrated data:", INTEGRATED_DATA, "\n") 
cat("  Analysis state:", ANALYSIS_STATE, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  Project:", PROJECT_NAME, "\n")
cat("  Date:", DATE_STAMP, "\n\n")

#' Load enrichment analysis results and transform them into performance metrics
#'
#' CRITICAL INTEGRATION: This is where enrichment becomes "performance"
#' - Loads enrichment metrics → becomes performance rankings
#' - Loads integrated data → enables gene-level model validation
#' - Preserves analysis state → maintains pipeline continuity
#'
#' @param enrichment_file Path to enrichment results CSV
#' @param integrated_file Path to integrated data CSV
#' @param state_file Path to analysis state RDS file
#' @return List containing loaded enrichment data and analysis state
load_performance_data <- function(enrichment_file, integrated_file, state_file) {

  cat("Loading performance analysis datasets from enrichment analysis...\n")

  # STEP 1: Load enrichment metrics as "performance metrics"
  # KEY TRANSFORMATION: enrichment scores become performance rankings
  # Each row = one motif with its biological relevance score (importance_score)
  if (!file.exists(enrichment_file)) {
    stop("Enrichment results file not found: ", enrichment_file)
  }
  enrichment_metrics <- read.csv(enrichment_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(enrichment_metrics), "motif enrichment results\n")

  # STEP 2: Load integrated data for gene-level model validation
  # Each row = one motif occurrence + gene expression + model prediction
  # This enables validation of model performance on genes containing ranked motifs
  if (!file.exists(integrated_file)) {
    stop("Integrated data file not found: ", integrated_file)
  }
  integrated_data <- read.csv(integrated_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(integrated_data), "integrated motif-gene associations\n")
  
  # Load analysis state
  analysis_state <- NULL
  if (file.exists(state_file)) {
    analysis_state <- readRDS(state_file)
    cat("  Loaded analysis state from enrichment run\n")
    cat("    Prediction data available:", analysis_state$pred_available, "\n")
    cat("    Original date stamp:", analysis_state$date_stamp, "\n")
  } else {
    cat("  Warning: Analysis state file not found, using defaults\n")
    analysis_state <- list(
      pred_available = TRUE,  # Default assumption
      date_stamp = DATE_STAMP,
      total_motifs = nrow(enrichment_metrics),
      total_associations = nrow(integrated_data)
    )
  }
  
  # Data validation and summary
  cat("\nData validation:\n")
  cat("  Unique motifs in enrichment results:", length(unique(enrichment_metrics$epm)), "\n")
  cat("  Unique genes in integrated data:", length(unique(integrated_data$gene_id)), "\n")
  cat("  Expression classes in integrated data:")
  if ("expr_binary" %in% colnames(integrated_data)) {
    expr_dist <- table(integrated_data$expr_binary)
    cat(" Low:", expr_dist["low"], " High:", expr_dist["high"], "\n")
  } else {
    cat(" [column not found]\n")
  }
  
  return(list(
    enrichment_metrics = enrichment_metrics,
    integrated_data = integrated_data,
    analysis_state = analysis_state,
    pred_available = analysis_state$pred_available
  ))
}

#' Generate comprehensive performance visualizations based on enrichment rankings
#'
#' VISUALIZATION LOGIC: All plots use enrichment-based performance metrics
#' - Performance metrics = enrichment scores from biological analysis
#' - Rankings based on importance_score (log2 fold-change enrichment)
#' - Model validation overlaid on biologically-ranked motifs
#'
#' @param performance_metrics Data frame with enrichment-based performance metrics
#' @param integrated_data Full integrated dataset with gene-level predictions
#' @param output_dir Directory for saving plots
#' @param pred_available Whether prediction data is available
#' @return List of ggplot objects
generate_performance_plots <- function(performance_metrics, integrated_data, output_dir, pred_available = TRUE) {

  cat("\nGenerating performance visualizations...\n")

  plot_list <- list()

  # PLOT 1: Distribution of biological enrichment scores by motif type
  # Shows how well p0m vs p1m motifs follow their biological expectations
  plot_list$importance_dist <- ggplot(performance_metrics,
                                     aes(x = factor(metacluster), y = importance_score,
                                         fill = factor(metacluster))) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
    scale_fill_manual(values = CUSTOM_COLORS[1:2],
                     labels = c("p0m (High-expr)", "p1m (Low-expr)")) +
    labs(title = "Motif Importance Scores by Metacluster",
         subtitle = "Log2 enrichment in expected expression class",
         x = "Metacluster", y = "Importance Score (log2)",
         fill = "Metacluster") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # PLOT 2: Top performing motifs - RANKED BY BIOLOGICAL ENRICHMENT
  # CRITICAL INSIGHT: This ranks motifs by biological relevance, NOT by significance filtering
  # Shows the top motifs based on how well they follow biological expectations
  top_motifs <- performance_metrics %>%
    filter(!is.na(importance_score)) %>%  # Only exclude motifs with calculation errors
    arrange(desc(importance_score)) %>%   # RANK by enrichment strength 
    slice_head(n = 20)                   # Take top 20 regardless of significance status
  
  # Debug info - investigate motif availability
  cat("  Debug: Total motifs in performance_metrics:", nrow(performance_metrics), "\n")
  cat("  Debug: Motifs with non-NA importance_score:", sum(!is.na(performance_metrics$importance_score)), "\n")
  cat("  Debug: Motifs with 1 occurrence:", sum(performance_metrics$total_occurrences == 1, na.rm = TRUE), "\n")
  cat("  Debug: Motifs with 2+ occurrences:", sum(performance_metrics$total_occurrences >= 2, na.rm = TRUE), "\n")
  cat("  Debug: Motifs with 5+ occurrences:", sum(performance_metrics$total_occurrences >= 5, na.rm = TRUE), "\n")
  cat("  Debug: Final motifs selected for plot:", nrow(top_motifs), "\n")
  
  if (nrow(top_motifs) > 0) {
    cat("  Debug: Occurrence range in top motifs:", 
        min(top_motifs$total_occurrences, na.rm = TRUE), "to", 
        max(top_motifs$total_occurrences, na.rm = TRUE), "\n")
  }
  
  if (nrow(top_motifs) > 0) {
    
    # Create enrichment analysis explanation text
    enrichment_explanation <- paste(
      "ENRICHMENT ANALYSIS INTERPRETATION:",
      "• Importance Score = log2 fold-change enrichment in expected expression class",
      "• p0m motifs (blue): Positive scores = enriched in HIGH expression genes",
      "• p1m motifs (red): Positive scores = enriched in LOW expression genes", 
      "• Higher absolute scores indicate stronger motif-expression associations",
      "• Negative scores indicate motifs depleted from their expected class",
      sep = "\n"
    )
    
    plot_list$top_performers <- ggplot(top_motifs, 
                                      aes(x = reorder(epm, importance_score), 
                                          y = importance_score,
                                          fill = factor(metacluster))) +
      geom_col(alpha = 0.8) +
      scale_fill_manual(values = CUSTOM_COLORS[1:2], 
                       labels = c("p0m (High-expr associated)", "p1m (Low-expr associated)")) +
      coord_flip() +
      labs(title = paste("Top", nrow(top_motifs), "Performing Motifs"),
           subtitle = "Ranked by log2 fold-change enrichment in expected expression class",
           x = "Motif (EPM)", y = "Importance Score (log2 fold-change)",
           fill = "Metacluster",
           caption = enrichment_explanation) +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.caption = element_text(hjust = 0, size = 9, color = "gray30", 
                                       margin = margin(t = 15), lineheight = 1.2))
  } else {
    cat("  Warning: No motifs available for top performers plot\n")
  }
  
  # 3. Dual Importance Scores Comparison
  plot_list$dual_importance <- ggplot(performance_metrics %>% filter(!is.na(imp_expr_score) & !is.na(imp_pred_score)), 
                                     aes(x = imp_expr_score, y = imp_pred_score, 
                                         color = factor(metacluster))) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = CUSTOM_COLORS[1:2], 
                      labels = c("p0m (High-expr)", "p1m (Low-expr)")) +
    labs(title = "Expression vs Prediction Importance Scores",
         subtitle = "Diagonal line shows perfect correlation",
         x = "Expression Importance Score (log2)", 
         y = "Prediction Importance Score (log2)",
         color = "Metacluster") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 4. Performance Distribution by Region
  if ("motif_region" %in% colnames(integrated_data)) {
    regional_performance <- integrated_data %>%
      group_by(motif_region, metacluster) %>%
      summarise(
        mean_performance = mean(pred_correct == "TRUE", na.rm = TRUE),
        count = n(),
        .groups = "drop"
      ) %>%
      filter(count >= 10)  # Only regions with sufficient data
    
    plot_list$regional_performance <- ggplot(regional_performance, 
                                            aes(x = motif_region, y = mean_performance,
                                                fill = factor(metacluster))) +
      geom_col(position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = CUSTOM_COLORS[1:2]) +
      labs(title = "Prediction Performance by Genomic Region",
           subtitle = "Average prediction accuracy",
           x = "Genomic Region", y = "Mean Prediction Accuracy",
           fill = "Metacluster") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")
  }
  
  if (pred_available) {
    # Confusion Matrix
    confusion_data <- integrated_data %>%
      filter(!is.na(pred_binary) & !is.na(expr_binary)) %>%  # Remove rows with NA predictions
      group_by(expr_binary, pred_binary) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(percentage = count / sum(count) * 100)
    
    plot_list$confusion_matrix <- ggplot(confusion_data, 
                                        aes(x = pred_binary, y = expr_binary, 
                                            fill = percentage)) +
      geom_tile(alpha = 0.8) +
      geom_text(aes(label = paste0(count, "\n(", round(percentage, 1), "%)")), 
               color = "white", size = 4) +
      scale_fill_gradient(low = "#a8c4e0", high = CUSTOM_COLORS[3]) +
      labs(title = "Prediction Confusion Matrix",
           subtitle = "Actual vs Predicted Expression Classes",
           x = "Predicted Class", y = "Actual Class",
           fill = "Percentage") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Probability Distribution by Expression Class
    prob_data <- integrated_data %>%
      filter(!is.na(pred_probs) & !is.na(expr_binary))
    
    plot_list$prob_distribution <- ggplot(prob_data, 
                                         aes(x = pred_probs, fill = expr_binary)) +
      geom_density(alpha = 0.6, adjust = 1.2) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
      scale_fill_manual(values = c("low" = CUSTOM_COLORS[2], "high" = CUSTOM_COLORS[1]), 
                       labels = c("low" = "Low Expression", "high" = "High Expression")) +
      labs(title = "Prediction Probability Distribution by Expression Class",
           subtitle = "Red line shows classification threshold (0.5). Colors: Low=pink, High=burgundy",
           x = "Prediction Probability", y = "Density",
           fill = "Actual Expression") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Threshold Analysis - Performance across different thresholds
    if (length(unique(prob_data$pred_probs)) > 10) {  # Only if we have continuous probabilities
      thresholds <- seq(0.1, 0.9, 0.05)
      threshold_results <- data.frame()
      
      for (thresh in thresholds) {
        temp_pred <- ifelse(prob_data$pred_probs >= thresh, "high", "low")
        temp_correct <- ifelse(prob_data$expr_binary == temp_pred, 1, 0)
        
        # Calculate metrics
        tp <- sum(temp_pred == "high" & prob_data$expr_binary == "high")
        fp <- sum(temp_pred == "high" & prob_data$expr_binary == "low")
        tn <- sum(temp_pred == "low" & prob_data$expr_binary == "low")
        fn <- sum(temp_pred == "low" & prob_data$expr_binary == "high")
        
        accuracy <- (tp + tn) / (tp + fp + tn + fn)
        precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
        recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
        f1 <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
        
        threshold_results <- rbind(threshold_results, 
                                 data.frame(threshold = thresh, 
                                           accuracy = accuracy,
                                           precision = precision,
                                           recall = recall,
                                           f1_score = f1))
      }
      
      threshold_long <- threshold_results %>%
        pivot_longer(cols = c(accuracy, precision, recall, f1_score), 
                    names_to = "metric", values_to = "value")
      
      plot_list$threshold_analysis <- ggplot(threshold_long, 
                                            aes(x = threshold, y = value, 
                                                color = metric, linetype = metric)) +
        geom_line(size = 1) +
        geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray", alpha = 0.7) +
        scale_color_manual(values = EXTENDED_COLORS[1:4]) +
        scale_linetype_manual(values = c("solid", "longdash", "dotted", "twodash")) +
        labs(title = "Model Performance Across Classification Thresholds",
             subtitle = "Gray line shows default threshold (0.5)",
             x = "Classification Threshold", y = "Performance Metric",
             color = "Metric", linetype = "Metric") +
        theme_minimal() +
        theme(legend.position = "bottom")
    }
  }
  
  # Save all plots
  for (plot_name in names(plot_list)) {
    if (!is.null(plot_list[[plot_name]])) {
      ggsave(
        filename = file.path(output_dir, paste0(DATE_STAMP, "_", PROJECT_NAME, "_", plot_name, ".pdf")),
        plot = plot_list[[plot_name]],
        width = 10, height = 8, units = "in"
      )
      cat("  Saved plot:", plot_name, "\n")
    }
  }
  
  return(plot_list)
}

#' Calculate model performance on genes containing biologically-ranked motifs
#' This evaluates how well the CNN model predicts expression on genes
#' that contain the motifs ranked by enrichment analysis.
#'
#'
#' @param integrated_data Integrated dataset with motif-gene-expression-prediction associations
#' @param pred_available Whether prediction data is available
#' @return List of model performance metrics
calculate_model_performance <- function(integrated_data, pred_available = TRUE) {

  cat("\nCalculating overall model performance...\n")

  if (!pred_available) {
    cat("  Skipping model performance (predictions not available)\n")
    return(NULL)
  }

  # Test model accuracy on genes containing enriched motifs
  # This validates that the model performs well on the genes we've identified
  # as containing biologically relevant motifs
  actual <- factor(integrated_data$expr_binary, levels = c("low", "high"))
  predicted <- factor(integrated_data$pred_binary, levels = c("low", "high"))

  # Standard classification performance metrics
  cm <- confusionMatrix(predicted, actual, positive = "high")
  
  # Extract key metrics
  performance_summary <- list(
    accuracy = cm$overall['Accuracy'],
    sensitivity = cm$byClass['Sensitivity'],  # TPR for "high"
    specificity = cm$byClass['Specificity'],  # TNR for "high"
    precision = cm$byClass['Pos Pred Value'],
    recall = cm$byClass['Sensitivity'],
    f1_score = cm$byClass['F1'],
    balanced_accuracy = cm$byClass['Balanced Accuracy']
  )
  
  # Calculate AUC if we have probability scores
  if ("pred_probs" %in% colnames(integrated_data)) {
    roc_obj <- roc(integrated_data$expr_binary, integrated_data$pred_probs, 
                   levels = c("low", "high"), direction = "<")
    performance_summary$auc <- auc(roc_obj)
  }
  
  cat("  Model Performance Summary:\n")
  cat("    Accuracy:", round(performance_summary$accuracy, 3), "\n")
  cat("    F1 Score:", round(performance_summary$f1_score, 3), "\n")
  cat("    Precision:", round(performance_summary$precision, 3), "\n")
  cat("    Recall:", round(performance_summary$recall, 3), "\n")
  if (!is.null(performance_summary$auc)) {
    cat("    AUC:", round(performance_summary$auc, 3), "\n")
  }
  
  return(performance_summary)
}

# MAIN EXECUTION: The Core Integration Logic
cat("Starting motif performance analysis pipeline...\n")

# STEP 1: Load enrichment analysis outputs
data_list <- load_performance_data(ENRICHMENT_RESULTS, INTEGRATED_DATA, ANALYSIS_STATE)

# STEP 2: CRITICAL TRANSFORMATION - Enrichment becomes Performance
# This is the key integration point where biological enrichment analysis
# becomes the foundation for performance evaluation
performance_metrics <- data_list$enrichment_metrics  # <-- ENRICHMENT SCORES = PERFORMANCE RANKINGS
integrated_data <- data_list$integrated_data         # <-- GENE-LEVEL DATA FOR MODEL VALIDATION

# STEP 3: Dual validation approach
# A) Technical validation: How well does the model predict on genes with these motifs?
model_performance <- calculate_model_performance(integrated_data, data_list$pred_available)

# B) Biological validation: Visualizations showing enrichment-based rankings
performance_plots <- generate_performance_plots(performance_metrics, integrated_data, OUTPUT_DIR, data_list$pred_available)

# STEP 4: Generate comprehensive summary statistics by biological relevance
# Creates metacluster-level summaries based on enrichment analysis results
metacluster_summary <- performance_metrics %>%
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

# STEP 5: Create final top performers list - RANKED BY BIOLOGICAL ENRICHMENT
# CRITICAL OUTPUT: This demonstrates the core principle of the pipeline
# Motifs are ranked by biological relevance (enrichment), not by model accuracy
top_performers_summary <- performance_metrics %>%
  filter(!is.na(imp_expr_score)) %>%           # Only exclude calculation errors
  arrange(desc(imp_expr_score)) %>%            # SORT BY BIOLOGICAL ENRICHMENT STRENGTH
  select(epm, metacluster, total_occurrences, imp_expr_score, imp_pred_score,
         TPR_expected, TPR_correct, chi_expr_pval, chi_pred_pval, significance_expr) %>%
  head(10)                                     # Top 10 by enrichment (not significance)

# Write output files
output_base <- file.path(OUTPUT_DIR, paste0(DATE_STAMP, "_", PROJECT_NAME))

# Performance metrics CSV
write.csv(performance_metrics, paste0(output_base, "_motif_performance.csv"), row.names = FALSE)
cat("\nOutput files generated:\n")
cat("  Motif performance:", basename(paste0(output_base, "_motif_performance.csv")), "\n")

# Summary statistics CSV
write.csv(metacluster_summary, paste0(output_base, "_metacluster_summary.csv"), row.names = FALSE)
cat("  Metacluster summary:", basename(paste0(output_base, "_metacluster_summary.csv")), "\n")

# Top performers CSV
write.csv(top_performers_summary, paste0(output_base, "_top_performers.csv"), row.names = FALSE)
cat("  Top performers:", basename(paste0(output_base, "_top_performers.csv")), "\n")

# Model performance summary
if (!is.null(model_performance)) {
  model_perf_df <- data.frame(
    metric = names(model_performance),
    value = unlist(model_performance)
  )
  write.csv(model_perf_df, paste0(output_base, "_model_performance.csv"), row.names = FALSE)
  cat("  Model performance:", basename(paste0(output_base, "_model_performance.csv")), "\n")
}

# Final summary
cat("\n", rep("=", 60), "\n", sep = "")
cat("MOTIF PERFORMANCE ANALYSIS COMPLETED\n")
cat(rep("=", 60), "\n\n", sep = "")

# ANALYSIS SUMMARY: Demonstrating the Integration Approach
cat("Analysis Summary (using enrichment results):\n")
cat("  Total motifs analyzed:", nrow(performance_metrics), "\n")
cat("  Significantly associated with expression:",
    sum(performance_metrics$significance_expr == "significant", na.rm = TRUE), "\n")
cat("  Significantly associated with predictions:",
    sum(performance_metrics$significance_pred == "significant", na.rm = TRUE), "\n")

# Show that top motif is ranked by biological enrichment, not significance
cat("  Top performing motif (expression):", performance_metrics$epm[1],
    "(imp_expr:", round(performance_metrics$imp_expr_score[1], 3),
    ", imp_pred:", round(performance_metrics$imp_pred_score[1], 3), ")\n")

cat("\nMetacluster Summary Table:\n")
print(metacluster_summary)

cat("\nTop 10 Performing Motifs:\n") 
print(top_performers_summary)

if (!is.null(model_performance)) {
  cat("\nOverall Model Performance:\n")
  cat("  Accuracy:", round(model_performance$accuracy, 3), "\n")
  cat("  F1 Score:", round(model_performance$f1_score, 3), "\n")
}

cat("\nOutput Directory:", OUTPUT_DIR, "\n")
cat("Files Generated: Performance metrics, summaries, and visualizations\n")

# PIPELINE SUMMARY: Explaining the Complete Integration Approach
cat("\n", rep("=", 60), "\n", sep = "")
cat("ENRICHMENT-PERFORMANCE INTEGRATION SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nPIPELINE LOGIC DEMONSTRATED:\n")
cat("1. Enrichment Analysis → Performance Rankings\n")
cat("   - importance_score (enrichment) = primary ranking metric\n")
cat("   - Motifs ranked by biological relevance, not model accuracy\n")
cat("2. Dual Validation Approach:\n")
cat("   - Biological: Do motifs follow p0m/p1m expression expectations?\n")
cat("   - Technical: Does model predict well on genes with these motifs?\n")
cat("3. No Significance Filtering:\n")
cat("   - All motifs ranked (preserves full spectrum)\n")
cat("   - Statistical significance noted but not used for filtering\n")
cat("4. Integrated Performance:\n")
cat("   - Combines biological relevance + technical validation\n")
cat("   - Ensures both biological meaning and predictive utility\n")

cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Performance analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Based on enrichment analysis from:", data_list$analysis_state$date_stamp, "\n")