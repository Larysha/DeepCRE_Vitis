#!/usr/bin/env Rscript
# Exploratory analysis of expression class distribution in clade shared genes
# Generates contingency tables showing low/high expression counts and percentages
#
# Usage:
#   Rscript exploratory_analysis.R [output_file]
#
# Example:
#   Rscript exploratory_analysis.R ../../../out/moca_results/mo_go/clade_expression_contingency.tsv
######################

library(dplyr)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default output file
OUTPUT_FILE <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  "../../../out/moca_results/mo_go/clade_expression_contingency.tsv"
}

# Paths
CLADE_BASE_DIR <- "../../../out/moca_results/mo_go/clade_shared_genes"

cat("=" , rep("=", 70), "\n", sep = "")
cat("CLADE EXPRESSION CLASS CONTINGENCY ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Parameters:\n")
cat("  Clade directory:", CLADE_BASE_DIR, "\n")
cat("  Output file:", OUTPUT_FILE, "\n\n")

# Check paths exist
if (!dir.exists(CLADE_BASE_DIR)) {
  stop("Clade directory not found: ", CLADE_BASE_DIR)
}

#' Read shared genes CSV with comment lines
#'
#' @param file_path Path to shared genes CSV
#' @return Data frame with gene_id and expression_class
read_shared_genes <- function(file_path) {

  # Read all lines
  all_lines <- readLines(file_path)

  # Find first non-comment line (header)
  non_comment_lines <- which(!grepl("^#", all_lines) & !grepl("^\"#", all_lines))

  if (length(non_comment_lines) == 0) {
    stop("No data found in file: ", file_path)
  }

  # First non-comment line is the header
  header_line <- non_comment_lines[1]
  skip_lines <- header_line - 1

  # Read CSV
  gene_data <- read.csv(file_path, skip = skip_lines, stringsAsFactors = FALSE)

  # Check required columns
  if (!"gene_id" %in% colnames(gene_data)) {
    stop("No 'gene_id' column found in file: ", file_path)
  }

  if (!"expression_class" %in% colnames(gene_data)) {
    stop("No 'expression_class' column found in file: ", file_path)
  }

  return(gene_data)
}

#' Analyze expression class distribution for a clade
#'
#' @param clade_name Name of the clade
#' @param gene_file Path to shared genes CSV
#' @return Data frame with contingency table row
analyze_clade <- function(clade_name, gene_file) {

  # Read gene data
  genes <- read_shared_genes(gene_file)

  # Count expression classes
  total_genes <- nrow(genes)

  if (total_genes == 0) {
    return(data.frame(
      clade = clade_name,
      total_genes = 0,
      count_low = 0,
      count_high = 0,
      percent_low = 0,
      percent_high = 0,
      stringsAsFactors = FALSE
    ))
  }

  # Count by expression class
  expr_counts <- table(genes$expression_class)

  count_low <- ifelse("low" %in% names(expr_counts), expr_counts["low"], 0)
  count_high <- ifelse("high" %in% names(expr_counts), expr_counts["high"], 0)

  # Calculate percentages
  percent_low <- (count_low / total_genes) * 100
  percent_high <- (count_high / total_genes) * 100

  return(data.frame(
    clade = clade_name,
    total_genes = total_genes,
    count_low = as.integer(count_low),
    count_high = as.integer(count_high),
    percent_low = round(percent_low, 2),
    percent_high = round(percent_high, 2),
    stringsAsFactors = FALSE
  ))
}

#' Process all clades in a metacluster directory
#'
#' @param clade_dir Path to clade directory (e.g., p0_clades)
#' @param metacluster_name Name for grouping (e.g., "p0m")
#' @return Data frame with contingency table rows
process_metacluster <- function(clade_dir, metacluster_name) {

  cat("\n", rep("-", 70), "\n", sep = "")
  cat("Processing", metacluster_name, "clades\n")
  cat(rep("-", 70), "\n", sep = "")

  # Find all *_shared_genes.csv files
  gene_files <- list.files(clade_dir,
                           pattern = "_shared_genes\\.csv$",
                           full.names = TRUE)

  if (length(gene_files) == 0) {
    cat("  Warning: No shared gene files found in", clade_dir, "\n")
    return(data.frame())
  }

  cat("  Found", length(gene_files), "clade gene lists\n\n")

  results <- list()

  for (gene_file in gene_files) {
    # Extract clade name from filename
    clade_name <- basename(gene_file)
    clade_name <- gsub("_shared_genes\\.csv$", "", clade_name)

    cat("  Processing:", clade_name, "\n")

    result <- tryCatch({
      analyze_clade(clade_name, gene_file)
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
      NULL
    })

    if (!is.null(result)) {
      results[[clade_name]] <- result
      cat("    Total:", result$total_genes,
          "| Low:", result$count_low,
          "(", result$percent_low, "%)",
          "| High:", result$count_high,
          "(", result$percent_high, "%)\n")
    }
  }

  if (length(results) == 0) {
    return(data.frame())
  }

  combined <- do.call(rbind, results)

  # Add metacluster column
  combined$metacluster <- metacluster_name

  # Reorder columns
  combined <- combined %>%
    select(metacluster, clade, total_genes, count_low, count_high,
           percent_low, percent_high)

  return(combined)
}

# MAIN EXECUTION
cat("Starting analysis...\n")

all_results <- list()

# Process p0 clades
p0_dir <- file.path(CLADE_BASE_DIR, "p0_clades")
if (dir.exists(p0_dir)) {
  p0_results <- process_metacluster(p0_dir, "p0m")
  if (nrow(p0_results) > 0) {
    all_results[["p0m"]] <- p0_results
  }
}

# Process p1 clades
p1_dir <- file.path(CLADE_BASE_DIR, "p1_clades")
if (dir.exists(p1_dir)) {
  p1_results <- process_metacluster(p1_dir, "p1m")
  if (nrow(p1_results) > 0) {
    all_results[["p1m"]] <- p1_results
  }
}

if (length(all_results) == 0) {
  stop("No clade results found")
}

# Combine all results
cat("\n", rep("=", 70), "\n", sep = "")
cat("COMBINING RESULTS\n")
cat(rep("=", 70), "\n\n", sep = "")

combined_table <- do.call(rbind, all_results)
rownames(combined_table) <- NULL

# Sort by metacluster and clade
combined_table <- combined_table %>%
  arrange(metacluster, clade)

cat("Total clades analyzed:", nrow(combined_table), "\n")

# Print summary statistics
cat("\nSummary by metacluster:\n")
summary_stats <- combined_table %>%
  group_by(metacluster) %>%
  summarise(
    n_clades = n(),
    total_genes = sum(total_genes),
    total_low = sum(count_low),
    total_high = sum(count_high),
    avg_percent_low = mean(percent_low),
    avg_percent_high = mean(percent_high),
    .groups = "drop"
  )

print(summary_stats, row.names = FALSE)

# Save results
cat("\nSaving results...\n")
cat("  Output file:", OUTPUT_FILE, "\n")

# Create output directory if needed
output_dir <- dirname(OUTPUT_FILE)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

write.table(combined_table,
            OUTPUT_FILE,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

cat("\nContingency table saved successfully!\n")

# Print the table to screen
cat("\n", rep("=", 70), "\n", sep = "")
cat("CONTINGENCY TABLE\n")
cat(rep("=", 70), "\n\n", sep = "")

print(combined_table, row.names = FALSE)

# Session info
cat("\n", rep("=", 70), "\n", sep = "")
cat("Session Information:\n")
cat("R version:", R.version.string, "\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
