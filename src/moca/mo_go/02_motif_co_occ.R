#!/usr/bin/env Rscript
# Motif Co-occurrence Analysis via Venn Diagrams and UpSet Plots
# Visualizes shared target genes between motifs within each metacluster
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_go")
######################

library(dplyr)
library(readr)
library(ggVennDiagram)
library(ggplot2)

# Source utility functions
source("../utils.R")

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Default paths
default_p0m_dir <- "../../../out/moca_results/mo_go/motif_gene_mapping/p0m_motifs/gene_lists"
default_p1m_dir <- "../../../out/moca_results/mo_go/motif_gene_mapping/p1m_motifs/gene_lists"
default_output_dir <- "../../../out/moca_results/mo_go/motif_gene_mapping/venn_diagrams"

# Parse arguments
P0M_DIR <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_p0m_dir
P1M_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_p1m_dir
OUTPUT_DIR <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_output_dir

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

cat("Motif Co-occurrence Analysis Parameters:\n")
cat("  p0m gene lists:", P0M_DIR, "\n")
cat("  p1m gene lists:", P1M_DIR, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n\n")

#' Read gene list CSV files with comment lines
#'
#' Skips the first 6 comment lines (starting with #) and reads gene data
#'
#' @param file_path Path to gene list CSV
#' @return Data frame with gene_id, expression_class, chromosome
read_gene_list <- function(file_path) {

  # Read all lines
  all_lines <- readLines(file_path)

  # Find where data starts (after comment lines and header)
  # Comment lines start with # in the gene_id column
  data_start <- which(!grepl("^#", all_lines) & !grepl("^\"#", all_lines))

  # Find header line
  header_line <- data_start[1] - 1

  # Read the file starting from header
  gene_data <- read.csv(file_path, skip = header_line - 1, stringsAsFactors = FALSE)

  # Remove any remaining comment rows AND duplicate header rows
  gene_data <- gene_data %>%
    filter(!grepl("^#", gene_id)) %>%
    filter(gene_id != "gene_id")  # Remove duplicate header rows

  return(gene_data)
}

#' Load all motif gene lists from a directory
#'
#' @param gene_list_dir Directory containing *_genes.csv files
#' @return List of gene sets, named by motif
load_motif_gene_sets <- function(gene_list_dir) {

  cat("\nLoading gene lists from:", gene_list_dir, "\n")

  # Find all gene list files
  gene_files <- list.files(gene_list_dir, pattern = "_genes\\.csv$", full.names = TRUE)

  if (length(gene_files) == 0) {
    stop("No gene list files found in ", gene_list_dir)
  }

  cat("  Found", length(gene_files), "gene list files\n")

  # Load each file
  gene_sets <- list()
  gene_counts <- data.frame(
    motif_name = character(),
    total_genes = integer(),
    high_expression_genes = integer(),
    low_expression_genes = integer(),
    stringsAsFactors = FALSE
  )

  for (file in gene_files) {
    # Extract motif name from filename
    motif_name <- gsub("_genes\\.csv$", "", basename(file))

    # Read gene list
    genes <- read_gene_list(file)

    # Store gene IDs
    gene_sets[[motif_name]] <- genes$gene_id

    # Count by expression class
    high_count <- sum(genes$expr_binary == "high", na.rm = TRUE)
    low_count <- sum(genes$expr_binary == "low", na.rm = TRUE)

    # Add to summary
    gene_counts <- rbind(gene_counts, data.frame(
      motif_name = motif_name,
      total_genes = nrow(genes),
      high_expression_genes = high_count,
      low_expression_genes = low_count
    ))
  }

  # Sort by total genes (descending)
  gene_counts <- gene_counts %>% arrange(desc(total_genes))

  cat("  Loaded gene sets for", length(gene_sets), "motifs\n")
  cat("  Gene count range:", min(gene_counts$total_genes), "-",
      max(gene_counts$total_genes), "\n")

  return(list(gene_sets = gene_sets, gene_counts = gene_counts))
}

#' Create Venn diagram for top 5 motifs by gene count
#'
#' @param gene_sets Named list of gene ID vectors
#' @param gene_counts Data frame with motif gene counts
#' @param metacluster "p0m" or "p1m"
#' @param output_dir Directory for output files
create_venn_diagram <- function(gene_sets, gene_counts, metacluster, output_dir) {

  cat("\nCreating Venn diagram for", metacluster, "motifs (top 5 by gene count)...\n")

  # Select top 5 motifs by gene count (for visualization clarity)
  top_n <- min(5, nrow(gene_counts))
  top_motifs <- gene_counts$motif_name[1:top_n]
  top_gene_sets <- gene_sets[top_motifs]

  cat("  Selected top", length(top_motifs), "motifs by gene count:\n")
  for (i in 1:length(top_motifs)) {
    cat("    ", i, ".", top_motifs[i], "(",
        length(top_gene_sets[[i]]), "genes )\n")
  }

  # Create Venn diagram
  venn_plot <- ggVennDiagram(top_gene_sets,
                             label = "count",
                             label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    ggtitle(paste0(metacluster, " Motifs: Shared Target Genes\n",
                   ifelse(metacluster == "p0m",
                          "(High Expression / Metacluster 0)",
                          "(Low Expression / Metacluster 1)"))) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  # Save plot
  output_file <- file.path(output_dir, paste0(metacluster, "_motif_gene_venn.pdf"))
  ggsave(output_file, venn_plot, width = 10, height = 8, dpi = 300)

  cat("  Saved Venn diagram:", basename(output_file), "\n")

  return(top_motifs)
}

# UpSet plot function removed - rendering issues
# Use Venn diagrams and summary statistics for co-occurrence analysis

#' Calculate and save overlap statistics for ALL motifs
#'
#' @param gene_sets Named list of gene ID vectors (ALL motifs)
#' @param top_motifs Vector of top motif names used in Venn (top 5)
#' @param metacluster "p0m" or "p1m"
#' @param output_dir Directory for output files
save_venn_summary <- function(gene_sets, top_motifs, metacluster, output_dir) {

  cat("\nGenerating summary statistics for", metacluster, "(ALL motifs)...\n")

  output_file <- file.path(output_dir, paste0(metacluster, "_venn_summary.txt"))

  sink(output_file)

  cat("=" , rep("=", 60), "\n", sep = "")
  cat("MOTIF CO-OCCURRENCE SUMMARY:", metacluster, "\n")
  cat(rep("=", 60), "\n\n", sep = "")

  cat("Metacluster:", ifelse(metacluster == "p0m",
                             "p0m (High Expression / MC0)",
                             "p1m (Low Expression / MC1)"), "\n\n")

  cat("Total Motifs Analyzed (ALL):", length(gene_sets), "\n")
  cat("Motifs in Venn Diagram (top 5 by gene count):", length(top_motifs), "\n\n")

  cat("Top Motifs Selected for Venn (top 5 by gene count):\n")
  for (i in 1:length(top_motifs)) {
    cat("  ", i, ".", top_motifs[i], ":", length(gene_sets[[top_motifs[i]]]), "genes\n")
  }

  # Get all motif names sorted by gene count
  all_motif_names <- names(gene_sets)

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("PAIRWISE OVERLAPS (ALL MOTIFS)\n")
  cat(rep("=", 60), "\n\n", sep = "")

  # Calculate all pairwise overlaps
  for (i in 1:(length(all_motif_names)-1)) {
    for (j in (i+1):length(all_motif_names)) {
      overlap <- length(intersect(gene_sets[[all_motif_names[i]]],
                                  gene_sets[[all_motif_names[j]]]))
      cat("  ", all_motif_names[i], "∩", all_motif_names[j], ":", overlap, "genes\n")
    }
  }

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("UNIQUE GENES PER MOTIF (ALL MOTIFS)\n")
  cat("(genes not shared with any other motif)\n")
  cat(rep("=", 60), "\n\n", sep = "")

  # Calculate genes unique to each motif (not shared with ANY other motif)
  for (motif in all_motif_names) {
    other_motifs <- setdiff(all_motif_names, motif)
    other_genes <- unique(unlist(gene_sets[other_motifs]))
    unique_genes <- setdiff(gene_sets[[motif]], other_genes)
    cat("  ", motif, ":", length(unique_genes), "unique genes\n")
  }

  # Overall statistics
  all_genes <- unique(unlist(gene_sets))
  cat("\nOverall Statistics:\n")
  cat("  Total unique genes across all motifs:", length(all_genes), "\n")
  cat("  Average genes per motif:",
      round(mean(sapply(gene_sets, length)), 1), "\n")
  cat("  Median genes per motif:",
      median(sapply(gene_sets, length)), "\n")

  sink()

  cat("  Saved summary:", basename(output_file), "\n")
}

#' Identify highly-regulated genes (≥3 motifs)
#'
#' @param gene_list_dir Directory containing gene list CSV files
#' @param metacluster "p0m" or "p1m"
#' @param output_dir Directory for output files
#' @param min_motifs Minimum number of motifs to be considered highly-regulated
identify_highly_regulated_genes <- function(gene_list_dir, metacluster, output_dir, min_motifs = 3) {

  cat("\nIdentifying highly-regulated genes for", metacluster, "(≥", min_motifs, "motifs)...\n")

  # Find all gene list files
  gene_files <- list.files(gene_list_dir, pattern = "_genes\\.csv$", full.names = TRUE)

  # Create master gene-motif mapping
  gene_motif_map <- data.frame()

  for (file in gene_files) {
    # Extract motif name and remove strand suffix (F/R)
    motif_name <- gsub("_genes\\.csv$", "", basename(file))
    motif_base <- gsub("[FR]$", "", motif_name)  # Remove strand suffix

    # Read gene list
    genes <- read_gene_list(file)

    # Add to mapping
    if (nrow(genes) > 0) {
      genes$motif_base <- motif_base
      gene_motif_map <- rbind(gene_motif_map,
                              genes %>% select(gene_id, motif_base, expr_binary, chr))
    }
  }

  # Count DISTINCT motifs per gene (strand-agnostic)
  highly_regulated <- gene_motif_map %>%
    group_by(gene_id) %>%
    summarise(
      motif_count = n_distinct(motif_base),
      motif_list = paste(unique(motif_base), collapse = ","),
      expression_class = first(expr_binary),
      chromosome = first(chr),
      .groups = "drop"
    ) %>%
    filter(motif_count >= min_motifs) %>%
    arrange(desc(motif_count))

  # Save output
  output_file <- file.path(output_dir, paste0(metacluster, "_highly_regulated_genes.csv"))
  write.csv(highly_regulated, output_file, row.names = FALSE)

  cat("  Found", nrow(highly_regulated), "genes with ≥", min_motifs, "motifs\n")
  cat("  Motif count range:", min(highly_regulated$motif_count), "-",
      max(highly_regulated$motif_count), "\n")
  cat("  Saved:", basename(output_file), "\n")

  return(highly_regulated)
}

#' Main analysis function for one metacluster
#'
#' @param gene_list_dir Directory with gene list CSV files
#' @param metacluster "p0m" or "p1m"
#' @param output_dir Output directory
analyze_metacluster <- function(gene_list_dir, metacluster, output_dir) {

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("ANALYZING METACLUSTER:", metacluster, "\n")
  cat(rep("=", 60), "\n", sep = "")

  # Load gene sets
  loaded_data <- load_motif_gene_sets(gene_list_dir)
  gene_sets <- loaded_data$gene_sets
  gene_counts <- loaded_data$gene_counts

  # Save gene count summary
  count_file <- file.path(output_dir, paste0(metacluster, "_motif_gene_counts.csv"))
  write.csv(gene_counts, count_file, row.names = FALSE)
  cat("\nSaved gene count summary:", basename(count_file), "\n")

  # Create Venn diagram (top 5 motifs by gene count only)
  top_motifs <- create_venn_diagram(gene_sets, gene_counts, metacluster, output_dir)

  # Save summary statistics (ALL motifs)
  save_venn_summary(gene_sets, top_motifs, metacluster, output_dir)

  # Identify highly-regulated genes (≥3 motifs)
  highly_regulated <- identify_highly_regulated_genes(gene_list_dir, metacluster,
                                                      output_dir, min_motifs = 3)

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("COMPLETED:", metacluster, "\n")
  cat(rep("=", 60), "\n\n", sep = "")

  return(list(
    gene_counts = gene_counts,
    highly_regulated = highly_regulated,
    top_motifs = top_motifs
  ))
}

# MAIN EXECUTION
cat("Starting Motif Co-occurrence Analysis...\n")
cat(rep("=", 70), "\n", sep = "")

# Analyze p0m metacluster
p0m_results <- analyze_metacluster(P0M_DIR, "p0m", OUTPUT_DIR)

# Analyze p1m metacluster
p1m_results <- analyze_metacluster(P1M_DIR, "p1m", OUTPUT_DIR)

# Final summary
cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Output Files Generated:\n")
cat("  Venn Diagrams:\n")
cat("    - p0m_motif_gene_venn.pdf\n")
cat("    - p1m_motif_gene_venn.pdf\n")
cat("  Summary Statistics:\n")
cat("    - p0m_venn_summary.txt\n")
cat("    - p1m_venn_summary.txt\n")
cat("  Gene Count Tables:\n")
cat("    - p0m_motif_gene_counts.csv\n")
cat("    - p1m_motif_gene_counts.csv\n")
cat("  Highly-Regulated Genes:\n")
cat("    - p0m_highly_regulated_genes.csv (", nrow(p0m_results$highly_regulated), " genes)\n", sep = "")
cat("    - p1m_highly_regulated_genes.csv (", nrow(p1m_results$highly_regulated), " genes)\n", sep = "")

cat("\nOutput Directory:", OUTPUT_DIR, "\n")

cat("\nKey Statistics:\n")
cat("  p0m: ", nrow(p0m_results$gene_counts), " motifs (ALL analyzed), ",
    nrow(p0m_results$highly_regulated), " highly-regulated genes (≥3 motifs)\n", sep = "")
cat("  p1m: ", nrow(p1m_results$gene_counts), " motifs (ALL analyzed), ",
    nrow(p1m_results$highly_regulated), " highly-regulated genes (≥3 motifs)\n", sep = "")

cat("\nSelection Criteria:\n")
cat("  - Venn diagrams: Top 5 motifs by GENE COUNT\n")
cat("  - Highly-regulated genes: ≥3 distinct motifs (strand-agnostic)\n")
cat("  - All analyses use ALL motifs (no filtering except Venn top 5)\n")

# Print session info
cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Script completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
