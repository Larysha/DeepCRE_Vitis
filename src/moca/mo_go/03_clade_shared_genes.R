#!/usr/bin/env Rscript
# Extract shared genes between motifs within hierarchical clustering clades
# Uses dendrogram structures from 01_motif_gene_mapper.R
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_go")
######################

library(dplyr)
library(readr)

# Source utility functions
source("../utils.R")

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Default paths
default_gene_list_dir_p0m <- "../../../out/moca_results/mo_go/motif_gene_mapping/p0m_motifs/gene_lists"
default_gene_list_dir_p1m <- "../../../out/moca_results/mo_go/motif_gene_mapping/p1m_motifs/gene_lists"
default_output_dir <- "../../../out/moca_results/mo_go/clade_shared_genes"

# Parse arguments
GENE_LIST_DIR_P0M <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_gene_list_dir_p0m
GENE_LIST_DIR_P1M <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_gene_list_dir_p1m
OUTPUT_DIR <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_output_dir

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

cat("Clade Shared Genes Extraction Parameters:\n")
cat("  p0m gene lists:", GENE_LIST_DIR_P0M, "\n")
cat("  p1m gene lists:", GENE_LIST_DIR_P1M, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n\n")


# Define clade structures from dendrogram analysis
p1m_clades <- list(
  p1m_l0_m6_blue = c("epm_vitis_ssr_p1m11", "epm_vitis_ssr_p1m14", 
                     "epm_vitis_ssr_p1m09", "epm_vitis_ssr_p1m04", 
                     "epm_vitis_ssr_p1m02", "epm_vitis_ssr_p1m03"),
  
  p1m_l2_m4_darkgreen = c("epm_vitis_ssr_p1m09", "epm_vitis_ssr_p1m04", 
                          "epm_vitis_ssr_p1m02", "epm_vitis_ssr_p1m03"),
  
  p1m_l0_m11_red = c("epm_vitis_ssr_p1m05", "epm_vitis_ssr_p1m13", 
                     "epm_vitis_ssr_p1m06", "epm_vitis_ssr_p1m10", 
                     "epm_vitis_ssr_p1m00", "epm_vitis_ssr_p1m15",
                     "epm_vitis_ssr_p1m01", "epm_vitis_ssr_p1m12", 
                     "epm_vitis_ssr_p1m07", "epm_vitis_ssr_p1m08", 
                     "epm_vitis_ssr_p1m16"),
  
  p1m_l1_m6_pink = c("epm_vitis_ssr_p1m05", "epm_vitis_ssr_p1m13", 
                     "epm_vitis_ssr_p1m06", "epm_vitis_ssr_p1m10", 
                     "epm_vitis_ssr_p1m00", "epm_vitis_ssr_p1m15"),
  
  p1m_l1_m5_lightgreen = c("epm_vitis_ssr_p1m01", "epm_vitis_ssr_p1m12", 
                           "epm_vitis_ssr_p1m07", "epm_vitis_ssr_p1m08", 
                           "epm_vitis_ssr_p1m16"),
  
  p1m_l2_m2_darkblue = c("epm_vitis_ssr_p1m06", "epm_vitis_ssr_p1m10")
)

p0m_clades <- list(
  p0m_l0_m4_red = c("epm_vitis_ssr_p0m06", "epm_vitis_ssr_p0m04", 
                    "epm_vitis_ssr_p0m02", "epm_vitis_ssr_p0m05"),
  
  p0m_l0_m8_grey = c("epm_vitis_ssr_p0m03", "epm_vitis_ssr_p0m11", 
                     "epm_vitis_ssr_p0m07", "epm_vitis_ssr_p0m01", 
                     "epm_vitis_ssr_p0m08", "epm_vitis_ssr_p0m10", 
                     "epm_vitis_ssr_p0m00", "epm_vitis_ssr_p0m09"),
  
  p0m_l1_m2_yellow = c("epm_vitis_ssr_p0m03", "epm_vitis_ssr_p0m11"),
  
  p0m_l1_m3_green = c("epm_vitis_ssr_p0m07", "epm_vitis_ssr_p0m01", 
                      "epm_vitis_ssr_p0m08"),
  
  p0m_l1_m3_pink = c("epm_vitis_ssr_p0m10", "epm_vitis_ssr_p0m00", 
                     "epm_vitis_ssr_p0m09")
)

#' Read gene list CSV with comment lines
#'
#' @param file_path Path to gene list CSV
#' @return Data frame with gene_id, expr_binary columns
read_gene_list <- function(file_path) {

  all_lines <- readLines(file_path)
  data_start <- which(!grepl("^#", all_lines) & !grepl("^\"#", all_lines))
  header_line <- data_start[1] - 1

  gene_data <- read.csv(file_path, skip = header_line - 1, stringsAsFactors = FALSE)

  # Remove comment rows and duplicate headers
  gene_data <- gene_data %>%
    filter(!grepl("^#", gene_id)) %>%
    filter(gene_id != "gene_id")

  return(gene_data)
}

#' Find shared genes across all motifs in a clade
#'
#' @param clade_motifs Vector of motif names in the clade
#' @param gene_list_dir Directory containing gene list CSV files
#' @return Data frame with shared genes
find_clade_shared_genes <- function(clade_motifs, gene_list_dir) {

  # Load gene sets for each motif in the clade
  gene_sets <- list()

  for (motif in clade_motifs) {
    gene_file <- file.path(gene_list_dir, paste0(motif, "_genes.csv"))

    if (!file.exists(gene_file)) {
      cat("  Warning: Gene file not found for", motif, "\n")
      next
    }

    genes <- read_gene_list(gene_file)
    gene_sets[[motif]] <- genes
  }

  if (length(gene_sets) == 0) {
    cat("  Error: No gene files found for this clade\n")
    return(data.frame())
  }

  # Find intersection of all gene sets
  shared_gene_ids <- gene_sets[[1]]$gene_id

  for (i in 2:length(gene_sets)) {
    shared_gene_ids <- intersect(shared_gene_ids, gene_sets[[i]]$gene_id)
  }

  # Get full gene information for shared genes (from first motif)
  shared_genes <- gene_sets[[1]] %>%
    filter(gene_id %in% shared_gene_ids) %>%
    select(gene_id, expr_binary) %>%
    rename(expression_class = expr_binary) %>%
    arrange(gene_id)

  return(shared_genes)
}

#' Create shared genes CSV for a clade
#'
#' @param clade_name Name of the clade
#' @param clade_motifs Vector of motif names
#' @param shared_genes Data frame with shared genes
#' @param output_dir Output directory
create_clade_csv <- function(clade_name, clade_motifs, shared_genes, output_dir) {

  # Create clade subdirectory
  clade_dir <- file.path(output_dir, clade_name)
  if (!dir.exists(clade_dir)) {
    dir.create(clade_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Prepare output file
  output_file <- file.path(clade_dir, paste0(clade_name, "_shared_genes.csv"))

  # Create complete output with header comments
  header1 <- paste0("# Clade: ", clade_name, " | Motifs: ",
                    paste(clade_motifs, collapse = ", "))

  header2 <- paste0("# Total motifs in clade: ", length(clade_motifs),
                    " | Shared genes: ", nrow(shared_genes))

  # Write complete file manually to avoid append warnings
  file_conn <- file(output_file, open = "w")

  # Write header comments
  writeLines(c(header1, header2), file_conn)

  # Write column headers
  writeLines(paste(colnames(shared_genes), collapse = ","), file_conn)

  # Write data rows
  if (nrow(shared_genes) > 0) {
    for (i in 1:nrow(shared_genes)) {
      row_text <- paste(shared_genes[i, ], collapse = ",")
      writeLines(row_text, file_conn)
    }
  }

  close(file_conn)

  cat("  Saved:", basename(output_file), "-", nrow(shared_genes), "shared genes\n")

  return(nrow(shared_genes))
}

#' Process all clades for one metacluster
#'
#' @param clades Named list of clade definitions
#' @param gene_list_dir Gene list directory
#' @param output_dir Output directory
#' @return Data frame with clade summary statistics
process_metacluster_clades <- function(clades, gene_list_dir, output_dir) {

  summary_data <- data.frame(
    clade_name = character(),
    n_motifs = integer(),
    n_shared_genes = integer(),
    motif_list = character(),
    stringsAsFactors = FALSE
  )

  for (clade_name in names(clades)) {
    cat("\nProcessing clade:", clade_name, "\n")

    clade_motifs <- clades[[clade_name]]
    cat("  Motifs:", paste(clade_motifs, collapse = ", "), "\n")

    # Find shared genes
    shared_genes <- find_clade_shared_genes(clade_motifs, gene_list_dir)

    # Create CSV
    n_shared <- create_clade_csv(clade_name, clade_motifs, shared_genes, output_dir)

    # Add to summary
    summary_data <- rbind(summary_data, data.frame(
      clade_name = clade_name,
      n_motifs = length(clade_motifs),
      n_shared_genes = n_shared,
      motif_list = paste(clade_motifs, collapse = ","),
      stringsAsFactors = FALSE
    ))
  }

  return(summary_data)
}

# MAIN EXECUTION
cat("Starting clade shared gene extraction...\n")
cat(rep("=", 70), "\n", sep = "")

# Process p0m clades
cat("\n", rep("=", 70), "\n", sep = "")
cat("PROCESSING p0m CLADES\n")
cat(rep("=", 70), "\n", sep = "")

p0m_summary <- process_metacluster_clades(p0m_clades, GENE_LIST_DIR_P0M, OUTPUT_DIR)

# Process p1m clades
cat("\n", rep("=", 70), "\n", sep = "")
cat("PROCESSING p1m CLADES\n")
cat(rep("=", 70), "\n", sep = "")

p1m_summary <- process_metacluster_clades(p1m_clades, GENE_LIST_DIR_P1M, OUTPUT_DIR)

# Combine summaries
all_summary <- rbind(p0m_summary, p1m_summary)

# Save overall summary
summary_file <- file.path(OUTPUT_DIR, "clade_summary.csv")
write.csv(all_summary, summary_file, row.names = FALSE, quote = FALSE)

cat("\n", rep("=", 70), "\n", sep = "")
cat("CLADE ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Summary Statistics:\n")
cat("  Total clades processed:", nrow(all_summary), "\n")
cat("  p0m clades:", nrow(p0m_summary), "\n")
cat("  p1m clades:", nrow(p1m_summary), "\n\n")

cat("Shared genes per clade:\n")
for (i in 1:nrow(all_summary)) {
  cat("  ", all_summary$clade_name[i], ":",
      all_summary$n_shared_genes[i], "genes (",
      all_summary$n_motifs[i], "motifs)\n")
}

cat("\nOutput Files:\n")
cat("  Summary: clade_summary.csv\n")
cat("  Clade directories: ", nrow(all_summary), "subdirectories\n")
cat("  Output directory:", OUTPUT_DIR, "\n")

# Print session info
cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Script completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
