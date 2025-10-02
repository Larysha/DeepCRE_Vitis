#!/usr/bin/env Rscript
# Generate per-motif BED files and gene lists from enrichment analysis
# Groups F and R strands together, separates metaclusters
# Creates IGV-compatible BED files with expression-based coloring
# Generates positional preference matrices and hierarchical clustering dendrograms
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_go")
######################

library(tidyr)
library(dplyr)
library(readr)
library(ape)

# Source utility functions
source("../utils.R")

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Default file paths
default_input_file <- "../../../out/moca_results/mo_proj/enrichment/20250930_vitis_motif_enrichment_integrated_data.csv"
default_output_dir <- "../../../out/moca_results/mo_go/motif_gene_mapping"

# Parse arguments
INPUT_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_input_file
OUTPUT_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_output_dir

cat("Motif-Gene Mapping Parameters:\n")
cat("  Input file:", INPUT_FILE, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n\n")

# Create output directory structure
p0m_bed_dir <- file.path(OUTPUT_DIR, "p0m_motifs", "bed")
p1m_bed_dir <- file.path(OUTPUT_DIR, "p1m_motifs", "bed")
p0m_gene_dir <- file.path(OUTPUT_DIR, "p0m_motifs", "gene_lists")
p1m_gene_dir <- file.path(OUTPUT_DIR, "p1m_motifs", "gene_lists")
p0m_matrix_dir <- file.path(OUTPUT_DIR, "p0m_motifs", "matrices")
p1m_matrix_dir <- file.path(OUTPUT_DIR, "p1m_motifs", "matrices")
p0m_dendro_dir <- file.path(OUTPUT_DIR, "p0m_motifs", "dendrograms")
p1m_dendro_dir <- file.path(OUTPUT_DIR, "p1m_motifs", "dendrograms")

for (dir in c(p0m_bed_dir, p1m_bed_dir, p0m_gene_dir, p1m_gene_dir,
              p0m_matrix_dir, p1m_matrix_dir, p0m_dendro_dir, p1m_dendro_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat("Created directory:", dir, "\n")
  }
}

#' Load enrichment integrated data
#'
#' @param input_file Path to integrated enrichment CSV
#' @return Data frame with motif-gene associations
load_integrated_data <- function(input_file) {

  cat("\nLoading integrated enrichment data...\n")

  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  data <- read.csv(input_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(data), "motif-gene associations\n")
  cat("  Unique genes:", length(unique(data$gene_id)), "\n")
  cat("  Unique motifs:", length(unique(data$motif)), "\n")

  return(data)
}

#' Create strand-agnostic motif identifier
#'
#' Removes F/R suffix to group forward and reverse strands together
#' @param motif_name Motif name (e.g., epm_vitis_ssr_p0m08F)
#' @return Strand-agnostic name (e.g., epm_vitis_ssr_p0m08)
create_motif_base_name <- function(motif_name) {
  # Remove F or R suffix from short_name format
  gsub("[FR]$", "", motif_name)
}

#' Determine metacluster from motif name
#'
#' @param motif_name Motif name
#' @return "p0m" or "p1m"
get_metacluster <- function(motif_name) {
  ifelse(grepl("p1m|_1_", motif_name), "p1m", "p0m")
}

#' Process motif data and prepare for output
#'
#' @param data Integrated enrichment data
#' @return List of data frames grouped by motif base name and metacluster
process_motif_data <- function(data) {

  cat("\nProcessing motif data for grouping...\n")

  # Use short_name if available (from dual naming system), otherwise use motif
  data$motif_to_use <- if ("short_name" %in% colnames(data) &&
                           sum(!is.na(data$short_name)) > 0) {
    ifelse(is.na(data$short_name), data$motif, data$short_name)
  } else {
    data$motif
  }

  # Create base motif name (without F/R strand suffix)
  data$motif_base <- create_motif_base_name(data$motif_to_use)

  # Determine metacluster
  data$metacluster <- get_metacluster(data$motif_to_use)

  # Add genomic coordinates if not present
  if (!"motif_genomic_start" %in% colnames(data)) {
    data$motif_genomic_start <- data$gene_start + data$pos_relative_to_gene - 1
  }

  if (!"motif_genomic_end" %in% colnames(data)) {
    data$motif_length <- data$mend - data$mstart + 1
    data$motif_genomic_end <- data$motif_genomic_start + data$motif_length - 1
  }

  # Ensure we have chr column
  if (!"chr" %in% colnames(data)) {
    stop("Missing 'chr' column in input data")
  }

  cat("  Unique base motifs (strand-agnostic):",
      length(unique(data$motif_base)), "\n")
  cat("  Metacluster distribution:\n")
  cat("    p0m motifs:", length(unique(data$motif_base[data$metacluster == "p0m"])), "\n")
  cat("    p1m motifs:", length(unique(data$motif_base[data$metacluster == "p1m"])), "\n")

  return(data)
}

#' Create BED file for a single motif
#'
#' @param motif_data Data frame with all occurrences of one motif
#' @param motif_base_name Base motif name (without strand)
#' @param output_file Path for output BED file
create_motif_bed <- function(motif_data, motif_base_name, output_file) {

  # Prepare BED format data
  # Color by expression: pink (255,105,180) for high, blue (65,105,225) for low
  bed_data <- motif_data %>%
    mutate(
      # BED format columns
      chrom = chr,
      chromStart = motif_genomic_start,
      chromEnd = motif_genomic_end,
      name = paste0(motif_base_name, "_", gene_id),  # Unique identifier
      score = score,
      # strand already exists
      thickStart = motif_genomic_start,
      thickEnd = motif_genomic_end,
      # Color based on expression: pink for high, blue for low
      itemRgb = ifelse(expr_binary == "high", "255,105,180", "65,105,225")
    ) %>%
    select(chrom, chromStart, chromEnd, name, score, strand,
           thickStart, thickEnd, itemRgb) %>%
    arrange(chrom, chromStart)

  # Write BED file (no header for BED format)
  write.table(bed_data, output_file, sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)

  return(nrow(bed_data))
}

#' Create gene list CSV for a single motif
#'
#' @param motif_data Data frame with all occurrences of one motif
#' @param motif_base_name Base motif name (without strand)
#' @param output_file Path for output CSV file
create_motif_gene_list <- function(motif_data, motif_base_name, output_file) {

  # Get unique genes with their expression class
  gene_list <- motif_data %>%
    select(gene_id, expr_binary, chr) %>%
    distinct() %>%
    arrange(chr, gene_id)

  # Add summary header information
  summary_info <- data.frame(
    gene_id = c(
      paste0("# Motif: ", motif_base_name),
      paste0("# Total unique genes: ", nrow(gene_list)),
      paste0("# High expression genes: ", sum(gene_list$expr_binary == "high")),
      paste0("# Low expression genes: ", sum(gene_list$expr_binary == "low")),
      "#",
      "gene_id"
    ),
    expr_binary = c(rep("", 5), "expression_class"),
    chr = c(rep("", 5), "chromosome")
  )

  # Combine summary and gene list
  output_data <- rbind(summary_info, gene_list)

  # Write CSV
  write.csv(output_data, output_file, row.names = FALSE, quote = FALSE)

  return(nrow(gene_list))
}

#' Generate all motif files (BED + gene lists)
#'
#' @param data Processed motif data
#' @param output_dirs Named list of output directories
generate_motif_files <- function(data, output_dirs) {

  cat("\nGenerating per-motif BED files and gene lists...\n")

  # Get unique combinations of base motif and metacluster
  motif_groups <- data %>%
    select(motif_base, metacluster) %>%
    distinct() %>%
    arrange(metacluster, motif_base)

  cat("  Processing", nrow(motif_groups), "unique motif-metacluster combinations\n\n")

  # Track statistics
  stats <- list(
    p0m_count = 0,
    p1m_count = 0,
    total_occurrences = 0,
    total_genes = 0
  )

  # Process each motif
  for (i in 1:nrow(motif_groups)) {
    motif_base <- motif_groups$motif_base[i]
    metacluster <- motif_groups$metacluster[i]

    # Get all data for this motif (both F and R strands)
    motif_data <- data %>%
      filter(motif_base == !!motif_base, metacluster == !!metacluster)

    # Determine output directories based on metacluster
    if (metacluster == "p0m") {
      bed_dir <- output_dirs$p0m_bed
      gene_dir <- output_dirs$p0m_gene
      stats$p0m_count <- stats$p0m_count + 1
    } else {
      bed_dir <- output_dirs$p1m_bed
      gene_dir <- output_dirs$p1m_gene
      stats$p1m_count <- stats$p1m_count + 1
    }

    # Create output file paths
    bed_file <- file.path(bed_dir, paste0(motif_base, ".bed"))
    gene_file <- file.path(gene_dir, paste0(motif_base, "_genes.csv"))

    # Generate files
    n_occurrences <- create_motif_bed(motif_data, motif_base, bed_file)
    n_genes <- create_motif_gene_list(motif_data, motif_base, gene_file)

    stats$total_occurrences <- stats$total_occurrences + n_occurrences
    stats$total_genes <- stats$total_genes + n_genes

    # Progress update
    if (i %% 10 == 0 || i == nrow(motif_groups)) {
      cat("  Processed", i, "/", nrow(motif_groups), "motifs\n")
    }
  }

  return(stats)
}

#' Create positional distance matrices for motifs
#'
#' Based on original moca_blue logic - creates matrices showing distance
#' to transcription borders for each motif-gene combination
#'
#' @param data Processed motif data with positional information
#' @param metacluster "p0m" or "p1m"
#' @param output_dir Directory for matrix outputs
create_distance_matrices <- function(data, metacluster, output_dir) {

  cat("\nCreating distance matrices for", metacluster, "motifs...\n")

  # Filter for this metacluster
  mc_data <- data %>% filter(metacluster == !!metacluster)

  if (nrow(mc_data) == 0) {
    cat("  No data for", metacluster, "- skipping\n")
    return(NULL)
  }

  # Create motif identifier with strand and region info (like original)
  mc_data <- mc_data %>%
    mutate(
      motif_id = paste0(motif_base, strand,
                       ifelse(classic_region_corrected == "upstream", "up", "do"))
    )

  # Create wide format: genes x motifs with distance values
  # Each cell contains distance to transcription border
  # Handle multiple occurrences by taking the minimum distance (closest to border)
  distance_matrix <- mc_data %>%
    select(gene_id, motif_id, dist_transc_border) %>%
    group_by(gene_id, motif_id) %>%
    summarise(dist_transc_border = min(dist_transc_border, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = motif_id, values_from = dist_transc_border)

  # Export distance matrix
  output_file <- file.path(output_dir,
                          paste0(metacluster, "_distance_matrix.csv"))
  write.csv(distance_matrix, output_file, row.names = FALSE)

  cat("  Created distance matrix:", basename(output_file), "\n")
  cat("    Dimensions:", nrow(distance_matrix), "genes x",
      ncol(distance_matrix) - 1, "motif patterns\n")

  return(distance_matrix)
}

#' Create binary presence/absence matrices for motifs
#'
#' Creates contingency tables showing motif presence (1) or absence (0) per gene
#'
#' @param data Processed motif data
#' @param metacluster "p0m" or "p1m"
#' @param output_dir Directory for matrix outputs
create_binary_matrices <- function(data, metacluster, output_dir) {

  cat("\nCreating binary matrices for", metacluster, "motifs...\n")

  # Filter for this metacluster
  mc_data <- data %>% filter(metacluster == !!metacluster)

  if (nrow(mc_data) == 0) {
    cat("  No data for", metacluster, "- skipping\n")
    return(NULL)
  }

  # Create contingency table: gene_id x motif_base
  contingency <- table(mc_data$gene_id, mc_data$motif_base)

  # Convert to data frame
  binary_df <- as.data.frame(contingency) %>%
    pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0)

  # Export binary matrix (numeric)
  binary_output <- file.path(output_dir,
                             paste0(metacluster, "_binary_matrix.csv"))
  write.csv(binary_df, binary_output, row.names = FALSE)

  # Create gene list version (replace counts with gene IDs like original)
  gene_list_df <- binary_df
  gene_list_df[, -1] <- lapply(gene_list_df[, -1], as.character)

  # Replace non-zero values with gene ID
  for (i in 2:ncol(gene_list_df)) {
    gene_list_df[gene_list_df[, i] != "0", i] <-
      gene_list_df$Var1[gene_list_df[, i] != "0"]
  }

  gene_list_output <- file.path(output_dir,
                                paste0(metacluster, "_gene_list_matrix.csv"))
  write.csv(gene_list_df, gene_list_output, row.names = FALSE)

  cat("  Created binary matrix:", basename(binary_output), "\n")
  cat("  Created gene list matrix:", basename(gene_list_output), "\n")
  cat("    Dimensions:", nrow(binary_df), "genes x",
      ncol(binary_df) - 1, "motifs\n")

  return(binary_df)
}

#' Generate hierarchical clustering dendrogram from binary matrix
#'
#' Creates dendrogram showing motif co-occurrence patterns
#'
#' @param binary_matrix Binary matrix (genes x motifs)
#' @param metacluster "p0m" or "p1m"
#' @param output_dir Directory for dendrogram outputs
create_dendrogram <- function(binary_matrix, metacluster, output_dir) {

  cat("\nGenerating dendrogram for", metacluster, "motifs...\n")

  if (is.null(binary_matrix) || nrow(binary_matrix) < 2) {
    cat("  Insufficient data for clustering - skipping\n")
    return(NULL)
  }

  # Convert to matrix (exclude gene_id column)
  mat <- as.matrix(binary_matrix[, -1])

  # Check if we have enough motifs for clustering
  if (ncol(mat) < 2) {
    cat("  Only", ncol(mat), "motif - cannot cluster\n")
    return(NULL)
  }

  # Compute correlation matrix
  corr <- cor(mat)

  # Compute dissimilarity matrix (1 - correlation)
  diss <- 1 - corr

  # Hierarchical clustering
  hc <- hclust(as.dist(diss))

  # Save as PDF
  pdf_file <- file.path(output_dir,
                       paste0(metacluster, "_motif_dendrogram.pdf"))
  pdf(pdf_file)
  plot(hc, hang = -1, main = paste(metacluster, "Motif Co-occurrence Dendrogram"),
       xlab = "Motifs", ylab = "Dissimilarity (1 - correlation)")
  dev.off()

  # Save as Newick tree format (for phylogenetic tools)
  phy <- as.phylo(hc)
  newick_file <- file.path(output_dir,
                          paste0(metacluster, "_motif_dendrogram.nwk"))
  write.tree(phy, file = newick_file)

  cat("  Created dendrogram PDF:", basename(pdf_file), "\n")
  cat("  Created Newick tree:", basename(newick_file), "\n")
  cat("    Clustered", ncol(mat), "motifs based on co-occurrence\n")

  return(hc)
}

# MAIN EXECUTION
cat("Starting motif-gene mapping file generation...\n")
cat(rep("=", 60), "\n", sep = "")

# Load data
integrated_data <- load_integrated_data(INPUT_FILE)

# Process and prepare data
processed_data <- process_motif_data(integrated_data)

# Set up output directories
output_dirs <- list(
  p0m_bed = p0m_bed_dir,
  p1m_bed = p1m_bed_dir,
  p0m_gene = p0m_gene_dir,
  p1m_gene = p1m_gene_dir
)

# Generate all motif files (BED + gene lists)
stats <- generate_motif_files(processed_data, output_dirs)

# Generate positional distance matrices for each metacluster
p0m_dist_matrix <- create_distance_matrices(processed_data, "p0m", p0m_matrix_dir)
p1m_dist_matrix <- create_distance_matrices(processed_data, "p1m", p1m_matrix_dir)

# Generate binary presence/absence matrices
p0m_binary_matrix <- create_binary_matrices(processed_data, "p0m", p0m_matrix_dir)
p1m_binary_matrix <- create_binary_matrices(processed_data, "p1m", p1m_matrix_dir)

# Generate hierarchical clustering dendrograms
p0m_dendrogram <- create_dendrogram(p0m_binary_matrix, "p0m", p0m_dendro_dir)
p1m_dendrogram <- create_dendrogram(p1m_binary_matrix, "p1m", p1m_dendro_dir)

# Summary report
cat("\n", rep("=", 60), "\n", sep = "")
cat("MOTIF-GENE MAPPING COMPLETED\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Output Summary:\n")
cat("  p0m motifs:", stats$p0m_count, "\n")
cat("  p1m motifs:", stats$p1m_count, "\n")
cat("  Total motif files:", stats$p0m_count + stats$p1m_count, "\n")
cat("  Total occurrences processed:", stats$total_occurrences, "\n")
cat("  Total unique genes:", stats$total_genes, "\n")

cat("\nOutput Structure:\n")
cat("  ", OUTPUT_DIR, "/\n")
cat("    p0m_motifs/\n")
cat("      bed/          - BED files for IGV (", stats$p0m_count, " files)\n", sep = "")
cat("      gene_lists/   - Gene CSV files (", stats$p0m_count, " files)\n", sep = "")
cat("      matrices/     - Distance & binary matrices\n")
cat("      dendrograms/  - Hierarchical clustering (PDF + Newick)\n")
cat("    p1m_motifs/\n")
cat("      bed/          - BED files for IGV (", stats$p1m_count, " files)\n", sep = "")
cat("      gene_lists/   - Gene CSV files (", stats$p1m_count, " files)\n", sep = "")
cat("      matrices/     - Distance & binary matrices\n")
cat("      dendrograms/  - Hierarchical clustering (PDF + Newick)\n")

cat("\nBED File Format:\n")
cat("  - Pink (255,105,180): High expression genes\n")
cat("  - Blue (65,105,225): Low expression genes\n")
cat("  - Load in IGV to visualize motif binding sites\n")

cat("\nGene List Format:\n")
cat("  - CSV files with unique genes per motif\n")
cat("  - Includes expression class and chromosome\n")
cat("  - Header comments with motif statistics\n")

cat("\nMatrix Outputs:\n")
cat("  - Distance matrices: Motif positional preferences (distance to TSS/TTS)\n")
cat("  - Binary matrices: Motif presence/absence per gene\n")
cat("  - Gene list matrices: Binary matrices with gene IDs instead of counts\n")

cat("\nDendrogram Outputs:\n")
cat("  - PDF: Visual hierarchical clustering of motif co-occurrence\n")
cat("  - Newick: Tree format for phylogenetic analysis tools\n")

cat("\nNext Steps:\n")
cat("  1. Load BED files in IGV to visualize motif locations\n")
cat("  2. Analyze dendrograms to identify co-occurring motif groups\n")
cat("  3. Use gene lists for GO enrichment analysis\n")
cat("  4. Explore distance matrices for positional preference patterns\n")

# Print session info
cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Script completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
