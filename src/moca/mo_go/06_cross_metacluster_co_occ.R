#!/usr/bin/env Rscript
# Cross-Metacluster Motif Co-occurrence Analysis
# Identifies genes with both P0 (activating) and P1 (repressing) motifs
# Uses Jaccard similarity and Pearson correlation with hierarchical clustering
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_go")
######################

library(dplyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dendextend)

# Source utility functions
source("../utils.R")

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Default paths
default_p0m_dir <- "../../../out/moca_results/mo_go/motif_gene_mapping/p0m_motifs/gene_lists"
default_p1m_dir <- "../../../out/moca_results/mo_go/motif_gene_mapping/p1m_motifs/gene_lists"
default_output_dir <- "../../../out/moca_results/mo_go/cross_metacluster_co_occ"

# Parse arguments
P0M_DIR <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_p0m_dir
P1M_DIR <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_p1m_dir
OUTPUT_DIR <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_output_dir

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

cat("Cross-Metacluster Co-occurrence Analysis Parameters:\n")
cat("  p0m gene lists:", P0M_DIR, "\n")
cat("  p1m gene lists:", P1M_DIR, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n\n")

#' Read gene list CSV files with comment lines
#'
#' @param file_path Path to gene list CSV
#' @return Data frame with gene_id, expression_class, chromosome
read_gene_list <- function(file_path) {

  # Read all lines
  all_lines <- readLines(file_path)

  # Find where data starts (after comment lines and header)
  data_start <- which(!grepl("^#", all_lines) & !grepl("^\"#", all_lines))

  # Find header line
  header_line <- data_start[1] - 1

  # Read the file starting from header
  gene_data <- read.csv(file_path, skip = header_line - 1, stringsAsFactors = FALSE)

  # Remove any remaining comment rows AND duplicate header rows
  gene_data <- gene_data %>%
    filter(!grepl("^#", gene_id)) %>%
    filter(gene_id != "gene_id")

  return(gene_data)
}

#' Load all motif gene lists from a directory
#'
#' @param gene_list_dir Directory containing *_genes.csv files
#' @param metacluster_prefix Prefix to add to motif names ("p0m" or "p1m")
#' @return List of gene sets, named by motif with prefix
load_motif_gene_sets <- function(gene_list_dir, metacluster_prefix) {

  cat("\nLoading gene lists from:", gene_list_dir, "\n")
  cat("  Metacluster prefix:", metacluster_prefix, "\n")

  # Find all gene list files
  gene_files <- list.files(gene_list_dir, pattern = "_genes\\.csv$", full.names = TRUE)

  if (length(gene_files) == 0) {
    stop("No gene list files found in ", gene_list_dir)
  }

  cat("  Found", length(gene_files), "gene list files\n")

  # Load each file
  gene_sets <- list()

  for (file in gene_files) {
    # Extract motif name from filename
    motif_name <- gsub("_genes\\.csv$", "", basename(file))

    # Add metacluster prefix for identification
    full_name <- paste0(metacluster_prefix, "_", motif_name)

    # Read gene list
    genes <- read_gene_list(file)

    # Store gene IDs
    gene_sets[[full_name]] <- genes$gene_id
  }

  cat("  Loaded gene sets for", length(gene_sets), "motifs\n")

  return(gene_sets)
}

#' Calculate Jaccard similarity between two gene sets
#'
#' @param set1 Vector of gene IDs
#' @param set2 Vector of gene IDs
#' @return Jaccard index (intersection/union)
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))

  if (union == 0) return(0)
  return(intersection / union)
}

#' Calculate cross-metacluster similarity matrix
#'
#' @param p0m_gene_sets List of P0 motif gene sets
#' @param p1m_gene_sets List of P1 motif gene sets
#' @return Matrix of Jaccard similarities (P0 motifs x P1 motifs)
calculate_cross_similarity_matrix <- function(p0m_gene_sets, p1m_gene_sets) {

  cat("\nCalculating cross-metacluster Jaccard similarity matrix...\n")

  p0m_motifs <- names(p0m_gene_sets)
  p1m_motifs <- names(p1m_gene_sets)

  # Initialize matrix
  similarity_matrix <- matrix(0,
                              nrow = length(p0m_motifs),
                              ncol = length(p1m_motifs),
                              dimnames = list(p0m_motifs, p1m_motifs))

  # Calculate Jaccard similarity for each P0-P1 pair
  for (i in 1:length(p0m_motifs)) {
    for (j in 1:length(p1m_motifs)) {
      similarity_matrix[i, j] <- jaccard_similarity(
        p0m_gene_sets[[p0m_motifs[i]]],
        p1m_gene_sets[[p1m_motifs[j]]]
      )
    }
  }

  cat("  Matrix dimensions:", nrow(similarity_matrix), "P0 motifs x",
      ncol(similarity_matrix), "P1 motifs\n")
  cat("  Jaccard range:", round(min(similarity_matrix), 4), "-",
      round(max(similarity_matrix), 4), "\n")

  return(similarity_matrix)
}

#' Create hierarchical clustering heatmap with Pearson correlation
#'
#' @param similarity_matrix Jaccard similarity matrix (P0 x P1)
#' @param output_dir Output directory
create_cross_metacluster_heatmap <- function(similarity_matrix, output_dir) {

  cat("\nCreating cross-metacluster co-occurrence heatmap...\n")

  # Clean motif names for display (remove prefixes for readability)
  rownames(similarity_matrix) <- gsub("^p0m_", "P0:", rownames(similarity_matrix))
  colnames(similarity_matrix) <- gsub("^p1m_", "P1:", colnames(similarity_matrix))

  # Create output file path
  output_file <- file.path(output_dir, "cross_metacluster_heatmap.pdf")

  # Create heatmap with hierarchical clustering - VERY LARGE
  pdf(output_file, width = 40, height = 32)

  pheatmap(
    similarity_matrix,
    color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 18,
    fontsize_row = 16,
    fontsize_col = 16,
    main = "Cross-Metacluster Motif Co-occurrence\n(P0 Activating vs P1 Repressing Motifs)",
    border_color = "grey90",
    angle_col = 45,
    legend = TRUE,
    legend_labels = "Jaccard\nSimilarity",
    cellwidth = 30,
    cellheight = 30
  )

  dev.off()

  cat("  Saved heatmap:", basename(output_file), "\n")
}

#' Create correlation-based dendrogram
#'
#' @param similarity_matrix Jaccard similarity matrix (P0 x P1)
#' @param output_dir Output directory
create_correlation_dendrogram <- function(similarity_matrix, output_dir) {

  cat("\nCreating Pearson correlation dendrogram...\n")

  # Create a symmetric matrix for correlation analysis
  # For asymmetric matrix (P0 x P1), we need to create a full similarity matrix

  n_p0 <- nrow(similarity_matrix)
  n_p1 <- ncol(similarity_matrix)
  n_total <- n_p0 + n_p1

  # Create full symmetric matrix
  full_matrix <- matrix(0, nrow = n_total, ncol = n_total)
  rownames(full_matrix) <- c(rownames(similarity_matrix), colnames(similarity_matrix))
  colnames(full_matrix) <- c(rownames(similarity_matrix), colnames(similarity_matrix))

  # Fill in the cross-metacluster similarities
  # Upper right: P0 x P1
  full_matrix[1:n_p0, (n_p0+1):n_total] <- similarity_matrix
  # Lower left: P1 x P0 (transpose)
  full_matrix[(n_p0+1):n_total, 1:n_p0] <- t(similarity_matrix)

  # Diagonal is 1 (self-similarity)
  diag(full_matrix) <- 1

  # Calculate Pearson correlation on this combined matrix
  cor_matrix <- cor(full_matrix, method = "pearson")

  # Convert to distance (1 - correlation for dissimilarity)
  cor_dist <- as.dist(1 - abs(cor_matrix))

  # Hierarchical clustering
  hclust_result <- hclust(cor_dist, method = "ward.D2")

  # Convert to dendrogram
  dend <- as.dendrogram(hclust_result)

  # Define colors: P0 = pink, P1 = teal
  p0_color <- "#E91E63"  # Pink
  p1_color <- "#176272"  # Teal

  # Color labels by metacluster
  motif_labels <- labels(dend)
  labels_colors <- ifelse(grepl("^p0m_", motif_labels), p0_color, p1_color)

  # Apply label colors only - keep branches black except terminal ones
  dend <- dend %>%
    set("labels_col", labels_colors) %>%
    set("labels_cex", 0.9) %>%
    set("branches_lwd", 2)  # Thicker branches

  # Color only the terminal (leaf) branches
  # This function colors just the last segment leading to each leaf
  dend <- dend %>%
    dendextend::assign_values_to_leaves_edgePar(value = labels_colors, edgePar = "col") %>%
    dendextend::assign_values_to_leaves_edgePar(value = 3, edgePar = "lwd")  # Thicker terminal branches

  # Save dendrogram
  output_file <- file.path(output_dir, "cross_metacluster_dendrogram.pdf")
  pdf(output_file, width = 22, height = 12)

  par(mar = c(12, 5, 4, 2))
  plot(dend,
       main = "Cross-Metacluster Motif Relationships\n(Pearson Correlation Clustering)",
       ylab = "Height (1 - |Pearson r|)",
       cex.main = 1.4,
       cex.axis = 1.1,
       cex.lab = 1.2)

  legend("topright",
         legend = c("P0 motifs (activating)", "P1 motifs (repressing)"),
         col = c(p0_color, p1_color),
         lwd = 4,
         cex = 1.2,
         bty = "n")

  dev.off()

  cat("  Saved dendrogram:", basename(output_file), "\n")

  return(list(hclust = hclust_result, cor_matrix = cor_matrix, full_matrix = full_matrix))
}

#' Identify significant cross-metacluster co-occurrences
#'
#' @param similarity_matrix Jaccard similarity matrix
#' @param p0m_gene_sets P0 motif gene sets
#' @param p1m_gene_sets P1 motif gene sets
#' @param output_dir Output directory
#' @param threshold Minimum Jaccard similarity threshold
identify_significant_pairs <- function(similarity_matrix, p0m_gene_sets,
                                      p1m_gene_sets, output_dir,
                                      threshold = 0.05) {

  cat("\nIdentifying significant cross-metacluster pairs...\n")
  cat("  Jaccard threshold:", threshold, "\n")

  # Create data frame of all pairs with their similarities
  pairs_df <- data.frame()

  for (p0_motif in rownames(similarity_matrix)) {
    for (p1_motif in colnames(similarity_matrix)) {
      jaccard <- similarity_matrix[p0_motif, p1_motif]

      if (jaccard >= threshold) {
        # Get original names (with prefixes)
        p0_original <- paste0("p0m_", gsub("^P0:", "", p0_motif))
        p1_original <- paste0("p1m_", gsub("^P1:", "", p1_motif))

        # Get shared genes
        shared_genes <- intersect(p0m_gene_sets[[p0_original]],
                                 p1m_gene_sets[[p1_original]])

        pairs_df <- rbind(pairs_df, data.frame(
          p0_motif = p0_motif,
          p1_motif = p1_motif,
          jaccard_similarity = round(jaccard, 4),
          p0_gene_count = length(p0m_gene_sets[[p0_original]]),
          p1_gene_count = length(p1m_gene_sets[[p1_original]]),
          shared_gene_count = length(shared_genes),
          shared_genes = paste(shared_genes, collapse = ","),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Sort by Jaccard similarity (descending)
  pairs_df <- pairs_df %>% arrange(desc(jaccard_similarity))

  cat("  Found", nrow(pairs_df), "significant pairs (Jaccard >=", threshold, ")\n")

  if (nrow(pairs_df) > 0) {
    cat("  Jaccard range:", min(pairs_df$jaccard_similarity), "-",
        max(pairs_df$jaccard_similarity), "\n")
    cat("  Shared gene range:", min(pairs_df$shared_gene_count), "-",
        max(pairs_df$shared_gene_count), "\n")
  }

  # Save results
  output_file <- file.path(output_dir, "significant_cross_metacluster_pairs.csv")
  write.csv(pairs_df, output_file, row.names = FALSE)

  cat("  Saved significant pairs:", basename(output_file), "\n")

  return(pairs_df)
}

#' Identify genes with mixed regulatory signals
#'
#' @param p0m_gene_sets P0 motif gene sets
#' @param p1m_gene_sets P1 motif gene sets
#' @param output_dir Output directory
identify_mixed_regulation_genes <- function(p0m_gene_sets, p1m_gene_sets, output_dir) {

  cat("\nIdentifying genes with mixed regulatory signals...\n")

  # Get all unique genes from each metacluster
  all_p0_genes <- unique(unlist(p0m_gene_sets))
  all_p1_genes <- unique(unlist(p1m_gene_sets))

  # Find genes present in both metaclusters
  mixed_genes <- intersect(all_p0_genes, all_p1_genes)

  cat("  Genes with P0 motifs only:", length(setdiff(all_p0_genes, all_p1_genes)), "\n")
  cat("  Genes with P1 motifs only:", length(setdiff(all_p1_genes, all_p0_genes)), "\n")
  cat("  Genes with BOTH P0 and P1 motifs:", length(mixed_genes), "\n")

  if (length(mixed_genes) > 0) {
    # For each mixed gene, count P0 and P1 motifs
    mixed_df <- data.frame()

    for (gene in mixed_genes) {
      # Count P0 motifs
      p0_motifs <- names(p0m_gene_sets)[sapply(p0m_gene_sets, function(x) gene %in% x)]
      p0_count <- length(p0_motifs)

      # Count P1 motifs
      p1_motifs <- names(p1m_gene_sets)[sapply(p1m_gene_sets, function(x) gene %in% x)]
      p1_count <- length(p1_motifs)

      mixed_df <- rbind(mixed_df, data.frame(
        gene_id = gene,
        p0_motif_count = p0_count,
        p1_motif_count = p1_count,
        total_motif_count = p0_count + p1_count,
        p0_motifs = paste(gsub("^p0m_", "", p0_motifs), collapse = ","),
        p1_motifs = paste(gsub("^p1m_", "", p1_motifs), collapse = ","),
        stringsAsFactors = FALSE
      ))
    }

    # Sort by total motif count
    mixed_df <- mixed_df %>% arrange(desc(total_motif_count))

    # Save results
    output_file <- file.path(output_dir, "genes_with_mixed_regulation.csv")
    write.csv(mixed_df, output_file, row.names = FALSE)

    cat("  Saved mixed regulation genes:", basename(output_file), "\n")
    cat("  Average P0 motifs per gene:", round(mean(mixed_df$p0_motif_count), 2), "\n")
    cat("  Average P1 motifs per gene:", round(mean(mixed_df$p1_motif_count), 2), "\n")

    # Create summary plot
    create_mixed_regulation_plot(mixed_df, output_dir)

    return(mixed_df)
  } else {
    cat("  No genes with mixed regulation found.\n")
    return(data.frame())
  }
}

#' Create plot showing distribution of mixed regulation
#'
#' @param mixed_df Data frame of genes with mixed regulation
#' @param output_dir Output directory
create_mixed_regulation_plot <- function(mixed_df, output_dir) {

  cat("\nCreating mixed regulation distribution plot...\n")

  # Create scatter plot with detailed annotations
  p <- ggplot(mixed_df, aes(x = p0_motif_count, y = p1_motif_count)) +
    geom_point(alpha = 0.7, size = 4, color = "#7D3C98") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
    labs(
      title = "Genes with Mixed Regulatory Signals",
      subtitle = "Each dot = one gene with both P0 (activating) and P1 (repressing) motifs",
      x = "Number of P0 Motifs (Activating)",
      y = "Number of P1 Motifs (Repressing)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95")
    ) +
    annotate("text", x = max(mixed_df$p0_motif_count) * 0.8,
             y = max(mixed_df$p1_motif_count) * 0.1,
             label = paste0("Total genes: ", nrow(mixed_df)),
             size = 5, color = "gray30", fontface = "bold") +
    annotate("text",
             x = max(mixed_df$p0_motif_count) * 0.05,
             y = max(mixed_df$p1_motif_count) * 0.95,
             label = "Diagonal (dashed line):\nEqual P0 and P1 motif counts",
             size = 4, color = "gray40", hjust = 0, vjust = 1) +
    annotate("text",
             x = max(mixed_df$p0_motif_count) * 0.95,
             y = max(mixed_df$p1_motif_count) * 0.25,
             label = "Below diagonal:\nMore P0 (activating) motifs",
             size = 4, color = "#E91E63", hjust = 1, fontface = "italic") +
    annotate("text",
             x = max(mixed_df$p0_motif_count) * 0.25,
             y = max(mixed_df$p1_motif_count) * 0.95,
             label = "Above diagonal:\nMore P1 (repressing) motifs",
             size = 4, color = "#176272", hjust = 0, vjust = 1, fontface = "italic")

  # Save plot with larger size
  output_file <- file.path(output_dir, "mixed_regulation_scatter.pdf")
  ggsave(output_file, p, width = 11, height = 9)

  cat("  Saved scatter plot:", basename(output_file), "\n")
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("CROSS-METACLUSTER CO-OCCURRENCE ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep = "")

# Load P0 and P1 gene sets
p0m_gene_sets <- load_motif_gene_sets(P0M_DIR, "p0m")
p1m_gene_sets <- load_motif_gene_sets(P1M_DIR, "p1m")

# Calculate cross-metacluster similarity matrix
similarity_matrix <- calculate_cross_similarity_matrix(p0m_gene_sets, p1m_gene_sets)

# Save similarity matrix
matrix_file <- file.path(OUTPUT_DIR, "jaccard_similarity_matrix.csv")
write.csv(similarity_matrix, matrix_file, row.names = TRUE)
cat("\nSaved similarity matrix:", basename(matrix_file), "\n")

# Create heatmap
create_cross_metacluster_heatmap(similarity_matrix, OUTPUT_DIR)

# Create correlation dendrogram
dend_results <- create_correlation_dendrogram(similarity_matrix, OUTPUT_DIR)

# Identify significant pairs
significant_pairs <- identify_significant_pairs(similarity_matrix,
                                                p0m_gene_sets,
                                                p1m_gene_sets,
                                                OUTPUT_DIR,
                                                threshold = 0.05)

# Identify genes with mixed regulation
mixed_genes <- identify_mixed_regulation_genes(p0m_gene_sets, p1m_gene_sets, OUTPUT_DIR)

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Output Files Generated:\n")
cat("  1. jaccard_similarity_matrix.csv\n")
cat("  2. cross_metacluster_heatmap.pdf\n")
cat("  3. cross_metacluster_dendrogram.pdf\n")
cat("  4. significant_cross_metacluster_pairs.csv\n")
cat("  5. genes_with_mixed_regulation.csv\n")
cat("  6. mixed_regulation_scatter.pdf\n")

cat("\nKey Statistics:\n")
cat("  P0 motifs analyzed:", length(p0m_gene_sets), "\n")
cat("  P1 motifs analyzed:", length(p1m_gene_sets), "\n")
cat("  Significant cross-pairs (Jaccard >= 0.05):", nrow(significant_pairs), "\n")
cat("  Genes with mixed regulation:", nrow(mixed_genes), "\n")

cat("\nBiological Interpretation:\n")
cat("  - Significant pairs: P0 and P1 motifs co-occurring in same genes\n")
cat("  - Mixed regulation: Genes with BOTH activating and repressing motifs\n")
cat("  - High Jaccard: Strong combinatorial regulation\n")
cat("  - Dendrogram: Shows which P0/P1 motifs co-occur most frequently\n")

cat("\nOutput Directory:", OUTPUT_DIR, "\n")

# Print session info
cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Script completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
