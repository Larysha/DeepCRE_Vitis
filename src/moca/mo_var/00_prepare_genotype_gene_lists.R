#!/usr/bin/env Rscript
# Preprocessing Script: Prepare Gene Lists for Genotype Variance Analysis
#
# PURPOSE:
# This script prepares gene lists for cross-genotype/variety motif variance analysis.
# It extracts differentially expressed (DE) and uniformly expressed (UE) gene sets
# from expression data and creates GFF annotation subsets for downstream BLAMM mapping.
#
# WORKFLOW OVERVIEW:
# 1. Load expression data for multiple genotypes/varieties
# 2. Identify differentially expressed genes (DE) across genotypes
# 3. Identify uniformly expressed genes (UE) across genotypes
# 4. Extract GFF annotations for DE and UE gene sets
# 5. Prepare coordinate files for FASTA extraction
#
# ORIGINAL PREPROCESSING STEPS (from mo_genotype_variance.v1.4.R comments):
# - Created list of differentially and uniformly expressed genes
# - Split gene_model.gff files with their annotations using:
#   while read id; do grep "$id" *_gene_models.gff >> Sol_genotypes_differentially_expressed_location2.txt; done < Sol_genotypes_differentially_expressed.txt
# - Created multifasta file using extract_range_to_fasta.sh
# - Merged all fastas into one file
# - Used BLAMM to find EPMs
# - Used occ_filter_v1.1.R to parse matches
#
# FUTURE INPUTS (when in silico perturbation data available):
# - Per-genotype expression predictions from CNN model
# - In silico perturbed sequences and predictions
# - Genotype-specific expression measurements
#
# REFERENCE:
# Methods from: https://www.nature.com/articles/s41467-024-47744-0
# (Section on genotype-specific motif conservation analysis)
#
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_var")
######################

library(tidyr)
library(dplyr)
library(readr)
library(stringr)

# Source utility functions
source("../utils.R")

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

######################
# FILE PATHS - PATTERN-BASED FINDING
######################

# Expression data for multiple genotypes/varieties
# FUTURE: Replace with actual multi-genotype expression data
# Expected format: gene_id, genotype, expression_level, class
default_expr_file <- "../../../vitis_data/tpm_counts/vitis_multi_genotype_expression.csv"  # PLACEHOLDER

# Gene annotation GFF
default_gff_file <- "../../../vitis_data/gene_models/vitis_vinifera_PN40024.gff3"

# Output directory
default_output_dir <- "../../../out/moca_results/mo_var/preprocessing"

# Parse arguments
EXPR_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_expr_file
GFF_FILE <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_gff_file
OUTPUT_DIR <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_output_dir

# Analysis parameters
DATE_STAMP <- format(Sys.Date(), "%Y%m%d")
PROJECT_NAME <- "vitis_genotype_variance_prep"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

cat("Genotype Variance Analysis - Preprocessing:\n")
cat("  Expression file:", EXPR_FILE, "\n")
cat("  GFF file:", GFF_FILE, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  Date:", DATE_STAMP, "\n\n")

######################
# SECTION 1: LOAD MULTI-GENOTYPE EXPRESSION DATA
######################
#
# DESCRIPTION:
# Load expression data for multiple Vitis varieties/genotypes
# This could come from:
# - Multi-variety RNA-seq experiments
# - In silico perturbation predictions from trained models
# - Cross-variety expression profiling studies
#
# EXPECTED DATA FORMAT (based on paper methodology):
# gene_id       | genotype    | cnn_probability | expression_class | variance_across_genotypes | differential_status
# Vitvi01g00010 | PN40024     | 0.62            | high            | 0.003                     | uniform
# Vitvi01g00010 | Cabernet    | 0.58            | high            | 0.003                     | uniform
# Vitvi01g00010 | Chardonnay  | 0.61            | high            | 0.003                     | uniform
# Vitvi01g00020 | PN40024     | 0.35            | low             | 0.012                     | differential
# Vitvi01g00020 | Cabernet    | 0.78            | high            | 0.012                     | differential
#
# WHERE:
# - cnn_probability: MSR model prediction (0-1 scale) from trained DeepCRE models
# - expression_class: "high" if probability >0.5, "low" if probability ≤0.5
# - variance_across_genotypes: variance of CNN probabilities across genotypes
# - differential_status: "differential" if variance >0.005, "uniform" otherwise
#
# NOTE: The paper uses CNN model predictions to classify genes, not direct RNA-seq data.
# This allows validation that predicted differential expression correlates with motif
# conservation patterns across genotypes.
#
######################

cat("\n=== SECTION 1: Load Multi-Genotype Expression Data ===\n")

# PSEUDOCODE - Replace with actual data loading when available
if (file.exists(EXPR_FILE)) {
  cat("  Loading multi-genotype expression data...\n")
  expression_data <- read.csv(EXPR_FILE, stringsAsFactors = FALSE)
} else {
  cat("  WARNING: Expression file not found - using PLACEHOLDER data structure\n")
  cat("  FUTURE: Replace with actual multi-genotype expression data\n")

  # PLACEHOLDER: Create example data structure for reference
  expression_data <- data.frame(
    gene_id = character(),
    genotype = character(),
    cnn_probability = numeric(),
    expression_class = character(),
    variance_across_genotypes = numeric(),
    differential_status = character(),
    stringsAsFactors = FALSE
  )

  cat("  Expected columns: gene_id, genotype, cnn_probability, expression_class,\n")
  cat("                    variance_across_genotypes, differential_status\n")
}

######################
# SECTION 2: IDENTIFY DIFFERENTIALLY EXPRESSED GENES
######################
#
# DESCRIPTION:
# Identify genes with differential expression across genotypes/varieties
#
# CRITERIA FOR DIFFERENTIAL EXPRESSION (based on paper methodology):
# - Variance in CNN probability predictions across genotypes >0.005
# - Expression class varies across genotypes (e.g., high in variety A, low in variety B)
# - Calculated from variance of MSR model predictions, not direct RNA-seq measurements
#
# BIOLOGICAL INTERPRETATION:
# DE genes may harbor genotype-specific regulatory elements
# Motif differences in these genes could explain expression divergence
# The paper uses this to validate that CNN predictions correlate with regulatory variation
#
######################

cat("\n=== SECTION 2: Identify Differentially Expressed Genes ===\n")

# PSEUDOCODE for identifying DE genes
# Based on paper methodology: variance >0.005 in CNN probability predictions
identify_differential_genes <- function(expr_data) {
  # Group by gene and calculate variance in CNN predictions across genotypes
  de_genes <- expr_data %>%
    group_by(gene_id) %>%
    summarise(
      n_genotypes = n_distinct(genotype),
      variance_cnn_prob = var(cnn_probability),
      n_classes = n_distinct(expression_class),
      is_differential = variance_cnn_prob > 0.005,  # Paper threshold
      .groups = "drop"
    ) %>%
    filter(is_differential == TRUE)

  return(de_genes$gene_id)
}

# WHEN DATA AVAILABLE:
# differential_genes <- identify_differential_genes(expression_data)
# cat("  Identified", length(differential_genes), "differentially expressed genes\n")

cat("  PLACEHOLDER: DE gene identification logic defined\n")
cat("  TODO: Run when multi-genotype expression data available\n")

######################
# SECTION 3: IDENTIFY UNIFORMLY EXPRESSED GENES
######################
#
# DESCRIPTION:
# Identify genes with uniform expression across genotypes/varieties
#
# CRITERIA FOR UNIFORM EXPRESSION (based on paper methodology):
# - Variance in CNN probability predictions across genotypes ≤0.005
# - Same expression class across all genotypes
# - Stable predicted expression regardless of genetic background
#
# BIOLOGICAL INTERPRETATION:
# UE genes serve as BACKGROUND/CONTROL set
# Expected to have conserved regulatory elements
# Used for comparison with DE genes to assess if motif variance
# is specifically enriched in DE genes
# Bootstrap resampling (100 genes, 1000 iterations) handles size imbalance
#
######################

cat("\n=== SECTION 3: Identify Uniformly Expressed Genes ===\n")

# PSEUDOCODE for identifying UE genes
# Based on paper methodology: variance ≤0.005 in CNN probability predictions
identify_uniform_genes <- function(expr_data) {
  # Group by gene and check variance in CNN predictions
  ue_genes <- expr_data %>%
    group_by(gene_id) %>%
    summarise(
      n_genotypes = n_distinct(genotype),
      variance_cnn_prob = var(cnn_probability),
      n_classes = n_distinct(expression_class),
      is_uniform = variance_cnn_prob <= 0.005,  # Paper threshold
      .groups = "drop"
    ) %>%
    filter(is_uniform == TRUE)

  return(ue_genes$gene_id)
}

# WHEN DATA AVAILABLE:
# uniform_genes <- identify_uniform_genes(expression_data)
# cat("  Identified", length(uniform_genes), "uniformly expressed genes\n")

cat("  PLACEHOLDER: UE gene identification logic defined\n")
cat("  TODO: Run when multi-genotype expression data available\n")

######################
# SECTION 4: EXTRACT GFF ANNOTATIONS FOR GENE SETS
######################
#
# DESCRIPTION:
# Extract genomic coordinates and annotations from GFF for DE and UE genes
# This replicates the original bash loop preprocessing:
#   while read id; do grep "$id" *_gene_models.gff >> output.txt; done < gene_list.txt
#
# OUTPUT FORMAT (following original):
# genotype:chr:feature:start-end:strand:attributes
# Used for downstream FASTA extraction
#
######################

cat("\n=== SECTION 4: Extract GFF Annotations for Gene Sets ===\n")

# Load GFF annotations
load_gff_annotations <- function(gff_file) {
  cat("  Loading GFF file:", gff_file, "\n")

  if (!file.exists(gff_file)) {
    cat("  WARNING: GFF file not found\n")
    return(NULL)
  }

  # Read GFF (skip comment lines starting with #)
  gff_data <- read.table(gff_file,
                         header = FALSE,
                         sep = "\t",
                         comment.char = "#",
                         quote = "",
                         stringsAsFactors = FALSE)

  colnames(gff_data) <- c("seqname", "source", "feature", "start", "end",
                          "score", "strand", "frame", "attributes")

  # Filter for gene features only
  gff_data <- gff_data %>%
    filter(feature == "gene")

  return(gff_data)
}

# Extract gene coordinates with flanking regions (±1000bp as in original)
extract_gene_coordinates <- function(gff_data, gene_list, flank_size = 1000, genotype_prefix = "") {
  # Filter GFF for genes in list
  gene_coords <- gff_data %>%
    filter(grepl(paste(gene_list, collapse = "|"), attributes)) %>%
    mutate(
      # Extract gene ID from attributes
      gene_id = str_extract(attributes, "ID=([^;]+)"),
      gene_id = str_remove(gene_id, "ID="),

      # Add flanking regions
      region_start = pmax(1, start - flank_size),  # Don't go below 1
      region_end = end + flank_size,

      # Create location string (format: genotype_chr:start-end)
      loc = paste0(
        ifelse(genotype_prefix != "", paste0(genotype_prefix, "_"), ""),
        seqname, ":", region_start, "-", region_end
      )
    ) %>%
    select(loc, gene_id, seqname, region_start, region_end, strand, attributes)

  return(gene_coords)
}

# WHEN DATA AVAILABLE:
# gff_data <- load_gff_annotations(GFF_FILE)
# de_coords <- extract_gene_coordinates(gff_data, differential_genes, genotype_prefix = "")
# ue_coords <- extract_gene_coordinates(gff_data, uniform_genes, genotype_prefix = "")

cat("  PLACEHOLDER: GFF extraction functions defined\n")
cat("  TODO: Run when gene lists available\n")

######################
# SECTION 5: PREPARE COORDINATE FILES FOR FASTA EXTRACTION
######################
#
# DESCRIPTION:
# Create coordinate files for extract_range_to_fasta.sh script
# Output format compatible with bedtools/genomic range tools
#
# NEXT STEPS (manual):
# 1. Run extract_range_to_fasta_genomic.sh with coordinate files
# 2. Merge FASTAs from multiple genotypes
# 3. Run BLAMM motif mapping on merged FASTA
# 4. Use filter_blamm_occ_00.R to parse BLAMM output
# 5. Proceed to genotype_variance_analysis.R
#
######################

cat("\n=== SECTION 5: Prepare Coordinate Files ===\n")

# Write coordinate files
write_coordinate_files <- function(de_coords, ue_coords, output_dir, date_stamp) {
  # Differentially expressed genes
  de_file <- file.path(output_dir, paste0(date_stamp, "_vitis_genotypes_differential_expressed_locations.txt"))

  # Format: genotype_chr:start-end [tab] gene_id [tab] strand
  de_output <- de_coords %>%
    mutate(output_line = paste(seqname, region_start, region_end, gene_id, strand, sep = "\t"))

  # PSEUDOCODE - write when data available
  # write.table(de_output$output_line, de_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Uniformly expressed genes
  ue_file <- file.path(output_dir, paste0(date_stamp, "_vitis_genotypes_uniformly_expressed_locations.txt"))

  ue_output <- ue_coords %>%
    mutate(output_line = paste(seqname, region_start, region_end, gene_id, strand, sep = "\t"))

  # PSEUDOCODE - write when data available
  # write.table(ue_output$output_line, ue_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

  cat("  Files to be created:\n")
  cat("    DE coordinates:", de_file, "\n")
  cat("    UE coordinates:", ue_file, "\n")

  return(list(de_file = de_file, ue_file = ue_file))
}

# WHEN DATA AVAILABLE:
# coord_files <- write_coordinate_files(de_coords, ue_coords, OUTPUT_DIR, DATE_STAMP)

cat("  PLACEHOLDER: Coordinate file writing logic defined\n")

######################
# SECTION 6: SUMMARY AND NEXT STEPS
######################

cat("\n=== PREPROCESSING COMPLETE ===\n")
cat("\nWORKFLOW OVERVIEW (aligned with paper methodology):\n")
cat("\nDATA PREPARATION PHASE:\n")
cat("1. Train DeepCRE MSR models on multiple Vitis varieties\n")
cat("2. Generate CNN probability predictions for each gene across varieties\n")
cat("3. Calculate variance in predictions (threshold: >0.005 for differential)\n")
cat("4. Classify genes as differential (DE) or uniform (UE) expression\n")
cat("\nMOTIF MAPPING PHASE:\n")
cat("5. Run this preprocessing script with CNN prediction data\n")
cat("6. Extract FASTA sequences using coordinate files:\n")
cat("   bash ../ref_seq/extract_range_to_fasta_genomic.sh <genome> <coord_file>\n")
cat("7. Merge FASTAs from multiple genotypes into single file\n")
cat("8. Run BLAMM motif mapping on merged FASTA:\n")
cat("   blamm -e 0.0001 -w 14 <motif_pwms> <merged_fasta> > occurrences.txt\n")
cat("9. Filter BLAMM occurrences:\n")
cat("   Rscript ../mo_proj/filter_blamm_occ_00.R <occurrences.txt>\n")
cat("10. Filter to EPM preferred positional ranges (from mo_range outputs)\n")
cat("\nANALYSIS PHASE:\n")
cat("11. Run genotype variance analysis:\n")
cat("    Rscript 01_genotype_variance_analysis.R\n")
cat("12. Statistical tests: Fisher's exact, Chi-squared (p <0.0001)\n")
cat("13. Bootstrap resampling: 100 genes, 1000 iterations\n")
cat("\nKEY PAPER INSIGHT:\n")
cat("This validates that CNN predictions correlate with regulatory variation -\n")
cat("genes with predicted differential expression should show enrichment for\n")
cat("mutated motifs, while uniformly expressed genes should have conserved motifs.\n")

cat("\nOUTPUT DIRECTORY:", OUTPUT_DIR, "\n")
cat("Session completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")