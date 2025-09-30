#!/usr/bin/env Rscript
# Genotype/Variety Variance Analysis: Motif Conservation and Mutation Patterns
#
# REFACTORED FROM: mo_genotype_variance.v1.4.R (original moca_blue project)
#
# PURPOSE:
# This script analyzes motif conservation and mutation patterns across multiple
# genotypes/varieties to understand the relationship between regulatory element
# variation and differential gene expression.
#
# KEY BIOLOGICAL QUESTIONS:
# 1. Are motifs more variable (mutated) in differentially expressed genes?
# 2. Are motifs more conserved in uniformly expressed genes?
# 3. Do specific motif types (p0m/p1m) show differential conservation patterns?
# 4. Can regulatory element variation explain expression divergence?
#
# ANALYSIS WORKFLOW:
# 1. Load motif occurrences mapped across multiple genotypes (from BLAMM)
# 2. Integrate with gene expression classifications (DE vs UE)
# 3. Count motif occurrences per gene per genotype
# 4. Classify motifs as "conserved" (present in all genotypes) or "mutated" (variable)
# 5. Statistical comparison of conservation patterns between DE and UE genes
# 6. Bootstrap resampling for robust statistical inference
# 7. Visualization of results
#
# REFERENCE:
# Methods adapted from: https://www.nature.com/articles/s41467-024-47744-0
# Original implementation: SMZ 2023-08-30
#
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_var")
######################

library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

# Source utility functions
source("../utils.R")

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

######################
# FILE PATHS - PATTERN-BASED FINDING
######################

# Motif occurrence data (filtered BLAMM output from preprocessing)
# Expected: Output from filter_blamm_occ_00.R run on multi-genotype merged FASTA
motif_pattern <- "../../../out/moca_results/mo_var/*_vitis_genotype_motif_occurrences_filtered.txt"
motif_files <- Sys.glob(motif_pattern)
default_motif_file <- if (length(motif_files) > 0) sort(motif_files, decreasing = TRUE)[1] else motif_pattern

# Differentially expressed gene annotations (from preprocessing)
de_pattern <- "../../../out/moca_results/mo_var/preprocessing/*_differential_expressed_locations.txt"
de_files <- Sys.glob(de_pattern)
default_de_file <- if (length(de_files) > 0) sort(de_files, decreasing = TRUE)[1] else de_pattern

# Uniformly expressed gene annotations (from preprocessing)
ue_pattern <- "../../../out/moca_results/mo_var/preprocessing/*_uniformly_expressed_locations.txt"
ue_files <- Sys.glob(ue_pattern)
default_ue_file <- if (length(ue_files) > 0) sort(ue_files, decreasing = TRUE)[1] else ue_pattern

# Output directory
default_output_dir <- "../../../out/moca_results/mo_var/analysis"

# Parse arguments
MOTIF_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_motif_file
DE_FILE <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_de_file
UE_FILE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_ue_file
OUTPUT_DIR <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_output_dir

# Analysis parameters
DATE_STAMP <- format(Sys.Date(), "%Y%m%d")
PROJECT_NAME <- "vitis_genotype_variance"
WORD_SIZE <- 14  # Minimum motif length (as in original)

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

cat("Genotype Variance Analysis:\n")
cat("  Motif occurrences:", MOTIF_FILE, "\n")
cat("  DE gene file:", DE_FILE, "\n")
cat("  UE gene file:", UE_FILE, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("  Minimum motif size:", WORD_SIZE, "bp\n")
cat("  Date:", DATE_STAMP, "\n\n")

######################################################################################
# SECTION 1: LOAD AND PARSE DIFFERENTIALLY EXPRESSED GENE ANNOTATIONS
######################################################################################
#
# DESCRIPTION:
# Load gene annotations for differentially expressed (DE) genes across genotypes
# Parse the GFF-style format created during preprocessing
#
# INPUT FORMAT (from original):
# genotype:chr:feature:start-end:strand:attributes
#
# PROCESSING STEPS:
# 1. Split composite location string
# 2. Extract gene IDs from attributes
# 3. Create location identifiers (loc) for merging with motif data
# 4. Add flanking regions (±1000bp)
#
######################################################################################

cat("\n=== SECTION 1: Load Differentially Expressed Gene Annotations ===\n")

load_gene_annotations <- function(annot_file, gene_type = "differential") {
  cat("  Loading", gene_type, "gene annotations from:", annot_file, "\n")

  if (!file.exists(annot_file)) {
    cat("  WARNING: File not found - returning empty data frame\n")
    return(data.frame(
      loc = character(),
      gene_id = character(),
      chromosome = character(),
      start = numeric(),
      end = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  # ORIGINAL LOGIC (lines 29-51 of original script):
  # Read tab-delimited file without headers
  gene_annot <- read.table(annot_file,
                           header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE)

  # STEP 1: Parse composite location string
  # Original uses strsplit with ":" separator
  split_columns <- strsplit(gene_annot$V1, ":")
  split_df <- as.data.frame(do.call(rbind, split_columns))
  colnames(split_df) <- paste("Column", 1:ncol(split_df), sep = "")

  # Combine split columns with remaining data
  gene_annot_parsed <- cbind(split_df, gene_annot[-1])

  # STEP 2: Filter for gene features only
  gene_annot_parsed <- gene_annot_parsed %>%
    filter(V3 == "gene")

  # STEP 3: Extract gene IDs from attributes column (V9)
  # Original regex: sub("^[^:]+:([^:]+).*", "\\1", V9)
  # Then: sub("^.*?=(.*?);.*", "\\1", gene_id)
  gene_annot_parsed <- gene_annot_parsed %>%
    mutate(
      gene_id = sub("^[^:]+:([^:]+).*", "\\1", V9),
      gene_id = sub("^.*?=(.*?);.*", "\\1", gene_id)
    )

  # STEP 4: Extract coordinates and add flanking regions
  # Select: chromosome (Column2), start (V4), end (V5), gene_id (extracted)
  gene_coords <- gene_annot_parsed[, c(2, 5, 6, 11)]
  colnames(gene_coords) <- c("chromosome", "start", "end", "gene_id")

  # Convert to numeric and add ±1000bp flanking regions
  gene_coords <- gene_coords %>%
    mutate(
      start = as.numeric(start) - 1000,
      end = as.numeric(end) + 1000
    )

  # STEP 5: Create location identifier for merging
  # Format: chr:start-end (used as key to merge with motif data)
  gene_coords$loc <- paste0(
    gene_coords$chromosome, ":",
    gene_coords$start, "-",
    gene_coords$end
  )
  gene_coords$loc <- gsub(" ", "-", gene_coords$loc)

  # Final output: loc and gene_id mapping
  gene_loc_id <- gene_coords[, c("loc", "gene_id")]
  gene_loc_id <- unique(gene_loc_id)

  cat("  Loaded", nrow(gene_loc_id), gene_type, "gene locations\n")

  return(gene_loc_id)
}

# LOAD DIFFERENTIALLY EXPRESSED GENES
# These are genes with variable expression across genotypes
# HYPOTHESIS: Should show enrichment for MUTATED motifs
gene_annot_diff_expr_locID <- load_gene_annotations(DE_FILE, "differentially expressed")

######################################################################################
# SECTION 2: LOAD AND PARSE UNIFORMLY EXPRESSED GENE ANNOTATIONS
######################################################################################
#
# DESCRIPTION:
# Load gene annotations for uniformly expressed (UE) genes across genotypes
# Same processing as DE genes (lines 58-80 of original)
#
# BIOLOGICAL ROLE:
# UE genes serve as CONTROL/BACKGROUND set
# HYPOTHESIS: Should show enrichment for CONSERVED motifs
#
######################################################################################

cat("\n=== SECTION 2: Load Uniformly Expressed Gene Annotations ===\n")

# LOAD UNIFORMLY EXPRESSED GENES
gene_annot_unif_expr_locID <- load_gene_annotations(UE_FILE, "uniformly expressed")

######################################################################################
# SECTION 3: LOAD AND FILTER MOTIF OCCURRENCE DATA
######################################################################################
#
# DESCRIPTION:
# Load motif occurrences from multi-genotype BLAMM mapping
# Apply quality filters:
# - Must be within 1500bp of gene start/end (promoter/terminator regions)
# - Motif length must be >= word_size (default 14bp)
#
# DATA STRUCTURE:
# loc | source | motif | mstart | mend | score | strand | gene_start | gene_end | ...
#
# KEY VARIABLES:
# - loc: Location identifier matching gene annotations (genotype_chr:start-end)
# - motif: EPM identifier (e.g., epm_vitis_ssr_p0m02F)
# - mstart/mend: Motif position within gene region
# - gene_start/gene_end: Gene boundaries
#
######################################################################################

cat("\n=== SECTION 3: Load and Filter Motif Occurrence Data ===\n")

load_motif_occurrences <- function(motif_file, word_size = 14) {
  cat("  Loading motif occurrences from:", motif_file, "\n")

  if (!file.exists(motif_file)) {
    cat("  WARNING: File not found - returning empty data frame\n")
    return(data.frame())
  }

  # ORIGINAL LOGIC (lines 82-94):
  # Read motif-gene matches (output from BLAMM + occ_filter)
  motif_gene_matches <- read.table(
    motif_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  cat("  Loaded", nrow(motif_gene_matches), "total motif occurrences\n")

  # FILTER 1: Proximity to gene boundaries (±1500bp)
  # Logic: motif must be within 1500bp of TSS OR TTS
  # mstart <= 1500 (near start) OR
  # abs(mstart - gene_length) <= 1500 (near end)
  motif_filtered <- motif_gene_matches[
    motif_gene_matches$mstart <= 1500 |
      abs(motif_gene_matches$mstart - (motif_gene_matches$gene_end - motif_gene_matches$gene_start) + 1) <= 1500,
  ]

  cat("  After proximity filter:", nrow(motif_filtered), "occurrences\n")

  # FILTER 2: Minimum motif length
  # Calculate motif length and filter
  motif_filtered$m_len <- abs(motif_filtered$mstart - motif_filtered$mend)
  motif_filtered <- motif_filtered %>%
    filter(m_len >= word_size)

  cat("  After length filter (>=", word_size, "bp):", nrow(motif_filtered), "occurrences\n")

  return(motif_filtered)
}

# LOAD FILTERED MOTIF OCCURRENCES
motif_gene_matches_filtered <- load_motif_occurrences(MOTIF_FILE, WORD_SIZE)

######################################################################################
# SECTION 4: MERGE MOTIF DATA WITH GENE ANNOTATIONS (DIFFERENTIAL)
######################################################################################
#
# DESCRIPTION:
# Merge motif occurrences with differentially expressed gene annotations
# This creates the primary dataset for variance analysis in DE genes
#
# ANALYSIS STEPS:
# 1. Merge motif data with DE gene locations (by "loc" identifier)
# 2. Extract genotype prefix from location string
# 3. Count unique genotypes in dataset
# 4. Count motif occurrences per gene
# 5. Classify motifs as "conserved" or "mutated" based on genotype coverage
#
# CONSERVATION LOGIC (line 106):
# - CONSERVED: n / n_genotypes = integer (motif present in ALL genotypes)
# - MUTATED: n / n_genotypes ≠ integer (motif present in SOME genotypes)
#
######################################################################################

cat("\n=== SECTION 4: Analyze Motif Variance in Differentially Expressed Genes ===\n")

analyze_motif_variance_de <- function(motif_data, gene_annotations) {
  cat("  Merging motif data with DE gene annotations...\n")

  # STEP 1: Merge by location identifier
  # ORIGINAL LOGIC (line 96-98):
  merged_df <- merge(
    motif_data,
    gene_annotations,
    by = "loc"
  )

  cat("  Merged dataset:", nrow(merged_df), "motif-gene associations\n")

  if (nrow(merged_df) == 0) {
    cat("  WARNING: No overlapping data - returning empty results\n")
    return(list(
      variance_counts = data.frame(),
      motif_counts = data.frame(),
      n_genotypes = 0
    ))
  }

  # STEP 2: Extract genotype prefix from location string
  # ORIGINAL LOGIC (line 99-101):
  # Format: genotype_chr:start-end -> extract "genotype"
  merged_df <- merged_df %>%
    mutate(loc_prefix = sub("^(.*?)_.*", "\\1", loc))

  genotypes_diff <- length(unique(merged_df$loc_prefix))
  cat("  Identified", genotypes_diff, "unique genotypes in DE genes\n")

  # STEP 3: Count motif occurrences per gene
  # ORIGINAL LOGIC (line 102-104):
  # Group by gene_id and count each motif
  motif_counts <- merged_df %>%
    group_by(gene_id) %>%
    count(motif)  # Creates column 'n' with count per gene-motif combination

  cat("  Total gene-motif combinations:", nrow(motif_counts), "\n")

  # STEP 4: Classify conservation status
  # ORIGINAL LOGIC (line 105-106):
  # If count divides evenly by number of genotypes, motif is CONSERVED
  # Otherwise, it's MUTATED (present in some but not all genotypes)
  #
  # MATHEMATICAL INTERPRETATION:
  # n=6, genotypes=3 -> n/3 = 2 (integer) -> CONSERVED (2 copies per genotype)
  # n=5, genotypes=3 -> n/3 = 1.67 (not integer) -> MUTATED (variable across genotypes)
  motif_counts <- motif_counts %>%
    mutate(
      variance = ifelse(
        n / genotypes_diff == floor(n / genotypes_diff),
        "conserved",
        "mutated"
      )
    )

  # STEP 5: Summarize variance distribution
  variance_counts <- table(motif_counts$variance)

  cat("  Variance distribution in DE genes:\n")
  cat("    Conserved motifs:", variance_counts["conserved"], "\n")
  cat("    Mutated motifs:", variance_counts["mutated"], "\n")

  return(list(
    variance_counts = as.data.frame(variance_counts),
    motif_counts = motif_counts,
    n_genotypes = genotypes_diff
  ))
}

# RUN ANALYSIS ON DIFFERENTIALLY EXPRESSED GENES
de_results <- analyze_motif_variance_de(motif_gene_matches_filtered, gene_annot_diff_expr_locID)

######################################################################################
# SECTION 5: MERGE MOTIF DATA WITH GENE ANNOTATIONS (UNIFORM)
######################################################################################
#
# DESCRIPTION:
# Analyze motif variance in uniformly expressed genes (CONTROL SET)
# Same logic as Section 4, but for UE genes
#
# BOOTSTRAP RESAMPLING (lines 122-141):
# Because UE gene set is much larger than DE gene set, we need to:
# 1. Randomly sample UE genes to match DE set size (for fair comparison)
# 2. Repeat sampling 100 times (bootstrap)
# 3. Calculate average conserved/mutated counts
# 4. Use average as representative UE baseline
#
######################################################################################

cat("\n=== SECTION 5: Analyze Motif Variance in Uniformly Expressed Genes ===\n")

analyze_motif_variance_ue <- function(motif_data, gene_annotations, sample_size = NULL, n_iterations = 100) {
  cat("  Merging motif data with UE gene annotations...\n")

  # STEP 1: Merge by location identifier
  merged_df <- merge(
    motif_data,
    gene_annotations,
    by = "loc"
  )

  cat("  Merged dataset:", nrow(merged_df), "motif-gene associations\n")

  if (nrow(merged_df) == 0) {
    cat("  WARNING: No overlapping data - returning empty results\n")
    return(list(
      variance_counts = data.frame(),
      motif_counts = data.frame(),
      n_genotypes = 0
    ))
  }

  # STEP 2: Extract genotype information
  merged_df <- merged_df %>%
    mutate(loc_prefix = sub("^(.*?)_.*", "\\1", loc))

  genotypes_unif <- length(unique(merged_df$loc_prefix))
  cat("  Identified", genotypes_unif, "unique genotypes in UE genes\n")

  # STEP 3: Count motif occurrences
  motif_counts <- merged_df %>%
    group_by(gene_id) %>%
    count(motif)

  cat("  Total gene-motif combinations:", nrow(motif_counts), "\n")

  # STEP 4: Classify conservation status
  motif_counts <- motif_counts %>%
    mutate(
      variance = ifelse(
        n / genotypes_unif == floor(n / genotypes_unif),
        "conserved",
        "mutated"
      )
    )

  # STEP 5: BOOTSTRAP RESAMPLING (if sample_size provided)
  # ORIGINAL LOGIC (lines 122-141):
  if (!is.null(sample_size) && sample_size < nrow(motif_counts)) {
    cat("  Performing bootstrap resampling...\n")
    cat("  Sample size:", sample_size, "| Iterations:", n_iterations, "\n")

    # Ungroup for sampling
    ungrouped_data <- motif_counts %>% ungroup()

    # Storage for bootstrap results
    conserved_counts <- numeric(n_iterations)
    mutated_counts <- numeric(n_iterations)

    # BOOTSTRAP LOOP
    for (i in 1:n_iterations) {
      # Random seed for reproducibility tracking
      iteration_seed <- sample.int(10^5, 1)
      set.seed(iteration_seed)

      # Random sample without replacement
      random_subset <- ungrouped_data %>%
        sample_n(size = sample_size, replace = FALSE)

      # Count variance categories
      variance_subset_counts <- table(random_subset$variance)
      conserved_counts[i] <- variance_subset_counts["conserved"]
      mutated_counts[i] <- variance_subset_counts["mutated"]
    }

    # Calculate bootstrap averages
    average_conserved <- mean(conserved_counts)
    average_mutated <- mean(mutated_counts)

    cat("  Bootstrap results:\n")
    cat("    Mean conserved:", round(average_conserved, 1), "\n")
    cat("    Mean mutated:", round(average_mutated, 1), "\n")

    # Create summary table
    variance_counts <- data.frame(
      Var1 = c("conserved", "mutated"),
      Freq = c(average_conserved, average_mutated)
    )
  } else {
    # No resampling - use full dataset
    variance_counts <- as.data.frame(table(motif_counts$variance))
  }

  cat("  Variance distribution in UE genes:\n")
  cat("    Conserved motifs:", variance_counts$Freq[variance_counts$Var1 == "conserved"], "\n")
  cat("    Mutated motifs:", variance_counts$Freq[variance_counts$Var1 == "mutated"], "\n")

  return(list(
    variance_counts = variance_counts,
    motif_counts = motif_counts,
    n_genotypes = genotypes_unif
  ))
}

# RUN ANALYSIS ON UNIFORMLY EXPRESSED GENES
# Use bootstrap with sample size matching DE set (for fair comparison)
de_sample_size <- nrow(de_results$motif_counts)
ue_results <- analyze_motif_variance_ue(
  motif_gene_matches_filtered,
  gene_annot_unif_expr_locID,
  sample_size = de_sample_size,
  n_iterations = 100
)

######################################################################################
# SECTION 6: STATISTICAL COMPARISON - CONSERVATION RATES
######################################################################################
#
# DESCRIPTION:
# Compare motif conservation rates between DE and UE genes
#
# KEY COMPARISON:
# - DE genes: Expected to have MORE mutated motifs (regulatory divergence)
# - UE genes: Expected to have MORE conserved motifs (regulatory stability)
#
# VISUALIZATION: Normalized bar plot (lines 143-164)
#
######################################################################################

cat("\n=== SECTION 6: Statistical Comparison of Conservation Rates ===\n")

# ORIGINAL LOGIC (lines 143-156):
# Merge DE and UE variance counts
# Normalize by total to get percentages
# Prepare for visualization

create_comparison_data <- function(de_variance, ue_variance) {
  # Ensure consistent naming
  colnames(de_variance) <- c("Variance", "Diff_Count")
  colnames(ue_variance) <- c("Variance", "Unif_Count")

  # Merge data frames
  combined_df <- merge(de_variance, ue_variance, by = "Variance", all = TRUE)

  # Handle NAs (in case one category is missing)
  combined_df[is.na(combined_df)] <- 0

  # Normalize to percentages
  combined_df$Normalized_Diff_Count <- combined_df$Diff_Count / sum(combined_df$Diff_Count)
  combined_df$Normalized_Unif_Count <- combined_df$Unif_Count / sum(combined_df$Unif_Count)

  # Reshape for plotting
  rearranged_df <- combined_df %>%
    pivot_longer(
      cols = starts_with("Normalized"),
      names_to = "Source",
      values_to = "Normalized"
    ) %>%
    mutate(
      Source = gsub("Normalized_", "", Source),
      Source = gsub("_Count", "", Source)
    )

  return(list(
    combined = combined_df,
    plot_data = rearranged_df
  ))
}

# PSEUDOCODE - run when data available
# comparison_data <- create_comparison_data(de_results$variance_counts, ue_results$variance_counts)

cat("  PLACEHOLDER: Statistical comparison logic defined\n")
cat("  TODO: Run when both DE and UE results available\n")

######################################################################################
# SECTION 7: BOOTSTRAP ANALYSIS - GENE-LEVEL CONSERVATION
######################################################################################
#
# DESCRIPTION:
# More detailed analysis: gene-level conservation patterns
# Identifies genes with ONLY conserved motifs vs genes with ANY mutated motifs
#
# BIOLOGICAL INTERPRETATION:
# - Genes with only conserved motifs: Stable regulatory architecture
# - Genes with mutated motifs: Variable regulatory control
#
# BOOTSTRAP STRATEGY (lines 181-221 for DE, 271-314 for UE):
# 1. Sample 100 random genes (repeated 1000 times)
# 2. Separate genes by motif variance status
# 3. Count genes with ONLY conserved vs ANY mutated motifs
# 4. Calculate error rates (proportion with mutations)
# 5. Generate confidence intervals from bootstrap distribution
#
######################################################################################

cat("\n=== SECTION 7: Gene-Level Conservation Bootstrap Analysis ===\n")

bootstrap_gene_conservation <- function(motif_counts_df, n_sample = 100, n_iterations = 1000, gene_type = "DE") {
  cat("  Running bootstrap analysis for", gene_type, "genes...\n")
  cat("  Sample size:", n_sample, "| Iterations:", n_iterations, "\n")

  # ORIGINAL LOGIC (lines 181-221):
  # Storage for results
  average_conserved_genes <- numeric(n_iterations)
  average_mutated_genes <- numeric(n_iterations)
  error_rates <- numeric(n_iterations)
  data_for_boxplots <- list()

  # BOOTSTRAP LOOP
  for (i in 1:n_iterations) {
    # STEP 1: Random sample of gene-motif combinations
    sampled_indices <- sample(nrow(motif_counts_df), n_sample)
    sampled_data <- motif_counts_df[sampled_indices, ]

    # STEP 2: Separate by variance status
    motifs_mutated <- sampled_data[sampled_data$variance == "mutated", ]
    motifs_conserved <- sampled_data[sampled_data$variance == "conserved", ]

    # STEP 3: Identify genes with ONLY conserved motifs
    # Anti-join: genes in conserved set that are NOT in mutated set
    genes_only_conserved <- anti_join(motifs_conserved, motifs_mutated, by = "gene_id")

    # STEP 4: Identify genes with ANY mutated motifs
    # Anti-join: genes in mutated set that are NOT in the conserved-only set
    genes_with_mutations <- anti_join(motifs_mutated, genes_only_conserved, by = "gene_id")

    # STEP 5: Count unique genes
    n_conserved <- length(unique(genes_only_conserved$gene_id))
    n_mutated <- length(unique(genes_with_mutations$gene_id))

    # STEP 6: Store results
    average_conserved_genes[i] <- n_conserved
    average_mutated_genes[i] <- n_mutated
    error_rates[i] <- n_mutated / (n_conserved + n_mutated)  # Proportion with mutations

    # Store for boxplot
    data_for_boxplots[[i]] <- c(n_conserved, n_mutated)
  }

  # CALCULATE STATISTICS
  avg_conserved <- mean(average_conserved_genes)
  avg_mutated <- mean(average_mutated_genes)
  avg_error_rate <- mean(error_rates)
  error_rate_sd <- sd(error_rates)

  # Calculate statistics for boxplot
  conserved_values <- sapply(data_for_boxplots, "[[", 1)
  mutated_values <- sapply(data_for_boxplots, "[[", 2)

  cat("  Results for", gene_type, "genes:\n")
  cat("    Mean genes with only conserved motifs:", round(avg_conserved, 1), "\n")
  cat("    Mean genes with mutated motifs:", round(avg_mutated, 1), "\n")
  cat("    Mean mutation rate:", round(avg_error_rate, 3), "\n")
  cat("    SD of mutation rate:", round(error_rate_sd, 3), "\n")

  return(list(
    mean_conserved = avg_conserved,
    mean_mutated = avg_mutated,
    mean_error_rate = avg_error_rate,
    sd_error_rate = error_rate_sd,
    conserved_distribution = conserved_values,
    mutated_distribution = mutated_values,
    boxplot_data = data_for_boxplots
  ))
}

# PSEUDOCODE - run when data available
# de_bootstrap <- bootstrap_gene_conservation(de_results$motif_counts, gene_type = "DE")
# ue_bootstrap <- bootstrap_gene_conservation(ue_results$motif_counts, gene_type = "UE")

cat("  PLACEHOLDER: Bootstrap analysis logic defined\n")
cat("  TODO: Run when motif count data available\n")

######################################################################################
# SECTION 8: CREATE GENE LISTS FOR DOWNSTREAM ANALYSIS
######################################################################################
#
# DESCRIPTION:
# Generate comprehensive gene lists classified by expression pattern and motif status
# These lists can be used for:
# - GO enrichment analysis
# - Pathway analysis
# - Further regulatory element investigation
#
# GENE CATEGORIES (lines 400-425):
# 1. genes_conserved_with_diff_expr: DE genes with only conserved motifs
# 2. genes_mutated_with_diff_expr: DE genes with mutated motifs
# 3. genes_conserved_with_unif_expr: UE genes with only conserved motifs
# 4. genes_mutated_with_unif_expr: UE genes with mutated motifs
#
######################################################################################

cat("\n=== SECTION 8: Create Classified Gene Lists ===\n")

create_gene_lists <- function(de_motif_counts, ue_motif_counts) {
  cat("  Classifying genes by expression and motif conservation status...\n")

  # ORIGINAL LOGIC (lines 400-425):
  # DIFFERENTIAL EXPRESSION GENES
  # Separate by variance
  de_motifs_mutated <- de_motif_counts[de_motif_counts$variance == "mutated", ]
  de_motifs_conserved <- de_motif_counts[de_motif_counts$variance == "conserved", ]

  # Identify exclusive sets
  genes_conserved_diff_expr <- anti_join(de_motifs_conserved, de_motifs_mutated, by = "gene_id")
  genes_mutated_diff_expr <- anti_join(de_motifs_mutated, genes_conserved_diff_expr, by = "gene_id")

  # UNIFORM EXPRESSION GENES
  ue_motifs_mutated <- ue_motif_counts[ue_motif_counts$variance == "mutated", ]
  ue_motifs_conserved <- ue_motif_counts[ue_motif_counts$variance == "conserved", ]

  genes_conserved_unif_expr <- anti_join(ue_motifs_conserved, ue_motifs_mutated, by = "gene_id")
  genes_mutated_unif_expr <- anti_join(ue_motifs_mutated, genes_conserved_unif_expr, by = "gene_id")

  # Add source labels
  genes_conserved_diff_expr$source_df <- "genes_conserved_with_diff_expr"
  genes_mutated_diff_expr$source_df <- "genes_mutated_with_diff_expr"
  genes_conserved_unif_expr$source_df <- "genes_conserved_with_unif_expr"
  genes_mutated_unif_expr$source_df <- "genes_mutated_with_unif_expr"

  # Combine into single data frame
  combined_genes <- bind_rows(
    genes_conserved_diff_expr,
    genes_mutated_diff_expr,
    genes_conserved_unif_expr,
    genes_mutated_unif_expr
  )

  cat("  Gene classification summary:\n")
  cat("    DE genes with conserved motifs:", nrow(genes_conserved_diff_expr), "\n")
  cat("    DE genes with mutated motifs:", nrow(genes_mutated_diff_expr), "\n")
  cat("    UE genes with conserved motifs:", nrow(genes_conserved_unif_expr), "\n")
  cat("    UE genes with mutated motifs:", nrow(genes_mutated_unif_expr), "\n")

  return(list(
    combined = combined_genes,
    de_conserved = genes_conserved_diff_expr,
    de_mutated = genes_mutated_diff_expr,
    ue_conserved = genes_conserved_unif_expr,
    ue_mutated = genes_mutated_unif_expr
  ))
}

# PSEUDOCODE - run when data available
# gene_lists <- create_gene_lists(de_results$motif_counts, ue_results$motif_counts)

cat("  PLACEHOLDER: Gene list creation logic defined\n")
cat("  TODO: Run when motif count data available\n")

######################################################################################
# SECTION 9: VISUALIZATION
######################################################################################
#
# DESCRIPTION:
# Generate publication-quality figures
#
# FIGURE 1: Normalized bar plot - Conservation rates (lines 158-164)
# FIGURE 2: Boxplot - Gene-level conservation distribution (lines 359-397)
#
######################################################################################

cat("\n=== SECTION 9: Generate Visualizations ===\n")

# FIGURE 1: Conservation rate comparison bar plot
plot_conservation_rates <- function(plot_data) {
  # ORIGINAL LOGIC (lines 158-164):
  p <- ggplot(plot_data, aes(x = Source, y = Normalized, fill = Variance)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
      x = "Gene Set",
      y = "Proportion",
      fill = "Motif Status",
      title = "Motif Conservation Rates: Differential vs Uniform Expression",
      subtitle = "Percentage of conserved vs mutated motifs across genotypes"
    ) +
    scale_fill_manual(values = c("conserved" = "lightgray", "mutated" = "darkgray")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}

# FIGURE 2: Gene-level conservation boxplot
plot_gene_conservation_boxplot <- function(de_bootstrap, ue_bootstrap) {
  # ORIGINAL LOGIC (lines 359-397):
  # Calculate boxplot statistics from bootstrap distributions

  # DE genes
  de_conserved_vals <- de_bootstrap$conserved_distribution
  de_mutated_vals <- de_bootstrap$mutated_distribution

  boxplot_de <- data.frame(
    Category = factor(c("DE_conserved", "DE_mutated")),
    Min = c(min(de_conserved_vals), min(de_mutated_vals)),
    Q1 = c(quantile(de_conserved_vals, 0.25), quantile(de_mutated_vals, 0.25)),
    Median = c(mean(de_conserved_vals), mean(de_mutated_vals)),
    Q3 = c(quantile(de_conserved_vals, 0.75), quantile(de_mutated_vals, 0.75)),
    Max = c(max(de_conserved_vals), max(de_mutated_vals))
  )

  # UE genes
  ue_conserved_vals <- ue_bootstrap$conserved_distribution
  ue_mutated_vals <- ue_bootstrap$mutated_distribution

  boxplot_ue <- data.frame(
    Category = factor(c("UE_conserved", "UE_mutated")),
    Min = c(min(ue_conserved_vals), min(ue_mutated_vals)),
    Q1 = c(quantile(ue_conserved_vals, 0.25), quantile(ue_mutated_vals, 0.25)),
    Median = c(mean(ue_conserved_vals), mean(ue_mutated_vals)),
    Q3 = c(quantile(ue_conserved_vals, 0.75), quantile(ue_mutated_vals, 0.75)),
    Max = c(max(ue_conserved_vals), max(ue_mutated_vals))
  )

  # Combine
  combined_boxplot_data <- bind_rows(boxplot_ue, boxplot_de)

  # Custom colors (from original)
  custom_colors <- c(
    "DE_conserved" = "cyan3",
    "DE_mutated" = "gold",
    "UE_conserved" = "cyan4",
    "UE_mutated" = "darkgoldenrod2"
  )

  p <- ggplot(combined_boxplot_data,
              aes(x = Category, ymin = Min, lower = Q1, middle = Median,
                  upper = Q3, ymax = Max, fill = Category)) +
    geom_boxplot(stat = "identity", width = 0.5) +
    scale_fill_manual(values = custom_colors) +
    labs(
      x = "Gene Category",
      y = "Number of Genes",
      title = "Gene-Level Motif Conservation Patterns",
      subtitle = "Bootstrap distribution (1000 iterations, 100 genes per sample)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.major.x = element_blank()
    )

  return(p)
}

cat("  PLACEHOLDER: Visualization functions defined\n")
cat("  TODO: Generate plots when analysis results available\n")

######################################################################################
# SECTION 10: WRITE OUTPUT FILES
######################################################################################

cat("\n=== SECTION 10: Write Output Files ===\n")

# Output file paths
output_base <- file.path(OUTPUT_DIR, paste0(DATE_STAMP, "_", PROJECT_NAME))

# PSEUDOCODE - write when data available:
# write.csv(gene_lists$combined, paste0(output_base, "_combined_gene_lists.csv"), row.names = FALSE)
# write.csv(comparison_data$combined, paste0(output_base, "_variance_comparison.csv"), row.names = FALSE)
# ggsave(paste0(output_base, "_conservation_rates.pdf"), plot = conservation_plot, width = 10, height = 8)
# ggsave(paste0(output_base, "_gene_conservation_boxplot.pdf"), plot = boxplot, width = 12, height = 8)

cat("  Output files will be saved to:", OUTPUT_DIR, "\n")
cat("  Base filename:", basename(output_base), "\n")

######################################################################################
# ANALYSIS COMPLETE
######################################################################################

cat("\n=== GENOTYPE VARIANCE ANALYSIS SETUP COMPLETE ===\n")
cat("\nTHIS SCRIPT IS A TEMPLATE/FRAMEWORK\n")
cat("All major analysis steps are defined with detailed annotations\n")
cat("Ready to be populated with actual multi-genotype data\n\n")

cat("KEY STEPS TO ACTIVATE THIS ANALYSIS:\n")
cat("1. Run preprocessing script (00_prepare_genotype_gene_lists.R)\n")
cat("2. Extract FASTA sequences for DE and UE genes\n")
cat("3. Run BLAMM on merged multi-genotype FASTA\n")
cat("4. Filter BLAMM output with filter_blamm_occ_00.R\n")
cat("5. Execute this script with actual data\n\n")

cat("OUTPUT DIRECTORY:", OUTPUT_DIR, "\n")
cat("Session completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")