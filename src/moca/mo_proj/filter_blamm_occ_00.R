#!/usr/bin/env Rscript
# Filters BLAMM genomic coordinate occurrence results using biologically meaningful criteria
# Enhanced from occ_filter_v1.1.R for genomic coordinate workflow
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca/mo_proj")
######################

library(tidyr)
library(dplyr)
library(readr)
library(magrittr)

# Source utility functions
source("../utils.R")

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Default file paths for genomic coordinate workflow
default_blamm_dir <- "../../../out/moca_results/mo_proj/genomic_vitis_vinifera_PN40024_1000bp_pt0.0001"
default_output_dir <- "../../../out/moca_results/mo_proj"
default_output_file <- "filtered_genomic_occurrences.txt"

# Parse arguments
INPUT_TYPE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "single"  # "chunks" or "single"
INPUT_PATH <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_blamm_dir
OUTPUT_DIR <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_output_dir
OUTPUT_FILE <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_output_file
ANNOTATION_FILE <- if (length(args) >= 5 && nzchar(args[5])) args[5] else "../../../vitis_data/gene_models/vitis_vinifera_PN40024.gff3"
FLANK_SIZE <- if (length(args) >= 6 && nzchar(args[6])) as.numeric(args[6]) else 1000

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

# Full output path
output_path <- file.path(OUTPUT_DIR, OUTPUT_FILE)

cat("BLAMM Genomic Occurrence Filtering Parameters:\n")
cat("  Input type:", INPUT_TYPE, "\n")
cat("  Input path:", INPUT_PATH, "\n")
cat("  Output file:", output_path, "\n")
cat("  Annotation file:", ANNOTATION_FILE, "\n")
cat("  Flank size:", FLANK_SIZE, "bp\n")
cat("  Filtering: 1000bp upstream + 500bp 5'UTR + 500bp 3'UTR + 1000bp downstream (matches model training)\n\n")

#' Load gene annotation data for coordinate conversion
#' 
#' @param annotation_file Path to GFF3 annotation file
#' @return Data frame with gene coordinates
load_gene_annotations <- function(annotation_file) {
  
  if (!file.exists(annotation_file)) {
    stop("Annotation file not found: ", annotation_file)
  }
  
  cat("Loading gene annotations from:", annotation_file, "\n")
  
  # Read GFF3 file, skip comment lines
  gff_data <- read_tsv(annotation_file, 
                     col_names = c("seqname", "source", "feature", "start", "end", 
                                   "score", "strand", "frame", "attribute"),
                     col_types = "ccciicccc",
                     comment = "#",
                     show_col_types = FALSE)

  # Debug: Check what got loaded
 cat("DEBUG - Loaded", nrow(gff_data), "total rows from GFF3\n")
 cat("DEBUG - Chromosomes present:", length(unique(gff_data$seqname)), "\n")
 cat("DEBUG - Chr12-16 rows loaded:", sum(grepl("^chr1[2-6]", gff_data$seqname)), "\n")
  
  colnames(gff_data) <- c("seqname", "source", "feature", "start", "end", 
                          "score", "strand", "frame", "attribute")
  
  # Filter for gene features
  genes <- gff_data %>%
    filter(feature == "gene") %>%
    mutate(
      # Extract gene ID from attributes
      gene_id = gsub(".*ID=([^;]+).*", "\\1", attribute)
    ) %>%
    select(gene_id, seqname, start, end, strand) %>%
    rename(chr = seqname, gene_start = start, gene_end = end, gene_strand = strand)
  
  # DEBUG: Check gene annotation loading
  cat("DEBUG - Gene annotations loaded:\n")
  cat("  Total genes:", nrow(genes), "\n")
  gene_chr_counts <- table(genes$chr)
  cat("  Genes per chromosome:\n")
  for (chr in names(gene_chr_counts)) {
    cat("    ", chr, ":", gene_chr_counts[chr], "\n")
  }
  
  # Check for chromosomes 12-16 specifically
  chr_12_16 <- genes$chr %in% c("chr12", "chr13", "chr14", "chr15", "chr16")
  cat("  Genes on chr12-16:", sum(chr_12_16), "\n")
  
  cat("Loaded", nrow(genes), "gene annotations\n")
  return(genes)
}

#' Filter BLAMM genomic occurrence data using model training criteria
#' 
#' This function processes BLAMM occurrence files from genomic coordinate mapping
#' and filters motifs based on regions used in model training:
#' - 1000bp upstream flank (promoter region)
#' - 500bp 5' UTR (beginning of gene)
#' - 500bp 3' UTR (end of gene)  
#' - 1000bp downstream flank (terminator region)
#' - Excludes coding regions (middle of gene not used in training)
#' 
#' @param file_paths Character vector of file paths to process
#' @param output_file Character path for output file
#' @param gene_annotations Data frame with gene coordinate information
#' @param flank_size Numeric size of flanking regions extracted
#' @return Invisible TRUE on success
filter_genomic_occurrences <- function(file_paths, output_file, gene_annotations, flank_size) {
  
  # Initialize counters
  combined_df <- data.frame()
  total_files <- length(file_paths)
  total_input_occurrences <- 0
  
  cat("Processing", total_files, "occurrence file(s)...\n")
  
  # Process each file
  for (i in seq_along(file_paths)) {
    filepath <- file_paths[i]
    
    cat("  Processing file", i, "of", total_files, ":", basename(filepath), "\n")
    
    # Read occurrence file
    occurrence_df <- read.table(filepath, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
    
    # Handle potential empty files
    if (nrow(occurrence_df) == 0) {
      cat("    Warning: Empty file, skipping\n")
      next
    }
    
    # Set column names based on BLAMM genomic output format
    colnames(occurrence_df) <- c("gene_id", "source", "motif", "mstart", "mend", 
                                "score", "strand", "V7", "V8")
    
    # DEBUG: Check BLAMM data loading for this file
    cat("    BLAMM data for file", i, ":\n")
    cat("      Total occurrences:", nrow(occurrence_df), "\n")
    
    # Extract chromosome from gene_id for debugging
    occurrence_df$debug_chr <- gsub(".*chr(\\d+)g.*", "chr\\1", occurrence_df$gene_id)
    blamm_chr_counts <- table(occurrence_df$debug_chr)
    cat("      Occurrences per chromosome:\n")
    for (chr in names(blamm_chr_counts)) {
      if (chr %in% c("chr12", "chr13", "chr14", "chr15", "chr16")) {
        cat("        ", chr, ":", blamm_chr_counts[chr], " *** \n")
      } else {
        cat("        ", chr, ":", blamm_chr_counts[chr], "\n")
      }
    }
    
    # Count input occurrences
    file_input_count <- nrow(occurrence_df)
    total_input_occurrences <- total_input_occurrences + file_input_count
    
    # Join with gene annotations to get genomic coordinates
    annotated_df <- occurrence_df %>%
      left_join(gene_annotations, by = "gene_id") %>%
      filter(!is.na(chr)) %>%  # Remove motifs for genes not in annotation
      mutate(
        # Convert motif coordinates: relative position within extracted sequence
        # Extracted sequence = [flank][gene][flank]
        # Gene starts at position (flank_size + 1) in extracted sequence
        
        # Calculate position relative to gene start
        pos_relative_to_gene = mstart - flank_size,
        
        # Calculate gene length
        gene_length = gene_end - gene_start + 1,
        
        # Classify motif location based on model training regions
        motif_region = case_when(
          pos_relative_to_gene <= 0 ~ "upstream_flank",              # 1000bp upstream flank
          pos_relative_to_gene > 0 & pos_relative_to_gene <= 500 ~ "five_prime_utr",  # 5' UTR (500bp into gene)
          pos_relative_to_gene > (gene_length - 500) & pos_relative_to_gene <= gene_length ~ "three_prime_utr",  # 3' UTR (last 500bp of gene)
          pos_relative_to_gene > gene_length & pos_relative_to_gene <= (gene_length + 1000) ~ "downstream_flank",   # 1000bp downstream flank
          TRUE ~ "coding_region"                                     # Middle of gene (exclude)
        )
      )
    
    # Apply biological filtering - keep regions used in model training
    # Keep: 1000bp upstream + 500bp 5'UTR + 500bp 3'UTR + 1000bp downstream
    filtered_df <- annotated_df %>%
      filter(
        motif_region %in% c("upstream_flank", "five_prime_utr", "three_prime_utr", "downstream_flank")
      )
    
    # Append to combined results
    combined_df <- rbind(combined_df, filtered_df)
    
    # Progress update
    file_filtered_count <- nrow(filtered_df)
    cat("    Input:", file_input_count, "occurrences, Filtered:", file_filtered_count, 
        "(", round(100 * file_filtered_count / file_input_count, 1), "%)\n")
    
    # Progress update for large file sets
    if (i %% 10 == 0 || i == total_files) {
      cat("    Processed", i, "/", total_files, "files. Current total rows:", nrow(combined_df), "\n")
    }
  }
  
  # Calculate filtering statistics
  total_output_occurrences <- nrow(combined_df)
  total_filtered_out <- total_input_occurrences - total_output_occurrences
  filter_percentage <- round(100 * total_output_occurrences / total_input_occurrences, 1)
  
  cat("\n=== FILTERING SUMMARY ===\n")
  cat("Total input occurrences:", total_input_occurrences, "\n")
  cat("Total filtered occurrences:", total_output_occurrences, "\n")
  cat("Occurrences removed:", total_filtered_out, "(", 100 - filter_percentage, "%)\n")
  cat("Occurrences retained:", filter_percentage, "%\n")
  
  # Region breakdown
  if (nrow(combined_df) > 0) {
    region_summary <- combined_df %>%
      count(motif_region, name = "count") %>%
      mutate(percentage = round(100 * count / sum(count), 1))
    
    cat("\nRetained motifs by region:\n")
    for (i in 1:nrow(region_summary)) {
      cat("  ", region_summary$motif_region[i], ":", region_summary$count[i], 
          "(", region_summary$percentage[i], "%)\n")
    }
  }
  
  # Write results to output file
  cat("\nWriting results to:", output_file, "\n")
  write.table(combined_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Genomic occurrence filtering completed successfully!\n")
  return(invisible(TRUE))
}

# Main execution logic
cat("Starting BLAMM genomic occurrence filtering...\n")

# Load gene annotations
gene_annotations <- load_gene_annotations(ANNOTATION_FILE)

if (INPUT_TYPE == "chunks") {
  # Process chunked files from directory
  if (!dir.exists(INPUT_PATH)) {
    stop("Chunks directory not found: ", INPUT_PATH, 
         "\nMake sure the directory exists and contains chunked occurrence files.")
  }
  
  # Find all chunk files (pattern matches "chunk_*" or "smallfile*")
  chunk_patterns <- c("^chunk_", "^smallfile")
  chunk_files <- c()
  
  for (pattern in chunk_patterns) {
    found_files <- list.files(INPUT_PATH, pattern = pattern, full.names = TRUE)
    chunk_files <- c(chunk_files, found_files)
  }
  
  if (length(chunk_files) == 0) {
    stop("No chunk files found in ", INPUT_PATH, 
         "\nLooking for files matching patterns: ", paste(chunk_patterns, collapse = ", "))
  }
  
  cat("Found", length(chunk_files), "chunk files to process.\n")
  file_paths <- sort(chunk_files)  # Process in order
  
} else if (INPUT_TYPE == "single") {
  # Process single large occurrence file
  occurrence_file <- file.path(INPUT_PATH, "occurrences.txt")
  
  if (!file.exists(occurrence_file)) {
    # Try alternative file patterns
    alt_files <- list.files(INPUT_PATH, pattern = "occurrences.*\\.txt$", full.names = TRUE)
    if (length(alt_files) > 0) {
      occurrence_file <- alt_files[1]
      cat("Using occurrence file:", occurrence_file, "\n")
    } else {
      stop("Occurrence file not found in ", INPUT_PATH, 
           "\nExpected: occurrences.txt or similar pattern")
    }
  }
  
  file_paths <- occurrence_file
  
} else {
  stop("Invalid INPUT_TYPE: ", INPUT_TYPE, ". Must be 'chunks' or 'single'")
}

# Execute filtering
filter_genomic_occurrences(file_paths, output_path, gene_annotations, FLANK_SIZE)

# Print session info for reproducibility
cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Script completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")