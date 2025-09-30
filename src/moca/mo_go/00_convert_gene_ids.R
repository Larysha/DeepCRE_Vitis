#!/usr/bin/env Rscript
# Convert GO gene IDs from standard format to match motif data naming convention
# Input: V_vinifera_ontology.gmt (GO format: Vitvi08g02097)
# Output: V_vinifera_ont_converted.gmt (Motif format: Vitvi05_01chr08g02097)
# Optimized GO Gene ID Conversion - uses hash table for faster lookups
# Processes full dataset efficiently by building lookup tables
# this script only needs to be run once per go file / genome annotation 
######################

library(stringr)

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# File paths
default_input_go <- "../../../vitis_data/gene_ontology/V_vinifera_ontology.gmt"
default_output_go <- "../../../vitis_data/gene_ontology/V_vinifera_ont_converted.gmt"
default_genome_annotation <- "../../../vitis_data/gene_models/vitis_vinifera_PN40024.gff3"

INPUT_GO <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_input_go
OUTPUT_GO <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_output_go
GENOME_ANNOTATION <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_genome_annotation

cat("Optimized GO Gene ID Conversion:\n")
cat("  Input GO file:", INPUT_GO, "\n")
cat("  Output GO file:", OUTPUT_GO, "\n")
cat("  Reference genome annotation:", GENOME_ANNOTATION, "\n\n")

#' Build hash table for fast gene ID lookups
build_gene_lookup_table <- function(gff_file) {
  cat("Building gene ID lookup table...\n")
  
  if (!file.exists(gff_file)) {
    stop("Genome annotation file not found: ", gff_file)
  }
  
  # Read gene features only using awk for efficiency
  cat("  Extracting gene features from GFF3...\n")
  temp_file <- tempfile()
  system(sprintf("grep '\tgene\t' %s > %s", gff_file, temp_file))
  
  gff_genes <- read.table(temp_file, sep = "\t", stringsAsFactors = FALSE, quote = "")
  unlink(temp_file)
  
  cat("  Found", nrow(gff_genes), "gene features\n")
  
  # Extract gene IDs from attributes
  gene_ids <- sapply(gff_genes$V9, function(attr) {
    id_match <- stringr::str_extract(attr, "ID=([^;]+)")
    if (!is.na(id_match)) {
      return(stringr::str_remove(id_match, "ID="))
    }
    return(NA)
  })
  
  valid_genes <- gene_ids[!is.na(gene_ids)]
  cat("  Extracted", length(valid_genes), "valid gene IDs\n")
  
  # Build lookup table: key = short format (e.g., "08g02097"), value = full format
  lookup_table <- list()
  
  for (gene_id in valid_genes) {
    # Extract the short pattern from full gene ID
    # Vitvi05_01chr08g02097 -> 08g02097
    short_pattern <- stringr::str_extract(gene_id, "chr(\\d+g\\d+)$")
    if (!is.na(short_pattern)) {
      short_key <- stringr::str_remove(short_pattern, "chr")
      lookup_table[[short_key]] <- gene_id
    }
  }
  
  cat("  Built lookup table with", length(lookup_table), "entries\n")
  cat("  Sample entries:", paste(head(names(lookup_table), 5), collapse = ", "), "\n")
  
  return(lookup_table)
}

#' Convert GO gene IDs using lookup table
convert_with_lookup <- function(go_gene_ids, lookup_table) {
  cat("Converting", length(go_gene_ids), "GO gene IDs using lookup table...\n")
  
  results <- sapply(go_gene_ids, function(go_id) {
    # Extract pattern from GO ID: Vitvi08g02097 -> 08g02097
    pattern <- stringr::str_extract(go_id, "(\\d+g\\d+)$")
    if (!is.na(pattern) && pattern %in% names(lookup_table)) {
      return(lookup_table[[pattern]])
    }
    return(NA)
  })
  
  successful <- sum(!is.na(results))
  cat("  Successful conversions:", successful, "/", length(go_gene_ids), 
      sprintf("(%.1f%%)", 100 * successful / length(go_gene_ids)), "\n")
  
  return(results)
}

#' Process GMT file with optimized conversion
process_gmt_optimized <- function(input_file, output_file, lookup_table) {
  cat("Processing GMT file...\n")
  
  gmt_lines <- readLines(input_file)
  cat("  Total categories:", length(gmt_lines), "\n")
  
  # Collect all unique GO genes first
  all_go_genes <- c()
  for (line in gmt_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      all_go_genes <- c(all_go_genes, parts[3:length(parts)])
    }
  }
  
  unique_go_genes <- unique(all_go_genes)
  cat("  Unique GO genes to convert:", length(unique_go_genes), "\n")
  
  # Build conversion mapping for all unique genes
  conversion_map <- convert_with_lookup(unique_go_genes, lookup_table)
  names(conversion_map) <- unique_go_genes
  
  # Process each category
  converted_lines <- c()
  stats <- list(processed = 0, converted_genes = 0, dropped_genes = 0)
  
  for (i in seq_along(gmt_lines)) {
    if (i %% 1000 == 0) cat("  Processing category", i, "/", length(gmt_lines), "\n")
    
    parts <- strsplit(gmt_lines[i], "\t")[[1]]
    if (length(parts) >= 3) {
      category_id <- parts[1]
      category_name <- parts[2]
      gene_ids <- parts[3:length(parts)]
      
      # Convert genes using pre-built mapping
      converted_genes <- conversion_map[gene_ids]
      valid_genes <- converted_genes[!is.na(converted_genes)]
      
      if (length(valid_genes) > 0) {
        new_line <- paste(c(category_id, category_name, valid_genes), collapse = "\t")
        converted_lines <- c(converted_lines, new_line)
        
        stats$processed <- stats$processed + 1
        stats$converted_genes <- stats$converted_genes + length(valid_genes)
        stats$dropped_genes <- stats$dropped_genes + (length(gene_ids) - length(valid_genes))
      }
    }
  }
  
  # Write output
  writeLines(converted_lines, output_file)
  
  cat("Conversion complete:\n")
  cat("  Input categories:", length(gmt_lines), "\n")
  cat("  Output categories:", stats$processed, "\n")  
  cat("  Genes converted:", stats$converted_genes, "\n")
  cat("  Genes dropped:", stats$dropped_genes, "\n")
  
  return(stats)
}

# Main execution
cat("Starting optimized GO gene ID conversion...\n")

# Build lookup table from genome annotation
lookup_table <- build_gene_lookup_table(GENOME_ANNOTATION)

# Process GMT file
stats <- process_gmt_optimized(INPUT_GO, OUTPUT_GO, lookup_table)

cat("\n", rep("=", 60), "\n", sep = "")
cat("OPTIMIZED GO CONVERSION COMPLETED\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Final Results:\n")
cat("  Categories retained:", stats$processed, "\n")
cat("  Total genes converted:", stats$converted_genes, "\n")
cat("  Conversion rate:", sprintf("%.1f%%", 100 * stats$converted_genes / (stats$converted_genes + stats$dropped_genes)), "\n")
cat("  Output file:", OUTPUT_GO, "\n")

cat("\nSession Information:\n")
cat("Conversion completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")