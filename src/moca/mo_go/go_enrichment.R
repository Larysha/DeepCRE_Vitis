#!/usr/bin/env Rscript
# GO Term Enrichment Analysis for Motif Performance Results
# This script performs functional enrichment analysis on top-performing motifs
# Extracted from motif_performance.R for modular analysis
######################

library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(grid)
library(circlize)
library(tibble)
library(RColorBrewer)

# Command line arguments with defaults
args <- commandArgs(trailingOnly = TRUE)

# Default file paths - find files by pattern and use mo_go output directory
# Find most recent motif performance files
perf_pattern <- "../../../out/moca_results/mo_proj/performance/*_vitis_motif_performance_motif_performance.csv"
perf_files <- Sys.glob(perf_pattern)
default_motif_performance <- if (length(perf_files) > 0) sort(perf_files, decreasing = TRUE)[1] else perf_pattern

# Find most recent integrated data files
integrated_pattern <- "../../../out/moca_results/mo_proj/*_vitis_motif_enrichment_integrated_data.csv"
integrated_files <- Sys.glob(integrated_pattern)
default_integrated_data <- if (length(integrated_files) > 0) sort(integrated_files, decreasing = TRUE)[1] else integrated_pattern

default_go_file <- "../../../vitis_data/gene_ontology/V_vinifera_ont_converted.gmt"  # Use converted file
default_output_dir <- "../../../out/moca_results/mo_go"

# Parse arguments
MOTIF_PERFORMANCE_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_motif_performance
INTEGRATED_DATA_FILE <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_integrated_data
GO_FILE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_go_file
OUTPUT_DIR <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_output_dir

DATE_STAMP <- format(Sys.Date(), "%Y%m%d")
PROJECT_NAME <- "vitis_go_enrichment"

#' Parse GMT file hierarchy for 3-level structure
#' 
#' @param gmt_file Path to GMT file with hierarchical structure
#' @return Data frame with hierarchy information
parse_gmt_hierarchy <- function(gmt_file) {
  
  cat("Parsing GMT hierarchy...\n")
  
  gmt_lines <- readLines(gmt_file)
  hierarchy_data <- data.frame()
  
  for (line in gmt_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      category_id <- parts[1]
      category_name <- parts[2]
      
      # Parse hierarchical levels
      id_parts <- strsplit(category_id, "\\.")[[1]]
      level1_id <- id_parts[1]
      level2_id <- if(length(id_parts) >= 2) paste(id_parts[1:2], collapse = ".") else NA
      level3_id <- if(length(id_parts) >= 3) paste(id_parts[1:3], collapse = ".") else NA
      
      # Extract category names at each level
      name_parts <- strsplit(category_name, "\\.")[[1]]
      level1_name <- name_parts[1]
      level2_name <- if(length(name_parts) >= 2) name_parts[2] else NA
      level3_name <- if(length(name_parts) >= 3) name_parts[3] else NA
      
      hierarchy_row <- data.frame(
        category_id = category_id,
        category_name = category_name,
        level1_id = level1_id,
        level1_name = level1_name,
        level2_id = level2_id,
        level2_name = level2_name,
        level3_id = level3_id,
        level3_name = level3_name,
        hierarchy_level = length(id_parts),
        stringsAsFactors = FALSE
      )
      
      hierarchy_data <- rbind(hierarchy_data, hierarchy_row)
    }
  }
  
  cat("  Parsed", nrow(hierarchy_data), "hierarchical categories\n")
  cat("  Level 1 categories:", length(unique(hierarchy_data$level1_id[!is.na(hierarchy_data$level1_id)])), "\n")
  cat("  Level 2 categories:", length(unique(hierarchy_data$level2_id[!is.na(hierarchy_data$level2_id)])), "\n")
  cat("  Level 3 categories:", length(unique(hierarchy_data$level3_id[!is.na(hierarchy_data$level3_id)])), "\n")
  
  return(hierarchy_data)
}

#' Load GO annotations from GMT or standard format
#' 
#' @param go_file Path to GO annotation file
#' @return Data frame with gene-category associations
load_go_annotations <- function(go_file) {
  
  cat("Loading GO annotations...\n")
  
  go_data <- NULL
  if (nzchar(go_file) && file.exists(go_file)) {
    # Check file format - GMT vs standard table
    if (grepl("\\.gmt$", go_file, ignore.case = TRUE)) {
      # Parse GMT format
      cat("  Loading GMT format GO annotations...\n")
      gmt_lines <- readLines(go_file)
      
      # Parse each line into gene-category mappings
      go_list <- list()
      for (line in gmt_lines) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 3) {
          category_id <- parts[1]
          category_name <- parts[2]
          gene_ids <- parts[3:length(parts)]
          
          # Create data frame for this category
          category_df <- data.frame(
            gene_id = gene_ids,
            category_id = category_id,
            category_name = category_name,
            stringsAsFactors = FALSE
          )
          go_list[[category_id]] <- category_df
        }
      }
      
      # Combine all categories
      go_data <- do.call(rbind, go_list)
      rownames(go_data) <- NULL
      
    } else {
      # Standard table format
      go_data <- read.table(go_file, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE, quote = "")
    }
    cat("  Loaded", nrow(go_data), "gene-category associations\n")
  } else {
    stop("GO file not found: ", go_file)
  }
  
  return(go_data)
}

#' Load motif performance results and integrated data
#' 
#' @param performance_file Path to motif performance CSV
#' @param integrated_file Path to integrated data CSV
#' @return List with performance metrics and integrated data
load_analysis_results <- function(performance_file, integrated_file) {
  
  cat("Loading analysis results...\n")
  
  # Load motif performance metrics
  if (!file.exists(performance_file)) {
    stop("Motif performance file not found: ", performance_file)
  }
  performance_metrics <- read.csv(performance_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(performance_metrics), "motif performance records\n")
  
  # Load integrated data
  if (!file.exists(integrated_file)) {
    stop("Integrated data file not found: ", integrated_file)
  }
  integrated_data <- read.csv(integrated_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(integrated_data), "integrated motif-gene associations\n")
  
  return(list(
    performance = performance_metrics,
    integrated = integrated_data
  ))
}

#' Calculate GO term enrichment for top motifs
#' 
#' @param performance_metrics Motif performance data
#' @param integrated_data Integrated motif-gene data with GO annotations
#' @return Data frame with enrichment statistics
calculate_go_enrichment <- function(performance_metrics, integrated_data) {
  
  cat("Calculating GO term enrichment...\n")
  
  # Calculate motif enrichment in functional categories
  go_enrichment <- integrated_data %>%
    filter(!is.na(category_name)) %>%
    group_by(epm, category_name) %>%
    summarise(motif_count = n(), .groups = "drop") %>%
    group_by(category_name) %>%
    filter(sum(motif_count) >= 5) %>%  # Categories with sufficient data
    ungroup() %>%
    left_join(
      integrated_data %>% 
        group_by(epm) %>% 
        summarise(total_motif_count = n(), .groups = "drop"),
      by = "epm"
    ) %>%
    mutate(
      enrichment_ratio = motif_count / total_motif_count,
      category_short = ifelse(nchar(category_name) > 30, 
                             paste0(substr(category_name, 1, 27), "..."), 
                             category_name)
    )
  
  cat("  Calculated enrichment for", length(unique(go_enrichment$category_name)), "GO categories\n")
  
  return(go_enrichment)
}

#' Generate GO enrichment visualization
#' 
#' @param go_enrichment GO enrichment data
#' @param performance_metrics Motif performance data  
#' @param output_dir Output directory
#' @param strand_type Type of strand analysis: "forward" or "reverse"
#' @return ggplot object
generate_go_enrichment_plot <- function(go_enrichment, performance_metrics, output_dir, strand_type = "forward") {
  
  cat("Generating GO enrichment plot for", strand_type, "strand motifs...\n")
  
  # Filter motifs based on strand type
  if (strand_type == "forward") {
    all_motifs <- performance_metrics %>%
      filter(!grepl("_R_", epm)) %>%  # Keep only forward strand (_F_) motifs
      pull(epm)
  } else if (strand_type == "reverse") {
    all_motifs <- performance_metrics %>%
      filter(grepl("_R_", epm)) %>%  # Keep only reverse strand (_R_) motifs
      pull(epm)
  } else {
    stop("Invalid strand_type. Must be 'forward' or 'reverse'")
  }
  
  cat("  Including", length(all_motifs), strand_type, "strand motifs in analysis\n")
  
  # Filter and prepare data, excluding non-informative categories
  excluded_patterns <- c("not assigned", "not annotated")
  
  excluded_count <- go_enrichment %>%
    filter(epm %in% all_motifs) %>%
    filter(grepl("not assigned|not annotated", tolower(category_name)) | 
           grepl("not assigned|not annotated", tolower(category_short))) %>%
    summarise(total_excluded = n()) %>%
    pull(total_excluded)
  
  # Use ALL detailed GO categories without aggregation
  go_plot_data <- go_enrichment %>%
    filter(epm %in% all_motifs) %>%
    filter(!grepl("not assigned|not annotated", tolower(category_name)) & 
           !grepl("not assigned|not annotated", tolower(category_short))) %>%
    # Filter categories with moderately stringent data quality threshold
    group_by(category_name) %>%
    filter(
      n() >= length(all_motifs) * 0.5 &  # Must have data for at least 50% of motifs
      max(enrichment_ratio, na.rm = TRUE) >= 0.003  # At least one motif must have enrichment >= 0.002
    ) %>%
    ungroup()
  
  if (nrow(go_plot_data) > 0) {
    # Parse GMT hierarchy and add Level 1 information for color annotation
    hierarchy <- parse_gmt_hierarchy(GO_FILE)
    
    # Add hierarchy information to the detailed GO categories 
    go_detailed_with_hierarchy <- go_plot_data %>%
      left_join(hierarchy %>% 
                select(category_name, category_id, level1_name, level1_id),
                by = "category_name") %>%
      filter(!is.na(level1_name))  # Only keep categories with Level 1 hierarchy info
    
    cat("  Detailed categories prepared:", nrow(go_detailed_with_hierarchy), "entries\n")
    cat("  Unique GO categories:", length(unique(go_detailed_with_hierarchy$category_name)), "\n")
    cat("  Level 1 groups represented:", length(unique(go_detailed_with_hierarchy$level1_name)), "\n")
    
    if (nrow(go_detailed_with_hierarchy) > 0) {
      # Create data matrix with shortened category names (last 2 hierarchical levels)
      heatmap_data_clean <- go_detailed_with_hierarchy %>%
        select(epm, category_name, enrichment_ratio) %>%
        mutate(
          category_name_short = sapply(category_name, function(name) {
            if (grepl("\\.", name)) {
              parts <- strsplit(name, "\\.")[[1]]
              if (length(parts) >= 2) {
                paste(parts[(length(parts)-1):length(parts)], collapse = ".")
              } else {
                name
              }
            } else {
              name
            }
          })
        )
      
      heatmap_matrix <- heatmap_data_clean %>%
        pivot_wider(names_from = "epm", values_from = "enrichment_ratio") %>%
        column_to_rownames("category_name") %>%
        as.matrix()
      
      # Ensure matrix is numeric
      heatmap_matrix <- apply(heatmap_matrix, c(1,2), as.numeric)
      
      # Create mapping from full category names to shortened names for display
      category_name_mapping <- heatmap_data_clean %>%
        select(category_name, category_name_short) %>%
        distinct() %>%
        deframe()
      
      # Replace NA values with 0
      heatmap_matrix[is.na(heatmap_matrix)] <- 0
      
      # Create mapping from detailed category names to Level 1 for annotation
      category_to_level1 <- go_detailed_with_hierarchy %>%
        select(category_name, level1_name) %>%
        distinct() %>%
        deframe()
      
      # Order rows by Level 1 groupings, with metabolism categories first
      metabolism_level1 <- c("Amino acid metabolism", "Carbohydrate metabolism", 
                            "Nucleotide metabolism", "Lipid metabolism", "Coenzyme metabolism")
      
      # Get detailed categories for each Level 1 group
      categories_by_level1 <- split(names(category_to_level1), category_to_level1)
      
      # Order detailed categories: metabolism Level 1 first, then others alphabetically
      row_order <- c()
      
      # First add metabolism-related detailed categories
      for (l1_cat in metabolism_level1) {
        if (l1_cat %in% names(categories_by_level1)) {
          row_order <- c(row_order, sort(categories_by_level1[[l1_cat]]))
        }
      }
      
      # Then add other Level 1 categories
      other_level1 <- setdiff(names(categories_by_level1), metabolism_level1)
      for (l1_cat in sort(other_level1)) {
        row_order <- c(row_order, sort(categories_by_level1[[l1_cat]]))
      }
      
      # Keep only rows that exist in the matrix
      row_order <- row_order[row_order %in% rownames(heatmap_matrix)]
      heatmap_matrix <- heatmap_matrix[row_order, , drop = FALSE]
      
      # Store original rownames for Level1 mapping before applying shortened names
      original_rownames <- rownames(heatmap_matrix)
      rownames(heatmap_matrix) <- category_name_mapping[original_rownames]
      
      # Create color function with better scale
      max_val <- max(heatmap_matrix, na.rm = TRUE)
      # Handle case where max_val might be -Inf, NA, or non-numeric
      if (!is.finite(max_val) || max_val <= 0) {
        max_val <- 1  # Default fallback value
      }
      col_fun <- colorRamp2(c(0, max_val/3, 2*max_val/3, max_val), 
                           c("#ffffff", "#fce4ec", "#f8bbd9", "#c62828"))
      
      # Create row annotation for Level 1 groupings
      row_level1_mapping <- category_to_level1[original_rownames]
      names(row_level1_mapping) <- rownames(heatmap_matrix)
      
      # Create a gradient color palette for Level 1 categories
      unique_level1 <- unique(row_level1_mapping)
      
      # Define gradients for different biological process groups
      metabolism_colors <- c(
        "Amino acid metabolism" = "#671436",    
        "Carbohydrate metabolism" = "#7d1e45",    
        "Nucleotide metabolism" = "#891b49",     
        "Lipid metabolism" = "#b31f5d",          
        "Coenzyme metabolism" = "#c52266"        
      )
      
      # Define gradients for cellular processes
      cellular_colors <- c(
        "Cell cycle organisation" = "#052b67",     
        "Cell wall organisation" = "#083988",      
        "Cellular respiration" = "#1d58b8",        
        "Cytoskeleton organisation" = "#2e71dc",  
        "Vesicle trafficking" = "#5992ee"          
      )
      
      # Define gradients for protein processes  
      protein_colors <- c(
        "Protein biosynthesis" = "#4A148C",        
        "Protein homeostasis" = "#6A1B9A",         
        "Protein modification" = "#7B1FA2",        
        "Protein translocation" = "#8E24AA"        
      )
      
      # Define gradients for information processes
      info_colors <- c(
        "Chromatin organisation" = "#1d778b",      
        "DNA damage response" = "#2494ad",         
        "RNA biosynthesis" = "#43bad5",            
        "RNA processing" = "#73d9f0"               
      )
      
      # Create final color mapping
      level1_colors <- c(metabolism_colors, cellular_colors, protein_colors, info_colors)
      
      # Add any remaining categories with custom colors
      remaining_level1 <- setdiff(unique_level1, names(level1_colors))
      if (length(remaining_level1) > 0) {
        remaining_colors <- c("#ad597e", "#008066", "#8d77ab", "#808000", "#657b9e")
        names(remaining_colors) <- remaining_level1[1:min(length(remaining_level1), length(remaining_colors))]
        level1_colors <- c(level1_colors, remaining_colors)
      }
      
      # Only keep colors for categories that actually appear
      level1_colors <- level1_colors[names(level1_colors) %in% unique_level1]
      
      row_anno <- rowAnnotation(
        "Biological Process Category" = row_level1_mapping,
        col = list("Biological Process Category" = level1_colors),
        annotation_name_side = "top",
        show_legend = TRUE,
        show_annotation_name = FALSE,
        annotation_width = unit(1.2, "cm")
      )
      
      # Create heatmap
      ht <- Heatmap(
        heatmap_matrix,
        name = "Enrichment\nRatio",
        col = col_fun,
        
        # Row settings - move row names to left side
        cluster_rows = FALSE,
        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),  # Same size as column names
        row_title = NULL,  # Remove any row title
        
        # Column settings - make blocks square and tilt column names
        cluster_columns = FALSE,
        show_column_names = TRUE,
        column_names_side = "bottom",
        column_names_rot = 45,  # Tilt motif names for better readability
        column_names_gp = gpar(fontsize = 8),
        column_names_max_height = unit(4, "cm"),  # Add space for tilted names
        column_title = NULL,  # Remove any column title
        width = unit(ncol(heatmap_matrix) * 0.5, "cm"),  # Square blocks (adjust for square ratio)
        height = unit(nrow(heatmap_matrix) * 0.5, "cm"),  # Match width for square blocks
        
        # Add row annotation
        left_annotation = row_anno,
        
        # Heatmap appearance
        rect_gp = gpar(col = "white", lwd = 0.5),
        
        
        # Legend
        heatmap_legend_param = list(
          title = "Enrichment\nRatio",
          at = round(seq(0, max_val, length.out = 5), 3),  # Round to 3 decimal places
          labels = sprintf("%.3f", round(seq(0, max_val, length.out = 5), 3)),  # Format labels to 3 decimal places
          labels_gp = gpar(fontsize = 8),
          title_gp = gpar(fontsize = 9)
        )
      )
      
      # Create title and subtitle
      title_text <- paste("Gene Ontology (GO) Enrichment -", str_to_title(strand_type), "Strand Motifs")
      subtitle_text <- paste0("Filtering: ≥50% motif coverage, max enrichment ≥0.003, showing last 2 GO term levels\n",
                             excluded_count, " 'not assigned' categories excluded (genes lacking GO annotations)\n",
                             "Analysis of ", length(all_motifs), " ", strand_type, " strand motifs")
      
      # Save the heatmap
      pdf_file <- file.path(output_dir, paste0(DATE_STAMP, "_", PROJECT_NAME, "_go_enrichment_", strand_type, ".pdf"))
      pdf(pdf_file, width = 12 + 1/2.54, height = max(6, nrow(heatmap_matrix) * 0.3 + 2.5) - 1/2.54)  # Add 1cm width, remove 1cm height
      
      # Draw the heatmap with title only
      draw(ht, 
           column_title = title_text,
           column_title_gp = gpar(fontsize = 11, fontface = "bold"),
           padding = unit(c(2, 4, 3, 2), "mm"),  # Reduce padding - top, right, bottom, left
           gap = unit(2, "mm"))  # Smaller gap between heatmap and legends
      
      # Add subtitle with proper spacing from bottom
      grid.text(subtitle_text, 
                x = 0.5, 
                y = unit(2, "cm"),  # Move up 2 cm from bottom with padding
                just = "centre",
                gp = gpar(fontsize = 9, col = "gray30"))
      
      dev.off()
      
      cat("  Saved", strand_type, "strand GO enrichment heatmap with Level 1 gradient annotation\n")
      return(ht)
    } else {
      cat("  Warning: No detailed GO categories with sufficient data across", strand_type, "strand motifs\n")
      return(NULL)
    }
  } else {
    cat("  Warning: No data available for", strand_type, "strand GO enrichment plot\n")
    return(NULL)
  }
}

#' Add GO annotations to integrated data
#' 
#' @param integrated_data Integrated motif-gene data
#' @param go_data GO annotation data
#' @return Integrated data with GO annotations
add_go_annotations <- function(integrated_data, go_data) {
  
  cat("Adding GO annotations to integrated data...\n")
  
  # For GMT format, gene IDs should match directly (assuming converted file)
  if ("category_name" %in% colnames(go_data)) {
    # GMT format - direct gene ID matching
    annotated_data <- integrated_data %>%
      left_join(go_data %>% 
               select(gene_id, category_name, category_id), 
               by = "gene_id", 
               relationship = "many-to-many")
  } else {
    # Standard format - clean gene IDs for matching  
    go_clean <- go_data %>%
      mutate(
        gene_id = tolower(IDENTIFIER),
        gene_id = gsub("[-'_]+", "", gene_id),
        gene_id = gsub("\\..*", "", gene_id)
      )
    
    integrated_data$gene_id_clean <- tolower(integrated_data$gene_id)
    integrated_data$gene_id_clean <- gsub("[-'_]+", "", integrated_data$gene_id_clean)
    
    annotated_data <- integrated_data %>%
      left_join(go_clean %>% select(gene_id_clean = gene_id, NAME, BINCODE), 
               by = "gene_id_clean")
  }
  
  # Report annotation statistics
  genes_with_go <- sum(!is.na(annotated_data$category_name))
  cat("  Added GO annotations to", genes_with_go, "gene associations\n")
  
  return(annotated_data)
}

#' Calculate and display unique gene count for heatmap categories
#' 
#' @param integrated_with_go Integrated data with GO annotations
#' @param go_plot_result Result from generate_go_enrichment_plot (contains heatmap matrix info)
#' @return None (prints to console)
calculate_heatmap_gene_count <- function(integrated_with_go, go_enrichment, performance_metrics) {
  
  # Filter to forward strand motifs only
  all_motifs <- performance_metrics %>%
    filter(!grepl("_R_", epm)) %>%
    pull(epm)
  
  # Get categories that passed filtering (similar logic to plot function)
  filtered_categories <- go_enrichment %>%
    filter(epm %in% all_motifs) %>%
    filter(!grepl("not assigned|not annotated", tolower(category_name))) %>%
    group_by(category_name) %>%
    filter(
      n() >= length(all_motifs) * 0.5 &
      max(enrichment_ratio, na.rm = TRUE) >= 0.003
    ) %>%
    pull(category_name) %>%
    unique()
  
  # Calculate unique genes for these categories
  unique_genes_in_heatmap <- integrated_with_go %>%
    filter(category_name %in% filtered_categories, !is.na(gene_id)) %>%
    pull(gene_id) %>%
    n_distinct()
  
  cat("  Unique genes captured in heatmap:", unique_genes_in_heatmap, "\n")
  cat("  GO categories displayed:", length(filtered_categories), "\n")
}

# Main execution
cat("Starting GO enrichment analysis...\n")

# Load GO annotations
go_data <- load_go_annotations(GO_FILE)

# Load analysis results
analysis_results <- load_analysis_results(MOTIF_PERFORMANCE_FILE, INTEGRATED_DATA_FILE)

# Add GO annotations to integrated data
integrated_with_go <- add_go_annotations(analysis_results$integrated, go_data)

# Calculate GO enrichment
go_enrichment <- calculate_go_enrichment(analysis_results$performance, integrated_with_go)

# Generate GO enrichment plots for both forward and reverse strand motifs
go_plot_forward <- generate_go_enrichment_plot(go_enrichment, analysis_results$performance, OUTPUT_DIR, strand_type = "forward")
go_plot_reverse <- generate_go_enrichment_plot(go_enrichment, analysis_results$performance, OUTPUT_DIR, strand_type = "reverse")

# Calculate and display gene count for heatmap
calculate_heatmap_gene_count(integrated_with_go, go_enrichment, analysis_results$performance)

# Save enrichment results
output_base <- file.path(OUTPUT_DIR, paste0(DATE_STAMP, "_", PROJECT_NAME))
write.csv(go_enrichment, paste0(output_base, "_enrichment_results.csv"), row.names = FALSE)

# Final summary
cat("\n", rep("=", 60), "\n", sep = "")
cat("GO ENRICHMENT ANALYSIS COMPLETED\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Results:\n")
cat("  GO categories analyzed:", length(unique(go_enrichment$category_name)), "\n")
cat("  Genes with GO annotations:", sum(!is.na(integrated_with_go$category_name)), "\n")
cat("  Output files: enrichment results CSV and plot PDF\n")

cat("\nOutput Directory:", OUTPUT_DIR, "\n")
cat("Files Generated: GO enrichment analysis complete\n")

cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")