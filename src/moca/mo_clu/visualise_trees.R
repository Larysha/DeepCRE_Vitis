# Visualize motif clustering trees from .nwk files
# This script creates publication-ready tree visualizations using Sandelin-Wassermann clustering
#
# Features:
# - Automatically finds and visualizes Sandelin-Wassermann phylogenetic trees
# - Works with forward strand representatives only (matching cluster_motifs.R output)
# - Generates metacluster-colored rectangular and circular tree layouts
# - Produces distance heatmaps showing motif similarity relationships
#
# Method: Sandelin-Wassermann (2004) similarity with hierarchical clustering
#
# Note: Strand-specific plots are NOT generated since all motifs are forward strand only
#
# Output files will be saved to: out/moca_results/mo_clu/trees/
######################

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)

# Set working directory to script location
script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
if (nchar(script_dir) > 0) {
  setwd(script_dir)
  cat("Working directory set to:", getwd(), "\n")
}

# Define output directory for tree visualizations
output_dir <- "../../../out/moca_results/mo_clu/trees"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Function to parse motif names and extract metadata
parse_motif_metadata <- function(motif_names) {
  metadata <- data.frame(
    motif_name = motif_names,
    species = NA,
    model = NA,
    pattern_id = NA,
    metacluster = NA,
    strand = NA,
    seqlet_count = NA,
    ic_score = NA,
    consensus = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(motif_names)) {
    name <- motif_names[i]
    
    # Parse naming convention: epm_species_model_pattern_metacluster_strand_seqlets_ic_consensus
    parts <- strsplit(name, "_")[[1]]
    
    if (length(parts) >= 8) {
      metadata$species[i] <- parts[2]
      metadata$model[i] <- parts[3]
      metadata$pattern_id[i] <- parts[4]
      metadata$metacluster[i] <- parts[5]
      metadata$strand[i] <- parts[6]
      metadata$seqlet_count[i] <- as.numeric(parts[7])
      metadata$ic_score[i] <- as.numeric(parts[8])
      
      # Consensus is everything after the 8th underscore
      if (length(parts) > 8) {
        metadata$consensus[i] <- parts[9]
      }
    }
  }
  
  # Create combined labels for visualization
  metadata$display_label <- paste0(
    "P", metadata$pattern_id, "_", 
    "MC", metadata$metacluster, "_", 
    metadata$strand, "\n",
    "n=", metadata$seqlet_count, " ",
    "IC=", metadata$ic_score
  )
  
  metadata$short_label <- paste0(
    "P", metadata$pattern_id, "_", 
    "MC", metadata$metacluster, "_", 
    metadata$strand
  )
  
  return(metadata)
}

# Function to create tree visualization
visualize_motif_tree <- function(tree_file, output_prefix = "motif_tree") {
  
  cat("Reading tree from:", tree_file, "\n")
  
  # Read the tree
  tree <- read.tree(tree_file)
  
  # Parse motif metadata
  metadata <- parse_motif_metadata(tree$tip.label)

  # Create color scheme for metaclusters
  metacluster_colors <- c("0" = "#804a5f", "1" = "#478fca")  # Dark Purple for MC0, Blue for MC1

  # Note: All motifs are forward strand only (from cluster_motifs.R filtering)
  # No need to filter further or create strand-specific plots
  cat("All motifs are forward strand representatives\n")
  cat("(Clustering script already filtered to forward strands only)\n")

  # Tree colored by metacluster (rectangular layout)
  p1 <- ggtree(tree, layout = "rectangular") +
    geom_tippoint(aes(color = metadata$metacluster[match(label, metadata$motif_name)]),
                  size = 4) +
    geom_tiplab(aes(label = metadata$short_label[match(label, metadata$motif_name)]),
                size = 4, hjust = -0.1, fontface = "bold") +
    scale_color_manual(values = metacluster_colors, name = "Metacluster",
                       labels = c("0" = "MC0 (High Expression)", "1" = "MC1 (Low Expression)")) +
    theme_tree2() +
    ggtitle("Motif Clustering Tree - Colored by Metacluster\n(Forward strand representatives only)") +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))

  # Circular tree with annotations
  p2 <- ggtree(tree, layout = "circular") +
    geom_tippoint(aes(color = metadata$metacluster[match(label, metadata$motif_name)]),
                  size = 3) +
    geom_tiplab(aes(label = metadata$short_label[match(label, metadata$motif_name)]),
                size = 3, hjust = -0.1, fontface = "bold") +
    scale_color_manual(values = metacluster_colors, name = "Metacluster",
                       labels = c("0" = "MC0 (High Expression)", "1" = "MC1 (Low Expression)")) +
    ggtitle("Circular Motif Tree - Colored by Metacluster\n(Forward strand representatives only)") +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))

  # Save plots to output directory
  mc_file <- file.path(output_dir, paste0(output_prefix, "_metacluster.png"))
  circular_file <- file.path(output_dir, paste0(output_prefix, "_circular.png"))

  ggsave(mc_file, p1, width = 15, height = 18, dpi = 300)
  ggsave(circular_file, p2, width = 14, height = 14, dpi = 300)

  cat("Saved tree visualizations:\n")
  cat("  ", mc_file, "\n")
  cat("  ", circular_file, "\n")
  
  # Print summary statistics
  cat("\nTree Summary:\n")
  cat("Total motifs:", length(tree$tip.label), "\n")
  cat("Metacluster distribution:\n")
  print(table(metadata$metacluster))
  cat("Strand (all should be F):\n")
  print(table(metadata$strand))
  cat("Seqlet count range:", min(metadata$seqlet_count, na.rm = TRUE), "-",
      max(metadata$seqlet_count, na.rm = TRUE), "\n")

  return(list(tree = tree, metadata = metadata, plots = list(p1, p2)))
}

# Function to create a heatmap of pairwise distances
create_distance_heatmap <- function(tree_file, output_prefix = "motif_distance") {
  
  # Read tree and calculate distances
  tree <- read.tree(tree_file)
  dist_matrix <- cophenetic(tree)
  metadata <- parse_motif_metadata(tree$tip.label)
  
  # Convert to long format for ggplot
  dist_df <- as.data.frame(as.table(dist_matrix))
  colnames(dist_df) <- c("Motif1", "Motif2", "Distance")
  
  # Add metadata
  dist_df$MC1 <- metadata$metacluster[match(dist_df$Motif1, metadata$motif_name)]
  dist_df$MC2 <- metadata$metacluster[match(dist_df$Motif2, metadata$motif_name)]
  
  # Create heatmap with larger text and better labels
  p_heatmap <- ggplot(dist_df, aes(x = Motif1, y = Motif2, fill = Distance)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Distance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, face = "bold"),
          axis.text.y = element_text(size = 8, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold")) +
    ggtitle("Motif Pairwise Distance Heatmap\n(Darker = More Similar, Lighter = More Different)") +
    labs(x = "Motifs", y = "Motifs")
  
  heatmap_file <- file.path(output_dir, paste0(output_prefix, "_heatmap.png"))
  ggsave(heatmap_file, p_heatmap, width = 16, height = 14, dpi = 300)
  cat("Saved distance heatmap:", heatmap_file, "\n")
  
  return(p_heatmap)
}

# Main execution
if (!interactive()) {
  # Command line usage
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    tree_file <- args[1]
    output_prefix <- if(length(args) > 1) args[2] else "motif_tree_visualization"
  } else {
    stop("Usage: Rscript visualize_motif_tree.R <tree_file.nwk> [output_prefix]")
  }
} else {
  # Interactive usage - automatically find tree file
  cat("=== Searching for Sandelin-Wassermann tree file ===\n")

  # Pattern to match forward_only tree files
  tree_pattern <- "*_cwm_clustering_forward_only_tree.nwk"

  # Try multiple possible directories
  possible_dirs <- c(
    ".",
    "../../../out/moca_results/mo_clu",
    "../../out/moca_results/mo_clu",
    "/home/rish/phd_2025/deepcre_vitis/vitis_cre/out/moca_results/mo_clu"
  )

  tree_files <- c()

  for (dir in possible_dirs) {
    if (dir.exists(dir)) {
      files <- list.files(dir, pattern = glob2rx(tree_pattern), full.names = TRUE)
      tree_files <- c(tree_files, files)
    }
  }

  if (length(tree_files) == 0) {
    stop("No tree files found matching pattern: ", tree_pattern,
         "\nSearched in directories: ", paste(possible_dirs, collapse = ", "))
  }

  # Use the most recent file (alphabetically last due to date prefix)
  tree_file <- sort(tree_files, decreasing = TRUE)[1]
  cat("Found tree file:", basename(tree_file), "\n")

  output_prefix <- "vitis_ssr_motif_tree_forward_only"
}

# Check if file exists
if (!file.exists(tree_file)) {
  stop("Tree file not found: ", tree_file)
}

###########################
# Generate visualizations
###########################
cat("\n=== Creating Motif Tree Visualizations ===\n")
cat("Method: Sandelin-Wassermann (2004)\n")
cat("Creating phylogenetic tree plots...\n")
results <- visualize_motif_tree(tree_file, output_prefix)

cat("Creating distance heatmap...\n")
heatmap_plot <- create_distance_heatmap(tree_file, output_prefix)

cat("\n=== Visualization Complete ===\n")
cat("All plots saved to:", output_dir, "\n\n")

# Print summary of generated files
cat("Generated files:\n")
cat("  - Metacluster tree (rectangular layout)\n")
cat("  - Circular tree\n")
cat("  - Distance heatmap\n")
cat("\nInterpretation tips:\n")
cat("  - Branch lengths represent motif dissimilarity (longer = more different)\n")
cat("  - Tight clusters indicate groups of similar motifs (potential TF families)\n")
cat("  - Metacluster colors show MC0 (high expression) vs MC1 (low expression)\n")
cat("  - Check consensus sequences for motifs that cluster together\n")
cat("\nAll visualizations complete!\n")