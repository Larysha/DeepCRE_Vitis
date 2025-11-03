# Cluster motifs based on Sandelin-Wassermann 2004 similarity metrics
# This script clusters motifs from JASPAR format files using Sandelin-Wassermann similarity
# and identifies groups of highly similar motifs for downstream analysis.
# Extended from original mo_cluster_v2.3.R
#
# IMPORTANT: This script filters to FORWARD STRAND MOTIFS ONLY to eliminate F/R redundancy.
# The Sandelin-Wassermann algorithm considers reverse complements during comparison,
# so this filtering does not lose biological information.
#
# Method: Sandelin-Wassermann (2004) similarity with hierarchical clustering
#
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre/src/moca")
######################

library(dplyr)
library(purrr)
library(universalmotif)
library(ape)

# Set working directory to script location
script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
if (nchar(script_dir) > 0) {
  setwd(script_dir)
  cat("Working directory set to:", getwd(), "\n")
}

# Source utility functions
source("../utils.R")

args <- commandArgs(trailingOnly = TRUE)

# Defaults - these should match the output from moca/mo_nom/get_matrices.R 
default_spec <- "vitis"
default_model <- "ssr"
default_input_dir <- "../../../out/moca_results/mo_nom"
default_output_dir <- "../../../out/moca_results/mo_clu"

SPEC <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_spec
MODEL <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_model
DATE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else format(Sys.Date(), "%Y%m%d")

MATRIX_TYPE <- "cwm"

if (length(args) >= 4 && nzchar(args[4])) {
  input_file <- args[4]
} else {
  # input file based on naming convention from get_matrices.R (mo_nom)
  # e.g., mo_nom/rdf5_epmvitisssr_cwm_scaled-motifs.jaspar (the seqlet normalised version)
  # Prefer scaled motifs for clustering as they have normalized information content
  pattern_scaled <- paste0("rdf5_.*", SPEC, MODEL, "_cwm_scaled-motifs\\.jaspar$")
  pattern_general <- paste0("rdf5_.*", SPEC, MODEL, "_cwm.*-motifs\\.jaspar$")

  # First try to find scaled motifs specifically
  jaspar_files <- list.files(default_input_dir, pattern = pattern_scaled, full.names = TRUE)

  if (length(jaspar_files) == 0) {
    # Fall back to any CWM motifs if scaled not found
    jaspar_files <- list.files(default_input_dir, pattern = pattern_general, full.names = TRUE)

    if (length(jaspar_files) == 0) {
      stop("No CWM JASPAR files found matching pattern: ", pattern_general,
           " in directory: ", default_input_dir)
    }
    if (length(jaspar_files) > 1) {
      cat("Multiple CWM JASPAR files found:\n")
      cat(paste(basename(jaspar_files), collapse = "\n"), "\n")
      stop("Multiple files match the pattern. Please specify one explicitly.")
    }
  }
  input_file <- jaspar_files[1]
}

dirpath_out <- if (length(args) >= 5 && nzchar(args[5])) args[5] else default_output_dir

# Validate inputs
if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

# Ensure output directory exists
if (!dir.exists(dirpath_out)) {
  dir.create(dirpath_out, recursive = TRUE)
  cat("Created output directory:", dirpath_out, "\n")
}

cat("Clustering CWM motifs from:", basename(input_file), "\n")
cat("Species:", SPEC, "| Model:", MODEL, "\n")
cat("Output directory:", dirpath_out, "\n")

#######################################################
# Load and validate motif data
#######################################################

load_and_validate_motifs <- function(file_path) {
  # Load JASPAR format CWM motifs (need to use the sequence patterns to cluster)
  # Returns a list of universalmotif objects with metadata

  cat("Loading CWM motifs from:", basename(file_path), "\n")

  tryCatch({
    motifs <- read_jaspar(file_path) # Using universalmotif package for JASPAR format
  }, error = function(e) {
    stop("Failed to read JASPAR file: ", e$message)
  })
  
  if (length(motifs) == 0) {
    stop("No motifs found in input file")
  }
  
  cat("Loaded", length(motifs), "motifs\n")
  
  # Validate motif structure
  for (i in seq_along(motifs)) {
    if (is.null(attr(motifs[[i]], "name"))) {
      stop("Motif ", i, " is missing a name attribute")
    }
    
    # Check matrix dimensions
    motif_matrix <- motifs[[i]]@motif
    if (nrow(motif_matrix) != 4) {
      stop("Motif ", attr(motifs[[i]], "name"), " does not have 4 rows (A,C,G,T)")
    }
    if (ncol(motif_matrix) < 3) {
      warning("Motif ", attr(motifs[[i]], "name"), " has fewer than 3 positions")
    }
  }
  
  return(motifs)
}

extract_motif_metadata <- function(motifs) {
  # Extract metadata from motif names and add computed statistics
  # Expects names in format: epm_<species>_<model>_<pattern>_<metacluster>_<strand>_<seqletcount>
  
  metadata <- map_dfr(motifs, function(motif) {
    motif_name <- attr(motif, "name")
    
    # Parse the systematic naming convention
    name_parts <- strsplit(motif_name, "_")[[1]]
    
    if (length(name_parts) >= 7) {
      # Standard format: epm_species_model_pattern_metacluster_strand_seqletcount
      # Extract from the end backwards
      seqlet_count <- as.numeric(name_parts[length(name_parts)])
      strand <- name_parts[length(name_parts) - 1]
      metacluster <- name_parts[length(name_parts) - 2]

      species <- name_parts[2]
      model <- name_parts[3]

      # Pattern ID is everything between model and metacluster
      # This handles cases where pattern_id might contain underscores
      if (length(name_parts) > 7) {
        pattern_id <- paste(name_parts[4:(length(name_parts) - 3)], collapse = "_")
      } else {
        pattern_id <- name_parts[4]
      }
    } else if (length(name_parts) == 6) {
      # Older format without seqlet count: epm_species_model_pattern_metacluster_strand
      seqlet_count <- NA
      species <- name_parts[2]
      model <- name_parts[3]
      pattern_id <- name_parts[4]
      metacluster <- name_parts[5]
      strand <- name_parts[6]
    } else {
      # Fallback for non-standard naming
      warning("Motif name does not follow expected convention: ", motif_name)
      seqlet_count <- NA
      species <- NA
      model <- NA
      pattern_id <- motif_name
      metacluster <- NA
      strand <- NA
    }
    
    # Compute motif statistics
    total_ic <- attr(motif, "icscore")
    consensus <- attr(motif, "consensus")
    
    tibble(
      motif_name = motif_name,
      species = species,
      model = model,
      pattern_id = pattern_id,
      metacluster = metacluster,
      strand = strand,
      seqlet_count = seqlet_count,
      total_ic = total_ic,
      consensus = consensus,
      motif_length = ncol(motif@motif)
    )
  })
  
  return(metadata)
}

#######################################################
# Enhanced motif preprocessing
#######################################################

enhance_motif_metadata <- function(motifs, metadata) {
  # Add computed statistics to motif objects and update names with IC and consensus
  # This follows the pattern from original script
  
  enhanced_motifs <- motifs
  
  for (i in seq_along(enhanced_motifs)) {
    motif_data <- metadata[i, ]
    
    # Add seqlet count as nsites (following original approach)
    if (!is.na(motif_data$seqlet_count)) {
      enhanced_motifs[[i]]["nsites"] <- motif_data$seqlet_count
    }
    
    # Create name with IC and consensus information
    if (!is.na(motif_data$total_ic) && !is.na(motif_data$consensus)) {
      total_ic_rounded <- round(motif_data$total_ic, 1)
      enhanced_name <- paste0(motif_data$motif_name, "_", 
                             total_ic_rounded, "_", 
                             motif_data$consensus)
      attr(enhanced_motifs[[i]], "name") <- enhanced_name
    }
  }
  
  return(enhanced_motifs)
}

#######################################################
# Motif comparison and similarity analysis
#######################################################

compute_similarity_matrix <- function(motifs, method = "SW") {
  # Compute pairwise similarity matrix using Sandelin-Wassermann method
  # Returns both the similarity matrix and distance matrix for clustering
  
  cat("Computing pairwise similarities using", method, "method...\n")
  cat("This may take a few minutes for large motif sets...\n")
  
  tryCatch({
    if (method == "SW") {
      # Sandelin-Wassermann method 
      similarity_matrix <- compare_motifs(motifs, method = "SW")
    } else {
      stop("Unsupported similarity method: ", method)
    }
  }, error = function(e) {
    stop("Failed to compute similarity matrix: ", e$message)
  })
  
  # Convert to data frame for easier manipulation
  sim_df <- as.data.frame(similarity_matrix)
  sim_matrix <- as.matrix(sim_df)
  
  # Create distance matrix (1 - similarity) for clustering
  dist_matrix <- 1 - sim_matrix
  
  return(list(
    similarity = sim_df,
    similarity_matrix = sim_matrix,
    distance = as.dist(dist_matrix)
  ))
}

identify_similarity_groups <- function(similarity_matrix, threshold_percentile = 0.95) {
  # Identify groups of highly similar motifs based on similarity threshold
  # Uses percentile-based thresholding as in original script
  # *** LOGIC EDIT *** — original script had diagonal handling where sum == 2
  # This version: unique motifs = zero highly similar partners (excluding self)
 
  sim_values <- as.vector(similarity_matrix)
  threshold <- quantile(sim_values, probs = threshold_percentile, na.rm = TRUE)
 
  cat("Using similarity threshold of", round(threshold, 4),
      "(", threshold_percentile * 100, "th percentile)\n")
 
  # Create binary matrix indicating high similarity
  high_sim_matrix <- similarity_matrix >= threshold
 
  # *** KEY CHANGE *** — explicitly exclude self-similarity from counts
  # since self-similarity is generally uninformative
  diag(high_sim_matrix) <- FALSE  # Remove diagonal before counting
  similarity_counts <- rowSums(high_sim_matrix)
 
  # Classify motifs based on number of similar partners
  # *** LOGIC *** — unique motifs have zero similar partners
  # This asks: "how many OTHER motifs is this similar to?"
  unique_motifs <- names(similarity_counts)[similarity_counts == 0]
  redundant_motifs <- names(similarity_counts)[similarity_counts > 0]
 
  cat("Found", length(unique_motifs), "motifs without highly similar counterparts\n")
  cat("Found", length(redundant_motifs), "motifs with highly similar counterparts\n")
 
  return(list(
    threshold = threshold,
    unique_motifs = unique_motifs,
    redundant_motifs = redundant_motifs,
    similarity_counts = similarity_counts
  ))
}

#######################################################
# Clustering and phylogenetic tree construction
#######################################################

build_phylogenetic_tree <- function(distance_matrix, method = "complete") {
  # Build phylogenetic tree from distance matrix
  # Handles edge length adjustment to ensure positive values
  
  cat("Building phylogenetic tree using", method, "linkage...\n")
  
  # Perform hierarchical clustering
  hc_result <- hclust(distance_matrix, method = method)
  
  # Convert to phylo object
  phylo_tree <- as.phylo(hc_result)
  
  # Ensure positive edge lengths 
  if (any(phylo_tree$edge.length <= 0)) {
    cat("Adjusting negative or zero edge lengths...\n")
    phylo_tree$edge.length <- phylo_tree$edge.length + 1
  }
  
  return(phylo_tree)
}


#######################################################
# Output generation and file writing
#######################################################

generate_output_files <- function(motifs, metadata, similarity_results, tree,
                                 output_dir, file_prefix) {
  # Generate comprehensive output files for the Sandelin-Wassermann clustering analysis

  # Create summary of motif characteristics
  motif_summary <- summarise_motifs(motifs)
  summary_file <- file.path(output_dir, paste0(file_prefix, "_motif_summary.csv"))
  write.csv(motif_summary, file = summary_file, row.names = FALSE)
  cat("Motif summary written to:", basename(summary_file), "\n")

  # Write similarity matrix
  sim_file <- file.path(output_dir, paste0(file_prefix, "_similarity_matrix.csv"))
  write.table(similarity_results$similarity, file = sim_file,
              sep = "\t", col.names = NA, quote = FALSE)
  cat("Similarity matrix written to:", basename(sim_file), "\n")

  # Write motif groupings
  unique_file <- file.path(output_dir, paste0(file_prefix, "_unique_motifs.csv"))
  redundant_file <- file.path(output_dir, paste0(file_prefix, "_redundant_motifs.csv"))

  write.table(similarity_results$groups$unique_motifs, file = unique_file,
              sep = "\t", col.names = "motif_name", quote = FALSE, row.names = FALSE)
  write.table(similarity_results$groups$redundant_motifs, file = redundant_file,
              sep = "\t", col.names = "motif_name", quote = FALSE, row.names = FALSE)

  cat("Unique motifs list written to:", basename(unique_file), "\n")
  cat("Redundant motifs list written to:", basename(redundant_file), "\n")

  # Write phylogenetic tree
  tree_file <- file.path(output_dir, paste0(file_prefix, "_tree.nwk"))
  write.tree(tree, file = tree_file)
  cat("Phylogenetic tree written to:", basename(tree_file), "\n")
  
  # Create comprehensive analysis report
  report_file <- file.path(output_dir, paste0(file_prefix, "_clustering_report.txt"))
  
  report_content <- c(
    paste("Motif Clustering Analysis Report"),
    paste("Method: Sandelin-Wassermann (2004)"),
    paste("Generated:", Sys.time()),
    paste("Input file:", basename(input_file)),
    paste("Species:", SPEC, "| Model:", MODEL, "| Matrix type:", toupper(MATRIX_TYPE)),
    "",
    paste("ANALYSIS SCOPE: Forward strand motifs only"),
    paste("Reverse strand motifs excluded to eliminate F/R pair redundancy"),
    paste("Sandelin-Wassermann algorithm considers reverse complements during comparison"),
    "",
    paste("CLUSTERING PARAMETERS:"),
    paste("  Method: Sandelin-Wassermann similarity"),
    paste("  Linkage: Complete linkage hierarchical clustering"),
    paste("  Similarity threshold: 95th percentile"),
    "",
    paste("RESULTS:"),
    paste("  Total motifs analysed:", length(motifs)),
    paste("  Similarity threshold (95th percentile):", round(similarity_results$groups$threshold, 4)),
    paste("  Unique motifs:", length(similarity_results$groups$unique_motifs)),
    paste("  Redundant motifs:", length(similarity_results$groups$redundant_motifs)),
    "",
    "Motif length distribution:",
    paste(capture.output(summary(metadata$motif_length)), collapse = "\n"),
    "",
    "Seqlet count distribution:",
    paste(capture.output(summary(metadata$seqlet_count, na.rm = TRUE)), collapse = "\n")
  )
  
  writeLines(report_content, report_file)
  cat("Analysis report written to:", basename(report_file), "\n")
  
  return(list(
    summary_file = summary_file,
    similarity_file = sim_file,
    unique_file = unique_file, 
    redundant_file = redundant_file,
    report_file = report_file
  ))
}

#######################################################
# Main execution pipeline
#######################################################

main <- function() {

  start_time <- Sys.time()

  tryCatch({
    # Load and validate motifs
    motifs <- load_and_validate_motifs(input_file)

    # Extract metadata
    metadata <- extract_motif_metadata(motifs)

    # Filter to forward strand motifs only
    # This eliminates F/R redundancy since SW algorithm considers reverse complements
    cat("\n=== Filtering to Forward Strand Motifs ===\n")
    cat("Total motifs loaded:", length(motifs), "\n")

    forward_indices <- which(metadata$strand == "F")
    reverse_count <- sum(metadata$strand == "R", na.rm = TRUE)

    if (length(forward_indices) == 0) {
      stop("No forward strand motifs found. Check motif naming convention.")
    }

    motifs <- motifs[forward_indices]
    metadata <- metadata[forward_indices, ]

    cat("Forward strand motifs retained:", length(motifs), "\n")
    cat("Reverse strand motifs excluded:", reverse_count, "\n")
    cat("This eliminates F/R pair redundancy from clustering analysis.\n")
    cat("Note: Sandelin-Wassermann will still consider reverse complements during comparison.\n")
    
    # Enhance motifs with additional metadata
    enhanced_motifs <- enhance_motif_metadata(motifs, metadata)
    
    # Compute similarity matrix
    similarity_results <- compute_similarity_matrix(enhanced_motifs, method = "SW")
    
    # Identify similarity groups
    similarity_groups <- identify_similarity_groups(similarity_results$similarity_matrix)
    similarity_results$groups <- similarity_groups

    # Build phylogenetic tree using Sandelin-Wassermann
    cat("\nBuilding phylogenetic tree...\n")
    phylo_tree <- build_phylogenetic_tree(similarity_results$distance)

    # Generate output files
    file_prefix <- paste0(DATE, "_", SPEC, "_", MODEL, "_cwm_clustering_forward_only")
    output_files <- generate_output_files(enhanced_motifs, metadata, similarity_results,
                                         phylo_tree, dirpath_out, file_prefix)
    
    # Print summary
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))
    
    cat("\n" , rep("=", 50), "\n")
    cat("CLUSTERING ANALYSIS COMPLETE\n")
    cat(rep("=", 50), "\n")
    cat("Runtime:", round(runtime, 2), "minutes\n")
    cat("Output directory:", dirpath_out, "\n")
    cat("Files generated:", length(output_files), "\n")
    
    return(invisible(list(
      motifs = enhanced_motifs,
      metadata = metadata,
      similarity = similarity_results,
      tree = phylo_tree,
      output_files = output_files
    )))
    
  }, error = function(e) {
    cat("ERROR: Analysis failed with message:", e$message, "\n")
    cat("Check input files and parameters\n")
    stop(e)
  })
}


run_main <- TRUE  # Set to FALSE to skip main() in interactive /CL mode
if (!interactive() || run_main) {
  main()
}