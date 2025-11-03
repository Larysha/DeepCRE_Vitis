# Cluster motifs based on Sandelin-Wassermann 2004 similarity metrics
# This script clusters PFM motifs from JASPAR format files using multiple similarity metrics
# and identifies groups of highly similar motifs for downstream analysis.
# Modified from cluster_motifs.R to work with PFM matrices
# DeepCRE/moca/mo_clu
######################

library(dplyr)
library(purrr)
library(grid)
library(TFBSTools)
library(motifStack)
library(universalmotif)
library(ape)
library(ggtree)

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
default_input_dir <- "../../out/moca_results/mo_nom"
default_output_dir <- "../../out/moca_results/mo_clu"

SPEC <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_spec
MODEL <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_model
DATE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else format(Sys.Date(), "%Y%m%d")

MATRIX_TYPE <- "pfm"

if (length(args) >= 4 && nzchar(args[4])) {
  input_file <- args[4]
} else {
  # input file based on naming convention from get_matrices.R (mo_nom)
  # e.g., ../../out/moca_results/mo_nom/rdf5_PFM_patternvitisssr_pfm-motifs.jaspar
  pattern <- paste0("rdf5_PFM_pattern.*", SPEC, MODEL, "_pfm-motifs\\.jaspar$")
  jaspar_files <- list.files(default_input_dir, pattern = pattern, full.names = TRUE)
  
  if (length(jaspar_files) == 0) {
    stop("No PFM JASPAR files found matching pattern: ", pattern,
         " in directory: ", default_input_dir)
  }
  if (length(jaspar_files) > 1) {
    cat("Multiple PFM JASPAR files found:\n")
    cat(paste(basename(jaspar_files), collapse = "\n"), "\n")
    stop("Multiple files match the pattern. Please specify one explicitly.")
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

cat("Clustering PFM motifs from:", basename(input_file), "\n")
cat("Species:", SPEC, "| Model:", MODEL, "\n")
cat("Output directory:", dirpath_out, "\n")

#######################################################
# Load and validate motif data
#######################################################

load_and_validate_motifs <- function(file_path) {
  # Load JASPAR format PFM motifs 
  # PFMs contain sequence frequency counts and should work directly with read_jaspar()
  # Returns a list of universalmotif objects with metadata

  cat("Loading PFM motifs from:", basename(file_path), "\n")

  tryCatch({
    motifs <- read_jaspar(file_path) # Using universalmotif package for JASPAR format
  }, error = function(e) {
    stop("Failed to read JASPAR file: ", e$message)
  })
  
  if (length(motifs) == 0) {
    stop("No motifs found in input file")
  }
  
  cat("Loaded", length(motifs), "motifs\n")
  
  # Basic filtering: remove very short motifs only
  valid_motifs <- list()
  removed_count <- 0
  
  for (i in seq_along(motifs)) {
    motif <- motifs[[i]]
    motif_name <- attr(motif, "name")
    motif_length <- ncol(motif@motif)
    
    # Only filter out extremely short motifs
    if (motif_length >= 3) {
      valid_motifs[[motif_name]] <- motif
    } else {
      removed_count <- removed_count + 1
      cat("Removed short motif:", motif_name, "(length:", motif_length, ")\n")
    }
  }
  
  if (removed_count > 0) {
    cat("Removed", removed_count, "motifs due to short length\n")
  }
  
  motifs <- valid_motifs
  
  if (length(motifs) == 0) {
    stop("No valid motifs remain after filtering")
  }
  
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
  }
  
  cat("Final motif count for clustering:", length(motifs), "\n")
  
  return(motifs)
}

extract_motif_metadata <- function(motifs) {
  # Extract metadata from motif names and add computed statistics
  # Expects names in format: epm_<species>_<model>_<pattern>_<metacluster>_<strand>_<seqletcount>
  
  metadata <- map_dfr(motifs, function(motif) {
    motif_name <- attr(motif, "name")
    
    # Parse the systematic naming convention
    name_parts <- strsplit(motif_name, "_")[[1]]
    
    if (length(name_parts) >= 6) {
      # Extract seqlet count from the name (last part)
      seqlet_count <- as.numeric(name_parts[length(name_parts)])
      
      species <- name_parts[2]
      model <- name_parts[3]
      pattern_id <- name_parts[4]
      metacluster <- name_parts[5]
      strand <- name_parts[6]
      
      # Handle cases where there might be additional parts in the name
      if (length(name_parts) > 6) {
        # Reconstruct pattern_id if it contains underscores
        pattern_id <- paste(name_parts[4:(length(name_parts)-2)], collapse = "_")
        metacluster <- name_parts[length(name_parts)-1]
        strand <- name_parts[length(name_parts)]
      }
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
  # Compute pairwise similarity matrix using specified method
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
  
  # Check for problematic values
  na_count <- sum(is.na(sim_matrix))
  inf_count <- sum(is.infinite(sim_matrix))
  
  if (na_count > 0) {
    cat("Warning: Found", na_count, "NA values in similarity matrix\n")
    cat("Replacing NA values with 0 (no similarity)...\n")
    sim_matrix[is.na(sim_matrix)] <- 0
    sim_df[is.na(sim_df)] <- 0
  }
  
  if (inf_count > 0) {
    cat("Warning: Found", inf_count, "infinite values in similarity matrix\n")
    cat("Replacing infinite values with 1 (maximum similarity)...\n")
    sim_matrix[is.infinite(sim_matrix) & sim_matrix > 0] <- 1
    sim_matrix[is.infinite(sim_matrix) & sim_matrix < 0] <- 0
    sim_df[is.infinite(sim_df) & sim_df > 0] <- 1
    sim_df[is.infinite(sim_df) & sim_df < 0] <- 0
  }
  
  # Create distance matrix (1 - similarity) for clustering
  dist_matrix <- 1 - sim_matrix
  
  # Ensure distance matrix values are valid
  dist_matrix[dist_matrix < 0] <- 0  # Similarity > 1 becomes distance 0
  dist_matrix[is.na(dist_matrix)] <- 1  # NA becomes maximum distance
  dist_matrix[is.infinite(dist_matrix)] <- 1  # Inf becomes maximum distance
  
  # Check diagonal should be 0 (self-distance)
  diag(dist_matrix) <- 0
  
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
  
  # Validate distance matrix before clustering
  if (any(is.na(distance_matrix))) {
    stop("Distance matrix contains NA values")
  }
  if (any(is.infinite(distance_matrix))) {
    stop("Distance matrix contains infinite values")
  }
  if (any(distance_matrix < 0)) {
    cat("Warning: Found negative distances, setting to 0...\n")
    distance_matrix[distance_matrix < 0] <- 0
  }
  
  # Check matrix properties
  cat("Distance matrix summary:\n")
  cat("  Range:", range(distance_matrix), "\n")
  cat("  NA count:", sum(is.na(distance_matrix)), "\n")
  cat("  Inf count:", sum(is.infinite(distance_matrix)), "\n")
  
  tryCatch({
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
    
  }, error = function(e) {
    cat("Clustering failed with error:", e$message, "\n")
    cat("Attempting alternative clustering method...\n")
    
    # Try with different method
    tryCatch({
      hc_result <- hclust(distance_matrix, method = "average")
      phylo_tree <- as.phylo(hc_result)
      
      if (any(phylo_tree$edge.length <= 0)) {
        phylo_tree$edge.length <- phylo_tree$edge.length + 1
      }
      
      cat("Successfully clustered using 'average' method\n")
      return(phylo_tree)
      
    }, error = function(e2) {
      cat("All clustering methods failed\n")
      return(NULL)
    })
  })
}

perform_motifstack_clustering <- function(motifs) {
  # Perform clustering using motifStack's clusterMotifs function
  # This provides an alternative clustering approach to compare with
  
  cat("Performing motifStack clustering (this may take several minutes)...\n")
  
  # Convert to pcm format for motifStack
  pcm_motifs <- convert_motifs(motifs, class = "motifStack-pcm")
  
  tryCatch({
    clustered_motifs <- clusterMotifs(pcm_motifs, method = "Smith-Waterman")
    phylo_tree <- as.phylo(clustered_motifs)
    
    return(list(
      clustered = clustered_motifs,
      tree = phylo_tree,
      ordered_motifs = pcm_motifs[clustered_motifs$order]
    ))
  }, error = function(e) {
    warning("motifStack clustering failed: ", e$message)
    return(NULL)
  })
}

#######################################################
# Output generation and file writing
#######################################################

generate_output_files <- function(motifs, metadata, similarity_results, trees, 
                                 output_dir, file_prefix) {
  # Generate comprehensive output files for the clustering analysis
  
  # Create summary of motif characteristics
  motif_summary <- summarise_motifs(motifs)
  summary_file <- file.path(output_dir, paste0(file_prefix, "_motif_summary.csv"))
  write.csv(motif_summary, file = summary_file, row.names = FALSE)
  cat("Motif summary written to:", basename(summary_file), "\n")
  
  # Write similarity matrix
  sim_file <- file.path(output_dir, paste0(file_prefix, "_similarity_matrix_SW.csv"))
  write.table(similarity_results$similarity, file = sim_file, 
              sep = "\t", col.names = NA, quote = FALSE)
  cat("Similarity matrix written to:", basename(sim_file), "\n")
  
  # Write motif groupings
  unique_file <- file.path(output_dir, paste0(file_prefix, "_unique_motifs_SW.csv"))
  redundant_file <- file.path(output_dir, paste0(file_prefix, "_redundant_motifs_SW.csv"))
  
  write.table(similarity_results$groups$unique_motifs, file = unique_file, 
              sep = "\t", col.names = "motif_name", quote = FALSE, row.names = FALSE)
  write.table(similarity_results$groups$redundant_motifs, file = redundant_file, 
              sep = "\t", col.names = "motif_name", quote = FALSE, row.names = FALSE)
  
  cat("Unique motifs list written to:", basename(unique_file), "\n")
  cat("Redundant motifs list written to:", basename(redundant_file), "\n")
  
  # Write phylogenetic trees
  if (!is.null(trees$sandelin_wassermann)) {
    sw_tree_file <- file.path(output_dir, paste0(file_prefix, "_tree_Sandelin_Wassermann.nwk"))
    write.tree(trees$sandelin_wassermann, file = sw_tree_file)
    cat("Sandelin-Wassermann tree written to:", basename(sw_tree_file), "\n")
  }
  
  if (!is.null(trees$motifstack)) {
    ms_tree_file <- file.path(output_dir, paste0(file_prefix, "_tree_Smith_Waterman.nwk"))
    write.tree(trees$motifstack, file = ms_tree_file)
    cat("Smith-Waterman tree written to:", basename(ms_tree_file), "\n")
  }
  
  # Create comprehensive analysis report
  report_file <- file.path(output_dir, paste0(file_prefix, "_clustering_report.txt"))
  
  report_content <- c(
    paste("Motif Clustering Analysis Report"),
    paste("Generated:", Sys.time()),
    paste("Input file:", basename(input_file)),
    paste("Species:", SPEC, "| Model:", MODEL, "| Matrix type:", toupper(MATRIX_TYPE)),
    "",
    paste("Total motifs analysed:", length(motifs)),
    paste("Similarity threshold (95th percentile):", round(similarity_results$groups$threshold, 4)),
    paste("Unique motifs:", length(similarity_results$groups$unique_motifs)),
    paste("Redundant motifs:", length(similarity_results$groups$redundant_motifs)),
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
    
    # Enhance motifs with additional metadata
    enhanced_motifs <- enhance_motif_metadata(motifs, metadata)
    
    # Compute similarity matrix
    similarity_results <- compute_similarity_matrix(enhanced_motifs, method = "SW")
    
    # Identify similarity groups
    similarity_groups <- identify_similarity_groups(similarity_results$similarity_matrix)
    similarity_results$groups <- similarity_groups
    
    # Build phylogenetic trees
    trees <- list()
    
    # Sandelin-Wassermann tree
    trees$sandelin_wassermann <- build_phylogenetic_tree(similarity_results$distance)
    
    # motifStack clustering (optional, can be slow)
    cat("\nPerforming alternative clustering with motifStack...\n")
    motifstack_result <- perform_motifstack_clustering(enhanced_motifs)
    if (!is.null(motifstack_result)) {
      trees$motifstack <- motifstack_result$tree
    }
    
    # Generate output files
    file_prefix <- paste0(DATE, "_", SPEC, "_", MODEL, "_pfm_clustering")
    output_files <- generate_output_files(enhanced_motifs, metadata, similarity_results, 
                                         trees, dirpath_out, file_prefix)
    
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
      trees = trees,
      output_files = output_files
    )))
    
  }, error = function(e) {
    cat("ERROR: Analysis failed with message:", e$message, "\n")
    cat("Check input files and parameters\n")
    stop(e)
  })
}

run_main <- TRUE  # Set to FALSE to skip main() in interactive/CL mode
if (!interactive() || run_main) {
  main()
}