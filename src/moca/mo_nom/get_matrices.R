# This script processes MoDisco HDF5 files to extract motif matrices (PFM, CWM) in JASPAR format
# and applies a systematic nomenclature for easy identification.
# Also generates summary statistics for CWM matrices of contribution scores
# handles PFMs from sequence field with adjusted background frequencies
#
# remember to setwd to project root:
# setwd("/home/rish/phd_2025/deepcre_vitis/vitis_cre")
######################

library(dplyr)
library(purrr)
library(rhdf5)

args <- commandArgs(trailingOnly = TRUE)

# Defaults
default_spec <- "vitis"
default_model <- "ssr"
default_input_dir <- "out/modisco/"
default_output_dir <- "out/moca_results/mo_nom"

# Assign positional arguments or use defaults
SPEC <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_spec
MODEL <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_model
DATE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else format(Sys.Date(), "%Y%m%d")

# Optional input file
if (length(args) >= 4 && nzchar(args[4])) {
  FILE1 <- args[4]
} else {
  modisco_files <- list.files(default_input_dir, pattern = ".*motifs.*\\.hdf5", full.names = TRUE)
  if (length(modisco_files) == 0) stop("No files matching pattern '.*motifs.*\\.hdf5' found in ", default_input_dir)
  if (length(modisco_files) > 1) stop("Multiple *motifs*.hdf5 files found. Please specify one explicitly.")
  FILE1 <- modisco_files[1]
}

dirpath_out <- if (length(args) >= 5 && nzchar(args[5])) args[5] else default_output_dir

# Ensure output directories exist
if (!dir.exists(dirpath_out)) dir.create(dirpath_out, recursive = TRUE)
if (!dir.exists("out/moca_results")) dir.create("out/moca_results", recursive = TRUE)

cat("Processing MoDisco file:", FILE1, "\n")
cat("Species:", SPEC, "Model:", MODEL, "Date:", DATE, "\n")

h5file <- H5Fopen(FILE1, "H5F_ACC_RDONLY")
metacluster_group <- h5read(h5file, "metacluster_idx_to_submetacluster_results")

# Optional sanity check: print metacluster names
cat("Found metaclusters:", paste(names(metacluster_group), collapse = ", "), "\n")

#######################################################
# Extract patterns from all metaclusters
#######################################################

patterns_list <- list()
for (metacluster_name in names(metacluster_group)) {
  metacluster <- metacluster_group[[metacluster_name]]
  if (!is.null(metacluster[["seqlets_to_patterns_result"]][["patterns"]])) {
    patterns_list[[metacluster_name]] <- metacluster[["seqlets_to_patterns_result"]][["patterns"]]
  } else {
    warning(paste("No patterns found in", metacluster_name))
  }
}

#######################################################
# Extract seqlet counts per pattern (needed for all matrix types)
#######################################################

get_seqlet_counts <- function(metacluster_group) {
  seqlet_counts <- list()
  
  for (metacluster_name in names(metacluster_group)) {
    seqlet_counts[[metacluster_name]] <- list()
    
    pattern_names <- metacluster_group[[metacluster_name]][["seqlets_to_patterns_result"]][["patterns"]][["all_pattern_names"]]
    
    for (pattern_name in pattern_names) {
      # Convert raw pattern names to character if needed
      pattern_key <- if (is.raw(pattern_name)) rawToChar(pattern_name) else as.character(pattern_name)
      
      seqlets <- metacluster_group[[metacluster_name]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["seqlets_and_alnmts"]][["seqlets"]]
      seqlet_counts[[metacluster_name]][[pattern_key]] <- length(seqlets)
    }
  }
  return(seqlet_counts)
}

seqlet_counts <- get_seqlet_counts(metacluster_group)

#######################################################
# Extract matrices for each format (PFM from sequence, CWM from contrib_scores)
#######################################################

extract_matrices_by_type <- function(patterns_list, matrix_type) {
  matrix_list <- list()
  
  for (metacluster_name in names(patterns_list)) {
    patterns <- patterns_list[[metacluster_name]]
    matrix_list[[metacluster_name]] <- list(fwd = list(), rev = list())
    
    pattern_names_raw <- patterns[["all_pattern_names"]]
    pattern_names <- sapply(pattern_names_raw, function(x) {
      if (is.raw(x)) rawToChar(x) else as.character(x)
    })
    
    for (pattern_name in pattern_names) {
      if (!is.null(patterns[[pattern_name]][[matrix_type]])) {
        matrix_fwd <- patterns[[pattern_name]][[matrix_type]][["fwd"]]
        matrix_rev <- patterns[[pattern_name]][[matrix_type]][["rev"]]
        
        matrix_list[[metacluster_name]][["fwd"]][[pattern_name]] <- matrix_fwd
        matrix_list[[metacluster_name]][["rev"]][[pattern_name]] <- matrix_rev
      } else {
        warning(paste("No", matrix_type, "found for", pattern_name, "in", metacluster_name))
      }
    }
  }
  return(matrix_list)
}

# Extract raw matrices for each type
sequence_matrices <- extract_matrices_by_type(patterns_list, "sequence")
cwm_matrices <- extract_matrices_by_type(patterns_list, "task0_contrib_scores")

#######################################################
# Process matrices according to their type
#######################################################

process_pfm_matrices <- function(seq_matrices, seqlet_counts) {
  # Convert sequence frequency counts to proper Position Frequency Matrices (PFMs)
  # This logic is unchanged from the original script
  processed <- list()
  
  for (metacluster_name in names(seq_matrices)) {
    processed[[metacluster_name]] <- list(fwd = list(), rev = list())
    
    for (strand in c("fwd", "rev")) {
      matrices <- seq_matrices[[metacluster_name]][[strand]]
      
      for (pattern_name in names(matrices)) {
        matrix <- matrices[[pattern_name]]
        count <- seqlet_counts[[metacluster_name]][[pattern_name]]
        
        # MULTIPLY BY SEQLET COUNT — converts normalised frequencies back to actual counts
        # This gives us the raw frequency matrix that motif analysis tools expect
        pfm <- matrix * count
        
        processed[[metacluster_name]][[strand]][[pattern_name]] <- pfm
      }
    }
  }
  return(processed)
}

process_cwm_matrices <- function(cwm_matrices, seqlet_counts) {
  # Process Contribution Weight Matrices — keep raw values initially
  # The original script had:
  # inflation where max_score = seqlet_count and abs value cwm 
  # This version preserves the raw contribution magnitudes and cwm signs for interpretation.
  processed <- list()
  
  for (metacluster_name in names(cwm_matrices)) {
    processed[[metacluster_name]] <- list(fwd = list(), rev = list())
    
    for (strand in c("fwd", "rev")) {
      matrices <- cwm_matrices[[metacluster_name]][[strand]]
      
      for (pattern_name in names(matrices)) {
        matrix <- matrices[[pattern_name]]
        # KEEP RAW CONTRIBUTION SCORES — preserves biological signal and directionality
        # Positive values = features that increase expression 
        # Negative values = features that decrease expression
        processed[[metacluster_name]][[strand]][[pattern_name]] <- matrix
      }
    }
  }
  return(processed)
}

# Process matrices for each type
cat("Processing PFMs...\n")
pfm_processed <- process_pfm_matrices(sequence_matrices, seqlet_counts)

cat("Processing CWMs...\n")
cwm_processed <- process_cwm_matrices(cwm_matrices, seqlet_counts)

#######################################################
# Apply nomenclature system
#######################################################

rename_motifs_with_seqlets <- function(motif_matrices, metacluster_name, strand, seqlet_counts) {
  # Apply systematic nomenclature to motifs
  # Format: epm_<species>_<model>_<pattern>_<metacluster>_<strand>_<seqletcount>
  renamed <- list()
  
  for (pattern_name in names(motif_matrices)) {
    count <- seqlet_counts[[metacluster_name]][[pattern_name]]
    
    # Extract pattern number (assuming format "pattern_X")
    pattern_num <- sub("pattern_", "", pattern_name)
    
    new_name <- paste0("epm_", SPEC, "_", MODEL, "_",
                       pattern_num, "_", 
                       sub("metacluster_", "", metacluster_name), "_", 
                       strand, "_", count)
    
    renamed[[new_name]] <- motif_matrices[[pattern_name]]
  }
  return(renamed)
}

apply_nomenclature <- function(processed_matrices, seqlet_counts) {
  # Apply nomenclature to all matrices in the processed structure
  result <- list()
  
  for (metacluster_name in names(processed_matrices)) {
    # Forward strand
    result[["F0"]] <- c(result[["F0"]], 
                        rename_motifs_with_seqlets(
                          processed_matrices[[metacluster_name]][["fwd"]], 
                          metacluster_name, "F", seqlet_counts))
    
    # Reverse strand  
    result[["R0"]] <- c(result[["R0"]], 
                        rename_motifs_with_seqlets(
                          processed_matrices[[metacluster_name]][["rev"]], 
                          metacluster_name, "R", seqlet_counts))
  }
  
  # Separate by metacluster for compatibility with original output format
  final_result <- list(
    matricesF0 = list(), matricesF1 = list(),
    matricesR0 = list(), matricesR1 = list()
  )
  
  # Sort motifs into appropriate lists based on naming
  all_motifs <- c(result[["F0"]], result[["R0"]])
  
  for (motif_name in names(all_motifs)) {
    if (grepl("_0_F_", motif_name)) {
      final_result[["matricesF0"]][[motif_name]] <- all_motifs[[motif_name]]
    } else if (grepl("_1_F_", motif_name)) {
      final_result[["matricesF1"]][[motif_name]] <- all_motifs[[motif_name]]
    } else if (grepl("_0_R_", motif_name)) {
      final_result[["matricesR0"]][[motif_name]] <- all_motifs[[motif_name]]
    } else if (grepl("_1_R_", motif_name)) {
      final_result[["matricesR1"]][[motif_name]] <- all_motifs[[motif_name]]
    }
  }
  
  return(final_result)
}

# Apply nomenclature to all matrix types
pfm_named <- apply_nomenclature(pfm_processed, seqlet_counts)
cwm_named_unscaled <- apply_nomenclature(cwm_processed, seqlet_counts)
cwm_named_scaled <- apply_nomenclature(cwm_processed, seqlet_counts)

# *** NOTE ON CWM SCALING ***
# up to now, this script preserves raw contribution scores which
# maintains regulatory directionality (positive = activating, negative = repressing)

####
# UPDATE ON CWM SCALING:
#
# It is neccessary to apply scaling for the motif clustering step
# so that motifs supported by more seqlets are given higher weight 
# focus on contribution magnitude rather than directionality
# allows CWMs to behave more like traditional motif matrices (expected by clustering algorithms like universalmotif)
#
# the CWM values are tailored to prioritise patterns with stronger "evidence" (more seqlets)
# and reflect what the deep learning model actually learned

#######################################################
# Apply CWM scaling logic (post-nomenclature)
#######################################################

apply_original_cwm_scaling <- function(named_matrices) {
  # Apply the CWM scaling logic after nomenclature
  # Modified to ensure compatibility with read_jaspar() function
  # abs(matrix) * seqlet_count / max(abs(matrix)) + minimum floor value
  
  all_matrices <- c(named_matrices[["matricesF0"]], named_matrices[["matricesF1"]], 
                    named_matrices[["matricesR0"]], named_matrices[["matricesR1"]])
  
  scaled_matrices <- list()
  
  for (i in seq_along(all_matrices)) {
    matrix <- all_matrices[[i]]
    motif_name <- names(all_matrices)[i]
    
    # Extract seqlet count from motif name
    seq_count <- as.numeric(sub(".*_([0-9]+)$", "\\1", motif_name))
    
    abs_matrix <- abs(matrix)
    max_val <- max(abs_matrix)
    
    if (max_val > 0) {
      scaled_matrix <- abs_matrix * seq_count / max_val
    } else {
      scaled_matrix <- abs_matrix
    }
    
    # ADD IMPROVEMENTS FOR read_jaspar() COMPATIBILITY:
    
    # 1. Add minimum floor value to avoid very small values
    # ensure no values fall in the problematic 0-1 range that's too small
    min_floor <- 2.0  # Minimum value for any matrix entry
    scaled_matrix <- pmax(scaled_matrix, min_floor)
    
    # 2. Add small random noise to break ties and avoid identical values
    #  helps with matrix singularity issues in clustering
    noise_factor <- 0.01
    set.seed(42)  # Reproducible noise - and the answer to life, the universe and everything
    noise <- matrix(runif(length(scaled_matrix), -noise_factor, noise_factor), 
                   nrow = nrow(scaled_matrix))
    scaled_matrix <- scaled_matrix + noise
    
    # 3. Ensure values are reasonably large integers (read_jaspar expects count-like values)
    scaled_matrix <- round(scaled_matrix * 10) / 10  # Round to 1 decimal place
    
    # 4. Final check: ensure all values are >= 1.0
    scaled_matrix <- pmax(scaled_matrix, 1.0)
    
    scaled_matrices[[motif_name]] <- scaled_matrix
  }
  
  # Reconstruct the original structure
  final_result <- list(
    matricesF0 = list(), matricesF1 = list(),
    matricesR0 = list(), matricesR1 = list()
  )
  
  for (motif_name in names(scaled_matrices)) {
    if (grepl("_0_F_", motif_name)) {
      final_result[["matricesF0"]][[motif_name]] <- scaled_matrices[[motif_name]]
    } else if (grepl("_1_F_", motif_name)) {
      final_result[["matricesF1"]][[motif_name]] <- scaled_matrices[[motif_name]]
    } else if (grepl("_0_R_", motif_name)) {
      final_result[["matricesR0"]][[motif_name]] <- scaled_matrices[[motif_name]]
    } else if (grepl("_1_R_", motif_name)) {
      final_result[["matricesR1"]][[motif_name]] <- scaled_matrices[[motif_name]]
    }
  }
  
  return(final_result)
}

cat("Applying original CWM scaling...\n")
cwm_named_scaled <- apply_original_cwm_scaling(cwm_named_scaled)

#######################################################
# Generate summary statistics
#######################################################

summarise_contrib_scores <- function(matrices, matrix_type = "CWM") {
  # Generate summary statistics for motif matrices - 
  # default is CWM since this represents the expression contribution motifs
  # Handles different matrix types appropriately
  all_matrices <- c(matrices[["matricesF0"]], matrices[["matricesF1"]], 
                    matrices[["matricesR0"]], matrices[["matricesR1"]])
  
  imap_dfr(all_matrices, function(mat, motif_name) {
    # Extract seqlet count from motif name (format: epm_<species>_<model>_<pattern>_<metacluster>_<strand>_<seqletcount>)
    seqlet_count <- as.numeric(sub(".*_([0-9]+)$", "\\1", motif_name))
    
    if (matrix_type == "CWM") {
      # For CWMs, we care about total contribution magnitude
      tibble(
        motif = motif_name,
        matrix_type = matrix_type,
        seqlet_count = seqlet_count,
        score_sum = sum(abs(mat)),  # Total absolute contribution
        score_max = max(mat),
        score_min = min(mat),
        score_range = max(mat) - min(mat)
      )
    } else {
      # For PFMs, we summarise frequency totals
      tibble(
        motif = motif_name,
        matrix_type = matrix_type,
        seqlet_count = seqlet_count,
        score_sum = sum(mat),  # Total frequency counts
        score_max = max(mat),
        score_min = min(mat),
        score_range = max(mat) - min(mat)
      )
    }
  })
}

# Generate summaries for all matrix types
cwm_summary <- summarise_contrib_scores(cwm_named_scaled, "CWM")
cwm_unscaled_summary <- summarise_contrib_scores(cwm_named_unscaled, "CWM_unscaled")
pfm_summary <- summarise_contrib_scores(pfm_named, "PFM")

# Combine all summaries
combined_summary <- bind_rows(cwm_summary, cwm_unscaled_summary, pfm_summary)

#######################################################
# Write output files
#######################################################

write_jaspar_motifs <- function(motif_matrices, out_file, matrix_type = "motif") {
  # Write motifs to JASPAR format file
  # Handles different matrix types and ensures proper formatting

  rows <- c('A', 'C', 'G', 'T')
  motif_strs <- c()
  
  # Combine all matrices from the four lists
  all_matrices <- c(motif_matrices[["matricesF0"]], motif_matrices[["matricesF1"]], 
                    motif_matrices[["matricesR0"]], motif_matrices[["matricesR1"]])
  
  for (motif_name in names(all_matrices)) {
    mat <- all_matrices[[motif_name]]
    
    # Start motif entry with header
    motif_str <- paste0(">", motif_name, "\n")
    
    # Format matrix values based on type
    for (j in 1:nrow(mat)) {
      # For PFM and CWM, use appropriate precision
      values <- as.character(round(mat[j, ], 4))
      motif_str <- paste0(motif_str, rows[j], " [", paste(values, collapse = "\t"), "]\n")
    }
    motif_strs <- c(motif_strs, motif_str)
  }
  
  writeLines(motif_strs, out_file)
}

# Write output files
cat("Writing output files...\n")

# Summary statistics - separate files for each matrix type
cwm_summary_file <- file.path(dirpath_out, paste0(DATE, "_", SPEC, MODEL, "_cwm_summary.csv"))
cwm_unscaled_summary_file <- file.path(dirpath_out, paste0(DATE, "_", SPEC, MODEL, "_cwm_unscaled_summary.csv"))
pfm_summary_file <- file.path(dirpath_out, paste0(DATE, "_", SPEC, MODEL, "_pfm_summary.csv"))

write.csv(cwm_summary, file = cwm_summary_file, row.names = FALSE)
write.csv(cwm_unscaled_summary, file = cwm_unscaled_summary_file, row.names = FALSE)
write.csv(pfm_summary, file = pfm_summary_file, row.names = FALSE)

cat("CWM summary statistics written to:", cwm_summary_file, "\n")
cat("CWM unscaled summary statistics written to:", cwm_unscaled_summary_file, "\n")
cat("PFM summary statistics written to:", pfm_summary_file, "\n")

# JASPAR format files with original naming convention
pfm_jaspar <- file.path(dirpath_out, paste0("rdf5_PFM_pattern", SPEC, MODEL, "_pfm-motifs.jaspar"))
cwm_jaspar_scaled <- file.path(dirpath_out, paste0("rdf5_epm", SPEC, MODEL, "_cwm_scaled-motifs.jaspar"))
cwm_jaspar_unscaled <- file.path(dirpath_out, paste0("rdf5_epm", SPEC, MODEL, "_cwm_unscaled-motifs.jaspar"))

write_jaspar_motifs(pfm_named, pfm_jaspar, "PFM")
write_jaspar_motifs(cwm_named_scaled, cwm_jaspar_scaled, "CWM")
write_jaspar_motifs(cwm_named_unscaled, cwm_jaspar_unscaled, "CWM")

cat("PFM motifs written to:", pfm_jaspar, "\n")
cat("CWM scaled motifs written to:", cwm_jaspar_scaled, "\n")
cat("CWM unscaled motifs written to:", cwm_jaspar_unscaled, "\n")

# Close HDF5 file
H5Fclose(h5file)

cat("Processing complete!\n")
cat("Found", nrow(cwm_summary), "CWM and", nrow(pfm_summary), "PFM motifs\n")
cat("Use CWM scaled file for clustering and CWM unscaled file for interpretation\n")