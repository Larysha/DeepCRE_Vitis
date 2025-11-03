# get_seqlets_per_pattern_fixed.R
# Extracts seqlet positional information from TF-MoDISco HDF5 output
# Fixed version based on working original script structure
######################

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(rhdf5)

# Set working directory to script location
script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
if (nchar(script_dir) > 0) {
  setwd(script_dir)
  cat("Working directory set to:", getwd(), "\n")
}

# Source utility functions
source("../utils.R")

args <- commandArgs(trailingOnly = TRUE)

default_spec <- "vitis"
default_model <- "ssr"
default_input_dir <- "../../../out/modisco"
default_output_dir <- "../../../out/moca_results/mo_range"

# Assign arguments or defaults
SPEC <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_spec
MODEL <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_model
DATE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else format(Sys.Date(), "%Y%m%d")

cat("Arguments: SPEC =", SPEC, ", MODEL =", MODEL, ", DATE =", DATE, "\n")

# Input file selection
if (length(args) >= 4 && nzchar(args[4])) {
  FILE1 <- args[4]
} else {
  cat("No input file provided, searching for HDF5 file...\n")
  FILE1 <- find_project_files(
    directory = default_input_dir,
    pattern = ".*motifs.*\\.hdf5$",
    description = "MoDISco HDF5 file",
    select_multiple = FALSE
  )
}

# Output directory
dirpath_out <- if (length(args) >= 5 && nzchar(args[5])) args[5] else default_output_dir

# Create output directory
if (!dir.exists(dirpath_out)) {
  dir.create(dirpath_out, recursive = TRUE)
  cat("Created output directory:", dirpath_out, "\n")
}

# Define output file name
output_filename <- paste0("rdf5_seqlet_pattern_", SPEC, "_", MODEL, "_", DATE, ".csv")
output_file <- file.path(dirpath_out, output_filename)

cat("Processing MoDISco file:", FILE1, "\n")
cat("Output will be saved to:", output_file, "\n")

# Function to safely extract seqlets from a pattern
extract_pattern_seqlets <- function(h5file, metacluster_name, pattern_name) {
  
  seqlet_path <- paste0("metacluster_idx_to_submetacluster_results/", 
                        metacluster_name, 
                        "/seqlets_to_patterns_result/patterns/", 
                        pattern_name, 
                        "/seqlets_and_alnmts/seqlets")
  
  seqlets <- tryCatch({
    h5read(h5file, seqlet_path)
  }, error = function(e) {
    cat("⚠ Failed to read seqlets for", pattern_name, "in", metacluster_name, "\n")
    return(NULL)
  })
  
  if (is.null(seqlets) || length(seqlets) == 0) {
    return(data.frame())
  }
  
  # Convert byte strings to character if needed
  if (is.raw(seqlets) || any(class(seqlets) %in% c("array", "matrix"))) {
    seqlets <- as.character(seqlets)
  }
  
  # Create data frame following original script structure
  seqlets_df <- data.frame(
    seqlets = seqlets,
    pattern = pattern_name,
    metacluster = metacluster_name,
    stringsAsFactors = FALSE
  )
  
  cat("  ✓ Extracted", nrow(seqlets_df), "seqlets for", pattern_name, "\n")
  return(seqlets_df)
}

# Function to process a single metacluster
process_metacluster <- function(h5file, metacluster_name) {
  cat("Processing", metacluster_name, "...\n")
  
  # Get patterns - following original script logic
  patterns_path <- paste0("metacluster_idx_to_submetacluster_results/", 
                         metacluster_name, 
                         "/seqlets_to_patterns_result/patterns")
  
  patterns <- tryCatch({
    h5read(h5file, patterns_path)
  }, error = function(e) {
    cat("⚠ Failed to access patterns for", metacluster_name, "\n")
    return(NULL)
  })
  
  if (is.null(patterns)) {
    return(data.frame())
  }
  
  # Get pattern names, excluding metadata (following original script)
  pattern_names <- names(patterns)
  pattern_names <- pattern_names[grepl("^pattern_\\d+$", pattern_names)]
  
  if (length(pattern_names) == 0) {
    cat("⚠ No valid patterns found in", metacluster_name, "\n")
    return(data.frame())
  }
  
  cat("  Found", length(pattern_names), "patterns\n")
  
  # Extract seqlets for each pattern
  all_seqlets <- data.frame()
  
  for (pattern_name in pattern_names) {
    pattern_seqlets <- extract_pattern_seqlets(h5file, metacluster_name, pattern_name)
    
    if (nrow(pattern_seqlets) > 0) {
      all_seqlets <- rbind(all_seqlets, pattern_seqlets)
    }
  }
  
  return(all_seqlets)
}

# Main execution function
main <- function() {
  # Validate input file
  if (!file.exists(FILE1)) {
    stop("Input HDF5 file does not exist: ", FILE1)
  }
  
  cat("Opening HDF5 file...\n")
  
  # Open HDF5 file
  h5file <- tryCatch({
    H5Fopen(FILE1, "H5F_ACC_RDONLY")
  }, error = function(e) {
    stop("Failed to open HDF5 file: ", e$message)
  })
  
  # Ensure file closure on exit
  on.exit({
    if (exists("h5file") && !is.null(h5file)) {
      tryCatch(H5Fclose(h5file), error = function(e) NULL)
    }
  }, add = TRUE)
  
  # Get metacluster names
  metacluster_group <- tryCatch({
    h5read(h5file, "metacluster_idx_to_submetacluster_results")
  }, error = function(e) {
    stop("Failed to read metacluster results: ", e$message)
  })
  
  metacluster_names <- names(metacluster_group)
  cat("Found metaclusters:", paste(metacluster_names, collapse = ", "), "\n")
  
  if (length(metacluster_names) == 0) {
    stop("No metaclusters found in HDF5 file")
  }
  
  # Process all metaclusters
  all_seqlets <- data.frame()
  
  for (metacluster_name in metacluster_names) {
    metacluster_seqlets <- process_metacluster(h5file, metacluster_name)
    
    if (nrow(metacluster_seqlets) > 0) {
      all_seqlets <- rbind(all_seqlets, metacluster_seqlets)
    }
  }
  
  if (nrow(all_seqlets) == 0) {
    stop("No seqlets extracted from any metacluster")
  }
  
  cat("Total seqlets extracted:", nrow(all_seqlets), "\n")
  
  # Parse seqlets data (following original script exactly)
  cat("Parsing seqlet data...\n")
  
  parsed_seqlets <- all_seqlets %>%
    mutate(example = NA, start = NA, end = NA, rc = NA) %>%
    separate(col = seqlets, 
             into = c("example", "start", "end", "rc"), 
             sep = "[,]", 
             fill = "right") %>%
    mutate(
      example = gsub("example:", "", example),
      start = as.integer(gsub("start:", "", start)),
      end = as.integer(gsub("end:", "", end)),
      rc = gsub("rc:", "", rc)
    ) %>%
    select(metacluster, pattern, example, start, end, rc)
  
  # Check for parsing issues
  na_count <- sum(is.na(parsed_seqlets$start) | is.na(parsed_seqlets$end))
  if (na_count > 0) {
    cat("⚠ Warning:", na_count, "seqlets have missing position data\n")
  }
  
  # Write main output
  write.csv(parsed_seqlets, file = output_file, row.names = FALSE)
  cat("✓ Wrote", nrow(parsed_seqlets), "seqlets to:", output_file, "\n")
  
  # Generate and write summary
  summary_df <- parsed_seqlets %>%
    group_by(metacluster, pattern) %>%
    summarise(
      seqlet_count = n(),
      mean_start = mean(start, na.rm = TRUE),
      mean_end = mean(end, na.rm = TRUE),
      mean_length = mean(end - start, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_file <- file.path(dirpath_out, 
                           paste0("seqlet_summary_", SPEC, "_", MODEL, "_", DATE, ".csv"))
  write.csv(summary_df, file = summary_file, row.names = FALSE)
  cat("✓ Wrote summary to:", summary_file, "\n")
  
  # Print final summary
  cat("\n=== EXTRACTION COMPLETE ===\n")
  cat("Total seqlets:", nrow(parsed_seqlets), "\n")
  cat("Metaclusters:", length(unique(parsed_seqlets$metacluster)), "\n")
  cat("Patterns:", length(unique(parsed_seqlets$pattern)), "\n")
  cat("===========================\n")
  
  return(invisible(list(
    seqlets = parsed_seqlets,
    summary = summary_df,
    output_file = output_file,
    summary_file = summary_file
  )))
}

# Execute the main function
cat("Starting seqlet extraction pipeline...\n")
result <- main()
cat("Pipeline completed successfully!\n")