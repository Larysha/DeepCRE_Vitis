# motif_ranges_meta_fixed.R
# Calculate motif positional characteristics for TSS and TTS regions
# Analyses seqlet position data to generate summary statistics for motifs
# DeepCRE/moca/mo_range
######################

library(dplyr)

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
default_date <- format(Sys.Date(), "%Y%m%d")
default_input_dir <- "../out/moca_results/mo_range"
default_output_dir <- "../out/moca_results/mo_range"
# Assign arguments or defaults
SPEC <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_spec
MODEL <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_model
DATE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_date

# Find the source HDF5 file dynamically
if (length(args) >= 4 && nzchar(args[4])) {
  SOURCE_FILE <- args[4]
} else {
  cat("Searching for HDF5 source file in ../out/modisco/...\n")
  source_file_path <- find_project_files(
    directory = "../out/modisco",
    pattern = ".*motifs.*\\.hdf5$",
    description = "MoDISco HDF5 source file",
    select_multiple = FALSE
  )
  SOURCE_FILE <- basename(source_file_path)  # Just keep the filename for metadata
}

cat("Arguments: SPEC =", SPEC, ", MODEL =", MODEL, ", DATE =", DATE, "\n")

# Input file selection
if (length(args) >= 5 && nzchar(args[5])) {
  seqlet_file <- args[5]
} else {
  # Construct expected filename
  seqlet_filename <- paste0("rdf5_seqlet_pattern_", SPEC, "_", MODEL, "_", DATE, ".csv")
  seqlet_file <- file.path(default_input_dir, seqlet_filename)
  
  if (!file.exists(seqlet_file)) {
    cat("Specific file not found, searching for seqlet pattern files...\n")
    seqlet_file <- find_project_files(
      directory = default_input_dir,
      pattern = paste0("rdf5_seqlet_pattern_", SPEC, "_", MODEL, ".*\\.csv$"),
      description = "seqlet pattern file",
      select_multiple = FALSE
    )
  }
}

# Output directory
output_dir <- if (length(args) >= 6 && nzchar(args[6])) args[6] else default_output_dir

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

cat("Processing seqlet file:", seqlet_file, "\n")
cat("Output directory:", output_dir, "\n")

#' Calculate motif positional characteristics for TSS and TTS regions
#' 
#' Analyses seqlet position data to generate summary statistics for motifs
#' in TSS-proximal (1-1500bp) and TTS-proximal (1520-3000bp) regions.
#' TTS coordinates are flipped to create a mirrored 5'->3' view for comparison.
calculate_motif_ranges <- function(seqlet_file, species_code, model_code, source_file, output_dir) {
  
  # Read and validate data
  cat("Reading seqlet data...\n")
  
  if (!file.exists(seqlet_file)) {
    stop("Input file not found: ", seqlet_file)
  }
  
  data <- tryCatch({
    read.csv(seqlet_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop("Failed to read input file: ", e$message)
  })
  
  cat("  ✓ Read", nrow(data), "seqlet records\n")
  
  # Validate required columns
  required_cols <- c("start", "end", "metacluster", "pattern")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  cat("  ✓ Data validation passed\n")
  
  # Create motif identifiers (following original script logic)
  data$motif <- ifelse(data$metacluster == "metacluster_0", "p0", "p1")
  data$motif <- paste0(data$motif, "m", sprintf("%02d", as.numeric(substring(data$pattern, 9))))
  
  cat("  ✓ Created motif identifiers for", length(unique(data$motif)), "unique motifs\n")
  
  # Split into TSS and TTS ranges
  tss_range <- data %>% filter(start >= 1 & end <= 1500) 
  tts_range <- data %>% filter(start >= 1520 & end <= 3000)
  
  cat("  ✓ TSS range:", nrow(tss_range), "seqlets\n")
  cat("  ✓ TTS range:", nrow(tts_range), "seqlets\n")
  
  # Apply coordinate transformation for TTS (flip to 5'->3' orientation)
  if (nrow(tts_range) > 0) {
    tts_range$start <- 3020 - tts_range$start
    tts_range$end <- 3020 - tts_range$end
  }
  
  # Create binned positions for mode calculation (10bp bins)
  if (nrow(tss_range) > 0) {
    tss_range$trunc_start <- floor(tss_range$start / 10)
  }
  if (nrow(tts_range) > 0) {
    tts_range$trunc_start <- floor(tts_range$start / 10)
  }
  
  # Custom mode function
  customMode <- function(x) {
    freq <- table(x)
    mode_val <- as.numeric(names(freq)[which.max(freq)])
    return(mode_val)
  }
  
  # Calculate summary statistics for both regions
  calculate_stats <- function(range_data, region_name) {
    if (nrow(range_data) == 0) {
      cat("⚠ Warning: No data found for", region_name, "region\n")
      return(data.frame())
    }
    
    cat("Processing", region_name, "statistics...\n")
    
    stats <- range_data %>%
      group_by(motif) %>%
      summarise(
        min = min(start),
        max = max(start),
        q10 = quantile(start, 0.1),
        median = median(start),
        q90 = quantile(start, 0.9),
        mode = customMode(trunc_start) * 10, # Convert back from binned values
        mean = mean(start),
        sd = sd(start),
        cv = sd(start) / mean(start) * 100,
        iqr = q90 - q10,
        number = n(),
        .groups = 'drop'
      )
    
    # Add metadata
    stats$Species <- species_code
    stats$Model <- model_code
    stats$source <- source_file
    stats$epm <- paste("epm", stats$Species, stats$Model, stats$motif, sep = "_")
    
    # Reorder columns (following original script column order)
    stats <- stats %>%
      select(epm, min, max, mean, median, mode, q10, q90, sd, cv, iqr, number, source)
    
    cat("  ✓", nrow(stats), "motifs processed for", region_name, "\n")
    return(stats)
  }
  
  # Process both regions
  tss_results <- calculate_stats(tss_range, "TSS")
  tts_results <- calculate_stats(tts_range, "TTS")
  
  # Write output files
  output_prefix <- file.path(output_dir, paste0(species_code, model_code))
  
  tss_output <- paste0(output_prefix, "-TSS_motif_ranges_q1q9.csv")
  tts_output <- paste0(output_prefix, "-TTS_motif_ranges_q1q9.csv")
  
  if (nrow(tss_results) > 0) {
    write.csv(tss_results, file = tss_output, row.names = FALSE)
    cat("✓ TSS results written to:", tss_output, "\n")
  }
  
  if (nrow(tts_results) > 0) {
    write.csv(tts_results, file = tts_output, row.names = FALSE)
    cat("✓ TTS results written to:", tts_output, "\n")
  }
  
  # Summary report
  cat("\n=== MOTIF RANGE ANALYSIS COMPLETE ===\n")
  cat("Input sequences:", nrow(data), "\n")
  cat("TSS motifs found:", nrow(tss_results), "\n")
  cat("TTS motifs found:", nrow(tts_results), "\n")
  cat("=====================================\n")
  
  # Return results invisibly
  invisible(list(tss = tss_results, tts = tts_results))
}

# Main execution
main <- function() {
  tryCatch({
    result <- calculate_motif_ranges(
      seqlet_file = seqlet_file,
      species_code = SPEC,
      model_code = MODEL,
      source_file = SOURCE_FILE,
      output_dir = output_dir
    )
    
    return(result)
  }, error = function(e) {
    cat("Analysis failed:", e$message, "\n")
    stop(e)
  })
}

# Execute the main function
cat("Starting motif range analysis pipeline...\n")
result <- main()
cat("Pipeline completed successfully!\n")