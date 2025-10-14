#!/usr/bin/env Rscript
# Split shared gene lists by expression class (high/low)
# Creates separate gene list files for GO enrichment analysis
#
# Usage:
#   Single file:
#     Rscript split_by_expression.R /path/to/shared_genes.csv
#
#   Batch processing (wildcard pattern):
#     Rscript split_by_expression.R "../../../out/moca_results/mo_go/clade_shared_genes/p1_clades/*_shared_genes.csv"
#
# Example:
#   Rscript split_by_expression.R ../../../out/moca_results/mo_go/clade_shared_genes/p1_clades/p1m_l1_m6_pink_shared_genes.csv
#
# Output:
#   - {original_name}_low.csv - Contains only low expression genes
#   - {original_name}_high.csv - Contains only high expression genes
#
# The output files are compatible with 04_go_enrichment.R
######################

library(dplyr)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for required input
if (length(args) < 1 || !nzchar(args[1])) {
  stop("\nERROR: Gene list file path is required!\n\n",
       "Usage: Rscript split_by_expression.R /path/to/shared_genes.csv\n",
       "   or: Rscript split_by_expression.R \"/path/to/*_shared_genes.csv\"\n\n",
       "Example:\n",
       "  Rscript split_by_expression.R ../../../out/moca_results/mo_go/clade_shared_genes/p1_clades/p1m_l1_m6_pink_shared_genes.csv\n")
}

INPUT_PATH <- args[1]

cat("=", rep("=", 70), "\n", sep = "")
cat("SPLIT GENE LISTS BY EXPRESSION CLASS\n")
cat(rep("=", 70), "\n\n", sep = "")

#' Expand glob pattern to list of files
#'
#' @param path String that may contain wildcards (*, ?)
#' @return Character vector of matching file paths
expand_file_list <- function(path) {
  # Check if path contains wildcard characters
  if (grepl("[*?]", path)) {
    # Expand glob pattern
    files <- Sys.glob(path)

    if (length(files) == 0) {
      stop("No files matched the pattern: ", path)
    }

    # Filter to only keep .csv files that contain "_shared_genes.csv"
    files <- files[grepl("_shared_genes\\.csv$", files)]

    if (length(files) == 0) {
      stop("No shared genes CSV files matched the pattern: ", path)
    }

    return(files)
  } else {
    # Single file - check it exists
    if (!file.exists(path)) {
      stop("File not found: ", path)
    }
    return(path)
  }
}

#' Read gene list with expression classes (handles comment lines)
#'
#' @param file_path Path to shared genes CSV
#' @return Data frame with gene_id and expression_class columns
read_shared_genes <- function(file_path) {
  # Read all lines
  all_lines <- readLines(file_path)

  # Find first non-comment line (this should be the header)
  non_comment_lines <- which(!grepl("^#", all_lines) & !grepl("^\"#", all_lines))

  if (length(non_comment_lines) == 0) {
    stop("No data found in file: ", file_path)
  }

  # First non-comment line is the header
  header_line <- non_comment_lines[1]

  # Number of lines to skip before header
  skip_lines <- header_line - 1

  # Read CSV starting from header
  gene_data <- read.csv(file_path, skip = skip_lines, stringsAsFactors = FALSE)

  # Check required columns
  if (!all(c("gene_id", "expression_class") %in% colnames(gene_data))) {
    stop("File must have 'gene_id' and 'expression_class' columns: ", file_path)
  }

  # Remove any remaining comment rows and duplicate headers
  gene_data <- gene_data %>%
    filter(!grepl("^#", gene_id)) %>%
    filter(gene_id != "gene_id")

  return(gene_data)
}

#' Process a single shared genes file
#'
#' Splits the file into separate high and low expression gene lists
#'
#' @param file_path Path to shared genes CSV
#' @return List with success status and file counts
process_single_file <- function(file_path) {

  cat("\n", rep("-", 70), "\n", sep = "")
  cat("Processing:", basename(file_path), "\n")
  cat(rep("-", 70), "\n", sep = "")

  tryCatch({

    # Read gene list with expression classes
    cat("Reading gene list...\n")
    gene_data <- read_shared_genes(file_path)

    # Count expression classes
    n_total <- nrow(gene_data)
    n_low <- sum(gene_data$expression_class == "low")
    n_high <- sum(gene_data$expression_class == "high")

    cat("  Total genes:", n_total, "\n")
    cat("  Low expression:", n_low, "\n")
    cat("  High expression:", n_high, "\n")

    # Check if we have both classes
    if (n_low == 0) {
      cat("  WARNING: No low expression genes. Skipping split.\n")
      return(list(success = FALSE, reason = "no_low"))
    }

    if (n_high == 0) {
      cat("  WARNING: No high expression genes. Skipping split.\n")
      return(list(success = FALSE, reason = "no_high"))
    }

    # Create output file paths
    base_path <- tools::file_path_sans_ext(file_path)
    output_low <- paste0(base_path, "_low.csv")
    output_high <- paste0(base_path, "_high.csv")

    # Split by expression class
    genes_low <- gene_data %>% filter(expression_class == "low")
    genes_high <- gene_data %>% filter(expression_class == "high")

    # Get original comment lines (metadata)
    all_lines <- readLines(file_path)
    comment_lines <- all_lines[grepl("^#", all_lines)]

    # Write low expression genes
    cat("Writing low expression genes...\n")

    # Create comment header for low expression file
    if (length(comment_lines) > 0) {
      # Update gene count in comments if present
      updated_comments_low <- sapply(comment_lines, function(line) {
        # Update shared gene count if line contains it
        if (grepl("Shared genes:", line)) {
          gsub("Shared genes: [0-9]+", paste0("Shared genes: ", n_low, " (low expression)"), line)
        } else {
          line
        }
      })
      writeLines(updated_comments_low, output_low)

      # Append data
      write.table(genes_low %>% select(gene_id),
                  output_low,
                  append = TRUE,
                  sep = ",",
                  row.names = FALSE,
                  col.names = TRUE,
                  quote = FALSE)
    } else {
      # No comments, just write data
      write.csv(genes_low %>% select(gene_id),
                output_low,
                row.names = FALSE,
                quote = FALSE)
    }
    cat("  Saved:", basename(output_low), "(", n_low, "genes )\n")

    # Write high expression genes
    cat("Writing high expression genes...\n")

    if (length(comment_lines) > 0) {
      # Update gene count in comments if present
      updated_comments_high <- sapply(comment_lines, function(line) {
        # Update shared gene count if line contains it
        if (grepl("Shared genes:", line)) {
          gsub("Shared genes: [0-9]+", paste0("Shared genes: ", n_high, " (high expression)"), line)
        } else {
          line
        }
      })
      writeLines(updated_comments_high, output_high)

      # Append data
      write.table(genes_high %>% select(gene_id),
                  output_high,
                  append = TRUE,
                  sep = ",",
                  row.names = FALSE,
                  col.names = TRUE,
                  quote = FALSE)
    } else {
      # No comments, just write data
      write.csv(genes_high %>% select(gene_id),
                output_high,
                row.names = FALSE,
                quote = FALSE)
    }
    cat("  Saved:", basename(output_high), "(", n_high, "genes )\n")

    return(list(
      success = TRUE,
      n_low = n_low,
      n_high = n_high,
      output_low = output_low,
      output_high = output_high
    ))

  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(list(success = FALSE, reason = "error", error = e$message))
  })
}

# Main execution
cat("Starting gene list splitting...\n\n")

# Expand input path (handles both single files and wildcards)
FILE_LIST <- expand_file_list(INPUT_PATH)

# Determine if batch mode
BATCH_MODE <- length(FILE_LIST) > 1

cat("Input mode:", if (BATCH_MODE) "BATCH" else "SINGLE FILE", "\n")
cat("Files to process:", length(FILE_LIST), "\n")
if (BATCH_MODE) {
  cat("  (Showing first 5):\n")
  for (f in head(FILE_LIST, 5)) {
    cat("    -", basename(f), "\n")
  }
  if (length(FILE_LIST) > 5) {
    cat("    ... and", length(FILE_LIST) - 5, "more\n")
  }
} else {
  cat("  File:", FILE_LIST, "\n")
}
cat("\n")

# Process all files
results_summary <- data.frame(
  file = character(),
  status = character(),
  n_low = integer(),
  n_high = integer(),
  stringsAsFactors = FALSE
)

for (file_path in FILE_LIST) {
  result <- process_single_file(file_path)

  # Track results
  results_summary <- rbind(results_summary, data.frame(
    file = basename(file_path),
    status = if (result$success) "SUCCESS" else result$reason,
    n_low = if (result$success) result$n_low else 0,
    n_high = if (result$success) result$n_high else 0,
    stringsAsFactors = FALSE
  ))
}

# Final summary
cat("\n", rep("=", 70), "\n", sep = "")
cat("PROCESSING COMPLETE\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Summary:\n")
cat("  Total files processed:", nrow(results_summary), "\n")
cat("  Successful:", sum(results_summary$status == "SUCCESS"), "\n")
cat("  Failed (no low):", sum(results_summary$status == "no_low"), "\n")
cat("  Failed (no high):", sum(results_summary$status == "no_high"), "\n")
cat("  Failed (errors):", sum(results_summary$status == "error"), "\n\n")

if (BATCH_MODE) {
  # Print detailed summary for batch mode
  cat("Detailed Results:\n")
  print(results_summary, row.names = FALSE)
}

cat("\nGene Totals:\n")
cat("  Total low expression genes:", sum(results_summary$n_low), "\n")
cat("  Total high expression genes:", sum(results_summary$n_high), "\n")

cat("\nNext Steps:\n")
cat("  Run GO enrichment on the split files using:\n")
cat("    Rscript 04_go_enrichment.R \"/path/to/*_low.csv\"\n")
cat("    Rscript 04_go_enrichment.R \"/path/to/*_high.csv\"\n")

cat("\nCompleted at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
