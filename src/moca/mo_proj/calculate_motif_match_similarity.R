#!/usr/bin/env Rscript

# Calculate nucleotide-level similarity between BLAMM motif matches and consensus sequences
# This script extracts actual matched sequences and calculates percent identity to the motif consensus
#
# Usage: Rscript calculate_motif_match_similarity.R
# Output: CSV with match similarity statistics for each motif occurrence

suppressPackageStartupMessages({
  library(Biostrings)
})

# File paths
ENRICHMENT_FILE <- "../../../out/moca_results/mo_proj/enrichment/20251103_vitis_motif_enrichment_integrated_data.csv"
FASTA_FILE <- "../../../out/moca_results/ref_seq/genomic_coords_vitis_vinifera_PN40024_1000bp_flanks.fa"
OUTPUT_FILE <- "../../../out/moca_results/mo_proj/enrichment/motif_match_similarity_analysis.csv"
SUMMARY_FILE <- "../../../out/moca_results/mo_proj/enrichment/motif_match_similarity_summary.csv"

cat("Loading enrichment data...\n")
enrichment_data <- read.csv(ENRICHMENT_FILE, stringsAsFactors = FALSE)

cat(sprintf("Loaded %d motif occurrences\n", nrow(enrichment_data)))

# Load FASTA sequences
cat("Loading FASTA sequences (this may take a moment)...\n")
fasta_sequences <- readDNAStringSet(FASTA_FILE)

# Create a lookup for sequence IDs
fasta_df <- data.frame(
  header = names(fasta_sequences),
  sequence = as.character(fasta_sequences),
  stringsAsFactors = FALSE
)
# Extract gene_id from FASTA header (first field before tab/space)
fasta_df$gene_id <- sub("\t.*", "", fasta_df$header)
fasta_df$gene_id <- sub(" .*", "", fasta_df$gene_id)

cat(sprintf("Loaded %d sequences from FASTA\n", nrow(fasta_df)))

# Function to calculate similarity between matched sequence and consensus
calculate_similarity <- function(matched_seq, consensus_seq) {
  # Handle NA or empty sequences
  if (is.na(matched_seq) || is.na(consensus_seq) ||
      nchar(matched_seq) == 0 || nchar(consensus_seq) == 0) {
    return(list(
      pct_identity = NA_real_,
      num_matches = NA_integer_,
      num_positions = NA_integer_
    ))
  }

  # Ensure same length
  if (nchar(matched_seq) != nchar(consensus_seq)) {
    return(list(
      pct_identity = NA_real_,
      num_matches = NA_integer_,
      num_positions = NA_integer_
    ))
  }

  # Convert to uppercase
  matched_seq <- toupper(matched_seq)
  consensus_seq <- toupper(consensus_seq)

  # IUPAC nucleotide codes
  iupac_match <- function(observed, consensus_code) {
    matches <- list(
      "A" = "A", "C" = "C", "G" = "G", "T" = "T",
      "R" = c("A", "G"), "Y" = c("C", "T"),
      "S" = c("G", "C"), "W" = c("A", "T"),
      "K" = c("G", "T"), "M" = c("A", "C"),
      "B" = c("C", "G", "T"), "D" = c("A", "G", "T"),
      "H" = c("A", "C", "T"), "V" = c("A", "C", "G"),
      "N" = c("A", "C", "G", "T")
    )

    allowed <- matches[[consensus_code]]
    if (is.null(allowed)) allowed <- c("A", "C", "G", "T")  # Treat unknown as N

    observed %in% allowed
  }

  # Split into characters
  matched_chars <- strsplit(matched_seq, "")[[1]]
  consensus_chars <- strsplit(consensus_seq, "")[[1]]

  # Check each position
  matches <- mapply(iupac_match, matched_chars, consensus_chars)

  num_matches <- sum(matches)
  num_positions <- length(matches)
  pct_identity <- (num_matches / num_positions) * 100

  list(
    pct_identity = pct_identity,
    num_matches = num_matches,
    num_positions = num_positions
  )
}

# Debug: Check gene IDs
cat("Sample gene_ids from enrichment data:\n")
cat(paste(head(unique(enrichment_data$gene_id), 3), collapse = ", "), "\n")
cat("Sample gene_ids from FASTA data:\n")
cat(paste(head(unique(fasta_df$gene_id), 3), collapse = ", "), "\n")

# Check for overlap
test_gene <- "Vitvi05_01chr01g00630"
cat(sprintf("Test gene '%s' in enrichment: %s\n", test_gene, test_gene %in% enrichment_data$gene_id))
cat(sprintf("Test gene '%s' in FASTA: %s\n", test_gene, test_gene %in% fasta_df$gene_id))

# Find common genes
common_genes <- intersect(enrichment_data$gene_id, fasta_df$gene_id)
cat(sprintf("Common genes between enrichment and FASTA: %d\n", length(common_genes)))

if (length(common_genes) == 0) {
  cat("ERROR: No common genes found! Checking for encoding issues...\n")
  cat(sprintf("Enrichment gene_id class: %s\n", class(enrichment_data$gene_id)))
  cat(sprintf("FASTA gene_id class: %s\n", class(fasta_df$gene_id)))
  quit(status = 1)
}

# Merge enrichment data with sequences
cat("Merging enrichment data with sequences...\n")

# Check if there's a 'sequence' column in enrichment_data already
if ("sequence" %in% colnames(enrichment_data)) {
  cat("WARNING: enrichment_data already has a 'sequence' column!\n")
  enrichment_data$sequence <- NULL  # Remove it
}

enrichment_with_seq <- merge(
  enrichment_data,
  fasta_df[, c("gene_id", "sequence")],
  by = "gene_id",
  all.x = TRUE,
  suffixes = c("", "_fasta")
)

# Debug: Check merge results and column names
cat(sprintf("After merge: %d rows\n", nrow(enrichment_with_seq)))
cat("Columns containing 'sequence':", paste(grep("sequence", colnames(enrichment_with_seq), value = TRUE), collapse = ", "), "\n")
if ("sequence" %in% colnames(enrichment_with_seq)) {
  cat(sprintf("Sequences (using 'sequence' column): %d non-NA\n", sum(!is.na(enrichment_with_seq$sequence))))
}
if ("sequence_fasta" %in% colnames(enrichment_with_seq)) {
  cat(sprintf("Sequences (using 'sequence_fasta' column): %d non-NA\n", sum(!is.na(enrichment_with_seq$sequence_fasta))))
  # Use the correct column name
  enrichment_with_seq$sequence <- enrichment_with_seq$sequence_fasta
}

cat("Extracting matched sequences from FASTA...\n")

# Extract actual matched sequences
enrichment_with_seq$actual_match <- rep(NA_character_, nrow(enrichment_with_seq))

for (i in 1:nrow(enrichment_with_seq)) {
  seq <- enrichment_with_seq$sequence[i]
  start <- enrichment_with_seq$mstart[i]
  end <- enrichment_with_seq$mend[i]

  # Check for NA/NULL values explicitly
  if (!is.null(seq) && !is.null(start) && !is.null(end) &&
      length(seq) > 0 && length(start) > 0 && length(end) > 0) {
    if (!is.na(seq) && !is.na(start) && !is.na(end)) {
      if (start <= nchar(seq) && end <= nchar(seq) && start > 0) {
        enrichment_with_seq$actual_match[i] <- substr(seq, start, end)
      }
    }
  }

  # Progress indicator for large datasets
  if (i %% 50000 == 0) {
    cat(sprintf("  Extracted %d/%d sequences\n", i, nrow(enrichment_with_seq)))
  }
}

cat("Calculating similarity metrics...\n")

# Calculate similarity for each match (with progress)
n_total <- nrow(enrichment_with_seq)
cat(sprintf("Processing %d matches", n_total))

if (n_total < 1000) {
  # Small dataset - process all at once
  similarity_list <- mapply(
    calculate_similarity,
    enrichment_with_seq$actual_match,
    enrichment_with_seq$motif_sequence,
    SIMPLIFY = FALSE
  )
  cat(" - Done!\n")
} else {
  # Large dataset - process in chunks with progress
  chunk_size <- 10000
  n_chunks <- ceiling(n_total / chunk_size)
  similarity_list <- vector("list", n_total)

  for (i in 1:n_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_total)

    similarity_list[start_idx:end_idx] <- mapply(
      calculate_similarity,
      enrichment_with_seq$actual_match[start_idx:end_idx],
      enrichment_with_seq$motif_sequence[start_idx:end_idx],
      SIMPLIFY = FALSE
    )

    if (i %% 10 == 0 || i == n_chunks) {
      cat(sprintf(" (%d/%d)", end_idx, n_total))
    }
  }
  cat(" - Done!\n")
}

# Extract metrics from list
enrichment_with_seq$pct_identity <- sapply(similarity_list, function(x) x$pct_identity)
enrichment_with_seq$num_matches <- sapply(similarity_list, function(x) x$num_matches)
enrichment_with_seq$num_positions <- sapply(similarity_list, function(x) x$num_positions)

# Remove large sequence column before saving
results <- enrichment_with_seq
results$sequence <- NULL

# Save detailed results
cat("Saving detailed results...\n")
write.csv(results, OUTPUT_FILE, row.names = FALSE)
cat(sprintf("Detailed results saved to: %s\n", OUTPUT_FILE))

# Generate summary statistics by motif
cat("\nGenerating summary statistics...\n")

# Filter out NAs
valid_results <- results[!is.na(results$pct_identity), ]

# Calculate summary by motif
motif_list <- unique(valid_results$motif)
summary_list <- lapply(motif_list, function(m) {
  subset_data <- valid_results[valid_results$motif == m, ]

  data.frame(
    motif = m,
    motif_sequence = unique(subset_data$motif_sequence)[1],
    motif_length = unique(subset_data$motif_length)[1],
    n_occurrences = nrow(subset_data),

    # Percent identity statistics
    mean_pct_identity = mean(subset_data$pct_identity, na.rm = TRUE),
    median_pct_identity = median(subset_data$pct_identity, na.rm = TRUE),
    sd_pct_identity = sd(subset_data$pct_identity, na.rm = TRUE),
    min_pct_identity = min(subset_data$pct_identity, na.rm = TRUE),
    max_pct_identity = max(subset_data$pct_identity, na.rm = TRUE),

    # Quartiles
    q25_pct_identity = quantile(subset_data$pct_identity, 0.25, na.rm = TRUE),
    q75_pct_identity = quantile(subset_data$pct_identity, 0.75, na.rm = TRUE),

    # Counts at different thresholds
    n_perfect_match = sum(subset_data$pct_identity == 100),
    n_above_90pct = sum(subset_data$pct_identity >= 90),
    n_above_80pct = sum(subset_data$pct_identity >= 80),
    n_above_70pct = sum(subset_data$pct_identity >= 70),

    # Proportions
    prop_perfect = sum(subset_data$pct_identity == 100) / nrow(subset_data),
    prop_above_90pct = sum(subset_data$pct_identity >= 90) / nrow(subset_data),
    prop_above_80pct = sum(subset_data$pct_identity >= 80) / nrow(subset_data),

    stringsAsFactors = FALSE
  )
})

summary_stats <- do.call(rbind, summary_list)
# Sort by mean_pct_identity (descending)
if (!is.null(summary_stats) && nrow(summary_stats) > 0) {
  summary_stats <- summary_stats[order(summary_stats$mean_pct_identity, decreasing = TRUE), ]
} else {
  cat("Warning: No summary statistics to report\n")
  quit(status = 1)
}

# Save summary
write.csv(summary_stats, SUMMARY_FILE, row.names = FALSE)
cat(sprintf("Summary statistics saved to: %s\n", SUMMARY_FILE))

# Print summary to console
cat("\n=== MOTIF MATCH SIMILARITY SUMMARY ===\n\n")

cat("Overall Statistics:\n")
cat(sprintf("Total matches analyzed: %d\n", nrow(valid_results)))
cat(sprintf("Mean percent identity: %.2f%%\n", mean(valid_results$pct_identity)))
cat(sprintf("Median percent identity: %.2f%%\n", median(valid_results$pct_identity)))
cat(sprintf("Perfect matches (100%%): %d (%.1f%%)\n",
            sum(valid_results$pct_identity == 100),
            100 * sum(valid_results$pct_identity == 100) / nrow(valid_results)))
cat(sprintf("Matches >= 90%% identity: %d (%.1f%%)\n",
            sum(valid_results$pct_identity >= 90),
            100 * sum(valid_results$pct_identity >= 90) / nrow(valid_results)))
cat(sprintf("Matches >= 80%% identity: %d (%.1f%%)\n",
            sum(valid_results$pct_identity >= 80),
            100 * sum(valid_results$pct_identity >= 80) / nrow(valid_results)))

cat("\n\nTop 10 Motifs by Mean Percent Identity:\n")
top10 <- head(summary_stats[, c("motif", "motif_sequence", "n_occurrences",
                                 "mean_pct_identity", "median_pct_identity",
                                 "prop_perfect")], 10)
print(top10, row.names = FALSE)

cat("\n\nBottom 10 Motifs by Mean Percent Identity:\n")
bottom10 <- tail(summary_stats[, c("motif", "motif_sequence", "n_occurrences",
                                    "mean_pct_identity", "median_pct_identity",
                                    "prop_above_80pct")], 10)
print(bottom10, row.names = FALSE)

# Distribution
cat("\n\nPercent Identity Distribution:\n")
breaks <- c(0, 50, 60, 70, 80, 90, 95, 100)
labels <- c("0-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-95%", "95-100%")
identity_bins <- cut(valid_results$pct_identity, breaks = breaks,
                     labels = labels, include.lowest = TRUE)
distribution <- table(identity_bins)
distribution_pct <- 100 * distribution / sum(distribution)

distribution_df <- data.frame(
  identity_range = names(distribution),
  count = as.numeric(distribution),
  percentage = as.numeric(distribution_pct)
)
print(distribution_df, row.names = FALSE)

cat("\n=== INTERPRETATION ===\n")
cat("The percent identity shows how similar each matched sequence is to the motif consensus.\n")
cat("Higher values indicate better matches to the consensus sequence.\n")
cat("IUPAC ambiguity codes (R, Y, S, W, etc.) count as matches if the nucleotide is one of the allowed options.\n")
cat("\nThe p-value threshold (0.0001) controls statistical significance,\n")
cat("but the actual sequences can vary from the consensus.\n")
cat("\nFor example:\n")
cat("  Consensus: WAAAAAAWATTTWW (W = A or T)\n")
cat("  Match:     AAAAAAAAATTTAA\n")
cat("  Identity:  100%% (all positions match the ambiguity code)\n")
cat("\n")

cat("\n=== FILES CREATED ===\n")
cat(sprintf("- %s (detailed results with all matches)\n", OUTPUT_FILE))
cat(sprintf("- %s (summary statistics by motif)\n", SUMMARY_FILE))
cat("\nDone!\n")
