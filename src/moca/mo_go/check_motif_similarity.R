#!/usr/bin/env Rscript
# Check if motifs from P0 (metacluster 0) match motifs from P1 (metacluster 1)
# This includes checking reverse complements

library(dplyr)
library(readr)

# Read the motif summary file
motif_file <- "../../../out/moca_results/mo_clu/20250918_vitis_ssr_cwm_clustering_motif_summary.csv"
motifs <- read_csv(motif_file, show_col_types = FALSE)

# Function to get reverse complement
reverse_complement <- function(seq) {
  # Define complement mapping
  complement_map <- c(
    "A" = "T", "T" = "A", "G" = "C", "C" = "G",
    "R" = "Y", "Y" = "R",  # R=A/G, Y=C/T
    "S" = "S",  # S=G/C (self-complement)
    "W" = "W",  # W=A/T (self-complement)
    "K" = "M", "M" = "K",  # K=G/T, M=A/C
    "B" = "V", "V" = "B",  # B=C/G/T, V=A/C/G
    "D" = "H", "H" = "D",  # D=A/G/T, H=A/C/T
    "N" = "N"   # N=any
  )

  # Split sequence, complement, and reverse
  seq_split <- strsplit(seq, "")[[1]]
  complemented <- sapply(seq_split, function(x) complement_map[x])
  rev_comp <- paste(rev(complemented), collapse = "")
  return(rev_comp)
}

# Parse motif names to extract pattern, metacluster, strand
parse_motif <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 7) {
    pattern <- parts[4]
    metacluster <- parts[5]
    strand <- parts[6]
    return(list(pattern = pattern, metacluster = metacluster, strand = strand))
  }
  return(NULL)
}

# Add parsed information
motifs$pattern <- sapply(motifs$name, function(x) {
  parsed <- parse_motif(x)
  if (!is.null(parsed)) parsed$pattern else NA
})

motifs$metacluster <- sapply(motifs$name, function(x) {
  parsed <- parse_motif(x)
  if (!is.null(parsed)) parsed$metacluster else NA
})

motifs$strand <- sapply(motifs$name, function(x) {
  parsed <- parse_motif(x)
  if (!is.null(parsed)) parsed$strand else NA
})

# Keep only forward strand for comparison (to avoid double-counting F and R versions)
motifs_forward <- motifs %>% filter(strand == "F")

# Separate P0 and P1 motifs
p0_motifs <- motifs_forward %>% filter(metacluster == "0")
p1_motifs <- motifs_forward %>% filter(metacluster == "1")

cat("=== Motif Similarity Analysis ===\n\n")
cat("P0 motifs (metacluster 0):", nrow(p0_motifs), "\n")
cat("P1 motifs (metacluster 1):", nrow(p1_motifs), "\n\n")

# Check for exact matches (forward)
cat("=== Checking for exact consensus matches ===\n")
exact_matches <- data.frame()

for (i in 1:nrow(p0_motifs)) {
  p0_consensus <- p0_motifs$consensus[i]
  p0_name <- p0_motifs$name[i]

  for (j in 1:nrow(p1_motifs)) {
    p1_consensus <- p1_motifs$consensus[j]
    p1_name <- p1_motifs$name[j]

    if (p0_consensus == p1_consensus) {
      exact_matches <- rbind(exact_matches, data.frame(
        p0_motif = p0_name,
        p1_motif = p1_name,
        consensus = p0_consensus,
        match_type = "exact",
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(exact_matches) > 0) {
  cat("\nFound", nrow(exact_matches), "exact matches:\n")
  print(exact_matches)
} else {
  cat("\nNo exact matches found between P0 and P1 motifs.\n")
}

# Check for reverse complement matches
cat("\n=== Checking for reverse complement matches ===\n")
rc_matches <- data.frame()

for (i in 1:nrow(p0_motifs)) {
  p0_consensus <- p0_motifs$consensus[i]
  p0_rc <- reverse_complement(p0_consensus)
  p0_name <- p0_motifs$name[i]

  for (j in 1:nrow(p1_motifs)) {
    p1_consensus <- p1_motifs$consensus[j]
    p1_name <- p1_motifs$name[j]

    if (p0_rc == p1_consensus) {
      rc_matches <- rbind(rc_matches, data.frame(
        p0_motif = p0_name,
        p1_motif = p1_name,
        p0_consensus = p0_consensus,
        p1_consensus = p1_consensus,
        match_type = "reverse_complement",
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(rc_matches) > 0) {
  cat("\nFound", nrow(rc_matches), "reverse complement matches:\n")
  print(rc_matches)
} else {
  cat("\nNo reverse complement matches found between P0 and P1 motifs.\n")
}

# Calculate similarity scores (simple: count matching positions including IUPAC)
calculate_similarity <- function(seq1, seq2) {
  if (nchar(seq1) != nchar(seq2)) return(0)

  seq1_split <- strsplit(seq1, "")[[1]]
  seq2_split <- strsplit(seq2, "")[[1]]

  matches <- sum(seq1_split == seq2_split)
  return(matches / nchar(seq1))
}

# Find high-similarity pairs (>80% similarity)
cat("\n=== Checking for high-similarity matches (>80%) ===\n")
high_sim_matches <- data.frame()

for (i in 1:nrow(p0_motifs)) {
  p0_consensus <- p0_motifs$consensus[i]
  p0_name <- p0_motifs$name[i]

  for (j in 1:nrow(p1_motifs)) {
    p1_consensus <- p1_motifs$consensus[j]
    p1_name <- p1_motifs$name[j]

    if (nchar(p0_consensus) == nchar(p1_consensus)) {
      sim <- calculate_similarity(p0_consensus, p1_consensus)

      if (sim > 0.8 && sim < 1.0) {  # Exclude exact matches already found
        high_sim_matches <- rbind(high_sim_matches, data.frame(
          p0_motif = p0_name,
          p1_motif = p1_name,
          p0_consensus = p0_consensus,
          p1_consensus = p1_consensus,
          similarity = round(sim, 3),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

if (nrow(high_sim_matches) > 0) {
  cat("\nFound", nrow(high_sim_matches), "high-similarity matches:\n")
  print(high_sim_matches)
} else {
  cat("\nNo high-similarity matches (>80%) found between P0 and P1 motifs.\n")
}

# Summary
cat("\n=== SUMMARY ===\n")
total_cross_matches <- nrow(exact_matches) + nrow(rc_matches) + nrow(high_sim_matches)

if (total_cross_matches == 0) {
  cat("\n✓ All motifs are UNIQUE and DISTINCT across P0 and P1 metaclusters.\n")
  cat("  No exact matches, reverse complements, or high-similarity matches found.\n")
  cat("  Each motif represents a distinct regulatory sequence specific to its metacluster.\n")
} else {
  cat("\nFound", total_cross_matches, "potential cross-metacluster matches:\n")
  cat("  - Exact matches:", nrow(exact_matches), "\n")
  cat("  - Reverse complement matches:", nrow(rc_matches), "\n")
  cat("  - High similarity (>80%):", nrow(high_sim_matches), "\n")
}

cat("\n")
