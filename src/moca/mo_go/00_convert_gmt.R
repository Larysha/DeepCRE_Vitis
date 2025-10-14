#!/usr/bin/env Rscript
# Convert blast2go annotation file to GMT format
# GMT format: GO_ID\tGO_NAME\tGene1\tGene2\t...

library(dplyr)
library(GO.db)

# Input and output file paths
ANNOT_FILE <- "../../../vitis_data/gene_ontology/PN40024_T2T_5.1_ref_blast2go.annot"
OUTPUT_FILE <- "../../../vitis_data/gene_ontology/blast2go_t2t_5.1.gmt"

cat("Converting blast2go annotation file to GMT format\n")
cat("=" , rep("=", 70), "\n", sep = "")

# Read the ANNOT file
cat("Reading annotation file...\n")
cat("  Input:", ANNOT_FILE, "\n")

# Read with flexible column handling (2 or 3 columns)
annot_raw <- readLines(ANNOT_FILE)

# Parse each line
annot_list <- lapply(annot_raw, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 2) {
    # Return gene_id and go_id (ignore description if present)
    data.frame(
      gene_id = parts[1],
      go_id = parts[2],
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
})

# Combine into data frame
annot <- do.call(rbind, annot_list[!sapply(annot_list, is.null)])

cat("  Loaded", nrow(annot), "gene-GO term associations\n")

# Filter to only GO terms (ignore EC numbers)
annot <- annot %>%
  filter(grepl("^GO:", go_id))

cat("  Filtered to", nrow(annot), "GO term associations (removed EC numbers)\n")

# Strip the protein suffix to get gene IDs (handles all transcript variants)
cat("Converting protein IDs to gene IDs...\n")
annot$gene_id <- gsub("_t[0-9]+_CDS[0-9]+\\.prot$", "", annot$gene_id)

cat("  Unique genes:", length(unique(annot$gene_id)), "\n")
cat("  Unique GO terms:", length(unique(annot$go_id)), "\n")

# Get GO term names from GO.db (only for unique GO terms)
cat("Fetching GO term names from GO.db...\n")
unique_go_terms <- unique(annot$go_id)
cat("  Querying", length(unique_go_terms), "unique GO terms...\n")

go_name_lookup <- sapply(unique_go_terms, function(x) {
  tryCatch({
    term <- Term(GOTERM[[x]])
    if (is.null(term) || is.na(term)) {
      "Unknown"
    } else {
      term
    }
  }, error = function(e) {
    "Unknown"
  })
})

# Create lookup data frame and merge
go_lookup_df <- data.frame(
  go_id = unique_go_terms,
  go_name = go_name_lookup,
  stringsAsFactors = FALSE
)

annot <- annot %>%
  left_join(go_lookup_df, by = "go_id")

# Count how many terms have names
n_with_names <- sum(annot$go_name != "Unknown")
cat("  Found names for", n_with_names, "of", length(unique(annot$go_id)),
    "GO terms (", round(100 * n_with_names / length(unique(annot$go_id)), 1), "%)\n")

# Create GMT format: GO_ID\tGO_NAME\tGene1\tGene2\t...
cat("Creating GMT format...\n")
gmt_data <- annot %>%
  group_by(go_id, go_name) %>%
  summarise(genes = paste(unique(gene_id), collapse = "\t"), .groups = "drop") %>%
  mutate(gmt_line = paste(go_id, go_name, genes, sep = "\t"))

cat("  Total GMT entries:", nrow(gmt_data), "\n")

# Calculate statistics
gene_counts <- annot %>%
  group_by(go_id) %>%
  summarise(n_genes = n_distinct(gene_id))

cat("  Genes per GO term: min =", min(gene_counts$n_genes),
    ", median =", median(gene_counts$n_genes),
    ", max =", max(gene_counts$n_genes), "\n")

# Write GMT
cat("Writing GMT file...\n")
cat("  Output:", OUTPUT_FILE, "\n")
writeLines(gmt_data$gmt_line, OUTPUT_FILE)

cat("\nConversion complete!\n")
cat("=" , rep("=", 70), "\n", sep = "")
