# Quick analysis to check if filtering rate (from filter_balmm_occ.R) matches gene length distribution

library(dplyr)

# Load your gene annotations
gff_data <- read.table("../../../vitis_data/gene_models/vitis_vinifera_PN40024.gff3", 
                       header = FALSE, sep = "\t", comment.char = "#", 
                       stringsAsFactors = FALSE)

colnames(gff_data) <- c("seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute")

# Calculate gene lengths and expected retention rates
gene_stats <- gff_data %>%
  filter(feature == "gene") %>%
  mutate(
    gene_length = end - start + 1,
    # Calculate how much of each gene would be retained
    retained_bp = pmin(gene_length, 1000),  # Max 1000bp from UTRs (500 each end)
    retention_rate = retained_bp / (gene_length + 2000),  # +2000 for flanks
    total_extracted = gene_length + 2000  # Total extracted sequence length
  )

# Summary statistics
cat("Gene Length Analysis:\n")
cat("Median gene length:", median(gene_stats$gene_length), "bp\n")
cat("Mean gene length:", round(mean(gene_stats$gene_length)), "bp\n")
cat("Expected retention rate:", round(mean(gene_stats$retention_rate) * 100, 1), "%\n")

# Length distribution
length_summary <- gene_stats %>%
  summarise(
    short_genes = sum(gene_length <= 1000),
    medium_genes = sum(gene_length > 1000 & gene_length <= 3000),
    long_genes = sum(gene_length > 3000),
    total_genes = n()
  )

cat("\nGene length distribution:\n")
cat("Short genes (≤1000bp):", length_summary$short_genes, "\n")
cat("Medium genes (1000-3000bp):", length_summary$medium_genes, "\n") 
cat("Long genes (>3000bp):", length_summary$long_genes, "\n")