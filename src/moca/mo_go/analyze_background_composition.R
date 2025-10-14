#!/usr/bin/env Rscript

# Quick analysis: Background gene composition in GO enrichment
# Shows how many target==0 (low expression) vs target==1 (high expression) genes
# have GO annotations

library(dplyr)

# File paths (same as in 04_go_enrichment.R)
GMT_FILE <- "../../../vitis_data/gene_ontology/blast2go_t2t_5.1.gmt"
TARGET_FILE <- "../../../vitis_data/tpm_counts/vitis_drought_leaf_targets.csv"

cat("=" , rep("=", 70), "\n", sep = "")
cat("BACKGROUND GENE SET COMPOSITION ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Analyzing background genes used for GO enrichment...\n\n")

# Load target classifications
cat("Loading target classifications...\n")
targets <- read.csv(TARGET_FILE, stringsAsFactors = FALSE)

cat("Target file summary:\n")
cat(sprintf("  Total genes: %d\n", nrow(targets)))
cat(sprintf("  target == 0 (low expression): %d (%.1f%%)\n",
            sum(targets$target == 0),
            100 * sum(targets$target == 0) / nrow(targets)))
cat(sprintf("  target == 1 (high expression): %d (%.1f%%)\n",
            sum(targets$target == 1),
            100 * sum(targets$target == 1) / nrow(targets)))
cat(sprintf("  target == 2 (excluded): %d (%.1f%%)\n\n",
            sum(targets$target == 2),
            100 * sum(targets$target == 2) / nrow(targets)))

# Load GMT file to get genes with GO annotations
cat("Loading GO annotations from GMT file...\n")
gmt_lines <- readLines(GMT_FILE)

# Extract all genes from GMT
all_go_genes <- character()
for (line in gmt_lines) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    genes <- parts[3:length(parts)]
    valid_genes <- genes[nzchar(genes)]
    all_go_genes <- c(all_go_genes, valid_genes)
  }
}

# Get unique genes with GO annotations
genes_with_go <- unique(all_go_genes)
cat(sprintf("  Total unique genes with GO annotations: %d\n\n", length(genes_with_go)))

# Filter targets to only include target==0 or target==1
valid_targets <- targets %>% filter(target %in% c(0, 1))

cat("After filtering to target==0 or target==1:\n")
cat(sprintf("  Genes included: %d\n\n", nrow(valid_targets)))

# Find which of these genes have GO annotations
valid_targets$has_go <- valid_targets$gene_id %in% genes_with_go

# Overall GO annotation coverage
cat("GO Annotation Coverage:\n")
cat(sprintf("  Genes with GO annotations: %d (%.1f%%)\n",
            sum(valid_targets$has_go),
            100 * sum(valid_targets$has_go) / nrow(valid_targets)))
cat(sprintf("  Genes without GO annotations: %d (%.1f%%)\n\n",
            sum(!valid_targets$has_go),
            100 * sum(!valid_targets$has_go) / nrow(valid_targets)))

# Breakdown by target class
cat("=" , rep("=", 70), "\n", sep = "")
cat("BACKGROUND COMPOSITION BY EXPRESSION CLASS\n")
cat(rep("=", 70), "\n\n", sep = "")

# Target == 0 (low expression)
target0_genes <- valid_targets %>% filter(target == 0)
target0_with_go <- sum(target0_genes$has_go)
target0_total <- nrow(target0_genes)

cat("Target == 0 (Low Expression Genes):\n")
cat(sprintf("  Total target==0 genes: %d\n", target0_total))
cat(sprintf("  With GO annotations: %d (%.1f%% of target==0)\n",
            target0_with_go,
            100 * target0_with_go / target0_total))
cat(sprintf("  Without GO annotations: %d (%.1f%% of target==0)\n\n",
            target0_total - target0_with_go,
            100 * (target0_total - target0_with_go) / target0_total))

# Target == 1 (high expression)
target1_genes <- valid_targets %>% filter(target == 1)
target1_with_go <- sum(target1_genes$has_go)
target1_total <- nrow(target1_genes)

cat("Target == 1 (High Expression Genes):\n")
cat(sprintf("  Total target==1 genes: %d\n", target1_total))
cat(sprintf("  With GO annotations: %d (%.1f%% of target==1)\n",
            target1_with_go,
            100 * target1_with_go / target1_total))
cat(sprintf("  Without GO annotations: %d (%.1f%% of target==1)\n\n",
            target1_total - target1_with_go,
            100 * (target1_total - target1_with_go) / target1_total))

# Background composition
background_genes <- valid_targets %>% filter(has_go)
total_background <- nrow(background_genes)

cat("=" , rep("=", 70), "\n", sep = "")
cat("FINAL BACKGROUND SET COMPOSITION\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("This is the background used for GO enrichment testing:\n\n")
cat(sprintf("Total background genes: %d\n\n", total_background))

cat("Composition by expression class:\n")
cat(sprintf("  Low expression (target==0): %d genes (%.1f%% of background)\n",
            target0_with_go,
            100 * target0_with_go / total_background))
cat(sprintf("  High expression (target==1): %d genes (%.1f%% of background)\n",
            target1_with_go,
            100 * target1_with_go / total_background))

cat("\n")
cat("Visual representation:\n")
cat("  [", paste(rep("■", round(50 * target0_with_go / total_background)), collapse = ""),
    paste(rep("□", round(50 * target1_with_go / total_background)), collapse = ""), "]\n", sep = "")
cat("   ■ = target==0 (low expression)\n")
cat("   □ = target==1 (high expression)\n\n")

# Compare proportions
cat("Proportion comparison:\n")
cat(sprintf("  In full dataset (target 0 or 1): %.1f%% low, %.1f%% high\n",
            100 * target0_total / nrow(valid_targets),
            100 * target1_total / nrow(valid_targets)))
cat(sprintf("  In background (with GO): %.1f%% low, %.1f%% high\n",
            100 * target0_with_go / total_background,
            100 * target1_with_go / total_background))

# Check if proportions are similar
prop_diff <- abs((target0_with_go / total_background) - (target0_total / nrow(valid_targets)))
if (prop_diff < 0.05) {
  cat("\n✓ Background proportions are well-balanced (similar to full dataset)\n")
} else {
  cat(sprintf("\n⚠ Background is slightly skewed (%.1f%% difference from full dataset)\n",
              100 * prop_diff))
}

cat("\n")
cat("=" , rep("=", 70), "\n", sep = "")
cat("INTERPRETATION\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("The background set contains genes that:\n")
cat("  1. Have target == 0 (low expression) OR target == 1 (high expression)\n")
cat("  2. Have GO annotations in the GMT file\n")
cat("  3. Are used as the reference for enrichment testing\n\n")

cat("When testing a gene list for GO enrichment:\n")
cat("  - Only genes in the background are considered\n")
cat("  - Enrichment is relative to this background composition\n")
cat("  - Target==2 genes are excluded from both gene lists and background\n\n")

cat("Analysis complete!\n")
