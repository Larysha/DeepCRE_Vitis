# Enhanced Saliency Map Analysis Script
# Generates importance profiles and sequence logos for top-contributing genes

###########################
# 0. Load Project Utilities
###########################

# Set working directory to script location
script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
if (nchar(script_dir) > 0) {
  setwd(script_dir)
  cat("Working directory set to:", getwd(), "\n")
}

source("../utils.R")

setup_project_environment(verbose = TRUE)
required_packages <- c("rhdf5", "dplyr", "tidyr", "ggplot2", "ggseqlogo")
if (!load_required_packages(required_packages)) {
  stop("Please install missing packages before continuing.")
}

# Close any existing HDF5 connections to prevent conflicts
rhdf5::h5closeAll()

args <- commandArgs(trailingOnly = TRUE)
# Define input parameters with sensible defaults
# Use list.files() to find files matching patterns instead of wildcards
h5_file <- if (length(args) >= 1 && !is.na(args[1])) {
  args[1]
} else {
  h5_files <- list.files("../out/shap", pattern = "\\.h5$", full.names = TRUE)
  if (length(h5_files) == 0) stop("No .h5 files found in ../out/shap/")
  if (length(h5_files) > 1) {
    cat("Multiple .h5 files found, using:", h5_files[1], "\n")
  }
  h5_files[1]
}

csv_file <- if (length(args) >= 2 && !is.na(args[2])) {
  args[2]
} else {
  csv_files <- list.files("../out/predictions", pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No .csv files found in ../out/predictions/")
  if (length(csv_files) > 1) {
    cat("Multiple .csv files found, using:", csv_files[1], "\n")
  }
  csv_files[1]
}

shap_meta_file <- if (length(args) >= 3 && !is.na(args[3])) {
  args[3]
} else {
  meta_files <- list.files("../out/shap", pattern = "shap_meta\\.csv$", full.names = TRUE)
  if (length(meta_files) == 0) stop("No shap_meta.csv files found in ../out/shap/")
  if (length(meta_files) > 1) {
    cat("Multiple shap_meta.csv files found, using:", meta_files[1], "\n")
  }
  meta_files[1]
} 
output_dir_plots <- if (length(args) >= 4) args[4] else "../out/moca_results/mo_imp/figures/saliency_top5"
output_dir_scores <- if (length(args) >= 5) args[5] else "../out/moca_results/mo_imp/importance_scores"

# Define analysis window - 
# changed from original 750-1500 so that midpoint coincides with TSS at 1000
window_start <- 500
window_end <- 1500

# Create output directories
dir.create(output_dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_scores, showWarnings = FALSE, recursive = TRUE)

cat("Starting saliency analysis...\n")
cat("Plots will be saved to:", output_dir_plots, "\n")
cat("Importance scores will be saved to:", output_dir_scores, "\n")

###########################
# 1. Load and Validate Data
###########################

# Load prediction results and filter for binary classification
prediction_df <- read.csv(csv_file)
prediction_df <- prediction_df %>% 
  filter(true_targets %in% c(0, 1))

cat("Loaded", nrow(prediction_df), "predictions\n")

# Load SHAP metadata to get gene ordering
shap_meta_df <- read.csv(shap_meta_file)
shap_genes <- shap_meta_df$gene_ids

cat("Loaded metadata for", length(shap_genes), "genes with SHAP scores\n")

# Load saliency scores from HDF5
if (!file.exists(h5_file)) {
  stop("HDF5 file not found: ", h5_file)
}

h5_data <- H5Fopen(h5_file)
saliency_scores <- h5read(h5_file, "contrib_scores")

# Validate data dimensions
# Expected format: [nucleotides, positions, genes] = [4, sequence_length, num_genes]
dim_sal <- dim(saliency_scores)
stopifnot("Saliency scores should be 3-dimensional" = length(dim_sal) == 3)
stopifnot("First dimension should be 4 (A,C,G,T)" = dim_sal[1] == 4)

num_genes <- dim_sal[3]
sequence_length <- dim_sal[2]

cat("Saliency data dimensions: [", paste(dim_sal, collapse = ", "), "]\n")

# Align prediction data with SHAP gene order
prediction_df <- prediction_df %>%
  filter(genes %in% shap_genes) %>%
  slice(match(shap_genes, genes))

cat("Retained", nrow(prediction_df), "genes with binary class labels and SHAP scores\n")

# Validate that we have the expected genes column
stopifnot("Prediction data must contain 'genes' column" = "genes" %in% names(prediction_df))

###########################
# 2. Calculate Position-wise Importance Scores
###########################

# Collapse saliency scores across nucleotides for each position
# This gives the total contribution of each position to the model's prediction
cat("Calculating position-wise importance scores...\n")

collapsed_matrix <- matrix(0, nrow = num_genes, ncol = sequence_length)
for (i in seq_len(num_genes)) {
  # Sum across nucleotide dimension (A, C, G, T) for each position
  collapsed_matrix[i, ] <- colSums(saliency_scores[, , i]) # means that negative and positive contributions can cancel out
}

# Convert to data frame with gene identifiers
collapsed_df <- as.data.frame(collapsed_matrix)
colnames(collapsed_df) <- paste0("pos_", seq_len(sequence_length))
collapsed_df <- cbind(genes = shap_genes[seq_len(nrow(collapsed_df))], collapsed_df)

# Calculate total importance per gene (useful for ranking)
collapsed_df <- collapsed_df %>%
  rowwise() %>%
  mutate(total_importance = sum(c_across(starts_with("pos_")))) %>%
  ungroup()

cat("Calculated importance scores for", nrow(collapsed_df), "genes across", sequence_length, "positions\n")

###########################
# 3. Save Full Importance Matrix
###########################

# Save the complete importance matrix for downstream analysis
importance_file <- file.path(output_dir_scores, "gene_position_importance_scores.csv")
write.csv(collapsed_df, importance_file, row.names = FALSE)
cat("Saved full importance matrix to:", importance_file, "\n")

# Also save a summary statistics file
summary_df <- collapsed_df %>%
  select(genes, total_importance) %>%
  arrange(desc(total_importance)) %>%
  mutate(rank = row_number(),
         percentile = round((1 - (rank - 1) / n()) * 100, 2)) # 2 decs

summary_file <- file.path(output_dir_scores, "gene_importance_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("Saved importance summary to:", summary_file, "\n")

###########################
# 4. Identify Top Contributing Genes
###########################

# Select top 5 genes by total importance for detailed visualization
top_genes <- summary_df %>%
  slice_head(n = 5) %>%
  pull(genes)

cat("Top 5 most important genes:\n")
for (i in seq_along(top_genes)) {
  gene <- top_genes[i]
  importance <- summary_df %>% filter(genes == gene) %>% pull(total_importance)
  cat(sprintf("  %d. %s (importance: %.4f)\n", i, gene, importance))
}

###########################
# 5. Generate Visualizations for Top Genes
###########################

cat("Generating visualizations for top genes...\n")

for (gene in top_genes) {
  cat("Processing", gene, "...\n")
  
  # Extract importance scores for the analysis window
  gene_row <- collapsed_df %>% filter(genes == gene)
  
  if (nrow(gene_row) == 0) {
    warning("Gene ", gene, " not found in importance data")
    next
  }
  
  # Extract window positions (adding 1 for R's 1-based indexing)
  window_cols <- paste0("pos_", window_start:window_end)
  pos_scores <- as.numeric(gene_row[1, window_cols])
  
  # Create line plot showing importance across positions
  df_plot <- data.frame(
    position = window_start:window_end,
    importance = pos_scores
  )
  
  p_line <- ggplot(df_plot, aes(x = position, y = importance)) +
    geom_line(size = 0.3, colour = "black") +
    geom_vline(xintercept = (window_start + window_end) / 2, 
               colour = "#723983", linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dotted", alpha = 0.5) +
    labs(title = paste("Importance Profile:", gene),
         subtitle = paste("Window:", window_start, "-", window_end),
         x = "Genomic Position",
         y = "Contribution Score") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  
  # Save line plot
  line_filename <- file.path(output_dir_plots, paste0(gene, "_importance_profile.png"))
  ggsave(filename = line_filename, plot = p_line, width = 8, height = 4, dpi = 300)
  
  # Generate sequence logo showing nucleotide-specific contributions
  gene_idx <- which(shap_genes == gene)
  
  if (length(gene_idx) == 0) {
    warning("Gene ", gene, " not found in SHAP gene list")
    next
  }
  
  # Extract saliency matrix for the analysis window
  saliency_slice <- saliency_scores[, window_start:window_end, gene_idx]
  rownames(saliency_slice) <- c("A", "C", "G", "T")
  
  # Create sequence logo plot
  p_logo <- ggseqlogo(saliency_slice, method = "custom", seq_type = "dna") +
    labs(title = paste("Nucleotide Contributions:", gene),
         subtitle = paste("Window:", window_start, "-", window_end)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  
  # Save logo plot
  logo_filename <- file.path(output_dir_plots, paste0(gene, "_nucleotide_logo.png"))
  ggsave(filename = logo_filename, plot = p_logo, width = 8, height = 3, dpi = 300)
}

###########################
# 6. Cleanup and Summary
###########################

# Close HDF5 connections
rhdf5::h5closeAll()

cat("\n=== Analysis Complete ===\n")
cat("Files generated:\n")
cat("  - Full importance matrix:", importance_file, "\n")
cat("  - Importance summary:", summary_file, "\n")
cat("  - Visualizations for", length(top_genes), "genes in:", output_dir_plots, "\n")
cat("  - Total plots created:", length(top_genes) * 2, "\n")

# Display session info for reproducibility
cat("\nSession Info:\n")
sessionInfo()