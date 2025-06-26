###########################
# 0. Setup
###########################
# 0. Setup
###########################
library(rhdf5)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggseqlogo)

rhdf5::h5closeAll()
args <- commandArgs(trailingOnly = TRUE)

# Set working directory to location of the script
this_file <- normalizePath(sys.frames()[[1]]$ofile)  # works when sourcing
# If using Rscript execution, fall back to commandArgs
if (is.null(this_file)) {
  this_file <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
}
setwd(dirname(this_file))


# Input arguments with fallbacks
h5_file <- if (length(args) >= 1) args[1] else "../results/shap/vitis_02_deepcre_interpret_250610_125951.h5"
csv_file <- if (length(args) >= 2) args[2] else "../results/predictions/vitis_02_SSR_deepcre_predict_250610_101946.csv"
shap_meta_file <- if (length(args) >= 3) args[3] else "../results/shap/vitis_02_deepcre_interpret_250610_125955_shap_meta.csv" 
output_dir <- if (length(args) >= 4) args[4] else "../results/modisco/figures/saliency_top5"

window_start <- 750
window_end <- 1500

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

###########################
# load prediction data and hdf saliency scores
###########################

prediction_df <- read.csv(csv_file)
prediction_df <- prediction_df %>% 
  filter(true_targets %in% c(0, 1))

shap_meta_df <- read.csv(shap_meta_file)
shap_genes <- shap_meta_df$gene_ids  # [ADDED] assign gene list for ordering

h5_data <- H5Fopen(h5_file)
saliency_scores <- h5read(h5_file, "contrib_scores")

prediction_df <- prediction_df %>%
  filter(genes %in% shap_genes) %>%
  slice(match(shap_genes, genes))  # ensure order matches shap_meta_df

cat("Retained", nrow(prediction_df), "genes with binary class labels\n")
stopifnot("genes" %in% names(prediction_df))

# check: [A,C,G,T] x sequence_length x num_genes
dim_sal <- dim(saliency_scores)
stopifnot(length(dim_sal) == 3)

num_genes <- dim_sal[3]
sequence_length <- dim_sal[2]

###########################
# collapse importance by position
###########################
# Create a matrix: [genes x positions]
collapsed_matrix <- matrix(0, nrow = num_genes, ncol = sequence_length)
for (i in seq_len(num_genes)) {
  collapsed_matrix[i, ] <- colSums(saliency_scores[, , i])
}

# Replace prediction_df$genes with shap_genes to match ordering
collapsed_df <- as.data.frame(collapsed_matrix)
collapsed_df <- cbind(genes = shap_genes[seq_len(nrow(collapsed_df))], collapsed_df)

collapsed_df <- collapsed_df %>%
  rowwise() %>%
  mutate(sum_importance = sum(c_across(!!paste0("V1"):!!paste0("V", sequence_length)))) %>%  # [MODIFIED] dynamic column range
  ungroup()

###########################
# get top 5 genes by importance
###########################
top_genes <- collapsed_df |>
  arrange(desc(sum_importance)) |>
  slice(1:5) |>
  pull(genes)

message("Top 5 genes:\n", paste(top_genes, collapse = "\n"))

###########################
# plotting top genes
###########################

for (gene in top_genes) {
  gene_row <- collapsed_df %>% filter(genes == gene)

  if (nrow(gene_row) == 0) next

  # Line plot
  pos_scores <- as.numeric(gene_row[1, (window_start + 1):(window_end + 1)])  # +1 for data.frame indexing
  df_plot <- data.frame(position = window_start:window_end,
                        importance = pos_scores)

  p_line <- ggplot(df_plot, aes(x = position, y = importance)) +
    geom_line(size = 0.2) +
    geom_vline(xintercept = (window_start + window_end)/2, colour = "#723983", linetype = "dashed") +
    labs(title = paste("Importance profile:", gene),
         x = "Position",
         y = "Contribution Score") +
    theme_minimal(base_size = 10)

  ggsave(filename = file.path(output_dir, paste0(gene, "_lineplot.png")),
         plot = p_line, width = 6, height = 3)

  # Saliency logo plot
  gene_idx <- which(shap_genes == gene)  # [MODIFIED] use shap_genes for index lookup
  saliency_slice <- saliency_scores[, window_start:window_end, gene_idx]
  rownames(saliency_slice) <- c("A", "C", "G", "T")

  p_logo <- ggseqlogo(saliency_slice, method = "custom", seq_type = "dna") +
    labs(title = paste("Saliency logo:", gene))

  ggsave(filename = file.path(output_dir, paste0(gene, "_logo.png")),
         plot = p_logo, width = 6, height = 3)
}

cat("Plots saved in: ", output_dir, "\n")

rhdf5::h5closeAll()
