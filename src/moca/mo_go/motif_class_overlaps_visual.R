# Motif-Clade Heatmap Visualisation
# Author: Larysha
# Description: Heatmap showing motif presence across P0/P1 clades with hierarchical clustering

# Load libraries
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Define motif membership per clade (hardcoded)
clades <- list(
  # P0 clades
  p0m_l0_m4_red = c("epm_vitis_ssr_p0m06", "epm_vitis_ssr_p0m04", 
                     "epm_vitis_ssr_p0m02", "epm_vitis_ssr_p0m05"),
  
  p0m_l1_m2_yellow = c("epm_vitis_ssr_p0m03", "epm_vitis_ssr_p0m11"),
  
  p0m_l1_m3_green = c("epm_vitis_ssr_p0m07", "epm_vitis_ssr_p0m01", 
                       "epm_vitis_ssr_p0m08"),
  
  p0m_l1_m3_pink = c("epm_vitis_ssr_p0m10", "epm_vitis_ssr_p0m00", 
                      "epm_vitis_ssr_p0m09"),
  
  # P1 clades
  p1m_l0_m6_blue = c("epm_vitis_ssr_p1m11", "epm_vitis_ssr_p1m14", 
                      "epm_vitis_ssr_p1m09", "epm_vitis_ssr_p1m04", 
                      "epm_vitis_ssr_p1m02", "epm_vitis_ssr_p1m03"),
  
  p1m_l2_m4_darkgreen = c("epm_vitis_ssr_p1m09", "epm_vitis_ssr_p1m04", 
                           "epm_vitis_ssr_p1m02", "epm_vitis_ssr_p1m03"),
  
  p1m_l0_m11_red = c("epm_vitis_ssr_p1m05", "epm_vitis_ssr_p1m13", 
                      "epm_vitis_ssr_p1m06", "epm_vitis_ssr_p1m10", 
                      "epm_vitis_ssr_p1m00", "epm_vitis_ssr_p1m15", 
                      "epm_vitis_ssr_p1m01", "epm_vitis_ssr_p1m12", 
                      "epm_vitis_ssr_p1m07", "epm_vitis_ssr_p1m08", 
                      "epm_vitis_ssr_p1m16"),
  
  p1m_l1_m6_pink = c("epm_vitis_ssr_p1m05", "epm_vitis_ssr_p1m13", 
                      "epm_vitis_ssr_p1m06", "epm_vitis_ssr_p1m10", 
                      "epm_vitis_ssr_p1m00", "epm_vitis_ssr_p1m15"),
  
  p1m_l1_m5_lightgreen = c("epm_vitis_ssr_p1m01", "epm_vitis_ssr_p1m12", 
                            "epm_vitis_ssr_p1m07", "epm_vitis_ssr_p1m08", 
                            "epm_vitis_ssr_p1m16"),
  
  p1m_l2_m2_darkblue = c("epm_vitis_ssr_p1m06", "epm_vitis_ssr_p1m10")
)

# Create output directory
dir.create("motif_overlap_plots", showWarnings = FALSE)

# ============================================================================
# Summary table: motif occurrence across clades
# ============================================================================

# Get all unique motifs
all_motifs <- unique(unlist(clades))

# Create occurrence matrix
occurrence_matrix <- data.frame(
  motif = all_motifs,
  stringsAsFactors = FALSE
)

for (clade_name in names(clades)) {
  occurrence_matrix[[clade_name]] <- as.integer(occurrence_matrix$motif %in% clades[[clade_name]])
}

# Add summary columns
occurrence_matrix$total_clades <- rowSums(occurrence_matrix[, -1])
occurrence_matrix$in_p0 <- rowSums(occurrence_matrix[, grepl("^p0m", names(occurrence_matrix))])
occurrence_matrix$in_p1 <- rowSums(occurrence_matrix[, grepl("^p1m", names(occurrence_matrix))])
occurrence_matrix$cross_boundary <- (occurrence_matrix$in_p0 > 0) & (occurrence_matrix$in_p1 > 0)

# Sort by total occurrences
occurrence_matrix <- occurrence_matrix[order(-occurrence_matrix$total_clades), ]

# Save summary table
write.csv(occurrence_matrix, 
          "motif_overlap_plots/motif_occurrence_summary.csv", 
          row.names = FALSE)

# Print summary statistics
cat("\nTotal unique motifs:", nrow(occurrence_matrix), "\n")
cat("Motifs appearing in P0 clades:", sum(occurrence_matrix$in_p0 > 0), "\n")
cat("Motifs appearing in P1 clades:", sum(occurrence_matrix$in_p1 > 0), "\n")
cat("Motifs crossing P0/P1 boundary:", sum(occurrence_matrix$cross_boundary), "\n\n")

# ============================================================================
# Heatmap of motif x clade occurrence
# ============================================================================

# Prepare matrix for heatmap (motifs in rows, clades in columns)
heatmap_matrix <- as.matrix(occurrence_matrix[, grepl("^p[01]m", names(occurrence_matrix))])
rownames(heatmap_matrix) <- occurrence_matrix$motif

# Simplify motif names by removing "epm_vitis_ssr_" prefix
simplified_names <- gsub("epm_vitis_ssr_", "", rownames(heatmap_matrix))
rownames(heatmap_matrix) <- simplified_names

# Reorder columns: P0 clades first, then P1
p0_cols <- grep("^p0m", colnames(heatmap_matrix))
p1_cols <- grep("^p1m", colnames(heatmap_matrix))
heatmap_matrix <- heatmap_matrix[, c(p0_cols, p1_cols)]

# Show ALL motifs (remove the filtering that was previously applied)
# No filtering - include all motifs

# Create annotation for column colors (clades - P0 vs P1)
col_annotation <- data.frame(
  Clade = ifelse(grepl("^p0m", colnames(heatmap_matrix)), "P0", "P1"),
  row.names = colnames(heatmap_matrix)
)

# Create annotation for row colors (motifs - P0 vs P1)
row_annotation <- data.frame(
  Motif = ifelse(grepl("^p0m", rownames(heatmap_matrix)), "P0", "P1"),
  row.names = rownames(heatmap_matrix)
)

# Define colors for annotations
annotation_colors <- list(
  Clade = c(P0 = "#052b67", P1 = "#176272"),
  Motif = c(P0 = "#052b67", P1 = "#176272")
)

# Create a modified matrix where we color based on motif phase
# For pheatmap, we need to create a matrix with values that will map to our colors
# 0 = absent (white), 1 = P0 motif present (blue), 2 = P1 motif present (teal)
colored_matrix <- heatmap_matrix
for (i in 1:nrow(colored_matrix)) {
  motif_name <- rownames(colored_matrix)[i]
  if (grepl("^p0m", motif_name)) {
    # P0 motif: set presence to 1
    colored_matrix[i, ] <- ifelse(colored_matrix[i, ] == 1, 1, 0)
  } else {
    # P1 motif: set presence to 2
    colored_matrix[i, ] <- ifelse(colored_matrix[i, ] == 1, 2, 0)
  }
}

# Create color palette: white (0), P0 blue (1), P1 teal (2)
color_palette <- c("white", "#052b67", "#176272")

# Create the heatmap
pdf("motif_overlap_plots/motif_clade_heatmap.pdf", width = 14, height = 12)
pheatmap(
  colored_matrix,
  color = color_palette,
  breaks = c(-0.5, 0.5, 1.5, 2.5),  # Define breaks for 0, 1, 2
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = col_annotation,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  annotation_legend = TRUE,
  legend = FALSE,
  main = "Motif Occurrence Across P0 and P1 Clades",
  fontsize = 9,
  fontsize_row = 8,
  fontsize_col = 9,
  border_color = "grey70",
  angle_col = 45
)
dev.off()

cat("\nHeatmap saved to: motif_overlap_plots/motif_clade_heatmap.pdf\n")
cat("Summary table saved to: motif_overlap_plots/motif_occurrence_summary.csv\n")