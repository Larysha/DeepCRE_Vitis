#!/usr/bin/env Rscript
# Professional Dendrogram Visualization with Clade Highlighting
# Generates publication-quality dendrograms with colored rectangles around specified clades

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)

# Color palette definitions
COLORS <- list(
  grey = "#CCCCCC",
  red = "#E74C3C",
  yellow = "#F39C12",
  green = "#27AE60",
  pink = "#E91E63",
  blue = "#3498DB",
  lightgreen = "#82E0AA",
  darkgreen = "#196F3D",
  darkblue = "#1F618D"
)

#' Get MRCA (Most Recent Common Ancestor) node for a set of tips
#'
#' @param tree phylo object
#' @param tip_labels character vector of tip labels
#' @return node number of MRCA
get_clade_mrca <- function(tree, tip_labels) {
  tip_indices <- which(tree$tip.label %in% tip_labels)
  if (length(tip_indices) < 2) {
    return(tip_indices[1])
  }
  return(getMRCA(tree, tip_indices))
}

#' Create highlighted dendrogram with colored rectangles
#'
#' @param tree_file path to .nwk file
#' @param clade_definitions list of lists with name, tips, color, alpha
#' @param output_file path for output PDF
#' @param title plot title
create_highlighted_dendrogram <- function(tree_file, clade_definitions, output_file, title) {

  cat("Processing:", tree_file, "\n")

  # Read tree
  tree <- read.tree(tree_file)

  # Create base ggtree plot
  p <- ggtree(tree, ladderize = TRUE, branch.length = "branch.length") +
    theme_tree2() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    ggtitle(title)

  # Add clade highlights in reverse order (so they layer correctly)
  for (clade in rev(clade_definitions)) {
    mrca_node <- get_clade_mrca(tree, clade$tips)

    if (!is.na(mrca_node) && mrca_node > 0) {
      p <- p + geom_hilight(
        node = mrca_node,
        fill = clade$color,
        alpha = clade$alpha,
        extend = 0.01
      )
    } else {
      cat("  Warning: Could not find MRCA for clade:", clade$name, "\n")
    }
  }

  # Add tip labels
  p <- p + geom_tiplab(size = 3, hjust = -0.1) +
    xlim(NA, max(node.depth.edgelength(tree)) * 1.3)

  # Save to PDF
  cat("  Saving to:", output_file, "\n")
  ggsave(
    filename = output_file,
    plot = p,
    width = 10,
    height = 8,
    units = "in",
    device = "pdf"
  )

  cat("  ✓ Completed\n\n")
}

# Define p0m clade highlights (back to front order)
p0m_clades <- list(
  list(
    name = "p0m_l0_m8_grey",
    tips = c("epm_vitis_ssr_p0m03", "epm_vitis_ssr_p0m11", "epm_vitis_ssr_p0m07",
             "epm_vitis_ssr_p0m01", "epm_vitis_ssr_p0m08", "epm_vitis_ssr_p0m10",
             "epm_vitis_ssr_p0m00", "epm_vitis_ssr_p0m09"),
    color = COLORS$grey,
    alpha = 0.3
  ),
  list(
    name = "p0m_l0_m4_red",
    tips = c("epm_vitis_ssr_p0m06", "epm_vitis_ssr_p0m04",
             "epm_vitis_ssr_p0m02", "epm_vitis_ssr_p0m05"),
    color = COLORS$red,
    alpha = 0.3
  ),
  list(
    name = "p0m_l1_m2_yellow",
    tips = c("epm_vitis_ssr_p0m03", "epm_vitis_ssr_p0m11"),
    color = COLORS$yellow,
    alpha = 0.3
  ),
  list(
    name = "p0m_l1_m3_green",
    tips = c("epm_vitis_ssr_p0m07", "epm_vitis_ssr_p0m01", "epm_vitis_ssr_p0m08"),
    color = COLORS$green,
    alpha = 0.3
  ),
  list(
    name = "p0m_l1_m3_pink",
    tips = c("epm_vitis_ssr_p0m10", "epm_vitis_ssr_p0m00", "epm_vitis_ssr_p0m09"),
    color = COLORS$pink,
    alpha = 0.3
  )
)

# Define p1m clade highlights (back to front order)
p1m_clades <- list(
  list(
    name = "p1m_l0_m11_red",
    tips = c("epm_vitis_ssr_p1m05", "epm_vitis_ssr_p1m13", "epm_vitis_ssr_p1m06",
             "epm_vitis_ssr_p1m10", "epm_vitis_ssr_p1m00", "epm_vitis_ssr_p1m15",
             "epm_vitis_ssr_p1m01", "epm_vitis_ssr_p1m12", "epm_vitis_ssr_p1m07",
             "epm_vitis_ssr_p1m08", "epm_vitis_ssr_p1m16"),
    color = COLORS$red,
    alpha = 0.3
  ),
  list(
    name = "p1m_l0_m6_blue",
    tips = c("epm_vitis_ssr_p1m11", "epm_vitis_ssr_p1m14", "epm_vitis_ssr_p1m09",
             "epm_vitis_ssr_p1m04", "epm_vitis_ssr_p1m02", "epm_vitis_ssr_p1m03"),
    color = COLORS$blue,
    alpha = 0.3
  ),
  list(
    name = "p1m_l1_m6_pink",
    tips = c("epm_vitis_ssr_p1m05", "epm_vitis_ssr_p1m13", "epm_vitis_ssr_p1m06",
             "epm_vitis_ssr_p1m10", "epm_vitis_ssr_p1m00", "epm_vitis_ssr_p1m15"),
    color = COLORS$pink,
    alpha = 0.3
  ),
  list(
    name = "p1m_l1_m5_lightgreen",
    tips = c("epm_vitis_ssr_p1m01", "epm_vitis_ssr_p1m12", "epm_vitis_ssr_p1m07",
             "epm_vitis_ssr_p1m08", "epm_vitis_ssr_p1m16"),
    color = COLORS$lightgreen,
    alpha = 0.3
  ),
  list(
    name = "p1m_l2_m4_darkgreen",
    tips = c("epm_vitis_ssr_p1m09", "epm_vitis_ssr_p1m04",
             "epm_vitis_ssr_p1m02", "epm_vitis_ssr_p1m03"),
    color = COLORS$darkgreen,
    alpha = 0.3
  ),
  list(
    name = "p1m_l2_m2_darkblue",
    tips = c("epm_vitis_ssr_p1m06", "epm_vitis_ssr_p1m10"),
    color = COLORS$darkblue,
    alpha = 0.3
  )
)

# File paths
base_dir <- "/home/rish/rish_phd/deepcre_vitis/vitis_cre/out/moca_results/mo_go/motif_gene_mapping"
output_dir <- file.path(base_dir, "motif_overlap_plots")

p0m_tree_file <- file.path(base_dir, "p0m_motifs/dendrograms/p0m_motif_dendrogram.nwk")
p1m_tree_file <- file.path(base_dir, "p1m_motifs/dendrograms/p1m_motif_dendrogram.nwk")

p0m_output <- file.path(output_dir, "p0m_motif_dendrogram_highlighted.pdf")
p1m_output <- file.path(output_dir, "p1m_motif_dendrogram_highlighted.pdf")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n\n")
}

# Generate dendrograms
cat(rep("=", 60), "\n")
cat("Generating Highlighted Dendrograms\n")
cat(rep("=", 60), "\n\n")

create_highlighted_dendrogram(
  tree_file = p0m_tree_file,
  clade_definitions = p0m_clades,
  output_file = p0m_output,
  title = "P0M Motif Dendrogram (High Expression Associated)"
)

create_highlighted_dendrogram(
  tree_file = p1m_tree_file,
  clade_definitions = p1m_clades,
  output_file = p1m_output,
  title = "P1M Motif Dendrogram (Low Expression Associated)"
)

cat(rep("=", 60), "\n")
cat("All dendrograms generated successfully!\n")
cat("Output location:", output_dir, "\n")
cat(rep("=", 60), "\n")
