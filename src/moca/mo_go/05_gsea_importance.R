#!/usr/bin/env Rscript
# Gene Set Enrichment Analysis (GSEA) using SHAP importance scores
# Uses clusterProfiler's gseGO() to identify GO terms enriched in high-importance genes
#
# Usage:
#   Rscript 06_gsea_importance.R [importance_file] [gmt_file] [pval] [qval] [minGS] [maxGS]
#
# Example:
#   Rscript 06_gsea_importance.R ../../../out/moca_results/mo_imp/importance_scores/gene_importance_summary.csv
######################

library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(GO.db)
library(DOSE)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parameters
IMPORTANCE_FILE <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  "../../../out/moca_results/mo_imp/importance_scores/gene_importance_summary.csv"
}

GMT_FILE <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  "../../../vitis_data/gene_ontology/blast2go_t2t_5.1.gmt"
}

P_VALUE_CUTOFF <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
Q_VALUE_CUTOFF <- if (length(args) >= 4) as.numeric(args[4]) else 0.2
MIN_GENE_SIZE <- if (length(args) >= 5) as.integer(args[5]) else 3
MAX_GENE_SIZE <- if (length(args) >= 6) as.integer(args[6]) else 500

# Output directory
OUTPUT_DIR <- "../../../out/moca_results/mo_go/gsea"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

cat("=" , rep("=", 70), "\n", sep = "")
cat("GENE SET ENRICHMENT ANALYSIS (GSEA)\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Parameters:\n")
cat("  Importance file:", IMPORTANCE_FILE, "\n")
cat("  GMT file:", GMT_FILE, "\n")
cat("  P-value cutoff:", P_VALUE_CUTOFF, "\n")
cat("  Q-value (FDR) cutoff:", Q_VALUE_CUTOFF, "\n")
cat("  Min gene set size:", MIN_GENE_SIZE, "\n")
cat("  Max gene set size:", MAX_GENE_SIZE, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n\n")

# Check files exist
if (!file.exists(IMPORTANCE_FILE)) {
  stop("Importance file not found: ", IMPORTANCE_FILE)
}

if (!file.exists(GMT_FILE)) {
  stop("GMT file not found: ", GMT_FILE)
}

#' Read gene importance scores and create ranked list
#'
#' @param file_path Path to gene importance CSV
#' @return Named numeric vector (gene IDs as names, importance as values)
read_gene_importance <- function(file_path) {

  cat("Loading gene importance scores...\n")

  # Read data
  importance_data <- read.csv(file_path, stringsAsFactors = FALSE)

  cat("  Total genes:", nrow(importance_data), "\n")

  # Check required columns
  if (!"genes" %in% colnames(importance_data)) {
    stop("Missing 'genes' column in importance file")
  }

  if (!"total_importance" %in% colnames(importance_data)) {
    stop("Missing 'total_importance' column in importance file")
  }

  # Create named vector
  gene_scores <- importance_data$total_importance
  names(gene_scores) <- importance_data$genes

  # Sort by importance (descending)
  gene_scores <- sort(gene_scores, decreasing = TRUE)

  cat("  Importance score range: [", min(gene_scores), ", ", max(gene_scores), "]\n")
  cat("  Top 5 genes:\n")
  for (i in 1:min(5, length(gene_scores))) {
    cat("    ", names(gene_scores)[i], ": ", gene_scores[i], "\n", sep = "")
  }
  cat("\n")

  return(gene_scores)
}

#' Read GMT file and create TERM2GENE and TERM2NAME data frames, split by ontology
#'
#' @param gmt_file Path to GMT file
#' @return List with term2gene, term2name, and data split by ontology (BP, MF, CC)
read_gmt_for_gsea <- function(gmt_file) {

  cat("Loading GO annotations from GMT file...\n")

  gmt_lines <- readLines(gmt_file)

  term2gene <- data.frame()
  term2name <- data.frame()

  for (line in gmt_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      term_id <- parts[1]
      term_name <- parts[2]
      genes <- parts[3:length(parts)]

      # Filter non-empty genes
      valid_genes <- genes[nzchar(genes)]

      if (length(valid_genes) > 0) {
        # Create term2gene mapping
        term2gene <- rbind(term2gene, data.frame(
          term = term_id,
          gene = valid_genes,
          stringsAsFactors = FALSE
        ))

        # Create term2name mapping (only once per term)
        if (!term_id %in% term2name$term) {
          term2name <- rbind(term2name, data.frame(
            term = term_id,
            name = term_name,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  # Get ontology for each GO term
  cat("  Annotating terms with ontology information...\n")
  term2name$ontology <- sapply(term2name$term, function(go_id) {
    tryCatch({
      Ontology(GOTERM[[go_id]])
    }, error = function(e) "Unknown")
  })

  # Split by ontology
  term2gene_bp <- term2gene %>% filter(term %in% term2name$term[term2name$ontology == "BP"])
  term2gene_mf <- term2gene %>% filter(term %in% term2name$term[term2name$ontology == "MF"])
  term2gene_cc <- term2gene %>% filter(term %in% term2name$term[term2name$ontology == "CC"])

  term2name_bp <- term2name %>% filter(ontology == "BP")
  term2name_mf <- term2name %>% filter(ontology == "MF")
  term2name_cc <- term2name %>% filter(ontology == "CC")

  cat("  Loaded", length(unique(term2gene$term)), "GO terms\n")
  cat("    BP:", nrow(term2name_bp), "terms\n")
  cat("    MF:", nrow(term2name_mf), "terms\n")
  cat("    CC:", nrow(term2name_cc), "terms\n")
  cat("  Total gene-term associations:", nrow(term2gene), "\n")
  cat("  Genes with GO annotation:", length(unique(term2gene$gene)), "\n\n")

  return(list(
    term2gene = term2gene,
    term2name = term2name,
    bp = list(term2gene = term2gene_bp, term2name = term2name_bp),
    mf = list(term2gene = term2gene_mf, term2name = term2name_mf),
    cc = list(term2gene = term2gene_cc, term2name = term2name_cc)
  ))
}

#' Perform GSEA for a single ontology
#'
#' @param gene_list Named vector of gene importance scores
#' @param term2gene TERM2GENE data frame
#' @param term2name TERM2NAME data frame
#' @param ontology Ontology name (BP, MF, or CC)
#' @return GSEA result object or NULL
perform_gsea_ontology <- function(gene_list, term2gene, term2name, ontology) {

  cat("  Running GSEA for", ontology, "...\n")

  if (nrow(term2gene) == 0) {
    cat("    No terms available\n")
    return(NULL)
  }

  # Run GSEA
  gsea_result <- tryCatch({
    GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      TERM2NAME = term2name[, c("term", "name")],
      pvalueCutoff = P_VALUE_CUTOFF,
      pAdjustMethod = "BH",
      minGSSize = MIN_GENE_SIZE,
      maxGSSize = MAX_GENE_SIZE,
      by = "fgsea"  # Fast GSEA algorithm
    )
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
    return(NULL)
  })

  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    cat("    No significant enrichment\n")
    return(NULL)
  }

  # Filter by q-value
  gsea_result@result <- gsea_result@result %>%
    filter(qvalue < Q_VALUE_CUTOFF)

  if (nrow(gsea_result@result) == 0) {
    cat("    No terms pass FDR threshold\n")
    return(NULL)
  }

  n_pos <- sum(gsea_result@result$NES > 0)
  n_neg <- sum(gsea_result@result$NES < 0)

  cat("    Found", nrow(gsea_result@result), "enriched terms\n")
  cat("      Positive NES (enriched in high importance):", n_pos, "\n")
  cat("      Negative NES (enriched in low importance):", n_neg, "\n")

  return(gsea_result)
}

#' Run GSEA for all ontologies
#'
#' @param gene_list Named vector of gene importance scores
#' @param gmt_data GMT data with term2gene/term2name split by ontology
#' @return List with GSEA results for BP, MF, CC
run_gsea_all_ontologies <- function(gene_list, gmt_data) {

  cat("Performing GSEA for all ontologies...\n")

  results <- list()

  # Biological Process
  results$bp <- perform_gsea_ontology(
    gene_list,
    gmt_data$bp$term2gene,
    gmt_data$bp$term2name,
    "BP"
  )

  # Molecular Function
  results$mf <- perform_gsea_ontology(
    gene_list,
    gmt_data$mf$term2gene,
    gmt_data$mf$term2name,
    "MF"
  )

  # Cellular Component
  results$cc <- perform_gsea_ontology(
    gene_list,
    gmt_data$cc$term2gene,
    gmt_data$cc$term2name,
    "CC"
  )

  cat("\n")

  return(results)
}

#' Save GSEA results to CSV files
#'
#' @param gsea_results List of GSEA results
#' @param output_dir Output directory
save_gsea_results <- function(gsea_results, output_dir) {

  cat("Saving GSEA results...\n")

  # Save individual ontology results
  ontologies <- c("bp" = "BP", "mf" = "MF", "cc" = "CC")

  all_results <- list()

  for (ont_code in names(ontologies)) {
    ont_name <- ontologies[ont_code]

    if (!is.null(gsea_results[[ont_code]])) {
      result_df <- as.data.frame(gsea_results[[ont_code]]@result)
      result_df$Ontology <- ont_name

      # Save individual file
      out_file <- file.path(output_dir, paste0("gsea_", ont_name, ".csv"))
      write.csv(result_df, out_file, row.names = FALSE)
      cat("  Saved", ont_name, "results:", basename(out_file), "(",
          nrow(result_df), "terms )\n")

      all_results[[ont_code]] <- result_df
    }
  }

  # Save combined results
  if (length(all_results) > 0) {
    combined_df <- bind_rows(all_results)
    combined_file <- file.path(output_dir, "gsea_all.csv")
    write.csv(combined_df, combined_file, row.names = FALSE)
    cat("  Saved combined results:", basename(combined_file), "(",
        nrow(combined_df), "terms )\n")
  }

  cat("\n")
}

#' Create dotplots for GSEA results
#'
#' @param gsea_results List of GSEA results
#' @param output_dir Output directory
create_gsea_dotplots <- function(gsea_results, output_dir) {

  cat("Creating GSEA dotplots...\n")

  ontology_colors <- c(
    "BP" = "#671436",
    "MF" = "#052b67",
    "CC" = "#176272"
  )

  # Combine all results for faceted plot
  all_results <- list()

  for (ont_code in c("bp", "mf", "cc")) {
    if (!is.null(gsea_results[[ont_code]])) {
      result_df <- as.data.frame(gsea_results[[ont_code]]@result)
      result_df$Ontology <- toupper(ont_code)
      all_results[[ont_code]] <- result_df
    }
  }

  if (length(all_results) == 0) {
    cat("  No results to plot\n")
    return()
  }

  combined_df <- bind_rows(all_results)

  # Select top 30 terms per ontology
  top_terms <- combined_df %>%
    group_by(Ontology) %>%
    arrange(pvalue) %>%
    slice_head(n = 30) %>%
    ungroup()

  if (nrow(top_terms) == 0) {
    cat("  No terms to plot\n")
    return()
  }

  # Create dotplot
  dotplot_file <- file.path(output_dir, "gsea_dotplot.pdf")

  # Calculate height based on number of terms
  n_terms <- nrow(top_terms)
  plot_height <- max(8, n_terms * 0.25)

  pdf(dotplot_file, width = 12, height = plot_height)

  p <- ggplot(top_terms, aes(x = NES, y = reorder(Description, NES))) +
    geom_point(aes(size = setSize, color = -log10(qvalue))) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(FDR)") +
    scale_size_continuous(name = "Gene Set Size", range = c(3, 10)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
    theme_minimal() +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = "",
      title = "GSEA: Gene Importance Enrichment by Ontology",
      subtitle = "NES > 0: enriched in high-importance genes | NES < 0: enriched in low-importance genes"
    ) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )

  print(p)
  dev.off()

  cat("  Saved dotplot:", basename(dotplot_file), "\n")
}

#' Create GSEA enrichment plots for top terms
#'
#' @param gsea_results List of GSEA results
#' @param gene_list Named vector of gene scores
#' @param output_dir Output directory
create_gsea_enrichment_plots <- function(gsea_results, gene_list, output_dir) {

  cat("Creating GSEA enrichment plots for top terms...\n")

  for (ont_code in c("bp", "mf", "cc")) {
    ont_name <- toupper(ont_code)

    if (is.null(gsea_results[[ont_code]])) {
      next
    }

    result_obj <- gsea_results[[ont_code]]
    result_df <- as.data.frame(result_obj@result)

    if (nrow(result_df) == 0) {
      next
    }

    # Get top 5 terms by p-value
    top_terms <- result_df %>%
      arrange(pvalue) %>%
      head(5)

    if (nrow(top_terms) == 0) {
      next
    }

    # Create enrichment plot file
    plot_file <- file.path(output_dir, paste0("gsea_enrichment_", ont_name, ".pdf"))

    pdf(plot_file, width = 12, height = 4 * nrow(top_terms))

    for (i in 1:nrow(top_terms)) {
      term_id <- top_terms$ID[i]

      tryCatch({
        p <- gseaplot2(result_obj,
                      geneSetID = term_id,
                      title = paste0(top_terms$Description[i],
                                   "\nNES=", round(top_terms$NES[i], 2),
                                   ", FDR=", format(top_terms$qvalue[i], digits=3)),
                      pvalue_table = TRUE,
                      ES_geom = "line")
        print(p)
      }, error = function(e) {
        cat("    Warning: Could not plot term", term_id, "\n")
      })
    }

    dev.off()

    cat("  Saved", ont_name, "enrichment plots:", basename(plot_file), "\n")
  }

  cat("\n")
}

#' Create ridge plot showing importance score distribution
#'
#' @param gsea_results List of GSEA results
#' @param gene_list Named vector of gene scores
#' @param output_dir Output directory
create_ridge_plot <- function(gsea_results, gene_list, output_dir) {

  cat("Creating ridge plot...\n")

  # Combine all results
  all_results <- list()

  for (ont_code in c("bp", "mf", "cc")) {
    if (!is.null(gsea_results[[ont_code]])) {
      result_df <- as.data.frame(gsea_results[[ont_code]]@result)
      result_df$Ontology <- toupper(ont_code)
      all_results[[ont_code]] <- result_df
    }
  }

  if (length(all_results) == 0) {
    cat("  No results for ridge plot\n")
    return()
  }

  combined_df <- bind_rows(all_results)

  # Select top 15 terms overall by p-value
  top_terms <- combined_df %>%
    arrange(pvalue) %>%
    head(15)

  if (nrow(top_terms) == 0) {
    cat("  No terms for ridge plot\n")
    return()
  }

  ridge_file <- file.path(output_dir, "gsea_ridgeplot.pdf")

  pdf(ridge_file, width = 14, height = max(8, nrow(top_terms) * 0.5))

  tryCatch({
    # Create ridge plot manually using ggplot2
    # Extract gene scores for each enriched term

    plot_data <- list()

    for (i in 1:nrow(top_terms)) {
      term_id <- top_terms$ID[i]
      term_desc <- top_terms$Description[i]
      term_ont <- top_terms$Ontology[i]

      # Get core enrichment genes
      core_genes <- strsplit(top_terms$core_enrichment[i], "/")[[1]]

      # Get scores for these genes
      if (length(core_genes) > 0 && any(core_genes %in% names(gene_list))) {
        gene_scores <- gene_list[core_genes[core_genes %in% names(gene_list)]]

        plot_data[[i]] <- data.frame(
          Term = paste0(term_desc, " [", term_ont, "]"),
          Score = gene_scores,
          Term_short = ifelse(nchar(term_desc) > 40,
                             paste0(substr(term_desc, 1, 37), "..."),
                             term_desc),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(plot_data) > 0) {
      plot_df <- bind_rows(plot_data)

      # Create ridge plot using ggridges
      if (requireNamespace("ggridges", quietly = TRUE)) {
        library(ggridges)

        p <- ggplot(plot_df, aes(x = Score, y = reorder(Term_short, Score, median),
                                fill = after_stat(x))) +
          geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
          scale_fill_gradient(low = "blue", high = "red",
                            name = "Importance\nScore") +
          theme_minimal() +
          labs(
            x = "Gene Importance Score",
            y = "",
            title = "Distribution of Gene Importance Scores in Enriched GO Terms",
            subtitle = paste0("Top ", nrow(top_terms), " enriched terms (by p-value)")
          ) +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 10),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 10)
          )

        print(p)
        cat("  Created ridge plot with ggridges\n")

      } else {
        # Fallback: use violin plot if ggridges not available
        p <- ggplot(plot_df, aes(x = Score, y = reorder(Term_short, Score, median))) +
          geom_violin(aes(fill = after_stat(x)), scale = "width") +
          geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.3) +
          scale_fill_gradient(low = "blue", high = "red",
                            name = "Importance\nScore") +
          theme_minimal() +
          labs(
            x = "Gene Importance Score",
            y = "",
            title = "Distribution of Gene Importance Scores in Enriched GO Terms",
            subtitle = paste0("Top ", nrow(top_terms), " enriched terms (by p-value)")
          ) +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 10),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 10)
          )

        print(p)
        cat("  Created violin plot (ggridges not available)\n")
      }
    } else {
      cat("  Warning: No gene scores found for plotting\n")
    }

  }, error = function(e) {
    cat("  Warning: Could not create ridge plot:", e$message, "\n")
  })

  dev.off()

  cat("  Saved ridge plot:", basename(ridge_file), "\n\n")
}

#' Print summary statistics
#'
#' @param gsea_results List of GSEA results
print_summary <- function(gsea_results) {

  cat(rep("=", 70), "\n", sep = "")
  cat("GSEA SUMMARY STATISTICS\n")
  cat(rep("=", 70), "\n\n", sep = "")

  for (ont_code in c("bp", "mf", "cc")) {
    ont_name <- toupper(ont_code)

    if (is.null(gsea_results[[ont_code]])) {
      cat(ont_name, ": No significant enrichment\n")
      next
    }

    result_df <- as.data.frame(gsea_results[[ont_code]]@result)

    if (nrow(result_df) == 0) {
      cat(ont_name, ": No significant enrichment\n")
      next
    }

    n_pos <- sum(result_df$NES > 0)
    n_neg <- sum(result_df$NES < 0)

    cat(ont_name, ":\n")
    cat("  Total enriched terms:", nrow(result_df), "\n")
    cat("  Positive NES (high-importance enriched):", n_pos, "\n")
    cat("  Negative NES (low-importance enriched):", n_neg, "\n")

    # Top 3 positive
    if (n_pos > 0) {
      cat("  Top 3 high-importance enriched:\n")
      top_pos <- result_df %>%
        filter(NES > 0) %>%
        arrange(pvalue) %>%
        head(3)

      for (i in 1:nrow(top_pos)) {
        cat("    ", i, ". ", top_pos$Description[i],
            " (NES=", round(top_pos$NES[i], 2), ", FDR=",
            format(top_pos$qvalue[i], digits=3), ")\n", sep = "")
      }
    }

    # Top 3 negative
    if (n_neg > 0) {
      cat("  Top 3 low-importance enriched:\n")
      top_neg <- result_df %>%
        filter(NES < 0) %>%
        arrange(pvalue) %>%
        head(3)

      for (i in 1:nrow(top_neg)) {
        cat("    ", i, ". ", top_neg$Description[i],
            " (NES=", round(top_neg$NES[i], 2), ", FDR=",
            format(top_neg$qvalue[i], digits=3), ")\n", sep = "")
      }
    }

    cat("\n")
  }
}

# MAIN EXECUTION
cat("Starting GSEA analysis...\n\n")

# 1. Load gene importance scores
gene_scores <- read_gene_importance(IMPORTANCE_FILE)

# 2. Load GMT annotations
gmt_data <- read_gmt_for_gsea(GMT_FILE)

# 3. Run GSEA for all ontologies
gsea_results <- run_gsea_all_ontologies(gene_scores, gmt_data)

# Check if any results
has_results <- any(sapply(gsea_results, function(x) !is.null(x)))

if (!has_results) {
  cat("No significant enrichment found at FDR <", Q_VALUE_CUTOFF, "\n")
  cat("Analysis complete.\n")
  quit(save = "no")
}

# 4. Save results
save_gsea_results(gsea_results, OUTPUT_DIR)

# 5. Create visualizations
cat("Creating visualizations...\n")
create_gsea_dotplots(gsea_results, OUTPUT_DIR)
create_gsea_enrichment_plots(gsea_results, gene_scores, OUTPUT_DIR)
create_ridge_plot(gsea_results, gene_scores, OUTPUT_DIR)

# 6. Print summary
print_summary(gsea_results)

# Session info
cat(rep("=", 70), "\n", sep = "")
cat("Session Information:\n")
cat("R version:", R.version.string, "\n")
cat("clusterProfiler version:", as.character(packageVersion("clusterProfiler")), "\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nAll outputs saved to:", OUTPUT_DIR, "\n")
