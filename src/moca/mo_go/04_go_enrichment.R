#!/usr/bin/env Rscript
# Simple GO enrichment analysis for gene lists using clusterProfiler
# Supports both single file and batch processing
#
# Usage from R terminal:
#   Single file:
#     Rscript 04_go_enrichment.R /path/to/your/gene_list.csv [gmt_file] [target_file] [pval] [qval] [minGS] [maxGS]
#
#   Batch processing (wildcard pattern):
#     Rscript 04_go_enrichment.R "../../../out/moca_results/mo_go/motif_gene_mapping/p0m_motifs/gene_lists/*.csv"
#
# Example:
#   Rscript 04_go_enrichment.R ../../../out/moca_results/mo_go/clade_shared_genes/p0m_l1_m2_yellow/p0m_l1_m2_yellow_shared_genes.csv
#
# The script will create an output directory for each file with:
#   - GO enrichment results: enrichment_results.csv, enrichment_barplot.pdf
#
# Note: Background is filtered to only include genes with target==0 or target==1 (excludes target==2)
#       to match the gene set used in motif-gene analysis
######################

library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(GO.db)
library(DOSE)  # Required for enrichplot to work properly

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for required input
if (length(args) < 1 || !nzchar(args[1])) {
  stop("\nERROR: Gene list file path is required!\n\n",
       "Usage: Rscript 04_go_enrichment.R /path/to/gene_list.csv\n",
       "   or: Rscript 04_go_enrichment.R \"/path/to/files/*.csv\"\n\n",
       "Example:\n",
       "  Rscript 04_go_enrichment.R ../../../out/moca_results/mo_go/clade_shared_genes/p0m_l1_m2_yellow/p0m_l1_m2_yellow_shared_genes.csv\n")
}

INPUT_PATH <- args[1]

# Optional parameters (can be overridden with additional arguments)
GMT_FILE <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  "../../../vitis_data/gene_ontology/blast2go_t2t_5.1.gmt"
}

TARGET_FILE <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  "../../../vitis_data/tpm_counts/vitis_drought_leaf_targets.csv"
}

P_VALUE_CUTOFF <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
Q_VALUE_CUTOFF <- if (length(args) >= 5) as.numeric(args[5]) else 0.2
MIN_GENE_SIZE <- if (length(args) >= 6) as.integer(args[6]) else 3
MAX_GENE_SIZE <- if (length(args) >= 7) as.integer(args[7]) else 500

cat("=" , rep("=", 70), "\n", sep = "")
cat("GO ENRICHMENT ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep = "")

#' Expand glob pattern to list of files
#'
#' This uses Sys.glob() which is R's built-in way to handle wildcards
#' like *.csv, similar to shell globbing
#'
#' @param path String that may contain wildcards (*, ?)
#' @return Character vector of matching file paths
expand_file_list <- function(path) {
  # Check if path contains wildcard characters
  if (grepl("[*?]", path)) {
    # Expand glob pattern
    files <- Sys.glob(path)

    if (length(files) == 0) {
      stop("No files matched the pattern: ", path)
    }

    # Filter to only keep .csv files that contain "_genes" (allows _genes.csv, _genes_high.csv, etc.)
    files <- files[grepl("_genes.*\\.csv$", files)]

    if (length(files) == 0) {
      stop("No gene list CSV files matched the pattern: ", path)
    }

    return(files)
  } else {
    # Single file - check it exists
    if (!file.exists(path)) {
      stop("File not found: ", path)
    }
    return(path)
  }
}

# Expand input path (handles both single files and wildcards)
FILE_LIST <- expand_file_list(INPUT_PATH)

# Determine if batch mode
BATCH_MODE <- length(FILE_LIST) > 1

cat("Input mode:", if (BATCH_MODE) "BATCH" else "SINGLE FILE", "\n")
cat("Files to process:", length(FILE_LIST), "\n")
if (BATCH_MODE) {
  cat("  (Showing first 5):\n")
  for (f in head(FILE_LIST, 5)) {
    cat("    -", basename(f), "\n")
  }
  if (length(FILE_LIST) > 5) {
    cat("    ... and", length(FILE_LIST) - 5, "more\n")
  }
} else {
  cat("  File:", FILE_LIST, "\n")
}
cat("\n")

cat("Parameters:\n")
cat("  GMT file:", GMT_FILE, "\n")
cat("  Target file:", TARGET_FILE, "\n")
cat("  P-value cutoff:", P_VALUE_CUTOFF, "\n")
cat("  Q-value (FDR) cutoff:", Q_VALUE_CUTOFF, "\n")
cat("  Min gene set size:", MIN_GENE_SIZE, "\n")
cat("  Max gene set size:", MAX_GENE_SIZE, "\n\n")

# Check GMT file exists
if (!file.exists(GMT_FILE)) {
  stop("GMT file not found: ", GMT_FILE)
}

# Check target file exists
if (!file.exists(TARGET_FILE)) {
  stop("Target file not found: ", TARGET_FILE)
}

#' Read gene list from CSV (handles comment lines)
#'
#' @param file_path Path to gene list CSV
#' @return Character vector of gene IDs
read_gene_list <- function(file_path) {
  # Read all lines
  all_lines <- readLines(file_path)

  # Find first non-comment line (this should be the header)
  non_comment_lines <- which(!grepl("^#", all_lines) & !grepl("^\"#", all_lines))

  if (length(non_comment_lines) == 0) {
    stop("No data found in file: ", file_path)
  }

  # First non-comment line is the header
  header_line <- non_comment_lines[1]

  # Number of lines to skip before header
  skip_lines <- header_line - 1

  # Read CSV starting from header
  gene_data <- read.csv(file_path, skip = skip_lines, stringsAsFactors = FALSE)

  # Check if gene_id column exists
  if (!"gene_id" %in% colnames(gene_data)) {
    stop("No 'gene_id' column found in file: ", file_path)
  }

  # Remove any remaining comment rows and duplicate headers
  gene_data <- gene_data %>%
    filter(!grepl("^#", gene_id)) %>%
    filter(gene_id != "gene_id")

  return(gene_data$gene_id)
}

#' Read GMT file and create TERM2GENE and TERM2NAME data frames, split by ontology
#'
#' Filters background to only include genes with target==0 or target==1 (excludes target==2)
#'
#' @param gmt_file Path to GMT file
#' @param target_file Path to target classification file
#' @return List with term2gene, term2name, and data split by ontology (BP, MF, CC)
read_gmt_for_enricher <- function(gmt_file, target_file) {
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

  # Load target classifications
  cat("  Loading target classifications to filter background...\n")
  targets <- read.csv(target_file, stringsAsFactors = FALSE)

  # Get genes with target==0 or target==1 (exclude target==2)
  valid_targets <- targets %>% filter(target %in% c(0, 1))
  valid_gene_set <- valid_targets$gene_id

  cat("  Target filtering:\n")
  cat("    Total genes in target file:", nrow(targets), "\n")
  cat("    Genes with target==0 or target==1:", length(valid_gene_set), "\n")
  cat("    Genes excluded (target==2):", sum(targets$target == 2), "\n")

  # Filter term2gene to only include valid genes (target==0 or target==1)
  n_before <- nrow(term2gene)
  term2gene <- term2gene %>% filter(gene %in% valid_gene_set)
  n_after <- nrow(term2gene)

  cat("  Background filtering applied:\n")
  cat("    Gene-term associations before filtering:", n_before, "\n")
  cat("    Gene-term associations after filtering:", n_after, "\n")
  cat("    Associations removed (target==2 genes):", n_before - n_after, "\n")

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

  # Get unique genes (now filtered to target==0 or target==1)
  unique_background <- unique(term2gene$gene)

  cat("  Loaded", length(unique(term2gene$term)), "GO terms\n")
  cat("    BP:", nrow(term2name_bp), "terms\n")
  cat("    MF:", nrow(term2name_mf), "terms\n")
  cat("    CC:", nrow(term2name_cc), "terms\n")
  cat("  Total gene-term associations:", nrow(term2gene), "\n")
  cat("  Genes with GO annotation (target==0 or 1 only):", length(unique_background), "\n")
  cat("  --> This is the filtered background for enrichment testing\n")
  cat("\n")

  return(list(
    term2gene = term2gene,
    term2name = term2name,
    bp = list(term2gene = term2gene_bp, term2name = term2name_bp),
    mf = list(term2gene = term2gene_mf, term2name = term2name_mf),
    cc = list(term2gene = term2gene_cc, term2name = term2name_cc)
  ))
}

#' Perform GO enrichment analysis separately for BP, MF, and CC
#'
#' @param gene_list Character vector of gene IDs
#' @param gmt_data List with term2gene and term2name split by ontology
#' @return List with enrichment results for BP, MF, CC
perform_enrichment_by_ontology <- function(gene_list, gmt_data) {
  results <- list()

  # Enrichment for Biological Process
  cat("  Running enrichment for Biological Process (BP)...\n")
  if (nrow(gmt_data$bp$term2gene) > 0) {
    results$bp <- tryCatch({
      enricher(
        gene = gene_list,
        TERM2GENE = gmt_data$bp$term2gene,
        TERM2NAME = gmt_data$bp$term2name[, c("term", "name")],
        pvalueCutoff = P_VALUE_CUTOFF,
        qvalueCutoff = Q_VALUE_CUTOFF,
        minGSSize = MIN_GENE_SIZE,
        maxGSSize = MAX_GENE_SIZE,
        pAdjustMethod = "BH"
      )
    }, error = function(e) NULL)

    if (!is.null(results$bp) && nrow(results$bp@result) > 0) {
      cat("    Found", nrow(results$bp@result), "enriched BP terms\n")
    } else {
      cat("    No significant BP enrichment\n")
      results$bp <- NULL
    }
  }

  # Enrichment for Molecular Function
  cat("  Running enrichment for Molecular Function (MF)...\n")
  if (nrow(gmt_data$mf$term2gene) > 0) {
    results$mf <- tryCatch({
      enricher(
        gene = gene_list,
        TERM2GENE = gmt_data$mf$term2gene,
        TERM2NAME = gmt_data$mf$term2name[, c("term", "name")],
        pvalueCutoff = P_VALUE_CUTOFF,
        qvalueCutoff = Q_VALUE_CUTOFF,
        minGSSize = MIN_GENE_SIZE,
        maxGSSize = MAX_GENE_SIZE,
        pAdjustMethod = "BH"
      )
    }, error = function(e) NULL)

    if (!is.null(results$mf) && nrow(results$mf@result) > 0) {
      cat("    Found", nrow(results$mf@result), "enriched MF terms\n")
    } else {
      cat("    No significant MF enrichment\n")
      results$mf <- NULL
    }
  }

  # Enrichment for Cellular Component
  cat("  Running enrichment for Cellular Component (CC)...\n")
  if (nrow(gmt_data$cc$term2gene) > 0) {
    results$cc <- tryCatch({
      enricher(
        gene = gene_list,
        TERM2GENE = gmt_data$cc$term2gene,
        TERM2NAME = gmt_data$cc$term2name[, c("term", "name")],
        pvalueCutoff = P_VALUE_CUTOFF,
        qvalueCutoff = Q_VALUE_CUTOFF,
        minGSSize = MIN_GENE_SIZE,
        maxGSSize = MAX_GENE_SIZE,
        pAdjustMethod = "BH"
      )
    }, error = function(e) NULL)

    if (!is.null(results$cc) && nrow(results$cc@result) > 0) {
      cat("    Found", nrow(results$cc@result), "enriched CC terms\n")
    } else {
      cat("    No significant CC enrichment\n")
      results$cc <- NULL
    }
  }

  # Check if we have any results
  if (is.null(results$bp) && is.null(results$mf) && is.null(results$cc)) {
    return(NULL)
  }

  return(results)
}


#' Create output directory for a gene list file
#'
#' @param input_file Path to input gene list file
#' @param batch_mode Logical, whether processing multiple files
#' @return Path to output directory
create_output_dir <- function(input_file, batch_mode = FALSE) {
  # Get directory and base name
  input_dir <- dirname(input_file)
  base_name <- tools::file_path_sans_ext(basename(input_file))

  # Remove "_genes" suffix if present (for cleaner output names)
  base_name <- gsub("_genes$", "", base_name)

  # Create output directory name
  # Batch mode: shorter naming (e.g., "epm_vitis_ssr_p0m00_go")
  # Single mode: descriptive naming (e.g., "p0m_l1_m2_yellow_enrichment")
  if (batch_mode) {
    output_dir <- file.path(input_dir, paste0(base_name, "_go"))
  } else {
    output_dir <- file.path(input_dir, paste0(base_name, "_enrichment"))
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  return(output_dir)
}

#' Process a single gene list file
#'
#' This is the core function that does the enrichment for one file.
#' By wrapping the processing in a function, we can easily call it
#' multiple times in a loop for batch processing.
#'
#' @param file_path Path to gene list CSV
#' @param gmt_data Pre-loaded GMT data (term2gene and term2name)
#' @param batch_mode Logical, whether this is part of batch processing
#' @return List with success status and results
process_single_file <- function(file_path, gmt_data, batch_mode = FALSE) {

  cat("\n", rep("-", 70), "\n", sep = "")
  cat("Processing:", basename(file_path), "\n")
  cat(rep("-", 70), "\n", sep = "")

  # Try to process, catch errors gracefully
  tryCatch({

    # Load gene list
    cat("Loading gene list...\n")
    gene_list <- read_gene_list(file_path)
    cat("  Loaded", length(gene_list), "genes\n")

    # Check overlap with background
    overlap <- intersect(gene_list, unique(gmt_data$term2gene$gene))
    cat("  Genes with GO annotations:", length(overlap),
        "(", round(100 * length(overlap) / length(gene_list), 1), "%)\n")

    if (length(overlap) < 3) {
      cat("  WARNING: Too few genes with GO annotations (< 3). Skipping.\n")
      return(list(success = FALSE, reason = "insufficient_overlap"))
    }

    # Perform enrichment by ontology
    cat("Performing GO enrichment analysis by ontology...\n")
    enrich_results <- perform_enrichment_by_ontology(gene_list, gmt_data)

    if (is.null(enrich_results)) {
      cat("  WARNING: No significant GO enrichment found. Skipping.\n")
      return(list(success = FALSE, reason = "no_enrichment"))
    }

    # Count total enriched terms
    n_total_terms <- sum(
      if(!is.null(enrich_results$bp)) nrow(enrich_results$bp@result) else 0,
      if(!is.null(enrich_results$mf)) nrow(enrich_results$mf@result) else 0,
      if(!is.null(enrich_results$cc)) nrow(enrich_results$cc@result) else 0
    )
    cat("  Total enriched terms across all ontologies:", n_total_terms, "\n")

    # Create output directory
    output_dir <- create_output_dir(file_path, batch_mode)

    # Combine all results into one dataframe for the combined CSV
    all_results <- bind_rows(
      if(!is.null(enrich_results$bp)) mutate(as.data.frame(enrich_results$bp@result), Ontology = "BP"),
      if(!is.null(enrich_results$mf)) mutate(as.data.frame(enrich_results$mf@result), Ontology = "MF"),
      if(!is.null(enrich_results$cc)) mutate(as.data.frame(enrich_results$cc@result), Ontology = "CC")
    )

    # Add calculated columns
    all_results <- all_results %>%
      mutate(
        GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
        BgRatio_numeric = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
        FoldEnrichment = GeneRatio_numeric / BgRatio_numeric,
        Count = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]))
      ) %>%
      arrange(Ontology, qvalue, pvalue)

    # Write separate CSV files for each ontology
    cat("Writing enrichment results...\n")

    if (!is.null(enrich_results$bp)) {
      bp_file <- file.path(output_dir, "enrichment_BP.csv")
      bp_results <- all_results %>% filter(Ontology == "BP")
      write.csv(bp_results, bp_file, row.names = FALSE)
      cat("  Saved BP results:", basename(bp_file), "(", nrow(bp_results), "terms )\n")
    }

    if (!is.null(enrich_results$mf)) {
      mf_file <- file.path(output_dir, "enrichment_MF.csv")
      mf_results <- all_results %>% filter(Ontology == "MF")
      write.csv(mf_results, mf_file, row.names = FALSE)
      cat("  Saved MF results:", basename(mf_file), "(", nrow(mf_results), "terms )\n")
    }

    if (!is.null(enrich_results$cc)) {
      cc_file <- file.path(output_dir, "enrichment_CC.csv")
      cc_results <- all_results %>% filter(Ontology == "CC")
      write.csv(cc_results, cc_file, row.names = FALSE)
      cat("  Saved CC results:", basename(cc_file), "(", nrow(cc_results), "terms )\n")
    }

    # Write combined results
    all_file <- file.path(output_dir, "enrichment_all.csv")
    write.csv(all_results, all_file, row.names = FALSE)
    cat("  Saved combined results:", basename(all_file), "(", nrow(all_results), "terms )\n")

    # Create visualizations
    cat("Creating visualizations...\n")

    # Define ontology colors
    ontology_colors <- c(
      "BP" = "#671436",  # Burgundy for Biological Process
      "MF" = "#052b67",  # Dark blue for Molecular Function
      "CC" = "#176272"   # Teal for Cellular Component
    )

    # 1. Create combined barplot
    cat("  Creating combined barplot...\n")
    barplot_file <- file.path(output_dir, "enrichment_barplot.pdf")

    # Select top terms per ontology for barplot
    top_per_ontology <- all_results %>%
      group_by(Ontology) %>%
      slice_min(order_by = qvalue, n = 10, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(Ontology, qvalue)

    n_terms_plot <- nrow(top_per_ontology)

    if (n_terms_plot > 0) {
      pdf(barplot_file, width = 12, height = max(8, n_terms_plot * 0.4))
      par(mar = c(5, 28, 4, 2))

      # Reverse order for plotting (most significant at top)
      plot_data <- rev(-log10(top_per_ontology$qvalue))
      plot_ontology <- rev(top_per_ontology$Ontology)
      plot_names <- rev(paste0("[", top_per_ontology$Ontology, "] ", top_per_ontology$Description))

      # Truncate long names
      plot_names <- sapply(plot_names, function(x) {
        if (nchar(x) > 65) paste0(substr(x, 1, 62), "...") else x
      })

      # Assign colors
      bar_colors <- sapply(plot_ontology, function(ont) ontology_colors[[ont]])

      # Create barplot
      barplot(
        plot_data,
        names.arg = plot_names,
        horiz = TRUE,
        las = 1,
        col = bar_colors,
        border = "black",
        xlab = "-log10(FDR)",
        main = paste0("GO Enrichment: ", tools::file_path_sans_ext(basename(file_path)),
                     "\n(Top ", n_terms_plot, " terms by ontology)"),
        cex.names = 0.65,
        cex.axis = 0.9,
        xlim = c(0, max(plot_data) * 1.1),
        cex.main = 0.9
      )

      # Add legend
      if (length(unique(plot_ontology)) > 1) {
        unique_onts <- unique(top_per_ontology$Ontology)
        legend("bottomright",
               legend = c("BP = Biological Process",
                         "MF = Molecular Function",
                         "CC = Cellular Component")[c("BP", "MF", "CC") %in% unique_onts],
               fill = ontology_colors[unique_onts],
               border = "black",
               cex = 0.7,
               bg = "white")
      }

      # Add FDR threshold line
      abline(v = -log10(Q_VALUE_CUTOFF), col = "red", lty = 2, lwd = 2)
      text(x = -log10(Q_VALUE_CUTOFF), y = par("usr")[4] * 0.95,
           labels = paste0("FDR=", Q_VALUE_CUTOFF), pos = 4, col = "red", cex = 0.8)

      dev.off()
      cat("    Saved:", basename(barplot_file), "\n")
    }

    # 2. Create faceted dotplot using ggplot2
    cat("  Creating faceted dotplot...\n")
    dotplot_file <- file.path(output_dir, "enrichment_dotplot_by_ontology.pdf")

    # Select top terms per ontology for dotplot
    top_for_dotplot <- all_results %>%
      group_by(Ontology) %>%
      slice_min(order_by = qvalue, n = 15, with_ties = FALSE) %>%
      ungroup()

    if (nrow(top_for_dotplot) >= 3) {
      # Wrap long descriptions (> 5 words) to multiple lines
      top_for_dotplot$Description_wrapped <- sapply(top_for_dotplot$Description, function(desc) {
        words <- strsplit(desc, " ")[[1]]
        if (length(words) > 5) {
          # Split approximately in half
          split_point <- ceiling(length(words) / 2)
          line1 <- paste(words[1:split_point], collapse = " ")
          line2 <- paste(words[(split_point + 1):length(words)], collapse = " ")
          return(paste0(line1, "\n", line2))
        } else {
          return(desc)
        }
      })

      # Calculate plot height based on number of terms per ontology
      max_terms_per_ont <- top_for_dotplot %>%
        count(Ontology) %>%
        pull(n) %>%
        max()

      calculated_height <- max(8, max_terms_per_ont * 0.35 * length(unique(top_for_dotplot$Ontology)))

      # Extract motif name from file path (remove _genes suffix if present)
      motif_name <- tools::file_path_sans_ext(basename(file_path))
      motif_name <- gsub("_genes$", "", motif_name)

      # Create ggplot2 dotplot
      p <- ggplot(top_for_dotplot, aes(x = GeneRatio_numeric, y = reorder(Description_wrapped, GeneRatio_numeric))) +
        geom_point(aes(size = Count, color = -log10(qvalue))) +
        scale_color_gradient(low = "blue", high = "red", name = "-log10(FDR)") +
        scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
        facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
        theme_minimal() +
        labs(
          x = "Gene Ratio",
          y = "",
          title = paste0("GO Enrichment: ", motif_name),
          subtitle = "BP = Biological Process, MF = Molecular Function, CC = Cellular Component"
        ) +
        theme(
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10)
        )

      ggsave(dotplot_file, plot = p, width = 12, height = calculated_height, limitsize = FALSE)
      cat("    Saved:", basename(dotplot_file), "\n")
    } else {
      cat("    Skipping dotplot (insufficient terms)\n")
    }

    return(list(
      success = TRUE,
      n_terms = n_total_terms,
      output_dir = output_dir
    ))

  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(list(success = FALSE, reason = "error", error = e$message))
  })
}
# Main execution
cat("Starting enrichment analysis...\n\n")

# Load GMT annotations ONCE (efficient for batch processing)
gmt_data <- read_gmt_for_enricher(GMT_FILE, TARGET_FILE)

# Process all files
# We use a for loop to iterate through FILE_LIST
# For single files, FILE_LIST has length 1, so loop runs once
# For batch mode, FILE_LIST has multiple files, loop runs for each

results_summary <- data.frame(
  file = character(),
  status = character(),
  n_terms = integer(),
  output_dir = character(),
  stringsAsFactors = FALSE
)

for (file_path in FILE_LIST) {
  result <- process_single_file(file_path, gmt_data, batch_mode = BATCH_MODE)

  # Track results
  results_summary <- rbind(results_summary, data.frame(
    file = basename(file_path),
    status = if (result$success) "SUCCESS" else result$reason,
    n_terms = if (result$success) result$n_terms else 0,
    output_dir = if (result$success) result$output_dir else "",
    stringsAsFactors = FALSE
  ))
}

# Final summary
cat("\n", rep("=", 70), "\n", sep = "")
cat("BATCH PROCESSING COMPLETE\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Summary:\n")
cat("  Total files processed:", nrow(results_summary), "\n")
cat("  Successful:", sum(results_summary$status == "SUCCESS"), "\n")
cat("  Failed (no enrichment):", sum(results_summary$status == "no_enrichment"), "\n")
cat("  Failed (insufficient data):", sum(results_summary$status == "insufficient_overlap"), "\n")
cat("  Failed (errors):", sum(results_summary$status == "error"), "\n\n")

if (BATCH_MODE) {
  # Print detailed summary for batch mode
  cat("Detailed Results:\n")
  print(results_summary, row.names = FALSE)
}

cat("\nEnrichment Totals:\n")
cat("  Total GO terms found:", sum(results_summary$n_terms), "\n")

cat("\nSession Information:\n")
cat("R version:", R.version.string, "\n")
cat("clusterProfiler version:", as.character(packageVersion("clusterProfiler")), "\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
