# Extract the contribution scores of EPMs
# take as input hdf5 files from modisco
# adapted from moca_blue by Dr. SM Zumkeller
# adaptations by L. Rothmann
# DeepCRE/moca_blue/mo_imp
######################

library(dplyr)
library(purrr)
library(rhdf5)

args <- commandArgs(trailingOnly = TRUE)

# Defaults
default_spec <- "vitis"
default_model <- "ssr"
default_input_dir <- "./results/modisco"
default_output_dir <- "./results/modisco/stats"

# Assign positional arguments or use defaults
SPEC <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_spec
MODEL <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_model
DATE <- if (length(args) >= 3 && nzchar(args[3])) args[3] else format(Sys.Date(), "%Y%m%d")

# example usage:
# Rscript get_cwms.R vitis ssr 20231001 ./results/modisco/vitis_02_modisco.hdf5 ./results/modisco/stats

# Optional input file
if (length(args) >= 4 && nzchar(args[4])) {
  FILE1 <- args[4]
} else {
  modisco_files <- list.files(default_input_dir, pattern = "*_modisco.hdf5$", full.names = TRUE)
  if (length(modisco_files) == 0) stop("No *_modisco.hdf5 files found in ", default_input_dir)
  if (length(modisco_files) > 1) stop("Multiple *_modisco.hdf5 files found. Please specify one explicitly.")
  FILE1 <- modisco_files[1]
}


dirpath_out <- if (length(args) >= 5 && nzchar(args[5])) args[5] else default_output_dir
file_path_out <- file.path(dirpath_out, paste0(DATE, "_", SPEC, MODEL, "_contrib_scores"))
h5file <- H5Fopen(FILE1, "H5F_ACC_RDONLY")

# Access the top-level group with metaclusters
metacluster_group <- h5read(h5file, "metacluster_idx_to_submetacluster_results")

# Optional sanity check: print metacluster names
cat("Found metaclusters:\n")
print(names(metacluster_group))

#######################################################
# loop through the metaclusters and generate a list of patterns
#######################################################

patterns_list <- list()
for (metacluster_name in names(metacluster_group)) {
  metacluster <- metacluster_group[[metacluster_name]]
  if (!is.null(metacluster[["seqlets_to_patterns_result"]][["patterns"]])) {
    patterns_list[[metacluster_name]] <- metacluster[["seqlets_to_patterns_result"]][["patterns"]]
  } else {
    warning(paste("No patterns found in", metacluster_name))
  }
}

#######################################################
# loop through the patterns and extract the contribution scores
#######################################################

cwm_list <- list()
metacluster_names <- names(metacluster_group)

for (metacluster_name in names(patterns_list)) {
  patterns <- patterns_list[[metacluster_name]]
  
  cwm_list[[metacluster_name]] <- list(fwd = list(), rev = list())
  
  pattern_names_raw <- patterns[["all_pattern_names"]]

  # Convert raw pattern names to character
  pattern_names <- sapply(pattern_names_raw, function(x) {
    if (is.raw(x)) rawToChar(x) else as.character(x)
  })
  
  for (pattern_name in pattern_names) {
    # Defensive check - sometimes patterns might be missing the contrib_scores group
    if (!is.null(patterns[[pattern_name]][["task0_contrib_scores"]])) {
      
      cwm_fwd <- patterns[[pattern_name]][["task0_contrib_scores"]][["fwd"]]
      cwm_rev <- patterns[[pattern_name]][["task0_contrib_scores"]][["rev"]]
      
      cwm_list[[metacluster_name]][["fwd"]][[pattern_name]] <- cwm_fwd
      cwm_list[[metacluster_name]][["rev"]][[pattern_name]] <- cwm_rev
    } else {
      warning(paste("No contrib scores found for", pattern_name, "in", metacluster_name))
    }
  }
}

matricesF0 <- cwm_list[["metacluster_0"]][["fwd"]]
matricesF1 <- cwm_list[["metacluster_1"]][["fwd"]]
matricesR0 <- cwm_list[["metacluster_0"]][["rev"]]
matricesR1 <- cwm_list[["metacluster_1"]][["rev"]]

###################################
# get seqlet counts per pattern
###################################

seqlet_counts <- list()

for (metacluster_name in metacluster_names) {
  seqlet_counts[[metacluster_name]] <- list()
  
  pattern_names <- metacluster_group[[metacluster_name]][["seqlets_to_patterns_result"]][["patterns"]][["all_pattern_names"]]

  for (pattern_name in pattern_names) {
    seqlets <- metacluster_group[[metacluster_name]][["seqlets_to_patterns_result"]][["patterns"]][[pattern_name]][["seqlets_and_alnmts"]][["seqlets"]]
    seqlet_counts[[metacluster_name]][[pattern_name]] <- length(seqlets)
  }
}

############################
# assign nomenclature
############################

rename_motifs_with_seqlets <- function(motif_list, metacluster_name, strand) {
  for (i in seq_along(motif_list)) {
    orig_name <- names(motif_list)[i] 
    pattern_name <- orig_name
    count <- seqlet_counts[[metacluster_name]][[pattern_name]]
    
    new_name <- paste0("epm_", SPEC, "_", MODEL, "_",
                       substring(pattern_name, 9),  # extract numeric part e.g. "0"
                       "_", metacluster_name, "_", strand,
                       "_", count)
    
    names(motif_list)[i] <- new_name
  }
  return(motif_list)
}

matricesF0 <- rename_motifs_with_seqlets(matricesF0, "metacluster_0", "F")
matricesF1 <- rename_motifs_with_seqlets(matricesF1, "metacluster_1", "F")
matricesR0 <- rename_motifs_with_seqlets(matricesR0, "metacluster_0", "R")
matricesR1 <- rename_motifs_with_seqlets(matricesR1, "metacluster_1", "R")


############################
# compute contribution score summaries
############################

summarise_contrib_scores <- function(matrices) {
  imap_dfr(matrices, function(mat, motif_name) {
    tibble(
      motif = motif_name,
      contrib_score_sum = sum(rowSums(mat)) + sum(colSums(mat)),
      contrib_score_max = max(mat),
      contrib_score_min = min(mat)
    )
  })
}

contrib_scores_F0 <- summarise_contrib_scores(matricesF0)
contrib_scores_F1 <- summarise_contrib_scores(matricesF1)
contrib_scores_R0 <- summarise_contrib_scores(matricesR0)
contrib_scores_R1 <- summarise_contrib_scores(matricesR1)


# Combine all results into a single data frame
contrib_score_table <- bind_rows(
  mutate(contrib_scores_F0, source = "F0"),
  mutate(contrib_scores_R0, source = "R0"),
  mutate(contrib_scores_F1, source = "F1"),
  mutate(contrib_scores_R1, source = "R1")
)

if (!dir.exists(dirpath_out)) dir.create(dirpath_out, recursive = TRUE)
write.csv(contrib_score_table, file = paste0(file_path_out, ".csv"), row.names = FALSE)

cat("Contribution scores written to:", paste0(file_path_out, ".csv"), "\n")