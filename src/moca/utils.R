# R Utilities for Vitis CRE moca Project
# This file contains reusable functions for project setup and common operations
# Source this file at the beginning of analysis scripts with: source("../utils.R")

#' Find and set project root directory
#' 
#' This function searches upward from the current location to find the project root,
#' identified by the presence of key project markers. It then sets the working
#' directory to the 'src' subdirectory for consistency across all scripts.
#' 
#' @param start_path Character. Starting path for the search (default: current working directory)
#' @param project_name Character. Name pattern to look for in project directory (default: "vitis_cre")
#' @param verbose Logical. Whether to print diagnostic messages (default: TRUE)
#' @return Character. Path to the project root directory
setup_project_environment <- function(start_path = getwd(), 
                                     project_name = "vitis_cre", 
                                     verbose = TRUE) {
  
  if (verbose) cat("Setting up project environment...\n")
  
  # Normalize the starting path to handle relative paths properly
  current_path <- normalizePath(start_path, mustWork = FALSE)
  
  # Define what we're looking for to identify the project root
  # These are common indicators that help us locate the correct directory
  project_indicators <- c(".git", "README.md", "requirements.txt", "train.log", ".gitignore")
  
  # Walk up the directory tree until we find our project or hit the filesystem root
  search_attempts <- 0
  max_attempts <- 10  # Prevent infinite loops on unusual filesystems
  
  while (current_path != dirname(current_path) && search_attempts < max_attempts) {
    search_attempts <- search_attempts + 1
    
    # Get contents of current directory (including hidden files like .git)
    dir_contents <- list.files(current_path, all.files = TRUE, include.dirs = TRUE)
    
    # Check if this looks like our project root
    has_project_name <- grepl(project_name, basename(current_path), ignore.case = TRUE)
    has_src_dir <- "src" %in% dir_contents
    has_git <- ".git" %in% dir_contents
    
    # We've found the project root if we have the project name and key indicators
    if (has_project_name && (has_src_dir || has_git)) {
      project_root <- current_path
      src_dir <- file.path(project_root, "src")
      
      # Verify that src directory exists, create if necessary
      if (!dir.exists(src_dir)) {
        if (verbose) cat("Creating src directory at:", src_dir, "\n")
        dir.create(src_dir, recursive = TRUE)
      }
      
      # Set working directory to src for consistency
      setwd(src_dir)
      if (verbose) {
        cat("✓ Project root found:", project_root, "\n")
        cat("✓ Working directory set to:", getwd(), "\n")
      }
      
      return(invisible(src_dir))
    }
    
    # Move up one directory level and continue searching
    current_path <- dirname(current_path)
  }
  
  # If we get here, we couldn't find the project root
  warning("Could not locate project root containing '", project_name, 
          "'. Using current directory: ", getwd())
  return(invisible(getwd()))
}

#' Load required packages with informative error handling
#' 
#' This function attempts to load a list of required packages, providing
#' clear error messages if packages are missing and suggestions for installation.
#' 
#' @param packages Character vector. Names of packages to load
#' @param verbose Logical. Whether to print loading messages (default: TRUE)
#' @return Logical. TRUE if all packages loaded successfully, FALSE otherwise
load_required_packages <- function(packages, verbose = TRUE) {
  
  if (verbose) cat("Loading required packages...\n")
  
  missing_packages <- c()
  loaded_packages <- c()
  
  for (pkg in packages) {
    # Check if package is installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    } else {
      # Load the package
      library(pkg, character.only = TRUE, quietly = !verbose)
      loaded_packages <- c(loaded_packages, pkg)
      if (verbose) cat("  ✓", pkg, "\n")
    }
  }
  
  # Report results
  if (length(missing_packages) > 0) {
    cat("\n⚠ Missing packages detected:\n")
    for (pkg in missing_packages) {
      cat("  ✗", pkg, "\n")
    }
    cat("\nTo install missing packages, run:\n")
    cat("install.packages(c(", paste0('"', missing_packages, '"', collapse = ", "), "))\n\n")
    
    # For Bioconductor packages, provide alternative installation instructions
    bioc_packages <- c("rhdf5", "GenomicRanges", "Biostrings", "rtracklayer")
    missing_bioc <- intersect(missing_packages, bioc_packages)
    if (length(missing_bioc) > 0) {
      cat("For Bioconductor packages, use:\n")
      cat("if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n")
      cat("BiocManager::install(c(", paste0('"', missing_bioc, '"', collapse = ", "), "))\n\n")
    }
    
    return(FALSE)
  }
  
  if (verbose) cat("All packages loaded successfully!\n\n")
  return(TRUE)
}

#' Find files matching a pattern in a directory
#' 
#' This function provides a more user-friendly interface for finding files
#' with pattern matching, including helpful error messages when files aren't found.
#' 
#' @param directory Character. Directory to search in
#' @param pattern Character. File pattern to match (regex)
#' @param description Character. Human-readable description of what we're looking for
#' @param select_multiple Logical. If TRUE and multiple files found, return all. If FALSE, return first and warn.
#' @return Character vector. Paths to matching files
find_project_files <- function(directory, pattern, description, select_multiple = FALSE) {
  
  # Check if directory exists
  if (!dir.exists(directory)) {
    stop("Directory not found: ", directory, 
         "\nMake sure you're running from the correct project location.")
  }
  
  # Find matching files
  matching_files <- list.files(directory, pattern = pattern, full.names = TRUE)
  
  # Handle different scenarios
  if (length(matching_files) == 0) {
    stop("No ", description, " found in ", directory, 
         "\nPattern searched: ", pattern)
  } else if (length(matching_files) == 1) {
    cat("Found", description, ":", matching_files[1], "\n")
    return(matching_files[1])
  } else {
    # Multiple files found
    if (select_multiple) {
      cat("Found", length(matching_files), description, "files:\n")
      for (file in matching_files) {
        cat("  -", file, "\n")
      }
      return(matching_files)
    } else {
      cat("Multiple", description, "files found. Using:", matching_files[1], "\n")
      cat("Other options:\n")
      for (file in matching_files[-1]) {
        cat("  -", file, "\n")
      }
      return(matching_files[1])
    }
  }
}

#' Create output directories with proper structure
#' 
#' Creates a standard directory structure for analysis outputs,
#' following consistent naming conventions across the project.
#' 
#' @param base_dir Character. Base directory for outputs (default: "../../out")
#' @param analysis_name Character. Name of the current analysis
#' @param subdirs Character vector. Subdirectories to create (default: c("figures", "tables", "data"))
#' @return Named character vector. Paths to created directories
setup_output_directories <- function(base_dir = "../../out", 
                                   analysis_name, 
                                   subdirs = c("figures", "tables", "data")) {
  
  if (missing(analysis_name)) {
    stop("analysis_name is required to create organized output structure")
  }
  
  # Create main analysis directory
  main_dir <- file.path(base_dir, analysis_name)
  dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create subdirectories
  created_dirs <- c(main = main_dir)
  for (subdir in subdirs) {
    full_path <- file.path(main_dir, subdir)
    dir.create(full_path, showWarnings = FALSE, recursive = TRUE)
    created_dirs[subdir] <- full_path
  }
  
  cat("Created output directories for analysis:", analysis_name, "\n")
  for (i in seq_along(created_dirs)) {
    cat("  ", names(created_dirs)[i], ":", created_dirs[i], "\n")
  }
  cat("\n")
  
  return(created_dirs)
}

#' Print session information for reproducibility
#' 
#' Prints comprehensive session information including R version,
#' loaded packages, and system details for research reproducibility.
print_session_info <- function() {
  cat("\n" , rep("=", 50), "\n", sep = "")
  cat("SESSION INFORMATION FOR REPRODUCIBILITY\n")
  cat(rep("=", 50), "\n\n", sep = "")
  
  cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
  cat("Working directory:", getwd(), "\n\n")
  cat("R Session Details:\n")
  print(sessionInfo())
}