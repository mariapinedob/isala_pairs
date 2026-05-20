#!/usr/bin/env Rscript 

# This script will calculate a denoised popANI-based SPM value given 
# inStrain mapping profiles of WGS/MGS samples to a reference assemblies.
#
# Run as: ./popani_denoiser.R comparisons.tsv output.tsv

suppressMessages(library(tidyverse))

# Function to compute denoised SPM 
denoised_spm <- function(profile_output, N = 5) {
  
  if (!dir.exists(profile_output)) stop("profile output not found")
  
  snvs <- 
    paste0(profile_output, "*_SNVs.tsv") %>%
    Sys.glob() %>%
    read_tsv(col_names = TRUE, col_types = cols())
  
  scaffolds <-
    paste0(profile_output, "*_scaffold_info.tsv") %>%
    Sys.glob() %>%
    read_tsv(col_names = TRUE, col_types = cols())
  
  genomes <-
    paste0(profile_output, "*_genome_info.tsv") %>%
    Sys.glob() %>%
    read_tsv(col_names = TRUE, col_types = cols())
  
  n_pop_snvs <- 
    snvs %>%
    filter(class == "pop_SNV" | class == "SNS") %>%
    left_join(select(scaffolds, scaffold, length), by = "scaffold") %>%
    filter(position >= {{N}}, position < (length - N)) %>%
    nrow()
  
  n_pop_snvs / genomes$length / genomes$breadth_minCov * 1e6
}

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Base directories (aligned with step 05)
dir_project <- "../../"
dir_profiles <- 
  file.path(
  dir_project,
  "results/strain_level_analyses/isala/ref_independent/profiles"
)

dout <- 
  file.path(
  dir_project,
  "results/strain_level_analyses/isala/ref_independent/processed"
)

dir.create(dout, recursive = TRUE, showWarnings = FALSE)

# Decide mode: file input vs automatic scan
if (length(args) >= 2) {
  
  # original behavior (use provided comparisons file)
  fin <- args[1]
  fout <- args[2]
  
  colnames <- c("comparison", "profile_output")
  comparisons <- read_tsv(fin, col_names = colnames, col_types = cols())
  
} else {
  
  # Automatic mode (scan profile directories)
  profiles <- list.dirs(dir_profiles, full.names = TRUE, recursive = FALSE)
  
  comparisons <- tibble(
    comparison = basename(profiles),
    profile_output = profiles
  )
  
  fout <- file.path(dout, "denoised_popspm_values.tsv")
}

# Compute denoised popSPM
comparisons$denoised_popspm <- 
  sapply(comparisons$profile_output, denoised_spm)

# Save results (only comparison + value)
comparisons %>% 
  select(comparison, denoised_popspm) %>%
  write_tsv(fout, col_names = FALSE)

message("finished computing denoised popSPM values")
