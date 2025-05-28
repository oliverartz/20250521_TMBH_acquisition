## -----------------------------------------------------------------------------
## Purpose of script: Make MAF files for neoantigen prediction pipeline and 
##                    signature deconvolution
##
## Author: Oliver Artz
## Date Created: May 26, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, here)

# set parameters ---------------------------------------------------------------
project_root <- here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# output dir for MAF files
output_dir <- paste0(project_dir, "/data/processed/006_maf/")

# load data --------------------------------------------------------------------
source(paste0(project_dir,"/analysis/scripts/utils/001_preprocess_utils.R"))

annotated_mutations_progression <- fread(paste0(project_dir, "/data/processed/005_annotated_mutations_progression.txt"))

impact_mod <- fread(paste0(project_dir, "/data/processed/004_impact_mod.txt"))

# wrangle ----------------------------------------------------------------------

# analysis ---------------------------------------------------------------------
# make output dir
system(paste0("mkdir -p ", output_dir))

# MAF files for OMNI data ------------------------------------------------------
# function to generate MAF files for a given grouping column
generate_maf_files <- function(group_column, data = annotated_mutations_progression, output_dir) {
  cat("Generating MAF files for:", group_column, "\n")
  
  # get unique groups
  all_samples <- unique(data[[group_column]])
  
  # generate MAF files
  for (sample in all_samples) {
    make_maf(sample, group_column, output_dir)
  }
  
  # verify file generation
  n_files <- list.files(path = paste0(output_dir, "/", group_column, "/"), pattern = "*.txt") %>% length()
  success <- n_files == length(all_samples)
  
  cat("Generated", n_files, "files for", length(all_samples), "samples -", 
      ifelse(success, "SUCCESS", "ERROR"), "\n\n")
  
  return(success)
}

# define grouping strategies
grouping_strategies <- c(
  "gh_request_id",
  "tmbh_status_compared_to_baseline", 
  "patient_status_compared_to_baseline"
)

# generate MAF files for all strategies
results <- sapply(grouping_strategies, function(x) generate_maf_files(x, output_dir = output_dir))

# MAF files for IMPACT data ----------------------------------------------------
# get all impact samples
all_impact_samples <- impact_mod$tumor_sample_barcode %>% unique()

results_impact <- sapply(all_impact_samples, function(x) make_maf_impact(x, sub_dir = "impact", output_dir = output_dir))

# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()