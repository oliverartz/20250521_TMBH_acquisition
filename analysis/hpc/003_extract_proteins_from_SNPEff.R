## -----------------------------------------------------------------------------
## Purpose of script: Extract protein sequences from SNPEff-annotated VCF files
##
## Author: Oliver Artz
## Date Created: May 29, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, vcfR, phylotools, here)

# set parameters ---------------------------------------------------------------
project_root <- here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/../", project_folder)

processed_data_dir <- paste0(project_dir, "/data/processed/hpc/")
snpeff_output_dir <- paste0(processed_data_dir, "snpeff/patient_status_compared_to_baseline")
fasta_output_dir <- paste0(processed_data_dir, "fasta/patient_status_compared_to_baseline/")

# load data --------------------------------------------------------------------
# get list of all fasta files in output dir
fasta_files <- list.files(snpeff_output_dir, pattern = "protein_seq.fasta", full.names = TRUE)

# read all fasta files into list
fasta <- lapply(fasta_files, read.fasta)

# get list of all vcf files in output dir
vcf_files <- list.files(snpeff_output_dir, pattern = "_annotated.vcf.tsv", full.names = TRUE)

## read all vcf files
vcf <- lapply(vcf_files, fread)

# wrangle ----------------------------------------------------------------------

# analysis ---------------------------------------------------------------------

# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()