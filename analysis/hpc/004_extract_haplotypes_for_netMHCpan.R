## -----------------------------------------------------------------------------
## Purpose of script: Extract HLA haplotypes from relevant patients and format
##                    for netMHCpan
##
## Author: Oliver Artz
## Date Created:
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table)

# set parameters ---------------------------------------------------------------
project_root <- here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/../", project_folder)


output_dir <- paste0(project_dir, "/data/processed/hpc/")

# load data --------------------------------------------------------------------
# all haplotypes
all_haplotypes <- fread("/data/ldiaz/artzo/assets/IMPACT/HLA/impact_hla_manifest_2025_01_16.txt")

# patient key
patient_key_mod <- fread(paste0(project_dir, "/data/processed/001_patient_key_mod.txt"))

# wrangle ----------------------------------------------------------------------
# filter haplotype data for relevant patients
hla_haplotypes_mod <- data.frame(patient = patient_key_mod$dmp_patient_id) %>% 
  left_join(all_haplotypes,
            by = c("patient")) %>% 
  distinct(patient, .keep_all = TRUE)

# analysis ---------------------------------------------------------------------

# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()