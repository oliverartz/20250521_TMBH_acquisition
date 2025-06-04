## -----------------------------------------------------------------------------
## Purpose of script: Clarify difference in number of acquired mutations between
##                    017 and 018
##
## Author: Oliver Artz
## Date Created: Jun 4, 2025
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
project_folder <- "20250521_TMBH_acquisition"
project_dir <- file.path(here::here(), project_folder)

# load data --------------------------------------------------------------------
# data from 017
annotated_mutations_progression <- fread(paste0(project_dir, "/data/processed/005_annotated_mutations_progression.txt"))

# data from 018
mut_mod_018 <- fread(file.path(project_dir, "/results/tables/018_mutations_with_corrected_vaf.txt"))

# patient metadata
patient_meta <- fread(paste0(project_dir, "/data/processed/003_patient_meta_mod.txt"))

# wrangle ----------------------------------------------------------------------
annotated_mutations_progression_mut <- annotated_mutations_progression %>% 
  mutate(pt_chrom_pos = paste0(patient_id, "_", chromosome, "_", position))

acquired_017 <- annotated_mutations_progression_mut %>% 
  filter(status_compared_to_baseline_impact == "Acquired") %>% 
  filter(patient_id %in% patient_meta$patient_id[patient_meta$acquired_tmbh == "Y"]) %>% 
  filter(variant_type %in% c("SNV", "Indel")) 

acquired_018 <- mut_mod_018 %>% 
  filter(status_compared_to_baseline_impact == "Acquired")

# analysis ---------------------------------------------------------------------
# Which mutations in 018 are not in 017?
mut_018 <- acquired_018 %>% 
  mutate(in_017 = case_when(pt_chrom_pos %in% acquired_017$pt_chrom_pos ~ "Yes",
                            TRUE ~ "No")) %>% 
  filter(in_017 == "No")

df_temp <- annotated_mutations_progression_mut %>% filter(pt_chrom_pos %in% mut_018$pt_chrom_pos)

df_temp_with_vaf <- df_temp %>% filter(!is.na(percentage))

# The missing mutations are Fusions and LGR, that had a reported VAF

# cleanup ----------------------------------------------------------------------
rm()