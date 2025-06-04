## -----------------------------------------------------------------------------
## Purpose of script: Calculate stats for clonality
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

p_load(tidyverse, data.table, janitor)

# set parameters ---------------------------------------------------------------
project_folder <- "20250521_TMBH_acquisition"
project_dir <- file.path(here::here(), project_folder)

# clonality threshold
subclonal_thresh <- 50
ultrasubclonal_thresh <- 20

# load data --------------------------------------------------------------------
mut_mod <- fread(file.path(project_dir, "/results/tables/018_mutations_with_corrected_vaf.txt"))

# patient key
patient_meta_mod <- fread(paste0(project_dir, "/data/processed/003_patient_meta_mod.txt"))

# annotated progression samples
annotated_mutations_progression <- fread(paste0(project_dir, 
                                                "/data/processed/005_annotated_mutations_progression.txt"))

# wrangle ----------------------------------------------------------------------
# annotate patient ID in mut_mod
idx <- match(mut_mod$gh_request_id, annotated_mutations_progression$gh_request_id)
mut_mod$patient_id <- annotated_mutations_progression$patient_id[idx]
mut_mod$acquired_tmbh <- annotated_mutations_progression$acquired_tmbh[idx]

# analysis ---------------------------------------------------------------------
# calculate fraction of subclonal acquired alterations
acquired <- mut_mod %>% 
  filter(status_compared_to_baseline_impact == "Acquired")

# filter for patients with acquired TMB-H
acquired <- acquired %>% 
  filter(acquired_tmbh == "Y")

# calculate fraction of subclonal
acquired <- acquired %>% 
  mutate(subclonal = case_when(percentage_corr < subclonal_thresh ~ "Subclonal",
                               TRUE ~ "Clonal"))

acquired$subclonal %>% tabyl()

# calculate fraction of ultrasubclonal within subclonal
ultrasubclonal <- acquired %>% 
  filter(subclonal == "Subclonal") %>% 
  mutate(ultrasubclonal = case_when(percentage_corr < ultrasubclonal_thresh ~ "Ultrasubclonal",
                                    TRUE ~ "Subclonal"))

ultrasubclonal$ultrasubclonal %>% tabyl()

# calculate number of subclonal mutations per patient
subclonal_per_patient <- acquired %>% 
  group_by(patient_id) %>% 
  summarize(total_mut = n(),
            n_subclonal = sum(subclonal == "Subclonal"),
            perc_subclonal = n_subclonal / total_mut * 100) %>% 
  ungroup()

summary(subclonal_per_patient$total_mut)
summary(subclonal_per_patient$perc_subclonal)

# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()