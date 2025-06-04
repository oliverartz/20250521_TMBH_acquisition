## -----------------------------------------------------------------------------
## Purpose of script: Stats for netMHCpan
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

p_load(tidyverse, data.table, janitor)

# set parameters ---------------------------------------------------------------
project_folder <- "20250521_TMBH_acquisition"
project_dir <- file.path(here::here(), project_folder)

threshold_weak_binder <- 2

# load data --------------------------------------------------------------------
source(file.path(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# netMHCpan results
patient_status_compared_to_baseline_phbr_vaf <- fread(file.path(project_dir, "data/processed/hpc/006_patient_status_compared_to_baseline_phbr_vaf.txt"))

# patient key
patient_meta_mod <- fread(paste0(project_dir, "/data/processed/003_patient_meta_mod.txt"))

# wrangle ----------------------------------------------------------------------
# get acquired mutations
acquired_mut <- patient_status_compared_to_baseline_phbr_vaf %>% 
  filter(status == "Acquired") %>% 
  mutate(binder_status = case_when(phbr < threshold_weak_binder ~ "Binder",
                                   TRUE ~ "Non-binder"))

# filter for patients with acquired TMB-H
acquired_mut <- acquired_mut %>% 
  filter(patient_id %in% patient_meta_mod$patient_id[patient_meta_mod$acquired_tmbh == "Y"])

# analysis ---------------------------------------------------------------------
# calculate proportion of binders
acquired_mut$binder_status %>% tabyl()


# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()