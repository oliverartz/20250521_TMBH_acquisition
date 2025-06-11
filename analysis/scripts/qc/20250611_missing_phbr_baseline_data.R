## -----------------------------------------------------------------------------
## Purpose of script: Why do pts 22 and 26 have no PHBR baseline data?
##
## Author: Oliver Artz
## Date Created: Jun 11, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, here, janitor)

# set parameters ---------------------------------------------------------------
project_name <- "20250521_TMBH_acquisition"
project_dir <- file.path(here(), project_name)

# load data --------------------------------------------------------------------
# acquired mutations
annotated_mutations_progression <- fread(file.path(project_dir, "data/processed/005_annotated_mutations_progression.txt"))

# wrangle ----------------------------------------------------------------------
# filter for patients in question
acquired_mut_mod <- annotated_mutations_progression %>% 
  filter(patient_id %in% c(22, 26)) %>% 
  filter(status_compared_to_baseline_impact == "Baseline")

# There are no 'baseline' mutations for these two patients
# They might have all been immunoedited

# analysis ---------------------------------------------------------------------

# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()