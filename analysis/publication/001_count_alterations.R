## -----------------------------------------------------------------------------
## Purpose of script: Count alterations
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
project_dir <- file.path(here::here(), "/", project_folder)

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# annotated mutations at progression
annotated_mutations_progression <- fread(paste0(project_dir, "/data/processed/005_annotated_mutations_progression.txt"))

# wrangle ----------------------------------------------------------------------
# annotated muts filter for acquired mutations
guardant_progression <- annotated_mutations_progression %>% 
  filter(status_compared_to_baseline_impact == "Acquired")

# analysis ---------------------------------------------------------------------
n_variant_type <- guardant_progression %>% 
  group_by(variant_type, acquired_tmbh) %>% 
  summarize(n_variant = n(), .groups = "drop") %>%
  group_by(acquired_tmbh) %>% 
  mutate(n_total = sum(n_variant))

mut_tmbh <- n_variant_type %>% filter(acquired_tmbh == "Y")

# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()