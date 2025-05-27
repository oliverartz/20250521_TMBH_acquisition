## -----------------------------------------------------------------------------
## Purpose of script: Plot fractions of variant types
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

p_load(tidyverse, data.table)

# set parameters ---------------------------------------------------------------
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# order of variant types
order_variant_type <- c("Missense", "Frameshift", "Nonsense", "In Frame", "Other")

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# IMPACT mutations
impact_mut_mod <- fread(paste0(project_dir, "/data/processed/004_impact_mod.txt"))

# wrangle ----------------------------------------------------------------------

# analysis ---------------------------------------------------------------------
# count variant types
count_variant_types <- impact_mut_mod %>% 
  count(variant_type, patient_id, acquired_tmbh) %>% 
  complete(variant_type, patient_id, fill = list(n = 0))

# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()