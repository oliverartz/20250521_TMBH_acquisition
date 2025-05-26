## -----------------------------------------------------------------------------
## Purpose of script: Plot MSI score for each patient
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

p_load(tidyverse, data.table, ggsci, ggpubr)

# set parameters ---------------------------------------------------------------
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# load data --------------------------------------------------------------------

# wrangle ----------------------------------------------------------------------

# analysis ---------------------------------------------------------------------

# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()