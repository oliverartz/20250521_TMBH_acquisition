## -----------------------------------------------------------------------------
## Purpose of script: Deconvolve SBS signatures for each IMPACT baseline sample
##
## Author: Oliver Artz
## Date Created: May 28, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, MutationalPatterns)

# set parameters ---------------------------------------------------------------
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# reference genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# signature etiology
sig_eti <- fread(paste0(project_dir, "/data/metadata/20230505_mut_signature_aetiology.csv"))

# patient info
patient_info <- read_excel(paste0(project_folder, "/data/20250123_Full Cohort Data Table_for Oliver_20240123.xlsx"))

# wrangle ----------------------------------------------------------------------

# analysis ---------------------------------------------------------------------

# get COSMIC signatures
signatures <- get_known_signatures(muttype = "snv",
                                   source = "COSMIC_v3.2",
                                   genome = "GRCh37")

# plot -------------------------------------------------------------------------

# export -----------------------------------------------------------------------

# cleanup ----------------------------------------------------------------------
rm()