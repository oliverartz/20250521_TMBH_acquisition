## -----------------------------------------------------------------------------
## Purpose of script: Acquired MHC binders stats
##
## Author: Oliver Artz
## Date Created: Jun 14, 2025
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

# thresholds to categorize binders
thresh_strong_binder <- 0.5
thresh_weak_binder <- 2

# load data --------------------------------------------------------------------
# PHBR data
phbr_vaf_mod <- fread(paste0(project_dir, "/results/tables/019_phbr_acquired_baseline_data.txt"))

# wrangle ----------------------------------------------------------------------
# calculate percent of mutations that are (weak) binders
phbr_vaf_mod <- phbr_vaf_mod %>% 
  mutate(binder_cat = case_when(phbr < thresh_weak_binder ~ "Binder",
                                TRUE ~ "Non binder"))

# get total mutations
total_mut <- phbr_vaf_mod %>% 
  group_by(patient_id, status) %>% 
  summarize(n_total = n()) %>% 
  ungroup() %>% 
  complete(patient_id, status, fill = list(n_total = 0))

# get binders
total_binders <- phbr_vaf_mod %>% 
  group_by(patient_id, status) %>% 
  summarize(n_binders = sum(binder_cat == "Binder")) %>% 
  ungroup() %>% 
  complete(patient_id, status, fill = list(n_binders = 0))

# merge and calculate percent of total
merge <- total_mut %>% 
  left_join(total_binders,
            by = c("patient_id", "status")) %>% 
  mutate(percent_binder = round(n_binders / n_total * 100, 0))

# get only acquired mutations
merge_acquired <- merge %>% 
  filter(status == "Acquired")

# summary stats
merge_acquired$percent_binder %>% summary()

# cleanup ----------------------------------------------------------------------
rm()