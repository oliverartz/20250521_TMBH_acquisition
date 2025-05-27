## -----------------------------------------------------------------------------
## Purpose of script: Preprocess patient key
##
## Author: Oliver Artz
## Date Created: May 21, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, readxl, here, janitor)

# set parameters ---------------------------------------------------------------
project_root <- here()
project_folder <- "20250521_TMBH_acquisition"

# load data --------------------------------------------------------------------
# patient key
patient_key <- read_xlsx(paste0(
  project_root, "/",
  project_folder, 
  "/data/metadata/20250123_Full Cohort Data Table_for Oliver_20240123.xlsx"))

# clinical IMPACT data
impact_clin <- fread("/Users/artzo/Documents/8_asset_files/IMPACT/clinical/20250206_data_clinical_sample.txt",
                     skip = 4)


# wrangle ----------------------------------------------------------------------
patient_key_mod <- patient_key %>% clean_names()

# reformat to computation-friendly
patient_key_mod <- patient_key_mod %>% 
  filter(!is.na(irb)) %>% 
  select(patient_id, guardant_id_baseline, guardant_id_progression, 
         acquired_high_tmb_y_yes_n_no, impact_id_baseline,
         baseline_t_tmb_mut_mb_18, baseline_p_tmb_mut_mb_19,
         progression_t_tmb_mut_mb_20, progression_p_tmb_mut_mb_21) %>% 
  mutate(dmp_patient_id = str_sub(impact_id_baseline, 1, 9))

# clean up baseline id
patient_key_mod <- patient_key_mod %>% 
  mutate(guardant_id_baseline = case_when(guardant_id_baseline == "A0540760 (Guardant CDx ONLY)" ~ "A0540760",
                                          TRUE ~ guardant_id_baseline))

# add IMPACT MSIScore
idx <- match(patient_key_mod$impact_id_baseline, impact_clin$SAMPLE_ID)
patient_key_mod$msi_score <- impact_clin$MSI_SCORE[idx]

# export -----------------------------------------------------------------------
write.table(patient_key_mod, 
            file = paste0(
              project_root, "/",
              project_folder, "/data/processed/001_patient_key_mod.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)
