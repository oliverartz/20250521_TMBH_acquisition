## -----------------------------------------------------------------------------
## Purpose of script: Preprocess IMPACT data
##
## Author: Oliver Artz
## Date Created: May 23, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, here)

# set parameters ---------------------------------------------------------------
project_root <- here()
project_folder <- "20250521_TMBH_acquisition"

# load data --------------------------------------------------------------------
# IMPACT mutations (cutoff 2024-12-06)
impact <- fread("/Users/artzo/Documents/8_asset_files/IMPACT/mutations/20241206_data_mutations_extended.txt.gz")

# patient meta data
patient_meta_mod <- fread(paste0(
  project_root, "/",
  project_folder, "/data/processed/003_patient_meta_mod.txt"))

# wrangle ----------------------------------------------------------------------
# clean names
impact_mod <- impact %>% clean_names() %>% remove_empty(which = "cols")

# keep relevant samples
samples <- patient_meta_mod %>% 
  filter(impact_include == "yes") %>% 
  select(impact_baseline, impact_progression) %>% 
  pivot_longer(cols = c(impact_baseline, impact_progression),
               names_to = "sample_status",
               values_to = "tumor_sample_barcode") 

impact_mod <- impact_mod %>% filter(tumor_sample_barcode %in% samples$tumor_sample_barcode)

# add mut_id
impact_mod <- impact_mod %>%
  mutate(mut_id_short = case_when(
    consequence == "frameshift_variant" ~ paste0(
      hugo_symbol,
      "_",
      str_remove(hgv_sp_short, "p.") %>%
        str_extract("\\D+\\d+"),
      "fs"
    ),
    TRUE ~ paste0(hugo_symbol, "_", str_remove(hgv_sp_short, "p."))
  ))

# add patient number
# make df with patient numbers
pt_num <- patient_meta_mod %>% 
  select(patient_id, impact_baseline, impact_progression) %>% 
  pivot_longer(cols = c(impact_baseline, impact_progression),
               names_to = "sample_status",
               values_to = "tumor_sample_barcode") %>% 
  filter(!tumor_sample_barcode == "None")

idx <- match(impact_mod$tumor_sample_barcode, pt_num$tumor_sample_barcode)
impact_mod$patient_id <- pt_num$patient_id[idx]

# export -----------------------------------------------------------------------
write.table(impact_mod, 
            file = paste0(project_root, "/",
                          project_folder, "/data/processed/004_impact_mod.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# cleanup ----------------------------------------------------------------------
rm()