## -----------------------------------------------------------------------------
## Purpose of script: Preprocess Guardant OMNI data
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

p_load(tidyverse, data.table, here, janitor)

# set parameters ---------------------------------------------------------------
project_root <- here()
project_folder <- "20250521_TMBH_acquisition"

# load data --------------------------------------------------------------------
guardant_1 <- fread(paste0(project_root, "/", project_folder, "/data/raw/guardant_omni/Yaeger G360 genes updated.csv 2025-01-08 .csv"))
guardant_2 <- fread(paste0(project_root, "/", project_folder, "/data/raw/guardant_omni/Yaeger G360 genes updated.csv 2025-01-28 .csv"))
guardant_3 <- fread(paste0(project_root, "/", project_folder, "/data/raw/guardant_omni/Yaeger G360 genes.csv 2024-12-10 .csv"))
guardant_4 <- fread(paste0(project_root, "/", project_folder, "/data/raw/guardant_omni/Yaeger OMNI.csv 2025-05-14 .csv"))

patient_meta <- fread(paste0(project_root, "/", project_folder, "/data/metadata/20250521_patient_meta.csv"))

# wrangle ----------------------------------------------------------------------
# make patient_meta_mod_long for merging with guardant ids
patient_meta_mod <- patient_meta %>% clean_names()

# add impact baseline version
patient_meta_mod <- patient_meta_mod %>% 
  mutate(baseline_impact_version = str_sub(impact_baseline, -3, -1))

# make long
patient_meta_mod_long <- patient_meta_mod %>%
  select(patient_id, guardant_baseline, guardant_progression, acquired_tmbh, baseline_impact_version) %>% 
  pivot_longer(cols = c(guardant_baseline, guardant_progression),
               names_to = "sample_status",
               values_to = "accession_id") %>% 
  mutate(sample_status = str_remove(sample_status, "guardant_"))


# merge guardant files
guardant <- rbind(guardant_1,
                  guardant_2,
                  guardant_3 %>% select(-c(`Patient ID`, `Sample Status`)),
                  guardant_4)

# fix column names
guardant_mod <- guardant %>% clean_names()

# add patient id
guardant_mod <- guardant_mod %>% 
  left_join(patient_meta_mod_long, 
            by = c("gh_request_id" = "accession_id")) %>% 
  dplyr::select(patient_id, everything())

# filter for somatic mutations
guardant_mod <- guardant_mod %>% 
  filter(grepl("somatic", somatic_status))

# select relevant columns
guardant_mod <- guardant_mod %>% 
  dplyr::select(-c(splice_effect, fusion_chrom_b, fusion_gene_b, fusion_position_a, 
                   fusion_position_b, direction_a, direction_b, downstream_gene, 
                   cosmic, db_snp, received_date, bloodcoll_date, reported_date, 
                   practice_name, physician_name))

# add mutation id
guardant_mod <- guardant_mod %>% 
  mutate(mut_id = paste0(gene, "_", position, "_", mut_cdna, "_", variant_type))

# make mut_id_short
guardant_mod <- guardant_mod %>% 
  mutate(mut_id_short = paste0(gene, "_", mut_aa))

# remove guardant samples we decided to exclude
samples_to_exclude <- c(
  "A0426213", # baseline sample from chemo patient - we have the chemo baseline sample instead
  "A0326263" # On-treatment sample from former patient 5
)

guardant_mod <- guardant_mod %>% 
  filter(!gh_request_id %in% samples_to_exclude)

# export -----------------------------------------------------------------------
write.table(guardant_mod, 
            file = paste0(
              project_root, "/",
              project_folder, "/data/processed/002_guardant_mod.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# cleanup ----------------------------------------------------------------------
rm()