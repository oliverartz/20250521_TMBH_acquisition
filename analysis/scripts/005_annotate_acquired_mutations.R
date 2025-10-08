## -----------------------------------------------------------------------------
## Purpose of script: Annotate baseline, acquired, and lost mutations
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
if (!require('pacman')) {
  install.packages('pacman', dependencies = TRUE)
}
library(pacman)

p_load(tidyverse, data.table, here)

# set parameters ---------------------------------------------------------------
project_root <- here()
project_folder <- "20250521_TMBH_acquisition"

# load data --------------------------------------------------------------------
# load function for annnotation
source(paste0(
  project_root,
  "/",
  project_folder,
  "/analysis/scripts/utils/001_preprocess_utils.R"
))

# Guardant OMNI mutations
guardant_mod <- fread(paste0(
  project_root,
  "/",
  project_folder,
  "/data/processed/002_guardant_mod.txt"
))

# IMPACT mutations
impact_mod <- fread(paste0(
  project_root,
  "/",
  project_folder,
  "/data/processed/004_impact_mod.txt"
))

# panel gene list
panel_gene_list <- fread(paste0(
  project_root,
  "/",
  project_folder,
  "/data/metadata/20250108_all_gene_list.csv"
))

# analysis ---------------------------------------------------------------------
# get all unique patient IDs
unique_patient_ids <- unique(guardant_mod$patient_id)

# run the function
annotated_mutations <- annotate_patient_mutations(
  guardant_mod,
  impact_mod,
  panel_gene_list,
  unique_patient_ids
)

# only keep unique mut_ids per patient
annotated_mutations_mod <- annotated_mutations %>%
  distinct(patient_id, mut_id_short, sample_status, .keep_all = TRUE)

# add columns for maf2vcf
annotated_mutations_mod <- annotated_mutations_mod %>%
  mutate(
    reference_allele = str_split_fixed(mut_nt, ">", 2)[, 1],
    tumor_seq_allele2 = str_split_fixed(mut_nt, ">", 2)[, 2],
    tumor_sample_barcode = gh_request_id
  )

# keep full data
annotated_mutations_all <- annotated_mutations_mod %>%
  mutate(
    patient_status_compared_to_baseline = paste0(
      patient_id,
      "_",
      status_compared_to_baseline_impact
    ),
    tmbh_status_compared_to_baseline = paste0(
      acquired_tmbh,
      "_",
      status_compared_to_baseline_impact
    )
  )

# filter for progression samples for each patient
annotated_mutations_progression <- annotated_mutations_all %>%
  filter(sample_status %in% c("progression"))

# export -----------------------------------------------------------------------
# these are only the progression samples. there will be no lost mutations in here, because lost mutations are in the IMPACT baseline
write.table(
  annotated_mutations_progression,
  file = paste0(
    project_root,
    "/",
    project_folder,
    "/data/processed/005_annotated_mutations_progression.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# this dataset contains all mutations from all sample types (i.e. it has the lost mutations from IMPACT)
write.table(
  annotated_mutations_all,
  file = paste0(
    project_root,
    "/",
    project_folder,
    "/data/processed/005_annotated_mutations_all.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# make reduced data set to upload with the publication
upload_df <- annotated_mutations_all |>
  select(
    patient_id,
    gene,
    variant_type,
    mut_aa,
    sample_status,
    status_compared_to_baseline_impact,
    percentage
  ) |>
  filter(!status_compared_to_baseline_impact == "Lost") |>
  arrange(patient_id) |>
  mutate(
    status_compared_to_baseline_impact = case_when(
      variant_type == "CNV" ~ NA,
      TRUE ~ status_compared_to_baseline_impact
    )
  )

write.table(
  upload_df,
  file = paste0(
    project_root,
    "/",
    project_folder,
    "/data/processed/005_progression_annotated.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# cleanup ----------------------------------------------------------------------
rm()
