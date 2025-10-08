## -----------------------------------------------------------------------------
## Purpose of script: Export mutation tables for publication
##
## Author: Oliver Artz
## Date Created: 2025-09-17
## Email: artzo@mskcc.org
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

p_load(tidyverse, data.table, here, janitor)

# set parameters ---------------------------------------------------------------
project_name <- "20250521_TMBH_acquisition"
project_dir <- file.path(here(), project_name)

# load data --------------------------------------------------------------------
# IMPACT mutations
impact <- fread(file.path(
  project_dir,
  "/data/processed/004_impact_mod.txt"
))

# Guardant mutations
guardant <- fread(file.path(
  project_dir,
  "/data/processed/002_guardant_mod.txt"
))

# patient meta data
patient_meta <- fread(file.path(
  project_dir,
  "/data/processed/003_patient_meta_mod.txt"
))

# wrangle ----------------------------------------------------------------------
# patient meta data ------------------------------------------------------------
patient_meta_mod <- patient_meta

# guardant data ----------------------------------------------------------------
guardant_mod <- guardant

# only keep essential columns
guardant_mod <- guardant |>
  select(patient_id, gene, variant_type, mut_aa, percentage, gh_request_id)

# annotate sample type
guardant_mod <- guardant_mod |>
  mutate(
    sample_status = case_when(
      gh_request_id %in% patient_meta_mod$guardant_baseline ~ "baseline",
      gh_request_id %in% patient_meta_mod$guardant_progression ~ "progression",
      TRUE ~ "other"
    )
  ) |>
  arrange(patient_id) |>
  select(-gh_request_id) |>
  mutate(technology = "Guardant")

# impact data ------------------------------------------------------------------
impact_mod <- impact

# only keep essential columns
impact_mod <- impact_mod |>
  select(
    patient_id,
    gene = hugo_symbol,
    variant_type,
    mut_aa = hgv_sp_short,
    tumor_sample_barcode
  ) |>
  mutate(mut_aa = str_remove(mut_aa, "p."))

# annotate sample type
impact_mod <- impact_mod |>
  mutate(
    sample_status = case_when(
      tumor_sample_barcode %in% patient_meta_mod$impact_baseline ~ "baseline",
      tumor_sample_barcode %in% patient_meta_mod$impact_progression ~
        "progression",
      TRUE ~ "other"
    )
  ) |>
  arrange(patient_id) |>
  select(-tumor_sample_barcode) |>
  mutate(percentage = NA) |>
  mutate(technology = "IMPACT")

# combine both data sets
merged_muts <- rbind(guardant_mod, impact_mod) |>
  arrange(patient_id)

# QC ---------------------------------------------------------------------------

# export -----------------------------------------------------------------------
fwrite(
  merged_muts,
  file = file.path(project_dir, "data/processed/006_impact_guardant_muts.txt"),
  sep = "\t"
)
