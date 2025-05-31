## -----------------------------------------------------------------------------
## Purpose of script: Plot acquired alteration types at progression. Guardant
##                    samples.
##
## Author: Oliver Artz
## Date Created: May 29, 2025
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

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# annotated 
annotated_mutations_progression <- fread(paste0(project_dir, "/data/processed/005_annotated_mutations_progression.txt"))

# patient metadata
patient_meta <- fread(paste0(project_dir, "/data/processed/003_patient_meta_mod.txt"))

# wrangle ----------------------------------------------------------------------
# annotated muts filter for acquired mutations
guardant_progression <- annotated_mutations_progression %>% 
  filter(status_compared_to_baseline_impact == "Acquired")

# analysis ---------------------------------------------------------------------
# number of total acquired mutations in progression samples --------------------

# variant types
n_variant_type <- guardant_progression %>% 
  group_by(variant_type, acquired_tmbh) %>% 
  summarize(n_variant = n(), .groups = "drop") %>%
  group_by(acquired_tmbh) %>% 
  mutate(n_total = sum(n_variant))

# variant subtypes -------------------------------------------------------------
## SNV
snv_subtype <- guardant_progression %>% 
  filter(variant_type == "SNV") %>% 
  mutate(subtype = case_when(is.na(mut_aa) ~ "Other",
                             str_sub(mut_aa, -1, -1) == "*" ~ "Nonsense",
                             str_detect(str_sub(mut_aa, -1, -1), "[A-Za-z]") ~ "Missense",
                             TRUE ~ "Other")) %>% 
  group_by(subtype, acquired_tmbh) %>% 
  count() %>% 
  mutate(type = "SNV")

## Indel
indel_subtype <- guardant_progression %>% 
  filter(variant_type == "Indel") %>% 
  group_by(indel_type, acquired_tmbh) %>% 
  count() %>% 
  mutate(type = "Indel") %>% 
  dplyr::rename("subtype" = "indel_type")

## CNV
## The CNV status in the detailed report has many subtypes. The category of the 
## subtype is sometimes not in line with the reported copy number, for example
## n = 5, but loh_deletion. After discussing with the team, we decided to 
## simplify the subtypes

cnv_subtype <- guardant_progression %>% 
  filter(variant_type == "CNV") %>% 
  mutate(subtype = case_when(amplification_type == "amplification" ~ "Amplification",
                             amplification_type == "aneuploid_amplification" ~ "Amplification",
                             amplification_type == "aneuploidy" ~ "Amplification",
                             amplification_type == "deletion" ~ "Deletion",
                             amplification_type == "focal" ~ "Amplification",
                             amplification_type == "focal_amplification" ~ "Amplification",
                             amplification_type == "homozygous_deletion" ~ "Deletion",
                             amplification_type == "loh_deletion" ~ "Deletion",
                             TRUE ~ "Other"
  )) 

cnv_subtype_count <- cnv_subtype %>% 
  filter(!subtype == "Other") %>% 
  group_by(subtype, acquired_tmbh) %>% 
  count() %>% 
  mutate(type = "CNV")

# merge counts
merge_subtype <- rbind(snv_subtype, indel_subtype, cnv_subtype_count)

# add type counts
merge_subtype <- merge_subtype %>% 
  left_join(n_variant_type,
            by = c("type" = "variant_type", "acquired_tmbh"))

# calculate fraction
merge_subtype <- merge_subtype %>% 
  mutate(frac_subtype = n / n_variant)

# refactor acquired tmb
merge_subtype <- merge_subtype %>% 
  mutate(acquired_tmbh = case_when(acquired_tmbh == "Y" ~ "Yes",
                                   acquired_tmbh == "N" ~ "No",
                                   TRUE ~ "Other"))


# plot -------------------------------------------------------------------------
# call plotting function
p_indel_n <- plot_variant_type("Indel", "No")
p_indel_y <- plot_variant_type("Indel", "Yes")

p_snv_n <- plot_variant_type("SNV", "No")
p_snv_y <- plot_variant_type("SNV", "Yes")

p_cnv_n <- plot_variant_type("CNV", "No")
p_cnv_y <- plot_variant_type("CNV", "Yes")

cowplot::plot_grid(p_snv_y, p_indel_y, p_cnv_y, 
                   p_snv_n, p_indel_n, p_cnv_n,
                   nrow = 2)

# export -----------------------------------------------------------------------
# plot
ggsave(paste0(project_dir, "/results/figures/017_acquired_alteration_types_progression.pdf"), 
       height = 6, width = 8)

# data
write.table(merge_subtype, 
            paste0(project_dir, "/results/tables/017_acquired_alteration_types_progression.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# cleanup ----------------------------------------------------------------------
rm()