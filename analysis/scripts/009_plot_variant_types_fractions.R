## -----------------------------------------------------------------------------
## Purpose of script: Plot fractions of variant types derived from IMPACT
##                    baseline samples
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

p_load(tidyverse, data.table, ggpubr, ggsci)

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

# patient key
patient_key_mod <- fread(paste0(project_dir, "/data/processed/001_patient_key_mod.txt"))

# wrangle ----------------------------------------------------------------------
# filter for relevant impact baseline samples
impact_mut_mod <- impact_mut_mod %>% 
  filter(tumor_sample_barcode %in% patient_key_mod$impact_id_baseline)

# add new categories for consequence
impact_mut_mod <- impact_mut_mod %>% 
  mutate(variant_type = case_when(consequence == "frameshift_variant" ~ "Frameshift",
                                  consequence %in% c("inframe_deletion", "inframe_insertion") ~ "In Frame",
                                  consequence %in% c("missense_variant", "missense_variant,splice_region_variant") ~ "Missense",
                                  consequence %in% c("stop_gained") ~ "Nonsense",
                                  consequence %in% c("splice_acceptor_variant,coding_sequence_variant,intron_variant",
                                                     "splice_donor_variant", 
                                                     "splice_donor_variant,coding_sequence_variant,intron_variant",
                                                     "stop_gained,splice_region_variant") ~ "Other",
                                  TRUE ~ "Not categorized")
  )

# analysis ---------------------------------------------------------------------
# count variant types
count_variant_types <- impact_mut_mod %>% 
  count(variant_type, patient_id) %>% 
  complete(variant_type, patient_id, fill = list(n = 0))

# add acquired high tmb info back in
idx <- match(count_variant_types$patient_id, impact_mut_mod$patient_id)
count_variant_types$acquired_high_tmb_y_yes_n_no <- impact_mut_mod$acquired_tmbh[idx]

# plot -------------------------------------------------------------------------
# fraction of variant types
fraction_variant_types <- count_variant_types %>% 
  group_by(patient_id) %>% 
  mutate(frac = n / sum(n)) %>% 
  ungroup()

# reorder variant types
fraction_variant_types$variant_type <- factor(fraction_variant_types$variant_type,
                                              levels = order_variant_type)

# plot
fraction_variant_types %>% 
  ggplot(aes(x = acquired_high_tmb_y_yes_n_no, y = frac, fill = acquired_high_tmb_y_yes_n_no)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(shape = 21, size = 3, alpha = 0.8) +
  stat_compare_means(label.y = 0.95, label = "p.format") +
  facet_grid(~ variant_type) +
  theme_jco() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none") +
  scale_fill_jco() +
  labs(
    #title = "Fraction of Variant Types",
    x = "Acquired High TMB",
    y = "Fraction of Variants",
    fill = "Acquired High TMB")

# export -----------------------------------------------------------------------
# plot
ggsave(paste0(project_folder, "/results/figures/009_variant_types_fractions.pdf"), height = 6, width = 6)

# export data
fraction_variant_types %>% 
  select(
    patient_id,
    acquired_high_tmb_y_yes_n_no,
    variant_type,
    frac
  ) %>% 
  fwrite(
    paste0(project_folder, "/results/tables/009_variant_types_fractions.txt"), 
    sep = "\t"
  )

# cleanup ----------------------------------------------------------------------
rm()