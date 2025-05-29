## -----------------------------------------------------------------------------
## Purpose of script: Plot mean number of variants for each patient at baseline
##                    and progression. Guardant samples.
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

# order of variant types
order_variant_types <- c("SNV", "Indel", "CNV", "Fusion", "LGR")

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# guardant data
guardant_mod <- fread(paste0(project_folder, "/data/processed/002_guardant_mod.txt"))

# patient metadata
patient_meta <- fread(paste0(project_folder, "/data/processed/003_patient_meta_mod.txt"))

# wrangle ----------------------------------------------------------------------
# refactor sample status
guardant_mod <- guardant_mod %>% 
  mutate(sample_status = case_when(sample_status == "baseline" ~ "Baseline",
                                   sample_status == "progression" ~ "Progression",
                                   TRUE ~ sample_status)) %>% 
  mutate(acquired_tmbh = case_when(acquired_tmbh == "N" ~ "No",
                                   acquired_tmbh == "Y" ~ "Yes",
                                   TRUE ~ acquired_tmbh))

# reorder variant types for plotting
guardant_mod$variant_type <- factor(guardant_mod$variant_type, levels = order_variant_types)

# analysis ---------------------------------------------------------------------
# count variants per time point and patient
n_sample_status_variant_type <- guardant_mod %>% 
  group_by(patient_id, sample_status, variant_type) %>% 
  summarize(n = n(), .groups = "drop") %>% 
  complete(patient_id, variant_type, sample_status, fill = list(n = 0)) %>% 
  mutate(plot_group = paste0(patient_id, "_", variant_type)) 

# annotate tmbh acquisition
idx <- match(n_sample_status_variant_type$patient_id, patient_meta$patient_id)
n_sample_status_variant_type$acquired_tmbh <- patient_meta$acquired_tmbh[idx]

# get average number of variants per sample
mean_sample_status_variant_type <- guardant_mod %>% 
  group_by(sample_status, variant_type, acquired_tmbh) %>% 
  summarize(n = n())

n_samples <- guardant_mod %>%
  distinct(paste0(patient_id, sample_status), .keep_all = TRUE) %>% 
  group_by(sample_status, acquired_tmbh) %>% 
  summarize(n_samples = n()) 

mean_sample_status_variant_type <- mean_sample_status_variant_type %>%
  mutate(merge = paste0(sample_status, acquired_tmbh)) %>% 
  left_join(n_samples %>% 
              mutate(merge = paste0(sample_status, acquired_tmbh)),
            by = "merge") %>% 
  mutate(mean_variant_type = n / n_samples) %>% 
  dplyr::rename("sample_status" = "sample_status.x",
                "acquired_tmbh" = "acquired_tmbh.x")

# plot -------------------------------------------------------------------------
# make df for plotting
df_plot <- mean_sample_status_variant_type

# plot
df_plot %>% 
  ggplot(aes(x = sample_status, y = mean_variant_type)) +
  geom_line(aes(color = variant_type, group = variant_type), alpha = 0.5, size = 1) +
  geom_point(shape = 21, size = 3, aes(fill = variant_type)) +
  geom_text(data = n_samples, aes(label = paste0("n=" , n_samples), y = -5)) +
  theme_jco() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  scale_color_jco() +
  scale_fill_jco() +
  labs(
    #title = "Mean number of variants per patient",
    x = "",
    y = "Mean number of variants",
    color = "Variant type",
    fill = "Variant type"
  ) +
  facet_wrap(~acquired_tmbh)

# N denotes the number of samples that were used to calculate the mean, not the
# number of data points in the plot

# export -----------------------------------------------------------------------
# plot
ggsave(paste0(project_folder, "/results/figures/014_mean_n_variants_per_patient.pdf"), height = 6, width = 6)

# data
df_plot %>% 
  select(-c(acquired_tmbh.y, merge, sample_status.y)) %>% 
  fwrite(
    paste0(project_folder, "/results/tables/014_mean_n_variants_per_patient.txt"), 
    sep = "\t"
  )

# cleanup ----------------------------------------------------------------------
rm()