## -----------------------------------------------------------------------------
## Purpose of script: Plot blood TMB vs tissue TMB
##
## Author: Oliver Artz
## Date Created: May 31, 2025
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
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# patient key
patient_key_mod <- fread(paste0(project_dir, "/data/processed/001_patient_key_mod.txt"))

# wrangle ----------------------------------------------------------------------
# make long format for plotting
patient_key_mod_long <- patient_key_mod %>% 
  mutate(baseline_t_tmb_mut_mb_18 = as.numeric(baseline_t_tmb_mut_mb_18), 
         baseline_p_tmb_mut_mb_19 = as.numeric(baseline_p_tmb_mut_mb_19),
         progression_t_tmb_mut_mb_20 = as.numeric(progression_t_tmb_mut_mb_20), 
         progression_p_tmb_mut_mb_21 = as.numeric(progression_p_tmb_mut_mb_21)) %>% 
  pivot_longer(cols = c(baseline_t_tmb_mut_mb_18, baseline_p_tmb_mut_mb_19,
                        progression_t_tmb_mut_mb_20, progression_p_tmb_mut_mb_21),
               names_to = "sample",
               values_to = "tmb") 

# add sample status
patient_key_mod_long <- patient_key_mod_long %>% 
  separate(sample, sep = "_", into = c("sample_status", "tissue")) %>% 
  mutate(tissue = case_when(tissue == "t" ~ "Tissue",
                            tissue == "p" ~ "Plasma",
                            TRUE ~ "Not defined"),
         sample_status = str_to_title(sample_status),
         acquired_high_tmb_y_yes_n_no = case_when(acquired_high_tmb_y_yes_n_no == "N" ~ "No",
                                                  acquired_high_tmb_y_yes_n_no == "Y" ~ "Yes",
                                                  TRUE ~ "Not defined"))

# analysis ---------------------------------------------------------------------

# plot -------------------------------------------------------------------------
# make df for plotting
df_plot <- patient_key_mod_long %>% 
  mutate(patient_status = paste0(patient_id, sample_status))

# plot
df_plot %>% 
  ggplot(aes(x = tissue, y = tmb)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_line(aes(group = patient_status,
                color = sample_status)) +
  geom_point(
    aes(fill = sample_status),
    size = 3, shape = 21, alpha = 0.8) +
  stat_compare_means(label.y = 140) +
  theme_jco() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  scale_fill_jco() +
  scale_color_jco() +
  labs(
    x = "",
    y = "TMB (mut/Mb)",
    fill = "Sample Status",
    color = "Sample Status",
    group = "Patient ID") +
  facet_wrap(~ acquired_high_tmb_y_yes_n_no)

# export -----------------------------------------------------------------------
ggsave(paste0(project_dir, "/results/figures/021_blood_tissue_tmb.pdf"), 
       width = 6, height = 6)

# cleanup ----------------------------------------------------------------------
rm()