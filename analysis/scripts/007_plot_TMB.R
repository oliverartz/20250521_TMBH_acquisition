## -----------------------------------------------------------------------------
## Purpose of script: Plot TMB for all patients
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

p_load(tidyverse, data.table, ggsci, ggpubr)

# set parameters ---------------------------------------------------------------
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# patient data
patient_key_mod <- fread(paste0(project_dir, "/data/processed/001_patient_key_mod.txt"))

# wrangle ----------------------------------------------------------------------
# refactor acquired_high_tmb_y_yes_n_no
patient_key_mod <- patient_key_mod %>% 
  filter(!is.na(acquired_high_tmb_y_yes_n_no)) %>% 
  filter(!impact_id_baseline == "None") %>% 
  mutate(acquired_high_tmb_y_yes_n_no = case_when(acquired_high_tmb_y_yes_n_no == "N" ~ "No",
                                                  acquired_high_tmb_y_yes_n_no == "Y" ~ "Yes",
                                                  TRUE ~ "Other"))

# plot -------------------------------------------------------------------------
# make df for plotting
df_plot <- patient_key_mod

# plot
df_plot %>% 
  ggplot(aes(x = acquired_high_tmb_y_yes_n_no, 
             y = baseline_t_tmb_mut_mb_18, 
             fill = acquired_high_tmb_y_yes_n_no)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(shape = 21, size = 3, alpha = 0.8) +
  stat_compare_means(label.x = 1.5) +
  theme_jco() +
  theme(legend.position = "none") +
  scale_fill_jco() +
  labs(
    #title = "TMB",
    x = "Acquired High TMB",
    y = "TMB (mut/Mb)",
    fill = "Acquired High TMB")

# export -----------------------------------------------------------------------
# plot
ggsave(paste0(project_folder, "/results/figures/007_TMB.pdf"), height = 6, width = 4)

# export data
df_plot %>%
  select(
    patient_id,
    acquired_high_tmb_y_yes_n_no,
    baseline_t_tmb_mut_mb_18
  ) %>%
  fwrite(
    paste0(project_folder, "/results/tables/007_TMB.txt"), 
    sep = "\t"
  )

# cleanup ----------------------------------------------------------------------
rm()
