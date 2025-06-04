## -----------------------------------------------------------------------------
## Purpose of script: Plot SBS per sample without etiologies that have no 
##                    exposures
##
## Author: Oliver Artz
## Date Created: Jun 4, 2025
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
project_folder <- "20250521_TMBH_acquisition"
project_dir <- file.path(here::here(), project_folder)

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# SBS signatures per sample
signatures_df <- fread(file.path(project_dir, "/results/tables/010_SBS_per_sample_exposures.txt"))

# wrangle ----------------------------------------------------------------------
# find signatures with no exposures
signatures_without_exp <- signatures_df %>% 
  group_by(category) %>% 
  summarize(total_exp = sum(rel_exposure)) %>% 
  filter(total_exp <= 0) %>% 
  pull(category)

# remove signatures with no exposures
sigs_filtered <- signatures_df %>% 
  filter(!category %in% signatures_without_exp)

# analysis ---------------------------------------------------------------------

# plot -------------------------------------------------------------------------
# make df for plotting
df_plot <- sigs_filtered %>% 
  group_by(status, category, acquired_tmbh) %>% 
  summarize(avg_cat_exp = mean(rel_exposure)) %>% 
  ungroup()

# plot
df_plot %>% 
  ggplot(aes(x = acquired_tmbh, y = avg_cat_exp, fill = acquired_tmbh)) +
  stat_compare_means(label = "p.format",
                     label.y = 0.17,
                     size = 3) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, shape = 21, alpha = 0.8) +
  facet_wrap(~ category) +
  theme_jco() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none") +
  scale_fill_jco() +
  labs(
    #title = "Mutational Signatures",
    #subtitle = "Average relative exposure per patient",
    x = "Acquired High TMB",
    y = "Relative Exposure",
    fill = "Acquired High TMB")

# export -----------------------------------------------------------------------
# plot
ggsave(paste0(project_dir, "/results/figures/022_sbs_per_sample_etiologies_with_exp.pdf"), height = 6, width = 8)

# cleanup ----------------------------------------------------------------------
rm()