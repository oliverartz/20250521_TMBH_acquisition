## -----------------------------------------------------------------------------
## Purpose of script: Plot comparison of PHBR by clonality category of acquired 
##                    mutations in patients with acquired TMB-H
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

p_load(tidyverse, data.table)

# set parameters ---------------------------------------------------------------
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# thresholds to categorize binders
thresh_strong_binder <- 0.5
thresh_weak_binder <- 2

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# PHBR and clonality data
phbr_vaf_mod <- fread(paste0(project_dir, "/results/tables/019_phbr_acquired_baseline_data.txt"))

# wrangle ----------------------------------------------------------------------

# analysis ---------------------------------------------------------------------

# plot -------------------------------------------------------------------------
# set seed for reproducibility
set.seed(2010)

# make df for plotting
df_plot <- phbr_vaf_mod %>% 
  filter(status == "Acquired")

# count mutations
n_mut <- df_plot %>% 
  group_by(clonality) %>% 
  summarize(n = n())

# plot
df_plot %>%
  ggplot(aes(x = clonality, y = phbr, fill = clonality)) +
  geom_hline(
    yintercept = c(thresh_strong_binder, thresh_weak_binder),
    linetype = "dashed", alpha = 0.5
  ) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.3, size = 3, shape = 21, alpha = 0.5) +
  geom_text(data = n_mut, aes(label = paste0("n=", n), y = -3)) +
  stat_compare_means(label.y = 38) +
  theme_jco() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    #title = "PHBR and Clonality",
    x = "Clonality",
    y = "Patient Harmonic-mean Best Rank"
  ) +
  scale_fill_jco()

# export -----------------------------------------------------------------------
# save plot
ggsave(paste0(project_dir, "/results/figures/020_phbr_clonality.pdf"),
       width = 4, height = 6)

# cleanup ----------------------------------------------------------------------
rm()