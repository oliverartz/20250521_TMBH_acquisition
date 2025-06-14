## -----------------------------------------------------------------------------
## Purpose of script: Plot Acquired and baseline PHBR per patient
##
## Author: Oliver Artz
## Date Created: Jun 11, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, here, janitor, ggsci)

# set parameters ---------------------------------------------------------------
project_name <- "20250521_TMBH_acquisition"
project_dir <- file.path(here(), project_name)

# thresholds to categorize binders
thresh_strong_binder <- 0.5
thresh_weak_binder <- 2

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# PHBR data
phbr_vaf_mod <- fread(paste0(project_dir, "/results/tables/019_phbr_acquired_baseline_data.txt"))

# wrangle ----------------------------------------------------------------------
phbr_vaf_mod <- phbr_vaf_mod %>% 
  mutate(patient_id = as.factor(patient_id))

# analysis ---------------------------------------------------------------------
# calculate per patient median and ranges, CI
patient_summary_stats <- phbr_vaf_mod %>% 
  group_by(patient_id, status) %>% 
  summarise(
    min = round(min(phbr), 2),
    q1 = round(quantile(phbr, 0.25), 2),
    median = round(median(phbr), 2),
    mean = round(mean(phbr), 2),
    q3 = round(quantile(phbr, 0.75), 2),
    max = round(max(phbr), 2),
    sd = round(sd(phbr), 2),
    n = n(),
    .groups = 'drop'
  ) %>% 
  complete(patient_id, status, fill = list(min = NA, q1 = NA, median = NA, 
                                           mean = NA, q3 = NA, max = NA, 
                                           sd = NA, n = 0)) %>% 
  ungroup()
  

# plot -------------------------------------------------------------------------
# make df for plotting
df_plot <- phbr_vaf_mod

# plot
df_plot %>% 
  ggplot(aes(x = status, y = phbr, fill = status)) +
  geom_hline(yintercept = c(thresh_strong_binder, thresh_weak_binder), 
             linetype = "dashed", alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, height = 0, shape = 21, size = 3, alpha = 0.8) +
  stat_compare_means(label = "p.format", label.y = max(df_plot$phbr) * 1.1) +
  theme_jco() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_jco() +
  facet_grid(~ patient_id) +
  labs(
    x = "Status compared to baseline",
    y = "Patient Harmonic-mean Best Rank"
  )
  

# export -----------------------------------------------------------------------
# plot
ggsave(paste0(project_dir, "/results/figures/023_phbr_acquired_baseline_per_patient.pdf"), 
       width = 8, height = 6)

# summary data
patient_summary_stats %>% 
  write_tsv(paste0(project_dir, "/results/tables/023_phbr_acquired_baseline_per_patient_summary_stats.txt"))

# cleanup ----------------------------------------------------------------------
rm()