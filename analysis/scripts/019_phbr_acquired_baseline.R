## -----------------------------------------------------------------------------
## Purpose of script: Plot PHBR in acquired vs baseline mutations in patients
##                    with acquired TMB-H
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

# thresholds to categorize clonality
subclonal_thresh <- 50
ultrasubclonal_thresh <- 20

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# patient_status_compared_to_baseline VAF and PHBR data
phbr_vaf <- fread(paste0(project_dir, "/data/processed/hpc/006_patient_status_compared_to_baseline_phbr_vaf.txt"))

# maxVAF data
vafmax_correct <- fread(paste0(project_folder, "/results/tables/018_mutations_with_corrected_vaf.txt"))

# patient key
patient_meta_mod <- fread(paste0(project_dir, "/data/processed/003_patient_meta_mod.txt"))

# wrangle ----------------------------------------------------------------------
# filter for patients with acquired TMB-H
phbr_vaf_mod <- phbr_vaf %>% 
  filter(patient_id %in% patient_meta_mod$patient_id[patient_meta_mod$acquired_tmbh == "Y"])

# add vafmax data
idx <- match(phbr_vaf_mod$pt_chrom_pos, vafmax_correct$pt_chrom_pos)
phbr_vaf_mod$percentage_corr <- vafmax_correct$percentage_corr[idx]

# categorize maxVAF corrected clonality
phbr_vaf_mod <- phbr_vaf_mod %>% 
  mutate(value = as.numeric(percentage_corr)) %>% 
  mutate(clonality = case_when(
    percentage_corr < ultrasubclonal_thresh ~ "Ultrasubclonal",
    percentage_corr < subclonal_thresh ~ "Subclonal",
    TRUE ~ "Clonal"))

# analysis ---------------------------------------------------------------------

# plot -------------------------------------------------------------------------
# set seed for reproducibility
set.seed(2010)

# make df for plotting
df_plot <- phbr_vaf_mod

# count mutations
n_mut <- df_plot %>% 
  group_by(status) %>% 
  summarize(n = n())

# plot
df_plot %>% 
  ggplot(aes(x = status, y = phbr, fill = status)) +
  geom_hline(yintercept = c(thresh_strong_binder, thresh_weak_binder), 
             linetype = "dashed", alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(shape = 21, alpha = 0.8, width = 0.2, size = 3) +
  geom_text(data = n_mut, aes(label = paste0("n=", n), y = -3)) +
  stat_compare_means(label.x = 1.5) +
  theme_jco() +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  labs(
    #title = "PHBR",
    x = "Status compared to baseline",
    y = "Patient Harmonic-mean Best Rank",
    fill = "Status compared to baseline") +
  scale_fill_jco()

# export -----------------------------------------------------------------------
# save plot
ggsave(paste0(project_dir, "/results/figures/019_phbr_acquired_baseline.pdf"),
       width = 4, height = 6)

# save data
write.table(df_plot, paste0(project_dir, "/results/tables/019_phbr_acquired_baseline_data.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
# cleanup ----------------------------------------------------------------------
rm()