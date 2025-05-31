## -----------------------------------------------------------------------------
## Purpose of script: Plot VAFs corrected by max. VAF for each sample. This 
##                    strategy was suggested by Guardant. In patients with 
##                    acquired TMB-H. Guardant samples at progression.
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

p_load(tidyverse, data.table, janitor, ggpubr, ggsci, grid)

# set parameters ---------------------------------------------------------------
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# clonality threshold
subclonal_thresh <- 0.5
ultrasubclonal_thresh <- 0.2

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# annotated progression samples
annotated_mutations_progression <- fread(paste0(project_dir, 
                                                "/data/processed/005_annotated_mutations_progression.txt"))

# patient key
patient_meta_mod <- fread(paste0(project_dir, "/data/processed/003_patient_meta_mod.txt"))

# wrangle ----------------------------------------------------------------------
# remove mutations without VAF data
mut_mod <- annotated_mutations_progression %>% 
  filter(!status_compared_to_baseline_impact == "Not on IMPACT") %>% 
  filter(!is.na(percentage))

# filter for patients with acquired TMB-H
mut_mod <- mut_mod %>% 
  filter(patient_id %in% patient_meta_mod$patient_id[patient_meta_mod$acquired_tmbh == "Y"])

# analysis ---------------------------------------------------------------------
# correct VAF by VAFmax for each sample ----------------------------------------
# correcting by the maximum VAF for each sample effectively corrects for tumor shed
# get VAFmax for each sample
vafmax <- mut_mod %>% 
  group_by(gh_request_id) %>% 
  summarize(vaf_max = max(percentage)) %>% 
  ungroup()

# add vaf_max to df
mut_mod <- mut_mod %>% 
  left_join(vafmax, by = "gh_request_id")

# calculate corrected vaf
mut_mod <- mut_mod %>% 
  mutate(percentage_corr = percentage / vaf_max * 100)

# add mut id
mut_mod <- mut_mod %>% 
  mutate(pt_chrom_pos = paste0(patient_id, "_", chromosome, "_", position))

# plot -------------------------------------------------------------------------
# add number of mutations to legend
# calculate mutation counts by status
status_counts <- mut_mod %>%
  count(status_compared_to_baseline_impact) %>%
  mutate(status_label = paste0(status_compared_to_baseline_impact, " (n=", n, ")"))

# create a named vector for mapping
status_labels <- setNames(status_counts$status_label, status_counts$status_compared_to_baseline_impact)

# plot histogram
p_hist <- mut_mod %>% 
  ggplot(aes(x = percentage_corr, fill = status_compared_to_baseline_impact, color = status_compared_to_baseline_impact)) +
  geom_histogram(alpha = 0.5, position = "dodge", binwidth = 0.1) +
  theme_jco() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.title = element_blank()
  ) +
  scale_fill_jco(labels = status_labels) +
  scale_color_jco(labels = status_labels) +
  scale_x_continuous(breaks = c(0, 10, 25, 50, 75, 100)) +
  labs(
    #title = "VAF distribution",
    x = "VAF / maxVAF (%)",
    y = "Count",
    fill = "Status compared to baseline",
    color = "Status compared to baseline"
  )

# zoom into low VAF region
p_zoom <- p_hist + 
  coord_cartesian(xlim = c(0, 2)) +
  geom_vline(xintercept = c(subclonal_thresh, ultrasubclonal_thresh), linetype = "dashed", alpha = 0.5) +
  theme(legend.position = "none",
        #panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        NULL
        ) +
  scale_x_continuous(breaks = c(0, ultrasubclonal_thresh, subclonal_thresh, 1, 2))

# convert the zoom plot to a grob
zoom_grob <- ggplotGrob(p_zoom)

# create the combined plot with the inset
p_combined <- p_hist +
  annotation_custom(
    grob = zoom_grob,
    xmin = 22,    
    xmax = 100,
    ymin = 5,
    ymax = 65
  )

# display the combined plot
p_combined

# export -----------------------------------------------------------------------
# plot
ggsave(paste0(project_dir, "/results/figures/018_maxVAF_histogram.pdf"), 
       plot = p_combined, width = 5, height = 6)

# data
mut_mod %>% 
  select(gh_request_id, percentage_corr, status_compared_to_baseline_impact, pt_chrom_pos) %>% 
  write.table(paste0(project_dir, "/results/tables/018_mutations_with_corrected_vaf.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

# cleanup ----------------------------------------------------------------------
rm()