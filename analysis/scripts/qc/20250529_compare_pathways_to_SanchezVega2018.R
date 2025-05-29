## -----------------------------------------------------------------------------
## Purpose of script: Compare our pathway annotation with Francisco's paper from
##                    2018 (10.1016/j.cell.2018.03.035)
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
# our pathway data
source(paste0(project_dir, "/analysis/scripts/015_pathways.R"))

# francisco's pathways
sanchezvega_pw <- fread("/Users/artzo/Documents/8_asset_files/pathways/20250519_SanchezVega_2018/processed/20250519_pathways.txt")

# wrangle ----------------------------------------------------------------------
compare_pathways <- mutations_with_pathways %>% 
  select(patient_id, gene, gs_name)

# annotate our results with published pathways
idx <- match(compare_pathways$gene, sanchezvega_pw$gene)
compare_pathways$sanchezvega_pw <- sanchezvega_pw$pathway[idx]
compare_pathways$sanchezvega_pw[is.na(compare_pathways$sanchezvega_pw)] <- "Not defined"

# analysis ---------------------------------------------------------------------
# prepare data for alluvial plot - include all genes
alluvial_data <- compare_pathways %>%
  # get unique gene mappings (remove patient duplicates)
  distinct(gene, gs_name, sanchezvega_pw) %>%
  # count genes for each pathway mapping (including "Not defined")
  count(gs_name, sanchezvega_pw, name = "n_genes")

# plot -------------------------------------------------------------------------
library(ggalluvial)

alluvial_data %>%
  # expand data to show individual gene flows
  uncount(n_genes) %>%
  ggplot(aes(axis1 = gs_name, axis2 = sanchezvega_pw)) +
  geom_alluvium(aes(fill = gs_name), alpha = 0.7, width = 0.2) +
  geom_stratum(alpha = 0.8, color = "black", width = 0.2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 3, angle = 0) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Gene Pathway Annotations: Our Analysis vs Sanchez-Vega 2018",
    x = ""
  ) +
  scale_x_discrete(
    limits = c("Our Pathways", "Sanchez-Vega 2018"),
    expand = c(0.1, 0.1)
  ) +
  scale_fill_brewer(type = "qual", palette = "Set2")

# export -----------------------------------------------------------------------
ggsave(paste0(project_dir, "/results/figures/qc/20250529_pathway_comparison_alluvial_complete.pdf"), width = 12, height = 8)

# cleanup ----------------------------------------------------------------------
rm()