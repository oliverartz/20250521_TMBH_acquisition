## -----------------------------------------------------------------------------
## Purpose of script: Annotate and plot pathways of acquired mutations. Guardant
##                    samples. Using Sanchez-Vega (2018) assignments.
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

# variant types to consider
variant_types_to_consider <- c("SNV", "Indel")

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# annotated 
annotated_mutations_progression <- fread(paste0(project_dir, "/data/processed/005_annotated_mutations_progression.txt"))

# patient metadata
patient_meta <- fread(paste0(project_dir, "/data/processed/003_patient_meta_mod.txt"))

# francisco's pathways
sanchezvega_pw <- fread("/Users/artzo/Documents/8_asset_files/pathways/20250519_SanchezVega_2018/processed/20250519_pathways.txt")

# wrangle ----------------------------------------------------------------------
# filter for acquired mutations
acquired_mutations <- annotated_mutations_progression %>%
  filter(!status_compared_to_baseline_impact == "Not on IMPACT") %>%
  filter(variant_type %in% variant_types_to_consider) %>% 
  filter(status_compared_to_baseline_impact == "Acquired")

# analysis ---------------------------------------------------------------------
# annotate mutations
idx <- match(acquired_mutations$gene, sanchezvega_pw$gene)
acquired_mutations$pathway <- sanchezvega_pw$pathway[idx]
acquired_mutations$pathway[is.na(acquired_mutations$pathway)] <- "Other"

# count pathways per patient
pathway_counts_per_patient <- acquired_mutations %>%
  group_by(patient_id, pathway) %>%
  summarize(n_mutations = n(), .groups = "drop")

# add patient metadata
idx <- match(pathway_counts_per_patient$patient_id, patient_meta$patient_id)
pathway_counts_per_patient$acquired_tmbh <- patient_meta$acquired_tmbh[idx]

# refactor acquired_tmbh and reorder pathways (reverse y-axis order)
pathway_counts_per_patient <- pathway_counts_per_patient %>% 
  mutate(
    acquired_tmbh = case_when(
      acquired_tmbh == "Y" ~ "Yes",
      acquired_tmbh == "N" ~ "No",
      TRUE ~ "Other"
    ),
    # Reverse the pathway order: "Other" first, then reverse alphabetical
    pathway = factor(pathway, 
                     levels = rev(c(sort(setdiff(unique(pathway), "Other")), "Other")))
  )

# calculate RTK-RAS fractions
rtkras_fractions <- pathway_counts_per_patient %>%
  group_by(patient_id, acquired_tmbh) %>%
  summarize(
    total_mutations = sum(n_mutations),
    rtkras_mutations = sum(n_mutations[pathway == "RTK-RAS"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rtkras_fraction = rtkras_mutations / total_mutations,
    rtkras_percentage = round(rtkras_fraction * 100, 0),
    rtkras_label = paste0(rtkras_percentage, "%")
  )

# create patient ordering based on RTK-RAS fraction
patient_order <- rtkras_fractions %>%
  arrange(acquired_tmbh, desc(rtkras_fraction)) %>%  # Order by RTK-RAS fraction (high to low)
  mutate(patient_id_ordered = factor(patient_id, levels = unique(patient_id)))

# add ordered factor to both datasets
pathway_counts_per_patient <- pathway_counts_per_patient %>%
  left_join(patient_order %>% select(patient_id, patient_id_ordered), by = "patient_id")

rtkras_fractions <- rtkras_fractions %>%
  left_join(patient_order %>% select(patient_id, patient_id_ordered), by = "patient_id")

# plot -------------------------------------------------------------------------
create_pathway_plot <- function(color_scheme = "viridis") {
  
  # Define color schemes
  if (color_scheme == "viridis") {
    tile_fill <- scale_fill_viridis_c(name = "# Mutations", trans = "sqrt",
                                      option = "plasma", na.value = "grey95")
    rtkras_fill <- scale_fill_viridis_c(name = "RTK-RAS\nFraction",
                                        option = "viridis", na.value = "grey95",
                                        labels = scales::percent)
    text_color <- "white"
    tile_color <- "white"
  } else if (color_scheme == "white_red") {
    tile_fill <- scale_fill_gradient(low = "white", high = "red", name = "# Mutations", 
                                     na.value = "grey95", trans = "sqrt")
    rtkras_fill <- scale_fill_gradient(low = "white", high = "red", name = "RTK-RAS\nFraction", 
                                       na.value = "grey95", labels = scales::percent)
    text_color <- "black"
    tile_color <- "black"
  }
  
  # tile plot
  p_tile <- pathway_counts_per_patient %>%
    ggplot(aes(x = patient_id_ordered, y = pathway, fill = n_mutations)) +
    geom_tile(color = tile_color, size = 0.5) +
    geom_text(aes(label = n_mutations), color = text_color, size = 3, fontface = "bold") +
    theme_jco() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "lightgrey"),
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    tile_fill +
    labs(x = NULL, y = "Pathway") +
    facet_grid(~ acquired_tmbh, scales = "free_x", space = "free_x")
  
  # RTK-RAS fraction plot
  p_rtkras <- rtkras_fractions %>%
    ggplot(aes(x = patient_id_ordered, y = "RTK-RAS\nFraction", fill = rtkras_fraction)) +
    geom_tile(color = tile_color, size = 0.5) +
    geom_text(aes(label = rtkras_label), color = text_color, size = 3, fontface = "bold") +
    theme_jco() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, linewidth = 1),
      legend.position = "none",
      strip.text = element_blank(),
      strip.background = element_blank(),
      plot.margin = margin(0, 5.5, 5.5, 5.5, "pt")
    ) +
    rtkras_fill +
    labs(x = "Patient ID", y = "") +
    facet_grid(~ acquired_tmbh, scales = "free_x", space = "free_x")
  
  # combine plots
  cowplot::plot_grid(p_tile, p_rtkras,
                     nrow = 2,
                     rel_heights = c(4.5, 1),
                     align = "v",
                     axis = "lr")
}

# export plots -----------------------------------------------------------------
# original viridis version
p_viridis <- create_pathway_plot("viridis")
ggsave(
  plot = p_viridis, 
  filename = paste0(project_folder, "/results/figures/016_pathways_SanchezVega2018.pdf"), 
  height = 6, width = 8)

# white-to-red version
p_white_red <- create_pathway_plot("white_red")
ggsave(
  plot = p_white_red, 
  filename = paste0(project_folder, "/results/figures/016_pathways_SanchezVega2018_white_red.pdf"), 
  height = 6, width = 8)

# export data ------------------------------------------------------------------
# pathway annotations for each mutation
df_export <- acquired_mutations %>% 
  select(patient_id,
         acquired_tmbh,
         variant_type,
         indel_type,
         gene,
         chromosome,
         position,
         mut_aa,
         mut_nt,
         status_compared_to_baseline_impact,
         pathway)

df_export %>% 
  fwrite(
    paste0(project_folder, "/results/tables/016_pathways_SanchezVega2018.txt"), 
    sep = "\t"
  )

 # cleanup ---------------------------------------------------------------------
rm()