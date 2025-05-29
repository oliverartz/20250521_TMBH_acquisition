## -----------------------------------------------------------------------------
## Purpose of script: Annotate and plot pathways of acquired mutations. Guardant
##                    samples.
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

p_load(tidyverse, data.table, msigdbr, ggsci, ggpubr)

# set parameters ---------------------------------------------------------------
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# variant types to consider
variant_types_to_consider <- c("SNV", "Indel")

# relevant pathways
relevant_pathways <- c("REACTOME_MAPK_FAMILY_SIGNALING_CASCADES",
                       "REACTOME_PI3K_AKT_ACTIVATION",
                       "REACTOME_DISEASES_OF_MISMATCH_REPAIR_MMR",
                       "REACTOME_DISEASES_OF_DNA_REPAIR")

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# annotated 
annotated_mutations_progression <- fread(paste0(project_dir, "/data/processed/005_annotated_mutations_progression.txt"))

# patient metadata
patient_meta <- fread(paste0(project_dir, "/data/processed/003_patient_meta_mod.txt"))

# wrangle ----------------------------------------------------------------------
# filter for acquired mutations and relevant variant types
acquired_mutations <- annotated_mutations_progression %>%
  filter(!status_compared_to_baseline_impact == "Not on IMPACT") %>%
  filter(variant_type %in% variant_types_to_consider) %>%
  filter(status_compared_to_baseline_impact == "Acquired")

# analysis ---------------------------------------------------------------------
# count gene mutations (filter for recurrent genes)
gene_counts <- acquired_mutations %>%
  group_by(gene) %>%
  summarize(n_mutations = n(), .groups = "drop")

# filter for recurrent genes
gene_counts <- gene_counts %>%
  filter(n_mutations > 2)

# get pathway annotations from MSigDB
pathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

# annotate genes with relevant pathways
gene_pathway_annotations <- pathways %>%
  filter(gene_symbol %in% gene_counts$gene,
         gs_name %in% relevant_pathways) %>%
  select(gene_symbol, gs_name)

# filter genes with multiple annotations
gene_pathway_annotations <- gene_pathway_annotations %>%
  filter(!(gene_symbol == "MSH6" & gs_name == "REACTOME_DISEASES_OF_DNA_REPAIR"),
         !(gene_symbol == "PIK3CA" & gs_name == "REACTOME_MAPK_FAMILY_SIGNALING_CASCADES"))

# create final gene annotations
genes_annotated <- gene_counts %>%
  left_join(gene_pathway_annotations, by = c("gene" = "gene_symbol")) %>%
  mutate(gs_name = case_when(
    is.na(gs_name) ~ "Other",
    gs_name == "REACTOME_MAPK_FAMILY_SIGNALING_CASCADES" ~ "MAPK",
    gs_name == "REACTOME_PI3K_AKT_ACTIVATION" ~ "PI3K", 
    gs_name == "REACTOME_DISEASES_OF_MISMATCH_REPAIR_MMR" ~ "Mismatch Repair",
    gs_name == "REACTOME_DISEASES_OF_DNA_REPAIR" ~ "DNA Repair",
    TRUE ~ gs_name
  ))

# annotate mutations with pathway information
mutations_with_pathways <- acquired_mutations %>%
  left_join(genes_annotated %>% select(gene, gs_name), by = "gene") %>%
  mutate(gs_name = case_when(
    is.na(gs_name) ~ "Other",
    TRUE ~ gs_name
  ))

# count pathway mutations per patient
pathway_counts_per_patient <- mutations_with_pathways %>%
  group_by(patient_id, gs_name) %>%
  summarize(n_mutations = n(), .groups = "drop")

# add patient metadata
pathway_counts_per_patient <- pathway_counts_per_patient %>%
  left_join(patient_meta %>% select(patient_id, acquired_tmbh), by = "patient_id") %>%
  mutate(
    acquired_tmbh = case_when(
      acquired_tmbh == "N" ~ "No",
      acquired_tmbh == "Y" ~ "Yes",
      TRUE ~ "Other"
    ),
    gs_name = factor(gs_name, levels = rev(c("MAPK", "PI3K", "Mismatch Repair", "DNA Repair", "Other")))
  )

# calculate MAPK fraction for all patients
all_patients <- pathway_counts_per_patient %>%
  distinct(patient_id, acquired_tmbh)

total_mutations_per_patient <- pathway_counts_per_patient %>%
  group_by(patient_id) %>%
  summarize(total_mutations = sum(n_mutations), .groups = "drop")

mapk_fractions <- all_patients %>%
  left_join(
    pathway_counts_per_patient %>% 
      filter(gs_name == "MAPK") %>%
      select(patient_id, mapk_mutations = n_mutations),
    by = "patient_id"
  ) %>%
  left_join(total_mutations_per_patient, by = "patient_id") %>%
  mutate(
    mapk_mutations = ifelse(is.na(mapk_mutations), 0, mapk_mutations),
    mapk_percentage = round((mapk_mutations / total_mutations) * 100, 0),
    mapk_percentage = ifelse(is.na(mapk_percentage), 0, mapk_percentage),
    mapk_label = paste0(mapk_percentage, "%")
  )

# calculate total mutations per patient for ordering
patient_order <- pathway_counts_per_patient %>%
  group_by(patient_id, acquired_tmbh) %>%
  summarize(total_mutations = sum(n_mutations), .groups = "drop") %>%
  arrange(acquired_tmbh, total_mutations) %>%
  mutate(patient_id_ordered = factor(patient_id, levels = unique(patient_id)))

# add ordered factor to the main data
pathway_counts_per_patient <- pathway_counts_per_patient %>%
  left_join(patient_order %>% select(patient_id, patient_id_ordered), by = "patient_id")

# also add to mapk_fractions for the text labels
mapk_fractions <- mapk_fractions %>%
  left_join(patient_order %>% select(patient_id, patient_id_ordered), by = "patient_id")

# plot -------------------------------------------------------------------------
# get JCO colors and reverse them manually
jco_colors <- pal_jco()(5)  # Get 5 colors for your 5 pathways
names(jco_colors) <- rev(c("Other", "DNA Repair", "Mismatch Repair", "PI3K", "MAPK"))

pathway_counts_per_patient %>%
  ggplot(aes(x = patient_id_ordered, y = n_mutations, fill = gs_name)) +
  geom_hline(yintercept = c(50, 100), linetype = "dashed", alpha = 0.3) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
    data = mapk_fractions,
    aes(label = mapk_label, y = -2, fill = NULL),
    size = 3) +
  theme_jco() +
  theme(
    legend.position = "right",
    panel.border = element_rect(fill = NA, linewidth = 1)) +
  labs(
    x = "Patient ID",
    y = "Number of Acquired Mutations",
    fill = "Pathway") +
  scale_fill_manual(values = jco_colors) + 
  facet_grid(~ acquired_tmbh, scales = "free_x", space = "free_x")

# export -----------------------------------------------------------------------
# plot
ggsave(paste0(project_folder, "/results/figures/015_pathways.pdf"), height = 6, width = 8)

 # data
pathway_counts_per_patient %>% 
  select(-patient_id_ordered) %>% 
  fwrite(
    paste0(project_folder, "/results/tables/015_pathways.txt"), 
    sep = "\t"
  )

# cleanup ----------------------------------------------------------------------
rm()