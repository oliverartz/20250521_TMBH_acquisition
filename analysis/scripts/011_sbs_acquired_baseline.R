## -----------------------------------------------------------------------------
## Purpose of script: Deconvolve SBS signatures pooled acquired vs baseline 
##                    mutations - Guardant progression samples compared to 
##                    baseline IMPACT for patients with acquired TMB-H.
##
## Author: Oliver Artz
## Date Created: May 28, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, MutationalPatterns, RColorBrewer)

# set parameters ---------------------------------------------------------------
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# reference genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"

# load data --------------------------------------------------------------------
source(paste0(project_dir, "/analysis/scripts/utils/002_plotting_utils.R"))

# signature etiology
sig_eti <- fread(paste0(project_dir, "/data/metadata/20230505_mut_signature_aetiology.csv"))

# patient info
patient_meta_mod <- fread(paste0(project_dir, "/data/processed/003_patient_meta_mod.txt"))

# wrangle ----------------------------------------------------------------------

# analysis ---------------------------------------------------------------------
# get COSMIC signatures
signatures <- get_known_signatures(muttype = "snv",
                                   source = "COSMIC_v3.2",
                                   genome = "GRCh37")

# make mutational matrix
vcf_files <- list.files(path = paste0(project_folder, "/data/processed/hpc/vcf/tmbh_status_compared_to_baseline"),
                        pattern = ".*vcf", 
                        full.names = TRUE)

# filter for patients who acquired tmb-h
sample_types <- c("Y_Acquired",
                  "Y_Baseline")

vcf_files <- vcf_files[basename(vcf_files) %in% paste0(sample_types, ".vcf")]

# define sample names
sample_names <- vcf_files %>% 
  basename() %>% 
  str_remove(".vcf")

# load VCF to GRangesList with SBS only
sbs_grl <- read_vcfs_as_granges(vcf_files, 
                                sample_names, 
                                ref_genome,
                                predefined_dbs_mbs = TRUE,
                                type = "snv")

# generate mutational matrix
mut_mat <- mut_matrix(vcf_list = sbs_grl, ref_genome = ref_genome)

# strict signature fitting by backwards removal --------------------------------
strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)

fit_res_strict <- strict_refit$fit_res

# etiology-based groups --------------------------------------------------------
# extract exposures and add etiology information
signatures_df <- fit_res_strict$contribution %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "SBS_signature") %>% 
  left_join(sig_eti, 
            by = "SBS_signature")

# make long for plotting
signatures_df <- signatures_df %>% 
  pivot_longer(names_to = "status",
               values_to = "exposure",
               colnames(mut_mat))

# calculate relative exposure
signatures_df <- signatures_df %>% 
  group_by(status) %>% 
  mutate(rel_exposure = exposure / sum(exposure)) %>% 
  ungroup()

# change "chemotherapy" for etiology
signatures_df$category[signatures_df$category == "Chemo"] <- "Chemotherapy"

# Cosine similarity ------------------------------------------------------------
# make data frame
cosine_similarities <- data.frame(status = colnames(mut_mat),
                                  cos_similiarity = NA)

# calculate cosine similarity
for (i in seq_along(colnames(mut_mat))) {
  cosine_similarities$cos_similiarity[i] <- round(cos_sim(fit_res_strict$reconstructed[ , i], mut_mat[ , i]), 3)
}

# count exposures from mut_mat
n_exposures <- data.frame(sample = colnames(mut_mat),
                          exposures = colSums(mut_mat)) %>% 
  left_join(cosine_similarities,
            by = c("sample" = "status"))

# add cosine to signatures_df
signatures_df <- signatures_df %>% 
  left_join(cosine_similarities,
            by = "status") %>% 
  mutate(cos_similiarity = round(cos_similiarity, 2))

# remove the N group and refactor
signatures_df <- signatures_df %>% 
  mutate(status = case_when(status == "Y_Acquired" ~ "Acquired",
                            status == "Y_Baseline" ~ "Baseline")) %>% 
  mutate(category = case_when(category == "Chemo" ~ "Chemotherapy",
                              SBS_signature == "SBS3" ~ "HRD",
                              TRUE ~ category))

# plot -------------------------------------------------------------------------
# define colors to match other plots
set3_colors <- brewer.pal(n = 11, name = "Set3")

etiology_colors <- c(
  "APOBEC" = set3_colors[1],
  "Chemotherapy" = set3_colors[2],
  "HRD" = set3_colors[4],
  "MMRd" = set3_colors[5],
  "ROS" = set3_colors[8],
  "UV" = set3_colors[11]
)

signatures_df %>%
  filter(rel_exposure > 0) %>% 
  ggplot(aes(x = status, y = rel_exposure, fill = category)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = SBS_signature), 
            position = position_stack(vjust = 0.5), 
            size = 3, 
            color = "black",
            alpha = 0.8) +
  geom_text(aes(label = cos_similiarity, y = -0.02),
            size = 3,
            alpha = 0.8) +
  theme_jco() +
  theme(panel.grid = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = etiology_colors) +
  labs(
    #title = "Relative exposure of SBS signatures",
    x = "",
    y = "Relative exposure",
    fill = "Etiology")

# export -----------------------------------------------------------------------
# plot
ggsave(paste0(project_folder, "/results/figures/011_sbs_acquired_baseline.pdf"), height = 6, width = 4)

# data
signatures_df %>% 
  fwrite(
    paste0(project_folder, "/results/tables/011_sbs_acquired_baseline.txt"), 
    sep = "\t"
  )

# cleanup ----------------------------------------------------------------------
rm()