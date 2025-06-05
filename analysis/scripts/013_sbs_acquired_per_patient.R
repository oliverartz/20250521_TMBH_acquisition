## -----------------------------------------------------------------------------
## Purpose of script: Deconvolve SBS signatures for acquired mutations in patients
##                    with and without acquired TMB-H. 
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

p_load(tidyverse, data.table, MutationalPatterns, RColorBrewer, ggpubr, ggsci)

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
# make key for guardant samples for annotation
guardant_sample_key <- patient_meta_mod %>% 
  select(patient_id, acquired_tmbh, guardant_baseline, guardant_progression) %>% 
  pivot_longer(col = -c(patient_id, acquired_tmbh),
               names_to = "sample_type",
               values_to = "gh_request_id")

# analysis ---------------------------------------------------------------------
# get COSMIC signatures
signatures <- get_known_signatures(muttype = "snv",
                                   source = "COSMIC_v3.2",
                                   genome = "GRCh37")

# make mutational matrix
vcf_files <- list.files(path = paste0(project_folder, "/data/processed/hpc/vcf/patient_status_compared_to_baseline"),
                        pattern = ".*vcf", 
                        full.names = TRUE)


# filter for acquired mutations only
vcf_files <- vcf_files[grepl("_Acquired\\.vcf$", basename(vcf_files))]

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

# refactor category
signatures_df <- signatures_df %>% 
  mutate(category = case_when(category == "Chemo" ~ "Chemotherapy",
                              SBS_signature == "SBS3" ~ "HRD",
                              TRUE ~ category))

# annotate samples with sample id and sample type
signatures_df$patient_id <- signatures_df$status %>% str_remove("_Acquired")
idx <- match(signatures_df$patient_id, guardant_sample_key$patient_id)
signatures_df$acquired_tmbh <- guardant_sample_key$acquired_tmbh[idx]

signatures_df <- signatures_df %>% 
  mutate(acquired_tmbh = case_when(acquired_tmbh == "N" ~ "No",
                                   acquired_tmbh == "Y" ~ "Yes",
                                   TRUE ~ "Other"))

# plot -------------------------------------------------------------------------
# make df for plotting
df_plot <- signatures_df %>% 
  filter(rel_exposure > 0)

# add TMZ to chemotherapy 
df_plot <- df_plot %>% 
  mutate(category = case_when(SBS_signature == "SBS11" ~ "Chemotherapy",
                              category == "DNA repair (not MMRd)" ~ "Other DNA repair",
                              TRUE ~ category))

# sort patient_id by cosine similarity within each group
df_plot <- df_plot %>%
  arrange(acquired_tmbh, desc(cos_similiarity)) %>%
  mutate(patient_id_ordered = factor(patient_id, levels = unique(patient_id)))

# define colors to match 012_sbs_signatures_samples
set3_colors <- brewer.pal(n = 11, name = "Set3")

etiology_colors <- c(
  "AID" = "#FBB4AE",
  "APOBEC" = set3_colors[1],
  "Chemotherapy" = set3_colors[2],
  "Clock" = set3_colors[3],
  "Other DNA repair" = "#E5D8BD",
  "HRD" = set3_colors[4],
  "MMRd" = set3_colors[5],
  "Platinum" = set3_colors[6],
  "POLE" = set3_colors[7],
  "ROS" = set3_colors[8],
  "Tobacco" = set3_colors[9],
  "Unknown" = set3_colors[10],
  "UV" = set3_colors[11])

# bar plot
p_bar <- df_plot %>%
  ggplot(aes(x = patient_id_ordered, y = rel_exposure, fill = category)) +
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
        panel.border = element_rect(fill = NA, linewidth = 1)) +
  scale_fill_manual(values = etiology_colors) +
  labs(
    #title = "Relative exposure of SBS signatures",
    x = "Patient ID",
    y = "Relative Exposure",
    fill = "Etiology") +
  facet_grid(~ acquired_tmbh, scales = "free_x", space = "free_x")

# box plot
df_plot <- signatures_df %>% 
  group_by(category) %>% 
  mutate(sum = sum(rel_exposure)) %>% 
  ungroup() %>% 
  filter(sum > 0)

p_box <- df_plot %>% 
  ggplot(aes(x = acquired_tmbh, y = rel_exposure, fill = acquired_tmbh)) +
  stat_compare_means(label = "p.format",
                     label.y = 0.55,
                     size = 3) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 3, shape = 21, alpha = 0.8) +
  facet_wrap(~ category) +
  ylim(0, 0.65) +
  theme_jco() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none") +
  scale_fill_jco() +
  labs(
    x = "Acquired High TMB",
    y = "Relative Exposure",
    fill = "Acquired High TMB")

# export -----------------------------------------------------------------------
# plot
ggsave(
  plot = p_bar,
  paste0(project_folder, "/results/figures/013_sbs_acquired_per_patient_barplot.pdf"), height = 6, width = 12)

ggsave(
  plot = p_box,
  paste0(project_folder, "/results/figures/013_sbs_acquired_per_patient_boxplot.pdf"), height = 6, width = 8)

# data
df_plot %>% 
  fwrite(
    paste0(project_folder, "/results/tables/013_sbs_acquired_per_patient.txt"), 
    sep = "\t"
  )

# cleanup ----------------------------------------------------------------------
rm()