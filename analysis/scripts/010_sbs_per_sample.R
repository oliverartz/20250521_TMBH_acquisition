## -----------------------------------------------------------------------------
## Purpose of script: Deconvolve SBS signatures for each IMPACT baseline sample
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

p_load(tidyverse, data.table, MutationalPatterns)

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
# make key for patient and impact samples
# impact_key <- patient_meta_mod %>% 
#   filter(impact_include == "yes") %>% 
#   select(patient_id, impact_baseline, impact_progression) %>% 
#   pivot_longer(col = -patient_id,
#                names_to = "sample_type",
#                values_to = "tumor_sample_barcode")

# analysis ---------------------------------------------------------------------
# get COSMIC signatures
signatures <- get_known_signatures(muttype = "snv",
                                   source = "COSMIC_v3.2",
                                   genome = "GRCh37")

# make mutational matrix
vcf_files <- list.files(path = paste0(project_folder, "/data/processed/hpc/vcf/impact"),
                        pattern = ".*vcf", 
                        full.names = TRUE)

# filter for IMPACT baseline samples
baseline_s<- vcf_files[grepl("impact_baseline", vcf_files)]amples <- patient_meta_mod %>%
  filter(impact_include == "yes", !is.na(impact_baseline)) %>%
  pull(impact_baseline)

# filter for IMPACT baseline samples
vcf_files <- vcf_files[basename(vcf_files) %in% paste0(baseline_samples, ".vcf")]

## QC: Which IMPACT sample is missing?
# df_qc <- patient_meta_mod %>% 
#   select(patient_id, impact_include, impact_baseline) %>% 
#   mutate(vcf_loaded = impact_baseline %in% str_remove(basename(vcf_files), ".vcf")) %>% 
#   filter(vcf_loaded == FALSE)


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

# add patient id
idx <- match(signatures_df$status, patient_meta_mod$impact_baseline)
signatures_df$patient_id <- impact_key$patient_id[idx]

# add acquired_tmbh info
idx <- match(signatures_df$patient_id, patient_meta_mod$patient_id)
signatures_df$acquired_tmbh <- patient_meta_mod$acquired_tmbh[idx]

signatures_df <- signatures_df %>% 
  mutate(acquired_tmbh = case_when(acquired_tmbh == "N" ~ "No",
                                   acquired_tmbh == "Y" ~ "Yes",
                                  TRUE ~ "Other"))

# refactor "chemotherapy" for etiology
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

# add patient id and tmb-h
idx <- match(n_exposures$sample, signatures_df$status)
n_exposures$patient_id <- signatures_df$patient_id[idx]
n_exposures$acquired_tmbh <- signatures_df$acquired_tmbh[idx]

# plot -------------------------------------------------------------------------
# make df for plotting
df_plot <- signatures_df %>% 
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
ggsave(paste0(project_folder, "/results/figures/010_SBS_per_sample.pdf"), height = 6, width = 8)

# data
signatures_df %>% 
  fwrite(
    paste0(project_folder, "/results/tables/010_SBS_per_sample_exposures.txt"), 
    sep = "\t"
  )

cosine_similarities %>% 
  fwrite(
    paste0(project_folder, "/results/tables/010_SBS_per_sample_cosine.txt"), 
    sep = "\t"
  )

# cleanup ----------------------------------------------------------------------
rm()