## -----------------------------------------------------------------------------
## Purpose of script: Run all scripts for the analysis
##
## Author: Oliver Artz
## Date Created: May 26, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# start timing -----------------------------------------------------------------
start_time <- Sys.time()
cat("Analysis started at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

# set parameters ---------------------------------------------------------------
project_root <- here::here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/", project_folder)

# load libraries ---------------------------------------------------------------
cat("Loading utilities...\n")
source(paste0(project_dir,"/analysis/scripts/utils/001_preprocess_utils.R"))

# preprocess data --------------------------------------------------------------
cat("Starting preprocessing...\n")
preprocess_start <- Sys.time()

source(paste0(project_dir, "/analysis/scripts/001_preprocess_patient_key.R"))
source(paste0(project_dir, "/analysis/scripts/002_preprocess_guardant_omni.R"))
source(paste0(project_dir, "/analysis/scripts/003_preprocess_patient_meta.R"))
source(paste0(project_dir, "/analysis/scripts/004_preprocess_impact.R"))
source(paste0(project_dir, "/analysis/scripts/005_annotate_acquired_mutations.R"))
source(paste0(project_dir, "/analysis/scripts/006_make_maf.R"))

preprocess_end <- Sys.time()
preprocess_duration <- as.numeric(difftime(preprocess_end, preprocess_start, units = "mins"))
cat("Preprocessing completed in", round(preprocess_duration, 2), "minutes\n")

# analysis ---------------------------------------------------------------------
cat("Starting analysis...\n")
analysis_start <- Sys.time()

source(paste0(project_dir, "/analysis/scripts/007_plot_TMB.R"))
source(paste0(project_dir, "/analysis/scripts/008_plot_MSI.R"))
source(paste0(project_dir, "/analysis/scripts/009_plot_variant_types_fractions.R"))

# input VCF files were generated on the hpc
# Require user input to confirm that HPC analyses are up-to-date
user_input <- readline(prompt = "The below analyses need files processed on the HPC. Are these files up-to-date? (y/n): ")

if (tolower(user_input) != "y" && tolower(user_input) != "yes") {
  cat("Please update the HPC files before continuing.\n")
  cat("Analysis stopped.\n")
  stop("HPC files need to be updated")
}

cat("Continuing with analysis...\n")

source(paste0(project_dir, "/analysis/scripts/010_sbs_per_sample.R"))
source(paste0(project_dir, "/analysis/scripts/011_sbs_acquired_baseline.R"))
source(paste0(project_dir, "/analysis/scripts/012_sbs_baseline_progression_samples.R"))
source(paste0(project_dir, "/analysis/scripts/013_sbs_acquired_per_patient.R"))
source(paste0(project_dir, "/analysis/scripts/014_mean_n_variants_per_patient.R"))
source(paste0(project_dir, "/analysis/scripts/015_pathways.R"))
source(paste0(project_dir, "/analysis/scripts/016_pathways_SanchezVega2018.R"))
source(paste0(project_dir, "/analysis/scripts/017_acquired_alteration_types_progression.R"))
source(paste0(project_dir, "/analysis/scripts/018_maxVAF.R"))
source(paste0(project_dir, "/analysis/scripts/019_phbr_acquired_baseline.R"))
source(paste0(project_dir, "/analysis/scripts/020_phbr_clonality.R"))
source(paste0(project_dir, "/analysis/scripts/021_blood_tissue_tmb.R"))

analysis_end <- Sys.time()
analysis_duration <- as.numeric(difftime(analysis_end, analysis_start, units = "mins"))
cat("Analysis completed in", round(analysis_duration, 2), "minutes\n")

# cleanup ----------------------------------------------------------------------
end_time <- Sys.time()
total_duration <- as.numeric(difftime(end_time, start_time, units = "mins"))

# upload results to shared drive
cat("Uploading results to shared drive...\n")
system(paste0(project_dir, "/analysis/results_to_onedrive.sh"), intern = FALSE)

cat("\n=== TIMING SUMMARY ===\n")
cat("Started at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Ended at:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total runtime:", round(total_duration, 2), "minutes\n")
cat("Preprocessing:", round(preprocess_duration, 2), "minutes\n")
cat("Analysis:", round(analysis_duration, 2), "minutes\n")

rm(start_time, end_time, preprocess_start, preprocess_end, analysis_start, analysis_end,
   total_duration, preprocess_duration, analysis_duration)