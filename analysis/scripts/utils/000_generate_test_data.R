## -----------------------------------------------------------------------------
## Purpose of script: Generate test data for CodeOcean
##
## Author: Oliver Artz
## Date Created: Aug 12, 2025
## Email: artzo@mskcc.org
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, here, janitor)

# set parameters ---------------------------------------------------------------
project_name <- "20250521_TMBH_acquisition"
project_dir <- file.path(here(), project_name)

# define directory or test data
tdir <- file.path(project_dir, "test_data")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)

# generate data ----------------------------------------------------------------
# helper functionn for creating dirs and saving files
wtsv <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); data.table::fwrite(x, path, sep = "\t") }

set.seed(2010)

# Patients
patients <- tibble(
  patient_id = 1:6,
  acquired_tmbh = c("Y","Y","Y","N","N","N"),
  impact_include = "yes",
  impact_baseline = sprintf("P-%07d-T01-IM7", 1000000 + 1:6),
  impact_progression = sprintf("P-%07d-T02-IM7", 1000000 + 1:6),
  guardant_baseline = sprintf("A%07d", 1000000 + 1:6),
  guardant_progression = sprintf("A%07d", 2000000 + 1:6)
)

# data/processed/003_patient_meta_mod.txt
patient_meta_mod <- patients
wtsv(patient_meta_mod, file.path(tdir, "processed/003_patient_meta_mod.txt"))

# data/processed/001_patient_key_mod.txt (minimal columns used downstream)
patient_key_mod <- patients %>%
  transmute(
    patient_id,
    dmp_patient_id = substr(impact_baseline, 1, 9),
    impact_baseline,
    impact_progression,
    guardant_id_baseline = guardant_baseline,
    guardant_id_progression = guardant_progression,
    acquired_high_tmb_y_yes_n_no = ifelse(acquired_tmbh=="Y","Yes","No"),
    baseline_t_tmb_mut_mb_18 = round(rnorm(n(), 6, 2), 2),
    baseline_p_tmb_mut_mb_19 = round(rnorm(n(), 6, 2), 2),
    progression_t_tmb_mut_mb_20 = round(rnorm(n(), 14, 6), 2),
    progression_p_tmb_mut_mb_21 = round(rnorm(n(), 14, 6), 2)
  )

wtsv(patient_key_mod, file.path(tdir, "processed/001_patient_key_mod.txt"))

# data/processed/005_annotated_mutations_progression.txt
genes <- c("KRAS","NRAS","BRAF","TP53","APC","PIK3CA","ERBB2","ATM","MLH1","MSH2")

mk <- function(pid, gh, tmbh, n=150) tibble(
  patient_id = pid,
  gh_request_id = gh,
  acquired_tmbh = tmbh,
  sample_status = "progression",
  analysis_sample_type = "guardant_progression",
  gene = sample(genes, n, TRUE),
  variant_type = sample(c("SNV","Indel"), n, TRUE, prob = c(0.8,0.2)),
  chromosome = as.character(sample(1:22, n, TRUE)),
  position = sample(1e5:9e5, n, TRUE),
  mut_aa = paste0("p.", sample(LETTERS, n, TRUE), sample(10:999, n, TRUE)),
  mut_nt = paste0("c.", sample(100:999, n, TRUE), sample(c("A>G","C>T","G>A","T>C"), n, TRUE)),
  status_compared_to_baseline_impact = sample(c("Acquired","Baseline","Lost"), n, TRUE, prob=c(0.6,0.3,0.1)),
  percentage = round(runif(n, 5, 70), 1)
)
annotated <- map_dfr(seq_len(nrow(patients)), \(i) mk(patients$patient_id[i], patients$guardant_progression[i], patients$acquired_tmbh[i]))

wtsv(annotated, file.path(tdir, "processed/005_annotated_mutations_progression.txt"))

# Minimal VCF writers
write_min_vcf <- function(path, sample, n=120) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  chr <- sample(c(1:22,"X"), n, TRUE)
  pos <- sample(1e5:9e5, n, TRUE)
  ref <- sample(c("A","C","G","T"), n, TRUE)
  
  # Fix: Generate alt bases that are different from corresponding ref bases
  bases <- c("A","C","G","T")
  alt <- character(n)
  for (i in seq_len(n)) {
    alt[i] <- sample(setdiff(bases, ref[i]), 1)
  }
  
  con <- file(path, "w")
  on.exit(close(con))
  writeLines(c(
    "##fileformat=VCFv4.2","##reference=GRCh37",
    paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", sample)
  ), con)
  for (i in seq_len(n)) {
    writeLines(paste(chr[i], pos[i], ".", ref[i], alt[i], 60, "PASS", ".", "GT", "0/1", sep="\t"), con)
  }
}

# data/processed/hpc/vcf/impact/*.vcf (used by 010)
walk(patients$impact_baseline, \(s) write_min_vcf(file.path(tdir, "processed/hpc/vcf/impact", paste0(s,".vcf")), s, 120))

# data/processed/hpc/vcf/gh_request_id/*.vcf (used by 012/013)
walk(patients$guardant_baseline, \(s) write_min_vcf(file.path(tdir, "processed/hpc/vcf/gh_request_id", paste0(s,".vcf")), s, 80))
walk(patients$guardant_progression, \(s) write_min_vcf(file.path(tdir, "processed/hpc/vcf/gh_request_id", paste0(s,".vcf")), s, 140))

# data/processed/hpc/vcf/tmbh_status_compared_to_baseline/Y_*.vcf (used by 011)
write_min_vcf(file.path(tdir, "processed/hpc/vcf/tmbh_status_compared_to_baseline/Y_Acquired.vcf"), "Y_Acquired", 400)
write_min_vcf(file.path(tdir, "processed/hpc/vcf/tmbh_status_compared_to_baseline/Y_Baseline.vcf"), "Y_Baseline", 200)

# data/processed/hpc/006_patient_status_compared_to_baseline_phbr_vaf.txt (used by 019/020)
phbr <- expand_grid(patient_id = patients$patient_id, status = c("Acquired","Baseline")) %>%
  mutate(phbr = round(runif(n(), 0.2, 30), 2), percentage = round(runif(n(), 5, 70), 1))

wtsv(phbr, file.path(tdir, "processed/hpc/006_patient_status_compared_to_baseline_phbr_vaf.txt"))
