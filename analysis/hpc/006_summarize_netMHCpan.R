## -----------------------------------------------------------------------------
## Purpose of script: Summarize netMHCpan output
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

p_load(tidyverse, data.table, here, janitor, DescTools)

# set parameters ---------------------------------------------------------------
project_root <- here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/../", project_folder)

# directory that has all the HPC processed data
output_dir <- paste0(project_dir, "/data/processed/hpc/")

# comparison group name
group_name <- "patient_status_compared_to_baseline"

# load data --------------------------------------------------------------------
# annotated mutations of progression samples
annotated_mutations_mod <- fread(paste0(project_dir, "/data/processed/005_annotated_mutations_progression.txt"))

# netMHCpan output -------------------------------------------------------------
# get file names
files <- list.files(path = paste0(output_dir, "/netMHCpan/", group_name, "/"),
                    pattern = ".xls",
                    full.names = TRUE)

# filter out empty files
files <- files[sapply(files, function(file) file.info(file)$size > 0)]

# load files
netMHCpan_output <- lapply(files, function(x) {
  
  haplotype_x <- fread(x) %>% dplyr::select(V4) %>% slice_head(n = 1) %>% pull()
  
  sample_x <- x %>% basename() %>% 
    str_remove(haplotype_x) %>% 
    str_remove("_.xls")
  
  # Read the file with proper handling of headers
  df <- fread(x, skip = 1) # Skip header and first description line
  
  # Check if columns exist before converting
  if("Pos" %in% names(df)) df$Pos <- as.integer(df$Pos)
  if("ID" %in% names(df)) df$ID <- as.character(df$ID)
  if("EL-score" %in% names(df)) df$`EL-score` <- as.numeric(df$`EL-score`)
  if("EL_Rank" %in% names(df)) df$EL_Rank <- as.numeric(df$EL_Rank) 
  if("BA-score" %in% names(df)) df$`BA-score` <- as.numeric(df$`BA-score`)
  if("BA_Rank" %in% names(df)) df$BA_Rank <- as.numeric(df$BA_Rank)
  if("Ave" %in% names(df)) df$Ave <- as.numeric(df$Ave)
  if("NB" %in% names(df)) df$NB <- as.numeric(df$NB)
  
  # Add sample and haplotype info
  df$haplotype <- haplotype_x
  df$sample <- sample_x
  
  return(df)
}) %>% 
  bind_rows()

# mutation identifier ----------------------------------------------------------
files <- list.files(path = paste0(output_dir, "/fasta/", group_name, "/"),
                    pattern = ".csv", full.names = TRUE)


mutation_id <- lapply(files, function(x){
  fread(x) %>% 
    mutate(CHROM = as.character(CHROM),
           pt_mut_id = paste0(sample, "_", peptide_ID))
}) %>% bind_rows()

# wrangle ----------------------------------------------------------------------
# fix colnames
netMHCpan_output_mod <- clean_names(netMHCpan_output)

# remove NA reads
netMHCpan_output_mod <- netMHCpan_output_mod %>% 
  filter(!is.na(el_rank))

# calculate PHBR for each mutation and patient ---------------------------------
# calculate best rank for each allele
best_EL_rank <- netMHCpan_output_mod %>% 
  group_by(id, haplotype, sample) %>% 
  summarize(best_el_rank = min(el_rank)) %>% 
  ungroup()

# calculate PHBR
phbr <- best_EL_rank %>% 
  group_by(id, sample) %>% 
  summarize(phbr = Hmean(best_el_rank, na.rm = TRUE)) %>% 
  ungroup()

# separate status and patient
phbr <- phbr %>% 
  separate(sample, into = c("patient_id", "status"))

# add peptide information to phbr
phbr <- phbr %>% 
  mutate(pt_mut_id = paste0(patient_id, "_", status, "_", id))

# combine PHBR and VAF data ----------------------------------------------------
phbr_vaf <- phbr

phbr_vaf <- phbr_vaf %>% 
  left_join(mutation_id %>% 
              dplyr::select(CHROM, POS, pt_mut_id),
            by = "pt_mut_id") %>% 
  mutate(pt_chrom_pos = paste0(patient_id, "_",
                               CHROM, "_",
                               POS))

annotated_mutations_mod <- annotated_mutations_mod %>% 
  mutate(pt_chrom_pos = paste0(patient_id, "_", 
                               chromosome, "_",
                               position))

phbr_vaf <- phbr_vaf %>% 
  left_join(annotated_mutations_mod %>% 
              dplyr::select(pt_chrom_pos, 
                            percentage),
            by = "pt_chrom_pos") %>% 
  filter(!is.na(percentage))

# QC ---------------------------------------------------------------------------
# What are the numbers for acquired and baseline?
status_counts <- phbr_vaf$status %>% tabyl()
print(status_counts)
cat("Total mutations analyzed:", nrow(phbr_vaf), "\n")
cat("Baseline mutations:", status_counts$n[status_counts$.== "Baseline"], "\n")
cat("Acquired mutations:", status_counts$n[status_counts$.== "Acquired"], "\n")

# Do we have the right patients?
unique_patients <- phbr_vaf$patient_id %>% 
  unique() %>% 
  as.numeric() %>% 
  sort()

cat("\nPatients with PHBR data:", paste(unique_patients, collapse = ", "), "\n")
cat("Number of patients with PHBR data:", length(unique_patients), "\n")

# export -----------------------------------------------------------------------
write.table(
  phbr_vaf,
  file = paste0(output_dir, "006_patient_status_compared_to_baseline_phbr_vaf.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

# cleanup ----------------------------------------------------------------------
rm()