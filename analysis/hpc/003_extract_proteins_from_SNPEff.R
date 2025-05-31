## -----------------------------------------------------------------------------
## Purpose of script: Extract protein sequences from SNPEff-annotated VCF files
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

p_load(tidyverse, data.table, vcfR, phylotools, here)

# set parameters ---------------------------------------------------------------
project_root <- here()
project_folder <- "20250521_TMBH_acquisition"
project_dir <- paste0(project_root, "/../", project_folder)

processed_data_dir <- paste0(project_dir, "/data/processed/hpc/")
snpeff_output_dir <- paste0(processed_data_dir, "snpeff/patient_status_compared_to_baseline")

# define final output dir for extracted protein seqs
fasta_output_dir <- paste0(processed_data_dir, "fasta/patient_status_compared_to_baseline/")

# load data --------------------------------------------------------------------
# get list of all fasta files in output dir
fasta_files <- list.files(snpeff_output_dir, pattern = "protein_seq.fasta", full.names = TRUE)

# read all fasta files into list
fasta <- lapply(fasta_files, read.fasta)

# get list of all vcf files in output dir
vcf_files <- list.files(snpeff_output_dir, pattern = "_annotated.vcf.tsv", full.names = TRUE)

## read all vcf files
vcf <- lapply(vcf_files, fread)

# analysis ---------------------------------------------------------------------
# make output dir
dir.create(fasta_output_dir, recursive = TRUE)

## define function to process each fasta and vcf file per sample ---------------
process_data <- function(fasta_file_path) {
  
  # get sample name from fasta file
  sample_name <- basename(fasta_file_path) %>% str_remove("_protein_seq.fasta")
  
  # get fasta
  fasta_file <- read.fasta(fasta_file_path)
  
  # split columns by white space
  fasta_df <- fasta_file %>%
    separate(seq.name, c("gene", "type", "pos", "Ref", "Alt", "HGVS.p"), sep = " ") %>% 
    dplyr::rename(sequence = "seq.text") %>% 
    suppressWarnings()
  
  # filter for mutated proteins
  fasta_df_mut <- fasta_df %>% 
    filter(type == "Variant") %>% 
    mutate(HGVS.p = gsub("HGVS.p:", "", HGVS.p))
  
  # make CHROM column
  fasta_df_mut$CHROM <- sapply(strsplit(fasta_df_mut$pos, ":"), function(x) x[1])
  
  # make START and END columns 
  fasta_df_mut$START <- sapply(strsplit(fasta_df_mut$pos, ":"), function(x) x[2])
  fasta_df_mut$START <- sapply(strsplit(fasta_df_mut$START, "-"), function(x) x[1])
  fasta_df_mut$END <- sapply(strsplit(fasta_df_mut$pos, "-"), function(x) x[2])
  
  # make mut_id column
  fasta_df_mut$mut_id <- paste(
    fasta_df_mut$gene,
    fasta_df_mut$HGVS.p,
    sep = "_"
  )
  
  # load corresponding VCF file
  vcf_file <- fread(paste0(snpeff_output_dir, "/", sample_name, "_annotated.vcf.tsv"))
  
  # add sample name column
  vcf <- vcf_file %>% mutate(sample = sample_name)
  
  # make mut_id column
  vcf$mut_id <- paste(
    vcf$`ANN[0].FEATUREID`,
    vcf$`ANN[0].HGVS_P`,
    sep = "_")
  
  # filter for mutations that have an impact on protein sequences
  vcf <- vcf %>% 
    filter(!`ANN[0].HGVS_P` == "")
  
  # There are mutations that are next to each other leading to the same HGVS_P
  # Since these mutations have the same outcome in terms of protein sequence,
  # I will filter the VCF file to only contain one of these mutations.
  
  vcf <- vcf %>%
    distinct(mut_id, .keep_all = TRUE)
  
  # add protein sequence to vcf file
  idx <- match(vcf$mut_id, fasta_df_mut$mut_id)
  vcf$sequence <- fasta_df_mut$sequence[idx]
  
  # extract position of mutation within peptide
  vcf <- vcf %>% 
    mutate(aa_position = str_extract(`ANN[0].HGVS_P`, "\\d+") %>% as.numeric())
  
  # find frameshift mutations
  vcf <- vcf %>% 
    mutate(frameshift = str_detect(`ANN[0].HGVS_P`, "fs"))
  
  # missense: restrict peptide size to window around the mutation
  # frameshift: take whole sequence
  vcf <- vcf %>% 
    mutate(mhcpan_input = case_when(
      frameshift == FALSE ~ ifelse(
        aa_position <= 7,
        str_sub(sequence, 1, aa_position + 7),
        str_sub(sequence, aa_position - 7, aa_position + 7)),
      frameshift == TRUE ~ str_sub(sequence, aa_position, -1)))
  
  # remove synonymous mutations
  # synonymous mutations have NA as sequence
  vcf <- vcf %>% 
    filter(!is.na(sequence))
  
  # remove nonsense mutations since they don't produce new AAs
  vcf <- vcf %>% 
    filter(!str_sub(`ANN[0].HGVS_P`, -1, -1) == "*")
  
  # remove everything after STOP codon (*)
  vcf <- vcf %>%
    mutate(mhcpan_input = str_remove(mhcpan_input, "\\*.*"))
  
  # determine length of peptide
  vcf <- vcf %>% 
    mutate(peptide_length = str_length(mhcpan_input))
  
  # if frameshift is very early in the sequence, get the sequence to be 15 AA long or start at the first AA
  vcf <- vcf %>%
    mutate(mhcpan_input = ifelse(
      peptide_length < 15 & frameshift == TRUE,
      ifelse(
        aa_position <= 7,
        str_sub(sequence, 1, 15),
        str_sub(sequence, aa_position - 7, -1)),
      mhcpan_input)) %>% 
    mutate(mhcpan_input = str_remove(mhcpan_input, "\\*.*")) %>% 
    mutate(peptide_length = str_length(mhcpan_input))
  
  # if frameshift produces a prompt stop codon so that sequence is too short, extend sequence upstream to 15 AA
  # df_temp <- vcf %>% 
  #   mutate(mhcpan_input = ifelse(
  #     peptide_length < 15 & frameshift == TRUE,
  #     str_sub(sequence, aa_position - (14 + peptide_length), aa_position - peptide_length),
  #     mhcpan_input)) %>% 
  #   mutate(new_peptide_length = nchar(mhcpan_input))
  
  
  # netMHCpan cuts peptide IDs after 15 characters
  # make peptide_ID based on row_number
  vcf <- vcf %>% 
    mutate(peptide_ID = row_number())
  
  # make input for netMHCpan
  mhc_pan_input <- vcf %>%
    na.omit() %>% 
    filter(!mhcpan_input == "")
  
  # export
  write_csv(mhc_pan_input, paste0(fasta_output_dir, sample_name, "_mhcpan_input.csv"))
  
  ## netMHCpan takes a fasta file as input
  # define the path for the output FASTA file
  output_file_path <- paste0(fasta_output_dir, sample_name, ".fasta")
  
  # start writing the file
  file_conn <- file(output_file_path, "w")
  
  # loop over every line in data frame
  apply(mhc_pan_input, 1, function(row) {
    writeLines(paste0(">", row['peptide_ID']), file_conn)
    writeLines(row['mhcpan_input'], file_conn)
  })
  
  # end writing the file
  close(file_conn)
}

# call function ----------------------------------------------------------------
# process all fasta files
lapply(fasta_files, process_data)

# cleanup ----------------------------------------------------------------------
rm()