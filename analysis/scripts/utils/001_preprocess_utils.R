## -----------------------------------------------------------------------------
## Purpose of script: Define functions used in data preprocessing steps
##
## Author: Oliver Artz
## Date Created: May 23, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# annotate_patient_mutations function ------------------------------------------
# Purpose:
# This function annotates Guardant360 data with mutation information from the IMPACT and 
# OMNI gene panels. It processes mutations for each patient in a given list of patient IDs, 
# determining whether mutations are Baseline, Acquired, or Lost based on their status in 
# the baseline IMPACT panel and OMNI panel.
#
# Parameters:
# - guardant_data: data frame containing mutation data from the Guardant360 assay.
#                  Must include columns such as patient_id, gene, and mut_id_short.
# - impact_data: data frame containing mutation data from the IMPACT gene panel.
#                Must include columns such as patient_id, hugo_symbol, mut_id_short, and more.
# - gene_panel_list: data frame mapping genes to specific gene panels (e.g., IMPACT, OMNI).
#                    Must include columns panel, gene, and gene_symbol.
# - patient_ids: vector of patient IDs to process and annotate.
#
# Returns:
# A combined data frame with annotated mutations for all patients. The output includes 
# columns indicating gene presence in the panels, mutation presence, and status relative to 
# the baseline IMPACT panel (e.g., Baseline, Acquired, Lost, Not on IMPACT, Not on OMNI).


annotate_patient_mutations <- function(guardant_data, impact_data, gene_panel_list, patient_ids) {
  
  # define a helper function for a single patient
  process_single_patient <- function(patient_id_of_interest) {
    
    # DEBUG ----
    if (FALSE) {
      guardant_data <- guardant_mod
      impact_data <- impact_mod
      gene_panel_list <- panel_gene_list
      patient_ids <- unique(guardant_mod$patient_id)
      patient_id_of_interest <- 14
    }
    
    #  END DEBUG ----
    
    # filter for the current patient in guardant_data
    guardant_patient_data <- guardant_data %>% filter(patient_id == patient_id_of_interest)
    
    # annotate gene presence in impact panel
    # important: this also labels mutations in baseline guardant samples as 
    # acquired if they are not present in the baseline IMPACT sample
    # this needs to be filtered later when making MAF and VCF files!
    guardant_patient_data <- guardant_patient_data %>%
      rowwise() %>%
      mutate(
        gene_in_impact = case_when(
          gene %in% gene_panel_list$gene[gene_panel_list$panel == baseline_impact_version] ~ "present_on_impact_ver",
          TRUE ~ "absent_on_impact_ver"),
        mut_in_impact_baseline = case_when(
          mut_id_short %in% impact_data$mut_id_short[impact_data$patient_id == patient_id] ~ "mut_on_impact",
          TRUE ~ "wt_on_impact")) %>%
      ungroup() %>%
      mutate(
        status_compared_to_baseline_impact = case_when(
          gene_in_impact == "present_on_impact_ver" & mut_in_impact_baseline == "mut_on_impact" ~ "Baseline",
          gene_in_impact == "present_on_impact_ver" & mut_in_impact_baseline == "wt_on_impact" ~ "Acquired",
          gene_in_impact == "absent_on_impact_ver" ~ "Not on IMPACT"),
        gene_in_omni = "present_on_omni")
    
    # process lost mutations based on impact_data
    impact_patient_data <- impact_data %>% filter(patient_id == patient_id)
    
    impact_patient_data <- impact_patient_data %>%
      rowwise() %>%
      mutate(
        gene_in_omni = case_when(
          hugo_symbol %in% gene_panel_list$gene_symbol[gene_panel_list$panel == "OMNI"] ~ "present_on_omni",
          TRUE ~ "absent_on_omni"),
        mut_in_omni = case_when(
          mut_id_short %in% guardant_patient_data$mut_id_short[guardant_patient_data$patient_id == patient_id] ~ "mut_on_omni",
          TRUE ~ "wt_on_omni"),
        status_compared_to_baseline_impact = case_when(
          gene_in_omni == "present_on_omni" & mut_in_omni == "mut_on_omni" ~ "Baseline",
          gene_in_omni == "present_on_omni" & mut_in_omni == "wt_on_omni" ~ "Lost",
          gene_in_omni == "absent_on_omni" ~ "Not on OMNI")) %>%
      ungroup() %>%
      filter(status_compared_to_baseline_impact == "Lost") %>%
      #dplyr::select(-gene) %>%
      dplyr::rename(
        gene = hugo_symbol,
        position = start_position,
        mut_aa = hgv_sp_short) %>% 
      mutate(mut_nt = paste0(reference_allele,">", tumor_seq_allele2))
    
    # fill missing columns and reformat
    missing_columns <- setdiff(colnames(guardant_patient_data), colnames(impact_patient_data))
    
    impact_patient_data <- impact_patient_data %>%
      dplyr::select(patient_id, gene, variant_type, chromosome, position, 
                    mut_aa, mut_id_short, status_compared_to_baseline_impact, 
                    gene_in_omni, mut_nt) %>%
      bind_cols(as_tibble(matrix(NA, nrow = nrow(.), 
                                 ncol = length(missing_columns), 
                                 dimnames = list(NULL, missing_columns)))) %>%
      mutate(
        sample_status = "IMPACT_baseline",
        analysis_sample_type = "IMPACT_baseline",
        gene_in_impact = "present_on_impact_ver",
        indel_type = case_when(
          variant_type == "INS" ~ "Insertion",
          variant_type == "DEL" ~ "Deletion",
          TRUE ~ variant_type),
        variant_type = case_when(
          variant_type == "SNP" ~ "SNV",
          variant_type %in% c("INS", "DEL") ~ "Indel",
          TRUE ~ variant_type
        ))
    
    # reorder column
    impact_patient_data <- impact_patient_data[, colnames(guardant_patient_data)]
    
    # combine processed data
    rbind(guardant_patient_data, impact_patient_data)
  }
  
  # apply function to all patients and combine results
  all_patient_results <- lapply(patient_ids, process_single_patient)
  do.call(rbind, all_patient_results)
}

# make MAF files for different grouping strategies -----------------------------
make_maf <- function(sample, sub_dir, output_dir){
  
  # DEBUG ----
  if (FALSE) {
    sub_dir <- "gh_request_id"
    sub_dir <- "tmbh_status_compared_to_baseline"
    output_dir <- paste0(project_dir, "/data/processed/006_maf/")
    sample <- "A0652291"
    sample <- "Y_Baseline"
  }
  
  # define columns to keep
  columns_to_keep <- c("chromosome", "position", "reference_allele", "tumor_seq_allele2", "tumor_sample_barcode")
  
  # define sample id
  sample_id <- sample
  
  # filter for SNV and Indels
  test_maf <- annotated_mutations_progression %>% 
    filter(variant_type %in% c("SNV", "Indel"))
  
  # filter for single sample 
  test_maf <- test_maf[get(sub_dir) == sample_id]
  
  # filter for columns to keep
  test_maf <- test_maf %>% 
    select(all_of(columns_to_keep))
  
  # rename position column to fit VCF convention
  test_maf <- test_maf %>% 
    dplyr::rename(start_Position = position)
  
  # make output dir
  output_dir <- paste0(output_dir, sub_dir, "/")
  
  # ensure directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # export with sample_id in filename
  output_file <- paste0(output_dir, sample_id, ".txt")

  # if output file already exists, remove it
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
cat("Deleting existing file:", sample, "\n")
cat("Writing new file:", sample, "\n")

  write.table(test_maf, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# make MAF files for IMPACT samples --------------------------------------------
make_maf_impact <- function(sample, sub_dir, output_dir){
  
  # DEBUG ----
  if (FALSE){
    sample <- "P-0005858-T03-IM7"
    sub_dir <- "impact"
  }
  
  # filter for sample of interest
  sample_impact <- impact_mod %>% 
    filter(tumor_sample_barcode == sample)
  
  # filter for SNV and Indel
  sample_impact <- sample_impact %>% 
    filter(variant_type %in% c("SNP", "INS", "DEL"))
  
  # define relevant columns
  columns_to_keep <- c("chromosome", "start_position", "reference_allele", "tumor_seq_allele2", "tumor_sample_barcode")
  
  # select relevant columns
  sample_impact <- sample_impact %>% 
    select(all_of(columns_to_keep)) %>% 
    dplyr::rename(start_Position = start_position)
  
  # make output dir
  output_dir <- paste0(output_dir, sub_dir, "/")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # export with sample_id in filename
  output_file <- paste0(output_dir, sample, ".txt")
  
  # if output file already exists, remove it
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  cat("Deleting existing file:", sample, "\n")
  cat("Writing new file:", sample, "\n")
  
  write.table(sample_impact, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
