## -----------------------------------------------------------------------------
## Purpose of script: Define functions used in data preprocessing steps
##
## Author: Oliver Artz
## Date Created: May 23, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# annotate_patient_mutations function
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


# DEBUG ----
if (FALSE) {
  guardant_data <- guardant_mod
  impact_data <- impact_mod
  gene_panel_list <- panel_gene_list
  patient_ids <- unique(guardant_mod$patient_id)
  patient_id_of_interest <- 14
}

# END DEBUG ----

annotate_patient_mutations <- function(guardant_data, impact_data, gene_panel_list, patient_ids) {
  
  # define a helper function for a single patient
  process_single_patient <- function(patient_id_of_interest) {
    
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