#!/bin/bash

# Manually sync intermediate files to HPC for further processing
hpc_dir="/data/ldiaz/artzo/analysis/20250521_TMBH_acquisition"

# MAF files
maf_dir="/Users/artzo/Documents/4_projects/20241016_TKI_TMB_increase/20250521_TMBH_acquisition/data/processed/006_maf"

rsync -avzPh ${maf_dir}/* artzo@lilac:${hpc_dir}/data/processed/006_maf

# patient key
patient_key_dir="/Users/artzo/Documents/4_projects/20241016_TKI_TMB_increase/20250521_TMBH_acquisition/data/processed/001_patient_key_mod.txt"

rsync -avzPh ${patient_key_dir} artzo@lilac:${hpc_dir}/data/processed/001_patient_key_mod.txt

# annotated mutations (progression samples)

mutations_file="/Users/artzo/Documents/4_projects/20241016_TKI_TMB_increase/20250521_TMBH_acquisition/data/processed/005_annotated_mutations_progression.txt"

rsync -avzPh ${mutations_file} artzo@lilac:${hpc_dir}/data/processed/005_annotated_mutations_progression.txt
