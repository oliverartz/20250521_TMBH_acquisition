#!/bin/bash

# Manually sync intermediate files to HPC for further processing
local_dir="/Users/artzo/Documents/4_projects/20241016_TKI_TMB_increase/20250521_TMBH_acquisition/"

# VCF files converted from MAF
vcf_dir="/data/ldiaz/artzo/analysis/20250521_TMBH_acquisition/data/processed/hpc/vcf"

mkdir -p ${local_dir}/data/processed/hpc/vcf
rsync -avzPh artzo@lilac:${vcf_dir}/* ${local_dir}/data/processed/hpc/vcf/

# Summarized netMHCpan results
netmhcpan_results="/data/ldiaz/artzo/analysis/20250521_TMBH_acquisition/data/processed/hpc/006_patient_status_compared_to_baseline_phbr_vaf.txt"

mkdir -p ${local_dir}/data/processed/hpc
rsync -avzPh artzo@lilac:${netmhcpan_results} ${local_dir}/data/processed/hpc/
