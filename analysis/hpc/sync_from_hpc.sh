#!/bin/bash

# Manually sync intermediate files to HPC for further processing
local_dir="/Users/artzo/Documents/4_projects/20241016_TKI_TMB_increase/20250521_TMBH_acquisition/"

# VCF files converted from MAF
vcf_dir="/data/ldiaz/artzo/analysis/20250521_TMBH_acquisition/data/processed/hpc/vcf"

mkdir -p ${local_dir}/data/processed/hpc/vcf
rsync -avzPh artzo@lilac:${vcf_dir}/* ${local_dir}/data/processed/hpc/vcf/