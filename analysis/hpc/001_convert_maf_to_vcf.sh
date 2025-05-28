#!/bin/bash

# Convert MAF to VCF using vcf2maf

# to run from interactive session
# bsub -n 2 -R "rusage[mem=4]" -W 4:00 -Is bash

source ~/.bashrc

# change to the correct working directory
cd /data/ldiaz/artzo/analysis/20250521_TMBH_acquisition

# define the path to the maf2vcf script
maf2vcf_path="/home/artzo/bin/vcf2maf/mskcc-vcf2maf-*/maf2vcf.pl"

# check if maf2vcf script exists
if ! ls ${maf2vcf_path} 1> /dev/null 2>&1; then
    echo "Error: maf2vcf script not found at ${maf2vcf_path}"
    exit 1
fi

# define the path to the reference genome
ref_fasta="/data/ldiaz/artzo/assets/vep/homo_sapiens/112_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz"

# check if reference genome exists
if [[ ! -f ${ref_fasta} ]]; then
    echo "Error: Reference genome not found at ${ref_fasta}"
    exit 1
fi

# define the directories to process
directories=("patient_status_compared_to_baseline" "gh_request_id" "tmbh_status_compared_to_baseline")

# loop through each directory
for dir in "${directories[@]}"; do
    echo "Processing directory: ${dir}"
    
    # define the path to the maf file
    input_dir="data/processed/006_maf/${dir}"
    output_dir="data/processed/hpc/vcf/${dir}"
    mkdir -p ${output_dir}
    
    # run maf2vcf for all files in the directory
    for input_maf in ${input_dir}/*.txt; do
        echo ${input_maf}
        perl ${maf2vcf_path} --input-maf ${input_maf} --output-dir ${output_dir} --ref-fasta ${ref_fasta}
    done
done