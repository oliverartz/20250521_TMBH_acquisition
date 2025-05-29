#!/bin/bash
#BSUB -J snpEff_job
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4]"
#BSUB -W 2:00

# Load required module
module load java/11.0.12

# Check if the input arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <sample> <input_vcf> <output_dir>"
    exit 1
fi

# Get input arguments
sample="$1"
input_vcf="$2"
output_dir="$3"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Run snpEff to get protein sequences
java -jar /home/artzo/bin/snpEff/snpEff.jar \
-fastaProt "$output_dir/${sample}_protein_seq.fasta" \
-canon \
-onlyProtein \
-strict GRCh37.75 "$input_vcf" > "$output_dir/${sample}_annotated.vcf"

# Convert VCF to TSV for easier parsing
java -jar /home/artzo/bin/snpEff/SnpSift.jar extractFields "$output_dir/${sample}_annotated.vcf" \
CHROM POS REF ALT ANN[0].FEATUREID ANN[0].HGVS_P > "$output_dir/${sample}_annotated.vcf.tsv"
