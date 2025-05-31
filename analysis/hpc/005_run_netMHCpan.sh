#!/bin/bash

# define parameters ------------------------------------------------------------
# change working directory
wkdir="/data/ldiaz/artzo/analysis/20250521_TMBH_acquisition"
cd $wkdir

# netMHCpan output directory
output_dir_netMHCpan="data/processed/hpc/netMHCpan/"

# path haplotype data
haplotype_data="data/processed/hpc/hla_haplotypes.txt"

# sample group name
group_name="patient_status_compared_to_baseline"

# load samples and haplotypes into bash array from hla_haplotypes.txt
mapfile -t samples < <(awk '{print $1}' ${haplotype_data})
mapfile -t haplotypes < <(awk '{print $2}' ${haplotype_data})

echo "Samples: ${samples[@]}"
echo "Haplotypes: ${haplotypes[@]}"

# run netMHCpan as array job ---------------------------------------------------
# get number of total samples to calculate progress
total_samples=${#samples[@]}

echo "Total samples: $total_samples"

# create array job script
bsub -J "netMHCpan[1-${total_samples}]" -n 4 -R "rusage[mem=4]" -W 2:00 <<EOF
#!/bin/bash
#BSUB -J netMHCpan[1-${total_samples}]
#BSUB -n 4
#BSUB -R "rusage[mem=4]"
#BSUB -W 2:00
#BSUB -o analysis/hpc/logs/netMHCpan.%J.log
#BSUB -e analysis/hpc/logs/netMHCpan.%J.err

# define path to netMHCpan
netMHCpan="/lila/data/ldiaz/artzo/bin/netMHCpan-4.1/"

# change to working directory
cd "${wkdir}"

# reload arrays from file inside job
mapfile -t samples < <(awk '{print \$1}' ${haplotype_data})
mapfile -t haplotypes < <(awk '{print \$2}' ${haplotype_data})

# define output directory
output_dir="${output_dir_netMHCpan}/${group_name}/"
mkdir -p \${output_dir}

# get sample and haplotype for this job
i=\$((LSB_JOBINDEX-1))
sample_name=\${samples[\$i]}
input_fasta="data/processed/hpc/fasta/${group_name}/\${sample_name}.fasta"
hla_haplotype=\${haplotypes[\$i]}

\${netMHCpan}/netMHCpan \\
-f \${input_fasta} \\
-a \${hla_haplotype} \\
-l 8,9,10,11 \\
-BA -s -xls \\
-xlsfile \${output_dir}\${sample_name}_\${hla_haplotype}.xls \\
> \${output_dir}\${sample_name}_\${hla_haplotype}.txt
EOF