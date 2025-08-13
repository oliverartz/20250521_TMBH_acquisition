# Acquired high tumor mutational burden and activity of immunotherapy after targeted therapy in microsatellite stable colorectal cancer

This repository contains the code and analysis pipeline for the publication "Acquired high tumor mutational burden and activity of immunotherapy after targeted therapy in microsatellite stable colorectal cancer" by Yeh et al.

## Overview

This study investigates the acquisition of high tumor mutational burden (TMB-H) in microsatellite stable (MSS) colorectal cancer patients following targeted therapy and its implications for subsequent immune checkpoint blockade response.

## Repository Structure

```
├── analysis/
│   ├── run_all.R                    # Master script to run entire analysis
│   ├── hpc/                         # High-performance computing scripts
│   ├── scripts/                     # Main analysis scripts
│   │   ├── utils/                   # Utility functions
│   │   └── 001-021 numbered scripts # Sequential analysis pipeline
│   └── publication/                 # Scripts to extract specific stats for main text
├── data/
│   ├── metadata/                    # Sample metadata and gene panels
│   ├── processed/                   # Processed data files
│   └── raw/                         # Raw sequencing data
└── results/
    ├── figures/                     # Generated figures
    └── tables/                      # Analysis results tables
└── test_data/                       # Synthetic, schema‑compatible inputs for demo
```

## Data availability and ethics

- Patient‑level raw and processed clinical data contain PHI and cannot be distributed.
- The repository ships with or supports synthetic “test_data” (no PHI) under `test_data/`.
- Real data access is governed by institutional approvals and data‑use agreements.

## Running the Analysis

### Dependencies

R packages (CRAN)
- pacman, tidyverse, data.table, here, janitor, readxl
- ggpubr, ggsci, RColorBrewer, viridis, msigdbr
- Optional (HPC helpers/QC): vcfR, phylotools, DescTools

Bioconductor
- MutationalPatterns, BSgenome.Hsapiens.UCSC.hg19

Install

```r
install.packages(c("pacman","tidyverse","data.table","here","janitor",
                   "readxl","ggpubr","ggsci","RColorBrewer","viridis",
                   "msigdbr","vcfR","phylotools","DescTools"))
install.packages("BiocManager")
BiocManager::install(c("MutationalPatterns","BSgenome.Hsapiens.UCSC.hg19"))
```

System/HPC tools (for neoantigen pipeline)
- snpEff (GRCh37.75 database)
- netMHCpan
- LSF or equivalent scheduler, Java
- Paths and submission scripts are in [analysis/hpc](analysis/hpc).

### Quick Start

```r
# Run complete analysis pipeline
source("analysis/run_all.R")
```

**Note**: The pipeline includes an interactive checkpoint before HPC-dependent analyses (scripts 010-021). User confirmation is required to proceed with analyses that depend on HPC-processed files.

### HPC Pipeline for Neoantigen Prediction

**HPC Path**: `/data/ldiaz/artzo/analysis/20250521_TMBH_acquisition`

#### 1. Initial Setup

Run locally:

```bash
bash analysis/hpc/sync_to_hpc.sh
```

#### 2. HPC Processing Pipeline

Run on HPC:

```bash

# navigate to project directory

# Step 1: Convert MAF to VCF format
bash analysis/hpc/001_convert_maf_to_vcf.sh

# Step 2: Annotate VCF files with SNPEff to extract protein sequences
bash analysis/hpc/002_vcf_anno_snpeff.sh

# Step 3: Extract protein sequences from SNPEff annotations
Rscript analysis/hpc/003_extract_proteins_from_SNPEff.R

# Step 4: Extract HLA haplotypes and format for netMHCpan
Rscript analysis/hpc/004_extract_haplotypes_for_netMHCpan.R

# Step 5: Run netMHCpan for neoantigen binding prediction
bash analysis/hpc/005_run_netMHCpan.sh

# Step 6: Summarize netMHCpan results and calculate PHBR scores
Rscript analysis/hpc/006_summarize_netMHCpan.R
```

#### 3. Retrieve Results

Run locally:

```bash
bash analysis/hpc/sync_from_hpc.sh
```

## Script-by-Script Analysis Summary

### **001_preprocess_patient_key.R**

- **Output**: `001_patient_key_mod.txt` (processed patient metadata file)
- **Samples Used**: Master patient key linking all sample types
- **Purpose**: Creates unified patient identifier mapping between IMPACT and Guardant360 samples

### **002_preprocess_guardant_omni.R**

- **Output**: `002_guardant_mod.txt` (processed Guardant360 data)
- **Samples Used**: All Guardant360 OMNI liquid biopsy samples (baseline + progression)
- **Purpose**: Processes and cleans Guardant360 sequencing data

### **003_preprocess_patient_meta.R**

- **Output**: `003_patient_meta_mod.txt` (processed patient metadata)
- **Samples Used**: Clinical metadata for all patients
- **Purpose**: Standardizes patient clinical annotations

### **004_preprocess_impact.R**

- **Output**: `004_impact_mod.txt` (processed IMPACT data)
- **Samples Used**: IMPACT tissue sequencing samples (baseline + progression where available)
- **Purpose**: Processes MSK-IMPACT tissue sequencing data

### **005_annotate_acquired_mutations.R**

- **Output**: `005_annotated_mutations_progression.txt` and `005_annotated_mutations_all.txt`
- **Samples Used**: 
  - **Baseline**: IMPACT baseline + Guardant360 baseline
  - **Progression**: Guardant360 progression samples
- **Purpose**: Identifies acquired, baseline, and lost mutations by comparing timepoints

### **006_make_maf.R**

- **Output**: MAF files for downstream analysis (multiple groupings)
- **Samples Used**: All processed mutation data (IMPACT + Guardant360)
- **Purpose**: Generates standardized MAF files for signature analysis and neoantigen prediction

### **007_plot_TMB.R**

- **Output Plot**: `007_TMB.pdf` - TMB distribution by acquisition status
- **Samples Used**: IMPACT
- **Comparison**: TMB-H acquiring vs non-acquiring patients
- **Data Source**: Baseline tissue TMB measurements

### **008_plot_MSI.R**

- **Output Plot**: `008_MSI.pdf` - MSI score distribution
- **Samples Used**: IMPACT
- **Comparison**: TMB-H acquiring vs non-acquiring patients
- **Purpose**: Validates microsatellite stable (MSS) status of cohort

### **009_plot_variant_types_fractions.R**

- **Output Plot**: `009_variant_types_fractions.pdf` - Variant type fractions
- **Samples Used**: **IMPACT baseline tissue samples only**
- **Analysis**: Fraction of missense, frameshift, nonsense, in-frame, and other variants
- **Comparison**: TMB-H acquiring vs non-acquiring patients

### **010_sbs_per_sample.R**

- **Output Plot**: `010_SBS_per_sample.pdf` - Mutational signatures per sample
- **Samples Used**: **IMPACT baseline tissue samples** 
- **Analysis**: COSMIC SBS signature deconvolution for individual samples
- **Comparison**: TMB-H acquiring vs non-acquiring patients

### **011_sbs_acquired_baseline.R**

- **Output Plot**: `011_sbs_acquired_baseline.pdf` - Signatures in acquired vs baseline mutations
- **Samples Used**: **TMB-H acquiring patients only**
- **Comparison**: 
  - **Acquired mutations**: From Guardant360 progression samples
  - **Baseline mutations**: From IMPACT baseline + Guardant360 baseline (pooled)
- **Analysis**: Signature differences between acquired and baseline mutation sets

### **012_sbs_baseline_progression_samples.R**

- **Output Plot**: `012_sbs_baseline_progression_samples.pdf` - Signatures in paired samples
- **Samples Used**: **Guardant360 paired baseline and progression samples** from TMB-H acquiring patients
- **Analysis**: Evolution of mutational signatures between baseline and progression timepoints
- **Sample Type**: Liquid biopsy only (Guardant360)

### **013_sbs_acquired_per_patient.R**

- **Output Plot**: `013_sbs_acquired_per_patient.pdf` - Signatures in acquired mutations per patient
- **Samples Used**: **Acquired mutations** from TMB-H acquiring and non-acquiring patients
- **Source**: Mutations identified as acquired in Guardant360 progression samples
- **Analysis**: Patient-level signature analysis of newly acquired variants

### **014_mean_n_variants_per_patient.R**

- **Output Plot**: `014_mean_n_variants_per_patient.pdf` - Mean variant counts over time
- **Samples Used**: **All Guardant360 samples** (baseline + progression)
- **Analysis**: Changes in mutation burden between baseline and progression
- **Comparison**: TMB-H acquiring vs non-acquiring patients

### **015_pathways.R**

- **Output Plot**: `015_pathways.pdf` - Pathway enrichment heatmap
- **Samples Used**: **Acquired mutations** from TMB-H acquiring patients (SNV + Indel only)
- **Source**: Guardant360 progression samples
- **Analysis**: MSigDB REACTOME pathway enrichment (MAPK, PI3K, DNA repair, MMR)

### **016_pathways_SanchezVega2018.R**

- **Output Plot**: `016_pathways_SanchezVega2018.pdf` - Cancer pathway analysis
- **Samples Used**: **Acquired mutations** from TMB-H acquiring patients (SNV + Indel only)
- **Source**: Guardant360 progression samples
- **Analysis**: Canonical cancer pathways from Sanchez-Vega et al. 2018

### **017_acquired_alteration_types_progression.R**

- **Output Plot**: `017_acquired_alteration_types_progression.pdf` - Alteration type distribution
- **Samples Used**: **Acquired mutations only** from Guardant360 progression samples
- **Analysis**: Types of alterations (missense, nonsense, frameshifts, CNV subtypes) in acquired mutations
- **Comparison**: TMB-H acquiring vs non-acquiring patients

### **018_maxVAF.R**

- **Output Plot**: `018_maxVAF_histogram.pdf` - VAF distribution with clonality inference
- **Samples Used**: **Guardant360 progression samples** from TMB-H acquiring patients with VAF data
- **Analysis**: VAF corrected by maximum VAF per sample (clonality proxy)
- **Comparison**: Acquired vs baseline mutations

### **019_phbr_acquired_baseline.R**

- **Output Plot**: `019_phbr_acquired_baseline.pdf` - Neoantigen presentation comparison
- **Samples Used**: **TMB-H acquiring patients** with HLA typing data
- **Comparison**: 
  - **Acquired mutations**: From Guardant360 progression
  - **Baseline mutations**: From IMPACT baseline + Guardant360 baseline (overlapping)
- **Analysis**: HLA-binding rank (PHBR) for neoantigen presentation potential
- **Requires**: HPC-processed neoantigen predictions

### **020_phbr_clonality.R**

- **Output Plot**: `020_phbr_clonality.pdf` - Neoantigen presentation by clonality
- **Samples Used**: **Acquired mutations** from TMB-H acquiring patients with clonality classification
- **Analysis**: PHBR stratified by clonal vs subclonal status (using VAF thresholds)

### **021_blood_tissue_tmb.R**

- **Output Plot**: `021_blood_tissue_tmb.pdf` - Platform concordance analysis
- **Samples Used**: **Patients with both IMPACT baseline and Guardant360 baseline/progression samples**
- **Comparison**: 
  - **Tissue TMB**: IMPACT baseline samples
  - **Plasma TMB**: Guardant360 baseline and progression samples
- **Analysis**: TMB concordance between tissue and liquid biopsy platforms
- **Faceting**: By TMB-H acquisition status

## Contact

Questions: artzo@mskcc.org