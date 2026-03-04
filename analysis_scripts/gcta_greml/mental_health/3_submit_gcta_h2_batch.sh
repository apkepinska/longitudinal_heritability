#!/bin/bash
# =============================================================================
# GCTA-GREML HERITABILITY - SUBMISSION SCRIPT
# Author: Ada Kepinska
# Date: 2025-02-26
#
# Description:
#   1. Creates the GCTA covariate file (extracted from the phenotype file)
#   2. Counts phenotypes in the list and submits 2_gcta_h2_batch.sh
#      as a SLURM array job with one task per phenotype.
#
# Usage:
#   bash 3_submit_gcta_h2_batch.sh
# =============================================================================

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR="/path/to/your/project"

# Phenotype file (FID, IID, then phenotype and covariate columns, tab-separated)
PHENO_FILE="${PROJECT_DIR}/results/GWAS/Mental_health/LGC_intercepts_GWAS_ready_EUR_plink2_final.txt"

# List of phenotype names to run (one per line)
PHENO_LIST="${PROJECT_DIR}/results/GWAS/Mental_health/phenotype_list_N3000.txt"

# Output directory for heritability results (must match path in 2_gcta_h2_batch.sh)
OUT_DIR="${PROJECT_DIR}/results/gcta/heritability"

# Path to the batch script
BATCH_SCRIPT="$(dirname "$0")/2_gcta_h2_batch.sh"

# =============================================================================
# SETUP
# =============================================================================

mkdir -p "${OUT_DIR}/logs"

# =============================================================================
# STEP 1: Extract covariate file from phenotype file
# =============================================================================
# GCTA requires covariates as a separate file (FID, IID, covariate columns).
# Covariates: sex, age at baseline, PCs 1-10 (V1-V10).
# Column positions assume the phenotype file was created by the GWAS pipeline
# (1=FID, 2=IID, 3=demo_sex_v2, 4=age_baseline, 5-14=V1-V10).
# Adjust column indices if your file has a different structure.

echo "Creating covariate file..."
awk 'BEGIN{OFS="\t"}
     NR==1{print "FID","IID","demo_sex_v2","age_baseline","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"}
     NR>1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' \
    "${PHENO_FILE}" > "${OUT_DIR}/covariates_gcta.txt"

echo "Covariate file written: ${OUT_DIR}/covariates_gcta.txt"

# =============================================================================
# STEP 2: Submit array job
# =============================================================================

N_PHENOS=$(wc -l < "${PHENO_LIST}")
echo "Submitting ${N_PHENOS} GCTA heritability jobs"

sbatch --array=1-${N_PHENOS} "${BATCH_SCRIPT}"
