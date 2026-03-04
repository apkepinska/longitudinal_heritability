#!/bin/bash
#SBATCH --job-name=gcta_h2_batch
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --partition=day
#SBATCH --output=logs/gcta_h2_%A_%a.out
#SBATCH --error=logs/gcta_h2_%A_%a.err
# =============================================================================
# GCTA-GREML HERITABILITY - BATCH (ALL PHENOTYPES)
# Author: Ada Kepinska
# Date: 2025-02-26
#
# Description:
#   SLURM array job. Each task estimates SNP heritability (h2) for one LGC
#   intercept phenotype using GCTA-GREML (REML). Phenotype names are read
#   from a list file and indexed by SLURM_ARRAY_TASK_ID.
#
# Submission:
#   Run 3_submit_gcta_h2_batch.sh — do not submit this script directly, as
#   the array size must first be computed from the phenotype list.
#
# Prerequisites:
#   - GRM created by 1_create_gcta_grm.sh
#   - Phenotype file formatted for plink2/GCTA (FID, IID, phenotype columns)
#   - Covariate file created by 3_submit_gcta_h2_batch.sh
#   - Phenotype list file (one name per line)
# =============================================================================

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR="/path/to/your/project"

# GRM prefix (output of 1_create_gcta_grm.sh, without file extension)
GRM_DIR="${PROJECT_DIR}/data/gcta/grm"

# Phenotype file (FID, IID, then phenotype and covariate columns, tab-separated)
# Created by 1_test_GWAS_single_phenotype.sh in the GWAS pipeline
PHENO_FILE="${PROJECT_DIR}/results/GWAS/Mental_health/LGC_intercepts_GWAS_ready_EUR_plink2_final.txt"

# List of phenotype names to run (one per line, e.g. filtered to N >= 3000)
PHENO_LIST="${PROJECT_DIR}/results/GWAS/Mental_health/phenotype_list_N3000.txt"

# Output directory for heritability results (.hsq files)
OUT_DIR="${PROJECT_DIR}/results/gcta/heritability"

# =============================================================================
# LOAD DEPENDENCIES
# =============================================================================

ml GCTA/1.94.1-gfbf-2022b

# =============================================================================
# RUN GCTA-GREML FOR THIS ARRAY TASK
# =============================================================================

PHENO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PHENO_LIST}")
echo "Processing phenotype: ${PHENO} (Task ID: ${SLURM_ARRAY_TASK_ID})"

# Find GCTA column number for this phenotype
# GCTA --mpheno uses 1-based indexing starting AFTER the FID and IID columns
COL_NUM=$(head -1 "${PHENO_FILE}" | tr '\t' '\n' | nl | grep -w "${PHENO}" | awk '{print $1-2}')
echo "Phenotype ${PHENO} is in column ${COL_NUM} (GCTA numbering)"

# Run REML heritability estimation
# Covariates: sex, age at baseline, PCs 1-10 (V1-V10) — extracted by submission script
echo "Running GCTA-GREML..."
gcta64 \
    --reml \
    --grm "${GRM_DIR}/eur_grm" \
    --pheno "${PHENO_FILE}" \
    --mpheno "${COL_NUM}" \
    --qcovar "${OUT_DIR}/covariates_gcta.txt" \
    --out "${OUT_DIR}/h2_${PHENO}" \
    --thread-num "${SLURM_CPUS_PER_TASK}"

echo "Completed: ${PHENO}"
echo "Results saved to: ${OUT_DIR}/h2_${PHENO}.hsq"
