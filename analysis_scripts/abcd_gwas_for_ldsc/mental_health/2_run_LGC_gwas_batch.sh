#!/bin/bash
#SBATCH --job-name=LGC_GWAS_EUR
#SBATCH --output=logs/gwas_%A_%a.out
#SBATCH --error=logs/gwas_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --partition=day
# =============================================================================
# GWAS PIPELINE FOR LGC INTERCEPT PHENOTYPES - EUR ONLY (BATCH)
# Author: Ada Kepinska
# Date: 2025-02-16
#
# Description:
#   SLURM array job. Each task runs a GWAS for one LGC intercept phenotype.
#   Phenotype names are read from a list file (one name per line), and the
#   task ID indexes into that list.
#
# Submission:
#   Run 3_submit_job_LGC_gwas_batch.sh — do not submit this script directly,
#   as it requires the array size to be computed from the phenotype list first.
#
# Prerequisites:
#   - Phenotype file created by 1_test_GWAS_single_phenotype.sh
#   - Phenotype list file (one phenotype name per line, N >= 3000 participants)
#   - QC'd EUR genotype files in plink binary format
# =============================================================================

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR="/path/to/your/project"

# QC'd genotype files (plink binary format, EUR ancestry)
GENO="${PROJECT_DIR}/data/genotypes/EUR/merged_chroms_rsIDs_EUR_final_QC"

# Phenotype file formatted for plink2
# Created by 1_test_GWAS_single_phenotype.sh
PHENO_FILE="${PROJECT_DIR}/results/GWAS/Mental_health/LGC_intercepts_GWAS_ready_EUR_plink2_final.txt"

# List of LGC intercept phenotype names to run GWAS on (one per line)
# Should contain only phenotypes with sufficient sample size (e.g. N >= 3000)
PHENO_LIST="${PROJECT_DIR}/results/GWAS/Mental_health/phenotype_list_N3000.txt"

# Output directory
OUT_DIR="${PROJECT_DIR}/results/GWAS/Mental_health"

# =============================================================================
# LOAD DEPENDENCIES
# =============================================================================

# Adjust module name for your HPC environment
ml PLINK/2_avx2_20221024

# =============================================================================
# RUN GWAS FOR THIS ARRAY TASK
# =============================================================================

# Get phenotype name for this array task (1-indexed)
PHENO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PHENO_LIST}")

echo "Running GWAS for: ${PHENO}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# Covariates: sex, age at baseline, PCs 1-10 (V1-V10)
# PCs variance-standardised to avoid numerical issues in plink2
plink2 \
    --bfile "${GENO}" \
    --pheno "${PHENO_FILE}" \
    --pheno-name "${PHENO}" \
    --covar "${PHENO_FILE}" \
    --covar-name demo_sex_v2 age_baseline V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 \
    --covar-variance-standardize V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 \
    --glm hide-covar \
    --out "${OUT_DIR}/${PHENO}"

# Compress output
gzip "${OUT_DIR}/${PHENO}".*.glm.linear

echo "Completed: ${PHENO}"
