#!/bin/bash
# =============================================================================
# GWAS PIPELINE FOR LGC INTERCEPT PHENOTYPES - EUR ONLY
# Prepare pipeline and test on single phenotype
# Author: Ada Kepinska
# Date: 2025-02-17
#
# Description:
#   Prepares a plink2-formatted phenotype file from LGC intercept output
#   and runs a test GWAS on a single phenotype to verify the pipeline.
#
# Prerequisites:
#   - LGC intercepts computed for EUR-ancestry participants (output of LGC pipeline)
#   - QC'd genotype files in plink binary format (.bed/.bim/.fam), EUR subset
#   - plink2 available as a module on your HPC (or adjust ml command below)
#
# Usage:
#   bash 1_test_GWAS_single_phenotype.sh
# =============================================================================

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR="/path/to/your/project"

# Input: LGC intercepts file (CSV, src_subject_id + phenotype columns)
# Produced by your LGC pipeline; must contain EUR participants only
LGC_INPUT_CSV="${PROJECT_DIR}/data/lgc_intercepts/LGC_intercepts_GWAS_ready_EUR.csv"

# QC'd genotype files (plink binary format, EUR ancestry, merged chromosomes)
# Must include ancestry PCs (V1-V10) and covariates in the .fam or separate file
GENO="${PROJECT_DIR}/data/genotypes/EUR/merged_chroms_rsIDs_EUR_final_QC"

# Output directory for GWAS results
OUT_DIR="${PROJECT_DIR}/results/GWAS/Mental_health"
mkdir -p "${OUT_DIR}"

# Formatted phenotype file (created by the R block below)
PHENO_FILE="${OUT_DIR}/LGC_intercepts_GWAS_ready_EUR_plink2_final.txt"

# Test phenotype — change to any column name in your LGC intercepts file
PHENO="upps_y_ss_negative_urgency_intercept"

# =============================================================================
# STEP 1: Format phenotype file for plink2
# =============================================================================
# plink2 requires:
#   - Columns FID and IID as the first two columns
#   - Tab-separated, NA for missing values

cd "${OUT_DIR}" || exit 1

# Load plink2 (adjust module name for your HPC environment)
ml PLINK/2_avx2_20221024

Rscript - <<EOF
df <- read.csv("${LGC_INPUT_CSV}",
               header    = TRUE,
               na.strings = c("", "NA"))

# Rename src_subject_id to IID and add FID (plink2 requires both)
colnames(df)[colnames(df) == "src_subject_id"] <- "IID"
df\$FID <- df\$IID

# Reorder so FID and IID are first
df <- df[, c("FID", "IID", setdiff(colnames(df), c("FID", "IID")))]

write.table(df,
            "${PHENO_FILE}",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE,
            na        = "NA")

cat("Phenotype file written:", nrow(df), "rows,", ncol(df), "columns\n")
EOF

# =============================================================================
# STEP 2: Run test GWAS on single phenotype
# =============================================================================
# Covariates: sex, age at baseline, PCs 1-10 (V1-V10)
# PCs are variance-standardized to avoid numerical issues in plink2
# --glm hide-covar: outputs only the phenotype association, not covariate stats

echo "Running test GWAS for: ${PHENO}"

plink2 \
    --bfile "${GENO}" \
    --pheno "${PHENO_FILE}" \
    --pheno-name "${PHENO}" \
    --covar "${PHENO_FILE}" \
    --covar-name demo_sex_v2 age_baseline V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 \
    --covar-variance-standardize V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 \
    --glm hide-covar \
    --out "${OUT_DIR}/test_${PHENO}"

echo "Test GWAS complete. Check ${OUT_DIR}/test_${PHENO}.*.glm.linear"
