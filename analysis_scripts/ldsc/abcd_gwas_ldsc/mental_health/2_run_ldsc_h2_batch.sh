#!/bin/bash
#SBATCH --job-name=LDSC_h2
#SBATCH --output=logs/ldsc_h2_%A_%a.out
#SBATCH --error=logs/ldsc_h2_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --partition=day
# =============================================================================
# LDSC HERITABILITY PIPELINE FOR LGC INTERCEPT PHENOTYPES - EUR ONLY
# Step 2: Munge sumstats and estimate h2
# Author: Ada Kepinska
# Date: 2025-02-28
#
# Description:
#   SLURM array job. Each task:
#     (a) Munges LDSC-formatted summary statistics with munge_sumstats.py
#     (b) Estimates SNP heritability with ldsc.py --h2
#
# Submission:
#   Run 3_submit_ldsc_h2_batch.sh — do not submit this script directly.
#
# Prerequisites:
#   - LDSC-formatted .txt.gz files from 1_reformat_plink2_for_ldsc.R
#   - LDSC installed (https://github.com/bulik/ldsc)
#   - EUR LD scores, HM3 SNP list, and LD weights (see 00_download_ldsc_setup.sh)
#   - ldsc conda environment active (Python 2.7 + dependencies)
# =============================================================================

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR="/path/to/your/project"

# LDSC-formatted summary statistics (output of 1_reformat_plink2_for_ldsc.R)
LDSC_INPUT="${PROJECT_DIR}/results/GWAS/Mental_health/LDSC_formatted"

# Output directory for munged sumstats and h2 results
LDSC_OUTPUT="${PROJECT_DIR}/results/LDSC/Mental_health"

# Phenotype list (one name per line)
PHENO_LIST="${PROJECT_DIR}/results/GWAS/Mental_health/phenotype_list_N3000.txt"

# LDSC software directory (cloned from https://github.com/bulik/ldsc)
LDSC_DIR="${PROJECT_DIR}/software/ldsc"

# LD reference files (downloaded/set up by 00_download_ldsc_setup.sh)
LD_SCORES="${PROJECT_DIR}/software/eur_w_ld_chr"          # Pre-computed EUR LD scores
W_LD="${PROJECT_DIR}/software/weights_hm3_no_hla"         # LD score weights
HM3_SNPS="${PROJECT_DIR}/software/w_hm3.snplist"          # HapMap3 SNP list

# =============================================================================
# SETUP
# =============================================================================

# Activate ldsc conda environment
# Adjust activation method for your HPC (module load miniconda, source activate, etc.)
module load miniconda
conda activate ldsc

mkdir -p "${LDSC_OUTPUT}/munged"
mkdir -p "${LDSC_OUTPUT}/h2"
mkdir -p "${LDSC_OUTPUT}/logs"

# =============================================================================
# GET PHENOTYPE FOR THIS ARRAY TASK
# =============================================================================

PHENO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PHENO_LIST}")

echo "============================================="
echo "Processing phenotype: ${PHENO}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "============================================="

# Check input file exists
INPUT_FILE="${LDSC_INPUT}/${PHENO}_ldsc.txt.gz"
if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: Input file not found: ${INPUT_FILE}"
    exit 1
fi

# =============================================================================
# STEP 2a: MUNGE SUMMARY STATISTICS
# =============================================================================

echo ""
echo "--- Munging summary statistics ---"

python "${LDSC_DIR}/munge_sumstats.py" \
    --sumstats "${INPUT_FILE}" \
    --out "${LDSC_OUTPUT}/munged/${PHENO}" \
    --merge-alleles "${HM3_SNPS}" \
    --signed-sumstats BETA,0

MUNGE_EXIT=$?
if [ ${MUNGE_EXIT} -ne 0 ]; then
    echo "ERROR: munge_sumstats.py failed for ${PHENO} (exit code ${MUNGE_EXIT})"
    exit 1
fi

MUNGED_FILE="${LDSC_OUTPUT}/munged/${PHENO}.sumstats.gz"
if [ ! -f "${MUNGED_FILE}" ]; then
    echo "ERROR: Munged file not created: ${MUNGED_FILE}"
    exit 1
fi

echo "Munging complete: ${MUNGED_FILE}"

# =============================================================================
# STEP 2b: ESTIMATE SNP HERITABILITY
# =============================================================================

echo ""
echo "--- Estimating h2 ---"

python "${LDSC_DIR}/ldsc.py" \
    --h2 "${MUNGED_FILE}" \
    --ref-ld-chr "${LD_SCORES}/" \
    --w-ld-chr "${W_LD}/weights." \
    --out "${LDSC_OUTPUT}/h2/${PHENO}_h2"

H2_EXIT=$?
if [ ${H2_EXIT} -ne 0 ]; then
    echo "ERROR: ldsc.py --h2 failed for ${PHENO} (exit code ${H2_EXIT})"
    exit 1
fi

echo ""
echo "--- Heritability results ---"
grep -A 5 "Heritability" "${LDSC_OUTPUT}/h2/${PHENO}_h2.log"

echo ""
echo "Done with ${PHENO} at $(date)"
