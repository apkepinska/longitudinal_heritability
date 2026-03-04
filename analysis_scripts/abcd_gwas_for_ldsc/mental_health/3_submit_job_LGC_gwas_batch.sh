#!/bin/bash
# =============================================================================
# GWAS PIPELINE FOR LGC INTERCEPT PHENOTYPES - EUR ONLY
# Submission script
# Author: Ada Kepinska
# Date: 2025-02-16
#
# Description:
#   Counts the number of phenotypes in the list file and submits
#   2_run_LGC_gwas_batch.sh as a SLURM array job with one task per phenotype.
#
# Usage:
#   bash 3_submit_job_LGC_gwas_batch.sh
# =============================================================================

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR="/path/to/your/project"

# List of phenotype names to run (one per line, N >= 3000 participants)
# Same file referenced in 2_run_LGC_gwas_batch.sh
PHENO_LIST="${PROJECT_DIR}/results/GWAS/Mental_health/phenotype_list_N3000.txt"

# Path to the batch GWAS script
BATCH_SCRIPT="$(dirname "$0")/2_run_LGC_gwas_batch.sh"

# Logs directory (must exist before sbatch runs)
LOGS_DIR="${PROJECT_DIR}/results/GWAS/Mental_health/logs"

# =============================================================================
# SUBMIT
# =============================================================================

# Create logs directory if it does not exist
mkdir -p "${LOGS_DIR}"

N_PHENOS=$(wc -l < "${PHENO_LIST}")
echo "Submitting ${N_PHENOS} GWAS jobs"

sbatch --array=1-${N_PHENOS} "${BATCH_SCRIPT}"
