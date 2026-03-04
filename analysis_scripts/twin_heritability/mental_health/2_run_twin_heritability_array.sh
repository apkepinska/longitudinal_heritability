#!/bin/bash
#SBATCH --job-name=twin_h2
#SBATCH --output=logs/h2_%A_%a.out
#SBATCH --error=logs/h2_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --partition=day
# =============================================================================
# TWIN HERITABILITY - SLURM ARRAY JOB
# Author: Ada Kepinska
# Date: 2025-02-28
#
# Each array task runs ACE + HLM for one phenotype × timepoint combination.
# The job list file maps SLURM_ARRAY_TASK_ID to phenotype + timepoint.
#
# Submission:
#   Run 3_submit_twin_heritability_array.sh — do not submit this script
#   directly, as the array size must first be computed from the job list.
# =============================================================================

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR="/path/to/your/project"

# Directory containing this script and 1_twin_heritability_single_phenotype.R
SCRIPT_DIR="$(dirname "$0")"

# Directory where per-phenotype results and the job list are stored
RESULTS_DIR="${PROJECT_DIR}/results/twin_heritability/cross_sectional"

# =============================================================================
# SETUP
# =============================================================================

# Load R (adjust module name for your HPC environment)
module load R/4.4.1

# Pass job list path to the R script via environment variable
export JOB_LIST_FILE="${RESULTS_DIR}/heritability_job_list.csv"

echo "=============================================="
echo "TWIN HERITABILITY - Task ${SLURM_ARRAY_TASK_ID}"
echo "Job list: ${JOB_LIST_FILE}"
echo "Start: $(date)"
echo "=============================================="

# =============================================================================
# RUN
# =============================================================================

Rscript "${SCRIPT_DIR}/1_twin_heritability_single_phenotype.R"

echo "=============================================="
echo "End: $(date)"
echo "=============================================="
