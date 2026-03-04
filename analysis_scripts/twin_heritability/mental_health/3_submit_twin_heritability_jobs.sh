#!/bin/bash
# =============================================================================
# SUBMIT TWIN HERITABILITY JOBS
# Author: Ada Kepinska
# Date: 2025-02-28
#
# Reads the power assessment output, selects phenotypes with good or adequate
# power, generates a job list (phenotype × timepoint), and submits array job.
#
# Prerequisites:
#   - twin_phenotype_availability_by_timepoint.csv from
#     Processing_scripts/Twin_heritability/2_check_twin_phenotype_availability_mh.R
# =============================================================================

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR="/path/to/your/project"

# Input: phenotype × timepoint availability file (from power check script)
AVAILABILITY_FILE="${PROJECT_DIR}/outputs/twin_identification/twin_phenotype_availability_by_timepoint.csv"

# Output directory for per-phenotype results and logs
RESULTS_DIR="${PROJECT_DIR}/results/twin_heritability/cross_sectional"

# Path to the array batch script (sibling of this script)
BATCH_SCRIPT="$(dirname "$0")/2_run_twin_heritability_array.sh"

# =============================================================================
# SETUP
# =============================================================================

mkdir -p "${RESULTS_DIR}/logs"

# Job list path (passed to array job via JOB_LIST_FILE env var)
JOB_LIST="${RESULTS_DIR}/heritability_job_list.csv"

echo "=============================================="
echo "TWIN HERITABILITY - JOB SUBMISSION"
echo "=============================================="
echo ""

# =============================================================================
# STEP 1: CHECK INPUT
# =============================================================================

if [ ! -f "${AVAILABILITY_FILE}" ]; then
    echo "ERROR: Availability file not found: ${AVAILABILITY_FILE}"
    echo "Run 2_check_twin_phenotype_availability_mh.R first."
    exit 1
fi

# =============================================================================
# STEP 2: GENERATE JOB LIST (filter to well-powered phenotypes)
# =============================================================================

echo "Generating job list from: ${AVAILABILITY_FILE}"

module load R/4.4.1

Rscript -e "
library(data.table)
dt <- fread('${AVAILABILITY_FILE}')

# Keep phenotype × timepoint combinations with good or adequate power
well_powered <- dt[power_level %in% c('good', 'adequate')]

job_list <- well_powered[, .(phenotype, timepoint, power_level,
                             mz_complete_pairs, dz_complete_pairs,
                             total_complete_pairs)]

# Sort for reproducibility
job_list <- job_list[order(timepoint, phenotype)]

fwrite(job_list, '${JOB_LIST}')

cat('Job list created:', nrow(job_list), 'phenotype x timepoint combinations\n')
cat('  Good power:    ', sum(job_list\$power_level == 'good'), '\n')
cat('  Adequate power:', sum(job_list\$power_level == 'adequate'), '\n')
cat('  By timepoint:\n')
print(job_list[, .N, by = timepoint])
"

# =============================================================================
# STEP 3: VALIDATE JOB LIST
# =============================================================================

N_JOBS=$(tail -n +2 "${JOB_LIST}" | wc -l)

if [ "${N_JOBS}" -eq 0 ]; then
    echo "ERROR: No well-powered phenotype x timepoint combinations found."
    exit 1
fi

echo ""
echo "Total jobs to submit: ${N_JOBS}"
echo "Job list saved to: ${JOB_LIST}"
echo ""

# =============================================================================
# STEP 4: CHECK FOR ALREADY-COMPLETED RESULTS
# =============================================================================

N_DONE=0
while IFS=, read -r pheno tp rest; do
    [ "$pheno" = "phenotype" ] && continue   # skip header
    if [ -f "${RESULTS_DIR}/${pheno}_${tp}_heritability.csv" ]; then
        N_DONE=$((N_DONE + 1))
    fi
done < "${JOB_LIST}"

if [ ${N_DONE} -gt 0 ]; then
    echo "NOTE: ${N_DONE}/${N_JOBS} jobs already have results."
    echo "  These will be overwritten. To skip completed jobs,"
    echo "  filter them from the job list manually before submitting."
    echo ""
fi

# =============================================================================
# STEP 5: SUBMIT ARRAY JOB
# =============================================================================

echo "Submitting SLURM array job (1-${N_JOBS})..."
sbatch --array=1-${N_JOBS} "${BATCH_SCRIPT}"

echo ""
echo "Done! Monitor with: squeue -u $(whoami) -n twin_h2"
echo "Results will be in: ${RESULTS_DIR}/"
