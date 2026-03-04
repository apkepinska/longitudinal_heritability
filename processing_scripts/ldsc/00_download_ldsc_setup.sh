#!/bin/bash
# =============================================================================
# DOWNLOAD LDSC AND SETUP
# Author: Ada Kepinska
# Date: 2025-02-27
#
# Description:
#   Downloads and sets up LDSC (LD Score Regression) for heritability and
#   genetic correlation analyses. Requires access to LD reference files
#   (EUR LD scores, HM3 SNP list, weights) from an external source such as
#   the Bulik-Sullivan lab or your HPC software repository.
#
# Dependencies:
#   - conda
#   - git
#   - Access to LD reference files (see REFERENCE FILES section below)
#
# Usage:
#   bash 00_download_ldsc_setup.sh
# =============================================================================

# =============================================================================
# USER CONFIGURATION — edit these paths before running
# =============================================================================

PROJECT_DIR="${HOME}/project/your_project_name"          # Root of your project
SOFTWARE_DIR="${PROJECT_DIR}/Software"                   # Where LDSC will be cloned
LD_REF_DIR="/path/to/your/ld_reference_files"           # Directory with EUR LD scores etc.
                                                          # e.g. /gpfs/.../software on your HPC

# SLURM resource allocation (adjust as needed for your cluster)
PARTITION="day"
MEM="15G"
TIME="2:00:00"
CPUS=2

# =============================================================================
# REFERENCE FILES NEEDED
# These files should be available from your HPC software repository or
# downloaded from: https://github.com/bulik/ldsc (see their README)
#   - eur_w_ld_chr.tar.bz2     (EUR LD scores)
#   - w_hm3.snplist             (HM3 SNP list, inside eur_w_ld_chr/)
#   - weights_hm3_no_hla/       (LD score weights directory)
# =============================================================================

# ---------------------------------------------------------------------------
# Step 1: Clone LDSC
# ---------------------------------------------------------------------------
mkdir -p "${SOFTWARE_DIR}"
cd "${SOFTWARE_DIR}" || exit 1

git clone https://github.com/bulik/ldsc.git

# ---------------------------------------------------------------------------
# Step 2: Copy and extract LD reference files
# ---------------------------------------------------------------------------
cp "${LD_REF_DIR}/eur_w_ld_chr.tar.bz2" .
tar -xjf eur_w_ld_chr.tar.bz2

# Copy HM3 SNP list
cp "${LD_REF_DIR}/eur_w_ld_chr/w_hm3.snplist" .

# Copy LD score weights
cp -r "${LD_REF_DIR}/ldsc_references/weights_hm3_no_hla" .

# ---------------------------------------------------------------------------
# Step 3: Set up conda environment
# ---------------------------------------------------------------------------
# Allocate interactive resources (SLURM — adjust flags for your cluster)
echo "Run the following to allocate an interactive session:"
echo "  salloc --partition=${PARTITION} --mem=${MEM} --time=${TIME} --cpus-per-task=${CPUS}"

# Create and activate LDSC conda environment (Python 2.7 required by LDSC)
conda create -n ldsc python=2.7 numpy scipy pandas bitarray -c conda-forge -y
conda activate ldsc

# ---------------------------------------------------------------------------
# Step 4: Test installation
# ---------------------------------------------------------------------------
cd "${SOFTWARE_DIR}/ldsc" || exit 1
python ldsc.py -h # check if help opens
