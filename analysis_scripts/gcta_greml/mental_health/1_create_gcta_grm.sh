#!/bin/bash
#SBATCH --job-name=gcta_grm
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=day
#SBATCH --output=gcta_grm_%j.out
#SBATCH --error=gcta_grm_%j.err
# =============================================================================
# GCTA GRM
# Create genetic relatedness matrix (GRM) for GCTA-GREML heritability analysis
# Author: Ada Kepinska
# Date: 2025-02-17
#
# Description:
#   Computes a whole-genome GRM from QC'd plink binary genotype files.
#   This GRM is used downstream by 2_gcta_h2_batch.sh for REML heritability
#   estimation across phenotypes.
#
# Usage:
#   sbatch 1_create_gcta_grm.sh
#
# Prerequisites:
#   - QC'd genotype files in plink binary format (.bed/.bim/.fam), EUR subset
#   - GCTA available as a module on your HPC (or adjust ml command below)
# =============================================================================

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR="/path/to/your/project"

# QC'd genotype files (plink binary format, EUR ancestry, merged chromosomes)
GENO="${PROJECT_DIR}/data/genotypes/EUR/merged_chroms_rsIDs_EUR_final_QC"

# Output directory for GRM files
OUT_DIR="${PROJECT_DIR}/data/gcta/grm"
mkdir -p "${OUT_DIR}"

# =============================================================================
# LOAD DEPENDENCIES
# =============================================================================

# Adjust module name for your HPC environment
ml GCTA/1.94.1-gfbf-2022b

# =============================================================================
# CREATE GRM
# =============================================================================
# Outputs: <OUT_DIR>/abcd_eur_grm.grm.bin, .grm.N.bin, .grm.id

gcta64 \
    --bfile "${GENO}" \
    --make-grm \
    --out "${OUT_DIR}/eur_grm"

echo "GRM created: ${OUT_DIR}/eur_grm"
