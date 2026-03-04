######################
# Create a list of sufficiently powered LGC phenotypes to run as GWASs
# Author: Ada Kepinska
# Date: 2025-02-16
#
# Description:
#   Reads the GWAS-ready LGC intercept file, counts non-missing participants
#   per phenotype, and outputs a list of phenotypes meeting the minimum
#   sample size threshold for GWAS power.
#
# Prerequisites:
#   - LGC_intercepts_GWAS_ready_EUR.csv produced by:
#     Latent_growth_curve/Mental_health/1_run_latent_growth_curves_mh_symptoms.R
######################

library(data.table)

# ============================================================================
# USER CONFIGURATION — update all paths before running
# ============================================================================

PROJECT_DIR <- "~/your_project_directory"

# Input: LGC intercept file (output of LGC pipeline, EUR participants only)
PHENO_FILE <- file.path(PROJECT_DIR, "results/LGC/Mental_health/LGC_intercepts_GWAS_ready_EUR.csv")

# Output: plain-text list of phenotype names passing the sample size filter
# (one name per line — used as input to GWAS and GCTA batch scripts)
OUT_FILE <- file.path(PROJECT_DIR, "results/GWAS/Mental_health/phenotype_list_N3000.txt")

# Minimum sample size (non-missing participants) to include a phenotype
N_MIN <- 3000

# ============================================================================
# FILTER
# ============================================================================

pheno_file     <- fread(PHENO_FILE)
intercept_cols <- grep("_intercept$", names(pheno_file), value = TRUE)

# Count non-missing participants per phenotype
sample_sizes <- sapply(intercept_cols, function(col) sum(!is.na(pheno_file[[col]])))

phenos_for_gwas <- names(sample_sizes[sample_sizes >= N_MIN])
excluded        <- names(sample_sizes[sample_sizes <  N_MIN])

cat("Total intercept phenotypes:", length(intercept_cols), "\n")
cat("Phenotypes with N >=", N_MIN, ":", length(phenos_for_gwas), "\n")
cat("Phenotypes excluded (N < ", N_MIN, "):", length(excluded), "\n")

if (length(excluded) > 0) {
  cat("\nExcluded phenotypes:\n")
  for (p in excluded) {
    cat("  ", p, ": N =", sample_sizes[p], "\n")
  }
}

# ============================================================================
# SAVE
# ============================================================================

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
writeLines(phenos_for_gwas, OUT_FILE)
cat("\nPhenotype list saved to:", OUT_FILE, "\n")
cat("Total phenotypes for GWAS:", length(phenos_for_gwas), "\n")
