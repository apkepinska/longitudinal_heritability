#!/usr/bin/env Rscript
# =============================================================================
# REFORMAT PLINK2 GWAS OUTPUT FOR LDSC
# Converts PLINK2 .glm.linear files to LDSC-compatible format
# Author: Ada Kepinska
# Date: 2025-02-27
#
# Description:
#   Reads compressed PLINK2 GWAS output (.glm.linear.gz) for each phenotype
#   in the list and reformats to the column structure expected by LDSC munge_sumstats.
#
# Column mapping:
#   PLINK2: #CHROM, POS, ID, REF, ALT, A1, TEST, OBS_CT, BETA, SE, T_STAT, P, ERRCODE
#   LDSC:   SNP, A1, A2, BETA, SE, P, N
#
# Prerequisites:
#   - GWAS results (.glm.linear.gz) from GWAS pipeline
#     (ABCD_GWAS_for_LDSC/Mental_health/2_run_LGC_gwas_batch.sh)
#   - Phenotype list file from:
#     Latent_growth_curve/Mental_health/2_create_phenotype_list_sample_size_filter_GWAS.R
# =============================================================================

library(data.table)

# ============================================================================
# USER CONFIGURATION — update all paths before running
# ============================================================================

PROJECT_DIR <- "~/your_project_directory"

# Directory containing PLINK2 GWAS output files (.glm.linear.gz)
GWAS_DIR <- file.path(PROJECT_DIR, "results/GWAS/Mental_health")

# Phenotype list (one name per line, produced by 2_create_phenotype_list_sample_size_filter_GWAS.R)
PHENO_LIST_FILE <- file.path(GWAS_DIR, "phenotype_list_N3000.txt")

# Output directory for LDSC-formatted summary statistics
LDSC_DIR <- file.path(GWAS_DIR, "LDSC_formatted")
dir.create(LDSC_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# REFORMAT EACH GWAS FILE
# ============================================================================

pheno_list <- readLines(PHENO_LIST_FILE)
cat("Found", length(pheno_list), "phenotypes to reformat\n\n")

for (pheno in pheno_list) {

  # Locate the GWAS file — PLINK2 inserts the phenotype name before .glm.linear
  # Try exact match first, then fall back to a glob
  gwas_file <- file.path(GWAS_DIR, paste0(pheno, ".glm.linear.gz"))

  if (!file.exists(gwas_file)) {
    matches <- Sys.glob(file.path(GWAS_DIR, paste0("*", pheno, "*.glm.linear.gz")))
    if (length(matches) == 0) {
      cat("WARNING: No GWAS file found for", pheno, "- skipping\n")
      next
    }
    gwas_file <- matches[1]
  }

  cat("Processing:", pheno, "from", basename(gwas_file), "\n")

  gwas <- fread(gwas_file)

  # Reformat to LDSC column structure
  ldsc_format <- data.table(
    SNP  = gwas$ID,
    A1   = gwas$A1,
    A2   = ifelse(gwas$A1 == gwas$REF, gwas$ALT, gwas$REF),
    BETA = gwas$BETA,
    SE   = gwas$SE,
    P    = gwas$P,
    N    = gwas$OBS_CT
  )

  # Remove rows with any missing values
  n_before   <- nrow(ldsc_format)
  ldsc_format <- ldsc_format[complete.cases(ldsc_format)]
  n_removed  <- n_before - nrow(ldsc_format)
  if (n_removed > 0) cat("  Removed", n_removed, "rows with missing values\n")

  # Write and compress
  out_file <- file.path(LDSC_DIR, paste0(pheno, "_ldsc.txt"))
  fwrite(ldsc_format, out_file, sep = "\t", quote = FALSE)
  system(paste("gzip -f", shQuote(out_file)))

  cat("  Written:", nrow(ldsc_format), "SNPs to", paste0(basename(out_file), ".gz"), "\n")
}

cat("\nDone! Reformatted files are in:", LDSC_DIR, "\n")
