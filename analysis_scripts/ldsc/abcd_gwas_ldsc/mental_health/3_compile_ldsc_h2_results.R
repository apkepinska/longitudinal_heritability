#!/usr/bin/env Rscript
# =============================================================================
# COMPILE LDSC HERITABILITY RESULTS
# Parses all h2 log files and creates a summary table
# Author: Ada Kepinska
# Date: 2025-02-28
#
# Prerequisites:
#   - h2 log files produced by 2_run_ldsc_h2_batch.sh
# =============================================================================

library(data.table)

# ============================================================================
# USER CONFIGURATION — update all paths before running
# ============================================================================

PROJECT_DIR <- "~/your_project_directory"

# Directory containing LDSC h2 .log files (output of 2_run_ldsc_h2_batch.sh)
H2_DIR  <- file.path(PROJECT_DIR, "results/LDSC/Mental_health/h2")

# Output directory for compiled results
OUT_DIR <- file.path(PROJECT_DIR, "results/LDSC/Mental_health")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# PARSE SINGLE LOG FILE
# ============================================================================

parse_h2_log <- function(log_file) {
  lines <- readLines(log_file)
  pheno <- sub("_h2\\.log$", "", basename(log_file))

  # Check for h2 estimate line
  h2_line <- grep("^Total Observed scale h2:", lines, value = TRUE)
  if (length(h2_line) == 0) {
    return(data.table(
      phenotype    = pheno,
      h2           = NA_real_, h2_se        = NA_real_,
      lambda_gc    = NA_real_, mean_chi2     = NA_real_,
      intercept    = NA_real_, intercept_se  = NA_real_,
      ratio        = NA_character_,
      status       = "FAILED"
    ))
  }

  # Parse h2 and SE: "Total Observed scale h2: 0.0283 (0.0954)"
  h2_match     <- regmatches(h2_line, regexec("h2: ([-.0-9]+) \\(([0-9.]+)\\)", h2_line))[[1]]
  h2           <- as.numeric(h2_match[2])
  h2_se        <- as.numeric(h2_match[3])

  # Lambda GC
  lambda_line  <- grep("^Lambda GC:", lines, value = TRUE)
  lambda_gc    <- as.numeric(sub("Lambda GC: ", "", lambda_line))

  # Mean Chi^2
  chi2_line    <- grep("^Mean Chi\\^2:", lines, value = TRUE)
  mean_chi2    <- as.numeric(sub("Mean Chi\\^2: ", "", chi2_line))

  # Intercept
  int_line     <- grep("^Intercept:", lines, value = TRUE)
  int_match    <- regmatches(int_line, regexec("Intercept: ([-.0-9]+) \\(([0-9.]+)\\)", int_line))[[1]]
  intercept    <- as.numeric(int_match[2])
  intercept_se <- as.numeric(int_match[3])

  # Ratio
  ratio_line   <- grep("^Ratio", lines, value = TRUE)
  ratio        <- sub("Ratio:? ?", "", ratio_line)

  data.table(
    phenotype    = pheno,
    h2           = h2,
    h2_se        = h2_se,
    lambda_gc    = lambda_gc,
    mean_chi2    = mean_chi2,
    intercept    = intercept,
    intercept_se = intercept_se,
    ratio        = ratio,
    status       = "OK"
  )
}

# ============================================================================
# PARSE ALL LOG FILES
# ============================================================================

log_files <- list.files(H2_DIR, pattern = "_h2\\.log$", full.names = TRUE)
cat("Found", length(log_files), "h2 log files\n")

results_list <- lapply(log_files, function(f) {
  tryCatch(parse_h2_log(f), error = function(e) {
    data.table(
      phenotype    = sub("_h2\\.log$", "", basename(f)),
      h2           = NA, h2_se        = NA,
      lambda_gc    = NA, mean_chi2     = NA,
      intercept    = NA, intercept_se  = NA,
      ratio        = NA,
      status       = paste("PARSE_ERROR:", e$message)
    )
  })
})

results <- rbindlist(results_list)

# Compute Z-test p-values for h2
results[, h2_z := h2 / h2_se]
results[, h2_p := 2 * pnorm(abs(h2_z), lower.tail = FALSE)]

# Sort largest to smallest h2
results <- results[order(-h2)]

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== LDSC h2 Summary ===\n")
cat("Total phenotypes:", nrow(results), "\n")
cat("Successful:      ", sum(results$status == "OK"), "\n")
cat("Failed:          ", sum(results$status != "OK"), "\n")

cat("\nSignificant h2 (p < 0.05):\n")
sig <- results[h2_p < 0.05 & status == "OK"]
if (nrow(sig) > 0) {
  print(sig[, .(phenotype, h2, h2_se, h2_p, mean_chi2, intercept)])
} else {
  cat("  None\n")
}

cat("\nAll results (successful):\n")
print(results[status == "OK", .(phenotype, h2, h2_se, h2_p, mean_chi2, intercept)], nrows = 100)

# ============================================================================
# SAVE
# ============================================================================

fwrite(results, file.path(OUT_DIR, "LDSC_h2_all_results.csv"))
cat("\nFull results saved to:", file.path(OUT_DIR, "LDSC_h2_all_results.csv"), "\n")

# Clean version for downstream use (successful runs only)
results_clean <- results[
  status == "OK",
  .(phenotype, h2, h2_se, h2_p, lambda_gc, mean_chi2, intercept, intercept_se)
]
fwrite(results_clean, file.path(OUT_DIR, "LDSC_h2_results_clean.csv"))
cat("Clean results saved to:", file.path(OUT_DIR, "LDSC_h2_results_clean.csv"), "\n")
