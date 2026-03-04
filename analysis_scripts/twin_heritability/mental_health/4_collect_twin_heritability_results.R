#!/usr/bin/env Rscript
# =============================================================================
# COLLECT TWIN HERITABILITY RESULTS
# Author: Ada Kepinska
# Date: 2025-02-28
#
# Run AFTER all SLURM array jobs complete.
# Collects individual per-phenotype result CSVs into a single summary table.
# Each result file has 4 rows: ACE, ACE_PC, HLM, HLM_PC
#
# Prerequisites:
#   - All jobs submitted by 3_submit_twin_heritability_jobs.sh have completed
# =============================================================================

library(data.table)

# =============================================================================
# USER CONFIGURATION — update path before running
# =============================================================================

PROJECT_DIR <- "~/your_project_directory"

# Directory containing per-phenotype _heritability.csv files
# (output of 1_twin_heritability_single_phenotype.R)
RESULTS_DIR <- file.path(PROJECT_DIR, "results/twin_heritability/cross_sectional")

# =============================================================================
# LOAD
# =============================================================================

cat("=============================================================\n")
cat("  COLLECTING TWIN HERITABILITY RESULTS\n")
cat("=============================================================\n\n")

result_files <- list.files(RESULTS_DIR, pattern = "_heritability\\.csv$", full.names = TRUE)
cat("Found", length(result_files), "result files\n\n")

if (length(result_files) == 0) stop("No result files found in: ", RESULTS_DIR)

all_results <- rbindlist(lapply(result_files, fread), fill = TRUE)
cat("Combined:", nrow(all_results), "rows (",
    nrow(all_results) / 4, "phenotype x timepoint combinations)\n\n")

# =============================================================================
# SUMMARIES
# =============================================================================

cat("--- STATUS BY MODEL ---\n")
print(all_results[, .N, by = .(model, status)])

cat("\n--- SUCCESSFUL RESULTS BY TIMEPOINT AND MODEL ---\n")
tp_summ <- all_results[grepl("success", status), .(
  n_phenotypes = .N,
  median_h2    = round(median(h2, na.rm = TRUE), 3),
  mean_h2      = round(mean(h2,   na.rm = TRUE), 3)
), by = .(timepoint, model)]
print(tp_summ[order(timepoint, model)])

cat("\n--- ACE vs HLM COMPARISON (base models, successful) ---\n")
ace_ok <- all_results[model == "ACE" & grepl("success", status), .(phenotype, timepoint, ace_h2 = h2)]
hlm_ok <- all_results[model == "HLM" & grepl("success", status), .(phenotype, timepoint, hlm_h2 = h2)]
both   <- merge(ace_ok, hlm_ok, by = c("phenotype", "timepoint"))
if (nrow(both) > 0) {
  cat("  N comparisons:", nrow(both), "\n")
  cat("  Correlation ACE h2 vs HLM h2:", round(cor(both$ace_h2, both$hlm_h2, use = "complete.obs"), 3), "\n")
  cat("  Mean diff (HLM - ACE):",        round(mean(both$hlm_h2 - both$ace_h2, na.rm = TRUE), 3), "\n")
}

cat("\n--- BASE vs PC SENSITIVITY ANALYSIS ---\n")
for (method in c("ACE", "HLM")) {
  base_ok <- all_results[model == method            & grepl("success", status), .(phenotype, timepoint, base_h2 = h2)]
  pc_ok   <- all_results[model == paste0(method, "_PC") & grepl("success", status), .(phenotype, timepoint, pc_h2 = h2)]
  comp    <- merge(base_ok, pc_ok, by = c("phenotype", "timepoint"))
  if (nrow(comp) > 0) {
    cat(sprintf("  %s: N=%d, corr=%.3f, mean diff (PC - base)=%.3f\n",
                method, nrow(comp),
                cor(comp$base_h2, comp$pc_h2, use = "complete.obs"),
                mean(comp$pc_h2 - comp$base_h2, na.rm = TRUE)))
  }
}

cat("\n--- TOP 20 MOST HERITABLE (HLM base, any timepoint) ---\n")
top_h2 <- all_results[model == "HLM" & grepl("success", status)][order(-h2)][1:min(20, .N)]
print(top_h2[, .(phenotype, timepoint, n_total_pairs,
                 h2 = round(h2, 3),
                 ci = paste0("[", round(h2_lo, 3), ", ", round(h2_hi, 3), "]"))])

# =============================================================================
# SAVE
# =============================================================================

out_file <- file.path(RESULTS_DIR, "all_heritability_results_combined.csv")
fwrite(all_results[order(phenotype, timepoint, model)], out_file)
cat("\nCombined results saved:   ", out_file, "\n")

summary_file <- file.path(RESULTS_DIR, "heritability_summary_successful.csv")
fwrite(all_results[grepl("success", status)][order(phenotype, timepoint, model)], summary_file)
cat("Successful results saved: ", summary_file, "\n")

cat("\nDone!\n")
