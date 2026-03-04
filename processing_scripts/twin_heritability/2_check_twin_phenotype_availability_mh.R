#!/usr/bin/env Rscript
# =============================================================================
# CHECK TWIN DATA AVAILABILITY FOR MENTAL HEALTH PHENOTYPES
# =============================================================================
# Author: Ada Kepinska
# Date: 2025-02-19
#
# Purpose: For each mental health phenotype at each timepoint, count how many
#          MZ and DZ twin PAIRS have complete data (both twins measured).
#          Assess whether we have enough power for ACE heritability estimation.
#
# Power thresholds (from Chen et al. 2025, Frontiers in Genetics):
#   - For questionnaire/sum-score phenotypes (low measurement error):
#     ~300 MZ + ~300 DZ pairs minimum (traditional SEM)
#   - For trial-level cognitive tasks (high measurement error):
#     Need thousands of pairs OR 100+ trials per person (HLM approach)
#
# Prerequisites:
#   - Run 1_find_and_crosscheck_ABCD_twins.R first to generate twin pair files
#   - ABCD mental health phenotype CSVs (one per timepoint), pre-QCed
#   - Data access via: https://nda.nih.gov/abcd/
#
# Expected phenotype file structure:
#   - One CSV per timepoint, containing src_subject_id + phenotype columns
#   - Organized into subfolders by timepoint (see TIMEPOINTS section below)
#   - Files should have refused/missing responses already handled
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

cat("=============================================================\n")
cat("  TWIN DATA AVAILABILITY CHECK FOR MENTAL HEALTH PHENOTYPES\n")
cat("=============================================================\n\n")

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR <- "~/your_project_directory"   # Root of your project

# Twin pair files (output from 1_find_and_crosscheck_ABCD_twins.R)
TWINS_DIR <- file.path(PROJECT_DIR, "outputs/twin_identification")
MZ_FILE   <- file.path(TWINS_DIR, "mz_twins_final.csv")
DZ_FILE   <- file.path(TWINS_DIR, "dz_twins_final.csv")
ALL_TWINS <- file.path(TWINS_DIR, "confirmed_twins_all.csv")

# Root directory for QC'd mental health phenotype CSVs
# Expected structure:
#   <phenotype_base_path>/
#     <timepoint_folder>/
#       <phenotype_file>.csv
# Adjust folder and file naming conventions to match your pipeline's output
PHENOTYPE_BASE_PATH <- file.path(PROJECT_DIR, "data/phenotypes/mental_health")

# Output directory
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs/twin_phenotype_availability")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# POWER THRESHOLDS
# =============================================================================

POWER_MIN_PAIRS_TOTAL <- 300   # Minimum total twin pairs for any analysis
POWER_MIN_MZ_PAIRS    <- 100   # Minimum MZ pairs
POWER_MIN_DZ_PAIRS    <- 100   # Minimum DZ pairs
POWER_GOOD_PAIRS      <- 600   # "Good power" threshold (traditional recommendation)

# =============================================================================
# TIMEPOINTS
# =============================================================================
# Update folder/file names to match your local ABCD phenome directory structure.
# ABCD uses the following event name conventions:
#   baseline_year_1_arm_1
#   1_year_follow_up_y_arm_1
#   2_year_follow_up_y_arm_1   etc.

timepoints <- data.frame(
  timepoint    = c("baseline", "year1", "year2", "year3", "year4"),
  time_code    = c(0, 1, 2, 3, 4),
  pheno_folder = c(
    "baseline_year_1_arm_1",
    "1_year_follow_up_y_arm_1",
    "2_year_follow_up_y_arm_1",
    "3_year_follow_up_y_arm_1",
    "4_year_follow_up_y_arm_1"
  ),
  # Adjust file naming to match your pipeline's output conventions
  pheno_file   = c(
    "mental_health_continuous_baseline_qced.csv",
    "mental_health_continuous_year1_qced.csv",
    "mental_health_continuous_year2_qced.csv",
    "mental_health_continuous_year3_qced.csv",
    "mental_health_continuous_year4_qced.csv"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# LOAD TWIN PAIRS
# =============================================================================

cat("Loading twin pair files...\n")

mz_pairs  <- fread(MZ_FILE)
dz_pairs  <- fread(DZ_FILE)
all_twins <- fread(ALL_TWINS)

cat("  MZ pairs:", nrow(mz_pairs), "\n")
cat("  DZ pairs:", nrow(dz_pairs), "\n")
cat("  Total twin pairs (incl Unknown zyg):", nrow(all_twins), "\n\n")

mz_ids       <- unique(c(mz_pairs$id1, mz_pairs$id2))
dz_ids       <- unique(c(dz_pairs$id1, dz_pairs$id2))
all_twin_ids <- unique(c(all_twins$id1, all_twins$id2))

cat("  MZ individuals:", length(mz_ids), "\n")
cat("  DZ individuals:", length(dz_ids), "\n")
cat("  All twin individuals:", length(all_twin_ids), "\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Count pairs where BOTH twins have non-NA data for a given phenotype
count_complete_pairs <- function(pair_dt, pheno_col, pheno_data) {
  if (!pheno_col %in% names(pheno_data)) return(NA)
  valid_ids <- pheno_data[!is.na(get(pheno_col)), src_subject_id]
  sum(pair_dt$id1 %in% valid_ids & pair_dt$id2 %in% valid_ids)
}

# Count pairs where AT LEAST ONE twin has data
count_any_twin_with_data <- function(pair_dt, pheno_col, pheno_data) {
  if (!pheno_col %in% names(pheno_data)) return(NA)
  valid_ids <- pheno_data[!is.na(get(pheno_col)), src_subject_id]
  sum(pair_dt$id1 %in% valid_ids | pair_dt$id2 %in% valid_ids)
}

# =============================================================================
# MAIN: CHECK EACH PHENOTYPE AT EACH TIMEPOINT
# =============================================================================

cat("=============================================================\n")
cat("Checking phenotype availability across timepoints...\n")
cat("=============================================================\n\n")

results_list <- list()

for (tp_idx in 1:nrow(timepoints)) {
  tp <- timepoints[tp_idx, ]
  
  file_path <- file.path(PHENOTYPE_BASE_PATH, tp$pheno_folder, tp$pheno_file)
  
  if (!file.exists(file_path)) {
    cat("  WARNING: File not found for", tp$timepoint, "- skipping\n")
    cat("  Expected:", file_path, "\n\n")
    next
  }
  
  cat("Loading", tp$timepoint, "...")
  pheno_data <- fread(file_path, na.strings = c("", "NA"))
  pheno_data[, src_subject_id := as.character(src_subject_id)]
  cat(" N =", nrow(pheno_data), "\n")
  
  twins_in_tp <- sum(all_twin_ids %in% pheno_data$src_subject_id)
  cat("  Twin individuals in this timepoint:", twins_in_tp, "/", length(all_twin_ids), "\n")
  
  # Phenotype columns (exclude metadata)
  exclude_cols <- c("src_subject_id", "eventname", "interview_date", "interview_age")
  pheno_cols   <- setdiff(names(pheno_data), exclude_cols)
  
  cat("  Phenotypes to check:", length(pheno_cols), "\n")
  
  for (pheno in pheno_cols) {
    mz_complete <- count_complete_pairs(mz_pairs, pheno, pheno_data)
    dz_complete <- count_complete_pairs(dz_pairs, pheno, pheno_data)
    mz_any      <- count_any_twin_with_data(mz_pairs, pheno, pheno_data)
    dz_any      <- count_any_twin_with_data(dz_pairs, pheno, pheno_data)
    
    total_complete <- mz_complete + dz_complete
    
    power_level <- "insufficient"
    if (!is.na(total_complete)) {
      if (total_complete >= POWER_GOOD_PAIRS &
          mz_complete >= POWER_MIN_MZ_PAIRS &
          dz_complete >= POWER_MIN_DZ_PAIRS) {
        power_level <- "good"
      } else if (total_complete >= POWER_MIN_PAIRS_TOTAL &
                 mz_complete >= POWER_MIN_MZ_PAIRS &
                 dz_complete >= POWER_MIN_DZ_PAIRS) {
        power_level <- "adequate"
      } else if (total_complete >= POWER_MIN_PAIRS_TOTAL) {
        power_level <- "marginal"
      }
    }
    
    results_list[[length(results_list) + 1]] <- data.table(
      phenotype            = pheno,
      timepoint            = tp$timepoint,
      time_code            = tp$time_code,
      mz_complete_pairs    = mz_complete,
      dz_complete_pairs    = dz_complete,
      total_complete_pairs = total_complete,
      mz_any_pairs         = mz_any,
      dz_any_pairs         = dz_any,
      power_level          = power_level
    )
  }
  
  cat("\n")
}

results <- rbindlist(results_list)

# =============================================================================
# SUMMARY BY TIMEPOINT
# =============================================================================

cat("=============================================================\n")
cat("RESULTS\n")
cat("=============================================================\n\n")

cat("--- POWER SUMMARY BY TIMEPOINT ---\n\n")

tp_summary <- results[, .(
  n_phenotypes       = .N,
  n_good             = sum(power_level == "good"),
  n_adequate         = sum(power_level == "adequate"),
  n_marginal         = sum(power_level == "marginal"),
  n_insufficient     = sum(power_level == "insufficient"),
  median_mz_pairs    = as.double(median(mz_complete_pairs, na.rm = TRUE)),
  median_dz_pairs    = as.double(median(dz_complete_pairs, na.rm = TRUE)),
  median_total_pairs = as.double(median(total_complete_pairs, na.rm = TRUE))
), by = timepoint]

print(tp_summary)

# =============================================================================
# SUMMARY BY PHENOTYPE (across timepoints)
# =============================================================================

cat("\n--- PHENOTYPES WITH GOOD POWER AT ALL AVAILABLE TIMEPOINTS ---\n\n")

pheno_summary <- results[, .(
  n_timepoints      = .N,
  n_good_power      = sum(power_level == "good"),
  n_adequate_plus   = sum(power_level %in% c("good", "adequate")),
  min_mz_pairs      = min(mz_complete_pairs, na.rm = TRUE),
  min_dz_pairs      = min(dz_complete_pairs, na.rm = TRUE),
  min_total_pairs   = min(total_complete_pairs, na.rm = TRUE),
  max_total_pairs   = max(total_complete_pairs, na.rm = TRUE)
), by = phenotype]

always_good <- pheno_summary[n_good_power == n_timepoints]
cat("Phenotypes with GOOD power at ALL timepoints:", nrow(always_good), "\n")
if (nrow(always_good) > 0 & nrow(always_good) <= 50) {
  print(always_good[order(-min_total_pairs)])
}

always_adequate <- pheno_summary[n_adequate_plus == n_timepoints]
cat("\nPhenotypes with at least ADEQUATE power at ALL timepoints:", nrow(always_adequate), "\n")

never_adequate <- pheno_summary[n_adequate_plus == 0]
cat("Phenotypes that NEVER reach adequate power:", nrow(never_adequate), "\n")
if (nrow(never_adequate) > 0 & nrow(never_adequate) <= 20) {
  print(never_adequate[order(max_total_pairs)])
}

# =============================================================================
# TOP PHENOTYPES
# =============================================================================

cat("\n--- TOP 30 PHENOTYPES BY MINIMUM COMPLETE PAIRS ACROSS TIMEPOINTS ---\n\n")
top_phenos <- pheno_summary[order(-min_total_pairs)][1:min(30, nrow(pheno_summary))]
print(top_phenos)

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("\n=============================================================\n")
cat("SAVING OUTPUTS\n")
cat("=============================================================\n\n")

out1 <- file.path(OUTPUT_DIR, "twin_phenotype_availability_by_timepoint_mental_health.csv")
fwrite(results, out1)
cat("  Full results:", out1, "\n")

# Well-powered phenotypes: adequate power at 3+ timepoints
well_powered <- pheno_summary[n_adequate_plus >= 3]
out2 <- file.path(OUTPUT_DIR, "twin_well_powered_phenotypes_mental_health.csv")
fwrite(well_powered[order(-min_total_pairs)], out2)
cat("  Well-powered phenotypes:", out2, "(", nrow(well_powered), "phenotypes)\n")

cat("\nDone!\n")
