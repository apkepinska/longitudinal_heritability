#!/usr/bin/env Rscript
# =============================================================================
# TWIN HERITABILITY - SINGLE PHENOTYPE × TIMEPOINT
# Author: Ada Kepinska
# Date: 2025-02-28
#
# Purpose: Estimate heritability for ONE phenotype at ONE timepoint using:
#   1. Conventional ACE model via mets::twinlm()
#   2. Hierarchical Linear Model (HLM) via brms (Bayesian)
#
# Each method is run twice: without and with genetic PCs as covariates.
# This yields 4 estimates per phenotype × timepoint:
#   - ACE       (base covariates: zygosity + sex + age + site)
#   - ACE + PCs (base + first 10 genetic PCs)
#   - HLM       (base covariates)
#   - HLM + PCs (base + first 10 genetic PCs)
#
# Rationale for sensitivity analysis with PCs:
#   - Chen et al. (2025): did NOT include PCs (standard twin practice)
#   - Smith et al. (2023, Behav Genet): DID include 10 PCs in ABCD
#     heritability models using FEMA with full sample + GRM
#   - Our twins span multiple ancestries (identified genetically across
#     full ABCD sample), so PCs are a reasonable sensitivity analysis
#
# Called by SLURM array job. Reads phenotype name and timepoint from a
# job list file using SLURM_ARRAY_TASK_ID.
#
# References:
#   - Chen et al. (2025) Frontiers in Genetics - HLM approach
#   - Smith et al. (2023) Behavior Genetics - FEMA/ABCD cognition
#   - Code: https://github.com/afni/apaper_heritability
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mets)
  library(brms)
})

options(mc.cores = parallel::detectCores())

# =============================================================================
# PARSE ARGUMENTS FROM JOB LIST
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  job_list_file <- Sys.getenv("JOB_LIST_FILE")
  task_id       <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

  if (job_list_file == "" | is.na(task_id)) {
    stop("Usage: Rscript 1_twin_heritability_single_phenotype.R <phenotype> <timepoint>\n",
         "  OR set JOB_LIST_FILE and SLURM_ARRAY_TASK_ID environment variables")
  }

  job_list       <- fread(job_list_file, header = TRUE)
  phenotype_name <- job_list$phenotype[task_id]
  timepoint_name <- job_list$timepoint[task_id]
  power_level    <- job_list$power_level[task_id]

} else {
  phenotype_name <- args[1]
  timepoint_name <- args[2]
  power_level    <- ifelse(length(args) >= 3, args[3], "unknown")
}

cat("=============================================================\n")
cat("  TWIN HERITABILITY ESTIMATION\n")
cat("  Phenotype:", phenotype_name, "\n")
cat("  Timepoint:", timepoint_name, "\n")
cat("  Power level:", power_level, "\n")
cat("  Models: ACE, ACE+PCs, HLM, HLM+PCs\n")
cat("  Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================\n\n")

# =============================================================================
# USER CONFIGURATION — update all paths before running
# =============================================================================

PROJECT_DIR <- "~/your_project_directory"

# QC'd mental health phenotype CSVs (one subfolder per timepoint)
PHENO_BASE_PATH <- file.path(PROJECT_DIR, "data/phenotypes/mental_health")

# Twin pair files (output of Processing_scripts/Twin_heritability/1_find_and_crosscheck_ABCD_twins.R)
TWIN_PATH <- file.path(PROJECT_DIR, "outputs/twin_identification")

# Output directory for per-phenotype heritability results
OUT_DIR <- file.path(PROJECT_DIR, "results/twin_heritability/cross_sectional")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Ancestry PCs from PC-AiR (must contain src_subject_id and V1–V50)
# Set to NULL or a non-existent path to skip PC-adjusted models
PC_PATH <- file.path(PROJECT_DIR, "data/ancestry/PCaAiR_50PCs.txt")

# Cleaned demographic files (abcd_p_demo, abcd_y_lt; baseline timepoint)
DEMO_PATH <- file.path(PROJECT_DIR, "data/demographics/baseline_year_1_arm_1")

# Number of PCs to include as covariates in PC-adjusted models
N_PCS <- 10

# =============================================================================
# TIMEPOINT DEFINITIONS
# =============================================================================
# Adjust pheno_folder and pheno_file names to match your pipeline's output.

tp_map <- data.table(
  timepoint    = c("baseline", "year1", "year2", "year3", "year4"),
  pheno_folder = c(
    "baseline_year_1_arm_1",
    "1_year_follow_up_y_arm_1",
    "2_year_follow_up_y_arm_1",
    "3_year_follow_up_y_arm_1",
    "4_year_follow_up_y_arm_1"
  ),
  pheno_file   = c(
    "mental_health_continuous_baseline_qced.csv",
    "mental_health_continuous_year1_qced.csv",
    "mental_health_continuous_year2_qced.csv",
    "mental_health_continuous_year3_qced.csv",
    "mental_health_continuous_year4_qced.csv"
  )
)

tp_info <- tp_map[timepoint == timepoint_name]
if (nrow(tp_info) == 0) stop("Unknown timepoint: ", timepoint_name)

# =============================================================================
# LOAD PHENOTYPE DATA
# =============================================================================

cat("Loading phenotype data...\n")

pheno_file_path <- file.path(PHENO_BASE_PATH, tp_info$pheno_folder, tp_info$pheno_file)

if (!file.exists(pheno_file_path)) stop("Phenotype file not found: ", pheno_file_path)

ph_data <- fread(pheno_file_path)
cat("  Loaded:", nrow(ph_data), "participants\n")

if (!phenotype_name %in% names(ph_data)) {
  stop("Phenotype '", phenotype_name, "' not found in data.")
}

# =============================================================================
# LOAD TWIN DATA
# =============================================================================

cat("Loading twin data...\n")

mz_pairs <- fread(file.path(TWIN_PATH, "mz_twins_final.csv"))
dz_pairs <- fread(file.path(TWIN_PATH, "dz_twins_final.csv"))

cat("  MZ pairs:", nrow(mz_pairs), ", DZ pairs:", nrow(dz_pairs), "\n")

# Long-format twin structure
mz_long <- rbind(
  data.table(src_subject_id = mz_pairs$id1, fam = paste0("MZ_", 1:nrow(mz_pairs)), zyg = "MZ", num = 1),
  data.table(src_subject_id = mz_pairs$id2, fam = paste0("MZ_", 1:nrow(mz_pairs)), zyg = "MZ", num = 2)
)

dz_long <- rbind(
  data.table(src_subject_id = dz_pairs$id1, fam = paste0("DZ_", 1:nrow(dz_pairs)), zyg = "DZ", num = 1),
  data.table(src_subject_id = dz_pairs$id2, fam = paste0("DZ_", 1:nrow(dz_pairs)), zyg = "DZ", num = 2)
)

twin_structure     <- rbind(mz_long, dz_long)
twin_structure[, MZ := ifelse(zyg == "MZ", 1, 0)]
twin_structure[, DZ := 1 - MZ]

# =============================================================================
# LOAD DEMOGRAPHICS AND GENETIC PCs
# =============================================================================

cat("Loading demographics...\n")

# Sex from abcd_p_demo (baseline); must contain: src_subject_id, demo_sex_v2
abcd_p_demo <- fread(file.path(DEMO_PATH, "abcd_p_demo_baseline_cleaned.csv"))
demo_sex    <- abcd_p_demo[, .(src_subject_id, sex = demo_sex_v2)]

# Age and site from abcd_y_lt (baseline); must contain: src_subject_id, interview_age, site_id_l
abcd_y_lt      <- fread(file.path(DEMO_PATH, "abcd_y_lt_baseline_cleaned.csv"))
demo_age_site  <- abcd_y_lt[, .(src_subject_id, age = interview_age, site = site_id_l)]

# Genetic PCs
pc_cols        <- paste0("PC", 1:N_PCS)
pcs_available  <- FALSE

if (!is.null(PC_PATH) && file.exists(PC_PATH)) {
  cat("Loading genetic PCs...\n")
  pc_data <- fread(PC_PATH)
  v_cols  <- paste0("V", 1:N_PCS)
  if (all(v_cols %in% names(pc_data))) {
    pc_data <- pc_data[, c("src_subject_id", v_cols), with = FALSE]
    setnames(pc_data, v_cols, pc_cols)
  } else {
    pc_data <- pc_data[, c(1, 2:(N_PCS + 1)), with = FALSE]
    setnames(pc_data, names(pc_data), c("src_subject_id", pc_cols))
  }
  pcs_available <- TRUE
  cat("  Loaded", N_PCS, "PCs for", nrow(pc_data), "individuals\n")
} else {
  cat("  WARNING: PC file not found or not set. Skipping PC-adjusted models.\n")
}

# =============================================================================
# MERGE DATA
# =============================================================================

cat("Merging data...\n")

twin_pheno <- merge(
  twin_structure,
  ph_data[, .(src_subject_id, pheno_value = get(phenotype_name))],
  by = "src_subject_id"
)
twin_pheno <- merge(twin_pheno, demo_sex,      by = "src_subject_id", all.x = TRUE)
twin_pheno <- merge(twin_pheno, demo_age_site, by = "src_subject_id", all.x = TRUE)

if (pcs_available) {
  twin_pheno <- merge(twin_pheno, pc_data, by = "src_subject_id", all.x = TRUE)
}

# =============================================================================
# FILTER TO COMPLETE PAIRS
# =============================================================================

cat("Filtering to complete pairs...\n")

twin_pheno <- twin_pheno[sex %in% c(1, 2)]
twin_pheno[, sex := factor(sex, levels = c(1, 2), labels = c("male", "female"))]

# Base model: complete on pheno + sex + age + site
complete_obs       <- twin_pheno[!is.na(pheno_value) & !is.na(sex) & !is.na(age) & !is.na(site)]
pair_counts        <- complete_obs[, .N, by = fam]
complete_pairs_base <- pair_counts[N == 2, fam]

# PC model: also complete on all PCs
if (pcs_available) {
  pc_na_check        <- rowSums(is.na(twin_pheno[, ..pc_cols])) == 0
  complete_obs_pc    <- twin_pheno[!is.na(pheno_value) & !is.na(sex) & !is.na(age) & !is.na(site) & pc_na_check]
  pair_counts_pc     <- complete_obs_pc[, .N, by = fam]
  complete_pairs_pc  <- pair_counts_pc[N == 2, fam]
} else {
  complete_pairs_pc  <- character(0)
}

twin_base <- twin_pheno[fam %in% complete_pairs_base]
twin_pc   <- if (pcs_available) twin_pheno[fam %in% complete_pairs_pc] else data.table()

n_mz_base    <- sum(grepl("^MZ_", complete_pairs_base))
n_dz_base    <- sum(grepl("^DZ_", complete_pairs_base))
n_total_base <- n_mz_base + n_dz_base

n_mz_pc      <- sum(grepl("^MZ_", complete_pairs_pc))
n_dz_pc      <- sum(grepl("^DZ_", complete_pairs_pc))
n_total_pc   <- n_mz_pc + n_dz_pc

cat("  Base model pairs:", n_total_base, "(MZ:", n_mz_base, ", DZ:", n_dz_base, ")\n")
if (pcs_available) cat("  PC model pairs:  ", n_total_pc,   "(MZ:", n_mz_pc,   ", DZ:", n_dz_pc,   ")\n")

# Helper: write an empty results file and exit cleanly
write_empty_and_quit <- function(reason) {
  empty_row <- data.table(
    phenotype     = phenotype_name, timepoint    = timepoint_name,
    power_level   = power_level,   model         = NA_character_,
    n_mz_pairs    = NA_integer_,   n_dz_pairs    = NA_integer_,
    n_total_pairs = NA_integer_,
    h2 = NA_real_, c2 = NA_real_, e2 = NA_real_,
    h2_se = NA_real_, c2_se = NA_real_, e2_se = NA_real_,
    h2_lo = NA_real_, h2_hi = NA_real_, c2_lo = NA_real_, c2_hi = NA_real_,
    e2_lo = NA_real_, e2_hi = NA_real_, rMZ = NA_real_, rDZ = NA_real_,
    status = reason
  )
  result <- rbindlist(lapply(c("ACE", "ACE_PC", "HLM", "HLM_PC"), function(m) {
    r <- copy(empty_row); r$model <- m; r
  }))
  fwrite(result, file.path(OUT_DIR, paste0(phenotype_name, "_", timepoint_name, "_heritability.csv")))
  quit(save = "no", status = 0)
}

if (n_total_base < 50) {
  cat("\nERROR: Too few complete pairs (", n_total_base, "). Skipping.\n")
  write_empty_and_quit("too_few_pairs")
}

pheno_var <- var(twin_base$pheno_value, na.rm = TRUE)
if (pheno_var < 1e-10) {
  cat("\nERROR: Near-zero variance. Skipping.\n")
  write_empty_and_quit("zero_variance")
}

pheno_min    <- min(twin_base$pheno_value, na.rm = TRUE)
use_lognormal <- pheno_min > 0
if (!use_lognormal) cat("Note: Phenotype has zero/negative values (min =", pheno_min, "). Using Gaussian.\n")

cat("\nPhenotype summary:\n"); print(summary(twin_base$pheno_value))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

run_ace <- function(dat, include_pcs = FALSE) {
  if (include_pcs && pcs_available) {
    pc_terms <- paste(pc_cols, collapse = " + ")
    f <- as.formula(paste0("pheno_value ~ zyg + sex + age + factor(site) + ", pc_terms))
  } else {
    f <- pheno_value ~ zyg + sex + age + factor(site)
  }

  out <- list(h2 = NA_real_, c2 = NA_real_, e2 = NA_real_,
              h2_se = NA_real_, c2_se = NA_real_, e2_se = NA_real_,
              status = "success")

  tryCatch({
    m <- twinlm(f, data = dat, DZ = "DZ", zyg = "zyg", id = "fam", type = "ace")
    s <- summary(m)
    out$h2 <- s$acde[1, 1]; out$c2 <- s$acde[2, 1]; out$e2 <- s$acde[3, 1]
    if (ncol(s$acde) >= 2) {
      out$h2_se <- s$acde[1, 2]; out$c2_se <- s$acde[2, 2]; out$e2_se <- s$acde[3, 2]
    }
  }, error = function(e) { out$status <<- paste0("error: ", conditionMessage(e)) })

  return(out)
}

run_hlm <- function(dat, include_pcs = FALSE) {
  if (include_pcs && pcs_available) {
    pc_terms  <- paste(pc_cols, collapse = " + ")
    fixed_str <- paste0("pheno_value ~ MZ + sex + s(age, by = sex) + ", pc_terms,
                        " + (0 + DZ | fam) + (0 + MZ | fam)")
  } else {
    fixed_str <- "pheno_value ~ MZ + sex + s(age, by = sex) + (0 + DZ | fam) + (0 + MZ | fam)"
  }

  brm_family <- if (use_lognormal) lognormal() else gaussian()

  out <- list(h2_median = NA_real_, c2_median = NA_real_, e2_median = NA_real_,
              h2_lo = NA_real_, h2_hi = NA_real_,
              c2_lo = NA_real_, c2_hi = NA_real_,
              e2_lo = NA_real_, e2_hi = NA_real_,
              rMZ = NA_real_, rDZ = NA_real_,
              n_divergent = NA_integer_,
              status = "success")

  tryCatch({
    m <- brm(
      bf(as.formula(fixed_str), sigma ~ 0 + zyg),
      data    = dat, family  = brm_family,
      chains  = 4,   iter    = 4000, warmup = 2000,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      silent  = 2,   refresh = 0
    )

    n_div           <- sum(subset(nuts_params(m), Parameter == "divergent__")$Value)
    out$n_divergent <- n_div
    if (n_div > 0)   cat("  WARNING:", n_div, "divergent transitions\n")
    if (n_div > 100) out$status <- paste0("success_with_", n_div, "_divergent")

    posterior <- as_draws_df(m)
    var_MZ    <- posterior$sd_fam__MZ^2
    var_DZ    <- posterior$sd_fam__DZ^2
    sigma_MZ  <- exp(posterior$b_sigma_zygMZ)^2
    sigma_DZ  <- exp(posterior$b_sigma_zygDZ)^2

    r_MZ <- var_MZ / (var_MZ + sigma_MZ)
    r_DZ <- var_DZ / (var_DZ + sigma_DZ)
    h2   <- 2 * (r_MZ - r_DZ)
    c2   <- 2 * r_DZ - r_MZ
    e2   <- 1 - r_MZ

    out$h2_median <- median(h2); out$h2_lo <- quantile(h2, 0.025); out$h2_hi <- quantile(h2, 0.975)
    out$c2_median <- median(c2); out$c2_lo <- quantile(c2, 0.025); out$c2_hi <- quantile(c2, 0.975)
    out$e2_median <- median(e2); out$e2_lo <- quantile(e2, 0.025); out$e2_hi <- quantile(e2, 0.975)
    out$rMZ       <- median(r_MZ)
    out$rDZ       <- median(r_DZ)
  }, error = function(e) { out$status <<- paste0("error: ", conditionMessage(e)) })

  return(out)
}

make_result_row <- function(model_name, ace_out = NULL, hlm_out = NULL, n_mz, n_dz, n_total) {
  if (!is.null(ace_out)) {
    return(data.table(
      model = model_name, n_mz_pairs = n_mz, n_dz_pairs = n_dz, n_total_pairs = n_total,
      h2 = ace_out$h2, c2 = ace_out$c2, e2 = ace_out$e2,
      h2_se = ace_out$h2_se, c2_se = ace_out$c2_se, e2_se = ace_out$e2_se,
      h2_lo = NA_real_, h2_hi = NA_real_, c2_lo = NA_real_, c2_hi = NA_real_,
      e2_lo = NA_real_, e2_hi = NA_real_, rMZ = NA_real_, rDZ = NA_real_,
      status = ace_out$status))
  } else {
    return(data.table(
      model = model_name, n_mz_pairs = n_mz, n_dz_pairs = n_dz, n_total_pairs = n_total,
      h2 = hlm_out$h2_median, c2 = hlm_out$c2_median, e2 = hlm_out$e2_median,
      h2_se = NA_real_, c2_se = NA_real_, e2_se = NA_real_,
      h2_lo = hlm_out$h2_lo, h2_hi = hlm_out$h2_hi,
      c2_lo = hlm_out$c2_lo, c2_hi = hlm_out$c2_hi,
      e2_lo = hlm_out$e2_lo, e2_hi = hlm_out$e2_hi,
      rMZ = hlm_out$rMZ, rDZ = hlm_out$rDZ,
      n_divergent = hlm_out$n_divergent,
      status = hlm_out$status))
  }
}

make_skipped_row <- function(model_name, n_mz, n_dz, n_total, reason) {
  data.table(
    model = model_name, n_mz_pairs = n_mz, n_dz_pairs = n_dz, n_total_pairs = n_total,
    h2 = NA_real_, c2 = NA_real_, e2 = NA_real_,
    h2_se = NA_real_, c2_se = NA_real_, e2_se = NA_real_,
    h2_lo = NA_real_, h2_hi = NA_real_, c2_lo = NA_real_, c2_hi = NA_real_,
    e2_lo = NA_real_, e2_hi = NA_real_, rMZ = NA_real_, rDZ = NA_real_,
    status = reason)
}

run_pc_ok <- pcs_available && n_total_pc >= 50

# =============================================================================
# RUN ALL 4 MODELS
# =============================================================================

results_list <- list()

# --- Model 1: ACE (base) ---
cat("\n=== MODEL 1/4: ACE (base) ===\n")
ace1 <- run_ace(twin_base, include_pcs = FALSE)
cat("  h2=", round(ace1$h2, 3), " c2=", round(ace1$c2, 3), " e2=", round(ace1$e2, 3), " [", ace1$status, "]\n")
results_list[[1]] <- make_result_row("ACE",    ace_out = ace1, n_mz = n_mz_base, n_dz = n_dz_base, n_total = n_total_base)

# --- Model 2: ACE + PCs ---
cat("\n=== MODEL 2/4: ACE + PCs ===\n")
if (run_pc_ok) {
  ace2 <- run_ace(twin_pc, include_pcs = TRUE)
  cat("  h2=", round(ace2$h2, 3), " c2=", round(ace2$c2, 3), " e2=", round(ace2$e2, 3), " [", ace2$status, "]\n")
  results_list[[2]] <- make_result_row("ACE_PC", ace_out = ace2, n_mz = n_mz_pc, n_dz = n_dz_pc, n_total = n_total_pc)
} else {
  cat("  Skipped (PCs unavailable or too few pairs:", n_total_pc, ")\n")
  results_list[[2]] <- make_skipped_row("ACE_PC", n_mz_pc, n_dz_pc, n_total_pc, "skipped_no_pcs")
}

# --- Model 3: HLM (base) ---
cat("\n=== MODEL 3/4: HLM (base) ===\n")
hlm1 <- run_hlm(twin_base, include_pcs = FALSE)
cat("  h2=", round(hlm1$h2_median, 3), " c2=", round(hlm1$c2_median, 3), " e2=", round(hlm1$e2_median, 3), " [", hlm1$status, "]\n")
results_list[[3]] <- make_result_row("HLM",    hlm_out = hlm1, n_mz = n_mz_base, n_dz = n_dz_base, n_total = n_total_base)

# --- Model 4: HLM + PCs ---
cat("\n=== MODEL 4/4: HLM + PCs ===\n")
if (run_pc_ok) {
  hlm2 <- run_hlm(twin_pc, include_pcs = TRUE)
  cat("  h2=", round(hlm2$h2_median, 3), " c2=", round(hlm2$c2_median, 3), " e2=", round(hlm2$e2_median, 3), " [", hlm2$status, "]\n")
  results_list[[4]] <- make_result_row("HLM_PC", hlm_out = hlm2, n_mz = n_mz_pc, n_dz = n_dz_pc, n_total = n_total_pc)
} else {
  cat("  Skipped (PCs unavailable or too few pairs:", n_total_pc, ")\n")
  results_list[[4]] <- make_skipped_row("HLM_PC", n_mz_pc, n_dz_pc, n_total_pc, "skipped_no_pcs")
}

# =============================================================================
# COMPILE AND SAVE
# =============================================================================

result <- rbindlist(results_list)
result[, phenotype   := phenotype_name]
result[, timepoint   := timepoint_name]
result[, power_level := power_level]
setcolorder(result, c("phenotype", "timepoint", "power_level", "model"))

out_file <- file.path(OUT_DIR, paste0(phenotype_name, "_", timepoint_name, "_heritability.csv"))
fwrite(result, out_file)
cat("\nResults saved:", out_file, "\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=============================================================\n")
cat("  COMPARISON: ALL MODELS\n")
cat("  Phenotype:", phenotype_name, " | Timepoint:", timepoint_name, "\n")
cat("=============================================================\n\n")

cat(sprintf("%-10s  %7s  %7s  %7s  %s\n", "Model", "h2", "c2", "e2", "Status"))
cat(sprintf("%-10s  %7s  %7s  %7s  %s\n", "-----", "---", "---", "---", "------"))
for (i in 1:nrow(result)) {
  r <- result[i]
  if (!is.na(r$h2)) {
    cat(sprintf("%-10s  %6.1f%%  %6.1f%%  %6.1f%%  %s\n",
                r$model, r$h2 * 100, r$c2 * 100, r$e2 * 100, r$status))
  } else {
    cat(sprintf("%-10s  %7s  %7s  %7s  %s\n", r$model, "NA", "NA", "NA", r$status))
  }
}

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
