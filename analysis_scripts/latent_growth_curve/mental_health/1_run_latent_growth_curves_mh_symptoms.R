######################
# Latent Growth Curve Modeling for ABCD mental health phenotypes
# EUR-ONLY VERSION: Fits LGC models on EUR participants only
# Extracts intercept (and slope) factor scores for GWAS/heritability
# Based on Deng et al. (2024) Molecular Psychiatry methodology
# https://www.nature.com/articles/s41380-024-02704-4
# Author: Ada Kepinska
# Date: 2025-02-16
#
# Prerequisites:
#   - QC'd mental health phenotype CSVs across 5 timepoints (ABCD release 5.1)
#   - EUR ancestry IDs derived from a PRS or ancestry inference pipeline
#   - Unrelated individual list from PC-AiR (GENESIS)
#   - Ancestry PCs (50 PCs from PC-AiR)
#   - Cleaned demographic files: abcd_p_demo, abcd_y_lt (baseline)
#   - ABCD data access: https://nda.nih.gov/abcd/
#
# Outputs (all saved to OUTPUT_DIR):
#   1. LGC_factor_scores_mental_health_EUR.csv       — all factor scores
#   2. LGC_phenotypes_GWAS_ready_EUR.csv             — all factor scores + covariates + PCs
#   3. LGC_intercepts_GWAS_ready_EUR.csv             — PRIMARY GWAS file (intercepts only)
#   4. ABCD_mental_health_phenotype_availability_EUR.csv
######################

library(data.table)
library(dplyr)
library(tidyverse)
library(lavaan)      # For LGC models
library(broom)

# ============================================================================
# USER CONFIGURATION — update all paths before running
# ============================================================================

PROJECT_DIR <- "~/your_project_directory"

# Root directory for QC'd mental health phenotype CSVs (one subfolder per timepoint)
PHENOTYPE_BASE_PATH <- file.path(PROJECT_DIR, "data/phenotypes/mental_health")

# Root directory for cleaned demographic files (abcd_p_demo, abcd_y_lt)
# Expected: <DEMO_PATH>/baseline_year_1_arm_1/abcd_p_demo_baseline_cleaned.csv
#           <DEMO_PATH>/baseline_year_1_arm_1/abcd_y_lt_baseline_cleaned.csv
DEMO_PATH <- file.path(PROJECT_DIR, "data/demographics")

# EUR ancestry IDs — one column named src_subject_id (or IID, renamed below)
# Typically derived from a PRS file or ancestry inference output for EUR participants
EUR_IDS_FILE <- file.path(PROJECT_DIR, "data/ancestry/EUR_sample_ids.txt")

# Unrelated individuals list from PC-AiR (GENESIS pipeline), no header, one ID per row
UNRELATEDS_FILE <- file.path(PROJECT_DIR, "data/ancestry/PCaAiR_unrelateds.txt")

# Ancestry PCs from PC-AiR (50 PCs); must contain src_subject_id and V1-V10 (at minimum)
PCS_FILE <- file.path(PROJECT_DIR, "data/ancestry/PCaAiR_50PCs.txt")

# Output directory
OUTPUT_DIR <- file.path(PROJECT_DIR, "results/LGC/Mental_health")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ============================================================================
# TIMEPOINT DEFINITIONS
# ============================================================================
# Adjust pheno_folder and pheno_file names to match your pipeline's output.
# ABCD uses the following event name conventions:
#   baseline_year_1_arm_1, 1_year_follow_up_y_arm_1, etc.

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
  pheno_file   = c(
    "mental_health_continuous_baseline_qced.csv",
    "mental_health_continuous_year1_qced.csv",
    "mental_health_continuous_year2_qced.csv",
    "mental_health_continuous_year3_qced.csv",
    "mental_health_continuous_year4_qced.csv"
  ),
  stringsAsFactors = FALSE
)

# ============================================================================
# PART 1: LOAD EUR SAMPLE IDs
# ============================================================================

cat("=== Loading EUR sample identifiers ===\n")

# EUR IDs (from ancestry inference or PRS file)
eur_raw <- fread(EUR_IDS_FILE)
# Rename IID to src_subject_id if needed (common in plink-format files)
if ("IID" %in% names(eur_raw) && !"src_subject_id" %in% names(eur_raw)) {
  names(eur_raw)[names(eur_raw) == "IID"] <- "src_subject_id"
}
EUR_ids <- unique(eur_raw$src_subject_id)
cat("EUR individuals identified:", length(EUR_ids), "\n")

# Unrelated individuals (PC-AiR output)
unrelateds <- fread(UNRELATEDS_FILE, header = FALSE)
names(unrelateds) <- "src_subject_id"
cat("Unrelated individuals:", nrow(unrelateds), "\n")

# Unrelated EUR IDs (intersection)
EUR_unrelated_ids <- intersect(EUR_ids, unrelateds$src_subject_id)
cat("Unrelated EUR individuals:", length(EUR_unrelated_ids), "\n")

# ============================================================================
# PART 2: LOAD AND MERGE PHENOTYPE DATA ACROSS TIMEPOINTS (EUR ONLY)
# ============================================================================

load_phenotype_timepoint <- function(tp_row, base_path, keep_ids) {
  file_path <- file.path(base_path, tp_row$pheno_folder, tp_row$pheno_file)

  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }

  cat("Loading:", tp_row$timepoint, "...")
  dat <- fread(file_path, header = TRUE, stringsAsFactors = FALSE)

  # Filter to EUR unrelated IDs immediately
  dat <- dat[dat$src_subject_id %in% keep_ids, ]

  dat$timepoint  <- tp_row$timepoint
  dat$time_code  <- tp_row$time_code

  cat(" N =", nrow(dat), "(EUR unrelated)\n")
  return(dat)
}

cat("\n=== Loading phenotype data across all timepoints (EUR only) ===\n")
pheno_list <- list()
for (i in 1:nrow(timepoints)) {
  pheno_list[[i]] <- load_phenotype_timepoint(timepoints[i, ], PHENOTYPE_BASE_PATH, EUR_unrelated_ids)
}

pheno_list <- pheno_list[!sapply(pheno_list, is.null)]

if (length(pheno_list) < 3) {
  stop("Need at least 3 timepoints for LGC modeling!")
}

cat("\nSuccessfully loaded", length(pheno_list), "timepoints for EUR sample\n")

# ============================================================================
# PART 3: IDENTIFY PHENOTYPES AVAILABLE ACROSS TIMEPOINTS
# ============================================================================

get_phenotype_cols <- function(dat) {
  exclude_cols <- c("src_subject_id", "timepoint", "time_code",
                    "eventname", "interview_date", "interview_age")
  cols <- names(dat)
  cols[!cols %in% exclude_cols]
}

pheno_cols_per_tp  <- lapply(pheno_list, get_phenotype_cols)

all_phenos         <- unique(unlist(pheno_cols_per_tp))
pheno_availability <- sapply(all_phenos, function(p) {
  sum(sapply(pheno_cols_per_tp, function(cols) p %in% cols))
})

phenos_for_lgc <- names(pheno_availability[pheno_availability >= 3])
cat("\nPhenotypes available at 3+ timepoints:", length(phenos_for_lgc), "\n")

# Save phenotype availability
pheno_summary <- data.frame(
  phenotype    = names(pheno_availability),
  n_timepoints = as.numeric(pheno_availability)
)
write.csv(pheno_summary,
          file.path(OUTPUT_DIR, "ABCD_mental_health_phenotype_availability_EUR.csv"),
          row.names = FALSE)

# ============================================================================
# PART 4: LOAD COVARIATES (FILTERED TO EUR)
# ============================================================================

cat("\n=== Loading covariates (EUR only) ===\n")

# Sex from baseline demographics (abcd_p_demo)
# Must contain: src_subject_id, demo_sex_v2
abcd_p_demo <- fread(file.path(DEMO_PATH, "baseline_year_1_arm_1/abcd_p_demo_baseline_cleaned.csv"))
demo_sex    <- abcd_p_demo[, c("src_subject_id", "demo_sex_v2")]

# Age at baseline from longitudinal tracking (abcd_y_lt)
# Must contain: src_subject_id, interview_age
abcd_y_lt_baseline <- fread(file.path(DEMO_PATH, "baseline_year_1_arm_1/abcd_y_lt_baseline_cleaned.csv"))
age_baseline        <- abcd_y_lt_baseline[, c("src_subject_id", "interview_age")]
names(age_baseline)[2] <- "age_baseline"

covariates <- merge(demo_sex, age_baseline, by = "src_subject_id", all = TRUE)
covariates <- covariates[demo_sex_v2 %in% c(1, 2), ]  # Keep coded sex values only

# Filter to EUR unrelated
covariates <- covariates[covariates$src_subject_id %in% EUR_unrelated_ids, ]

cat("Covariates prepared for", nrow(covariates), "EUR unrelated participants\n")

# ============================================================================
# PART 5: HELPER FUNCTION TO CREATE WIDE FORMAT DATA
# ============================================================================

create_wide_data <- function(pheno_name, pheno_list) {
  extracted <- lapply(seq_along(pheno_list), function(i) {
    dat <- pheno_list[[i]]
    if (pheno_name %in% names(dat)) {
      out <- dat[, c("src_subject_id", pheno_name), with = FALSE]
      names(out)[2] <- paste0(pheno_name, "_t", i - 1)
      return(out)
    } else {
      return(NULL)
    }
  })

  extracted <- extracted[!sapply(extracted, is.null)]
  if (length(extracted) < 3) return(NULL)

  wide_dat <- extracted[[1]]
  for (j in 2:length(extracted)) {
    wide_dat <- merge(wide_dat, extracted[[j]], by = "src_subject_id", all = TRUE)
  }

  return(wide_dat)
}

# ============================================================================
# PART 6: MODEL QUALITY CHECK FUNCTIONS
# ============================================================================

check_converged <- function(fit) {
  if (is.null(fit)) return(FALSE)

  tryCatch({
    converged   <- lavInspect(fit, "converged")
    if (!converged) return(FALSE)
    post_check  <- lavInspect(fit, "post.check")
    return(converged && post_check)
  }, error = function(e) {
    return(FALSE)
  })
}

check_heywood <- function(fit) {
  if (is.null(fit)) {
    return(list(has_heywood = NA, negative_variances = character(0), boundary = FALSE))
  }

  if (!check_converged(fit)) {
    return(list(has_heywood = NA, negative_variances = character(0), boundary = FALSE,
                note = "model did not converge"))
  }

  tryCatch({
    params    <- parameterEstimates(fit)
    variances <- params[params$op == "~~" & params$lhs == params$rhs, ]
    lv_names  <- lavNames(fit, "lv")
    lv_variances <- variances[variances$lhs %in% lv_names, ]

    neg_lv_var   <- any(lv_variances$est < 0, na.rm = TRUE)
    boundary_var <- any(lv_variances$est < 0.001 & lv_variances$est >= 0, na.rm = TRUE)

    if (neg_lv_var) {
      neg_vars <- lv_variances[lv_variances$est < 0, "lhs"]
      return(list(has_heywood = TRUE, negative_variances = neg_vars, boundary = FALSE))
    } else if (boundary_var) {
      bound_vars <- lv_variances[lv_variances$est < 0.001 & lv_variances$est >= 0, "lhs"]
      return(list(has_heywood = FALSE, negative_variances = character(0),
                  boundary = TRUE, boundary_variances = bound_vars))
    } else {
      return(list(has_heywood = FALSE, negative_variances = character(0), boundary = FALSE))
    }
  }, error = function(e) {
    return(list(has_heywood = NA, negative_variances = character(0), boundary = FALSE,
                note = paste("error:", e$message)))
  })
}

has_heywood_safe <- function(heywood_result) {
  if (is.null(heywood_result))             return(NA)
  if (!is.list(heywood_result))            return(NA)
  if (is.null(heywood_result$has_heywood)) return(NA)
  return(isTRUE(heywood_result$has_heywood))
}

# ============================================================================
# PART 6B: COMPREHENSIVE MODEL EXTRACTION FUNCTION
# Extracts fit indices AND parameter estimates from a single model.
# Used to produce Table 1-style output (cf. Deng et al. 2024, Table 1).
# ============================================================================

extract_model_info <- function(fit, model_name) {
  result <- list(
    model        = model_name,
    converged    = FALSE,
    # Fit indices
    aic = NA, bic = NA, cfi = NA, tli = NA, rmsea = NA, srmr = NA,
    # Intercept
    intercept_est = NA, intercept_se = NA, intercept_pval = NA, intercept_var = NA,
    # Slope
    slope_est = NA, slope_se = NA, slope_pval = NA, slope_var = NA,
    # Quadratic
    quad_est = NA, quad_se = NA, quad_pval = NA, quad_var = NA,
    # Covariances
    int_slope_cov = NA, int_slope_cov_pval = NA,
    int_quad_cov  = NA, slope_quad_cov = NA
  )

  if (is.null(fit)) return(result)

  converged <- tryCatch({
    lavInspect(fit, "converged") && lavInspect(fit, "post.check")
  }, error = function(e) FALSE)

  result$converged <- converged
  if (!converged) return(result)

  fit_measures <- tryCatch({
    fitMeasures(fit, c("aic", "bic", "cfi", "tli", "rmsea", "srmr"))
  }, error = function(e) {
    c(aic = NA, bic = NA, cfi = NA, tli = NA, rmsea = NA, srmr = NA)
  })

  result$aic   <- as.numeric(fit_measures["aic"])
  result$bic   <- as.numeric(fit_measures["bic"])
  result$cfi   <- as.numeric(fit_measures["cfi"])
  result$tli   <- as.numeric(fit_measures["tli"])
  result$rmsea <- as.numeric(fit_measures["rmsea"])
  result$srmr  <- as.numeric(fit_measures["srmr"])

  params <- tryCatch(parameterEstimates(fit), error = function(e) NULL)

  if (!is.null(params)) {
    get_param <- function(lhs_val, op_val, rhs_val = NULL) {
      if (is.null(rhs_val)) {
        row <- params[params$lhs == lhs_val & params$op == op_val, ]
      } else {
        row <- params[params$lhs == lhs_val & params$op == op_val & params$rhs == rhs_val, ]
      }
      if (nrow(row) == 0) return(list(est = NA, se = NA, pvalue = NA))
      return(list(est = row$est[1], se = row$se[1], pvalue = row$pvalue[1]))
    }

    # Latent means
    int_mean   <- get_param("intercept", "~1")
    result$intercept_est  <- int_mean$est
    result$intercept_se   <- int_mean$se
    result$intercept_pval <- int_mean$pvalue

    slope_mean  <- get_param("slope", "~1")
    result$slope_est  <- slope_mean$est
    result$slope_se   <- slope_mean$se
    result$slope_pval <- slope_mean$pvalue

    quad_mean   <- get_param("quad", "~1")
    result$quad_est  <- quad_mean$est
    result$quad_se   <- quad_mean$se
    result$quad_pval <- quad_mean$pvalue

    # Latent variances
    result$intercept_var <- get_param("intercept", "~~", "intercept")$est
    result$slope_var     <- get_param("slope",     "~~", "slope")$est
    result$quad_var      <- get_param("quad",      "~~", "quad")$est

    # Covariances
    int_slope <- get_param("intercept", "~~", "slope")
    result$int_slope_cov      <- int_slope$est
    result$int_slope_cov_pval <- int_slope$pvalue
    result$int_quad_cov       <- get_param("intercept", "~~", "quad")$est
    result$slope_quad_cov     <- get_param("slope",     "~~", "quad")$est
  }

  return(result)
}

# ============================================================================
# PART 7: MAIN LGC MODEL FITTING FUNCTION
# ============================================================================

fit_lgc_model <- function(pheno_name, pheno_list, covariates, test_quadratic = TRUE) {

  cat("\n--- Fitting LGC for:", pheno_name, "---\n")

  skip_result <- list(
    phenotype            = pheno_name,
    status               = "skipped",
    skip_reason          = NULL,
    n_timepoints_initial = NA,
    n_timepoints_valid   = NA,
    n_participants       = NA,
    empty_timepoints     = NULL,
    selected_model       = NULL,
    factor_scores        = NULL,
    fit_indices          = NULL,
    heywood_info         = NULL,
    convergence          = NULL,
    selection_reason     = NULL,
    model_info_ng        = NULL,
    model_info_lin       = NULL,
    model_info_quad      = NULL
  )

  wide_dat <- create_wide_data(pheno_name, pheno_list)

  if (is.null(wide_dat)) {
    cat("  Skipping - insufficient timepoints in source data\n")
    skip_result$skip_reason <- "insufficient timepoints in source data"
    return(skip_result)
  }

  pheno_escaped <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", pheno_name)
  tp_cols <- grep(paste0("^", pheno_escaped, "_t"), names(wide_dat), value = TRUE)
  n_tp    <- length(tp_cols)
  skip_result$n_timepoints_initial <- n_tp

  if (n_tp < 3) {
    cat("  Skipping - only", n_tp, "timepoints available\n")
    skip_result$skip_reason <- paste0("only ", n_tp, " timepoints available")
    return(skip_result)
  }

  cat("  Timepoints available:", n_tp, "\n")

  # Merge with covariates
  analysis_dat <- merge(wide_dat, covariates, by = "src_subject_id", all.x = TRUE)
  analysis_dat <- analysis_dat[!is.na(analysis_dat$age_baseline), ]

  # Remove rows with no phenotype data at any timepoint
  has_any_data <- rowSums(!is.na(analysis_dat[, ..tp_cols])) > 0
  analysis_dat <- analysis_dat[has_any_data, ]

  for (col in tp_cols) {
    analysis_dat[[col]] <- as.numeric(analysis_dat[[col]])
  }
  analysis_dat$age_baseline <- as.numeric(analysis_dat$age_baseline)

  # --- Timepoint diagnostics ---
  tp_diagnostics <- lapply(tp_cols, function(col) {
    vals       <- analysis_dat[[col]]
    n_total    <- length(vals)
    n_valid    <- sum(!is.na(vals))
    n_missing  <- sum(is.na(vals))
    val_var    <- var(vals, na.rm = TRUE)
    has_var    <- !is.na(val_var) && val_var > 0
    is_usable  <- n_valid >= 50 && has_var

    list(
      column      = col,
      n_total     = n_total,
      n_valid     = n_valid,
      n_missing   = n_missing,
      pct_missing = round(100 * n_missing / n_total, 1),
      variance    = ifelse(is.na(val_var), NA, round(val_var, 4)),
      has_variance = has_var,
      is_usable   = is_usable,
      fail_reason = if (is_usable) NA
                    else if (n_valid == 0)   "no observations"
                    else if (n_valid < 50)   paste0("too few observations (n=", n_valid, ")")
                    else if (!has_var)       "no variance (constant value)"
                    else                     "unknown"
    )
  })
  names(tp_diagnostics) <- tp_cols

  tp_valid <- sapply(tp_diagnostics, function(x) x$is_usable)
  skip_result$n_timepoints_valid   <- sum(tp_valid)
  skip_result$timepoint_diagnostics <- tp_diagnostics

  if (sum(tp_valid) < 3) {
    cat("  Timepoint diagnostics:\n")
    for (tp in names(tp_diagnostics)) {
      diag   <- tp_diagnostics[[tp]]
      status <- if (diag$is_usable) "OK" else paste0("INVALID: ", diag$fail_reason)
      cat("    ", tp, ": n=", diag$n_valid, "/", diag$n_total,
          " (", diag$pct_missing, "% missing), var=",
          ifelse(is.na(diag$variance), "NA", diag$variance),
          " [", status, "]\n", sep = "")
    }
    cat("  Skipping - only", sum(tp_valid), "valid timepoints (need 3+)\n")
    skip_result$skip_reason  <- paste0("only ", sum(tp_valid), " valid timepoints: ",
                                       paste(sapply(tp_diagnostics[!tp_valid], function(x) x$fail_reason),
                                             collapse = "; "))
    skip_result$empty_timepoints <- tp_cols[!tp_valid]
    return(skip_result)
  }

  if (sum(tp_valid) < length(tp_cols)) {
    empty_tps <- tp_cols[!tp_valid]
    skip_result$empty_timepoints <- empty_tps

    cat("  Removing invalid timepoints:\n")
    for (tp in empty_tps) {
      diag <- tp_diagnostics[[tp]]
      cat("    ", tp, ": ", diag$fail_reason,
          " (n=", diag$n_valid, ", var=",
          ifelse(is.na(diag$variance), "NA", diag$variance), ")\n", sep = "")
    }

    tp_cols    <- tp_cols[tp_valid]
    n_tp       <- length(tp_cols)
    tp_nums    <- as.numeric(gsub(paste0("^", pheno_escaped, "_t"), "", tp_cols))
    time_codes <- tp_nums
    cat("  Using timepoints:", paste(tp_cols, collapse = ", "), "\n")
    cat("  Time codes:", paste(time_codes, collapse = ", "), "\n")
  } else {
    time_codes <- 0:(n_tp - 1)
  }

  skip_result$n_participants <- nrow(analysis_dat)
  cat("  Sample size (EUR unrelated):", nrow(analysis_dat), "\n")

  if (nrow(analysis_dat) < 100) {
    cat("  Skipping - sample size too small (<100)\n")
    skip_result$skip_reason <- paste0("sample size too small (n=", nrow(analysis_dat), ", need 100+)")
    return(skip_result)
  }

  time_codes_quad <- time_codes^2

  # Build model syntax
  int_loadings   <- paste(paste0("1*",              tp_cols), collapse = " + ")
  slope_loadings <- paste(paste0(time_codes,     "*", tp_cols), collapse = " + ")
  quad_loadings  <- paste(paste0(time_codes_quad, "*", tp_cols), collapse = " + ")

  # --- Model 1: No-growth ---
  no_growth_syntax <- paste0('
    intercept =~ ', int_loadings, '
    intercept ~ age_baseline
  ')

  fit_ng <- tryCatch({
    growth(no_growth_syntax, data = analysis_dat, missing = "fiml")
  }, error = function(e) {
    cat("  No-growth model failed:", e$message, "\n")
    return(NULL)
  })

  # --- Model 2: Linear growth ---
  linear_syntax <- paste0('
    intercept =~ ', int_loadings, '
    slope     =~ ', slope_loadings, '
    intercept ~~ slope
    intercept ~ age_baseline
    slope     ~ age_baseline
  ')

  fit_lin <- tryCatch({
    growth(linear_syntax, data = analysis_dat, missing = "fiml")
  }, error = function(e) {
    cat("  Linear growth model failed:", e$message, "\n")
    return(NULL)
  })

  # --- Model 3: Quadratic growth ---
  fit_quad <- NULL
  if (test_quadratic && n_tp >= 4) {
    quad_syntax <- paste0('
      intercept =~ ', int_loadings, '
      slope     =~ ', slope_loadings, '
      quad      =~ ', quad_loadings, '
      intercept ~~ slope
      intercept ~~ quad
      slope     ~~ quad
      intercept ~ age_baseline
      slope     ~ age_baseline
      quad      ~ age_baseline
    ')

    fit_quad <- tryCatch({
      growth(quad_syntax, data = analysis_dat, missing = "fiml")
    }, error = function(e) {
      cat("  Quadratic model failed:", e$message, "\n")
      return(NULL)
    })
  }

  # --- Extract comprehensive info from all models (Table 1 style) ---
  model_info_ng   <- extract_model_info(fit_ng,   "no_growth")
  model_info_lin  <- extract_model_info(fit_lin,  "linear")
  model_info_quad <- extract_model_info(fit_quad, "quadratic")

  # --- Assemble results list ---
  results <- list(
    phenotype            = pheno_name,
    status               = "fitted",
    skip_reason          = NULL,
    n_timepoints_initial = skip_result$n_timepoints_initial,
    n_timepoints_valid   = n_tp,
    n_participants       = nrow(analysis_dat),
    empty_timepoints     = skip_result$empty_timepoints,
    selected_model       = NULL,
    factor_scores        = NULL,
    fit_indices          = NULL,
    heywood_info         = NULL,
    convergence          = NULL,
    selection_reason     = NULL,
    fit_ng               = fit_ng,
    fit_lin              = fit_lin,
    fit_quad             = fit_quad,
    model_info_ng        = model_info_ng,
    model_info_lin       = model_info_lin,
    model_info_quad      = model_info_quad
  )

  # --- Convergence and Heywood checks ---
  conv_ng   <- check_converged(fit_ng)
  conv_lin  <- check_converged(fit_lin)
  conv_quad <- check_converged(fit_quad)

  fit_ng_measures   <- if (conv_ng)   tryCatch(fitMeasures(fit_ng,   c("cfi", "tli", "rmsea", "srmr", "aic", "bic")), error = function(e) NULL) else NULL
  fit_lin_measures  <- if (conv_lin)  tryCatch(fitMeasures(fit_lin,  c("cfi", "tli", "rmsea", "srmr", "aic", "bic")), error = function(e) NULL) else NULL
  fit_quad_measures <- if (conv_quad) tryCatch(fitMeasures(fit_quad, c("cfi", "tli", "rmsea", "srmr", "aic", "bic")), error = function(e) NULL) else NULL

  heywood_ng   <- if (conv_ng)   check_heywood(fit_ng)   else list(has_heywood = NA)
  heywood_lin  <- if (conv_lin)  check_heywood(fit_lin)  else list(has_heywood = NA)
  heywood_quad <- if (conv_quad) check_heywood(fit_quad) else list(has_heywood = NA)

  # Print fit indices
  fmt_fit <- function(m, label, heywood) {
    if (!is.null(m)) {
      hw <- if (has_heywood_safe(heywood)) " [HEYWOOD]" else ""
      cat("  ", label, ": CFI=", round(m["cfi"], 3),
          ", TLI=", round(m["tli"], 3),
          ", RMSEA=", round(m["rmsea"], 3), hw, "\n", sep = "")
    } else {
      cat("  ", label, ": DID NOT CONVERGE\n", sep = "")
    }
  }
  fmt_fit(fit_ng_measures,   "No-growth", heywood_ng)
  fmt_fit(fit_lin_measures,  "Linear   ", heywood_lin)
  if (!is.null(fit_quad)) fmt_fit(fit_quad_measures, "Quadratic", heywood_quad)

  results$heywood_info    <- list(no_growth = heywood_ng, linear = heywood_lin, quadratic = heywood_quad)
  results$convergence     <- list(no_growth = conv_ng,    linear = conv_lin,    quadratic = conv_quad)
  results$all_fit_indices <- list(no_growth = fit_ng_measures, linear = fit_lin_measures, quadratic = fit_quad_measures)

  results$lrt_results <- list(
    ng_vs_lin_chisq   = NA, ng_vs_lin_df   = NA, ng_vs_lin_pval   = NA,
    lin_vs_quad_chisq = NA, lin_vs_quad_df = NA, lin_vs_quad_pval = NA
  )

  # --- Model selection via nested LRT ---
  selected_fit        <- NULL
  selected_model_name <- NULL
  selected_measures   <- NULL
  selection_reason    <- NULL

  if (conv_ng && conv_lin) {
    comp_ng_lin <- tryCatch(anova(fit_ng, fit_lin), error = function(e) NULL)

    if (!is.null(comp_ng_lin)) {
      results$lrt_results$ng_vs_lin_chisq <- comp_ng_lin[["Chisq diff"]][2]
      results$lrt_results$ng_vs_lin_df    <- comp_ng_lin[["Df diff"]][2]
      results$lrt_results$ng_vs_lin_pval  <- comp_ng_lin[["Pr(>Chisq)"]][2]
      p_lin <- comp_ng_lin[["Pr(>Chisq)"]][2]

      if (!is.na(p_lin) && p_lin < 0.05) {
        if (has_heywood_safe(heywood_lin)) {
          selected_fit        <- fit_ng
          selected_model_name <- "no_growth"
          selected_measures   <- fit_ng_measures
          selection_reason    <- "linear has Heywood cases"
          cat("  Selected: NO-GROWTH (linear has Heywood)\n")
        } else {
          selected_fit        <- fit_lin
          selected_model_name <- "linear"
          selected_measures   <- fit_lin_measures

          if (!is.null(fit_quad) && conv_quad) {
            if (has_heywood_safe(heywood_quad)) {
              selection_reason <- "quadratic has Heywood cases"
              cat("  Selected: LINEAR (quadratic has Heywood)\n")
            } else {
              comp_lin_quad <- tryCatch(anova(fit_lin, fit_quad), error = function(e) NULL)

              if (!is.null(comp_lin_quad)) {
                results$lrt_results$lin_vs_quad_chisq <- comp_lin_quad[["Chisq diff"]][2]
                results$lrt_results$lin_vs_quad_df    <- comp_lin_quad[["Df diff"]][2]
                results$lrt_results$lin_vs_quad_pval  <- comp_lin_quad[["Pr(>Chisq)"]][2]
                p_quad <- comp_lin_quad[["Pr(>Chisq)"]][2]

                if (!is.na(p_quad) && p_quad < 0.05) {
                  selected_fit        <- fit_quad
                  selected_model_name <- "quadratic"
                  selected_measures   <- fit_quad_measures
                  selection_reason    <- "LRT significant, no Heywood"
                  cat("  Selected: QUADRATIC (p =", round(p_quad, 4), ")\n")
                } else {
                  selection_reason <- "quadratic LRT not significant"
                  cat("  Selected: LINEAR (quadratic p =", ifelse(is.na(p_quad), "NA", round(p_quad, 4)), ")\n")
                }
              } else {
                selection_reason <- "quadratic comparison failed"
                cat("  Selected: LINEAR (quadratic comparison failed)\n")
              }
            }
          } else if (!is.null(fit_quad) && !conv_quad) {
            selection_reason <- "quadratic did not converge"
            cat("  Selected: LINEAR (quadratic did not converge)\n")
          } else {
            selection_reason <- "no quadratic model"
            cat("  Selected: LINEAR (p =", round(p_lin, 4), ")\n")
          }
        }
      } else {
        selected_fit        <- fit_ng
        selected_model_name <- "no_growth"
        selected_measures   <- fit_ng_measures
        selection_reason    <- "linear LRT not significant"
        cat("  Selected: NO-GROWTH (p =", ifelse(is.na(p_lin), "NA", round(p_lin, 4)), ")\n")
      }
    } else {
      if (!has_heywood_safe(heywood_lin)) {
        selected_fit        <- fit_lin
        selected_model_name <- "linear"
        selected_measures   <- fit_lin_measures
        selection_reason    <- "LRT failed, using linear"
        cat("  Selected: LINEAR (LRT comparison failed)\n")
      } else {
        selected_fit        <- fit_ng
        selected_model_name <- "no_growth"
        selected_measures   <- fit_ng_measures
        selection_reason    <- "LRT failed, linear has Heywood"
        cat("  Selected: NO-GROWTH (LRT failed, linear has Heywood)\n")
      }
    }
  } else if (conv_ng && !conv_lin) {
    selected_fit        <- fit_ng
    selected_model_name <- "no_growth"
    selected_measures   <- fit_ng_measures
    selection_reason    <- "only no-growth converged"
    cat("  Selected: NO-GROWTH (only model that converged)\n")
  } else if (!conv_ng && conv_lin) {
    if (!has_heywood_safe(heywood_lin)) {
      selected_fit        <- fit_lin
      selected_model_name <- "linear"
      selected_measures   <- fit_lin_measures
      selection_reason    <- "only linear converged"
      cat("  Selected: LINEAR (no-growth did not converge)\n")
    } else {
      cat("  SKIPPING: no-growth failed, linear has Heywood\n")
      results$status      <- "failed"
      results$skip_reason <- "no-growth failed, linear has Heywood"
    }
  } else {
    cat("  SKIPPING: no models converged\n")
    results$status      <- "failed"
    results$skip_reason <- "no models converged"
  }

  results$selected_model   <- selected_model_name
  results$fit_indices      <- selected_measures
  results$selection_reason <- selection_reason

  if (!is.null(selected_model_name)) results$status <- "success"

  # --- Extract factor scores ---
  if (!is.null(selected_fit)) {
    fs <- tryCatch(lavPredict(selected_fit), error = function(e) NULL)

    if (!is.null(fs)) {
      used_cases <- lavInspect(selected_fit, "case.idx")

      score_df <- data.frame(
        src_subject_id = analysis_dat$src_subject_id[used_cases],
        intercept      = fs[, "intercept"]
      )

      if (selected_model_name %in% c("linear", "quadratic") && "slope" %in% colnames(fs)) {
        score_df$slope <- fs[, "slope"]
      }
      if (selected_model_name == "quadratic" && "quad" %in% colnames(fs)) {
        score_df$quad <- fs[, "quad"]
      }

      results$factor_scores <- score_df
      cat("  Factor scores extracted for", nrow(score_df), "EUR participants\n")
    }
  }

  return(results)
}

# ============================================================================
# PART 8: RUN LGC MODELS FOR ALL PHENOTYPES
# ============================================================================

# Uncomment to test on a subset first:
# phenos_for_lgc <- phenos_for_lgc[1:10]

all_lgc_results   <- list()
all_factor_scores <- list()
all_phenotype_log <- list()

for (pheno in phenos_for_lgc) {
  result <- fit_lgc_model(pheno, pheno_list, covariates, test_quadratic = TRUE)

  all_phenotype_log[[pheno]] <- result

  if (!is.null(result) && result$status == "success" && !is.null(result$factor_scores)) {
    all_lgc_results[[pheno]] <- result

    scores <- result$factor_scores
    names(scores)[names(scores) == "intercept"] <- paste0(pheno, "_intercept")
    if ("slope" %in% names(scores)) names(scores)[names(scores) == "slope"] <- paste0(pheno, "_slope")
    if ("quad"  %in% names(scores)) names(scores)[names(scores) == "quad"]  <- paste0(pheno, "_quad")

    all_factor_scores[[pheno]] <- scores
  }
}

# ============================================================================
# PART 9: COMPREHENSIVE SUMMARY
# ============================================================================

cat("\n")
cat("################################################################################\n")
cat("#                         COMPREHENSIVE LGC SUMMARY                            #\n")
cat("################################################################################\n")

status_counts <- table(sapply(all_phenotype_log, function(x) x$status))
cat("\n=== OVERALL STATUS ===\n")
cat("Total phenotypes attempted:", length(phenos_for_lgc), "\n")
print(status_counts)

skipped <- all_phenotype_log[sapply(all_phenotype_log, function(x) x$status == "skipped")]
if (length(skipped) > 0) {
  cat("\n=== SKIPPED PHENOTYPES (", length(skipped), ") ===\n")
  skipped_df <- do.call(rbind, lapply(names(skipped), function(pheno) {
    res <- skipped[[pheno]]
    data.frame(
      phenotype            = pheno,
      reason               = ifelse(is.null(res$skip_reason), "unknown", res$skip_reason),
      n_timepoints_initial = ifelse(is.na(res$n_timepoints_initial), NA, res$n_timepoints_initial),
      n_timepoints_valid   = ifelse(is.na(res$n_timepoints_valid),   NA, res$n_timepoints_valid),
      n_participants       = ifelse(is.na(res$n_participants),        NA, res$n_participants)
    )
  }))
  cat("\n  Skip reasons summary:\n")
  print(table(skipped_df$reason))
}

failed <- all_phenotype_log[sapply(all_phenotype_log, function(x) x$status == "failed")]
if (length(failed) > 0) {
  cat("\n=== FAILED PHENOTYPES (", length(failed), ") ===\n")
  for (pheno in names(failed)) {
    cat("  ", pheno, ": ", failed[[pheno]]$skip_reason, "\n")
  }
}

successful <- all_phenotype_log[sapply(all_phenotype_log, function(x) x$status == "success")]
cat("\n=== SUCCESSFUL PHENOTYPES (", length(successful), ") ===\n")

if (length(successful) > 0) {
  model_types  <- sapply(successful, function(x) x$selected_model)
  cat("\nModel selection breakdown:\n")
  print(table(model_types))

  sample_sizes <- sapply(successful, function(x) x$n_participants)
  cat("\nSample size range:", min(sample_sizes), "-", max(sample_sizes),
      "(median:", median(sample_sizes), ")\n")
}

cat("\n################################################################################\n")

# ============================================================================
# PART 10: MERGE AND SAVE FACTOR SCORES
# ============================================================================

if (length(all_factor_scores) > 0) {
  merged_scores <- all_factor_scores[[1]]
  if (length(all_factor_scores) > 1) {
    for (i in 2:length(all_factor_scores)) {
      merged_scores <- merge(merged_scores, all_factor_scores[[i]],
                             by = "src_subject_id", all = TRUE)
    }
  }

  cat("\nFinal merged dataset (EUR only):\n")
  cat("  Participants:", nrow(merged_scores), "\n")
  cat("  Variables:",    ncol(merged_scores), "\n")

  fwrite(merged_scores,
         file.path(OUTPUT_DIR, "LGC_factor_scores_mental_health_EUR.csv"),
         row.names = FALSE)
  cat("Saved:", file.path(OUTPUT_DIR, "LGC_factor_scores_mental_health_EUR.csv"), "\n")
}

# ============================================================================
# PART 11: PREPARE GWAS-READY FILES
# ============================================================================

cat("\n========================================\n")
cat("PREPARING GWAS-READY FILE FOR EUR\n")
cat("========================================\n")

# Load ancestry PCs; keep only the first 10 for GWAS covariates
ABCD_PCs <- fread(PCS_FILE)
PCs_EUR  <- ABCD_PCs[ABCD_PCs$src_subject_id %in% EUR_ids, ]

gwas_dat <- merge(merged_scores, covariates, by = "src_subject_id", all.x = TRUE)
gwas_dat <- merge(gwas_dat, PCs_EUR[, c("src_subject_id", paste0("V", 1:10))],
                  by = "src_subject_id")

# Keep only unrelated EUR participants with valid sex coding
gwas_dat <- gwas_dat[gwas_dat$src_subject_id %in% EUR_unrelated_ids, ]
gwas_dat <- gwas_dat[gwas_dat$demo_sex_v2 %in% c(1, 2), ]

cat("Final EUR GWAS sample size:", nrow(gwas_dat), "\n")

fwrite(gwas_dat,
       file.path(OUTPUT_DIR, "LGC_phenotypes_GWAS_ready_EUR.csv"),
       row.names = FALSE)
cat("Saved:", file.path(OUTPUT_DIR, "LGC_phenotypes_GWAS_ready_EUR.csv"), "\n")

# Intercept-only file (primary GWAS input)
all_cols       <- names(gwas_dat)
intercept_cols <- grep("_intercept$", all_cols, value = TRUE)
covariate_cols <- c("src_subject_id", "demo_sex_v2", "age_baseline", paste0("V", 1:10))

gwas_intercepts_only <- gwas_dat[, c(covariate_cols, intercept_cols)]

cat("Intercept phenotypes for GWAS:", length(intercept_cols), "\n")

fwrite(gwas_intercepts_only,
       file.path(OUTPUT_DIR, "LGC_intercepts_GWAS_ready_EUR.csv"),
       row.names = FALSE)
cat("Saved:", file.path(OUTPUT_DIR, "LGC_intercepts_GWAS_ready_EUR.csv"), "\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("COMPLETE!\n")
cat("========================================\n")

cat("\n=== OUTPUT FILES ===\n")
cat("1. LGC_factor_scores_mental_health_EUR.csv\n")
cat("   - All factor scores (intercept / slope / quadratic)\n\n")
cat("2. LGC_phenotypes_GWAS_ready_EUR.csv\n")
cat("   - All factor scores + covariates + 10 PCs\n\n")
cat("3. LGC_intercepts_GWAS_ready_EUR.csv\n")
cat("   - PRIMARY GWAS FILE (intercepts only + covariates + 10 PCs)\n\n")
cat("4. ABCD_mental_health_phenotype_availability_EUR.csv\n")
cat("   - Phenotype availability across timepoints\n\n")

cat("=== MODEL SELECTION SUMMARY ===\n")
if (length(successful) > 0) {
  print(table(sapply(successful, function(x) x$selected_model)))
}

cat("\n=== GWAS RECOMMENDATIONS ===\n")
cat("1. PRIMARY: Use intercepts for GWAS (", length(intercept_cols), "phenotypes)\n")
cat("2. Covariates: genetic sex, age_baseline, PC1-PC10\n")
cat("3. Tool: PLINK2 --glm for autosomes\n")
cat("4. For heritability: GCTA-GREML on intercepts\n")
