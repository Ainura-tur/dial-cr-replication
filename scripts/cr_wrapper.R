#================================================================
# DIAL_NEURIPS_CR_wrapper.R  —  v3.1  (Phase 1 corrections applied)
#
# CHANGES IN THIS VERSION (relative to v3)
# ─────────────────────────────────────────────────────────────────
# [P1-WRAPPER-1]  .empirical_studies() Banerjee canonical updated
#                 from c(-0.8, -0.6) to c(-0.6, -0.4) to match the
#                 corrected run_banerjee_analysis_fixed() in
#                 Fin_Empirical5_clean.R.  This ensures scenario4
#                 (rxu_range sensitivity) uses the same canonical
#                 domain as the Table 5/6 main results.
#                 OPEN: PDF page-15 comment says "-0.5 upper bound".
#                 Run DIAL_Phase1_runner.R to compare and choose.
#
# [P1-WRAPPER-2]  .empirical_studies() Burde canonical KEPT as
#                 c(0.0, 0.4).  The Phase 1 sensitivity check
#                 c(-0.4, 0) is available via run_burde_sensitivity_fixed()
#                 and is run by DIAL_Phase1_runner.R, but it is NOT the
#                 canonical domain for scenario4 until the structural
#                 justification for the negative sign is documented.
#
# [P1-WRAPPER-3]  .empirical_studies() Di Tella canonical updated from
#                 c(0, 0.25) to c(0.15, 0.25) to match Table 4 of the
#                 NeurIPS submission. This is the within-sign canonical
#                 elicited per Section 5.4.
#
# SOURCE ORDER IS MANDATORY — do not change:
#   source("Fin_sim3_clean.R")
#   source("Fin_Empirical5_clean.R")
#   source("DIAL_NEURIPS_CR_wrapper.R")
#
# PRIOR BUG HISTORY (unchanged from v3):
#  [B1-orig] check_compatibility() → check_compatibility_simple()
#  [B2-orig] r_grid ignored rxu_range in analyze_single_sample()
#  [B3-orig] standardize_result() never called
#  [NEW-B1]  seed= removed from CIhybrid() call in wrapper
#  [NEW-B2]  standardize_cr_result() rewritten (row-safe)
#  [NEW-B3]  Default rxu_range = c(-0.4, 0.4) for simulations
#  [BUG-A]   Zero_in_CI_MCUB rename key corrected
#  [BUG-B]   enrich_s6() added to fill r_xz, r_zu, beta_IV
#================================================================

suppressPackageStartupMessages({
  library(MASS); library(ivreg); library(lmtest); library(sandwich)
  library(haven); library(dplyr); library(tidyr); library(Matrix)
  library(fGarch)
})

#================================================================
# SECTION 0 — Source guard and session diagnostics
#================================================================

check_sources_loaded <- function() {
  sim_fns <- c("generate_data", "monte_carlo_simulation",
               "analyze_single_sample", "summarize_sim_results", "fs")
  emp_fns <- c("ci_simple_union", "CIhybrid", "CIproj_p", "CIcon",
               "CIcon_TNbounds", "g_xu_safe",
               "local_compute_gradient_safe",
               "estimate_cov_corr_boot",
               "check_compatibility_simple",
               "pvalue_mcub_zero_fast",
               "run_ditella_analysis_fixed",
               "run_burde_analysis_fixed",
               "run_burde_sensitivity_fixed",    # [P1-WRAPPER-2] new function
               "run_galiani_analysis_fixed",
               "run_banerjee_analysis_fixed")
  
  miss_s <- sim_fns[!vapply(sim_fns, exists, logical(1))]
  miss_e <- emp_fns[!vapply(emp_fns, exists, logical(1))]
  
  if (length(miss_s))
    warning("Fin_sim3_clean.R not fully loaded. Missing: ",
            paste(miss_s, collapse = ", "))
  if (length(miss_e))
    warning("Fin_Empirical5_clean.R not fully loaded. Missing: ",
            paste(miss_e, collapse = ", "))
  
  ci_args <- tryCatch(names(formals(CIhybrid)), error = function(e) character(0))
  if ("seed" %in% ci_args) {
    cat("  [OK] CIhybrid() has seed parameter (Fin_Empirical5_clean.R version active)\n")
  } else {
    warning("[NEW-B1] Active CIhybrid() lacks 'seed' parameter. ",
            "Source Fin_Empirical5_clean.R AFTER Fin_sim3_clean.R.")
  }
  
  # [P1-WRAPPER-1] confirm Banerjee canonical is corrected
  if (exists("run_banerjee_analysis_fixed")) {
    fn_body <- deparse(body(run_banerjee_analysis_fixed))
    if (any(grepl("-0.8, -0.6", fn_body, fixed = TRUE))) {
      warning("[P1-WRAPPER-1] run_banerjee_analysis_fixed() still uses old ",
              "rxu_range = c(-0.8, -0.6). Re-source Fin_Empirical5_clean.R.")
    } else if (any(grepl("-0.6, -0.4", fn_body, fixed = TRUE))) {
      cat("  [OK] run_banerjee_analysis_fixed() uses corrected c(-0.6, -0.4)\n")
    }
  }
  
  invisible(length(miss_s) + length(miss_e) == 0)
}

#================================================================
# SECTION 1 — Unified CR test entry point
#================================================================

cr_test_unified <- function(
    df          = NULL,
    data_fn     = NULL,
    N           = 800L,
    seed        = 123L,
    rxu_range   = c(-0.4, 0.4),   # [NEW-B3] symmetric default for simulations
    alpha       = 0.05,
    method      = c("both", "mcub", "simple"),
    label       = "unnamed",
    i           = 1L,
    B_mcub      = 500L,
    Blarge_mcub = 5000L,
    verbose     = TRUE
) {
  method <- match.arg(method)
  
  if (is.null(df) && is.null(data_fn))
    stop("cr_test_unified: supply either df or data_fn.")
  
  if (is.null(df)) {
    set.seed(seed)
    raw <- data_fn(N = N, seed = seed)
    if (!all(c("x","y","z") %in% names(raw)))
      stop("data_fn must return a data.frame with columns x, y, z.")
    df <- raw[, c("x","y","z"), drop = FALSE]
    if ("u" %in% names(raw)) df$u <- raw$u
  }
  
  df <- df[complete.cases(df[, c("x","y","z")]), , drop = FALSE]
  n  <- nrow(df)
  if (n < 30) { warning("cr_test_unified: n = ", n, " < 30 — skipping."); return(NULL) }
  
  if (verbose)
    cat(sprintf("  [%s] n = %d | rxu in [%.2f, %.2f]\n",
                label, n, rxu_range[1], rxu_range[2]))
  
  rho_xy   <- cor(df$x, df$y, use = "complete.obs")
  rho_xz   <- cor(df$x, df$z, use = "complete.obs")
  rho_yz   <- cor(df$y, df$z, use = "complete.obs")
  deltahat <- c(rho_xy, rho_xz, rho_yz)
  
  r_zu    <- if ("u" %in% names(df))
    cor(df$z, df$u, use = "complete.obs") else NA_real_
  beta_IV <- cov(df$z, df$y, use = "complete.obs") /
    cov(df$z, df$x, use = "complete.obs")
  
  deltaSigma <- estimate_cov_corr_boot(df$x, df$y, df$z, B = 800L, seed = 1000L + i)
  
  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = 50L)
  g      <- function(delta) g_xu_safe(r_grid, delta)
  
  A  <- t(sapply(r_grid, function(r_xu)
    as.numeric(local_compute_gradient_safe(r_xu, deltahat[1], deltahat[2], deltahat[3]))))
  Al <- A; Au <- A
  
  eta <- 0.001; alphac <- 0.8 * alpha; tol <- 1e-3; tol_r <- 1e-3
  
  CI_s <- tryCatch({
    res <- ci_simple_union(df$x, df$y, df$z,
                           rxu_range = rxu_range, alpha = alpha,
                           cov_method = "bootstrap", seed = seed + i)
    list(CI = res$CI, plug_in = res$plug_in, ok = TRUE)
  }, error = function(e) {
    if (verbose) cat("    simple CI failed:", conditionMessage(e), "\n")
    list(CI = c(NA, NA), plug_in = c(NA, NA), ok = FALSE)
  })
  
  CI_b   <- CI_s; p_zero <- NA_real_
  
  if (method %in% c("both","mcub") && isTRUE(CI_s$ok)) {
    CI_b <- tryCatch({
      set.seed(2000L + i)
      res <- CIhybrid(deltahat, deltaSigma, Al, Au,
                      alpha = alpha, alphac = alphac, eta = eta,
                      B = B_mcub, Blarge = Blarge_mcub,
                      tol = tol, tol_r = tol_r, index = NULL, g = g)
      list(CI = res$CI_h, plug_in = res$CI_c, ok = TRUE)
    }, error = function(e) {
      if (verbose) cat("    MCUB failed:", conditionMessage(e), "\n")
      CI_s
    })
    
    if (isTRUE(CI_b$ok)) {
      p_zero <- tryCatch({
        set.seed(3000L + i)
        pvalue_mcub_zero_fast(deltahat, deltaSigma, Al, Au, g,
                              eta = eta, tol = tol, tol_r = tol_r,
                              B_fast = 300L, Blarge_fast = 3000L)
      }, error = function(e) NA_real_)
    }
  }
  
  data.frame(
    Study        = label,
    Instrument   = paste0("Z", i),
    n            = n,
    rxu_lo       = rxu_range[1],
    rxu_hi       = rxu_range[2],
    r_xz         = round(rho_xz,  3),
    r_zu         = round(r_zu,    3),
    beta_IV      = round(beta_IV, 3),
    plug_in      = sprintf("[%.3f, %.3f]", CI_b$plug_in[1], CI_b$plug_in[2]),
    CI_MCUB      = sprintf("[%.3f, %.3f]", CI_b$CI[1], CI_b$CI[2]),
    CI_simple    = sprintf("[%.3f, %.3f]", max(-1, CI_s$CI[1]), min(1, CI_s$CI[2])),
    width_MCUB   = round(diff(CI_b$CI), 3),
    width_simple = round(diff(pmax(-1, pmin(1, CI_s$CI))), 3),
    Zero_MCUB    = ifelse(!is.na(CI_b$CI[1]) & CI_b$CI[1] <= 0 & CI_b$CI[2] >= 0,
                          "\u2713", "\u00d7"),
    Zero_simple  = ifelse(!is.na(CI_s$CI[1]) & CI_s$CI[1] <= 0 & CI_s$CI[2] >= 0,
                          "\u2713", "\u00d7"),
    p_zero       = round(p_zero, 3),
    verdict      = ifelse(!is.na(CI_b$CI[1]) & CI_b$CI[1] <= 0 & CI_b$CI[2] >= 0,
                          "Valid", "Invalid"),
    stringsAsFactors = FALSE
  )
}

#================================================================
# SECTION 2 — Output standardisation
#================================================================

standardize_cr_result <- function(res, study_name) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  res$Study <- study_name
  
  # [BUG-A] rename covers both column name variants
  rename_map <- c(CI_Bei          = "CI_MCUB",
                  Zero_in_CI_MCUB = "Zero_MCUB",
                  Zero_in_CI      = "Zero_MCUB")
  for (old_nm in names(rename_map)) {
    new_nm <- rename_map[[old_nm]]
    if (old_nm %in% names(res) && !new_nm %in% names(res))
      names(res)[names(res) == old_nm] <- new_nm
  }
  
  if (!"Instrument" %in% names(res) && "Z" %in% names(res))
    res$Instrument <- res$Z
  
  unified_cols <- c("Study","Instrument","n","rxu_lo","rxu_hi",
                    "r_xz","r_zu","beta_IV","plug_in",
                    "CI_MCUB","CI_simple","width_MCUB","width_simple",
                    "Zero_MCUB","Zero_simple","p_zero","verdict")
  for (col in unified_cols) {
    if (!col %in% names(res))
      res[[col]] <- if (col %in% c("rxu_lo","rxu_hi","width_MCUB",
                                   "width_simple","r_xz","r_zu","beta_IV"))
        NA_real_ else NA_character_
  }
  
  needs_verdict <- is.na(res$verdict) | res$verdict %in% c("NA","")
  if (any(needs_verdict, na.rm = TRUE))
    res$verdict[needs_verdict] <- ifelse(
      !is.na(res$Zero_MCUB[needs_verdict]) &
        res$Zero_MCUB[needs_verdict] == "\u2713",
      "Valid", "Invalid")
  
  res[, unified_cols, drop = FALSE]
}

#================================================================
# SECTION 3 — Scenario 1: d_vec validity spectrum
#================================================================

scenario1_dvec_spectrum <- function(
    d_vec     = c(1.25, 1.10, 0.95, 0.83, 0.70, 0.62,
                  0.50, 0.45, 0.35, 0.20, 0.05, -0.05),
    N         = 800L,
    rxu_range = c(-0.4, 0.4),
    alpha     = 0.05,
    k         = 1,
    seed      = 9058L,
    verbose   = TRUE
) {
  cat("\n=== Scenario 1: d_vec validity spectrum ===\n")
  cat(sprintf("  N = %d | rxu_range = [%.2f, %.2f] | seed = %d\n",
              N, rxu_range[1], rxu_range[2], seed))
  
  sim_data <- generate_data(N = N, seed = seed)
  rows     <- vector("list", length(d_vec))
  
  for (j in seq_along(d_vec)) {
    d    <- d_vec[j]
    df   <- sim_data
    df$z <- df$x - k * d * df$R + df$A0
    
    rows[[j]] <- cr_test_unified(df = df, rxu_range = rxu_range, alpha = alpha,
                                 label = sprintf("sim_d=%.2f", d), i = j, verbose = verbose)
    if (!is.null(rows[[j]])) rows[[j]]$d <- d
  }
  
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  rownames(out) <- NULL
  
  if (!is.null(out) && nrow(out) > 1) {
    flips <- which(diff(as.integer(out$verdict == "Valid")) != 0)
    if (length(flips)) {
      cat(sprintf(
        "  Verdict flips between d = %.2f (%s, r_zu = %.3f) and d = %.2f (%s, r_zu = %.3f)\n",
        out$d[flips[1]], out$verdict[flips[1]], out$r_zu[flips[1]],
        out$d[flips[1]+1], out$verdict[flips[1]+1], out$r_zu[flips[1]+1]))
    } else {
      cat("  No verdict flip detected across d_vec.\n")
    }
  }
  out
}

#================================================================
# SECTION 4 — Scenario 2: Monte Carlo power analysis
#================================================================

scenario2_mc_power <- function(
    d_boundary = 0.62,
    N_vec      = c(500L, 800L, 2000L),
    n_sim      = 50L,
    rxu_range  = c(-0.4, 0.4),
    alpha      = 0.05,
    k          = 1,
    seed0      = 9058L,
    verbose    = TRUE
) {
  cat("\n=== Scenario 2: Monte Carlo power analysis ===\n")
  cat(sprintf("  d = %.2f | n_sim = %d | rxu_range = [%.2f, %.2f]\n",
              d_boundary, n_sim, rxu_range[1], rxu_range[2]))
  
  power_rows <- vector("list", length(N_vec))
  
  for (ni in seq_along(N_vec)) {
    N <- N_vec[ni]; cat(sprintf("  N = %d ...\n", N))
    sim_rows <- vector("list", n_sim)
    pb       <- txtProgressBar(min = 0, max = n_sim, style = 3)
    
    for (s in seq_len(n_sim)) {
      sim_data <- generate_data(N = N, seed = seed0 + s)
      df       <- sim_data; df$z <- df$x - k * d_boundary * df$R + df$A0
      sim_rows[[s]] <- tryCatch(
        cr_test_unified(df = df, rxu_range = rxu_range, alpha = alpha,
                        label = "power_sim", i = s, verbose = FALSE),
        error = function(e) NULL)
      setTxtProgressBar(pb, s)
    }
    close(pb)
    
    all_s <- do.call(rbind, sim_rows[!vapply(sim_rows, is.null, logical(1))])
    if (is.null(all_s) || nrow(all_s) == 0) {
      power_rows[[ni]] <- data.frame(N = N, n_sim = 0, d = d_boundary,
                                     reject_MCUB = NA, reject_simple = NA,
                                     mean_width_MCUB = NA, mean_width_simple = NA,
                                     mean_p_zero = NA, stringsAsFactors = FALSE)
      next
    }
    
    power_rows[[ni]] <- data.frame(
      N                 = N,
      n_sim             = nrow(all_s),
      d                 = d_boundary,
      reject_MCUB       = round(mean(all_s$Zero_MCUB   == "\u00d7", na.rm = TRUE), 3),
      reject_simple     = round(mean(all_s$Zero_simple  == "\u00d7", na.rm = TRUE), 3),
      mean_width_MCUB   = round(mean(all_s$width_MCUB,   na.rm = TRUE), 3),
      mean_width_simple = round(mean(all_s$width_simple, na.rm = TRUE), 3),
      mean_p_zero       = round(mean(all_s$p_zero,       na.rm = TRUE), 4),
      stringsAsFactors  = FALSE
    )
  }
  
  out <- do.call(rbind, power_rows)
  cat("\n  Power table (rejection rate at alpha =", alpha, "):\n")
  print(out, digits = 3)
  out
}

#================================================================
# SECTION 5 — Scenario 3: Non-normality stress test
#================================================================

scenario3_nonnormal_stress <- function(
    d_vec     = c(0.62, 0.45, 0.20),
    N         = 800L,
    rxu_range = c(-0.4, 0.4),
    alpha     = 0.05,
    k         = 1,
    seed      = 9058L,
    verbose   = TRUE
) {
  cat("\n=== Scenario 3: Non-normality stress test ===\n")
  sim_data <- generate_data(N = N, seed = seed)
  rows <- vector("list", length(d_vec) * 2); idx <- 1L
  
  for (d in d_vec) {
    df <- sim_data; df$z <- df$x - k * d * df$R + df$A0
    
    for (cov_m in c("bootstrap","diag")) {
      res_s <- tryCatch(
        ci_simple_union(df$x, df$y, df$z, rxu_range = rxu_range,
                        alpha = alpha, cov_method = cov_m, seed = seed),
        error = function(e) NULL)
      
      width_s <- if (!is.null(res_s)) round(diff(pmax(-1, pmin(1, res_s$CI))), 4) else NA_real_
      zero_s  <- if (!is.null(res_s)) res_s$CI[1] <= 0 & res_s$CI[2] >= 0 else NA
      
      rows[[idx]] <- data.frame(d = d, cov_method = cov_m,
                                CI_lo = if (!is.null(res_s)) round(res_s$CI[1], 3) else NA_real_,
                                CI_hi = if (!is.null(res_s)) round(res_s$CI[2], 3) else NA_real_,
                                width_simple = width_s,
                                zero_in_CI = ifelse(isTRUE(zero_s), "\u2713", "\u00d7"),
                                stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
  
  out  <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  wide <- out %>%
    tidyr::pivot_wider(id_cols = "d", names_from = "cov_method",
                       values_from = c("width_simple","zero_in_CI")) %>%
    dplyr::mutate(width_inflation = round(width_simple_diag / width_simple_bootstrap, 3))
  
  cat("  Width inflation factor (diag / bootstrap):\n")
  print(wide, digits = 3)
  list(full = out, summary = wide)
}

#================================================================
# SECTION 6 — Scenario 4: rxu_range sensitivity (empirical)
#================================================================

.residualise_df <- function(data, Y, X, Z, controls) {
  vars <- c(Y, X, Z, controls)
  d    <- data[complete.cases(data[, vars]), vars, drop = FALSE]
  fml  <- function(lhs) as.formula(paste(lhs, "~", paste(controls, collapse = "+")))
  y    <- resid(lm(fml(Y), data = d))
  x    <- resid(lm(fml(X), data = d))
  z0   <- resid(lm(fml(Z), data = d))
  z    <- predict(lm(z0 ~ x + y))
  data.frame(x = x, y = y, z = z)
}

.empirical_studies <- function() {
  base <- paste0("https://github.com/ratbekd/Orientation_paper/",
                 "raw/refs/heads/main/")
  list(
    list(
      name     = "Di Tella & Schargrodsky",
      loader   = function() {
        d <- haven::read_dta(paste0(base,
                                    "JPE%20-%20Di%20Tella%20and%20Schargrodsky%20-%20",
                                    "CriminalRecidivismAfterPrisonAndElectronicMonitoring%20.dta"))
        d[d$offendersPerJudge > 9, ]
      },
      Y = "recidivism", X = "electronicMonitoring",
      Z = "judgeAlreadyUsedEM",
      Zs = c("judgeAlreadyUsedEM","percJudgeSentToEM"),
      controls = c("mostSeriousCrime","age","ageSquared","argentine",
                   "numberPreviousImprisonments","judicialDistrict",
                   "yearOfImprisonment"),
      # [P1-WRAPPER-3] Di Tella canonical: (0.15, 0.25) for the primary
      # instrument (judgeAlreadyUsedEM). The companion percJudgeSentToEM
      # uses (0, 0.25) at run time, set per-instrument by the Phase 1
      # runner Dial_phase1_runner.R; the value below applies only to the
      # scenario4 sensitivity sweep, which is anchored on the primary
      # instrument's canonical.
      canonical = c(0.15, 0.25)
    ),
    list(
      name   = "Burde & Linden",
      loader = function() haven::read_dta(paste0(base, "afgan.dta")),
      Y = "testscore", X = "enrolled", Z = "buildschool",
      Zs = "buildschool",
      controls = c("headchild","age","yrsvill","farsi","tajik","farmers",
                   "agehead","educhead","nhh","land","sheep",
                   "distschool","chagcharan"),
      canonical = c(0.0, 0.4)    # [P1-WRAPPER-2] unchanged pending structural justification
    ),
    list(
      name   = "Galiani et al.",
      loader = function() {
        d <- haven::read_dta(paste0(base, "Crime.dta"))
        d[d$cohort > 1957 & d$cohort < 1963, ]
      },
      Y = "crimerate", X = "sm", Z = "highnumber",
      Zs = "highnumber",
      controls = c("cohort","draftnumber","navy"),
      canonical = c(-0.2, 0.0)
    ),
    list(
      name   = "Banerjee et al.",
      loader = function() {
        d <- haven::read_dta(paste0(base, "yld_sett_aug03.dta"))
        d[d$year >= 1965 & d$phwht <= 1, ]
      },
      Y = "phwht", X = "p_nland", Z = "instru",
      Zs = "instru",
      controls = c("alt","totrain","so_black","so_red","so_all",
                   "lat","coastal","brule1","year"),
      # [P1-WRAPPER-1] CORRECTED from c(-0.8, -0.6) to c(-0.6, -0.4)
      # to match run_banerjee_analysis_fixed() and the sign-fix hypothesis.
      # NOTE: PDF page-15 says "-0.5 upper bound". Run DIAL_Phase1_runner.R
      # to compare and finalise before updating Table 6.
      canonical = c(-0.6, -0.4)
    )
  )
}

scenario4_rxu_sensitivity <- function(alpha = 0.05, verbose = TRUE) {
  cat("\n=== Scenario 4: rxu_range sensitivity ===\n")
  studies   <- .empirical_studies()
  wide_canonical <- list(
    "Di Tella & Schargrodsky" = c(0.00,  0.80),
    "Burde & Linden"          = c(0.00,  0.80),
    "Galiani et al."          = c(-0.80, 0.00),
    "Banerjee et al."         = c(-0.80, 0.00)
  )
  out_rows  <- list()
  
  for (st in studies) {
    if (verbose) cat(sprintf("  Loading %s ...\n", st$name))
    raw <- tryCatch(st$loader(), error = function(e) {
      cat("    Data load failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(raw)) next
    
    df <- tryCatch(
      .residualise_df(raw, st$Y, st$X, st$Z, st$controls),
      error = function(e) { cat("    Residualise failed.\n"); NULL })
    if (is.null(df)) next
    
    for (rng_name in c("canonical","wide_canonical")) {
      rng <- if (rng_name == "canonical") st$canonical else wide_canonical[[st$name]]
      res <- tryCatch(
        cr_test_unified(df = df, rxu_range = rng, alpha = alpha,
                        label = st$name, i = 1L, verbose = FALSE),
        error = function(e) NULL)
      if (!is.null(res)) {
        res$range_type <- rng_name
        out_rows[[length(out_rows) + 1]] <- res
      }
    }
  }
  
  out <- do.call(rbind, out_rows[!vapply(out_rows, is.null, logical(1))])
  rownames(out) <- NULL
  
  if (!is.null(out) && nrow(out) > 0) {
    sensitivity <- out %>%
      dplyr::select(Study, range_type, rxu_lo, rxu_hi, CI_MCUB, verdict) %>%
      tidyr::pivot_wider(id_cols = "Study", names_from = "range_type",
                         values_from = c("CI_MCUB","verdict")) %>%
      dplyr::mutate(robust = verdict_canonical == verdict_symmetric)
    cat("\n  Sensitivity summary (robust = verdict unchanged):\n")
    print(sensitivity)
  }
  out
}

#================================================================
# SECTION 7 — Scenario 5: MCUB vs. simple union comparator
#================================================================

scenario5_ci_comparator <- function(
    df_list   = NULL,
    N         = 800L,
    d_vec     = c(0.62, 0.45, 0.20),
    rxu_range = c(-0.4, 0.4),
    alpha     = 0.05,
    k         = 1,
    seed      = 9058L,
    verbose   = TRUE
) {
  cat("\n=== Scenario 5: MCUB vs. simple union comparator ===\n")
  
  if (is.null(df_list)) {
    sim_data <- generate_data(N = N, seed = seed)
    df_list  <- stats::setNames(
      lapply(d_vec, function(d) { df <- sim_data; df$z <- df$x - k*d*df$R + df$A0; df }),
      sprintf("d=%.2f", d_vec))
  }
  
  rows <- vector("list", length(df_list))
  for (j in seq_along(df_list)) {
    res <- tryCatch(
      cr_test_unified(df = df_list[[j]], rxu_range = rxu_range, alpha = alpha,
                      label = names(df_list)[j], i = j, verbose = verbose),
      error = function(e) NULL)
    if (!is.null(res)) rows[[j]] <- res
  }
  
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  rownames(out) <- NULL
  
  if (!is.null(out) && nrow(out) > 0) {
    out <- out %>%
      dplyr::mutate(width_ratio = round(width_MCUB / width_simple, 3),
                    agree       = Zero_MCUB == Zero_simple)
    cat("\n  CI comparator summary:\n")
    print(out[, c("Study","r_xz","r_zu","CI_MCUB","CI_simple",
                  "width_MCUB","width_simple","width_ratio",
                  "Zero_MCUB","Zero_simple","agree","p_zero")])
  }
  out
}

#================================================================
# SECTION 8 — Scenario 6: Multi-study unified empirical table
#================================================================

scenario6_multistudy_table <- function(alpha = 0.05, enrich = TRUE, verbose = TRUE) {
  cat("\n=== Scenario 6: Multi-study unified table ===\n")
  
  study_runners <- list(
    "Di Tella & Schargrodsky" = function()
      tryCatch(run_ditella_analysis_fixed(alpha = alpha),
               error = function(e) { message(conditionMessage(e)); NULL }),
    "Burde & Linden" = function()
      tryCatch(run_burde_analysis_fixed(),
               error = function(e) { message(conditionMessage(e)); NULL }),
    "Galiani et al." = function()
      tryCatch(run_galiani_analysis_fixed(),
               error = function(e) { message(conditionMessage(e)); NULL }),
    "Banerjee et al." = function()
      tryCatch(run_banerjee_analysis_fixed(),
               error = function(e) { message(conditionMessage(e)); NULL })
  )
  
  all_std <- list()
  for (study_name in names(study_runners)) {
    if (verbose) cat(sprintf("  Running %s...\n", study_name))
    raw <- study_runners[[study_name]]()
    if (!is.null(raw) && nrow(raw) > 0) {
      std <- standardize_cr_result(raw, study_name)
      if (!is.null(std)) all_std[[study_name]] <- std
    }
  }
  
  all_std <- all_std[!vapply(all_std, is.null, logical(1))]
  if (!length(all_std)) { warning("scenario6: all studies failed."); return(NULL) }
  
  combined <- dplyr::bind_rows(all_std)
  keep     <- c("Study","Instrument","n","r_xz","r_zu","beta_IV","plug_in",
                "CI_MCUB","CI_simple","Zero_MCUB","p_zero","verdict")
  combined <- combined[, intersect(keep, names(combined)), drop = FALSE]
  rownames(combined) <- NULL
  
  if (enrich) {
    if (verbose) cat("\n  Enriching with r_xz, r_zu, beta_IV...\n")
    combined <- enrich_s6(combined, verbose = verbose)
  }
  
  if (verbose) { cat("\n  Unified empirical results:\n"); print(combined) }
  combined
}

#================================================================
# SECTION 9 — enrich_s6()
#================================================================

enrich_s6 <- function(combined, verbose = TRUE) {
  studies <- .empirical_studies()
  
  for (st in studies) {
    if (verbose) cat(sprintf("    Enriching %s ...\n", st$name))
    raw <- tryCatch(st$loader(), error = function(e) NULL)
    if (is.null(raw)) {
      cat("      Data load failed — skipping enrichment for", st$name, "\n"); next
    }
    
    for (Z in st$Zs) {
      vars <- c(st$Y, st$X, Z, st$controls)
      d    <- raw[complete.cases(raw[, vars]), vars, drop = FALSE]
      fml  <- function(lhs) as.formula(paste(lhs,"~",paste(st$controls,collapse="+")))
      
      y    <- resid(lm(fml(st$Y), data = d))
      x    <- resid(lm(fml(st$X), data = d))
      z0   <- resid(lm(fml(Z),    data = d))
      zp   <- predict(lm(z0 ~ x + y))
      u_proxy <- resid(lm(y ~ x))
      
      rho_xz  <- round(cor(x, z0),       3)
      rho_zu  <- round(cor(zp, u_proxy), 3)
      beta_iv <- round(cov(zp, y) / cov(zp, x), 3)
      
      mask <- combined$Study == st$name & combined$Instrument == Z
      if (any(mask)) {
        combined$r_xz[mask]    <- rho_xz
        combined$r_zu[mask]    <- rho_zu
        combined$beta_IV[mask] <- beta_iv
      }
    }
  }
  combined
}

#================================================================
# SECTION 10 — Output utilities
#================================================================

save_results <- function(results, prefix = "cr_neurips", outdir = ".") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  map <- list(s1 = "s1_dvec_spectrum", s2 = "s2_mc_power",
              s3 = "s3_nonnormal",     s4 = "s4_rxu_sensitivity",
              s5 = "s5_ci_comparator", s6 = "s6_empirical")
  for (tag in names(map)) {
    obj <- results[[tag]]
    if (is.null(obj)) next
    if (is.list(obj) && !is.data.frame(obj)) obj <- obj$summary
    if (is.data.frame(obj) && nrow(obj) > 0) {
      path <- file.path(outdir, paste0(prefix, "_", map[[tag]], ".csv"))
      write.csv(obj, path, row.names = FALSE)
      cat(sprintf("  Saved %s\n", path))
    }
  }
  invisible(NULL)
}

latex_table <- function(df, caption = "CR test results", label = "tab:cr_results") {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      col <- gsub("\u2713", "\\\\checkmark", col)
      col <- gsub("\u00d7", "$\\\\times$",   col)
    }
    col
  })
  header <- paste(names(df), collapse = " & ")
  rows   <- apply(df, 1, function(r) paste(r, collapse = " & "))
  cat("\\begin{table}[h]\n")
  cat(sprintf("\\caption{%s}\n\\label{%s}\n", caption, label))
  cat("\\centering\n")
  cat(sprintf("\\begin{tabular}{%s}\n", paste(rep("l", ncol(df)), collapse = "")))
  cat("\\toprule\n")
  cat(header, "\\\\\n\\midrule\n")
  cat(paste(rows, collapse = " \\\\\n"), "\\\\\n")
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
  invisible(df)
}

print_summary <- function(results) {
  cat("\n", strrep("=", 64), "\n", sep = "")
  cat("DIAL / CR TEST  —  RESULTS SUMMARY\n")
  cat(strrep("=", 64), "\n", sep = "")
  
  if (!is.null(results$s1)) {
    s1 <- results$s1
    flips <- which(diff(as.integer(s1$verdict == "Valid")) != 0)
    if (length(flips))
      cat(sprintf("S1  Detection boundary: d in (%.2f, %.2f) | r_zu in (%.3f, %.3f)\n",
                  s1$d[flips[1]], s1$d[flips[1]+1],
                  s1$r_zu[flips[1]], s1$r_zu[flips[1]+1]))
  }
  
  if (!is.null(results$s2)) {
    s2 <- results$s2
    for (k in seq_len(nrow(s2)))
      cat(sprintf("S2  N = %4d: MCUB power = %.0f%% | simple power = %.0f%%\n",
                  s2$N[k], 100*s2$reject_MCUB[k], 100*s2$reject_simple[k]))
  }
  
  if (!is.null(results$s3)) {
    s3 <- if (is.list(results$s3) && !is.data.frame(results$s3))
      results$s3$summary else results$s3
    cat(sprintf("S3  Width inflation (diag/boot): %.2f – %.2f\n",
                min(s3$width_inflation, na.rm=TRUE), max(s3$width_inflation, na.rm=TRUE)))
  }
  
  if (!is.null(results$s4)) {
    s4  <- results$s4
    rob <- s4 %>%
      dplyr::filter(!is.na(range_type)) %>%
      tidyr::pivot_wider(id_cols = "Study", names_from = "range_type", values_from = "verdict") %>%
      dplyr::mutate(robust = canonical == symmetric)
    for (k in seq_len(nrow(rob)))
      cat(sprintf("S4  %-22s robust = %s  (%s canonical | %s symmetric)\n",
                  rob$Study[k], ifelse(rob$robust[k], "YES", "NO "),
                  rob$canonical[k], rob$symmetric[k]))
  }
  
  if (!is.null(results$s5)) {
    s5 <- results$s5
    cat(sprintf("S5  Width ratio range: %.3f – %.3f (MCUB narrower)\n",
                min(s5$width_ratio, na.rm=TRUE), max(s5$width_ratio, na.rm=TRUE)))
  }
  
  if (!is.null(results$s6)) {
    s6 <- results$s6
    valid   <- s6$Instrument[s6$verdict == "Valid"]
    invalid <- s6$Instrument[s6$verdict == "Invalid"]
    cat(sprintf("S6  Valid  : %s\n", paste(valid,   collapse = ", ")))
    cat(sprintf("S6  Invalid: %s\n", paste(invalid, collapse = ", ")))
  }
  
  cat(strrep("=", 64), "\n", sep = "")
  invisible(NULL)
}

#================================================================
# SECTION 11 — Master runner
#================================================================

run_all_scenarios <- function(
    run     = c("s1","s2","s3","s4","s5","s6"),
    alpha   = 0.05,
    seed    = 9058L,
    N       = 800L,
    n_sim   = 20L,
    save    = FALSE,
    outdir  = "cr_results",
    verbose = TRUE
) {
  check_sources_loaded()
  results <- list()
  
  safe <- function(tag, expr) {
    if (!tag %in% run) return(invisible(NULL))
    cat(sprintf("\n%s\n", strrep("=", 60)))
    tryCatch(expr,
             error   = function(e) { message("  ERROR in ", tag, ": ", conditionMessage(e)); NULL },
             warning = function(w) { message("  WARN  in ", tag, ": ", conditionMessage(w)); NULL })
  }
  
  results$s1 <- safe("s1", scenario1_dvec_spectrum(N = N, seed = seed,
                                                   rxu_range = c(-0.4, 0.4), verbose = verbose))
  results$s2 <- safe("s2", scenario2_mc_power(N_vec = c(500L, 800L, 2000L), n_sim = n_sim,
                                              rxu_range = c(-0.4, 0.4), seed0 = seed, verbose = verbose))
  results$s3 <- safe("s3", scenario3_nonnormal_stress(N = N, seed = seed,
                                                      rxu_range = c(-0.4, 0.4), verbose = verbose))
  results$s4 <- safe("s4", scenario4_rxu_sensitivity(alpha = alpha, verbose = verbose))
  results$s5 <- safe("s5", scenario5_ci_comparator(N = N, seed = seed,
                                                   rxu_range = c(-0.4, 0.4), verbose = verbose))
  results$s6 <- safe("s6", scenario6_multistudy_table(alpha = alpha, enrich = TRUE, verbose = verbose))
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  print_summary(results)
  
  if (save) save_results(results, outdir = outdir)
  
  cat("\nAccess results via: results$s1 ... results$s6\n\n")
  invisible(results)
}

#================================================================
# SECTION 12 — Quick-start reference
#
# DIAL_SKIP_AUTORUN <- TRUE
.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- sub('--file=', '', args[grep('--file=', args)])
    if (length(f) && nchar(f)) dirname(normalizePath(f))
    else getwd()
  }
)
.repo_root <- normalizePath(file.path(.script_dir, '..'))

source(file.path(.repo_root, 'foundation', 'Fin_sim3_clean.R'))
source(file.path(.repo_root, 'foundation', 'Fin_Empirical5_clean.R'))
source("DIAL_NEURIPS_CR_wrapper.R")
check_sources_loaded()
#
# # Phase 1 verification (run first after applying corrections):
source("DIAL_Phase1_runner.R")
#
# # Full pipeline:
res <- run_all_scenarios(n_sim = 50, save = TRUE)
#
# # Individual scenarios:
res_s4 <- scenario4_rxu_sensitivity()   # confirms Banerjee verdict flip
res_s6 <- scenario6_multistudy_table(enrich = TRUE)
#
# # LaTeX tables:
latex_table(res_s6, caption = "CR test results: empirical studies",
            label = "tab:cr_empirical")
#================================================================
