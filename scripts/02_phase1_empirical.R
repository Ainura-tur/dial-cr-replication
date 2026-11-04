# ==============================================================================
# 02_phase1_empirical.R  -  Phase 1 runner for tab:empirical_full and
#                           tab:domain_sens
# ==============================================================================
#
# PURPOSE
#   Calls the four "phase 1" empirical pipeline functions in
#   Fin_Empirical5_clean.R with the within-sign canonical domains, plus
#   the within-sign widening pairs needed for the domain-sensitivity
#   table. Writes a flat phase1_results.csv digest.
#
#   Tables and figures are referenced by LaTeX label, not by number, so
#   a float move in the paper cannot make these comments stale.
#
#   Phase 1 covers the four .dta-loaded studies:
#     - Di Tella & Schargrodsky (2013)   -- 2 instruments
#     - Galiani et al. (2011)            -- 1 instrument
#     - Banerjee et al. (2007)           -- 1 instrument
#     - Burde & Linden (2013)            -- 1 instrument + cross-sign diagnostic
#
#   Phase 2 (MovieLens, Papke, NHANES, Card) is run separately by
#   03_phase2_empirical.R.
#
# CANONICAL DOMAINS (tab:empirical_full)
#   Di Tella judgeAlreadyUsedEM   c( 0.15,  0.25)
#   Di Tella percJudgeSentToEM    c( 0.00,  0.25)   # starts at 0, not 0.15
#   Galiani  highnumber           c(-0.20,  0.00)
#   Banerjee instru               c(-0.60, -0.40)
#   Burde    buildschool          c( 0.00,  0.40)
#
# WITHIN-SIGN WIDENING (tab:domain_sens, robustness column)
#   Di Tella                      c( 0.00,  0.80)
#   Galiani                       c(-0.80,  0.00)
#   Banerjee                      c(-0.80,  0.00)
#   Burde                         c( 0.00,  0.80)
#
# REQUIRES (source order is mandatory)
#   DIAL_SKIP_AUTORUN <- TRUE
#   source("../foundation/Fin_sim3_clean.R")
#   source("../foundation/Fin_Empirical5_clean.R")
#   source("02_phase1_empirical.R")    # this file
#
# Paths assume cwd = <repo>/scripts/ (the standard entry point).
# When running via 00_setup.R, foundation functions are already loaded
# and this script runs as-is.
#
# BOOTSTRAP VALIDITY RULE (prop:finite_sample)
#   This script sources ../foundation/bootstrap_rule.R, which wraps
#   check_compatibility_simple(). Each MCUB call then also runs the
#   tau_V rule of Algorithm 1: B_rule = 500 bootstrap resamples, one MCUB
#   set per resample, and
#     p_hat_B = share of resamples whose MCUB set contains zero.
#   The instrument is Valid if p_hat_B >= tau_V = 0.95. The verdict is
#   strict (zero_tol = 0), like the single-sample column Zero_in_CI.
#   Each resample takes 1-4 s on one core, so the rule adds roughly
#   10-30 core-minutes per row; dial_boot_plan() runs the resamples in
#   parallel (up to 7 workers). Resample results are cached in
#   results/boot_rule/, so an interrupted run resumes where it stopped.
#   To skip the rule: options(dial.boot_rule = FALSE) before sourcing.
#
# RAW-DATA BOOTSTRAP
#   Every bootstrap resample redraws the rows of the analysis data and
#   redoes the residualise-and-project step: the rule's resamples, the
#   covariance bootstrap inside each resample's MCUB set, and the
#   covariance bootstrap behind the full-sample MCUB set (CI_MCUB). The
#   Phase 1 functions read their .dta files with read_dta(); the shim
#   below passes each call through unchanged and registers the analysis
#   specifications that go with the file, so the wrapper can find the raw
#   rows. options(dial.boot_unit = "rows") before sourcing restores the
#   earlier behaviour (resample the prepared (x, y, z) rows).
#
# OUTPUT
#   phase1_results.csv  -- one row per (study x instrument x domain)
#     columns p_hat_B, n_boot_ok and Valid_tau2 hold the bootstrap rule;
#     Zero_in_CI is the single-sample check (zero in the full-sample set);
#     cov_unit and boot_unit record the resampling unit ("raw")
#   Returned (invisibly) as a data frame for further analysis.
#
# RUN
#   In R, after the source() steps above:
#     res <- run_phase1_all()
# ==============================================================================

# Source guard -- emit a clear error if dependencies are not loaded.
.required_phase1_fns <- c(
  "run_ditella_analysis_fixed",
  "run_galiani_analysis_fixed",
  "run_banerjee_analysis_fixed",
  "run_burde_analysis_fixed",
  "run_burde_sensitivity_fixed",
  "check_compatibility_simple"
)

.missing_p1 <- setdiff(.required_phase1_fns, ls(envir = .GlobalEnv))
if (length(.missing_p1) > 0) {
  msg <- paste(
    "02_phase1_empirical.R: required functions not loaded:",
    paste(" -", .missing_p1, collapse = "\n"),
    "",
    "Source the dependencies first, in this exact order",
    "(assumes cwd = <repo>/scripts/):",
    "",
    "    DIAL_SKIP_AUTORUN <- TRUE",
    "    source(\"../foundation/Fin_sim3_clean.R\")",
    "    source(\"../foundation/Fin_Empirical5_clean.R\")",
    "    source(\"02_phase1_empirical.R\")",
    "",
    "Or simply run 00_setup.R once per R session before any script.",
    "",
    "The DIAL_SKIP_AUTORUN flag suppresses the auto-run blocks at the",
    "bottom of Fin_sim3_clean.R and Fin_Empirical5_clean.R so they",
    "load the function definitions without immediately starting a",
    "Monte Carlo loop.",
    sep = "\n"
  )
  stop(msg, call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
})

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;    filter    <- dplyr::filter
mutate    <- dplyr::mutate;    slice     <- dplyr::slice
recode    <- dplyr::recode;    rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange   <- dplyr::arrange;   count     <- dplyr::count
lag       <- dplyr::lag

# ---- Bootstrap validity rule (see header) ---------------------------------

.br_file <- file.path("..", "foundation", "bootstrap_rule.R")
if (!file.exists(.br_file) && exists("FOUNDATION_DIR"))
  .br_file <- file.path(FOUNDATION_DIR, "bootstrap_rule.R")
if (!file.exists(.br_file))
  stop("02_phase1_empirical.R: foundation/bootstrap_rule.R not found ",
       "(expected cwd = <repo>/scripts/).", call. = FALSE)
options(dial.zero_tol = 0)          # strict verdict in the empirical tables
source(.br_file)
dial_boot_plan()
rm(.br_file)

# Raw data for the raw-data bootstrap (the default dial.boot_unit = "raw";
# see foundation/bootstrap_rule.R). The Phase 1 functions read their .dta
# files with read_dta(). This shim passes each call through unchanged and
# registers the analysis specifications that go with the file, so every
# bootstrap can resample the raw rows and redo the residualise-and-project
# step.
.p1_raw_specs <- list(
  list(match = "Di%20Tella|Di Tella|Di_Tella",
       keep  = function(d) d[which(d$offendersPerJudge > 9), , drop = FALSE],
       Y = "recidivism", X = "electronicMonitoring",
       Z = c("judgeAlreadyUsedEM", "percJudgeSentToEM"),
       controls = c("mostSeriousCrime", "age", "ageSquared", "argentine",
                    "numberPreviousImprisonments", "judicialDistrict",
                    "yearOfImprisonment")),
  list(match = "Crime\\.dta",
       keep  = function(d) d[which(d$cohort > 1957 & d$cohort < 1963), , drop = FALSE],
       Y = "crimerate", X = "sm", Z = "highnumber",
       controls = c("cohort", "draftnumber", "navy")),
  list(match = "yld_sett_aug03\\.dta",
       keep  = function(d) d[which(d$year >= 1965 & d$phwht <= 1), , drop = FALSE],
       Y = "phwht", X = "p_nland", Z = "instru",
       controls = c("alt", "totrain", "so_black", "so_red", "so_all",
                    "lat", "coastal", "brule1", "year")),
  list(match = "afgan\\.dta",
       keep  = function(d) d,
       Y = "testscore", X = "enrolled", Z = "buildschool",
       controls = c("headchild", "age", "yrsvill", "farsi", "tajik", "farmers",
                    "agehead", "educhead", "nhh", "land", "sheep",
                    "distschool", "chagcharan")))
if (identical(getOption("dial.boot_unit"), "raw")) {
  .read_dta_prev <- if (exists("read_dta", envir = .GlobalEnv, inherits = FALSE))
    get("read_dta", envir = .GlobalEnv) else haven::read_dta
  read_dta <- function(file, ...) {
    d <- .read_dta_prev(file, ...)
    for (s in .p1_raw_specs) if (grepl(s$match, file)) {
      dd <- as.data.frame(s$keep(as.data.frame(d)))
      options(dial.boot_raw = lapply(s$Z, function(z)
        dial_raw_spec(dd, s$Y, s$X, z, s$controls)))
    }
    d
  }
}

# ---- Helper: standardise a single result row ------------------------------

.standardise_p1 <- function(res, study, instrument, domain_label, rxu_lo, rxu_hi) {
  if (is.null(res) || nrow(res) == 0L) {
    return(data.frame(
      study = study, instrument = instrument,
      domain_label = domain_label,
      rxu_lo = rxu_lo, rxu_hi = rxu_hi,
      n = NA_integer_,
      r_xz = NA_real_, beta_IV = NA_real_,
      plug_in = NA_character_, CI_MCUB = NA_character_,
      Zero_in_CI = NA_character_, p_zero = NA_real_,
      cov_unit = NA_character_,
      p_hat_B = NA_real_, n_boot_ok = NA_integer_,
      boot_unit = NA_character_, Valid_tau2 = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  # Fin_Empirical5_clean.R returns 1 row per instrument; pick first.
  # CI_MCUB is the MCUB set (CI_Bei); plug_in is the plug-in identified
  # set. The paper reports CI_MCUB.
  r <- res[1L, , drop = FALSE]
  data.frame(
    study        = study,
    instrument   = instrument,
    domain_label = domain_label,
    rxu_lo       = rxu_lo,
    rxu_hi       = rxu_hi,
    n            = if ("n" %in% names(r)) r$n else NA_integer_,
    r_xz         = if ("r_xz" %in% names(r)) r$r_xz else NA_real_,
    beta_IV      = if ("beta_IV" %in% names(r)) r$beta_IV else NA_real_,
    plug_in      = if ("plug_in" %in% names(r)) r$plug_in else NA_character_,
    CI_MCUB      = if ("CI_Bei"  %in% names(r)) r$CI_Bei  else NA_character_,
    Zero_in_CI   = if ("Zero_in_CI" %in% names(r)) r$Zero_in_CI else NA_character_,
    p_zero       = if ("p_zero" %in% names(r)) r$p_zero else NA_real_,
    # Resampling unit of the full-sample covariance ("raw" = preparation
    # redone on every draw; see foundation/bootstrap_rule.R)
    cov_unit     = if ("cov_unit" %in% names(r)) r$cov_unit else NA_character_,
    # Bootstrap validity rule (bootstrap_rule.R): verdict of Algorithm 1
    p_hat_B      = if ("p_hat_B" %in% names(r)) r$p_hat_B else NA_real_,
    n_boot_ok    = if ("n_boot_ok" %in% names(r)) r$n_boot_ok else NA_integer_,
    boot_unit    = if ("boot_unit" %in% names(r)) r$boot_unit else NA_character_,
    Valid_tau2   = if ("Valid_tau2" %in% names(r)) r$Valid_tau2 else NA_character_,
    stringsAsFactors = FALSE
  )
}

# ---- Per-study runners ----------------------------------------------------

.run_ditella_pair <- function() {
  cat("\n--- Di Tella & Schargrodsky (2013) ---\n")
  out <- list()

  # Canonical domains DIFFER by instrument:
  #   judgeAlreadyUsedEM : (0.15, 0.25)  -- judge-history measure pinned
  #                       to the upper part of the positive half by
  #                       sign elicitation (sec:domain).
  #   percJudgeSentToEM  : (0.00, 0.25)  -- judge-preference measure
  #                       has wider plausible support starting at zero.
  #
  # Earlier versions of this runner applied (0.15, 0.25) uniformly to
  # both instruments; that is corrected here. The wide canonical
  # (within-sign widening) remains (0, 0.8) for both.

  cat("  judgeAlreadyUsedEM canonical: (0.15, 0.25)\n")
  res_jaue <- tryCatch(
    run_ditella_analysis_fixed(
      instruments = c("judgeAlreadyUsedEM"),
      rxu_range   = c(0.15, 0.25)
    ),
    error = function(e) { cat("    FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(res_jaue)) {
    out[[length(out) + 1L]] <- .standardise_p1(
      res_jaue, study = "Di Tella & Schargrodsky",
      instrument = "judgeAlreadyUsedEM",
      domain_label = "canonical", rxu_lo = 0.15, rxu_hi = 0.25
    )
  }

  cat("  percJudgeSentToEM canonical: (0, 0.25)\n")
  res_pjse <- tryCatch(
    run_ditella_analysis_fixed(
      instruments = c("percJudgeSentToEM"),
      rxu_range   = c(0.0, 0.25)
    ),
    error = function(e) { cat("    FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(res_pjse)) {
    out[[length(out) + 1L]] <- .standardise_p1(
      res_pjse, study = "Di Tella & Schargrodsky",
      instrument = "percJudgeSentToEM",
      domain_label = "canonical", rxu_lo = 0.0, rxu_hi = 0.25
    )
  }

  # Within-sign widening: (0, 0.8) for both
  cat("  Wide canonical (both): (0, 0.8)\n")
  res_wide <- tryCatch(
    run_ditella_analysis_fixed(
      instruments = c("judgeAlreadyUsedEM", "percJudgeSentToEM"),
      rxu_range   = c(0.0, 0.8)
    ),
    error = function(e) { cat("    FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(res_wide)) {
    instr_names <- c("judgeAlreadyUsedEM", "percJudgeSentToEM")
    for (j in seq_along(instr_names)) {
      out[[length(out) + 1L]] <- .standardise_p1(
        res_wide[j, , drop = FALSE],
        study = "Di Tella & Schargrodsky", instrument = instr_names[j],
        domain_label = "wide_canonical", rxu_lo = 0.0, rxu_hi = 0.8
      )
    }
  }

  do.call(rbind, out)
}

.run_galiani <- function() {
  cat("\n--- Galiani et al. (2011) ---\n")
  out <- list()
  cat("  Canonical: (-0.2, 0.0)\n")
  res1 <- tryCatch(run_galiani_analysis_fixed(rxu_range = c(-0.2, 0.0)),
                   error = function(e) NULL)
  out[[1L]] <- .standardise_p1(res1, "Galiani et al.", "highnumber",
                               "canonical", -0.2, 0.0)

  cat("  Wide canonical: (-0.8, 0.0)\n")
  res2 <- tryCatch(run_galiani_analysis_fixed(rxu_range = c(-0.8, 0.0)),
                   error = function(e) NULL)
  out[[2L]] <- .standardise_p1(res2, "Galiani et al.", "highnumber",
                               "wide_canonical", -0.8, 0.0)
  do.call(rbind, out)
}

.run_banerjee <- function() {
  cat("\n--- Banerjee et al. (2007) ---\n")
  out <- list()
  cat("  Canonical: (-0.6, -0.4)\n")
  res1 <- tryCatch(run_banerjee_analysis_fixed(rxu_range = c(-0.6, -0.4)),
                   error = function(e) NULL)
  out[[1L]] <- .standardise_p1(res1, "Banerjee et al.", "instru",
                               "canonical", -0.6, -0.4)

  cat("  Wide canonical: (-0.8, 0.0)\n")
  res2 <- tryCatch(run_banerjee_analysis_fixed(rxu_range = c(-0.8, 0.0)),
                   error = function(e) NULL)
  out[[2L]] <- .standardise_p1(res2, "Banerjee et al.", "instru",
                               "wide_canonical", -0.8, 0.0)
  do.call(rbind, out)
}

.run_burde <- function() {
  cat("\n--- Burde & Linden (2013) ---\n")
  out <- list()
  cat("  Canonical: (0, 0.4)\n")
  res1 <- tryCatch(run_burde_analysis_fixed(rxu_range = c(0.0, 0.4)),
                   error = function(e) NULL)
  out[[1L]] <- .standardise_p1(res1, "Burde & Linden", "buildschool",
                               "canonical", 0.0, 0.4)

  cat("  Wide canonical: (0, 0.8)\n")
  res2 <- tryCatch(run_burde_analysis_fixed(rxu_range = c(0.0, 0.8)),
                   error = function(e) NULL)
  out[[2L]] <- .standardise_p1(res2, "Burde & Linden", "buildschool",
                               "wide_canonical", 0.0, 0.8)

  # Burde cross-sign sensitivity check.
  # NOTE (cor:domain, part c): extending the domain across the sign
  # boundary is NOT a within-sign robustness check. It can re-admit
  # rho_DU = 0 and collapse coverage for a valid instrument. We report
  # it here for completeness as a cross-sign diagnostic, NOT as a
  # verdict change.
  #
  # DOWNSTREAM FILTER: 08_figures.R must EXCLUDE rows where
  #   domain_label == "cross_sign_sensitivity"
  # from fig:domain_sens (domain-sensitivity bars) and
  # fig:diagnostic_space. Only "canonical" and "wide_canonical" rows
  # belong in those figures.
  cat("  Cross-sign sensitivity (NOT a robustness check per cor:domain): (-0.4, 0)\n")
  res3 <- tryCatch(run_burde_sensitivity_fixed(),
                   error = function(e) NULL)
  out[[3L]] <- .standardise_p1(res3, "Burde & Linden", "buildschool",
                               "cross_sign_sensitivity", -0.4, 0.0)
  do.call(rbind, out)
}

# ---- Public entry point ---------------------------------------------------

run_phase1_all <- function(out_csv = "phase1_results.csv", verbose = TRUE) {
  t0 <- Sys.time()

  rows <- list(
    .run_ditella_pair(),
    .run_galiani(),
    .run_banerjee(),
    .run_burde()
  )
  combined <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  rownames(combined) <- NULL

  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("PHASE 1 SUMMARY\n")
  cat(strrep("=", 78), "\n", sep = "")
  print(combined, row.names = FALSE)

  tryCatch({
    utils::write.csv(combined, out_csv, row.names = FALSE)
    cat(sprintf("\nWrote %d rows to %s\n", nrow(combined), out_csv))
  }, error = function(e) cat("\nCSV write failed:", conditionMessage(e), "\n"))

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("Elapsed: %.1f s\n", elapsed))

  invisible(combined)
}

# ---- Auto-run guard (mirrors Fin_*_clean.R convention) --------------------

.dial_autorun_p1 <- !exists("DIAL_SKIP_AUTORUN") || !isTRUE(DIAL_SKIP_AUTORUN)
if (.dial_autorun_p1) {
  res_phase1 <- run_phase1_all()
} else {
  cat("[02_phase1_empirical.R] DIAL_SKIP_AUTORUN = TRUE -- ",
      "call run_phase1_all() manually.\n", sep = "")
}
