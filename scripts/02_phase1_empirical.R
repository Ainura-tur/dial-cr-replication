# ==============================================================================
# Dial_phase1_runner.R  -  Verification runner for Table 4 / Table 5 (Phase 1)
# ==============================================================================
#
# PURPOSE
#   Calls the four "phase 1" empirical pipeline functions in
#   Fin_Empirical5_clean.R with the corrected within-sign canonical
#   domains, plus the within-sign widening pairs needed for Table 5
#   (domain sensitivity). Writes a flat phase1_results.csv digest.
#
#   Phase 1 covers the four .dta-loaded studies:
#     - Di Tella & Schargrodsky (2013)   -- 2 instruments
#     - Galiani et al. (2011)            -- 1 instrument
#     - Banerjee et al. (2007)           -- 1 instrument
#     - Burde & Linden (2013)            -- 1 instrument + sensitivity
#
#   Phase 2 (MovieLens, Papke, NHANES, Card) is run separately by
#   DIAL_NeurIPS_fin_empirical5_phase2a.R.
#
# CANONICAL DOMAINS (matches Table 4 of the paper)
#   Di Tella judgeAlreadyUsedEM   c( 0.15,  0.25)
#   Di Tella percJudgeSentToEM    c( 0.15,  0.25)
#   Galiani  highnumber           c(-0.20,  0.00)
#   Banerjee instru               c(-0.60, -0.40)
#   Burde    buildschool          c( 0.00,  0.40)
#
# WITHIN-SIGN WIDENING (Table 5 robust column)
#   Di Tella                      c( 0.00,  0.80)
#   Galiani                       c(-0.80,  0.00)
#   Banerjee                      c(-0.80,  0.00)
#   Burde                         c( 0.00,  0.80)
#
# REQUIRES (source order is mandatory)
#   DIAL_SKIP_AUTORUN <- TRUE
#   source("Fin_sim3_clean.R")
#   source("Fin_Empirical5_clean.R")
#   source("Dial_phase1_runner.R")     # this file
#
# OUTPUT
#   phase1_results.csv  -- one row per (study x instrument x domain)
#   Returned (invisibly) as a tibble for further analysis.
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
    "Dial_phase1_runner.R: required functions not loaded:",
    paste(" -", .missing_p1, collapse = "\n"),
    "",
    "Source the dependencies first, in this exact order:",
    "",
    "    DIAL_SKIP_AUTORUN <- TRUE",
    "    source(\"Fin_sim3_clean.R\")",
    "    source(\"Fin_Empirical5_clean.R\")",
    "    source(\"Dial_phase1_runner.R\")",
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
      stringsAsFactors = FALSE
    ))
  }
  # Fin_Empirical5_clean.R returns 1 row per instrument; pick first.
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
  #                       Section 5.4 sign elicitation.
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
  # NOTE (Corollary 4): extending the domain across the sign boundary
  # is NOT a "within-sign robustness check" in the corollary's sense.
  # By Corollary 4(c) it can re-admit rho_DU = 0 and collapse coverage
  # for a valid instrument. We report it here for completeness as a
  # cross-sign sensitivity diagnostic, NOT as a verdict change.
  cat("  Cross-sign sensitivity (NOT a robustness check per Cor. 4): (-0.4, 0)\n")
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
  cat("[Dial_phase1_runner.R] DIAL_SKIP_AUTORUN = TRUE -- ",
      "call run_phase1_all() manually.\n", sep = "")
}