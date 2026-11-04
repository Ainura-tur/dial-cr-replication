# ==============================================================================
# 04_synthetic_benchmark.R   --   tab:synthetic (Synthetic Benchmark, CR vs HSIC)
# ------------------------------------------------------------------------------
# Tables are referenced by LaTeX label, not by number.
#
# Reproduces tab:synthetic under the Jones-Pewsey DGP with the per-row
# correctly-elicited domain protocol (cor:domain, practical protocol):
#
#   For each instrument with target cor(z, u) = c:
#     - if c < 0:  draw from a lambda < 0 DGP (so true rho_xu < 0),
#                  test under D = (-0.8, 0)
#     - if c > 0:  draw from a lambda > 0 DGP (so true rho_xu > 0),
#                  test under D = (0, 0.8)
#     - if c = 0:  test the valid instrument under D = (0, 0.8)
#                  (single-sign convention; results are identical under
#                   D = (-0.8, 0) up to sampling noise).
#
# This is single-sweep per row using the elicited sign that the
# practitioner would commit to given knowledge of the DGP's true rho_xu.
# Avoids the dual-sweep AND-rule failure: by cor:domain part (c), testing
# under the wrong sign collapses coverage even for valid instruments,
# so two sweeps cannot both pass for any nontrivial cor(z, u).
#
# DGP (sourced from foundation/jp_dgp.R)
#   x ~ J-P (sinh-arcsinh transform of normal, kappa=0.9 + Skew-Normal(1.2))
#   u = lambda (x - mean(x)) + epsilon + v_perp_x
#   y = beta * x + u, beta = 2
#   z_j: 11 instruments at known cor(z, u) in {-0.20, -0.16, ..., 0.20},
#         cor(z, x) ~ 0.25
#   z_het: one extra valid instrument with heteroskedastic errors, drawn
#         from its own sample (see "EXTRA ROW" below). cor(z, u) = 0 but
#         Var(u | z) depends on z, so HSIC should reject it while the CR
#         test accepts it.
#
# VERDICT
#   Tolerance-based zero-in-MCUB on CI_Bei.
#   Two-sided tolerance: zero is "in MCUB" iff
#     (ci_lo <= ZERO_TOL) AND (ci_hi >= -ZERO_TOL)
#   so the test handles positive-D and negative-D MCUBs symmetrically.
#   ZERO_TOL = 0.05.
#
#   Bootstrap validity rule (prop:finite_sample).
#   The DIAL column applies the tau_V rule of Algorithm 1 through
#   foundation/bootstrap_rule.R: B_RULE = 500 bootstrap resamples, one
#   MCUB set per resample, and
#     p_hat_B = share of resamples whose set contains zero (within ZERO_TOL)
#   A (valid) iff p_hat_B >= TAU2 = 0.95, else C. The single-sample check
#   above is kept in the column DIAL_single. The resamples run in parallel
#   (dial_boot_plan) and are cached in results/boot_rule/.
#
#   The verdict is read from the Valid_tau2 column of bootstrap_rule.R.
#   bootstrap_rule.R computes Valid_tau2 from the unrounded share but
#   returns p_hat_B rounded to three decimals, so comparing the rounded
#   p_hat_B with TAU2 could pass a share just below the threshold
#   (474/499 = 0.9499 prints as 0.950). The output also records n_boot_ok,
#   the number of resamples whose MCUB set could be computed, and
#   boot_unit, the resampling unit ("raw" or "rows").
#
#   Raw-data bootstrap: every bootstrap resample redraws the
#   simulated (y, x, z) rows and redoes the centring and projection -- the
#   rule's resamples, the covariance bootstrap inside each resample's MCUB
#   set, and the covariance bootstrap behind the full-sample set (ci_lo,
#   ci_hi). options(dial.boot_unit = "rows") before sourcing restores the
#   earlier behaviour (resample the prepared rows).
#
# SAMPLE SIZE
#   Every simulated sample has N = 800 rows, the sample size behind
#   tab:synthetic. The script does not source Fin_sim3_clean.R (it needs
#   nothing from that file, and its top-level code sets N = 800), and it
#   stops if a sourced file changes N.
#
# COMPARISON BASELINE
#   HSIC test of independence between z and the IV residual e_hat = y - beta_IV * x.
#
# DEPENDENCIES (the script sources these itself, in this order)
#   foundation/Fin_Empirical5_clean.R   (with DIAL_SKIP_AUTORUN <- TRUE)
#   foundation/jp_dgp.R
#   foundation/bootstrap_rule.R
#
# Foundation path:
#   This script lives in <repo>/scripts/ and the foundation files live
#   in <repo>/foundation/. The relative path "../foundation" assumes the
#   working directory is <repo>/scripts/ (the standard entry point).
#   The resolution block below also supports being sourced from <repo>/.
#
# OUTPUT
#   table1_results.csv  -- one row per instrument (file name kept for
#                          continuity with earlier runs); n_boot_ok and
#                          boot_unit are the last two columns
#   Console: LaTeX rows for tab:synthetic in the paper's column order,
#            ready to paste between the header \midrule and \bottomrule,
#            and the table-note lines for clipped MCUB endpoints and for
#            the sample correlation of Z_het with U.
#
# RUN
#   Rscript scripts/04_synthetic_benchmark.R
# ==============================================================================

rm(list = ls())

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org/")
}
install_if_missing(c("dHSIC", "ivreg", "dplyr"))
suppressPackageStartupMessages({
  library(dHSIC); library(ivreg)
})

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;    filter    <- dplyr::filter
mutate    <- dplyr::mutate;    slice     <- dplyr::slice
recode    <- dplyr::recode;    rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange   <- dplyr::arrange;   count     <- dplyr::count
lag       <- dplyr::lag

# ==============================================================================
# CONFIG
# ==============================================================================

N             <- 800L          # rows per simulated sample (see SAMPLE SIZE)
TARGET_COR_XU <- 0.30          # |rho_xu| target; sign set per-DGP
COR_ZX        <- 0.25
ALPHA         <- 0.10
ZERO_TOL      <- 0.05
RXU_POS       <- c(0, 0.80)    # used when lambda > 0 (cor_zu >= 0)
RXU_NEG       <- c(-0.80, 0)   # used when lambda < 0 (cor_zu  < 0)
MASTER_SEED   <- 1234L
B_RULE        <- 500L          # bootstrap resamples for the tau_V rule
TAU2          <- 0.95          # validity threshold tau_V

# ==============================================================================
# FOUNDATION SOURCE
# ==============================================================================

# Resolve the foundation directory robustly:
#   - If 00_setup.R has already set DIAL_FOUNDATION_DIR, use that.
#   - Otherwise try ../foundation (cwd = <repo>/scripts) then foundation
#     (cwd = <repo>) in that order.
.resolve_foundation_dir <- function() {
  if (exists("DIAL_FOUNDATION_DIR", envir = .GlobalEnv))
    return(get("DIAL_FOUNDATION_DIR", envir = .GlobalEnv))
  for (candidate in c("../foundation", "foundation")) {
    if (dir.exists(candidate))
      return(normalizePath(candidate))
  }
  stop("Could not locate the 'foundation' directory. ",
       "Set DIAL_FOUNDATION_DIR or run from <repo>/scripts/ or <repo>/.")
}
FOUNDATION_DIR <- .resolve_foundation_dir()

# Fin_sim3_clean.R is not sourced: nothing here needs it, and its
# top-level code sets N = 800 (see SAMPLE SIZE).
FOUNDATION_FILES <- normalizePath(file.path(
  FOUNDATION_DIR,
  c("Fin_Empirical5_clean.R",
    "jp_dgp.R")
))

.N_config <- N    # checked after sourcing (see SAMPLE SIZE)
assign("DIAL_SKIP_AUTORUN", TRUE, envir = .GlobalEnv)

for (.f in FOUNDATION_FILES) {
  if (!file.exists(.f)) stop("Foundation file not found: ", .f)
  invisible(capture.output(source(.f, local = FALSE)))
}

stopifnot(exists("check_compatibility_simple"),
          exists("draw_jp_sample"),
          exists("calibrate_lambda"),
          exists("JP_INSTRUMENT_GRID"))

# Bootstrap validity rule: wraps check_compatibility_simple() (see VERDICT)
options(dial.boot_B = B_RULE, dial.tau2 = TAU2, dial.zero_tol = ZERO_TOL)
source(file.path(FOUNDATION_DIR, "bootstrap_rule.R"))
dial_boot_plan()

# No sourced file may change the sample size (see SAMPLE SIZE)
if (!identical(N, .N_config))
  stop("A sourced file changed N from ", .N_config, " to ", N, ".")

cat("================================================================================\n")
cat("DIAL synthetic benchmark (tab:synthetic): J-P DGP, per-row correctly-elicited domain\n")
cat("================================================================================\n")
cat(sprintf("  n = %d, |target cor(x,u)| = %.2f, cor(z,x) = %.2f\n",
            N, TARGET_COR_XU, COR_ZX))
cat(sprintf("  D for lambda>0: (%.2f, %.2f)   D for lambda<0: (%.2f, %.2f)\n",
            RXU_POS[1], RXU_POS[2], RXU_NEG[1], RXU_NEG[2]))
cat(sprintf("  alpha = %.2f, ZERO_TOL = %.2f, master seed = %d\n",
            ALPHA, ZERO_TOL, MASTER_SEED))
cat(sprintf("  Instrument grid: cor(z,u) in {%s}\n",
            paste(sprintf("%+.2f", JP_INSTRUMENT_GRID), collapse = ", ")))
cat("\n")

# ==============================================================================
# CALIBRATE LAMBDA FOR BOTH SIGNS
# ==============================================================================

cat("Calibrating lambda > 0 for cor(x, u) ~ +", TARGET_COR_XU, "...\n", sep = "")
cal_pos <- calibrate_lambda(n = N, target_cor_xu = TARGET_COR_XU,
                            cor_zx_default = COR_ZX,
                            seed = MASTER_SEED, sign = "positive",
                            verbose = FALSE)
cat(sprintf("  lambda_+ = %+.4f, realized cor(x, u) = %+.4f\n",
            cal_pos$lambda, cal_pos$realized))

cat("Calibrating lambda < 0 for cor(x, u) ~ -", TARGET_COR_XU, "...\n", sep = "")
cal_neg <- calibrate_lambda(n = N, target_cor_xu = TARGET_COR_XU,
                            cor_zx_default = COR_ZX,
                            seed = MASTER_SEED, sign = "negative",
                            verbose = FALSE)
cat(sprintf("  lambda_- = %+.4f, realized cor(x, u) = %+.4f\n\n",
            cal_neg$lambda, cal_neg$realized))

# ==============================================================================
# DRAW MASTER SAMPLES (one per sign)
# ==============================================================================

cat("Drawing master samples (one per sign of lambda)...\n")
sim_pos <- draw_jp_sample(n = N, lambda = cal_pos$lambda,
                          seed = MASTER_SEED + 1L,
                          cor_zu_grid = JP_INSTRUMENT_GRID,
                          cor_zx_default = COR_ZX)
sim_neg <- draw_jp_sample(n = N, lambda = cal_neg$lambda,
                          seed = MASTER_SEED + 2L,
                          cor_zu_grid = JP_INSTRUMENT_GRID,
                          cor_zx_default = COR_ZX)

cat(sprintf("  lambda_+ sample: realized cor(x, u) = %+.4f\n",
            sim_pos$realized_cor_xu))
cat(sprintf("  lambda_- sample: realized cor(x, u) = %+.4f\n\n",
            sim_neg$realized_cor_xu))

# ==============================================================================
# Per-instrument pipeline
# ==============================================================================

# One value from the check_compatibility_simple() output, or `na` if the
# column is missing or NA (the rule's columns exist only when the rule ran).
.cr_col <- function(cr, name, na) {
  if (is.null(cr) || !NROW(cr) || is.null(cr[[name]])) return(na)
  v <- cr[[name]][1]
  if (is.na(v)) na else v
}

run_one_instrument <- function(j) {
  cor_target <- JP_INSTRUMENT_GRID[j]

  # Pick the sign-appropriate sample and domain.
  if (cor_target < 0) {
    sim     <- sim_neg
    rxu     <- RXU_NEG
    sign    <- "neg"
    lambda  <- cal_neg$lambda
  } else {
    sim     <- sim_pos
    rxu     <- RXU_POS
    sign    <- "pos"
    lambda  <- cal_pos$lambda
  }

  z       <- sim$z[, j]
  df_full <- data.frame(y = sim$y, x = sim$x, z = z)

  iv_mod  <- ivreg::ivreg(y ~ x | z, data = df_full)
  beta_iv <- as.numeric(stats::coef(iv_mod)["x"])

  y_r <- sim$y - mean(sim$y)
  x_r <- sim$x - mean(sim$x)
  z0  <- z - mean(z)
  z_p <- stats::predict(stats::lm(z0 ~ x_r + y_r))
  df_cr <- data.frame(x = x_r, y = y_r, z = z_p)
  # Raw data for the raw-data bootstrap (bootstrap_rule.R)
  df_cr <- dial_attach_raw(df_cr, data.frame(Y = sim$y, X = sim$x, Z = z),
                           "Y", "X", "Z")

  cr <- tryCatch({
    invisible(capture.output(
      out <- check_compatibility_simple(df_cr, i = j + 100L,
                                        alpha = ALPHA,
                                        rxu_range = rxu)
    ))
    out
  }, error = function(e) NULL)

  if (is.null(cr) || nrow(cr) == 0L) {
    return(data.frame(
      Instrument    = sprintf("z_%+0.2f", cor_target),
      cor_zu_target = cor_target,
      cor_zu_real   = sim$realized_cor_zu[j],
      DGP_sign      = sign,
      lambda        = round(lambda, 3),
      D_lo          = rxu[1], D_hi = rxu[2],
      beta_IV       = round(beta_iv, 3),
      ci_lo         = NA_real_, ci_hi = NA_real_,
      MCUB          = NA_character_,
      p_hat_B       = NA_real_,
      DIAL          = "?",
      DIAL_single   = "?",
      HSIC_p        = NA_real_,
      HSIC_verdict  = "?",
      HSIC_agrees   = "?",
      n_boot_ok     = NA_integer_,
      boot_unit     = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  ci_bei <- as.character(cr$CI_Bei[1])
  m <- regmatches(ci_bei,
                  regexec("\\[\\s*(-?[0-9.]+)\\s*,\\s*(-?[0-9.]+)\\s*\\]",
                          ci_bei))[[1]]
  ci_lo <- if (length(m) >= 3L) as.numeric(m[2L]) else NA_real_
  ci_hi <- if (length(m) >= 3L) as.numeric(m[3L]) else NA_real_

  # Single-sample check. Two-sided tolerance: handles positive and
  # negative D symmetrically.
  zero_in_tol <- !is.na(ci_lo) && !is.na(ci_hi) &&
    (ci_lo <= ZERO_TOL) && (ci_hi >= -ZERO_TOL)
  dial_single <- if (is.na(ci_lo) || is.na(ci_hi)) "?"
  else if (zero_in_tol) "A" else "C"
  # Bootstrap validity rule (Algorithm 1, Stage 2). The verdict comes from
  # Valid_tau2, which uses the unrounded share; p_hat_B is rounded to three
  # decimals for display only (see VERDICT in the header).
  p_hat_B    <- as.numeric(.cr_col(cr, "p_hat_B", NA_real_))
  n_boot_ok  <- as.integer(.cr_col(cr, "n_boot_ok", NA_integer_))
  boot_unit  <- as.character(.cr_col(cr, "boot_unit", NA_character_))
  valid_tau2 <- as.character(.cr_col(cr, "Valid_tau2", NA_character_))
  cr_verdict <- if (is.na(valid_tau2)) dial_single
  else if (valid_tau2 == "Valid") "A" else "C"

  e_hat <- sim$y - beta_iv * sim$x
  hsic_res <- tryCatch(
    dHSIC::dhsic.test(list(z, e_hat), method = "gamma", B = 100),
    error = function(e) NULL
  )
  hsic_p <- if (is.null(hsic_res)) NA_real_ else as.numeric(hsic_res$p.value)
  hsic_verdict <- if (is.na(hsic_p)) "?"
  else if (hsic_p >= 0.05) "A" else "C"
  hsic_agrees  <- if (cr_verdict == hsic_verdict) "Yes" else "No"

  cat(sprintf(
    "  z_%+0.2f  sign=%s  D=[%+.2f,%+.2f]  beta_IV=%+.3f  CI_Bei=[%+.3f,%+.3f]  p_hat_B=%.3f (%d/%d ok)  -> %s (single %s)   HSIC=%.3f (%s)\n",
    cor_target, sign, rxu[1], rxu[2], beta_iv, ci_lo, ci_hi, p_hat_B,
    n_boot_ok, B_RULE, cr_verdict, dial_single, hsic_p, hsic_verdict
  ))

  # Guard against negative zero in display (R's sprintf("%+0.3f", -1e-17)
  # otherwise prints "-0.000"). Snap values within 1e-9 of zero to +0.
  .snap0 <- function(x) ifelse(abs(x) < 1e-9, 0, x)

  data.frame(
    Instrument    = sprintf("z_%+0.2f", cor_target),
    cor_zu_target = .snap0(cor_target),
    cor_zu_real   = .snap0(round(sim$realized_cor_zu[j], 3)),
    DGP_sign      = sign,
    lambda        = round(lambda, 3),
    D_lo          = rxu[1], D_hi = rxu[2],
    beta_IV       = round(beta_iv, 3),
    ci_lo         = round(ci_lo, 3),
    ci_hi         = round(ci_hi, 3),
    MCUB          = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi),
    p_hat_B       = round(p_hat_B, 3),
    DIAL          = cr_verdict,
    DIAL_single   = dial_single,
    HSIC_p        = round(hsic_p, 3),
    HSIC_verdict  = hsic_verdict,
    HSIC_agrees   = hsic_agrees,
    n_boot_ok     = n_boot_ok,
    boot_unit     = boot_unit,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# EXTRA ROW: valid instrument with heteroskedastic errors (z_het)
# ------------------------------------------------------------------------------
# The eleven grid instruments test whether HSIC detects linear leakage. This
# row tests the other direction of the CR-vs-HSIC comparison: a valid
# instrument that is not independent of the structural error.
#
#   z  ~ J-P (same transform as x), standardised
#   u0 ~ Uniform(-1, 1) + J-P sinh-arcsinh draw, standardised, independent of z
#   u  = s(z) * u0,  s(z) = sqrt((1 + z^2) / 2),  so E[s(z)^2] = 1
#   v  = r * u0 + sqrt(1 - r^2) * w,  w ~ N(0, 1)   (first-stage error)
#   x  = pi * z + v,  pi = cor_zx / sqrt(1 - cor_zx^2)  so cor(z, x) = cor_zx
#   y  = beta * x + u,  beta = 2
#
# E[u | z] = s(z) E[u0] = 0, so z is valid (cor(z, u) = 0), but
# Var(u | z) = s(z)^2 changes with z, so z and u are dependent. The
# endogeneity loading r is set by uniroot so that cor(x, u) = TARGET_COR_XU,
# which is positive, so the elicited domain is RXU_POS.
# ==============================================================================

# >>> HET-ROW
draw_het_sample <- function(n, seed, cor_zx = COR_ZX,
                            target_cor_xu = TARGET_COR_XU, beta = 2) {
  set.seed(seed)
  .std <- function(a) (a - mean(a)) / stats::sd(a)
  z  <- .std(5 * jp_h_transform(stats::rnorm(n), zeta = 0, kappa = 0.9) +
               jp_skew_normal(n, zeta = 0, kappa = 1.2))
  u0 <- .std(stats::runif(n, -1, 1) + jp_skew_normal(n, zeta = 0, kappa = 1.2))
  w  <- stats::rnorm(n)
  u  <- sqrt((1 + z^2) / 2) * u0
  pi_z   <- cor_zx / sqrt(1 - cor_zx^2)
  make_x <- function(r) pi_z * z + r * u0 + sqrt(1 - r^2) * w
  r <- stats::uniroot(function(r) stats::cor(make_x(r), u) - target_cor_xu,
                      c(0, 0.95), tol = 1e-8)$root
  x <- make_x(r)
  list(x = x, y = beta * x + u, z = z, u = u, r = r,
       realized_cor_zu = stats::cor(z, u),
       realized_cor_xu = stats::cor(x, u),
       realized_cor_zx = stats::cor(z, x))
}

run_het_instrument <- function() {
  sim <- draw_het_sample(N, seed = MASTER_SEED + 3L)
  rxu <- RXU_POS

  df_full <- data.frame(y = sim$y, x = sim$x, z = sim$z)
  iv_mod  <- ivreg::ivreg(y ~ x | z, data = df_full)
  beta_iv <- as.numeric(stats::coef(iv_mod)["x"])

  y_r <- sim$y - mean(sim$y)
  x_r <- sim$x - mean(sim$x)
  z0  <- sim$z - mean(sim$z)
  z_p <- stats::predict(stats::lm(z0 ~ x_r + y_r))
  df_cr <- data.frame(x = x_r, y = y_r, z = z_p)
  # Raw data for the raw-data bootstrap (bootstrap_rule.R)
  df_cr <- dial_attach_raw(df_cr, data.frame(Y = sim$y, X = sim$x, Z = sim$z),
                           "Y", "X", "Z")

  cr <- tryCatch({
    invisible(capture.output(
      out <- check_compatibility_simple(df_cr, i = 200L,
                                        alpha = ALPHA,
                                        rxu_range = rxu)
    ))
    out
  }, error = function(e) NULL)

  ci_lo <- ci_hi <- NA_real_
  if (!is.null(cr) && nrow(cr) > 0L) {
    m <- regmatches(as.character(cr$CI_Bei[1]),
                    regexec("\\[\\s*(-?[0-9.]+)\\s*,\\s*(-?[0-9.]+)\\s*\\]",
                            as.character(cr$CI_Bei[1])))[[1]]
    if (length(m) >= 3L) { ci_lo <- as.numeric(m[2L]); ci_hi <- as.numeric(m[3L]) }
  }
  zero_in_tol <- !is.na(ci_lo) && !is.na(ci_hi) &&
    (ci_lo <= ZERO_TOL) && (ci_hi >= -ZERO_TOL)
  dial_single <- if (is.na(ci_lo) || is.na(ci_hi)) "?"
  else if (zero_in_tol) "A" else "C"
  # Bootstrap validity rule (Algorithm 1, Stage 2); verdict from Valid_tau2,
  # as in run_one_instrument()
  p_hat_B    <- as.numeric(.cr_col(cr, "p_hat_B", NA_real_))
  n_boot_ok  <- as.integer(.cr_col(cr, "n_boot_ok", NA_integer_))
  boot_unit  <- as.character(.cr_col(cr, "boot_unit", NA_character_))
  valid_tau2 <- as.character(.cr_col(cr, "Valid_tau2", NA_character_))
  cr_verdict <- if (is.na(valid_tau2)) dial_single
  else if (valid_tau2 == "Valid") "A" else "C"

  e_hat <- sim$y - beta_iv * sim$x
  hsic_res <- tryCatch(
    dHSIC::dhsic.test(list(sim$z, e_hat), method = "gamma", B = 100),
    error = function(e) NULL
  )
  hsic_p <- if (is.null(hsic_res)) NA_real_ else as.numeric(hsic_res$p.value)
  hsic_verdict <- if (is.na(hsic_p)) "?"
  else if (hsic_p >= 0.05) "A" else "C"
  hsic_agrees  <- if (cr_verdict == hsic_verdict) "Yes" else "No"

  cat(sprintf(
    "  z_het    r=%.3f  cor(z,u)=%+.4f  cor(x,u)=%+.4f  cor(z,x)=%+.4f  beta_IV=%+.3f  CI_Bei=[%+.3f,%+.3f]  p_hat_B=%.3f (%d/%d ok)  -> %s (single %s)   HSIC=%.3g (%s)\n",
    sim$r, sim$realized_cor_zu, sim$realized_cor_xu, sim$realized_cor_zx,
    beta_iv, ci_lo, ci_hi, p_hat_B, n_boot_ok, B_RULE,
    cr_verdict, dial_single, hsic_p, hsic_verdict
  ))

  data.frame(
    Instrument    = "z_het",
    cor_zu_target = 0,
    cor_zu_real   = round(sim$realized_cor_zu, 3),
    DGP_sign      = "pos",
    lambda        = NA_real_,
    D_lo          = rxu[1], D_hi = rxu[2],
    beta_IV       = round(beta_iv, 3),
    ci_lo         = round(ci_lo, 3),
    ci_hi         = round(ci_hi, 3),
    MCUB          = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi),
    p_hat_B       = round(p_hat_B, 3),
    DIAL          = cr_verdict,
    DIAL_single   = dial_single,
    HSIC_p        = round(hsic_p, 3),
    HSIC_verdict  = hsic_verdict,
    HSIC_agrees   = hsic_agrees,
    n_boot_ok     = n_boot_ok,
    boot_unit     = boot_unit,
    stringsAsFactors = FALSE
  )
}
# <<< HET-ROW

# ==============================================================================
# LATEX ROWS FOR tab:synthetic
# ------------------------------------------------------------------------------
# Same columns and format as the table in the paper:
#   Instrument & rho_ZU & D & beta_IV & C_MCUB & p_hat_B & HSIC p & DIAL &
#   HSIC agrees?
# Paste the rows between the header \midrule and \bottomrule; the z_het row
# comes after its own \midrule. MCUB endpoints outside [-1, 1] are clipped
# to +/-1.000 and marked with \dagger, and the printed note line lists the
# unclipped values. The verdicts use the unclipped endpoints.
# ==============================================================================

print_tab_synthetic <- function(res) {
  snap0 <- function(v) ifelse(abs(v) < 1e-9, 0, v)
  tex_name <- function(s) {
    if (s == "z_het") "$Z_{\\text{het}}$"
    else sprintf("$Z_{%s}$", sub("^z_", "", s))
  }
  tex_rho <- function(v) {
    if (is.na(v)) "n/a"
    else if (snap0(v) == 0) "$\\phantom{+}0.00$"
    else sprintf("$%+.2f$", v)
  }
  tex_num <- function(v) if (is.na(v)) "n/a" else sprintf("$%.3f$", snap0(v))
  tex_set <- function(lo, hi) {
    if (is.na(lo) || is.na(hi)) return("n/a")
    mark <- if (lo < -1 || hi > 1) "^{\\dagger}" else ""
    sprintf("$[%+.3f,\\,%+.3f]%s$", snap0(max(lo, -1)), snap0(min(hi, 1)), mark)
  }
  tex_hsic <- function(p) {
    if (is.na(p)) "n/a"
    else if (p < 0.001) "$<0.001$"
    else sprintf("$%.3f$", p)
  }

  clipped <- character(0)
  for (i in seq_len(nrow(res))) {
    r  <- res[i, ]
    nm <- tex_name(r$Instrument)
    if (r$Instrument == "z_het") cat("\\midrule\n")
    cat(sprintf("%s & %s & $(%g, %g)$ & %s & %s & %s & %s & %s & %s \\\\\n",
                nm, tex_rho(r$cor_zu_target),
                snap0(r$D_lo), snap0(r$D_hi),
                tex_num(r$beta_IV), tex_set(r$ci_lo, r$ci_hi),
                tex_num(r$p_hat_B), tex_hsic(r$HSIC_p),
                r$DIAL, r$HSIC_agrees))
    if (!is.na(r$ci_lo) && r$ci_lo < -1)
      clipped <- c(clipped, sprintf("$%+.3f$ (%s)", r$ci_lo, nm))
    if (!is.na(r$ci_hi) && r$ci_hi > 1)
      clipped <- c(clipped, sprintf("$%+.3f$ (%s)", r$ci_hi, nm))
  }

  if (length(clipped)) {
    k   <- length(clipped)
    lst <- if (k == 1L) clipped
    else if (k == 2L) paste(clipped, collapse = " and ")
    else paste0(paste(clipped[-k], collapse = ", "), ", and ", clipped[k])
    cat("\nTable note (replaces the \\item for clipped endpoints):\n")
    cat(sprintf(paste0("\\item[$\\dagger$] MCUB endpoints of %s are clipped ",
                       "to $[-1, 1]$; the verdicts are unaffected.\n"), lst))
  }

  # Unlike the grid instruments, Z_het does not have an exact sample
  # correlation with U, so its beta_IV moves with the sample correlation.
  het <- res[res$Instrument == "z_het", , drop = FALSE]
  if (nrow(het) == 1L && !is.na(het$cor_zu_real) && abs(het$cor_zu_real) >= 0.0005)
    cat(sprintf(paste0("\\item $Z_{\\text{het}}$ has population $\\rho_{ZU} = 0$; ",
                       "in this sample, its correlation with $U$ is $%.3f$, ",
                       "which accounts for $\\hat\\beta_{\\text{IV}} = %.3f$.\n"),
                het$cor_zu_real, het$beta_IV))
  invisible(clipped)
}

# ==============================================================================
# RUN ALL
# ==============================================================================

t0 <- Sys.time()
res_list <- vector("list", length(JP_INSTRUMENT_GRID))
for (j in seq_along(JP_INSTRUMENT_GRID)) {
  res_list[[j]] <- run_one_instrument(j)
}
res_list[[length(res_list) + 1L]] <- run_het_instrument()
res <- do.call(rbind, res_list)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ==============================================================================
# REPORT
# ==============================================================================

cat("\n================================================================================\n")
cat("tab:synthetic results\n")
cat("================================================================================\n\n")
print(res, row.names = FALSE)
cat(sprintf("\nElapsed: %.1f s\n", elapsed))

if (any(!is.na(res$n_boot_ok) & res$n_boot_ok < B_RULE))
  cat(sprintf("\nNOTE: some MCUB sets failed in the bootstrap rule (n_boot_ok < %d);\n",
              B_RULE),
      "p_hat_B is the share among the resamples that succeeded.\n")
if (any(!is.na(res$boot_unit) & res$boot_unit != "raw"))
  cat("\nWARNING: some rows resampled the prepared rows, not the raw data.\n")

cat("\nLaTeX rows for tab:synthetic (paste between the header \\midrule and \\bottomrule):\n\n")
print_tab_synthetic(res)

utils::write.csv(res, "table1_results.csv", row.names = FALSE)
cat("\nWrote table1_results.csv\n")
cat("DONE.\n")
