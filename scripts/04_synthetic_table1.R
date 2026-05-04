# ==============================================================================
# DIAL_NeurIPS_table1_MCUB_sets.R   --   Table 1 (Synthetic Benchmark)
# ------------------------------------------------------------------------------
# Reproduces Table 1 under the Jones-Pewsey DGP with the per-row
# correctly-elicited domain protocol (Corollary 4):
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
# Avoids the dual-sweep AND-rule failure: by Corollary 4(c), testing
# under the wrong sign collapses coverage even for valid instruments,
# so two sweeps cannot both pass for any nontrivial cor(z, u).
#
# DGP (sourced from DIAL_NeurIPS_JP_dgp.R)
#   x ~ J-P (sinh-arcsinh transform of normal, kappa=0.9 + Skew-Normal(1.2))
#   u = lambda (x - mean(x)) + epsilon + v_perp_x
#   y = beta * x + u, beta = 2
#   z_j: 11 instruments at known cor(z, u) in {-0.20, -0.16, ..., 0.20},
#         cor(z, x) ~ 0.25
#
# VERDICT
#   Tolerance-based zero-in-MCUB on CI_Bei.
#   Two-sided tolerance: zero is "in MCUB" iff
#     (ci_lo <= ZERO_TOL) AND (ci_hi >= -ZERO_TOL)
#   so the test handles positive-D and negative-D MCUBs symmetrically.
#   ZERO_TOL = 0.05.
#
# COMPARISON BASELINE
#   HSIC test of independence between z and the IV residual e_hat = y - beta_IV * x.
#
# DEPENDENCIES (mandatory source order)
#   DIAL_SKIP_AUTORUN <- TRUE
#   source("Fin_sim3_clean.R")
#   source("Fin_Empirical5_clean.R")
#   source("DIAL_NeurIPS_JP_dgp.R")
#   then this script.
#
# OUTPUT
#   table1_results.csv  -- one row per instrument
#   Console: paste-ready Table 1.
#
# RUN
#   Rscript DIAL_NeurIPS_table1_MCUB_sets.R
# ==============================================================================

rm(list = ls())

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org/")
}
install_if_missing(c("dHSIC", "ivreg"))
suppressPackageStartupMessages({
  library(dHSIC); library(ivreg)
})

# ==============================================================================
# CONFIG
# ==============================================================================

N             <- 20000L
TARGET_COR_XU <- 0.30          # |rho_xu| target; sign set per-DGP
COR_ZX        <- 0.25
ALPHA         <- 0.10
ZERO_TOL      <- 0.05
RXU_POS       <- c(0, 0.80)    # used when lambda > 0 (cor_zu >= 0)
RXU_NEG       <- c(-0.80, 0)   # used when lambda < 0 (cor_zu  < 0)
MASTER_SEED   <- 1234L

# ==============================================================================
# FOUNDATION SOURCE
# ==============================================================================

FOUNDATION_FILES <- normalizePath(c("Fin_sim3_clean.R",
                                    "Fin_Empirical5_clean.R",
                                    "DIAL_NeurIPS_JP_dgp.R"))

assign("DIAL_SKIP_AUTORUN", TRUE, envir = .GlobalEnv)
for (.f in FOUNDATION_FILES) {
  if (!file.exists(.f)) stop("Foundation file not found: ", .f)
  invisible(capture.output(source(.f, local = FALSE)))
}

stopifnot(exists("check_compatibility_simple"),
          exists("draw_jp_sample"),
          exists("calibrate_lambda"),
          exists("JP_INSTRUMENT_GRID"))

cat("================================================================================\n")
cat("DIAL Table 1: J-P synthetic benchmark (per-row correctly-elicited domain)\n")
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
      DIAL          = "?",
      HSIC_p        = NA_real_,
      HSIC_verdict  = "?",
      HSIC_agrees   = "?",
      stringsAsFactors = FALSE
    ))
  }
  
  ci_bei <- as.character(cr$CI_Bei[1])
  m <- regmatches(ci_bei,
                  regexec("\\[\\s*(-?[0-9.]+)\\s*,\\s*(-?[0-9.]+)\\s*\\]",
                          ci_bei))[[1]]
  ci_lo <- if (length(m) >= 3L) as.numeric(m[2L]) else NA_real_
  ci_hi <- if (length(m) >= 3L) as.numeric(m[3L]) else NA_real_
  
  # Two-sided tolerance: handles positive and negative D symmetrically.
  zero_in_tol <- !is.na(ci_lo) && !is.na(ci_hi) &&
    (ci_lo <= ZERO_TOL) && (ci_hi >= -ZERO_TOL)
  cr_verdict  <- if (is.na(ci_lo) || is.na(ci_hi)) "?"
  else if (zero_in_tol) "A" else "C"
  
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
    "  z_%+0.2f  sign=%s  D=[%+.2f,%+.2f]  beta_IV=%+.3f  CI_Bei=[%+.3f,%+.3f]  -> %s   HSIC=%.3f (%s)\n",
    cor_target, sign, rxu[1], rxu[2], beta_iv, ci_lo, ci_hi,
    cr_verdict, hsic_p, hsic_verdict
  ))
  
  data.frame(
    Instrument    = sprintf("z_%+0.2f", cor_target),
    cor_zu_target = cor_target,
    cor_zu_real   = round(sim$realized_cor_zu[j], 3),
    DGP_sign      = sign,
    lambda        = round(lambda, 3),
    D_lo          = rxu[1], D_hi = rxu[2],
    beta_IV       = round(beta_iv, 3),
    ci_lo         = round(ci_lo, 3),
    ci_hi         = round(ci_hi, 3),
    MCUB          = sprintf("[%+.3f, %+.3f]", ci_lo, ci_hi),
    DIAL          = cr_verdict,
    HSIC_p        = round(hsic_p, 3),
    HSIC_verdict  = hsic_verdict,
    HSIC_agrees   = hsic_agrees,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# RUN ALL
# ==============================================================================

t0 <- Sys.time()
res_list <- vector("list", length(JP_INSTRUMENT_GRID))
for (j in seq_along(JP_INSTRUMENT_GRID)) {
  res_list[[j]] <- run_one_instrument(j)
}
res <- do.call(rbind, res_list)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ==============================================================================
# REPORT
# ==============================================================================

cat("\n================================================================================\n")
cat("TABLE 1 (paste-ready)\n")
cat("================================================================================\n\n")
print(res, row.names = FALSE)
cat(sprintf("\nElapsed: %.1f s\n", elapsed))

cat("\nLaTeX rows for tab:synthetic:\n\n")
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  hsic_str <- if (is.na(r$HSIC_p)) "n/a"
  else if (r$HSIC_p < 0.001) "$<0.001$"
  else sprintf("%.3f", r$HSIC_p)
  cat(sprintf("$%s$ & %+.2f & %+.3f & %+.3f & [%+.2f, %+.2f] & %s & %s & %s & %s \\\\\n",
              gsub("_", "\\\\_", r$Instrument),
              r$cor_zu_target, r$cor_zu_real, r$beta_IV,
              r$D_lo, r$D_hi, r$MCUB, r$DIAL,
              hsic_str, r$HSIC_agrees))
}

utils::write.csv(res, "table1_results.csv", row.names = FALSE)
cat("\nWrote table1_results.csv\n")
cat("DONE.\n")