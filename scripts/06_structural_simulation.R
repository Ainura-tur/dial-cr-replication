# ==============================================================================
# 06_structural_simulation.R  -  Structural simulation for tab:sim
# ==============================================================================
#
# PURPOSE
#   Reproduces the five rows of tab:sim with the DGP of app:dgp. The rows
#   cover every DIAL outcome, and two of them reach Scenario B by
#   different routes:
#     A          valid instrument, homogeneous effect, confounder observed
#     B (het.)   valid instrument, effects differ by latent type, and the
#                low-effect type complies more, so LATE < ATE
#     B (conf.)  the sample of row A, but the analyst does not observe the
#                confounder C, so DML is biased
#     C          exclusion violation (gamma = 0.30)
#     C (weak)   weak instrument (first-stage F about 7)
#
# DESIGN NOTES
#   - Every row has its own IV and DML estimates.
#   - Ratio rule: symmetric, B when |log R| >= log(tau_R); IV and DML
#     estimates of opposite sign give "B (sign)" (the sign gate).
#   - Each row reports the implied endogeneity r*, the zero of the CR map
#     g(rho, r). On the prepared data, r* = cor(x, y - b_IV x) and
#     b_OLS - b_IV = r* sd(y - b_IV x) / sd(x). r* therefore absorbs both
#     residual confounding and effect heterogeneity; the two B rows have
#     similar r* for different reasons.
#   - Every DML fit starts from set.seed(dml_seed), so its estimate does
#     not depend on what ran before it (bootstrap rule in parallel,
#     sequential or cached).
#   - Everything the script needs is in Fin_Empirical5_clean.R and
#     bootstrap_rule.R.
#   - IV and the first-stage F come from lm() via Frisch-Waugh. They equal
#     the 2SLS estimate and the weak-instrument F of ivreg().
#
# DGP (app:dgp)
#     Y = beta_i D + gamma Z + phi' H + U
#     D = pi_i Z + lambda' H + V
#     V = 0.40 C + eps_V,  U = 0.40 C + eps_U,  eps ~ N(0, 0.6^2)
#     H ~ N(0, I_5), C ~ N(0, 1), Z ~ N(0, 1); phi, lambda ~ N(0, 1),
#     drawn once under COEF_SEED = 2026.
#   Latent types High and Low have equal probability.
#     Row A, B (conf.), C, C (weak): beta = 0.425 for both types.
#     Row B (het.): beta_H = 0.80, beta_L = 0.05 (ATE = 0.425) and
#       pi_H = 0.05, pi_L = 0.60. IV targets
#       LATE = E[beta pi] / E[pi] = 0.108. DML (partially linear model)
#       targets E[beta w] / E[w] with w = pi^2 + 0.6^2, which is 0.301.
#   Controls seen by the analyst: (H, C) in every row except B (conf.),
#   where only H is seen. Then DML targets
#   0.425 + Cov(V, U) / Var(D | H) = 0.425 + 0.16 / 0.5825 = 0.700.
#
# CR TEST
#   Domain D = (0, 0.8); alpha = 0.10; ZERO_TOL = 0.05 (footnote fn:tol).
#   The CR test residualises on the same controls as DML.
#
# BOOTSTRAP VALIDITY RULE (prop:finite_sample)
#   foundation/bootstrap_rule.R wraps check_compatibility_simple(): B = 500
#   resamples of the raw rows, one MCUB set per resample; the instrument
#   passes when p_hat_B >= TAU_V = 0.95 (column Valid_tau2, computed from
#   the unrounded share). Results are cached in results/boot_rule/, so a
#   rerun with the same data and seeds reuses them.
#
# OUTPUT
#   - console table and LaTeX rows for tab:sim
#   - results/sim_dial_v7_results.rds
#   - results/coef_vectors.rds
#
# RUNTIME  about 40 minutes on 7 cores when the rule results of rows C and
#   C (weak) are already cached; longer from a fresh clone, because the
#   rule then runs for all five rows.
#
# RUN
#   Rscript scripts/06_structural_simulation.R
# ==============================================================================

rm(list = ls())
gc()

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org/")
}
install_if_missing(c("dplyr", "MASS", "DoubleML", "mlr3", "mlr3learners",
                     "ranger", "data.table", "future.apply"))

suppressPackageStartupMessages({
  library(dplyr); library(MASS); library(DoubleML); library(mlr3)
  library(mlr3learners); library(ranger); library(data.table)
})

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;    filter    <- dplyr::filter
mutate    <- dplyr::mutate;    slice     <- dplyr::slice
recode    <- dplyr::recode;    rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange   <- dplyr::arrange;   count     <- dplyr::count
lag       <- dplyr::lag

# ==============================================================================
# REPO-ROOT DETECTION  (works under Rscript and in interactive sessions)
# ==============================================================================

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) && nchar(f)) dirname(normalizePath(f)) else getwd()
  }
)
if (!length(.script_dir) || is.na(.script_dir[1])) .script_dir <- getwd()
.repo_root <- normalizePath(file.path(.script_dir, ".."))

FOUNDATION_DIR <- file.path(.repo_root, "foundation")
RESULTS_DIR    <- file.path(.repo_root, "results")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ==============================================================================
# FOUNDATION: the CR test (Fin_Empirical5_clean.R) and the bootstrap rule
# ==============================================================================

assign("DIAL_SKIP_AUTORUN", TRUE, envir = .GlobalEnv)
.f <- file.path(FOUNDATION_DIR, "Fin_Empirical5_clean.R")
if (!file.exists(.f)) stop("Foundation file not found: ", .f,
                           "\n  Expected layout: <repo_root>/foundation/<name>",
                           "\n  .repo_root currently = ", .repo_root)
invisible(capture.output(source(.f, local = FALSE)))
rm(.f)

.missing <- setdiff(c("check_compatibility_simple", "CIhybrid",
                      "estimate_cov_corr_boot"), ls(envir = .GlobalEnv))
if (length(.missing))
  stop("Foundation source did not provide: ", paste(.missing, collapse = ", "))

B_RULE   <- 500L            # bootstrap resamples for the validity rule
TAU_V    <- 0.95            # validity threshold (tau_V in the paper)
TAU_R    <- 1.5             # ratio threshold tau_R
ZERO_TOL <- 0.05            # equivalence margin (footnote fn:tol)
ALPHA    <- 0.10            # MCUB level in the simulations
RXU      <- c(0, 0.8)       # maintained domain D for rho_DU
options(dial.boot_B = B_RULE, dial.tau2 = TAU_V)
source(file.path(FOUNDATION_DIR, "bootstrap_rule.R"))
dial_boot_plan()
options(dial.zero_tol = ZERO_TOL)   # same margin in the bootstrap rule

# ==============================================================================
# COEFFICIENT VECTORS phi, lambda  (drawn once, fixed across rows)
# ==============================================================================

COEF_SEED <- 2026L
H_DIM     <- 5L

cat("================================================================================\n")
cat(sprintf("DIAL Structural Simulation  (tab:sim; DGP from app:dgp; COEF_SEED = %d)\n",
            COEF_SEED))
cat("================================================================================\n")
cat(sprintf("  repo_root = %s\n\n", .repo_root))

set.seed(COEF_SEED)
phi    <- rnorm(H_DIM)
lambda <- rnorm(H_DIM)
saveRDS(list(phi = phi, lambda = lambda, seed = COEF_SEED, H_dim = H_DIM),
        file.path(RESULTS_DIR, "coef_vectors.rds"))
cat("phi    =", sprintf("%+.4f", phi),    "\n")
cat("lambda =", sprintf("%+.4f", lambda), "\n\n")

# ==============================================================================
# DATA-GENERATING PROCESS
# ==============================================================================

SIGMA_EPS    <- 0.6
CONF_LOADING <- 0.40   # cor(V, U) = 0.16 / 0.52 = 0.31 before conditioning on C

generate_data <- function(n, beta_h, beta_l, pi_h, pi_l, gamma_z, seed,
                          phi_vec = phi, lambda_vec = lambda,
                          sigma_eps = SIGMA_EPS) {
  set.seed(seed)
  H <- matrix(rnorm(n * H_DIM), nrow = n, ncol = H_DIM)
  colnames(H) <- paste0("H", seq_len(H_DIM))
  C <- rnorm(n)
  Z <- rnorm(n)
  type   <- sample(c("High", "Low"), size = n, replace = TRUE)
  beta_i <- ifelse(type == "High", beta_h, beta_l)
  pi_i   <- ifelse(type == "High", pi_h,   pi_l)
  eps_V  <- rnorm(n, sd = sigma_eps)
  eps_U  <- rnorm(n, sd = sigma_eps)
  V <- CONF_LOADING * C + eps_V
  U <- CONF_LOADING * C + eps_U
  D <- pi_i * Z + as.numeric(H %*% lambda_vec) + V
  Y <- beta_i * D + gamma_z * Z + as.numeric(H %*% phi_vec) + U
  out <- as.data.frame(H)
  out$C <- C; out$Z <- Z; out$D <- D; out$Y <- Y
  out$type <- type; out$beta <- beta_i; out$pi <- pi_i
  out
}

# Population targets of IV (LATE) and of DML (partially linear model) for a
# row, given the type parameters and whether the analyst sees C.
row_targets <- function(beta_h, beta_l, pi_h, pi_l, observe_C) {
  b  <- c(beta_h, beta_l); p <- c(pi_h, pi_l)
  late <- sum(b * p) / sum(p)
  if (observe_C) {
    w   <- p^2 + SIGMA_EPS^2                   # Var(D | H, C, type)
    dml <- sum(b * w) / sum(w)
  } else {                                     # homogeneous rows only
    vdh <- mean(p^2) + CONF_LOADING^2 + SIGMA_EPS^2
    dml <- mean(b) + CONF_LOADING^2 / vdh
  }
  c(LATE = late, ATE = mean(b), DML = dml)
}

# ==============================================================================
# ESTIMATORS
# ==============================================================================

# Implied endogeneity r*: the zero of the CR map
#   g(rho, r) = rho_xz r - (rho_xy rho_xz - rho_yz) sqrt((1 - r^2) / (1 - rho_xy^2)),
# rho = (rho_xy, rho_xz, rho_yz) of the prepared data. Equals
# cor(x, y - b_IV x).
implied_rdu <- function(rho) {
  k <- (rho[1] * rho[2] - rho[3]) / (rho[2] * sqrt(1 - rho[1]^2))
  k / sqrt(1 + k^2)
}

# DML: partially linear model, random forests, 5-fold cross-fitting.
# (Defined only if absent, so a test harness can supply its own.)
if (!exists("dml_plr", mode = "function")) {
  dml_plr <- function(df, Y, D, controls, num.trees = 300L, n_folds = 5L) {
    if (requireNamespace("lgr", quietly = TRUE))
      lgr::get_logger("mlr3")$set_threshold("warn")
    dd <- DoubleML::DoubleMLData$new(
      data = data.table::as.data.table(df[, c(Y, D, controls)]),
      y_col = Y, d_cols = D, x_cols = controls)
    dm <- DoubleML::DoubleMLPLR$new(
      dd,
      ml_l = mlr3::lrn("regr.ranger", num.trees = num.trees, min.node.size = 5),
      ml_m = mlr3::lrn("regr.ranger", num.trees = num.trees, min.node.size = 5),
      n_folds = n_folds)
    dm$fit()
    as.numeric(dm$coef)
  }
}

run_dial <- function(df, label, seed, controls, targets) {
  n        <- nrow(df)
  ctrl_str <- paste(controls, collapse = " + ")
  fml      <- function(lhs) stats::as.formula(paste(lhs, "~", ctrl_str))

  # --- DML, from a fixed RNG state ---
  set.seed(seed + 7L)
  beta_dml <- tryCatch(dml_plr(df, "Y", "D", controls), error = function(e) {
    cat(sprintf("    DML FAILED for %s: %s\n", label, conditionMessage(e)))
    NA_real_
  })

  # --- IV, OLS and first-stage F (Frisch-Waugh) ---
  y_r <- stats::resid(stats::lm(fml("Y"), data = df))
  x_r <- stats::resid(stats::lm(fml("D"), data = df))
  z0  <- stats::resid(stats::lm(fml("Z"), data = df))
  beta_iv  <- sum(z0 * y_r) / sum(z0 * x_r)
  beta_ols <- sum(x_r * y_r) / sum(x_r^2)
  fs_full  <- stats::lm(stats::as.formula(paste("D ~ Z +", ctrl_str)), data = df)
  first_F  <- stats::anova(stats::lm(fml("D"), data = df), fs_full)$F[2]

  # --- CR test on the prepared frame, with the bootstrap rule ---
  z_p   <- stats::predict(stats::lm(z0 ~ x_r + y_r))
  df_cr <- dial_attach_raw(data.frame(x = x_r, y = y_r, z = z_p),
                           df, "Y", "D", "Z", controls)
  rho    <- c(stats::cor(x_r, y_r), stats::cor(x_r, z_p), stats::cor(y_r, z_p))
  r_star <- implied_rdu(rho)
  u_hat  <- y_r - beta_iv * x_r

  ident_set <- NA_character_; p_hat_B <- NA_real_; n_boot_ok <- NA_integer_
  valid_tau2 <- NA_character_
  tryCatch({
    cr <- check_compatibility_simple(df_cr, i = as.integer(seed), alpha = ALPHA,
                                     rxu_range = RXU)
    if (!is.null(cr) && nrow(cr) > 0L) {
      ident_set  <- as.character(cr$CI_Bei[1])
      if (!is.null(cr$p_hat_B))    p_hat_B    <- as.numeric(cr$p_hat_B[1])
      if (!is.null(cr$Valid_tau2)) valid_tau2 <- as.character(cr$Valid_tau2[1])
      if (!is.null(cr$n_boot_ok))  n_boot_ok  <- as.integer(cr$n_boot_ok[1])
    }
  }, error = function(e) cat(sprintf("    CR ERROR: %s\n", conditionMessage(e))))

  ci <- suppressWarnings(as.numeric(regmatches(
    ident_set, gregexpr("-?[0-9]+\\.?[0-9]*", ident_set))[[1]]))
  ci_lo <- if (length(ci) >= 2) ci[1] else NA_real_
  ci_hi <- if (length(ci) >= 2) ci[2] else NA_real_
  zero_in_tol <- !is.na(ci_lo) && ci_lo <= ZERO_TOL && ci_hi >= -ZERO_TOL

  # --- DIAL classification (Algorithm 1) ---
  valid <- if (!is.na(valid_tau2)) valid_tau2 == "Valid" else zero_in_tol
  ratio <- abs(beta_iv) / abs(beta_dml)
  classify <- function(valid) {
    if (is.na(first_F) || first_F < 10)        return("C (weak)")
    if (is.na(valid) || is.na(ci_lo))          return("?")
    if (!valid)                                return("C (invalid)")
    if (is.na(beta_dml))                       return("?")
    if (sign(beta_iv) != sign(beta_dml))       return("B (sign)")
    if (abs(log(ratio)) >= log(TAU_R))         return("B")
    "A"
  }
  dial        <- classify(valid)
  dial_single <- classify(zero_in_tol)

  cat(sprintf(
    "  %-26s n=%5d DML=%6.3f IV=%6.3f R=%5.2f F=%7.0f %17s p_hat_B=%5.3f (%s ok) r*=%+.3f -> %s (single %s)\n",
    label, n, beta_dml, beta_iv, ratio, first_F, ident_set, p_hat_B,
    as.character(n_boot_ok), r_star, dial, dial_single))
  cat(sprintf(
    "  %-26s targets: LATE=%.3f ATE=%.3f DML=%.3f | OLS-IV=%+.4f, r*sd(u)/sd(x)=%+.4f\n",
    "", targets["LATE"], targets["ATE"], targets["DML"],
    beta_ols - beta_iv, r_star * stats::sd(u_hat) / stats::sd(x_r)))

  list(label = label, n = n, controls = controls,
       beta_dml = beta_dml, beta_iv = beta_iv, beta_ols = beta_ols,
       ratio = ratio, first_F = first_F, ident_set = ident_set,
       ci_lo = ci_lo, ci_hi = ci_hi, p_hat_B = p_hat_B, n_boot_ok = n_boot_ok,
       r_star = r_star, dial = dial, dial_single = dial_single,
       targets = targets)
}

# ==============================================================================
# PARAMETERS AND DATA
# ==============================================================================

n             <- 50000L
beta_hom      <- 0.425                 # = ATE of the two-type design
beta_high     <- 0.80;  beta_low <- 0.05
pi_hom        <- 0.25
pi_high       <- 0.05;  pi_low   <- 0.60   # low-effect type complies more
gamma_invalid <- 0.30
pi_weak       <- sqrt(5 * (CONF_LOADING^2 + SIGMA_EPS^2) / n)   # F about 7

H_ONLY <- paste0("H", seq_len(H_DIM))
H_AND_C <- c(H_ONLY, "C")

cat(sprintf("Parameters: n = %d, beta = %.3f (homogeneous rows); beta_H = %.2f, beta_L = %.2f\n",
            n, beta_hom, beta_high, beta_low))
cat(sprintf("            pi = %.2f; pi_H = %.2f, pi_L = %.2f; gamma = %.2f; pi_weak = %.4f\n\n",
            pi_hom, pi_high, pi_low, gamma_invalid, pi_weak))

df_hom   <- generate_data(n, beta_hom, beta_hom, pi_hom, pi_hom, 0, seed = 101L)
df_het   <- generate_data(n, beta_high, beta_low, pi_high, pi_low, 0, seed = 102L)
df_inv   <- generate_data(n, beta_hom, beta_hom, pi_hom, pi_hom,
                          gamma_invalid, seed = 201L)
df_weak  <- generate_data(n, beta_hom, beta_hom, pi_weak, pi_weak, 0, seed = 401L)

# ==============================================================================
# tab:sim ROWS
# ==============================================================================

cat("================================================================================\n")
cat("tab:sim: structural simulation -- every DIAL outcome from one DGP\n")
cat("================================================================================\n\n")

res <- list(
  run_dial(df_hom,  "A: homogeneous",          111L, H_AND_C,
           row_targets(beta_hom, beta_hom, pi_hom, pi_hom, TRUE)),
  run_dial(df_het,  "B: heterogeneous effects", 112L, H_AND_C,
           row_targets(beta_high, beta_low, pi_high, pi_low, TRUE)),
  run_dial(df_hom,  "B: C unobserved",          113L, H_ONLY,
           row_targets(beta_hom, beta_hom, pi_hom, pi_hom, FALSE)),
  run_dial(df_inv,  "C: exclusion violation",   211L, H_AND_C,
           row_targets(beta_hom, beta_hom, pi_hom, pi_hom, TRUE)),
  run_dial(df_weak, "C: weak instrument",       411L, H_AND_C,
           row_targets(beta_hom, beta_hom, pi_weak, pi_weak, TRUE))
)

# ==============================================================================
# SUMMARY AND LaTeX ROWS
# ==============================================================================

cat("\n================================================================================\n")
cat("MASTER SUMMARY (tab:sim)\n")
cat("================================================================================\n\n")
cat(sprintf("%-26s %7s %7s %6s %8s %17s %7s %7s %12s\n",
            "Analysis", "DML", "IV", "R", "F", "Ident.set", "p_hat_B", "r*", "DIAL"))
cat(strrep("-", 104), "\n")
for (r in res)
  cat(sprintf("%-26s %7.3f %7.3f %6.2f %8.0f %17s %7.3f %+7.3f %12s\n",
              r$label, r$beta_dml, r$beta_iv, r$ratio, r$first_F,
              r$ident_set, r$p_hat_B, r$r_star, r$dial))

fmt_end <- function(v) sprintf("%+.3f", max(-1, min(1, v)))
cat("\nLaTeX rows for tab:sim (endpoints clipped to [-1, 1]):\n\n")
for (r in res)
  cat(sprintf("%s & $%.3f$ & $%.3f$ & $%.2f$ & $%s$ & $[%s,\\,%s]$ & $%.3f$ & $%+.3f$ & %s \\\\\n",
              r$label, r$beta_dml, r$beta_iv, r$ratio,
              formatC(round(r$first_F), format = "d", big.mark = "{,}"),
              fmt_end(r$ci_lo), fmt_end(r$ci_hi), r$p_hat_B, r$r_star, r$dial))

saveRDS(res, file.path(RESULTS_DIR, "sim_dial_v7_results.rds"))
cat(sprintf("\nCoef seed = %d (phi, lambda in results/coef_vectors.rds)\n", COEF_SEED))
cat("\nDONE.\n")
