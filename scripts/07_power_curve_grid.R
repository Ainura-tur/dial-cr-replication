# =============================================================================
# 07_power_curve_grid.R
#
# Lean power-curve MC (based on run5_LEAN + extend).
#
# Per-rep cost is dominated by a single CIhybrid call. Relative to the full
# version, three expensive operations are removed:
#   1. pvalue_mcub_zero_fast  -- 18 x CIhybrid at Blarge 3000 per rep.
#   2. ci_simple_union        -- 800-rep bootstrap always called first.
#   3. Duplicate estimate_cov_corr_boot at B = 800 inside ci_simple_union.
#
# Workers source only Fin_Empirical5_clean.R, avoiding the heavyweight
# top-level execution in Fin_sim3_clean.R (stargazer, fGarch, ivreg calls).
#
# Grid: 35 base delta points + 8 extension points (indices 36-43) that fill
# in the transition zone for n=1107 and n=50000. Both sets share the same
# cache directory (results/cp_grid_lean/) and can be run in any order.
#
# Usage:
#   Rscript scripts/07_power_curve_grid.R
#
# Runtime: ~75 minutes on 7 cores.
# =============================================================================

rm(list = ls())

# ---- Repo-root detection ----------------------------------------------------
.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) && nchar(f)) dirname(normalizePath(f))
    else getwd()
  }
)
.repo_root <- normalizePath(file.path(.script_dir, ".."))

EMP5_PATH <- file.path(.repo_root, "foundation", "Fin_Empirical5_clean.R")
stopifnot(file.exists(EMP5_PATH))

# ---- Dependencies -----------------------------------------------------------
install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org/")
}
install_if_missing(c("MASS", "future.apply",
                     "ivreg", "lmtest", "sandwich", "haven",
                     "dplyr", "Matrix", "numDeriv"))

suppressPackageStartupMessages({ library(MASS); library(future.apply) })

# ---- Source EMP5 on main process -------------------------------------------
DIAL_SKIP_AUTORUN <- TRUE
source(EMP5_PATH)

required_fns <- c("CIhybrid", "g_xu_safe", "local_compute_gradient_safe",
                  "estimate_cov_corr_boot")
missing_fns  <- required_fns[!vapply(required_fns, exists, logical(1))]
if (length(missing_fns))
  stop("Functions still missing after sourcing EMP5: ",
       paste(missing_fns, collapse = ", "))
cat("Fin_Empirical5_clean.R loaded OK.\n")

# ---- Parallel backend -------------------------------------------------------
N_WORKERS <- 7L
plan(multisession, workers = N_WORKERS)
cat(sprintf("Parallel plan: multisession with %d workers\n\n", N_WORKERS))

# =============================================================================
# Configuration
# =============================================================================

N_GRID <- c(1107L, 1534L, 5000L, 10000L, 50000L)

# 35-point base grid
DELTA_GRID <- local({
  low_anchors  <- c(0.001, 0.005, 0.010, 0.020, 0.030, 0.045, 0.055, 0.060)
  pre_trans    <- seq(0.07,  0.10,  length.out = 4)
  trans_dense  <- seq(0.105, 0.155, length.out = 11)
  post_trans   <- seq(0.165, 0.20,  length.out = 5)
  high_anchors <- c(0.22, 0.25, 0.28, 0.30, 0.32, 0.35, 0.40)
  sort(unique(c(low_anchors, pre_trans, trans_dense, post_trans, high_anchors)))
})
stopifnot(length(DELTA_GRID) == 35L)

# 8-point extension: fills in transition zone for n=1107 and n=50000
EXT_DELTAS  <- c(0.131, 0.132, 0.133, 0.134, 0.167, 0.169, 0.171, 0.173)
EXT_INDICES <- 36L:43L
stopifnot(length(EXT_DELTAS) == length(EXT_INDICES))

# R_MC schedule: more reps at small n where curves are plotted
R_MC_BY_N <- c(
  "1107"  = 250L,
  "1534"  = 250L,
  "5000"  = 150L,
  "10000" = 100L,
  "50000" =  50L
)

# CIhybrid tuning
B_MCUB <- 199L   # per-bisection draws
BLARGE <- 499L   # c_bd and final coverage draws
B_BOOT <- 200L   # bootstrap B for estimate_cov_corr_boot

OUT_DIR     <- file.path(.repo_root, "results", "cp_grid_lean")
SUMMARY_CSV <- file.path(.repo_root, "results", "cp_grid_summary_lean.csv")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

ALPHA     <- 0.05
TAU2      <- 0.95
RHO_DU    <- 0.50
COR_ZD    <- 0.25
BETA_TRUE <- 0.30
RXU_RANGE <- c(0.0, 0.80)

total_base <- length(N_GRID) * length(DELTA_GRID)
total_ext  <- length(N_GRID) * length(EXT_DELTAS)
cat(sprintf("N grid         : %s\n", paste(N_GRID, collapse = ", ")))
cat(sprintf("Base grid      : %d points in [%.4f, %.4f]\n",
            length(DELTA_GRID), min(DELTA_GRID), max(DELTA_GRID)))
cat(sprintf("Extension      : %d points (%s)\n",
            length(EXT_DELTAS), paste(EXT_DELTAS, collapse = ", ")))
cat(sprintf("Total cells    : %d base + %d ext = %d\n",
            total_base, total_ext, total_base + total_ext))
cat(sprintf("R_MC schedule  : %s\n",
            paste(sprintf("n=%s:%d", names(R_MC_BY_N), R_MC_BY_N),
                  collapse = " | ")))
cat(sprintf("B_mcub / Blarge / B_boot : %d / %d / %d\n",
            B_MCUB, BLARGE, B_BOOT))
cat(sprintf("Output dir     : %s\n\n", OUT_DIR))

# =============================================================================
# DGP helpers (defined locally; exported to workers via future.globals)
# =============================================================================

construct_instrument <- function(Dtilde, Utilde, rho, a, b, seed_w = 1L) {
  set.seed(seed_w)
  n     <- length(Dtilde)
  V     <- (Utilde - rho * Dtilde) / sqrt(1 - rho^2)
  W_raw <- rnorm(n)
  W     <- resid(lm(W_raw ~ Dtilde + V))
  W     <- W / sd(W)
  c1    <- a
  c2    <- (b - a * rho) / sqrt(1 - rho^2)
  c3_sq <- 1 - c1^2 - c2^2
  if (c3_sq <= 0) stop("Infeasible Cholesky weights.")
  c1 * Dtilde + c2 * V + sqrt(c3_sq) * W
}

simulate_one <- function(n, delta, rep_seed) {
  set.seed(rep_seed)
  Sigma <- matrix(c(1, RHO_DU, RHO_DU, 1), 2, 2)
  DU    <- MASS::mvrnorm(n, c(0, 0), Sigma)
  Z     <- construct_instrument(DU[, 1], DU[, 2],
                                rho    = RHO_DU,
                                a      = COR_ZD,
                                b      = delta,
                                seed_w = rep_seed + 1L)
  data.frame(x = DU[, 1],
             y = BETA_TRUE * DU[, 1] + DU[, 2],
             z = Z)
}

# =============================================================================
# Lean one_rep
#
# Steps:
#   1. simulate_one
#   2. sample correlations
#   3. estimate_cov_corr_boot at B_BOOT
#   4. build r_grid / g / Al / Au
#   5. single CIhybrid call at B_MCUB / BLARGE
#   6. extract zero_in
#
# NOT called: pvalue_mcub_zero_fast, ci_simple_union.
# =============================================================================

one_rep <- function(rep_seed, n_val, delta,
                    rxu_range, alpha, b_mcub, blarge, b_boot,
                    emp5_path) {

  empty <- function(msg = NA_character_)
    list(zero_in = NA, ci_lo = NA_real_, ci_hi = NA_real_,
         plug_lo = NA_real_, plug_hi = NA_real_, err = msg)

  # Source EMP5 once per worker session
  if (!exists("CIhybrid", envir = globalenv(), inherits = FALSE)) {
    err <- tryCatch({
      assign("DIAL_SKIP_AUTORUN", TRUE, envir = globalenv())
      source(emp5_path, local = FALSE)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(err))
      return(empty(paste("worker EMP5 source failed:", err)))
    if (!exists("CIhybrid", envir = globalenv(), inherits = FALSE))
      return(empty("CIhybrid still missing after sourcing EMP5"))
  }

  # 1. Simulate
  df <- tryCatch(simulate_one(n_val, delta, rep_seed), error = function(e) NULL)
  if (is.null(df)) return(empty("simulate_one failed"))

  # 2. Sample correlations
  rho_xy   <- cor(df$x, df$y)
  rho_xz   <- cor(df$x, df$z)
  rho_yz   <- cor(df$y, df$z)
  deltahat <- c(rho_xy, rho_xz, rho_yz)
  if (any(!is.finite(deltahat))) return(empty("non-finite correlations"))

  # 3. Bootstrap covariance
  deltaSigma <- tryCatch(
    estimate_cov_corr_boot(df$x, df$y, df$z, B = b_boot, seed = rep_seed + 2L),
    error = function(e) NULL
  )
  if (is.null(deltaSigma)) return(empty("estimate_cov_corr_boot failed"))

  # 4. Grid, g map, Jacobian rows
  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = 50L)
  g      <- function(delta_mat) g_xu_safe(r_grid, delta_mat)
  A      <- t(vapply(r_grid, function(r_xu)
    as.numeric(local_compute_gradient_safe(
      r_xu, deltahat[1], deltahat[2], deltahat[3])),
    numeric(3L)))
  Al <- A; Au <- A

  eta    <- 0.001
  alphac <- 0.8 * alpha
  tol    <- 1e-3
  tol_r  <- 1e-3

  # 5. Single CIhybrid call
  res <- tryCatch(
    CIhybrid(deltahat, deltaSigma, Al, Au,
             alpha  = alpha, alphac = alphac, eta = eta,
             B      = b_mcub, Blarge = blarge,
             tol    = tol, tol_r  = tol_r,
             index  = NULL, g = g, seed = rep_seed + 3L),
    error = function(e) list(.err = conditionMessage(e))
  )
  if (!is.null(res$.err)) return(empty(res$.err))

  ci   <- res$CI_h
  plug <- res$CI_c

  list(
    zero_in = !is.na(ci[1]) && ci[1] <= 0 && 0 <= ci[2],
    ci_lo   = ci[1],   ci_hi   = ci[2],
    plug_lo = plug[1], plug_hi = plug[2],
    err     = NA_character_
  )
}

WORKER_GLOBALS <- c("simulate_one", "construct_instrument",
                    "RHO_DU", "COR_ZD", "BETA_TRUE")

# =============================================================================
# Helper: run a set of (delta, d_idx) pairs for all n values
# =============================================================================

sc  <- function(x, fill = NA_real_)
  if (is.null(x) || length(x) == 0) fill else x[1]
scc <- function(x, fill = NA_character_)
  if (is.null(x) || length(x) == 0) fill else as.character(x[1])

run_grid <- function(deltas, indices, label = "") {
  cell_counter   <- 0L
  total_cells    <- length(N_GRID) * length(deltas)
  t_start        <- Sys.time()

  for (n_val in N_GRID) {
    R_MC <- R_MC_BY_N[[as.character(n_val)]]

    for (k in seq_along(deltas)) {
      delta        <- deltas[k]
      d_idx        <- indices[k]
      cell_counter <- cell_counter + 1L

      fname <- sprintf("n%05d_d%02d.rds", n_val, d_idx)
      fpath <- file.path(OUT_DIR, fname)

      if (file.exists(fpath)) {
        cat(sprintf("[%3d/%3d]%s n=%d d=%.4f  (cached)\n",
                    cell_counter, total_cells, label, n_val, delta))
        next
      }

      t_cell        <- Sys.time()
      seed_base     <- 20000L * n_val + d_idx
      rep_seeds     <- seed_base + seq_len(R_MC)

      rep_results <- future_lapply(
        rep_seeds, one_rep,
        n_val     = n_val,
        delta     = delta,
        rxu_range = RXU_RANGE,
        alpha     = ALPHA,
        b_mcub    = B_MCUB,
        blarge    = BLARGE,
        b_boot    = B_BOOT,
        emp5_path = normalizePath(EMP5_PATH),
        future.seed     = TRUE,
        future.packages = c("MASS", "ivreg", "lmtest", "sandwich",
                            "haven", "dplyr", "Matrix", "numDeriv"),
        future.globals  = WORKER_GLOBALS
      )

      zero_in <- vapply(rep_results, function(r) sc(r$zero_in, NA), logical(1))
      ci_lo   <- vapply(rep_results, function(r) sc(r$ci_lo),       numeric(1))
      ci_hi   <- vapply(rep_results, function(r) sc(r$ci_hi),       numeric(1))
      plug_lo <- vapply(rep_results, function(r) sc(r$plug_lo),     numeric(1))
      plug_hi <- vapply(rep_results, function(r) sc(r$plug_hi),     numeric(1))
      errs    <- vapply(rep_results, function(r) scc(r$err),        character(1))
      n_errs  <- sum(!is.na(errs))

      cell_result <- list(
        n            = n_val,      delta       = delta,
        delta_index  = d_idx,      cp_raw      = mean(zero_in,           na.rm = TRUE),
        plug_lo_mean = mean(plug_lo, na.rm = TRUE),
        plug_hi_mean = mean(plug_hi, na.rm = TRUE),
        plug_w_mean  = mean(plug_hi - plug_lo, na.rm = TRUE),
        ci_lo_mean   = mean(ci_lo,  na.rm = TRUE),
        ci_hi_mean   = mean(ci_hi,  na.rm = TRUE),
        ci_w_mean    = mean(ci_hi  - ci_lo,    na.rm = TRUE),
        n_valid_reps = sum(!is.na(zero_in)),
        n_errors     = n_errs,
        first_error  = if (n_errs > 0) scc(errs[!is.na(errs)][1]) else NA_character_,
        R_MC         = R_MC,       B_mcub      = B_MCUB,
        Blarge       = BLARGE,     B_boot      = B_BOOT,
        elapsed_sec  = as.numeric(difftime(Sys.time(), t_cell, units = "secs"))
      )
      saveRDS(cell_result, fpath)

      err_tag <- if (n_errs > 0) sprintf(" [%d errors]", n_errs) else ""
      cat(sprintf(
        "[%3d/%3d]%s n=%5d d=%.4f  cp=%.3f  plug_w=%.3f  CI_w=%.3f  (%.0fs, R=%d)%s\n",
        cell_counter, total_cells, label, n_val, delta,
        cell_result$cp_raw, cell_result$plug_w_mean, cell_result$ci_w_mean,
        cell_result$elapsed_sec, R_MC, err_tag))
      if (n_errs > 0)
        cat(sprintf("    first error: %s\n", cell_result$first_error))
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  cat(sprintf("\n%s grid done in %.1f minutes.\n\n", label, elapsed))
}

# =============================================================================
# Run both grids
# =============================================================================

cat("=== Base grid (35 points) ===\n")
run_grid(DELTA_GRID, seq_along(DELTA_GRID), label = " [base]")

cat("=== Extension grid (8 points) ===\n")
run_grid(EXT_DELTAS, EXT_INDICES, label = " [ext]")

# =============================================================================
# Aggregate all cached cells to summary CSV
# =============================================================================

ALL_INDICES <- 1L:43L
all_rows <- list(); skipped <- 0L

for (n_val in N_GRID) {
  for (d_idx in ALL_INDICES) {
    fpath <- file.path(OUT_DIR, sprintf("n%05d_d%02d.rds", n_val, d_idx))
    if (!file.exists(fpath)) next
    x <- tryCatch(readRDS(fpath), error = function(e) NULL)
    if (is.null(x)) { skipped <- skipped + 1L; next }
    row <- tryCatch(data.frame(
      n           = sc(x$n),
      delta       = sc(x$delta),
      delta_index = sc(x$delta_index),
      cp_raw      = sc(x$cp_raw),
      plug_w_mean = sc(x$plug_w_mean),
      ci_w_mean   = sc(x$ci_w_mean),
      ci_lo_mean  = sc(x$ci_lo_mean),
      ci_hi_mean  = sc(x$ci_hi_mean),
      n_valid     = sc(x$n_valid_reps, 0L),
      n_errors    = sc(x$n_errors,     0L)
    ), error = function(e) NULL)
    if (is.null(row)) { skipped <- skipped + 1L; next }
    all_rows[[length(all_rows) + 1L]] <- row
  }
}
if (skipped > 0) cat(sprintf("Skipped %d malformed rds files.\n", skipped))

if (length(all_rows) > 0) {
  summary_df <- do.call(rbind, all_rows)
  summary_df <- summary_df[order(summary_df$n, summary_df$delta), ]
  utils::write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)
  cat(sprintf("Wrote %d rows to %s\n\n", nrow(summary_df), SUMMARY_CSV))

  cat(strrep("=", 78), "\n")
  cat("cp_raw pivoted (rows = delta, columns = n) -- transition zone only\n")
  cat(strrep("=", 78), "\n")
  trans <- summary_df[summary_df$delta >= 0.125 & summary_df$delta <= 0.18, ]
  piv   <- reshape(trans[, c("n", "delta", "cp_raw")],
                   idvar = "delta", timevar = "n", direction = "wide")
  names(piv) <- gsub("cp_raw\\.", "n=", names(piv))
  piv[, -1]  <- round(piv[, -1], 3)
  print(piv, row.names = FALSE)
} else {
  cat("No complete cells found -- check errors above.\n")
}
