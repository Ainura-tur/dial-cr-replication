# ==============================================================================
# DIAL_NeurIPS_MonteCarlo.R   --   Table 6 (Monte Carlo Validation)  v2-LEAN
# ------------------------------------------------------------------------------
# Changes from v1:
#
#   [FAST-1]  Bypass check_compatibility_simple entirely.
#             That function calls (per rep):
#               estimate_cov_corr_boot  B = 800     (1 call)
#               ci_simple_union         B_boot = 800 (1 internal call)
#               CIhybrid                B=500  Blarge=5000
#               pvalue_mcub_zero_fast   B_fast=300 Blarge_fast=3000
#                                        x up to 18 CIhybrid calls
#             Lean one_rep calls:
#               estimate_cov_corr_boot  B = B_BOOT (200)   (1 call)
#               CIhybrid                B = B_MCUB Blarge = BLARGE
#             pvalue_mcub_zero_fast and ci_simple_union are NOT called.
#             Verdict uses CI_h endpoints only.
#
#   [FAST-2]  Workers source only Fin_Empirical5_clean.R and
#             DIAL_NeurIPS_JP_dgp.R; Fin_sim3_clean.R is NOT sourced.
#             Fin_sim3_clean.R has ~90 lines of top-level execution
#             (ivreg calls, heavy package loads) that run on every
#             worker source(), adding latency and fragile dependencies.
#
#   [FAST-3]  capture.output wrapper removed from the CR call.
#             check_compatibility_simple printed per-instrument progress
#             via cat(); those are gone since we call CIhybrid directly.
#
#   [FAST-4]  Stale-cache guard schema bumped to 8 so cells from the
#             old check_compatibility_simple schema are not reused
#             (they stored different fields).
#
# Expected speedup: 10-20x per rep (same argument as run5_LEAN).
#
# DEPENDENCIES (main session only -- workers load a subset):
#   DIAL_SKIP_AUTORUN <- TRUE
#   source("Fin_sim3_clean.R")        # only needed for main session
#   source("Fin_Empirical5_clean.R")
#   source("DIAL_NeurIPS_JP_dgp.R")
#   then this script.
#
# RUNTIME (estimated)
#   R = 50  sanity run,  7 cores  ~  5-20 minutes
#   R = 1000 full Table 6, 7 cores  ~  30-60 minutes  (was 4-8 hours)
#
# CACHE
#   table6_grid_lean/   (new dir; incompatible with v1 table6_grid/ cells)
#   If you want to keep old cells, point CACHE_DIR back to table6_grid/
#   and set SCHEMA_VERSION <- 7L so the stale-cache guard passes.
#
# OUTPUT
#   table6_results.csv
#   table6_grid_lean/*.rds
#   table6_timing.csv
# ==============================================================================

rm(list = ls())

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org/")
}
install_if_missing(c("future.apply", "ivreg", "lmtest", "sandwich",
                     "haven", "dplyr", "Matrix", "numDeriv"))
suppressPackageStartupMessages({ library(future.apply) })

# ==============================================================================
# CONFIG
# ==============================================================================

N             <- 5000L
R_PER_INSTR   <- 1000L        # set to 50L for a sanity check
TARGET_COR_XU <- 0.30
COR_ZX        <- 0.25
ZERO_TOL      <- 0.05         # verdict: 0 in CI iff lo <= ZERO_TOL & hi >= -ZERO_TOL
ALPHA         <- 0.10
RXU_POS       <- c( 0.00,  0.80)
RXU_NEG       <- c(-0.80,  0.00)
N_WORKERS     <- 7L
MASTER_SEED   <- 4242L

# [FAST-1] CIhybrid tuning (one call per rep)
B_MCUB  <- 199L   # per-bisection draw count
BLARGE  <- 499L   # c_bd and final coverage draws
B_BOOT  <- 200L   # estimate_cov_corr_boot bootstrap size

SCHEMA_VERSION <- 8L   # [FAST-4] bump so v1 cells are not reused

CACHE_DIR   <- file.path(normalizePath("."), "table6_grid_lean")
SUMMARY_CSV <- file.path(normalizePath("."), "table6_results.csv")
TIMING_LOG  <- file.path(normalizePath("."), "table6_timing.csv")
if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)

# [FAST-2] Workers do NOT source Fin_sim3_clean.R
WORKER_FILES <- normalizePath(c("Fin_Empirical5_clean.R",
                                "DIAL_NeurIPS_JP_dgp.R"))

# Main session still needs sim3 for calibrate_lambda / JP_INSTRUMENT_GRID
MAIN_FILES <- normalizePath(c("Fin_sim3_clean.R",
                              "Fin_Empirical5_clean.R",
                              "DIAL_NeurIPS_JP_dgp.R"))

cat("================================================================================\n")
cat("DIAL Table 6: Monte Carlo Validation (J-P DGP, per-row sign elicitation) v2-LEAN\n")
cat("================================================================================\n")
cat(sprintf("  n = %d,  R per instrument = %d\n", N, R_PER_INSTR))
cat(sprintf("  |target cor(x,u)| = %.2f,  cor(z,x) = %.2f\n",
            TARGET_COR_XU, COR_ZX))
cat(sprintf("  D pos = [%.2f, %.2f]   D neg = [%.2f, %.2f]\n",
            RXU_POS[1], RXU_POS[2], RXU_NEG[1], RXU_NEG[2]))
cat(sprintf("  alpha = %.2f,  ZERO_TOL = %.2f,  master seed = %d\n",
            ALPHA, ZERO_TOL, MASTER_SEED))
cat(sprintf("  B_mcub = %d,  Blarge = %d,  B_boot = %d\n",
            B_MCUB, BLARGE, B_BOOT))
cat(sprintf("  Workers = %d,  cache dir: %s\n\n", N_WORKERS, CACHE_DIR))

# ==============================================================================
# Load foundation in main session
# ==============================================================================

load_foundation <- function(files) {
  assign("DIAL_SKIP_AUTORUN", TRUE, envir = .GlobalEnv)
  for (f in files) {
    if (!file.exists(f)) stop("Foundation file not found: ", f)
    invisible(capture.output(source(f, local = FALSE)))
  }
  required <- c("CIhybrid", "g_xu_safe", "local_compute_gradient_safe",
                "estimate_cov_corr_boot",
                "draw_jp_sample", "calibrate_lambda", "JP_INSTRUMENT_GRID")
  miss <- setdiff(required, ls(envir = .GlobalEnv))
  if (length(miss))
    stop("Foundation source did not provide: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

cat("Sourcing foundation in main session... ")
load_foundation(MAIN_FILES)
cat("OK\n\n")

# ==============================================================================
# Calibrate lambda for both signs
# ==============================================================================

cat("Calibrating lambda...\n")
cal_pos <- calibrate_lambda(n = N, target_cor_xu =  TARGET_COR_XU,
                            cor_zx_default = COR_ZX,
                            seed = MASTER_SEED, sign = "positive", verbose = FALSE)
cal_neg <- calibrate_lambda(n = N, target_cor_xu =  TARGET_COR_XU,
                            cor_zx_default = COR_ZX,
                            seed = MASTER_SEED, sign = "negative", verbose = FALSE)
LAMBDA_POS <- cal_pos$lambda
LAMBDA_NEG <- cal_neg$lambda
cat(sprintf("  lambda_+ = %+.4f  (realized cor(x,u) = %+.4f)\n",
            LAMBDA_POS, cal_pos$realized))
cat(sprintf("  lambda_- = %+.4f  (realized cor(x,u) = %+.4f)\n\n",
            LAMBDA_NEG, cal_neg$realized))

# ==============================================================================
# Instrument grid + truth labels
# ==============================================================================

GRID <- JP_INSTRUMENT_GRID

truth_label <- function(cor_zu) {
  r <- abs(cor_zu)
  if (r <= 0.04 + 1e-9) return("A")
  if (r >= 0.12 - 1e-9) return("C")
  "Boundary"
}

INSTR_GRID <- lapply(seq_along(GRID), function(j) {
  list(idx    = j,
       name   = sprintf("z_%+0.2f", GRID[j]),
       cor_zu = GRID[j],
       truth  = truth_label(GRID[j]))
})

# ==============================================================================
# [FAST-1] Lean one_rep: CIhybrid only, no pvalue_mcub_zero_fast
# ==============================================================================

one_rep <- function(rep_seed, instr,
                    worker_files,
                    lambda_pos, lambda_neg,
                    n_sample, cor_zx_default,
                    alpha, rxu_pos, rxu_neg,
                    zero_tol, b_mcub, blarge, b_boot) {
  
  empty <- function(msg = NA_character_)
    list(ci_lo = NA_real_, ci_hi = NA_real_,
         verdict = NA_character_, zero_in = NA, err = msg)
  
  # [FAST-2] Source only EMP5 + JP_DGP on workers (no Fin_sim3_clean.R)
  if (!exists("CIhybrid", envir = globalenv(), inherits = FALSE)) {
    err <- tryCatch({
      assign("DIAL_SKIP_AUTORUN", TRUE, envir = globalenv())
      for (f in worker_files) source(f, local = FALSE)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(err)) return(empty(paste("worker setup failed:", err)))
    if (!exists("CIhybrid", envir = globalenv(), inherits = FALSE))
      return(empty("CIhybrid missing after worker setup"))
  }
  
  # Pick sign-appropriate lambda and domain
  if (instr$cor_zu < 0) {
    lam <- lambda_neg; rxu <- rxu_neg
  } else {
    lam <- lambda_pos; rxu <- rxu_pos
  }
  
  # 1. Draw J-P sample
  sim <- tryCatch(
    draw_jp_sample(n = n_sample, lambda = lam, seed = rep_seed,
                   cor_zu_grid    = c(instr$cor_zu),
                   cor_zx_default = cor_zx_default),
    error = function(e) NULL
  )
  if (is.null(sim)) return(empty("draw_jp_sample failed"))
  
  # 2. Residualise (identity residualisation: no controls)
  x_r <- sim$x - mean(sim$x)
  y_r <- sim$y - mean(sim$y)
  z0  <- sim$z[, 1L] - mean(sim$z[, 1L])
  z_p <- stats::predict(stats::lm(z0 ~ x_r + y_r))
  
  # 3. Sample correlations
  rho_xy   <- cor(x_r, y_r)
  rho_xz   <- cor(x_r, z_p)
  rho_yz   <- cor(y_r, z_p)
  deltahat <- c(rho_xy, rho_xz, rho_yz)
  if (any(!is.finite(deltahat))) return(empty("non-finite correlations"))
  
  # 4. Bootstrap covariance -- one call at reduced B [FAST-1]
  deltaSigma <- tryCatch(
    estimate_cov_corr_boot(x_r, y_r, z_p, B = b_boot, seed = rep_seed + 2L),
    error = function(e) NULL
  )
  if (is.null(deltaSigma)) return(empty("estimate_cov_corr_boot failed"))
  
  # 5. Grid, g, Jacobian
  r_grid <- seq(rxu[1], rxu[2], length.out = 50L)
  g      <- function(dm) g_xu_safe(r_grid, dm)
  A      <- t(vapply(r_grid, function(r_xu)
    as.numeric(local_compute_gradient_safe(
      r_xu, deltahat[1], deltahat[2], deltahat[3])),
    numeric(3L)))
  Al <- A; Au <- A
  
  # 6. Single CIhybrid call [FAST-1]
  eta    <- 0.001
  alphac <- 0.8 * alpha
  tol    <- 1e-3
  tol_r  <- 1e-3
  
  res <- tryCatch(
    CIhybrid(deltahat, deltaSigma, Al, Au,
             alpha  = alpha, alphac = alphac, eta = eta,
             B      = b_mcub, Blarge = blarge,
             tol    = tol, tol_r = tol_r,
             index  = NULL, g = g,
             seed   = rep_seed + 3L),
    error = function(e) list(.err = conditionMessage(e))
  )
  if (!is.null(res$.err)) return(empty(res$.err))
  
  lo <- res$CI_h[1]; hi <- res$CI_h[2]
  
  # Verdict: relaxed zero-coverage (ZERO_TOL neighbourhood)
  zin     <- (lo <= zero_tol) && (hi >= -zero_tol)
  verdict <- if (zin) "A" else "C"
  
  list(ci_lo = lo, ci_hi = hi, verdict = verdict,
       zero_in = zin, err = NA_character_)
}

# ==============================================================================
# Stale-cache guard
# ==============================================================================

.check_cache <- function() {
  fs <- list.files(CACHE_DIR, pattern = "\\.rds$", full.names = TRUE)
  if (!length(fs)) return(invisible(NULL))
  x <- tryCatch(readRDS(fs[1]), error = function(e) NULL)
  if (is.null(x)) return(invisible(NULL))
  if (is.null(x$schema_version) || x$schema_version < SCHEMA_VERSION) {
    cat("\n*** WARNING: stale cache detected (schema", x$schema_version,
        "< required", SCHEMA_VERSION, ")\n")
    cat("DELETE", CACHE_DIR, "before continuing.\n\n")
    Sys.sleep(5)
  }
}
.check_cache()

# ==============================================================================
# Main loop
# ==============================================================================

plan(multisession, workers = N_WORKERS)

t_start      <- Sys.time()
all_summary  <- list()
timing_rows  <- list()

for (j in seq_along(INSTR_GRID)) {
  instr <- INSTR_GRID[[j]]
  cache_path <- file.path(CACHE_DIR,
                          sprintf("instr%02d_%s.rds", j,
                                  gsub("[^A-Za-z0-9_+-]", "_", instr$name)))
  
  if (file.exists(cache_path)) {
    cat(sprintf("[%2d/%2d] %s  (cached)\n", j, length(INSTR_GRID), instr$name))
    cell <- readRDS(cache_path)
  } else {
    sign_label <- if (instr$cor_zu < 0) "neg" else "pos"
    cat(sprintf("[%2d/%2d] %s  cor_zu = %+0.2f  truth = %-8s  sign = %s\n",
                j, length(INSTR_GRID), instr$name, instr$cor_zu,
                instr$truth, sign_label))
    rep_seeds <- MASTER_SEED + 1000L * j + seq_len(R_PER_INSTR)
    t_cell    <- Sys.time()
    
    rep_results <- future_lapply(
      rep_seeds, one_rep,
      instr          = instr,
      worker_files   = WORKER_FILES,   # [FAST-2]
      lambda_pos     = LAMBDA_POS,
      lambda_neg     = LAMBDA_NEG,
      n_sample       = N,
      cor_zx_default = COR_ZX,
      alpha          = ALPHA,
      rxu_pos        = RXU_POS,
      rxu_neg        = RXU_NEG,
      zero_tol       = ZERO_TOL,
      b_mcub         = B_MCUB,        # [FAST-1]
      blarge         = BLARGE,
      b_boot         = B_BOOT,
      future.seed      = TRUE,
      future.packages  = c("MASS", "ivreg", "lmtest", "sandwich",
                           "haven", "dplyr", "Matrix", "numDeriv")
    )
    
    pluck_num <- function(field)
      vapply(rep_results,
             function(r) { v <- r[[field]]; if (is.null(v)) NA_real_ else as.numeric(v[1]) },
             numeric(1))
    pluck_chr <- function(field)
      vapply(rep_results,
             function(r) { v <- r[[field]]; if (is.null(v)) NA_character_ else as.character(v[1]) },
             character(1))
    
    ci_los   <- pluck_num("ci_lo")
    ci_his   <- pluck_num("ci_hi")
    verdicts <- pluck_chr("verdict")
    zero_in  <- vapply(rep_results,
                       function(r) if (is.null(r$zero_in)) NA else as.logical(r$zero_in[1]),
                       logical(1))
    err_vec  <- pluck_chr("err")
    err_n    <- sum(!is.na(err_vec))
    
    cell <- list(
      schema_version = SCHEMA_VERSION,
      instrument     = instr$name,
      cor_zu         = instr$cor_zu,
      truth          = instr$truth,
      sign_label     = sign_label,
      lambda_pos     = LAMBDA_POS,
      lambda_neg     = LAMBDA_NEG,
      ci_los         = ci_los,
      ci_his         = ci_his,
      verdicts       = verdicts,
      zero_in        = zero_in,
      n_failed       = err_n,
      R_completed    = sum(!is.na(zero_in)),
      elapsed_sec    = as.numeric(difftime(Sys.time(), t_cell, units = "secs"))
    )
    saveRDS(cell, cache_path)
    cat(sprintf("    elapsed: %.1f min  (%d failures)\n",
                cell$elapsed_sec / 60, err_n))
  }
  
  p_hat         <- mean(cell$zero_in, na.rm = TRUE)
  mean_ci_lo    <- mean(cell$ci_los,  na.rm = TRUE)
  mean_ci_hi    <- mean(cell$ci_his,  na.rm = TRUE)
  modal_verdict <- {
    tab <- table(cell$verdicts, useNA = "no")
    if (!length(tab)) NA_character_ else names(tab)[which.max(tab)]
  }
  
  if (cell$truth %in% c("A", "C")) {
    misclass     <- sum(cell$verdicts != cell$truth & !is.na(cell$verdicts))
    misclass_str <- sprintf("%d/%d", misclass, cell$R_completed)
    bound_holds  <- (misclass / cell$R_completed) <=
      exp(-2 * cell$R_completed * (1 - 0.95)^2)
    bound_str    <- if (bound_holds) "Yes" else "No"
  } else {
    misclass_str <- "---"; bound_str <- "---"
  }
  
  cat(sprintf("    p_hat = %.3f  CI mean = [%+.3f, %+.3f]  modal = %-3s  misclass = %s  bound = %s\n",
              p_hat, mean_ci_lo, mean_ci_hi,
              modal_verdict, misclass_str, bound_str))
  
  all_summary[[j]] <- data.frame(
    Instrument  = cell$instrument,
    cor_zu      = cell$cor_zu,
    sign        = cell$sign_label,
    p_hat       = round(p_hat, 4),
    CI_lo_mean  = round(mean_ci_lo, 4),
    CI_hi_mean  = round(mean_ci_hi, 4),
    DIAL        = modal_verdict,
    Misclass    = misclass_str,
    Bound_holds = bound_str,
    R_completed = cell$R_completed,
    elapsed_min = round(cell$elapsed_sec / 60, 2),
    stringsAsFactors = FALSE
  )
  timing_rows[[j]] <- data.frame(
    Instrument  = cell$instrument,
    cor_zu      = cell$cor_zu,
    R_completed = cell$R_completed,
    elapsed_min = round(cell$elapsed_sec / 60, 2),
    stringsAsFactors = FALSE
  )
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

# ==============================================================================
# Report
# ==============================================================================

summary_df <- do.call(rbind, all_summary)
rownames(summary_df) <- NULL

cat("\n================================================================================\n")
cat("TABLE 6\n")
cat("================================================================================\n\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nTotal elapsed: %.2f min (%.2f hours)\n",
            elapsed_total, elapsed_total / 60))

utils::write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)
utils::write.csv(do.call(rbind, timing_rows), TIMING_LOG, row.names = FALSE)
cat(sprintf("Wrote %s\n", SUMMARY_CSV))
cat(sprintf("Wrote %s\n", TIMING_LOG))

cat("\nLaTeX rows for tab:mc:\n\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf(
    "$z_{%+0.2f}$ & %+0.2f & %s & %.3f & [%+.3f, %+.3f] & %s & %s & %s \\\\\n",
    r$cor_zu, r$cor_zu, r$sign,
    r$p_hat, r$CI_lo_mean, r$CI_hi_mean,
    r$DIAL, r$Misclass, r$Bound_holds
  ))
}

cat("\nDONE.\n")