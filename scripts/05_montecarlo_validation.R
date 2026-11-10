# ==============================================================================
# 05_montecarlo_validation.R   --   Monte Carlo study of single MCUB sets
#                                   and the MCUB-share rule (app:mc)
# ------------------------------------------------------------------------------
# Tables are referenced by LaTeX label, not by number.
#
# NOTE  The paper mentions this study in app:mc but does not tabulate it;
#       tab:mc below names this script's own table, and prop:finite_sample
#       and "bootstrap validity rule" refer to the MCUB-share rule.
#
# PURPOSE
#   For each of the 11 instruments of JP_INSTRUMENT_GRID (target cor(z, u)
#   from -0.20 to 0.20 in steps of 0.04), draws R_PER_INSTR = 1,000 samples
#   of n = 800 from the Jones-Pewsey DGP of app:dgp_jp (foundation/jp_dgp.R),
#   computes one MCUB set per replication, and scores the verdict (A if
#   zero lies in the set within ZERO_TOL, else C) against the truth label.
#   Per-row sign elicitation: an instrument with cor(z, u) < 0 is drawn
#   with lambda < 0 and tested under D = (-0.8, 0); otherwise lambda > 0
#   and D = (0, 0.8).
#
# IMPLEMENTATION
#   - Each replication calls CIhybrid() once (B_MCUB = 500 draws per
#     bisection step, BLARGE = 499 coverage draws). The covariance of the
#     correlations comes from B_BOOT = 200 raw-data bootstrap draws
#     (dial_cov_corr_boot_raw() in foundation/bootstrap_rule.R): each draw
#     resamples the simulated (y, x, z) rows and redoes the preparation
#     (demean, then project the instrument on the demeaned (x, y)).
#     check_compatibility_simple(), pvalue_mcub_zero_fast() and
#     ci_simple_union() are not called; the verdict uses the CIhybrid
#     endpoints only.
#   - The main session and the workers source Fin_Empirical5_clean.R,
#     jp_dgp.R and bootstrap_rule.R. Fin_sim3_clean.R is not needed, and
#     its top-level code sets N = 800; the script stops if a sourced file
#     changes N.
#   - Truth labels: at n = 800 a single MCUB set accepts |rho_ZU| = 0.12
#     in 98-99% of replications, so the frontier of a single set lies
#     between 0.12 and 0.16. truth_label() scores 0.08 and 0.12 as
#     boundary and starts C at 0.16. The summary scores each cell against
#     the current truth_label(), so a change to the labels takes effect
#     without recomputing cells. The LaTeX rows print the zero row as
#     \phantom{+}0.00, as in the paper.
#
# ERROR-RATE BOUND
#   The script reports the exact one-sided 95% (Clopper-Pearson) upper
#   confidence bound on each scored cell's misclassification rate and
#   compares it with the target pi* = 0.05. With 0 errors in R = 1,000
#   replications the bound is 1 - 0.05^(1/1000) = 0.0030.
#
# PATHS
#   All foundation files are resolved against a repo root detected from
#   the script's own location, so the script runs from the repo root or
#   from scripts/.
#
# RUNTIME (n = 800, 7 workers)
#   R = 50   sanity run           ~ 10 minutes
#   R = 1000 full tab:mc          ~ 3 hours (the eleven cells took 5 to 29
#                                   minutes each, 174 minutes in all)
#
# CACHE
#   results/table6_grid_lean/, one .rds per cell. Each cell records
#   SCHEMA_VERSION and n_sample, and a cached cell is reused only if both
#   are current. (The cache and CSV names keep an earlier table number.)
#
# OUTPUT
#   results/table6_results.csv
#   results/table6_grid_lean/*.rds
#   results/table6_timing.csv
#
# RUN
#   Rscript scripts/05_montecarlo_validation.R
# ==============================================================================

rm(list = ls())

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org/")
}
install_if_missing(c("future.apply", "ivreg", "lmtest", "sandwich",
                     "haven", "dplyr", "Matrix", "numDeriv"))
suppressPackageStartupMessages({ library(future.apply) })

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
    if (length(f) && nchar(f)) dirname(normalizePath(f))
    else getwd()
  }
)
.repo_root <- normalizePath(file.path(.script_dir, ".."))

FOUNDATION_DIR <- file.path(.repo_root, "foundation")
RESULTS_DIR    <- file.path(.repo_root, "results")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ==============================================================================
# CONFIG
# ==============================================================================

N             <- 800L         # rows per simulated sample
R_PER_INSTR   <- 1000L        # set to 50L for a sanity check
TARGET_COR_XU <- 0.30
COR_ZX        <- 0.25
ZERO_TOL      <- 0.05         # verdict: 0 in CI iff lo <= ZERO_TOL & hi >= -ZERO_TOL
ALPHA         <- 0.10
PI_STAR       <- 0.05         # target misclassification rate
RXU_POS       <- c( 0.00,  0.80)
RXU_NEG       <- c(-0.80,  0.00)
N_WORKERS     <- min(7L, max(1L, parallel::detectCores(logical = FALSE) - 1L))
MASTER_SEED   <- 4242L

# CIhybrid tuning (one call per rep)
B_MCUB  <- 500L   # per-bisection draw count
BLARGE  <- 499L   # c_bd and final coverage draws
B_BOOT  <- 200L   # covariance bootstrap size (raw-data draws)

SCHEMA_VERSION <- 10L  # cache schema; cells with an older schema are recomputed

CACHE_DIR   <- file.path(RESULTS_DIR, "table6_grid_lean")
SUMMARY_CSV <- file.path(RESULTS_DIR, "table6_results.csv")
TIMING_LOG  <- file.path(RESULTS_DIR, "table6_timing.csv")
if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)

# Workers source Fin_Empirical5_clean.R, jp_dgp.R and bootstrap_rule.R
# (for dial_cov_corr_boot_raw()); Fin_sim3_clean.R is not needed.
WORKER_FILES <- file.path(FOUNDATION_DIR,
                          c("Fin_Empirical5_clean.R",
                            "jp_dgp.R",
                            "bootstrap_rule.R"))

# The main session sources the same files. calibrate_lambda() and
# JP_INSTRUMENT_GRID come from jp_dgp.R; Fin_sim3_clean.R is
# not needed, and its top-level code would reset N to 800.
MAIN_FILES   <- WORKER_FILES

for (.f in c(WORKER_FILES, MAIN_FILES)) {
  if (!file.exists(.f))
    stop("Foundation file not found: ", .f,
         "\n  Expected layout: <repo_root>/foundation/<name>",
         "\n  .repo_root currently = ", .repo_root,
         "\n  If this is wrong, override by setting .repo_root before sourcing.")
}
rm(.f)

cat("================================================================================\n")
cat("DIAL Monte Carlo validation (tab:mc): J-P DGP, per-row sign elicitation\n")
cat("================================================================================\n")
cat(sprintf("  repo_root       = %s\n", .repo_root))
cat(sprintf("  n = %d,  R per instrument = %d\n", N, R_PER_INSTR))
cat(sprintf("  |target cor(x,u)| = %.2f,  cor(z,x) = %.2f\n",
            TARGET_COR_XU, COR_ZX))
cat(sprintf("  D pos = [%.2f, %.2f]   D neg = [%.2f, %.2f]\n",
            RXU_POS[1], RXU_POS[2], RXU_NEG[1], RXU_NEG[2]))
cat(sprintf("  alpha = %.2f,  ZERO_TOL = %.2f,  pi* = %.2f,  master seed = %d\n",
            ALPHA, ZERO_TOL, PI_STAR, MASTER_SEED))
cat(sprintf("  B_mcub = %d,  Blarge = %d,  B_boot = %d\n",
            B_MCUB, BLARGE, B_BOOT))
cat(sprintf("  Workers = %d,  cache dir: %s\n\n", N_WORKERS, CACHE_DIR))

# ==============================================================================
# Load foundation in main session (absolute paths; no base_dir prefix)
# ==============================================================================

load_foundation <- function(files) {
  assign("DIAL_SKIP_AUTORUN", TRUE, envir = .GlobalEnv)
  for (f in files) {
    if (!file.exists(f)) stop("Foundation file not found: ", f)
    invisible(capture.output(source(f, local = FALSE)))
  }
  required <- c("CIhybrid", "g_xu_safe", "local_compute_gradient_safe",
                "estimate_cov_corr_boot", "dial_cov_corr_boot_raw",
                "draw_jp_sample", "calibrate_lambda", "JP_INSTRUMENT_GRID")
  miss <- setdiff(required, ls(envir = .GlobalEnv))
  if (length(miss))
    stop("Foundation source did not provide: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

cat("Sourcing foundation in main session... ")
.N_config <- N    # checked after sourcing
load_foundation(MAIN_FILES)
if (!identical(N, .N_config))
  stop("A sourced file changed N from ", .N_config, " to ", N, ".")
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

# At n = 800 a single MCUB set does not detect |rho_ZU| = 0.12 (it
# accepts in 98-99% of replications), so 0.08 and 0.12 form the boundary
# region and C starts at 0.16 (the grey region of ass:class).
truth_label <- function(cor_zu) {
  r <- abs(cor_zu)
  if (r <= 0.04 + 1e-9) return("A")
  if (r >= 0.16 - 1e-9) return("C")
  "Boundary"
}

INSTR_GRID <- lapply(seq_along(GRID), function(j) {
  list(idx    = j,
       name   = sprintf("z_%+0.2f", GRID[j]),
       cor_zu = GRID[j],
       truth  = truth_label(GRID[j]))
})

# ==============================================================================
# Lean one_rep: CIhybrid only, no pvalue_mcub_zero_fast
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

  # Workers source WORKER_FILES only (not Fin_sim3_clean.R)
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

  # 4. Bootstrap covariance -- one call at reduced B. Each
  #    draw resamples the simulated rows and redoes step 2 (demean, then
  #    project the instrument on the demeaned (x, y)).
  deltaSigma <- tryCatch(
    dial_cov_corr_boot_raw(
      dial_raw_spec(data.frame(Y = sim$y, X = sim$x, Z = sim$z[, 1L]),
                    "Y", "X", "Z"),
      B = b_boot, seed = rep_seed + 2L),
    error = function(e) NULL
  )
  if (is.null(deltaSigma)) return(empty("dial_cov_corr_boot_raw failed"))

  # 5. Grid, g, Jacobian
  r_grid <- seq(rxu[1], rxu[2], length.out = 50L)
  g      <- function(dm) g_xu_safe(r_grid, dm)
  A      <- t(vapply(r_grid, function(r_xu)
    as.numeric(local_compute_gradient_safe(
      r_xu, deltahat[1], deltahat[2], deltahat[3])),
    numeric(3L)))
  Al <- A; Au <- A

  # 6. Single CIhybrid call
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
# Exact one-sided 95% upper confidence bound on an error rate
# (Clopper-Pearson): x errors in R trials.
# ==============================================================================

ucb95 <- function(x, R) {
  if (R <= 0) return(NA_real_)
  if (x >= R) return(1)
  stats::qbeta(0.95, x + 1, R - x)
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
    cat("\n*** NOTE: stale cache detected (schema", x$schema_version,
        "< required", SCHEMA_VERSION, ")\n")
    cat("Stale cells in", CACHE_DIR, "will be recomputed and overwritten.\n\n")
  }
}
.check_cache()

# A cached cell is reused only if it was written with the current
# schema and the current sample size (a cell without n_sample used n = 800).
.cell_reusable <- function(cell) {
  if (is.null(cell) || is.null(cell$schema_version) ||
      cell$schema_version < SCHEMA_VERSION) return(FALSE)
  cell_n <- if (is.null(cell$n_sample)) 800L else as.integer(cell$n_sample)
  identical(cell_n, N)
}

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

  # Reuse a cached cell only if its schema and n are current.
  cell <- if (file.exists(cache_path))
    tryCatch(readRDS(cache_path), error = function(e) NULL) else NULL
  if (.cell_reusable(cell)) {
    cat(sprintf("[%2d/%2d] %s  (cached)\n", j, length(INSTR_GRID), instr$name))
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
      worker_files   = WORKER_FILES,
      lambda_pos     = LAMBDA_POS,
      lambda_neg     = LAMBDA_NEG,
      n_sample       = N,
      cor_zx_default = COR_ZX,
      alpha          = ALPHA,
      rxu_pos        = RXU_POS,
      rxu_neg        = RXU_NEG,
      zero_tol       = ZERO_TOL,
      b_mcub         = B_MCUB,
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
      n_sample       = N,
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

  # Score against the current truth label, so a change to
  # truth_label() takes effect without recomputing cached cells.
  truth <- instr$truth
  if (truth %in% c("A", "C")) {
    misclass     <- sum(cell$verdicts != truth & !is.na(cell$verdicts))
    misclass_str <- sprintf("%d/%d", misclass, cell$R_completed)
    ucb          <- ucb95(misclass, cell$R_completed)
    ucb_str      <- sprintf("%.3f", ucb)
    below_str    <- if (!is.na(ucb) && ucb < PI_STAR) "Yes" else "No"
  } else {
    misclass_str <- "---"; ucb_str <- "---"; below_str <- "---"
  }

  cat(sprintf("    p_hat = %.3f  CI mean = [%+.3f, %+.3f]  modal = %-3s  misclass = %s  UCB95 = %s\n",
              p_hat, mean_ci_lo, mean_ci_hi,
              modal_verdict, misclass_str, ucb_str))

  all_summary[[j]] <- data.frame(
    Instrument  = cell$instrument,
    cor_zu      = cell$cor_zu,
    sign        = cell$sign_label,
    truth       = truth,
    p_hat       = round(p_hat, 4),
    CI_lo_mean  = round(mean_ci_lo, 4),
    CI_hi_mean  = round(mean_ci_hi, 4),
    DIAL        = modal_verdict,
    Misclass    = misclass_str,
    UCB95       = ucb_str,
    Below_pi    = below_str,
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
cat("tab:mc\n")
cat("================================================================================\n\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nTotal elapsed: %.2f min (%.2f hours)\n",
            elapsed_total, elapsed_total / 60))

utils::write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)
utils::write.csv(do.call(rbind, timing_rows), TIMING_LOG, row.names = FALSE)
cat(sprintf("Wrote %s\n", SUMMARY_CSV))
cat(sprintf("Wrote %s\n", TIMING_LOG))

# Columns: Instrument & rho_ZU & D & p-bar & DIAL & Misclass. & 95% UCB
cat("\nLaTeX rows for tab:mc:\n\n")
for (i in seq_len(nrow(summary_df))) {
  r     <- summary_df[i, ]
  d_str <- if (r$sign == "neg") "$(-0.8,\\,0)$" else "$(0,\\,0.8)$"
  dial  <- if (r$truth == "Boundary") "Boundary" else sprintf("\\textbf{%s}", r$DIAL)
  mis   <- if (r$Misclass == "---") "---" else sprintf("$%s$", r$Misclass)
  ucb   <- if (r$UCB95 == "---") "---" else sprintf("$%s$", r$UCB95)
  rho   <- if (abs(r$cor_zu) < 1e-9) "$\\phantom{+}0.00$" else sprintf("$%+0.2f$", r$cor_zu)
  cat(sprintf(
    "$Z_{%+0.2f}$ & %s & %s & $%.3f$ & %s & %s & %s \\\\\n",
    r$cor_zu, rho, d_str, r$p_hat, dial, mis, ucb
  ))
}

cat("\nDONE.\n")
