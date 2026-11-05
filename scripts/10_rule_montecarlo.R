# ==============================================================================
# 10_rule_montecarlo.R  -  Monte Carlo of the bootstrap validity rule
#                          (tab:rule_mc, app:mc)
# ==============================================================================
#
# PURPOSE
#   tab:mc (05_montecarlo_validation.R) and fig:power_curves (07) report
#   how often a single MCUB set contains zero. DIAL's validity gate is the
#   bootstrap rule p_hat_B >= tau_V. This script runs the rule itself on
#   simulated data: for each cell rho_ZU it draws R_REPS samples from the
#   J-P DGP (n = 800), applies the rule with B_RULE resamples, and records
#   p_hat_B and the verdict.
#
# DESIGN
#   Same DGP, domains, alpha, margin and seeds as 05, so the first R_REPS
#   samples of each cell are the samples of tab:mc. The script reads 05's
#   cached single-set verdicts for those samples (run 05 first; without its
#   cache that column is NA), so the two decisions are compared on
#   identical data.
#     Cells   rho_ZU = 0 (valid), 0.12 (boundary), 0.20 (invalid)
#     Rule    B_RULE = 200 resamples (the paper's tables use 500); each
#             resample computes a full MCUB set (500 bisection draws, 499
#             coverage draws, covariance from 200 raw-data draws: the
#             defaults of foundation/bootstrap_rule.R)
#     Verdict p_hat_B >= tau_V = 0.95, with the margin ZERO_TOL = 0.05 (the
#             simulation convention) and with strict containment, both from
#             the same resamples
#
# CACHE AND RESUMING
#   results/boot_rule/  every replication's resampled sets (bootstrap_rule.R)
#   results/rule_mc/    the per-cell results, saved after each replication
#   An interrupted run resumes at the last finished replication.
#
# RUNTIME  about 8.5 hours in all on 7 cores (504 minutes for the three cells).
#
# OUTPUT
#   results/rule_mc_results.csv and LaTeX rows for tab:rule_mc
#
# RUN
#   Rscript scripts/10_rule_montecarlo.R
# ==============================================================================

rm(list = ls())

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org/")
}
install_if_missing(c("future.apply", "MASS", "Matrix", "dplyr", "haven",
                     "lmtest", "sandwich", "ivreg"))
suppressPackageStartupMessages({ library(future.apply) })

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

N             <- 800L
R_REPS        <- 100L          # replications per cell
B_RULE        <- 200L          # bootstrap resamples per replication
CELLS         <- c(0.00, 0.12, 0.20)
TARGET_COR_XU <- 0.30
COR_ZX        <- 0.25
RXU_POS       <- c( 0.00, 0.80)
RXU_NEG       <- c(-0.80, 0.00)
ALPHA         <- 0.10
ZERO_TOL      <- 0.05
TAU_V         <- 0.95
PI_STAR       <- 0.05
MASTER_SEED   <- 4242L         # as in 05
N_WORKERS     <- min(7L, max(1L, parallel::detectCores(logical = FALSE) - 1L))

# ==============================================================================
# PATHS AND FOUNDATION
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
.repo_root     <- normalizePath(file.path(.script_dir, ".."))
FOUNDATION_DIR <- file.path(.repo_root, "foundation")
RESULTS_DIR    <- file.path(.repo_root, "results")
MC_DIR         <- file.path(RESULTS_DIR, "rule_mc")
T6_CACHE       <- file.path(RESULTS_DIR, "table6_grid_lean")   # 05's cells
for (dd in c(RESULTS_DIR, MC_DIR)) if (!dir.exists(dd)) dir.create(dd, recursive = TRUE)

.N_config <- N
assign("DIAL_SKIP_AUTORUN", TRUE, envir = .GlobalEnv)
for (f in c("Fin_Empirical5_clean.R", "jp_dgp.R", "bootstrap_rule.R")) {
  p <- file.path(FOUNDATION_DIR, f)
  if (!file.exists(p)) stop("Foundation file not found: ", p)
  invisible(capture.output(source(p, local = FALSE)))
}
if (!identical(N, .N_config)) stop("A sourced file changed N from ", .N_config, " to ", N, ".")
need <- c("draw_jp_sample", "calibrate_lambda", "JP_INSTRUMENT_GRID",
          "mcub_bootstrap_rule", "dial_raw_spec", "dial_prepare_xyz")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Foundation did not provide: ", paste(miss, collapse = ", "))
options(dial.boot_B = B_RULE, dial.tau2 = TAU_V, dial.zero_tol = ZERO_TOL)

cat("================================================================================\n")
cat("Monte Carlo of the bootstrap validity rule (tab:rule_mc)\n")
cat("================================================================================\n")
cat(sprintf("  repo_root = %s\n", .repo_root))
cat(sprintf("  n = %d, R = %d per cell, B_rule = %d, cells rho_ZU = %s\n",
            N, R_REPS, B_RULE, paste(sprintf("%+.2f", CELLS), collapse = ", ")))
cat(sprintf("  alpha = %.2f, tau_V = %.2f, ZERO_TOL = %.2f, workers = %d\n\n",
            ALPHA, TAU_V, ZERO_TOL, N_WORKERS))

cal_pos <- calibrate_lambda(n = N, target_cor_xu = TARGET_COR_XU,
                            cor_zx_default = COR_ZX, seed = MASTER_SEED,
                            sign = "positive", verbose = FALSE)
cal_neg <- calibrate_lambda(n = N, target_cor_xu = TARGET_COR_XU,
                            cor_zx_default = COR_ZX, seed = MASTER_SEED,
                            sign = "negative", verbose = FALSE)
cat(sprintf("  lambda_+ = %+.4f, lambda_- = %+.4f\n\n", cal_pos$lambda, cal_neg$lambda))

truth_label <- function(cor_zu) {            # as in 05
  r <- abs(cor_zu)
  if (r <= 0.04 + 1e-9) return("A")
  if (r >= 0.16 - 1e-9) return("C")
  "Boundary"
}
ucb95 <- function(x, R) {
  if (R <= 0) return(NA_real_)
  if (x >= R) return(1)
  stats::qbeta(0.95, x + 1, R - x)
}

# Single-set verdicts of 05 for the same samples, if its cache is present
single_set <- function(j, cor_zu, R) {
  f <- file.path(T6_CACHE, sprintf("instr%02d_%s.rds", j,
                                   gsub("[^A-Za-z0-9_+-]", "_", sprintf("z_%+0.2f", cor_zu))))
  if (!file.exists(f)) return(rep(NA, R))
  cell <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(cell) || length(cell$zero_in) < R) return(rep(NA, R))
  as.logical(cell$zero_in[seq_len(R)])
}

# ==============================================================================
# MAIN LOOP
# ==============================================================================

future::plan(future::multisession, workers = N_WORKERS)   # used by the rule
GRID    <- JP_INSTRUMENT_GRID
t_start <- Sys.time()
summary_rows <- list()

for (cz in CELLS) {
  j <- which.min(abs(GRID - cz))
  if (abs(GRID[j] - cz) > 1e-9) stop("rho_ZU = ", cz, " is not on JP_INSTRUMENT_GRID")
  neg <- cz < 0
  lam <- if (neg) cal_neg$lambda else cal_pos$lambda
  rxu <- if (neg) RXU_NEG else RXU_POS
  rep_seeds <- MASTER_SEED + 1000L * j + seq_len(R_REPS)     # as in 05

  cell_file <- file.path(MC_DIR, sprintf("cell_%+.2f.rds", cz))
  key  <- list(N = N, B_RULE = B_RULE, ALPHA = ALPHA, lambda = lam, rxu = rxu)
  cell <- if (file.exists(cell_file)) tryCatch(readRDS(cell_file), error = function(e) NULL)
  if (is.null(cell) || !identical(cell$key, key))
    cell <- list(cor_zu = cz, key = key, reps = vector("list", R_REPS))
  if (length(cell$reps) < R_REPS) length(cell$reps) <- R_REPS

  cat(sprintf("[cell %+.2f] truth = %s, domain (%.1f, %.1f)\n",
              cz, truth_label(cz), rxu[1], rxu[2]))
  t_cell <- Sys.time(); n_new <- 0L
  for (r in seq_len(R_REPS)) {
    if (!is.null(cell$reps[[r]])) next
    sim <- draw_jp_sample(n = N, lambda = lam, seed = rep_seeds[r],
                          cor_zu_grid = c(cz), cor_zx_default = COR_ZX)
    d_raw <- data.frame(Y = sim$y, X = sim$x, Z = sim$z[, 1L])
    spec  <- dial_raw_spec(d_raw, "Y", "X", "Z")
    df    <- dial_prepare_xyz(spec$data, "Y", "X", "Z")   # demean, then project
    br <- mcub_bootstrap_rule(df, rxu_range = rxu, alpha = ALPHA,
                              seed = rep_seeds[r], B_rule = B_RULE,
                              zero_tol = ZERO_TOL, raw = spec, verbose = FALSE)
    ok <- stats::complete.cases(br$ends)
    p_strict <- if (any(ok)) mean(br$ends[ok, "lo"] <= 0 & br$ends[ok, "hi"] >= 0) else NA_real_
    cell$reps[[r]] <- list(seed = rep_seeds[r], p_tol = br$p_hat,
                           p_strict = p_strict, n_ok = sum(ok))
    saveRDS(cell, cell_file)
    n_new <- n_new + 1L
    done  <- sum(!vapply(cell$reps, is.null, logical(1)))
    el    <- as.numeric(difftime(Sys.time(), t_cell, units = "mins"))
    cat(sprintf("  rep %3d/%d  p_hat_B = %.3f (strict %.3f) -> %s   %.1f min, about %.0f min left in this cell\n",
                r, R_REPS, br$p_hat, p_strict,
                if (isTRUE(br$p_hat >= TAU_V)) "Valid" else "Invalid",
                el, el / n_new * (R_REPS - done)))
  }

  reps   <- cell$reps[!vapply(cell$reps, is.null, logical(1))]
  p_tol  <- vapply(reps, `[[`, numeric(1), "p_tol")
  p_str  <- vapply(reps, `[[`, numeric(1), "p_strict")
  v_tol  <- p_tol >= TAU_V
  v_str  <- p_str >= TAU_V
  single <- single_set(j, GRID[j], R_REPS)
  truth  <- truth_label(cz)
  R_done <- length(reps)
  mis <- if (truth == "A") sum(!v_tol) else if (truth == "C") sum(v_tol) else NA_integer_
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    cor_zu = cz, truth = truth, R = R_done,
    valid_rule = mean(v_tol), valid_rule_strict = mean(v_str),
    valid_single_set = mean(single, na.rm = TRUE),
    p_hat_mean = mean(p_tol), p_hat_q10 = stats::quantile(p_tol, 0.10, names = FALSE),
    p_hat_q90 = stats::quantile(p_tol, 0.90, names = FALSE),
    misclass = if (is.na(mis)) NA_character_ else sprintf("%d/%d", mis, R_done),
    ucb95 = if (is.na(mis)) NA_real_ else ucb95(mis, R_done),
    stringsAsFactors = FALSE)
  cat(sprintf("  => Valid by the rule: %.3f (strict %.3f); single set: %.3f; mean p_hat_B %.3f\n\n",
              mean(v_tol), mean(v_str), mean(single, na.rm = TRUE), mean(p_tol)))
}
future::plan(future::sequential)

out <- do.call(rbind, summary_rows)
utils::write.csv(out, file.path(RESULTS_DIR, "rule_mc_results.csv"), row.names = FALSE)

cat("================================================================================\n")
cat("tab:rule_mc\n")
cat("================================================================================\n")
print(out, digits = 3, row.names = FALSE)
cat(sprintf("\nTotal elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

cat("\nLaTeX rows (rho_ZU & truth & single set & rule & rule, strict & mean p_hat_B & misclass. & 95% UCB):\n\n")
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  rho <- if (abs(r$cor_zu) < 1e-9) "$\\phantom{+}0.00$" else sprintf("$%+.2f$", r$cor_zu)
  ss  <- if (is.finite(r$valid_single_set)) sprintf("$%.2f$", r$valid_single_set) else "---"
  cat(sprintf("%s & %s & %s & $%.2f$ & $%.2f$ & $%.3f$ & %s & %s \\\\\n",
              rho, if (r$truth == "Boundary") "Boundary" else r$truth,
              ss, r$valid_rule, r$valid_rule_strict, r$p_hat_mean,
              if (is.na(r$misclass)) "---" else sprintf("$%s$", r$misclass),
              if (is.na(r$ucb95)) "---" else sprintf("$%.3f$", r$ucb95)))
}
cat(sprintf("\nWrote %s\n", file.path(RESULTS_DIR, "rule_mc_results.csv")))
cat("\nDONE.\n")
