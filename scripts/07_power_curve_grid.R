# =============================================================================
# DIAL_NeurIPS_cp_power_curve_run4_FULL.R
#
# Full power-curve run, designed using validate run3e results.
#
# DELTA GRID (35 points) -- designed from validate observations:
#   * Validate showed both n=1107 and n=50000 transition between delta=0.10
#     and delta=0.15 (curves bracket the same cell). Plug-in upper bound hits
#     zero at delta ~ 0.085 (population-driven). Predicted threshold gap
#     between n=1107 and n=50000 is ~0.02 in delta-space.
#   * Dense grid of 18 points across [0.07, 0.18] to resolve the transition.
#   * 9 anchor points at smaller delta to confirm CP=1 plateau.
#   * 8 anchor points at larger delta to confirm CP=0 plateau.
#
# N GRID: 5 sample sizes from 1107 to 50000, log-spaced.
# R_MC schedule: more reps at small n where we'll plot, fewer at large n.
#
# RUNTIME: estimated 15-25 hours wall time on 7 workers.
#   - n=1107  cells: ~28s/rep, 500 reps / 7 workers ~ 33min/cell × 35 cells = 19hr
#   - n=50000 cells: ~28s/rep, 100 reps / 7 workers ~ 7min/cell × 35 cells = 4hr
#   (small n dominates because we use more reps there)
#
# KEY FIX FROM run3e: existence-based worker setup gate (not flag-based).
# =============================================================================

rm(list = ls())

# ---- Paths --------------------------------------------------------------
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

SIM3_PATH    <- file.path(.repo_root, 'foundation', 'Fin_sim3_clean.R')
EMP5_PATH    <- file.path(.repo_root, 'foundation', 'Fin_Empirical5_clean.R')
WRAPPER_PATH <- file.path(.repo_root, 'scripts',    'cr_wrapper.R')

stopifnot(file.exists(WRAPPER_PATH), file.exists(SIM3_PATH),
          file.exists(EMP5_PATH))

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org/")
}
install_if_missing(c("MASS", "future.apply", "dplyr", "tidyr",
                     "numDeriv", "Matrix", "ivreg", "lmtest", "sandwich",
                     "haven", "fGarch"))

suppressPackageStartupMessages({
  library(MASS); library(future.apply)
})

# ---- Master setup -------------------------------------------------------
source_safe <- function(path, envir = globalenv()) {
  lines <- readLines(path, warn = FALSE)
  cut_at <- grep("SECTION 12|Quick-?start reference", lines, perl = TRUE)
  if (length(cut_at)) {
    eval(parse(text = paste(lines[seq_len(cut_at[1] - 1L)], collapse = "\n")),
         envir = envir)
  } else source(path, local = FALSE)
  invisible(NULL)
}

DIAL_SKIP_AUTORUN <- TRUE
source(SIM3_PATH); source(EMP5_PATH); source_safe(WRAPPER_PATH)
stopifnot(exists("cr_test_unified"))
cat("Master loaded.\n")

plan(sequential)
plan(multisession, workers = 7L)
cat("Parallel plan: multisession with 7 workers\n")

# =============================================================================
# Configuration
# =============================================================================
N_GRID <- c(1107L, 1534L, 5000L, 10000L, 50000L)

DELTA_GRID <- local({
  low_anchors  <- c(0.001, 0.005, 0.010, 0.020, 0.030, 0.045, 0.055, 0.060)  # 8
  pre_trans    <- seq(0.07, 0.10, length.out = 4)                            # 4
  trans_dense  <- seq(0.105, 0.155, length.out = 11)                         # 11
  post_trans   <- seq(0.165, 0.20, length.out = 5)                           # 5
  high_anchors <- c(0.22, 0.25, 0.28, 0.30, 0.32, 0.35, 0.40)                # 7
  sort(unique(c(low_anchors, pre_trans, trans_dense, post_trans, high_anchors)))
})
stopifnot(length(DELTA_GRID) == 35)

R_MC_BY_N <- c(
  "1107"  = 500L,
  "1534"  = 500L,
  "5000"  = 300L,
  "10000" = 200L,
  "50000" = 100L
)

OUT_DIR     <- file.path(normalizePath("."), "cp_grid_full")
SUMMARY_CSV <- file.path(normalizePath("."), "cp_grid_summary_full.csv")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

TAU2 <- 0.95; ALPHA <- 0.05
RHO_DU <- 0.50; COR_ZD <- 0.25; BETA_TRUE <- 0.30
RXU_RANGE <- c(0, 0.80)
B_MCUB <- 199L; BLARGE_MCUB <- 1999L

cat(sprintf("\nN grid     : %s\n", paste(N_GRID, collapse = ", ")))
cat(sprintf("delta grid : %d points in [%.4f, %.4f]\n",
            length(DELTA_GRID), min(DELTA_GRID), max(DELTA_GRID)))
cat(sprintf("Total cells: %d\n", length(N_GRID) * length(DELTA_GRID)))
total_reps <- sum(R_MC_BY_N) * length(DELTA_GRID)
cat(sprintf("Total reps : %d\n", total_reps))
cat(sprintf("Output dir : %s\n", OUT_DIR))
cat(sprintf("Cache      : preserved between runs (caching is safe in this script)\n\n"))

# ---- DGP ----------------------------------------------------------------
construct_instrument <- function(Dtilde, Utilde, rho, a, b, seed_w = 1L) {
  set.seed(seed_w); n <- length(Dtilde)
  V <- (Utilde - rho * Dtilde) / sqrt(1 - rho^2)
  W_raw <- rnorm(n); W <- resid(lm(W_raw ~ Dtilde + V)); W <- W / sd(W)
  c1 <- a; c2 <- (b - a * rho) / sqrt(1 - rho^2)
  c3_sq <- 1 - c1^2 - c2^2; if (c3_sq <= 0) stop("Infeasible.")
  c1 * Dtilde + c2 * V + sqrt(c3_sq) * W
}
simulate_one <- function(n, delta, rep_seed) {
  set.seed(rep_seed)
  Sigma <- matrix(c(1, RHO_DU, RHO_DU, 1), 2, 2)
  DU <- MASS::mvrnorm(n, c(0, 0), Sigma)
  Z <- construct_instrument(DU[,1], DU[,2], rho = RHO_DU,
                            a = COR_ZD, b = delta, seed_w = rep_seed + 1L)
  data.frame(x = DU[,1], y = BETA_TRUE * DU[,1] + DU[,2], z = Z)
}
parse_iv <- function(s) {
  if (is.null(s) || length(s) == 0 || is.na(s)[1]) return(c(NA_real_, NA_real_))
  nums <- as.numeric(strsplit(gsub("[\\[\\]]", "", as.character(s),
                                   perl = TRUE), "[,;]\\s*")[[1]])
  if (length(nums) >= 2) nums[1:2] else c(NA_real_, NA_real_)
}
sc <- function(x, fill = NA_real_) if (is.null(x) || length(x) == 0) fill else x[1]
scc <- function(x, fill = NA_character_) if (is.null(x) || length(x) == 0) fill else as.character(x[1])

# ---- one_rep with EXISTENCE-BASED setup gate ---------------------------
one_rep <- function(rep_seed, n_val, delta, rxu_range, alpha,
                    b_mcub, blarge_mcub, sim3, emp5, wrap) {
  empty_result <- function(err_msg = NA_character_) {
    list(p = NA_real_, zero_in = NA,
         plug_lo = NA_real_, plug_hi = NA_real_,
         ci_lo   = NA_real_, ci_hi   = NA_real_,
         err     = err_msg)
  }
  
  if (!exists("cr_test_unified", envir = globalenv(), inherits = FALSE)) {
    setup_err <- tryCatch({
      assign("DIAL_SKIP_AUTORUN", TRUE, envir = globalenv())
      source(sim3, local = FALSE)
      source(emp5, local = FALSE)
      lines <- readLines(wrap, warn = FALSE)
      cut_at <- grep("SECTION 12|Quick-?start reference", lines, perl = TRUE)
      if (length(cut_at)) {
        eval(parse(text = paste(lines[seq_len(cut_at[1] - 1L)], collapse = "\n")),
             envir = globalenv())
      } else {
        source(wrap, local = FALSE)
      }
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(setup_err))
      return(empty_result(paste("worker setup failed:", setup_err)))
    if (!exists("cr_test_unified", envir = globalenv(), inherits = FALSE))
      return(empty_result("crtu still missing after setup"))
  }
  
  df <- tryCatch(simulate_one(n_val, delta, rep_seed),
                 error = function(e) NULL)
  if (is.null(df)) return(empty_result("sim error"))
  
  cr_fn <- get("cr_test_unified", envir = globalenv())
  
  res <- tryCatch(
    cr_fn(df = df, rxu_range = rxu_range, alpha = alpha,
          method = "mcub", label = "cp", i = rep_seed,
          B_mcub = b_mcub, Blarge_mcub = blarge_mcub, verbose = FALSE),
    error = function(e) list(.cr_err = conditionMessage(e))
  )
  if (!is.null(res$.cr_err)) return(empty_result(res$.cr_err))
  
  pin <- parse_iv(res$plug_in); ci <- parse_iv(res$CI_MCUB)
  list(
    p       = as.numeric(res$p_zero),
    zero_in = !is.na(ci[1]) && ci[1] <= 0 && 0 <= ci[2],
    plug_lo = pin[1], plug_hi = pin[2],
    ci_lo   = ci[1],  ci_hi   = ci[2],
    err     = NA_character_
  )
}

WORKER_GLOBALS <- c("simulate_one", "construct_instrument", "parse_iv",
                    "RHO_DU", "COR_ZD", "BETA_TRUE")

# ---- Main loop ----------------------------------------------------------
cell_counter <- 0L
total_cells  <- length(N_GRID) * length(DELTA_GRID)
t_start_global <- Sys.time()

for (n_val in N_GRID) {
  R_MC <- R_MC_BY_N[[as.character(n_val)]]
  for (d_idx in seq_along(DELTA_GRID)) {
    delta <- DELTA_GRID[d_idx]
    cell_counter <- cell_counter + 1L
    
    fname <- sprintf("n%05d_d%02d.rds", n_val, d_idx)
    fpath <- file.path(OUT_DIR, fname)
    if (file.exists(fpath)) {
      cat(sprintf("[%3d/%3d] n=%d d=%.4f  (cached)\n",
                  cell_counter, total_cells, n_val, delta))
      next
    }
    
    t_cell_start <- Sys.time()
    cell_seed_base <- 10000L * n_val + d_idx
    rep_seeds <- cell_seed_base + seq_len(R_MC)
    
    rep_results <- future_lapply(
      rep_seeds, one_rep,
      n_val = n_val, delta = delta,
      rxu_range = RXU_RANGE, alpha = ALPHA,
      b_mcub = B_MCUB, blarge_mcub = BLARGE_MCUB,
      sim3 = normalizePath(SIM3_PATH),
      emp5 = normalizePath(EMP5_PATH),
      wrap = normalizePath(WRAPPER_PATH),
      future.seed = TRUE,
      future.packages = c("MASS", "dplyr", "tidyr", "numDeriv", "Matrix",
                          "ivreg", "lmtest", "sandwich", "haven", "fGarch"),
      future.globals = WORKER_GLOBALS
    )
    
    p_vals  <- vapply(rep_results, function(r) sc(r$p),       numeric(1))
    zero_in <- vapply(rep_results, function(r) sc(r$zero_in, NA), logical(1))
    plug_lo <- vapply(rep_results, function(r) sc(r$plug_lo), numeric(1))
    plug_hi <- vapply(rep_results, function(r) sc(r$plug_hi), numeric(1))
    ci_lo   <- vapply(rep_results, function(r) sc(r$ci_lo),   numeric(1))
    ci_hi   <- vapply(rep_results, function(r) sc(r$ci_hi),   numeric(1))
    errs    <- vapply(rep_results, function(r) scc(r$err),    character(1))
    n_errs  <- sum(!is.na(errs))
    
    cell_result <- list(
      n = n_val, delta = delta, delta_index = d_idx,
      cp_raw          = sc(mean(zero_in,         na.rm = TRUE)),
      cp_decision     = sc(mean(p_vals >= TAU2,  na.rm = TRUE)),
      plug_lo_mean    = sc(mean(plug_lo,         na.rm = TRUE)),
      plug_hi_mean    = sc(mean(plug_hi,         na.rm = TRUE)),
      plug_width_mean = sc(mean(plug_hi-plug_lo, na.rm = TRUE)),
      ci_lo_mean      = sc(mean(ci_lo,           na.rm = TRUE)),
      ci_hi_mean      = sc(mean(ci_hi,           na.rm = TRUE)),
      ci_width_mean   = sc(mean(ci_hi  -ci_lo,   na.rm = TRUE)),
      n_valid_reps    = sum(!is.na(p_vals)),
      n_errors        = n_errs,
      first_error     = if (n_errs > 0) scc(errs[!is.na(errs)][1]) else NA_character_,
      R_MC = R_MC, B_mcub = B_MCUB,
      elapsed_sec = as.numeric(difftime(Sys.time(), t_cell_start, units = "secs"))
    )
    saveRDS(cell_result, fpath)
    
    err_tag <- if (n_errs > 0) sprintf(" [%d errors]", n_errs) else ""
    cat(sprintf(
      "[%3d/%3d] n=%5d d=%.4f  cp=%.3f  plug_w=%.3f  CI_w=%.3f  (%.0fs, R=%d)%s\n",
      cell_counter, total_cells, n_val, delta,
      cell_result$cp_raw, cell_result$plug_width_mean,
      cell_result$ci_width_mean, cell_result$elapsed_sec, R_MC, err_tag))
    if (n_errs > 0)
      cat(sprintf("    first error: %s\n", cell_result$first_error))
  }
}

t_total <- as.numeric(difftime(Sys.time(), t_start_global, units = "mins"))
cat(sprintf("\nDone in %.1f minutes (%.2f hours).\n", t_total, t_total/60))

# ---- Aggregate ----------------------------------------------------------
all_rows <- list(); skipped <- 0L
for (n_val in N_GRID) {
  for (d_idx in seq_along(DELTA_GRID)) {
    fpath <- file.path(OUT_DIR, sprintf("n%05d_d%02d.rds", n_val, d_idx))
    if (!file.exists(fpath)) next
    x <- tryCatch(readRDS(fpath), error = function(e) NULL)
    if (is.null(x)) { skipped <- skipped + 1L; next }
    row <- tryCatch(data.frame(
      n=sc(x$n), delta=sc(x$delta), delta_index=sc(x$delta_index),
      cp_raw=sc(x$cp_raw), cp_decision=sc(x$cp_decision),
      plug_lo_mean=sc(x$plug_lo_mean), plug_hi_mean=sc(x$plug_hi_mean),
      plug_width_mean=sc(x$plug_width_mean),
      ci_lo_mean=sc(x$ci_lo_mean), ci_hi_mean=sc(x$ci_hi_mean),
      ci_width_mean=sc(x$ci_width_mean),
      n_valid=sc(x$n_valid_reps, 0L), n_errors=sc(x$n_errors, 0L)
    ), error = function(e) NULL)
    if (is.null(row)) { skipped <- skipped + 1L; next }
    all_rows[[length(all_rows) + 1L]] <- row
  }
}
if (skipped > 0) cat(sprintf("Skipped %d malformed rds files.\n", skipped))

if (length(all_rows) > 0) {
  summary_df <- do.call(rbind, all_rows)
  utils::write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)
  cat(sprintf("Wrote %d rows to %s\n\n", nrow(summary_df), SUMMARY_CSV))
  
  cat(strrep("=", 78), "\n", sep = "")
  cat("cp_raw pivoted (delta x n)\n", strrep("=", 78), "\n", sep = "")
  piv <- reshape(summary_df[, c("n","delta","cp_raw")],
                 idvar = "delta", timevar = "n", direction = "wide")
  names(piv) <- gsub("cp_raw\\.", "n=", names(piv))
  piv[, -1] <- round(piv[, -1], 3)
  print(piv, row.names = FALSE)
  
  err_total <- sum(summary_df$n_errors, na.rm = TRUE)
  if (err_total > 0)
    cat(sprintf("\n*** %d errors across all cells -- check live log ***\n", err_total))
}
