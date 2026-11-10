# ==============================================================================
# 09_implied_endogeneity.R  -  implied endogeneity r* and intervals for R
#                              (tab:locality, sec:space, sec:empirical)
# ==============================================================================
#
# PURPOSE
#   For each of the 13 empirical instruments, under its canonical domain:
#     r*    the zero of the CR map g(rho, r), that is, the residual
#           endogeneity rho_DU at which the data are consistent with a
#           valid instrument. On the prepared data,
#             r* = cor(x, y - b_IV x)   and
#             b_OLS - b_IV = r* sd(y - b_IV x) / sd(x),
#           so for an instrument that passes the validity gate, r* is the
#           residual confounding that would explain the whole IV-OLS gap if
#           effects were homogeneous. Reported with a 95% percentile
#           interval from B_RSTAR bootstrap resamples of the raw rows (the
#           residualise-and-project step is redone on every resample).
#   For the six instruments that pass both gates (tab:empirical_full):
#     b_IV, b_DML, R = |b_IV| / |b_DML|, the sign check, and a 95%
#     percentile interval for R from B_R bootstrap resamples. IV and DML
#     are refitted on each resample. DML folds are assigned to the
#     ORIGINAL rows, so copies of a row never sit in both the training and
#     the test fold.
#
# DML SPECIFICATION (as in the paper)
#   Partially linear model, random forests (ranger, 500 trees, minimum node
#   size 5), 5-fold cross-fitting. Point estimates use 3 repetitions of the
#   sample split; bootstrap refits use 1. Each instrument has its own seed,
#   so the estimates do not depend on the order of the runs. MovieLens
#   uses seed 43, the seed of 03_phase2_empirical.R, which reproduces its
#   estimate. The paper reports these seeded values: R = 3.13, 1.37, 0.62,
#   2.44, 1.17 and 1.84 for judgeAlreadyUsedEM, percJudgeSentToEM,
#   num_genres, highnumber, fatheduc and instru.
#
# DATA
#   Phase 1 .dta files are read from data/ at the repository root through
#   dial_data_path() in foundation/Fin_Empirical5_clean.R. MovieLens needs
#   ML_DIR (the folder that holds u.data and u.item), Papke and Card the
#   wooldridge package, and NHANES the nhanes_mercury.csv written by
#   01_nhanes_build.R.
#
# OUTPUT
#   results/implied_endogeneity.csv, and LaTeX rows for tab:locality
#
# RUNTIME  about 20 minutes on 7 cores (the DML refits); r* alone takes
#   seconds (set B_R <- 0L to skip the intervals for R).
#
# RUN (from any folder; 11 needs this output, so run 09 first)
#   Rscript 09_implied_endogeneity.R
# or through run_all.R, which checks the setup first: Rscript run_all.R 09 11
# Set DIAL_REQUIRE_DOUBLEML=1 to stop instead of using the DML fallback.
# ==============================================================================

rm(list = ls())

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) try(install.packages(new, repos = "https://cloud.r-project.org/"))
}
install_if_missing(c("dplyr", "haven", "Matrix", "MASS", "data.table",
                     "DoubleML", "mlr3", "mlr3learners", "ranger",
                     "future.apply", "wooldridge"))

suppressPackageStartupMessages({
  library(dplyr); library(haven); library(Matrix); library(data.table)
  library(future.apply)
})

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;    filter    <- dplyr::filter
mutate    <- dplyr::mutate;    slice     <- dplyr::slice
recode    <- dplyr::recode;    rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange   <- dplyr::arrange;   count     <- dplyr::count
lag       <- dplyr::lag

# ==============================================================================
# SETTINGS
# ==============================================================================

# MovieLens (ml-100k/) and NHANES (nhanes_mercury.csv) are found by
# dial_specs.R; see its header, or set DIAL_ML_DIR and DIAL_NHANES_CSV.
B_RSTAR    <- 2000L                  # bootstrap resamples for r*
B_R        <- 200L                   # bootstrap resamples for R (DML refits)
N_FOLDS    <- 5L
N_TREES    <- 500L
TAU_R      <- 1.5
# Parallel workers for the bootstrap of R (multisession: works on Windows
# too). Override with the environment variable DIAL_WORKERS.
N_WORKERS  <- as.integer(Sys.getenv("DIAL_WORKERS",
                min(7L, max(1L, parallel::detectCores(logical = FALSE) - 1L))))

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
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

assign("DIAL_SKIP_AUTORUN", TRUE, envir = .GlobalEnv)
invisible(capture.output(source(file.path(FOUNDATION_DIR, "Fin_Empirical5_clean.R"))))
invisible(capture.output(source(file.path(FOUNDATION_DIR, "bootstrap_rule.R"))))
# 03 provides the MovieLens reader (.load_ml100k); sourcing it runs nothing.
.f03 <- file.path(.script_dir, "03_phase2_empirical.R")
if (file.exists(.f03)) {
  .wd <- setwd(.script_dir)
  invisible(capture.output(source(.f03)))
  setwd(.wd)
}
for (fn in c("dial_raw_spec", "dial_prep_engine", "dial_prep_rho"))
  if (!exists(fn, mode = "function")) stop("bootstrap_rule.R did not provide ", fn)

# ==============================================================================
# ESTIMATORS
# ==============================================================================

# Zero of g(rho, r) = rho_xz r - (rho_xy rho_xz - rho_yz) sqrt((1 - r^2) / (1 - rho_xy^2)),
# rho = (rho_xy, rho_xz, rho_yz) of the prepared data; equals cor(x, y - b_IV x).
implied_rdu <- function(rho) {
  if (any(!is.finite(rho))) return(NA_real_)
  k <- (rho[1] * rho[2] - rho[3]) / (rho[2] * sqrt(1 - rho[1]^2))
  k / sqrt(1 + k^2)
}

# IV, OLS and first-stage F with controls (Frisch-Waugh; same numbers as
# 2SLS and the weak-instrument F of ivreg).
iv_ols_F <- function(d, Y, X, Z, controls) {
  rhs <- if (length(controls)) paste(controls, collapse = " + ") else "1"
  res <- function(v) stats::resid(stats::lm(stats::as.formula(paste(v, "~", rhs)), data = d))
  y <- res(Y); x <- res(X); z <- res(Z)
  full <- stats::lm(stats::as.formula(paste(X, "~", Z, "+", rhs)), data = d)
  rest <- stats::lm(stats::as.formula(paste(X, "~", rhs)), data = d)
  c(iv  = sum(z * y) / sum(z * x),
    ols = sum(x * y) / sum(x^2),
    F   = stats::anova(rest, full)$F[2])
}

# DML (partially linear model). With 'folds' (one fold number per row), the
# cross-fitting uses that split; otherwise DoubleML draws n_rep splits.
# Defined only if absent, so a test harness can supply its own. Without
# DoubleML, foundation/dml_fallback.R supplies the same estimator.
source(file.path(FOUNDATION_DIR, "dml_fallback.R"))
DML_ENGINE <- dial_dml_engine_report()     # stops here if DIAL_REQUIRE_DOUBLEML=1
                                           # and DoubleML does not load
if (!exists("dml_plr", mode = "function") && DML_ENGINE != "DoubleML") {
  dml_plr <- function(d, Y, X, controls, folds = NULL, n_rep = 1L,
                      num.trees = N_TREES)
    dial_dml_plr_rf(d, Y, X, controls, folds = folds, n_rep = n_rep,
                    n_folds = N_FOLDS, num.trees = num.trees)$coef
}
if (!exists("dml_plr", mode = "function")) {
  dml_plr <- function(d, Y, X, controls, folds = NULL, n_rep = 1L,
                      num.trees = N_TREES) {
    if (!requireNamespace("mlr3learners", quietly = TRUE))   # registers regr.ranger
      stop("mlr3learners is required for the DML estimate")
    if (requireNamespace("lgr", quietly = TRUE))
      lgr::get_logger("mlr3")$set_threshold("warn")
    dd <- DoubleML::DoubleMLData$new(
      data = data.table::as.data.table(d[, c(Y, X, controls)]),
      y_col = Y, d_cols = X, x_cols = controls)
    lrn_rf <- function() mlr3::lrn("regr.ranger", num.trees = num.trees,
                                   min.node.size = 5)
    if (is.null(folds)) {
      dm <- DoubleML::DoubleMLPLR$new(dd, ml_l = lrn_rf(), ml_m = lrn_rf(),
                                      n_folds = N_FOLDS, n_rep = n_rep)
    } else {
      K <- max(folds)
      smpls <- list(list(
        train_ids = lapply(seq_len(K), function(k) which(folds != k)),
        test_ids  = lapply(seq_len(K), function(k) which(folds == k))))
      dm <- DoubleML::DoubleMLPLR$new(dd, ml_l = lrn_rf(), ml_m = lrn_rf(),
                                      n_folds = K, draw_sample_splitting = FALSE)
      dm$set_sample_splitting(smpls)
    }
    dm$fit()
    as.numeric(dm$coef)
  }
}

# Bootstrap of R: IV and DML refitted on each resample of the rows.
boot_ratio <- function(d, Y, X, Z, controls, B, seed) {
  n <- nrow(d)
  set.seed(seed)
  fold_orig <- sample(rep_len(seq_len(N_FOLDS), n))   # one fold per original row
  idx_list  <- lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE))
  one <- function(b) {
    idx <- idx_list[[b]]
    db  <- d[idx, , drop = FALSE]; rownames(db) <- NULL
    iv  <- iv_ols_F(db, Y, X, Z, controls)[["iv"]]
    set.seed(seed + b)
    dml <- tryCatch(dml_plr(db, Y, X, controls, folds = fold_orig[idx], n_rep = 1L),
                    error = function(e) NA_real_)
    c(iv = iv, dml = dml)
  }
  out <- if (N_WORKERS > 1L) {
    oplan <- future::plan(future::multisession, workers = N_WORKERS)
    on.exit(future::plan(oplan), add = TRUE)
    future.apply::future_lapply(seq_len(B), one, future.seed = TRUE,
      future.packages = if (dial_dml_available())
        c("DoubleML", "mlr3", "mlr3learners", "data.table") else "ranger")
  } else lapply(seq_len(B), one)
  do.call(rbind, out)
}

# ==============================================================================
# DATA: the 13 instruments, shared with 11_diagnostics.R and 12_rstar_gate.R
# ==============================================================================

source(file.path(.script_dir, "dial_specs.R"))
loaders <- LOADERS
dial_data_report()

# ==============================================================================
# MAIN LOOP
# ==============================================================================

cat("================================================================================\n")
cat("Implied endogeneity r* and bootstrap intervals for R (09_implied_endogeneity.R)\n")
cat("================================================================================\n")
cat(sprintf("  B_RSTAR = %d, B_R = %d, workers = %d, DML engine = %s\n\n",
            B_RSTAR, B_R, N_WORKERS, DML_ENGINE))

data_cache <- list()
rows <- list()
t0 <- Sys.time()

# DIAL_ONLY=card,fatheduc runs only the named instruments (study or
# instrument names, case-insensitive) and merges them into the existing CSV.
.only <- tolower(trimws(strsplit(Sys.getenv("DIAL_ONLY", ""), ",")[[1]]))
.only <- .only[nzchar(.only)]
for (s in SPECS) {
  if (length(.only) && !(tolower(s$instrument) %in% .only || tolower(s$study) %in% .only)) next
  cat(sprintf("--- %s %s (D = (%.2f, %.2f)) ---\n", s$study, s$instrument, s$D[1], s$D[2]))
  d_all <- data_cache[[s$data]]
  if (is.null(d_all)) {
    d_all <- tryCatch(loaders[[s$data]](), error = function(e) {
      cat("    DATA NOT AVAILABLE:", conditionMessage(e), "\n"); NULL })
    if (!is.null(d_all)) data_cache[[s$data]] <- d_all
  }
  if (is.null(d_all)) next

  spec <- dial_raw_spec(d_all, s$Y, s$X, s$Z, s$controls)   # complete cases
  d    <- spec$data
  eng  <- dial_prep_engine(spec)
  rho  <- dial_prep_rho(eng)
  r_star <- implied_rdu(rho)

  # Check the identity on the prepared data
  prep  <- dial_prepare_xyz(d, s$Y, s$X, s$Z, s$controls)
  b_ivp <- sum(prep$z * prep$y) / sum(prep$z * prep$x)
  r_chk <- stats::cor(prep$x, prep$y - b_ivp * prep$x)

  set.seed(s$seed)
  rs <- vapply(seq_len(B_RSTAR), function(b) {
    idx <- sample.int(eng$n, eng$n, replace = TRUE)
    implied_rdu(dial_prep_rho(eng, tabulate(idx, nbins = eng$n)))
  }, numeric(1))
  r_ci <- stats::quantile(rs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

  est <- iv_ols_F(d, s$Y, s$X, s$Z, s$controls)
  out <- data.frame(study = s$study, instrument = s$instrument,
                    D_lo = s$D[1], D_hi = s$D[2], n = nrow(d),
                    F_first = unname(est["F"]), b_IV = unname(est["iv"]),
                    b_OLS = unname(est["ols"]),
                    r_star = r_star, r_lo = r_ci[1], r_hi = r_ci[2],
                    r_in_D = r_star >= s$D[1] & r_star <= s$D[2],
                    b_DML = NA_real_, R = NA_real_, R_lo = NA_real_, R_hi = NA_real_,
                    share_B = NA_real_, share_sign = NA_real_, n_boot = NA_integer_,
                    dml_engine = DML_ENGINE, stringsAsFactors = FALSE)
  cat(sprintf("    n = %d  F = %.1f  b_IV = %.4f  b_OLS = %.4f  r* = %+.3f [%+.3f, %+.3f] (check %+.3f)  r* in D: %s\n",
              out$n, out$F_first, out$b_IV, out$b_OLS, r_star, r_ci[1], r_ci[2],
              r_chk, out$r_in_D))

  if (isTRUE(s$valid)) {
    set.seed(s$seed)
    b_dml <- dml_plr(d, s$Y, s$X, s$controls, folds = NULL, n_rep = 3L)
    out$b_DML <- b_dml
    out$R     <- abs(out$b_IV) / abs(b_dml)
    cat(sprintf("    b_DML = %.4f  R = %.2f  same sign: %s\n",
                b_dml, out$R, sign(out$b_IV) == sign(b_dml)))
    if (B_R > 0L) {
      t1 <- Sys.time()
      bt <- boot_ratio(d, s$Y, s$X, s$Z, s$controls, B_R, seed = s$seed * 100L)
      ok <- stats::complete.cases(bt)
      Rb <- abs(bt[ok, "iv"]) / abs(bt[ok, "dml"])
      same_b <- sign(bt[ok, "iv"]) == sign(bt[ok, "dml"])
      out$R_lo <- stats::quantile(Rb, 0.025, names = FALSE)
      out$R_hi <- stats::quantile(Rb, 0.975, names = FALSE)
      out$share_B    <- mean(!same_b | abs(log(Rb)) >= log(TAU_R))
      out$share_sign <- mean(!same_b)
      out$n_boot     <- sum(ok)
      cat(sprintf("    R 95%% interval [%.2f, %.2f]; share of resamples in B = %.3f (sign flips %.3f); %d/%d ok, %.1f min\n",
                  out$R_lo, out$R_hi, out$share_B, out$share_sign, sum(ok), B_R,
                  as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    }
  }
  rows[[length(rows) + 1L]] <- out
}

res <- do.call(rbind, rows)
.csv <- file.path(RESULTS_DIR, "implied_endogeneity.csv")
if (length(.only) && file.exists(.csv)) {          # merge into the earlier run
  old <- utils::read.csv(.csv, stringsAsFactors = FALSE)
  old <- old[!(old$instrument %in% res$instrument), , drop = FALSE]
  if (!"dml_engine" %in% names(old)) old$dml_engine <- rep(NA_character_, nrow(old))
  res <- rbind(old[, names(res)], res)
  ord <- vapply(SPECS, `[[`, "", "instrument")
  res <- res[order(match(res$instrument, ord)), ]
}
utils::write.csv(res, .csv, row.names = FALSE)

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")
print(res[, c("study", "instrument", "n", "b_IV", "b_DML", "R", "R_lo", "R_hi",
              "share_B", "r_star", "r_lo", "r_hi", "r_in_D")], digits = 3, row.names = FALSE)

cat("\nLaTeX rows (instrument & D & r* & 95% interval & R & 95% interval):\n\n")
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  Rtxt <- if (is.na(r$R)) "---" else sprintf("$%.2f$", r$R)
  Rci  <- if (is.na(r$R_lo)) "---" else sprintf("$[%.2f,\\,%.2f]$", r$R_lo, r$R_hi)
  cat(sprintf("%s \\texttt{%s} & $(%g,\\,%g)$ & $%+.3f$ & $[%+.3f,\\,%+.3f]$ & %s & %s \\\\\n",
              r$study, gsub("_", "\\\\_", r$instrument), r$D_lo, r$D_hi,
              r$r_star, r$r_lo, r$r_hi, Rtxt, Rci))
}
cat(sprintf("\nWrote %s\n", file.path(RESULTS_DIR, "implied_endogeneity.csv")))
cat(sprintf("Total elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("\nDONE.\n")
