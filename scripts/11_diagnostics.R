# ==============================================================================
# 11_diagnostics.R  -  scale of the MCUB sets, breakdown distances,
#                      the two Di Tella LATEs, and the ratio labels
#                      (tab:locality, sec:empirical, app:repro)
# ==============================================================================
#
# PURPOSE
#   Quantities that the paper reports next to the validity verdicts. None
#   of them reruns an MCUB set; the sets come from the outputs of 02 and 03.
#
#   1. Scale of the reported sets. The CR test runs on a prepared frame in
#      which the residualised instrument z0 is replaced by its fitted value
#      z_hat from a regression on the residualised (x, y). The projection
#      keeps Cov(z, x), Cov(z, y) and hence Cov(z, u) for every u = y - b x,
#      but sd(z_hat) = m sd(z0), with m the multiple correlation of z0 with
#      (x, y). So rho(z_hat, u) = rho(z0, u) / m: the sets of 02 and 03 are on
#      the z_hat scale, and multiplying their endpoints by m puts them on the
#      scale of rho_ZU. The zero of g(rho, .), r*, and every verdict (does
#      the set contain zero?) are the same on both scales. The script checks
#      g(rho0, r) = m g(rho_hat, r) on a grid of r.
#   2. Breakdown distance: how far r* lies from the domain D (zero if it lies
#      inside). For an instrument that fails the gate, this is the size of
#      the heterogeneity term in r* (prop:implied(c)) that would make it
#      valid. For one that passes, the margin is the distance from r* to the
#      nearer edge of D.
#   3. Share of the B_RSTAR bootstrap values of r* that lie in the closed
#      domain, and the MCUB margin in r units: the MCUB set adds about e to
#      each end of the plug-in set, which corresponds to e / |dg/dr| at r*.
#   4. For the MCUB-share rule (app:mc, rem:mcub_margin): its frontier from
#      each instrument's own sigma_hat, and the Gaussian prediction of its
#      share. The gate's own frontier is in
#      12_rstar_gate.R.
#   5. Bootstrap intervals for the IV estimates of the instruments that pass,
#      and for the difference of the two Di Tella IV estimates
#      (the same resamples for both).
#   6. DML for the instruments that pass: standard error, and the spread
#      of R across N_SPLITS sample splits.
#   7. Three-way ratio labels from the bootstrap shares of 09: A if at most
#      LABEL_LO of the resamples land in Scenario B, B if at least LABEL_HI
#      do, and "inconclusive" otherwise; "uninformative" when the intervals
#      of both estimates cover zero.
#
# INPUTS
#   scripts/phase1_results.csv, scripts/phase2_results.csv  (02, 03)
#   results/implied_endogeneity.csv                         (09)
#   the data of dial_specs.R
#
# DML
#   DoubleML when it is installed, otherwise foundation/dml_fallback.R (the
#   same estimator); the output records which one ran.
#
# OUTPUT
#   results/diagnostics.csv            one row per instrument
#   results/ditella_late_diff.csv      the two Di Tella IVs and their difference
#   results/dml_split_variation.csv    R for each sample split
#   LaTeX snippets on the console
#
# RUNTIME  a few minutes for 1-5; the DML fits of 6 take about 15 minutes
#   on 2 cores (set N_SPLITS <- 0L to skip them).
#
# RUN (after 09; from any folder)
#   Rscript 11_diagnostics.R
# or through run_all.R, which checks the setup first: Rscript run_all.R 09 11
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(haven); library(data.table)
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

# MovieLens and NHANES files are found by dial_specs.R (see its header).
B_RSTAR    <- 2000L      # as in 09, with the same seeds
B_IV       <- 2000L      # bootstrap resamples for the IV intervals
N_SPLITS   <- as.integer(Sys.getenv("DIAL_N_SPLITS", "10"))  # DML splits for the spread of R
N_FOLDS    <- 5L
N_TREES    <- 500L
TAU_R      <- 1.5
LABEL_LO   <- 0.10       # A if the B share is at most this
LABEL_HI   <- 0.90       # B if the B share is at least this
# Frontier settings (tab:frontier)
TAU_V <- 0.95; B_RULE <- 500L; ALPHA0 <- 0.05; PI_STAR <- 0.05

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
source(file.path(FOUNDATION_DIR, "dml_fallback.R"))
.f03 <- file.path(.script_dir, "03_phase2_empirical.R")
if (file.exists(.f03)) {
  .wd <- setwd(.script_dir)
  tryCatch(invisible(capture.output(source(.f03))),
           error = function(e) cat("NOTE: could not source 03 (MovieLens reader):",
                                   conditionMessage(e), "\n"))
  setwd(.wd)
}
source(file.path(.script_dir, "dial_specs.R"))
dial_data_report()

DML_ENGINE <- dial_dml_engine_report()     # stops if DIAL_REQUIRE_DOUBLEML=1
                                           # and DoubleML does not load
dml_fit <- function(d, Y, X, controls, n_rep, seed) {
  set.seed(seed)
  if (DML_ENGINE == "DoubleML") {
    if (requireNamespace("lgr", quietly = TRUE))
      lgr::get_logger("mlr3")$set_threshold("warn")
    dd <- DoubleML::DoubleMLData$new(
      data = data.table::as.data.table(d[, c(Y, X, controls)]),
      y_col = Y, d_cols = X, x_cols = controls)
    lrn <- function() mlr3::lrn("regr.ranger", num.trees = N_TREES, min.node.size = 5)
    dm <- DoubleML::DoubleMLPLR$new(dd, ml_l = lrn(), ml_m = lrn(),
                                    n_folds = N_FOLDS, n_rep = n_rep)
    dm$fit()
    return(list(coef = as.numeric(dm$coef), se = as.numeric(dm$se)))
  }
  f <- dial_dml_plr_rf(d, Y, X, controls, n_rep = n_rep, n_folds = N_FOLDS,
                       num.trees = N_TREES)
  list(coef = f$coef, se = f$se)
}
options(dial.ranger_threads = max(1L, parallel::detectCores() - 1L))

# ==============================================================================
# HELPERS
# ==============================================================================

g_map <- function(rho, r)            # eq:cr_mapping, rho = (xy, xz, yz)
  rho[2] * r - (rho[1] * rho[2] - rho[3]) * sqrt((1 - r^2) / (1 - rho[1]^2))
g_slope <- function(rho, r)          # dg/dr
  rho[2] + (rho[1] * rho[2] - rho[3]) * r / sqrt((1 - r^2) * (1 - rho[1]^2))
implied_rdu <- function(rho) {
  if (any(!is.finite(rho))) return(NA_real_)
  k <- (rho[1] * rho[2] - rho[3]) / (rho[2] * sqrt(1 - rho[1]^2))
  k / sqrt(1 + k^2)
}

# Prepared correlations with z_hat (as dial_prep_rho), the raw ones with z0,
# and m = sd(z_hat) / sd(z0), for the resample with row counts w.
prep_both <- function(eng, w = rep(1, eng$n)) {
  G   <- crossprod(eng$A, eng$A * w)
  iw  <- 3L + seq_len(eng$k)
  Gvw <- G[1:3, iw, drop = FALSE]
  qw  <- qr(G[iw, iw, drop = FALSE], tol = 1e-10)
  cf  <- qr.coef(qw, t(Gvw)); cf[is.na(cf)] <- 0
  R   <- G[1:3, 1:3] - Gvw %*% cf
  syy <- R[1, 1]; sxx <- R[2, 2]; szz <- R[3, 3]
  sxy <- R[1, 2]; syz <- R[1, 3]; sxz <- R[2, 3]
  det <- sxx * syy - sxy^2
  q   <- (syy * sxz^2 - 2 * sxy * sxz * syz + sxx * syz^2) / det
  list(hat = c(sxy / sqrt(sxx * syy), sxz / sqrt(sxx * q), syz / sqrt(syy * q)),
       raw = c(sxy / sqrt(sxx * syy), sxz / sqrt(sxx * szz), syz / sqrt(syy * szz)),
       m   = sqrt(q / szz),
       # IV slope on the standardised columns of the engine
       b_iv_std = syz / sxz)
}

parse_set <- function(txt) {
  v <- suppressWarnings(as.numeric(regmatches(txt, gregexpr("-?[0-9]+\\.?[0-9]*", txt))[[1]]))
  if (length(v) >= 2) v[1:2] else c(NA_real_, NA_real_)
}
dist_to <- function(r, D) if (r >= D[1] && r <= D[2]) 0 else min(abs(r - D))
margin_in <- function(r, D)          # > 0 inside D, < 0 outside
  if (r >= D[1] && r <= D[2]) min(r - D[1], D[2] - r) else -min(abs(r - D))
kappa_star <- function(n, sigma)     # eq:frontier_approx
  sigma / sqrt(n) * (stats::qnorm(1 - ALPHA0 / 2) -
    stats::qnorm(TAU_V - stats::qnorm(1 - PI_STAR) * sqrt(TAU_V * (1 - TAU_V) / B_RULE)))

# MCUB sets of the canonical domain from 02 and 03
p1 <- utils::read.csv(file.path(.script_dir, "phase1_results.csv"), stringsAsFactors = FALSE)
p2 <- utils::read.csv(file.path(.script_dir, "phase2_results.csv"), stringsAsFactors = FALSE)
canon_sets <- function(s) {
  if (s$phase == 1L) {
    r <- p1[p1$instrument == s$instrument & p1$domain_label == "canonical", ]
    mc <- r$CI_MCUB
  } else {
    r <- p2[p2$Instrument == s$instrument & p2$domain_label == "canonical", ]
    mc <- r$CI_Bei
  }
  if (!nrow(r)) return(NULL)
  list(plug = parse_set(r$plug_in[1]), mcub = parse_set(mc[1]),
       p_hat = as.numeric(r$p_hat_B[1]))
}
# 11 uses the bootstrap shares of 09, so 09 must run first, with the same
# DML engine.
.ie_csv <- file.path(RESULTS_DIR, "implied_endogeneity.csv")
if (!file.exists(.ie_csv))
  stop("results/implied_endogeneity.csv not found: run 09_implied_endogeneity.R first.",
       call. = FALSE)
ie <- utils::read.csv(.ie_csv, stringsAsFactors = FALSE)
.ie_eng <- if ("dml_engine" %in% names(ie)) unique(stats::na.omit(ie$dml_engine[!is.na(ie$R)])) else character(0)
cat(sprintf("09 output: %s (written %s; DML engine %s)\n", .ie_csv,
            format(file.mtime(.ie_csv), "%Y-%m-%d %H:%M"),
            if (length(.ie_eng)) paste(.ie_eng, collapse = " + ") else "not recorded"))
if (!identical(.ie_eng, DML_ENGINE)) {
  msg <- sprintf("09 used DML engine '%s' but this run uses '%s'. Rerun 09 with the same R.",
                 if (length(.ie_eng)) paste(.ie_eng, collapse = " + ") else "not recorded", DML_ENGINE)
  if (dial_require_doubleml()) stop(msg, call. = FALSE)
  cat("WARNING:", msg, "\n")
}
.need <- vapply(Filter(function(s) isTRUE(s$valid), SPECS), `[[`, "", "instrument")
.miss <- setdiff(.need, ie$instrument[!is.na(ie$R)])
if (length(.miss))
  cat("WARNING: no R from 09 for", paste(.miss, collapse = ", "),
      "(their labels will be NA; check 09's data report).\n")

# ==============================================================================
# 1-4: scale, breakdown, r* shares, frontier (all 13 instruments)
# ==============================================================================

cat("================================================================================\n")
cat("Diagnostics (11_diagnostics.R)\n")
cat("================================================================================\n\n")

rows <- list(); engines <- list(); datas <- list()
for (s in SPECS) {
  key <- s$instrument
  cat(sprintf("--- %s %s (D = [%.2f, %.2f]) ---\n", s$study, key, s$D[1], s$D[2]))
  cs <- canon_sets(s)
  d_all <- dial_load_spec_data(s)
  out <- data.frame(study = s$study, instrument = key, D_lo = s$D[1], D_hi = s$D[2],
                    valid = s$valid, n = NA_integer_, m = NA_real_, g_check = NA_real_,
                    r_star = NA_real_, r_star_raw = NA_real_,
                    breakdown = NA_real_, margin = NA_real_, share_r_in_D = NA_real_,
                    plug_lo = cs$plug[1], plug_hi = cs$plug[2],
                    mcub_lo = cs$mcub[1], mcub_hi = cs$mcub[2], p_hat_B = cs$p_hat,
                    mcub_lo_resc = NA_real_, mcub_hi_resc = NA_real_,
                    plug_lo_resc = NA_real_, plug_hi_resc = NA_real_,
                    d_signed = NA_real_, d_signed_resc = NA_real_,
                    e_mcub = NA_real_, slope_rstar = NA_real_, margin_r_units = NA_real_,
                    sigma_hat = NA_real_, kappa = NA_real_, kappa_star = NA_real_,
                    share_gauss = NA_real_, stringsAsFactors = FALSE)
  out$d_signed <- min(-cs$mcub[1], cs$mcub[2])
  # MCUB margin: how far the MCUB set extends beyond the plug-in set
  e_ends <- c(cs$plug[1] - cs$mcub[1], cs$mcub[2] - cs$plug[2])
  out$e_mcub <- mean(e_ends)
  # sigma_hat at the endpoint nearest zero (app:frontier)
  near <- which.min(abs(cs$plug))
  out$kappa <- if (cs$plug[1] <= 0 && cs$plug[2] >= 0) 0 else min(abs(cs$plug))
  if (!is.null(d_all)) {
    spec <- dial_raw_spec(d_all, s$Y, s$X, s$Z, s$controls)
    eng  <- dial_prep_engine(spec)
    engines[[key]] <- eng; datas[[key]] <- spec$data
    pb   <- prep_both(eng)
    stopifnot(max(abs(pb$hat - dial_prep_rho(eng))) < 1e-10)
    out$n <- eng$n
    out$m <- pb$m
    rg <- seq(-0.95, 0.95, by = 0.05)
    out$g_check <- max(abs(g_map(pb$raw, rg) - pb$m * g_map(pb$hat, rg)))
    out$r_star     <- implied_rdu(pb$hat)
    out$r_star_raw <- implied_rdu(pb$raw)
    out$breakdown  <- dist_to(out$r_star, s$D)
    out$margin     <- margin_in(out$r_star, s$D)
    set.seed(s$seed)                  # the seed of 09, so the draws match
    rs <- vapply(seq_len(B_RSTAR), function(b) {
      idx <- sample.int(eng$n, eng$n, replace = TRUE)
      implied_rdu(dial_prep_rho(eng, tabulate(idx, nbins = eng$n)))
    }, numeric(1))
    out$share_r_in_D <- mean(rs >= s$D[1] & rs <= s$D[2], na.rm = TRUE)
    out$slope_rstar  <- abs(g_slope(pb$hat, out$r_star))
    out$margin_r_units <- out$e_mcub / out$slope_rstar
    clip <- function(v) pmax(-1, pmin(1, v))
    out$mcub_lo_resc <- clip(pb$m * cs$mcub[1]); out$mcub_hi_resc <- clip(pb$m * cs$mcub[2])
    out$plug_lo_resc <- clip(pb$m * cs$plug[1]); out$plug_hi_resc <- clip(pb$m * cs$plug[2])
    out$d_signed_resc <- min(-out$mcub_lo_resc, out$mcub_hi_resc)
  } else {
    out$n <- if (s$phase == 1L) NA_integer_ else
      as.integer(p2$n[p2$Instrument == key & p2$domain_label == "canonical"][1])
  }
  if (!s$valid) {
    nn <- if (is.na(out$n)) ie$n[ie$instrument == key][1] else out$n
    out$sigma_hat <- sqrt(nn) * abs(cs$plug[near] - cs$mcub[near]) / stats::qnorm(0.975)
    out$kappa_star <- kappa_star(nn, out$sigma_hat)
    out$share_gauss <- stats::pnorm(stats::qnorm(0.975) - sqrt(nn) * out$kappa / out$sigma_hat)
  }
  cat(sprintf("    n = %s  m = %.3f  check %.1e  r* = %+.3f (raw %+.3f)  breakdown %.3f  margin %+.3f  share r* in D %.3f\n",
              format(out$n), out$m, out$g_check, out$r_star, out$r_star_raw,
              out$breakdown, out$margin, out$share_r_in_D))
  cat(sprintf("    MCUB [%+.3f, %+.3f] -> rescaled [%+.3f, %+.3f]; d %+.3f -> %+.3f; e = %.3f, |g'(r*)| = %.3f, margin in r units %.3f\n",
              out$mcub_lo, out$mcub_hi, out$mcub_lo_resc, out$mcub_hi_resc,
              out$d_signed, out$d_signed_resc, out$e_mcub, out$slope_rstar,
              out$margin_r_units))
  if (!s$valid)
    cat(sprintf("    sigma_hat = %.2f  kappa = %.3f  kappa*_N = %.3f  kappa/kappa* = %.1f  Gaussian share = %.3f\n",
                out$sigma_hat, out$kappa, out$kappa_star, out$kappa / out$kappa_star,
                out$share_gauss))
  rows[[key]] <- out
}
res <- do.call(rbind, rows); rownames(res) <- NULL

cat(sprintf("\nMedian sigma_hat of the %d instruments that fail a gate (z_hat scale): %.2f (range %.2f to %.2f)\n",
            sum(!is.na(res$sigma_hat)), stats::median(res$sigma_hat, na.rm = TRUE),
            min(res$sigma_hat, na.rm = TRUE), max(res$sigma_hat, na.rm = TRUE)))
cat(sprintf("Median sigma_hat on the rho_ZU scale (m sigma_hat), where m is known: %.2f\n",
            stats::median(res$m * res$sigma_hat, na.rm = TRUE)))

# ==============================================================================
# 5: IV intervals; the two Di Tella instruments
# ==============================================================================

iv_boot <- function(key, seed) {
  eng <- engines[[key]]; if (is.null(eng)) return(NULL)
  d <- datas[[key]]; s <- Filter(function(z) z$instrument == key, SPECS)[[1]]
  sy <- stats::sd(d[[s$Y]]); sx <- stats::sd(d[[s$X]])
  b0 <- prep_both(eng)$b_iv_std * sy / sx
  set.seed(seed)
  idx_list <- lapply(seq_len(B_IV), function(b) sample.int(eng$n, eng$n, replace = TRUE))
  bb <- vapply(idx_list, function(idx)
    prep_both(eng, tabulate(idx, nbins = eng$n))$b_iv_std * sy / sx, numeric(1))
  list(b = b0, draws = bb, idx = idx_list)
}

cat("\n--- IV bootstrap intervals (valid instruments) ---\n")
iv_ci <- list()
for (s in Filter(function(z) z$valid, SPECS)) {
  ib <- iv_boot(s$instrument, s$seed + 7L)
  if (is.null(ib)) { cat(sprintf("    %s: data not available\n", s$instrument)); next }
  ci <- stats::quantile(ib$draws, c(0.025, 0.975), names = FALSE)
  iv_ci[[s$instrument]] <- c(b = ib$b, lo = ci[1], hi = ci[2])
  cat(sprintf("    %-20s b_IV = %+.4f  95%% [%+.4f, %+.4f]\n", s$instrument, ib$b, ci[1], ci[2]))
}

# The two Di Tella instruments on the same resamples of the same rows
cat("\n--- Di Tella: difference of the two IV estimates ---\n")
dt <- NULL
s1 <- SPECS[[1]]; s2 <- SPECS[[2]]
d_dt <- dial_load_spec_data(s1)
if (!is.null(d_dt)) {
  vars <- unique(c(s1$Y, s1$X, s1$Z, s2$Z, s1$controls))
  dd <- d_dt[stats::complete.cases(d_dt[, vars]), vars]
  e1 <- dial_prep_engine(list(data = dd, Y = s1$Y, X = s1$X, Z = s1$Z, controls = s1$controls))
  e2 <- dial_prep_engine(list(data = dd, Y = s2$Y, X = s2$X, Z = s2$Z, controls = s2$controls))
  sy <- stats::sd(dd[[s1$Y]]); sx <- stats::sd(dd[[s1$X]])
  biv <- function(e, w) prep_both(e, w)$b_iv_std * sy / sx
  b1 <- biv(e1, rep(1, e1$n)); b2 <- biv(e2, rep(1, e2$n))
  set.seed(2026L)
  dr <- t(vapply(seq_len(B_IV), function(b) {
    w <- tabulate(sample.int(e1$n, e1$n, replace = TRUE), nbins = e1$n)
    c(biv(e1, w), biv(e2, w))
  }, numeric(2)))
  diff_b <- dr[, 1] - dr[, 2]
  ci_d   <- stats::quantile(diff_b, c(0.025, 0.975), names = FALSE)
  # bootstrap-t style: share of resamples on the other side of zero, two-sided
  p_d <- 2 * min(mean(diff_b >= 0), mean(diff_b <= 0))
  dt <- data.frame(n = e1$n, b_judgeAlreadyUsedEM = b1, b_percJudgeSentToEM = b2,
                   diff = b1 - b2, diff_lo = ci_d[1], diff_hi = ci_d[2],
                   se_diff = stats::sd(diff_b), p_boot = p_d,
                   cor_draws = stats::cor(dr[, 1], dr[, 2]), B = B_IV)
  cat(sprintf("    n = %d  b1 = %+.4f  b2 = %+.4f  diff = %+.4f  95%% [%+.4f, %+.4f]  se %.4f  p = %.3f  cor(b1, b2) = %.2f\n",
              dt$n, b1, b2, dt$diff, ci_d[1], ci_d[2], dt$se_diff, p_d, dt$cor_draws))
  utils::write.csv(dt, file.path(RESULTS_DIR, "ditella_late_diff.csv"), row.names = FALSE)
}

# ==============================================================================
# 6: DML standard errors and the spread of R across sample splits
# ==============================================================================

split_rows <- list()
res$b_IV <- NA_real_; res$iv_lo <- NA_real_; res$iv_hi <- NA_real_
res$b_DML_fit <- NA_real_; res$se_DML <- NA_real_
res$R_split_min <- NA_real_; res$R_split_max <- NA_real_; res$R_split_med <- NA_real_
for (s in Filter(function(z) z$valid, SPECS)) {
  i <- which(res$instrument == s$instrument)
  if (!is.null(iv_ci[[s$instrument]])) {
    res$b_IV[i] <- iv_ci[[s$instrument]]["b"]
    res$iv_lo[i] <- iv_ci[[s$instrument]]["lo"]; res$iv_hi[i] <- iv_ci[[s$instrument]]["hi"]
  }
  d <- datas[[s$instrument]]
  if (is.null(d) || N_SPLITS < 1L) next
  t1 <- Sys.time()
  f0 <- dml_fit(d, s$Y, s$X, s$controls, n_rep = 3L, seed = s$seed)
  res$b_DML_fit[i] <- f0$coef; res$se_DML[i] <- f0$se
  th <- vapply(seq_len(N_SPLITS), function(k)
    dml_fit(d, s$Y, s$X, s$controls, n_rep = 1L, seed = s$seed + 100L * k)$coef, numeric(1))
  Rk <- abs(res$b_IV[i]) / abs(th)
  split_rows[[s$instrument]] <- data.frame(instrument = s$instrument, split = seq_len(N_SPLITS),
                                           b_DML = th, R = Rk, engine = DML_ENGINE)
  res$R_split_min[i] <- min(Rk); res$R_split_max[i] <- max(Rk); res$R_split_med[i] <- stats::median(Rk)
  cat(sprintf("    %-20s DML (3 splits) %+.4f (se %.4f); R over %d splits: %.2f to %.2f (median %.2f); %.1f min\n",
              s$instrument, f0$coef, f0$se, N_SPLITS, min(Rk), max(Rk), stats::median(Rk),
              as.numeric(difftime(Sys.time(), t1, units = "mins"))))
}
if (length(split_rows))
  utils::write.csv(do.call(rbind, split_rows),
                   file.path(RESULTS_DIR, "dml_split_variation.csv"), row.names = FALSE)

# ==============================================================================
# 7: three-way ratio labels
# ==============================================================================

res <- merge(res, ie[, c("instrument", "b_DML", "R", "R_lo", "R_hi", "share_B", "share_sign")],
             by = "instrument", all.x = TRUE, sort = FALSE)
lab <- function(i) {
  r <- res[i, ]
  if (!r$valid) return(NA_character_)
  dml_lo <- r$b_DML_fit - 1.96 * r$se_DML; dml_hi <- r$b_DML_fit + 1.96 * r$se_DML
  if (all(is.finite(c(r$iv_lo, r$iv_hi, dml_lo, dml_hi))) &&
      r$iv_lo <= 0 && r$iv_hi >= 0 && dml_lo <= 0 && dml_hi >= 0)
    return("uninformative")
  if (is.na(r$share_B)) return(NA_character_)
  if (r$share_B >= LABEL_HI) return(if (isTRUE(r$share_sign >= 0.5)) "B (sign)" else "B")
  if (r$share_B <= LABEL_LO) return("A")
  "inconclusive"
}
res$label3 <- vapply(seq_len(nrow(res)), lab, character(1))
res$dml_engine <- DML_ENGINE
ord <- vapply(SPECS, `[[`, "", "instrument")
res <- res[match(ord, res$instrument), ]
utils::write.csv(res, file.path(RESULTS_DIR, "diagnostics.csv"), row.names = FALSE)

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")
print(res[, c("study", "instrument", "n", "m", "r_star", "breakdown", "margin",
              "share_r_in_D", "p_hat_B", "margin_r_units", "mcub_lo_resc", "mcub_hi_resc",
              "share_B", "label3")], digits = 3, row.names = FALSE)

cat("\nLaTeX rows (instrument & D & r* & breakdown & share r* in D & p_hat_B & label):\n\n")
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  cat(sprintf("%s \\texttt{%s} & $[%g,\\,%g]$ & $%+.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & %s \\\\\n",
              r$study, gsub("_", "\\\\_", r$instrument), r$D_lo, r$D_hi, r$r_star,
              r$breakdown, r$share_r_in_D, r$p_hat_B,
              ifelse(is.na(r$label3), "C", r$label3)))
}
cat("\nRescaled MCUB sets (rho_ZU scale):\n\n")
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  cat(sprintf("%-20s m = %5.3f  [%+.3f, %+.3f]  (z_hat scale [%+.3f, %+.3f])\n",
              r$instrument, r$m, r$mcub_lo_resc, r$mcub_hi_resc, r$mcub_lo, r$mcub_hi))
}
cat(sprintf("\nWrote %s\n", file.path(RESULTS_DIR, "diagnostics.csv")))
cat("\nDONE.\n")
