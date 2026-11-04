# ==============================================================================
# bootstrap_rule.R  --  DIAL bootstrap validity rule (Proposition 2, tau_V)
# ------------------------------------------------------------------------------
# WHAT IT DOES
#   For an instrument with analysis data (outcome Y, treatment X, instrument
#   Z, controls), the rule draws B_rule bootstrap resamples of the rows. On
#   each resample b it redoes the residualise-and-project step, computes a
#   new MCUB set C^b and records
#       I_b = 1{ 0 in C^b }        (zero within +/- zero_tol of the set)
#   The bootstrap validity fraction is p_hat_B = mean(I_b). DIAL calls the
#   instrument Valid when p_hat_B >= tau_V (default 0.95), else Invalid.
#   (The code calls this threshold tau2: option dial.tau2, column Valid_tau2.)
#
#   Each resample repeats the whole MCUB computation: its own correlations,
#   its own bootstrap covariance of the correlations, its own Jacobian rows
#   and its own CIhybrid() call. By default each resample uses the settings
#   of scripts 05 and 07: B = 500 draws per bisection step, B_large = 499
#   coverage draws and a covariance from B_cov = 200 bootstrap draws. The
#   full-sample set from check_compatibility_simple() keeps its own
#   settings (B = 500, B_large = 5,000, B_cov = 800).
#
#   THE PREPARATION STEP IS REDONE ON EVERY RESAMPLE. The pipeline
#   residualises Y, X and Z on the controls and replaces the residualised
#   instrument by its fitted values on the residualised (X, Y). Those
#   fitted values depend on estimated coefficients, so a bootstrap that
#   resamples the prepared (x, y, z) rows holds the coefficients fixed and
#   understates the sampling variability of the correlations. Here every
#   bootstrap resample -- each of the B_rule resamples, every draw of each
#   resample's covariance bootstrap, and every draw of the covariance
#   bootstrap behind the full-sample MCUB set -- redraws the rows of the
#   analysis data and redoes the residualise-and-project step.
#   dial_prep_engine() makes this cheap: it turns the step into one
#   weighted cross-product per draw (see "Residualise-and-project").
#
#   Cost: 1-5 s per resample on one core, depending on the data, so
#   roughly 10-40 core-minutes per instrument at B_rule = 500. The
#   resamples use faster copies of the MCUB routines that return
#   bit-identical sets (see "Faster copies" below).
#
# USAGE
#   Source this file AFTER foundation/Fin_Empirical5_clean.R:
#     source("../foundation/bootstrap_rule.R")
#     dial_boot_plan()      # optional: parallel workers via future.apply
#   Sourcing wraps check_compatibility_simple(), so every call (including
#   the calls inside run_*_analysis_fixed()) uses the raw-data covariance
#   for its full-sample MCUB set and returns five extra columns:
#     cov_unit    resampling unit of the full-sample covariance ("raw"/"rows")
#     p_hat_B     bootstrap validity fraction
#     n_boot_ok   resamples whose MCUB set could be computed
#     boot_unit   resampling unit of the rule ("raw" or "rows"; see below)
#     Valid_tau2  "Valid" if p_hat_B >= tau_V, else "Invalid"
#   The rule uses seed 5000 + i, where i is the instrument index passed to
#   check_compatibility_simple(). The full-sample covariance keeps the seed
#   and size of the original call (seed 1000 + i, 800 draws).
#
# RESAMPLING UNIT  (option dial.boot_unit)
#   "raw" (default): resample the rows of the analysis data (outcome,
#     treatment, instrument, controls) and redo the residualise-and-project
#     step on each resample, in the rule and in every covariance bootstrap
#     (see THE PREPARATION STEP above). The wrapper needs the raw data:
#     attach it with dial_attach_raw(df, data, Y, X, Z, controls), or list
#     candidate specifications from dial_raw_spec() in
#     options(dial.boot_raw = list(...)); the wrapper uses the one whose
#     preparation reproduces df. Without raw data, the call falls back to
#     "rows" and says so.
#   "rows": the earlier behaviour, kept for comparison. Resample the rows
#     of the prepared frame (x, y, z), as estimate_cov_corr_boot() does, so
#     the coefficients of the projection stay fixed across resamples.
#
# OPTIONS (set before or after sourcing; read at call time)
#   dial.boot_rule       TRUE  -> run the rule; FALSE -> plain CR test only
#   dial.boot_B          B_rule, number of resamples              (500)
#   dial.tau2            threshold tau_V                           (0.95)
#   dial.zero_tol        tolerance on the set endpoints            (0)
#                        0 = strict verdict (empirical tables);
#                        scripts 04 and 06 set 0.05 (ZERO_TOL)
#   dial.boot_mcub       list(B = 500, Blarge = 499, B_cov = 200)
#                        MCUB settings for each resample
#   dial.boot_cache_dir  cache folder              (<repo>/results/boot_rule)
#                        NA switches caching off
#   dial.boot_workers    parallel workers; set by dial_boot_plan()     (1)
#   dial.boot_unit       "raw" or "rows" (see RESAMPLING UNIT)      ("raw")
#   dial.boot_raw        candidate raw-data specifications (see above)
#   dial.emp5_path       Fin_Empirical5_clean.R, sourced by parallel workers
#
# CACHE
#   Each call saves the B_rule x 2 matrix of resample set endpoints, so an
#   interrupted run resumes where it stopped, and p_hat_B can be recomputed
#   for another zero_tol without new resamples. The file name records n,
#   the domain, alpha, B_rule, the seed, the MCUB settings, the sample
#   correlations and the resampling unit ("v2_..._uraw" for the raw unit).
#   Resample b always uses seed (5000 + i) * 10000 + b, so a run with more
#   resamples reuses a cached run with fewer (it computes only the new
#   ones), and a run with fewer takes the first rows of a larger one. Delete
#   the folder after editing the MCUB code.
# ==============================================================================

# This file attaches no packages; it calls stats::, utils::, MASS::, future::
# and future.apply:: explicitly.

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;    filter    <- dplyr::filter
mutate    <- dplyr::mutate;    slice     <- dplyr::slice
recode    <- dplyr::recode;    rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange   <- dplyr::arrange;   count     <- dplyr::count
lag       <- dplyr::lag

# ---- Paths and defaults -----------------------------------------------------
.dial_br_file <- local({
  f <- NA_character_
  for (fr in rev(sys.frames())) {
    if (exists("ofile", envir = fr, inherits = FALSE)) {
      cand <- get("ofile", envir = fr, inherits = FALSE)
      if (is.character(cand) && length(cand) == 1L && file.exists(cand)) {
        f <- normalizePath(cand); break
      }
    }
  }
  f
})
.dial_br_dir <- if (is.na(.dial_br_file)) getwd() else dirname(.dial_br_file)
if (is.null(getOption("dial.boot_rule_path")))
  options(dial.boot_rule_path = if (is.na(.dial_br_file))
    normalizePath(file.path(.dial_br_dir, "bootstrap_rule.R"), mustWork = FALSE)
    else .dial_br_file)

if (is.null(getOption("dial.emp5_path"))) {
  .p <- file.path(.dial_br_dir, "Fin_Empirical5_clean.R")
  if (!file.exists(.p) && exists("FOUNDATION_DIR", envir = globalenv()))
    .p <- file.path(get("FOUNDATION_DIR", envir = globalenv()), "Fin_Empirical5_clean.R")
  options(dial.emp5_path = normalizePath(.p, mustWork = FALSE))
  rm(.p)
}
if (is.null(getOption("dial.boot_cache_dir")))
  options(dial.boot_cache_dir = file.path(dirname(.dial_br_dir), "results", "boot_rule"))
if (is.null(getOption("dial.boot_rule")))  options(dial.boot_rule = TRUE)
if (is.null(getOption("dial.boot_B")))     options(dial.boot_B = 500L)
if (is.null(getOption("dial.tau2")))       options(dial.tau2 = 0.95)
if (is.null(getOption("dial.zero_tol")))   options(dial.zero_tol = 0)
if (is.null(getOption("dial.boot_mcub")))
  options(dial.boot_mcub = list(B = 500L, Blarge = 499L, B_cov = 200L))
if (is.null(getOption("dial.boot_unit")))  options(dial.boot_unit = "raw")

# ---- Raw-data specifications (for dial.boot_unit = "raw") ---------------------
# The residualise-and-project step of the pipeline: residualise Y, X and Z on
# the controls (on an intercept if there are none), then replace the
# residualised instrument by its fitted values on the residualised (X, Y).
dial_prepare_xyz <- function(data, Y, X, Z, controls = character(0)) {
  rhs <- if (length(controls)) paste(controls, collapse = "+") else "1"
  fml <- function(lhs) stats::as.formula(paste(lhs, "~", rhs))
  y  <- stats::resid(stats::lm(fml(Y), data = data))
  x  <- stats::resid(stats::lm(fml(X), data = data))
  z0 <- stats::resid(stats::lm(fml(Z), data = data))
  z  <- stats::predict(stats::lm(z0 ~ x + y))
  data.frame(x = as.numeric(x), y = as.numeric(y), z = as.numeric(z))
}
dial_raw_spec <- function(data, Y, X, Z, controls = character(0)) {
  vars <- c(Y, X, Z, controls)
  d <- as.data.frame(data)[, vars, drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  rownames(d) <- NULL
  list(data = d, Y = Y, X = X, Z = Z, controls = controls)
}
dial_attach_raw <- function(df, data, Y, X, Z, controls = character(0)) {
  attr(df, "dial_raw") <- dial_raw_spec(data, Y, X, Z, controls)
  df
}
# Pick the specification whose preparation reproduces the prepared frame df.
.dial_match_raw <- function(df, specs) {
  if (is.null(specs)) return(NULL)
  if (!is.null(specs$data)) specs <- list(specs)
  r0 <- c(stats::cor(df$x, df$y), stats::cor(df$x, df$z), stats::cor(df$y, df$z))
  for (s in specs) {
    p <- tryCatch(dial_prepare_xyz(s$data, s$Y, s$X, s$Z, s$controls),
                  error = function(e) NULL)
    if (is.null(p) || nrow(p) != nrow(df)) next
    r1 <- c(stats::cor(p$x, p$y), stats::cor(p$x, p$z), stats::cor(p$y, p$z))
    if (max(abs(r1 - r0)) < 1e-8) return(s)
  }
  NULL
}

# ---- Residualise-and-project for resampled rows -------------------------------
# A bootstrap resample is a vector of row counts w (w_i = number of times row
# i is drawn). Running lm() on the resampled rows is the same as weighted
# least squares with weights w, so the prepared correlations of a resample
# follow from the weighted cross-products of V = (Y, X, Z) and the control
# matrix W (with an intercept):
#   G = A' diag(w) A,  A = [V, W],
#   R = G_VV - G_VW G_WW^- G_WV   (cross-products of the residualised y, x, z0).
# The projection of z0 on (x, y) then gives, with s = entries of R,
#   cor(x, y) = s_xy / sqrt(s_xx s_yy),
#   cor(x, z) = s_xz0 / sqrt(s_xx q),  cor(y, z) = s_yz0 / sqrt(s_yy q),
#   q = c' S^-1 c,  c = (s_xz0, s_yz0),  S = [[s_xx, s_xy], [s_xy, s_yy]],
# because the fitted values keep the covariances of z0 with x and y, and q
# is their cross-product with themselves. The residuals have weighted mean
# zero (W has an intercept), so the fitted values have no intercept term.
# One draw therefore costs one cross-product instead of four lm() fits. The
# columns are centred and scaled once; with the intercept in W this leaves
# every column space unchanged and keeps G well conditioned. A control that
# is constant within a resample is dropped, as lm() drops aliased columns.
dial_prep_engine <- function(spec) {
  d <- as.data.frame(spec$data)
  for (v in names(d))
    if (inherits(d[[v]], "haven_labelled")) d[[v]] <- as.vector(unclass(d[[v]]))
  rhs <- if (length(spec$controls)) paste(spec$controls, collapse = " + ") else "1"
  W <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data = d)
  V <- cbind(y = as.numeric(d[[spec$Y]]), x = as.numeric(d[[spec$X]]),
             z = as.numeric(d[[spec$Z]]))
  if (nrow(W) != nrow(V)) stop("dial_prep_engine: the control matrix dropped rows")
  std <- function(M) {
    m <- colMeans(M); s <- apply(M, 2, stats::sd)
    s[!is.finite(s) | s == 0] <- 1
    sweep(sweep(M, 2, m, "-"), 2, s, "/")
  }
  int <- colnames(W) == "(Intercept)"
  if (any(!int)) W[, !int] <- std(W[, !int, drop = FALSE])
  A <- cbind(std(V), W)
  dimnames(A) <- NULL
  list(A = A, n = nrow(A), k = ncol(W))
}

# Prepared correlations (rho_xy, rho_xz, rho_yz) of the resample with row
# counts w, in the order of check_compatibility_simple(); NA if degenerate.
dial_prep_rho <- function(eng, w = rep(1, eng$n)) {
  G   <- crossprod(eng$A, eng$A * w)
  iw  <- 3L + seq_len(eng$k)
  Gvw <- G[1:3, iw, drop = FALSE]
  qw  <- qr(G[iw, iw, drop = FALSE], tol = 1e-10)
  cf  <- qr.coef(qw, t(Gvw)); cf[is.na(cf)] <- 0
  R   <- G[1:3, 1:3] - Gvw %*% cf
  syy <- R[1, 1]; sxx <- R[2, 2]; sxy <- R[1, 2]; syz <- R[1, 3]; sxz <- R[2, 3]
  det <- sxx * syy - sxy^2
  if (!is.finite(det) || det <= 0 || sxx <= 0 || syy <= 0) return(rep(NA_real_, 3))
  q <- (syy * sxz^2 - 2 * sxy * sxz * syz + sxx * syz^2) / det
  if (!is.finite(q) || q <= 0) return(rep(NA_real_, 3))
  c(sxy / sqrt(sxx * syy), sxz / sqrt(sxx * q), syz / sqrt(syy * q))
}

# Bootstrap covariance of the prepared correlations, with the preparation
# redone on every draw. Same seed handling, draw scheme and output as
# estimate_cov_corr_boot() in Fin_Empirical5_clean.R. 'base' holds the rows
# of an outer resample: the draws then resample that resample (the nested
# covariance of each resample of the bootstrap rule).
dial_cov_corr_boot_raw <- function(eng, B = 800, seed = 123, base = NULL) {
  if (!is.null(eng$data)) eng <- dial_prep_engine(eng)   # a spec was passed
  set.seed(seed)
  n  <- eng$n
  Rb <- matrix(NA_real_, nrow = B, ncol = 3)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (!is.null(base)) idx <- base[idx]
    Rb[b, ] <- dial_prep_rho(eng, tabulate(idx, nbins = n))
  }
  S <- stats::cov(Rb, use = "complete.obs")
  S <- (S + t(S)) / 2
  if (requireNamespace("Matrix", quietly = TRUE))
    S <- as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
  S
}

# The engine for a prepared frame d = (x, y, z), from a raw specification
# or an engine, or NULL if there is none or it does not reproduce d (same
# rows and the same three correlations).
.dial_engine_for <- function(d, raw) {
  if (is.null(raw)) return(NULL)
  eng <- if (!is.null(raw$A)) raw else
    tryCatch(dial_prep_engine(raw), error = function(e) NULL)
  if (is.null(eng) || eng$n != nrow(d)) return(NULL)
  r0 <- c(stats::cor(d$y, d$x), stats::cor(d$x, d$z), stats::cor(d$y, d$z))
  r1 <- dial_prep_rho(eng)
  if (any(!is.finite(r1)) || max(abs(r1 - r0)) > 1e-6) return(NULL)
  eng
}

# ---- Faster copies of the MCUB routines --------------------------------------
# Profiling shows that row-wise apply(M, 1, min) and apply(M, 1, which.min)
# and repeated evaluations of g() take most of the time in CIhybrid(). The
# copies below are the routines of Fin_Empirical5_clean.R with three changes
# that do not alter any number:
#   1. apply(M, 1, min / max / which.min) is replaced by max.col(..., "first"),
#      which compares exactly and breaks ties at the first column, as
#      which.min() does;
#   2. g(deltastar) is computed once per CIproj_p() call instead of twice;
#   3. CIcon_TNbounds() reuses those g() rows instead of recomputing them.
# The copies return bit-identical sets (checked against the originals with
# identical() on resampled data from the Di Tella, Galiani, Banerjee and
# Burde applications at alpha = 0.05 and 0.10, and on a full-sample call
# with B_large = 5,000). If a matrix contains NA or a standard error is zero,
# the copy hands the call to the original routine, so those cases also
# behave the same.
.dial_rowwhichmin <- function(M) max.col(-M, ties.method = "first")
.dial_rowmin <- function(M) M[cbind(seq_len(nrow(M)), max.col(-M, ties.method = "first"))]
.dial_rowmax <- function(M) M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]

.dial_CIhybrid_fast <- function(deltahat, deltaSigma, Al, Au,
                                alpha, alphac, eta,
                                B, Blarge,
                                tol, tol_r, index, g,
                                seed = 123) {
  set.seed(seed)
  k_dim <- length(deltahat)
  deltastar_demean_large <- MASS::mvrnorm(Blarge, mu = rep(0, k_dim), Sigma = deltaSigma)
  # Same expression as the original (including its vector recycling).
  deviation <- .dial_rowmax(abs(deltastar_demean_large / sqrt(diag(deltaSigma))))
  c_bd      <- quantile(deviation, 1 - eta)
  sigma_delta <- sqrt(diag(deltaSigma))
  Lambda      <- rbind(Al, Au)
  lambdaSigma <- Lambda %*% deltaSigma %*% t(Lambda)
  lambdasigma <- sqrt(diag(lambdaSigma))
  kk      <- nrow(Al)
  sigma_l <- lambdasigma[1:kk]
  sigma_u <- lambdasigma[(kk + 1):(2 * kk)]
  corr_all <- diag(1 / lambdasigma) %*% lambdaSigma %*% diag(1 / lambdasigma)
  corr_l   <- corr_all[1:kk, 1:kk]
  corr_u   <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]
  corr_m   <- corr_all[1:kk, (kk + 1):(2 * kk)]
  lb     <- deltahat - sigma_delta * c_bd
  ub     <- deltahat + sigma_delta * c_bd
  delta1 <- deltahat
  if (!is.null(index) && length(index) > 0) {
    lb[index] <- 0; ub[index] <- 0; delta1[index] <- 0
  }
  cl <- 0; cu <- qnorm(1 - alpha / 2); c <- (cl + cu) / 2
  delta_fea <- list(); c_fea <- numeric()
  obj_large <- function(delta, c_check) {
    (alpha - .dial_CIproj_p_fast(c_check, delta, alphac, c_bd,
                                 deltastar_demean_large,
                                 Al, Au, sigma_l, sigma_u,
                                 deltaSigma, corr_m, corr_l, corr_u,
                                 eta, tol_r, g)) * 100
  }
  k <- 1
  while ((cu - cl) > tol) {
    set.seed(seed + k); k <- k + 1
    deltastar_demean <- MASS::mvrnorm(B, mu = rep(0, k_dim), Sigma = deltaSigma)
    obj <- function(delta) {
      (alpha - .dial_CIproj_p_fast(c, delta, alphac, c_bd, deltastar_demean,
                                   Al, Au, sigma_l, sigma_u, deltaSigma,
                                   corr_m, corr_l, corr_u, eta, tol_r, g)) * 100
    }
    p1 <- obj_large(delta1, c)
    if (!is.finite(p1)) {
      cat("Warning: p1 is non-finite. delta1 =", round(delta1, 3), " c =", round(c, 3), "\n")
    }
    if (p1 >= 0) {
      f_optim <- tryCatch({
        optim(par = delta1, fn = obj, method = "L-BFGS-B",
              lower = lb, upper = ub, control = list(maxit = 1000))
      }, error = function(e) NULL)
      if (!is.null(f_optim)) delta1 <- f_optim$par
      p1 <- obj_large(delta1, c)
      if (p1 >= 0) cu <- c else {
        cl <- c
        delta_fea[[length(delta_fea) + 1]] <- delta1
        c_fea <- c(c, c_fea)
      }
    } else {
      cl <- c
      delta_fea[[length(delta_fea) + 1]] <- delta1
      c_fea <- c(c, c_fea)
    }
    c <- (cl + cu) / 2
  }
  cl <- c; cu <- qnorm(1 - alpha / 2)
  while ((cu - cl) > tol) {
    c  <- (cl + cu) / 2
    p  <- obj_large(delta1, c)
    if (p >= 0) cu <- c else cl <- c
  }
  lambdahat_l <- g(deltahat); lambdahat_u <- lambdahat_l
  c_LF  <- qnorm(1 - eta)
  CI_p  <- c(min(lambdahat_l - c * sigma_l), max(lambdahat_u + c * sigma_u))
  CI_c  <- .dial_CIcon_fast(deltahat, deltaSigma, Al, Au, c_LF, alphac, tol, tol_r, g)
  CI_h  <- c(min(CI_p[1], CI_c[1]), max(CI_p[2], CI_c[2]))
  list(CI_h = CI_h, CI_c = CI_c, CI_p = CI_p, delta1 = delta1, c = c)
}

.dial_CIproj_p_fast <- function(c, delta, alphac, c_bd,
                                deltastar_demean,
                                Al, Au,
                                sigma_l, sigma_u,
                                deltaSigma,
                                corr_m, corr_l, corr_u,
                                eta, tol_r, g) {
  lambda_l <- g(matrix(delta, nrow = 1))
  lambda_u <- lambda_l                                  # same g() value
  lb <- min(lambda_l); ub <- max(lambda_u); mb <- (lb + ub) / 2
  deltastar    <- sweep(deltastar_demean, 2, delta, "+")
  lambdastar_l <- g(deltastar)
  lambdastar_u <- lambdastar_l                          # same g() values
  if (anyNA(lambdastar_l) || anyNA(lambda_l) ||         # rare: use the original
      any(!is.finite(sigma_l)) || any(sigma_l == 0) ||
      any(!is.finite(sigma_u)) || any(sigma_u == 0))
    return(CIproj_p(c, delta, alphac, c_bd, deltastar_demean, Al, Au,
                    sigma_l, sigma_u, deltaSigma, corr_m, corr_l, corr_u,
                    eta, tol_r, g))
  delta_sigma <- sqrt(diag(deltaSigma))
  maxdev      <- .dial_rowmax(
    abs(deltastar_demean / matrix(delta_sigma, nrow = nrow(deltastar_demean),
                                  ncol = length(delta_sigma), byrow = TRUE)))
  ind_Delta <- as.numeric(maxdev <= c_bd)
  SL <- matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE)
  SU <- matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE)
  make_Tstar <- function(bound) {
    pmax(.dial_rowmin((lambdastar_l - bound) / SL),
         .dial_rowmin((bound - lambdastar_u) / SU))
  }
  Tstar_l <- make_Tstar(lb); Tstar_m <- make_Tstar(mb); Tstar_u <- make_Tstar(ub)
  ind_proj_l <- as.numeric(Tstar_l <= c)
  ind_proj_m <- as.numeric(Tstar_m <= c)
  ind_proj_u <- as.numeric(Tstar_u <= c)
  cLF <- qnorm(1 - eta)
  cond_block <- function(bound, select) {
    if (sum(select) == 0) return(rep(0, nrow(deltastar)))
    sel <- select == 1
    sub_deltastar <- deltastar[sel, , drop = FALSE]
    sub_l_l       <- lambdastar_l[sel, , drop = FALSE]
    sub_l_u       <- lambdastar_u[sel, , drop = FALSE]
    Tc1 <- .dial_rowmin((sub_l_l - bound) / SL[sel, , drop = FALSE])
    Tc2 <- .dial_rowmin((bound - sub_l_u) / SU[sel, , drop = FALSE])
    Tc  <- pmax(Tc1, Tc2)
    th_bounds <- .dial_CIcon_TNbounds_fast(bound, sub_deltastar, Al, Au, sigma_l, sigma_u,
                                           corr_m, corr_l, corr_u, cLF, cLF, tol_r, g,
                                           lam = sub_l_l)
    th1 <- th_bounds[[1]]; th2 <- th_bounds[[2]]
    phi     <- (pnorm(Tc) - pnorm(th1)) / (pnorm(th2) - pnorm(th1))
    sub_ind <- as.numeric(phi <= 1 - alphac)
    ind_c      <- rep(0, nrow(deltastar))
    ind_c[sel] <- sub_ind
    ind_c
  }
  selectl  <- (1 - ind_proj_l) * ind_Delta; ind_c_l  <- cond_block(lb, selectl)
  selectm  <- (1 - ind_proj_m) * ind_Delta; ind_c_m  <- cond_block(mb, selectm)
  selectu  <- (1 - ind_proj_u) * ind_Delta; ind_c_u  <- cond_block(ub, selectu)
  ind_l <- 1 - (1 - ind_c_l) * (1 - ind_proj_l)
  ind_m <- 1 - (1 - ind_c_m) * (1 - ind_proj_m)
  ind_u <- 1 - (1 - ind_c_u) * (1 - ind_proj_u)
  p1 <- mean((1 - ind_l * ind_m) * ind_Delta, na.rm = TRUE)
  p2 <- mean((1 - ind_u * ind_m) * ind_Delta, na.rm = TRUE)
  p  <- max(p1, p2, na.rm = TRUE) + eta
  if (!is.finite(p)) p <- NA_real_
  p
}

.dial_CIcon_fast <- function(deltahat, deltaSigma, Al, Au,
                             c_LF, alphac, tol, tol_r, g) {
  kk <- nrow(Al)
  Lambda       <- rbind(Al, Au)
  lambdaSigma  <- Lambda %*% deltaSigma %*% t(Lambda)
  lambda_sigma <- sqrt(diag(lambdaSigma))
  corr_all <- diag(1 / lambda_sigma) %*% lambdaSigma %*% diag(1 / lambda_sigma)
  corr_m   <- corr_all[1:kk, (kk + 1):(2 * kk)]
  corr_l   <- corr_all[1:kk, 1:kk]
  corr_u   <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]
  sigma_l  <- lambda_sigma[1:kk]
  sigma_u  <- lambda_sigma[(kk + 1):(2 * kk)]
  gval  <- as.numeric(g(matrix(deltahat, nrow = 1)))
  lb    <- gval - sigma_l * c_LF; ub <- gval + sigma_u * c_LF
  lb    <- min(lb); ub <- max(ub); mid <- (lb + ub) / 2
  lb_pt <- min(gval); ub_pt <- max(gval)
  rej   <- 1; theta <- lb
  while (rej == 1 && theta <= min(mid, lb_pt)) {
    Tcl <- min((gval - theta) / sigma_l); Tcu <- min((theta - gval) / sigma_u)
    Tc  <- max(Tcl, Tcu)
    th_bounds <- .dial_CIcon_TNbounds_fast(theta, matrix(deltahat, nrow = 1),
                                           Al, Au, sigma_l, sigma_u,
                                           corr_m, corr_l, corr_u,
                                           c_LF, c_LF, tol_r, g)
    th_1 <- th_bounds[[1]]; th_2 <- th_bounds[[2]]
    t   <- qnorm((1 - alphac) * pnorm(th_2) + alphac * pnorm(th_1))
    rej <- as.numeric(Tc > t); theta <- theta + tol
  }
  CI_lower <- theta
  rej   <- 1; theta <- ub
  while (rej == 1 && theta >= max(mid, ub_pt)) {
    Tcl <- min((gval - theta) / sigma_l); Tcu <- min((theta - gval) / sigma_u)
    Tc  <- max(Tcl, Tcu)
    th_bounds <- .dial_CIcon_TNbounds_fast(theta, matrix(deltahat, nrow = 1),
                                           Al, Au, sigma_l, sigma_u,
                                           corr_m, corr_l, corr_u,
                                           c_LF, c_LF, tol_r, g)
    th_1 <- th_bounds[[1]]; th_2 <- th_bounds[[2]]
    t   <- qnorm((1 - alphac) * pnorm(th_2) + alphac * pnorm(th_1))
    rej <- as.numeric(Tc > t); theta <- theta - tol
  }
  CI_upper <- theta
  c(CI_lower, CI_upper)
}

.dial_CIcon_TNbounds_fast <- function(theta, deltahat,
                                      Al, Au,
                                      sigma_l, sigma_u,
                                      corr_m, corr_l, corr_u,
                                      cLFl, cLFu, tol_r, g, lam = NULL) {
  B   <- nrow(deltahat)
  k_l <- nrow(Al); k_u <- nrow(Au)
  # lam = g(deltahat) when the caller has it already (same values)
  lambdahat_l <- if (is.null(lam)) g(deltahat) else lam
  lambdahat_u <- lambdahat_l
  Tlb <- sweep(lambdahat_l - theta, 2, sigma_l, "/")
  Tub <- sweep(theta - lambdahat_u, 2, sigma_u, "/")
  if (anyNA(Tlb) || anyNA(Tub))                         # rare: use the original
    return(CIcon_TNbounds(theta, deltahat, Al, Au, sigma_l, sigma_u,
                          corr_m, corr_l, corr_u, cLFl, cLFu, tol_r, g))
  bl <- .dial_rowwhichmin(Tlb); Tl <- Tlb[cbind(seq_len(B), bl)]
  bu <- .dial_rowwhichmin(Tub); Tu <- Tub[cbind(seq_len(B), bu)]
  corr_m_blb <- corr_m[bl, , drop = FALSE]
  corr_m_bub <- t(corr_m[, bu, drop = FALSE])
  tTS_l  <- matrix(1e10, nrow = B, ncol = k_u)
  tTS_ll <- (1 + corr_m_blb)^(-1) * (Tub + corr_m_blb * Tl)
  maskTS_l <- (1 + corr_m_blb) > tol_r
  tTS_l[maskTS_l] <- tTS_ll[maskTS_l]
  th_1_1 <- .dial_rowmin(tTS_l) * (Tl >= Tu)
  th_1_1[th_1_1 == 1e10] <- -1e10
  th_1_1 <- th_1_1 * (Tl >= Tu)
  tLFl <- rep(cLFl, B)
  CL_bl   <- corr_l[bl, , drop = FALSE]
  Tlb_min <- matrix(Tlb[cbind(1:B, bl)], nrow = B, ncol = k_l)
  Tl_rep  <- matrix(Tl, nrow = B, ncol = k_l)
  taux_1  <- (1 - CL_bl)^(-1) * (Tlb_min - CL_bl * Tl_rep)
  tB_l2   <- matrix(1e10, nrow = B, ncol = k_l)
  maskB_l <- (1 > CL_bl + tol_r)
  tB_l2[maskB_l] <- taux_1[maskB_l]
  th_2_1 <- pmin(tLFl, .dial_rowmin(tB_l2)) * (Tl >= Tu)
  tTS_u  <- matrix(1e10, nrow = B, ncol = k_l)
  tTS_u1 <- (1 + corr_m_bub)^(-1) * (Tlb + corr_m_bub * Tu)
  maskTS_u <- (1 + corr_m_bub) > tol_r
  tTS_u[maskTS_u] <- tTS_u1[maskTS_u]
  th_1_2 <- .dial_rowmin(tTS_u) * (Tl < Tu)
  th_1_2[th_1_2 == 1e10] <- -1e10
  th_1_2 <- th_1_2 * (Tl < Tu)
  tLFu <- rep(cLFu, B)
  CU_bu   <- corr_u[bu, , drop = FALSE]
  Tub_min <- matrix(Tub[cbind(1:B, bu)], nrow = B, ncol = k_u)
  Tu_rep  <- matrix(Tu, nrow = B, ncol = k_u)
  taux_2  <- (1 - CU_bu)^(-1) * (Tub_min - CU_bu * Tu_rep)
  tB_u2   <- matrix(1e10, nrow = B, ncol = k_u)
  maskB_u <- (1 > CU_bu + tol_r)
  tB_u2[maskB_u] <- taux_2[maskB_u]
  th_2_2 <- pmin(tLFu, .dial_rowmin(tB_u2)) * (Tl < Tu)
  th_1 <- th_1_1 + th_1_2
  th_2 <- th_2_1 + th_2_2
  list(th_1 = th_1, th_2 = th_2)
}

# ---- One MCUB set on one (resampled) data set -------------------------------
# Same steps as the MCUB call inside check_compatibility_simple(). Returns
# c(lower, upper), or c(NA, NA) if the computation fails.
# .dial_mcub_set():  from prepared (x, y, z) rows ("rows" unit)
# .dial_mcub_from(): from the correlations and their covariance ("raw" unit)
.dial_mcub_set <- function(x, y, z, rxu_range, alpha,
                           B = 500L, Blarge = 499L, B_cov = 200L, seed = 1L,
                           fast = TRUE) {
  deltahat <- c(stats::cor(y, x), stats::cor(x, z), stats::cor(y, z))
  if (any(!is.finite(deltahat))) return(c(NA_real_, NA_real_))
  S <- tryCatch(estimate_cov_corr_boot(x, y, z, B = B_cov, seed = seed),
                error = function(e) NULL)
  .dial_mcub_from(deltahat, S, rxu_range, alpha, B = B, Blarge = Blarge,
                  seed = seed, fast = fast)
}
.dial_mcub_from <- function(deltahat, S, rxu_range, alpha,
                            B = 500L, Blarge = 499L, seed = 1L, fast = TRUE) {
  if (any(!is.finite(deltahat))) return(c(NA_real_, NA_real_))
  if (is.null(S) || any(!is.finite(S))) return(c(NA_real_, NA_real_))
  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = 50L)
  g      <- function(delta) g_xu_safe(r_grid, delta)
  A <- t(vapply(r_grid, function(r)
    as.numeric(local_compute_gradient_safe(r, deltahat[1], deltahat[2], deltahat[3])),
    numeric(3L)))
  if (any(!is.finite(A))) return(c(NA_real_, NA_real_))
  ci_fun <- if (fast) .dial_CIhybrid_fast else CIhybrid
  res <- NULL
  invisible(utils::capture.output(
    res <- tryCatch(
      ci_fun(deltahat, S, A, A,
             alpha = alpha, alphac = 0.8 * alpha, eta = 0.001,
             B = B, Blarge = Blarge, tol = 1e-3, tol_r = 1e-3,
             index = NULL, g = g, seed = seed + 1L),
      error = function(e) NULL)))
  if (is.null(res)) c(NA_real_, NA_real_) else as.numeric(res$CI_h)
}

# ---- One bootstrap resample (runs in a parallel worker or in this session) --
.dial_boot_one <- function(b, d, seed, rxu_range, alpha, emp5, rule_path,
                           mcub, mcub_fun, raw = NULL) {
  # A fresh parallel worker starts empty: load the MCUB functions and the
  # fast copies above into its global environment.
  if (!exists("CIhybrid", mode = "function")) {
    assign("DIAL_SKIP_AUTORUN", TRUE, envir = globalenv())
    invisible(utils::capture.output(suppressMessages(suppressWarnings(
      source(emp5, local = globalenv())))))
  }
  if (!exists(".dial_CIhybrid_fast", mode = "function"))
    invisible(utils::capture.output(suppressMessages(suppressWarnings(
      source(rule_path, local = globalenv())))))
  s_b <- as.integer(seed) * 10000L + as.integer(b)
  set.seed(s_b, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  if (is.null(raw)) {                      # unit "rows": prepared frame
    n   <- nrow(d)
    idx <- sample.int(n, n, replace = TRUE)
    db  <- d[idx, , drop = FALSE]
    return(mcub_fun(db$x, db$y, db$z, rxu_range, alpha,
                    B = mcub$B, Blarge = mcub$Blarge, B_cov = mcub$B_cov,
                    seed = 1000000000L + s_b))
  }
  # Unit "raw" ('raw' is a dial_prep_engine()): redraw the rows of the
  # analysis data, redo the residualise-and-project step, and redo it again
  # on every draw of this resample's own covariance bootstrap.
  n   <- raw$n
  idx <- sample.int(n, n, replace = TRUE)
  deltahat <- dial_prep_rho(raw, tabulate(idx, nbins = n))
  if (any(!is.finite(deltahat))) return(c(NA_real_, NA_real_))
  S <- tryCatch(dial_cov_corr_boot_raw(raw, B = mcub$B_cov,
                                       seed = 1000000000L + s_b, base = idx),
                error = function(e) NULL)
  .dial_mcub_from(deltahat, S, rxu_range, alpha,
                  B = mcub$B, Blarge = mcub$Blarge, seed = 1000000000L + s_b)
}

# ---- The bootstrap rule -----------------------------------------------------
mcub_bootstrap_rule <- function(df, rxu_range, alpha,
                                seed      = 1L,
                                B_rule    = getOption("dial.boot_B", 500L),
                                tau2      = getOption("dial.tau2", 0.95),
                                zero_tol  = getOption("dial.zero_tol", 0),
                                mcub      = getOption("dial.boot_mcub"),
                                cache_dir = getOption("dial.boot_cache_dir", NA),
                                unit      = getOption("dial.boot_unit", "raw"),
                                raw       = attr(df, "dial_raw"),
                                verbose   = TRUE) {
  if (is.null(mcub)) mcub <- list(B = 500L, Blarge = 499L, B_cov = 200L)
  unit <- match.arg(unit, c("raw", "rows"))
  d <- as.data.frame(df)
  d <- d[stats::complete.cases(d[, c("x", "y", "z")]), c("x", "y", "z")]
  rownames(d) <- NULL
  n   <- nrow(d)
  rho <- c(stats::cor(d$y, d$x), stats::cor(d$x, d$z), stats::cor(d$y, d$z))
  # Unit "raw": the engine that redoes the preparation (raw may be a spec
  # from dial_raw_spec() or an engine already).
  eng <- NULL
  if (unit == "raw") {
    eng <- .dial_engine_for(d, raw)
    if (is.null(eng)) {
      warning("bootstrap rule: no raw data that reproduce this frame; ",
              "resampling the prepared rows instead.", call. = FALSE)
      unit <- "rows"
    }
  }
  key <- sprintf("%s_n%d_D%+.3f%+.3f_a%.3f_B%d_s%d_m%d-%d-%d_r%s",
                 if (unit == "raw") "v2" else "v1",
                 n, rxu_range[1], rxu_range[2], alpha, as.integer(B_rule),
                 as.integer(seed), as.integer(mcub$B), as.integer(mcub$Blarge),
                 as.integer(mcub$B_cov), paste(sprintf("%+.6f", rho), collapse = ""))
  if (unit == "raw") key <- paste0(key, "_uraw")
  use_cache  <- length(cache_dir) == 1L && !is.na(cache_dir) && nzchar(cache_dir)
  cache_file <- if (use_cache) file.path(cache_dir, paste0(key, ".rds")) else NULL

  t0 <- Sys.time()
  # Resample b has its own seed (seed * 10000 + b), so the first B' resamples
  # of a cached run with the same settings are the first B' resamples of
  # this run. Reuse them: take the first B_rule rows of a larger run, or
  # extend the largest smaller run.
  prev <- NULL
  if (use_cache && !file.exists(cache_file) && dir.exists(cache_dir)) {
    noB  <- function(f) sub("_B[0-9]+_", "_B#_", f)
    fs   <- list.files(cache_dir, pattern = "\\.rds$")
    fs   <- fs[noB(fs) == noB(paste0(key, ".rds"))]
    Bs   <- as.integer(sub(".*_B([0-9]+)_s.*", "\\1", fs))
    if (length(fs)) {
      big <- which(Bs >= B_rule)
      j   <- if (length(big)) big[which.min(Bs[big])] else which.max(Bs)
      prev <- tryCatch(readRDS(file.path(cache_dir, fs[j])), error = function(e) NULL)
      if (!is.null(prev) && nrow(prev) > B_rule) prev <- prev[seq_len(B_rule), , drop = FALSE]
    }
  }
  if (use_cache && file.exists(cache_file)) {
    ends <- readRDS(cache_file)
    from_cache <- TRUE
  } else if (!is.null(prev) && nrow(prev) == B_rule) {
    ends <- prev
    from_cache <- TRUE
    saveRDS(ends, cache_file)
  } else {
    from_cache <- FALSE
    b_todo <- seq.int(if (is.null(prev)) 1L else nrow(prev) + 1L, B_rule)
    emp5      <- getOption("dial.emp5_path")
    rule_path <- getOption("dial.boot_rule_path")
    # Parallel workers: the caller's own future plan if it has more than one
    # worker; otherwise a temporary multisession plan with
    # getOption("dial.boot_workers") workers (set by dial_boot_plan()). The
    # temporary plan ends when this call returns, so code that runs later
    # (for example DoubleML, whose mlr3 resampling follows the future plan)
    # keeps running sequentially, as before.
    use_par <- FALSE
    if (requireNamespace("future.apply", quietly = TRUE)) {
      if (future::nbrOfWorkers() > 1) {
        use_par <- TRUE
      } else {
        w <- as.integer(getOption("dial.boot_workers", 1L))
        if (!is.na(w) && w > 1L) {
          oplan <- future::plan(future::multisession, workers = w)
          on.exit(future::plan(oplan), add = TRUE)
          use_par <- TRUE
        }
      }
    }
    if (use_par) {
      out <- future.apply::future_lapply(
        b_todo, .dial_boot_one,
        d = d, seed = seed, rxu_range = rxu_range, alpha = alpha,
        emp5 = emp5, rule_path = rule_path, mcub = mcub,
        mcub_fun = .dial_mcub_set, raw = eng, future.seed = TRUE)
    } else {
      out <- lapply(b_todo, .dial_boot_one,
                    d = d, seed = seed, rxu_range = rxu_range, alpha = alpha,
                    emp5 = emp5, rule_path = rule_path, mcub = mcub,
                    mcub_fun = .dial_mcub_set, raw = eng)
    }
    ends <- rbind(prev, matrix(unlist(out), ncol = 2L, byrow = TRUE,
                               dimnames = list(NULL, c("lo", "hi"))))
    if (use_cache) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(ends, cache_file)
    }
  }
  ok    <- stats::complete.cases(ends)
  I     <- ifelse(ok, as.numeric(ends[, "lo"] <= zero_tol & ends[, "hi"] >= -zero_tol),
                  NA_real_)
  p_hat <- if (any(ok)) mean(I[ok]) else NA_real_
  n_reused <- if (is.null(prev)) 0L else nrow(prev)
  if (verbose)
    cat(sprintf("    bootstrap rule: p_hat_B = %.3f (%d/%d resamples of the %s, zero_tol = %.2f)%s, %.1f min\n",
                p_hat, sum(ok), B_rule,
                if (unit == "raw") "raw data" else "prepared rows", zero_tol,
                if (from_cache) " [cache]" else if (n_reused > 0L)
                  sprintf(" [%d reused from cache]", n_reused) else "",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(p_hat = p_hat, n_ok = sum(ok), valid = !is.na(p_hat) && p_hat >= tau2,
       tau2 = tau2, zero_tol = zero_tol, unit = unit, I = I, ends = ends)
}

# ---- Parallel workers -------------------------------------------------------
# Sets how many parallel workers the resamples use (default: physical cores
# minus one, at most 7, as in scripts 05 and 07). It does not change the
# session's future plan; see mcub_bootstrap_rule().
dial_boot_plan <- function(workers = NULL) {
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    message("future.apply is not installed: the bootstrap rule runs sequentially.")
    options(dial.boot_workers = 1L)
    return(invisible(FALSE))
  }
  if (is.null(workers)) {
    k <- parallel::detectCores(logical = FALSE)
    if (is.na(k)) k <- parallel::detectCores()
    workers <- max(1L, min(7L, as.integer(k) - 1L))
  }
  options(dial.boot_workers = as.integer(workers))
  invisible(TRUE)
}

# ---- Wrap check_compatibility_simple() --------------------------------------
.dial_install_boot_wrapper <- function(env = globalenv()) {
  if (!exists("check_compatibility_simple", envir = env, inherits = FALSE)) {
    warning("bootstrap_rule.R: check_compatibility_simple() not found. ",
            "Source Fin_Empirical5_clean.R first, then this file.")
    return(invisible(FALSE))
  }
  cur <- get("check_compatibility_simple", envir = env, inherits = FALSE)
  if (isTRUE(attr(cur, "dial_boot_wrapper"))) return(invisible(TRUE))
  assign(".dial_cc_simple_orig", cur, envir = env)
  wrapper <- function(df, i, alpha = 0.05, rxu_range = c(0, 0.8), ...) {
    unit <- match.arg(getOption("dial.boot_unit", "raw"), c("raw", "rows"))
    eng  <- NULL
    if (unit == "raw") {
      raw <- attr(df, "dial_raw")
      if (is.null(raw)) raw <- .dial_match_raw(df, getOption("dial.boot_raw"))
      dd  <- as.data.frame(df)
      dd  <- dd[stats::complete.cases(dd[, c("x", "y", "z")]), c("x", "y", "z")]
      eng <- .dial_engine_for(dd, raw)
      if (is.null(eng)) {
        warning("check_compatibility_simple(): no raw data that reproduce this ",
                "frame; its bootstraps resample the prepared rows.", call. = FALSE)
        unit <- "rows"
      }
    }
    # Unit "raw": while the original function runs, its covariance bootstraps
    # (estimate_cov_corr_boot() for the MCUB set and the bootstrap inside
    # ci_simple_union()) redraw the raw rows and redo the preparation, with
    # the same seeds and sizes. Both are restored when the call ends.
    if (!is.null(eng)) {
      home <- function(name) {
        e <- environment(.dial_cc_simple_orig)
        while (!identical(e, emptyenv())) {
          if (exists(name, envir = e, inherits = FALSE)) return(e)
          e <- parent.env(e)
        }
        NULL
      }
      e_cov <- home("estimate_cov_corr_boot"); e_csu <- home("ci_simple_union")
      if (!is.null(e_cov)) {
        cov_old <- get("estimate_cov_corr_boot", envir = e_cov, inherits = FALSE)
        assign("estimate_cov_corr_boot", function(x, y, z, B = 800, seed = 123) {
          if (length(x) != eng$n) return(cov_old(x, y, z, B = B, seed = seed))
          dial_cov_corr_boot_raw(eng, B = B, seed = seed)
        }, envir = e_cov)
        on.exit(assign("estimate_cov_corr_boot", cov_old, envir = e_cov), add = TRUE)
      }
      if (!is.null(e_csu)) {
        csu_old <- get("ci_simple_union", envir = e_csu, inherits = FALSE)
        assign("ci_simple_union", function(x, y, z, ..., Sigma_rho = NULL,
                                           cov_method = c("bootstrap", "diag"),
                                           B_boot = 800, seed = 123) {
          cov_method <- match.arg(cov_method)
          if (is.null(Sigma_rho) && cov_method == "bootstrap" && length(x) == eng$n)
            Sigma_rho <- dial_cov_corr_boot_raw(eng, B = B_boot, seed = seed)
          csu_old(x, y, z, ..., Sigma_rho = Sigma_rho, cov_method = cov_method,
                  B_boot = B_boot, seed = seed)
        }, envir = e_csu)
        on.exit(assign("ci_simple_union", csu_old, envir = e_csu), add = TRUE)
      }
    }
    out <- .dial_cc_simple_orig(df, i, alpha = alpha, rxu_range = rxu_range, ...)
    if (is.null(out) || !NROW(out)) return(out)
    out$cov_unit <- unit
    if (!isTRUE(getOption("dial.boot_rule", TRUE))) return(out)
    br <- mcub_bootstrap_rule(df, rxu_range = rxu_range, alpha = alpha,
                              seed = 5000L + as.integer(i), unit = unit, raw = eng)
    out$p_hat_B    <- round(br$p_hat, 3)
    out$n_boot_ok  <- br$n_ok
    out$boot_unit  <- br$unit
    out$Valid_tau2 <- if (is.na(br$p_hat)) NA_character_ else
      if (br$valid) "Valid" else "Invalid"
    out
  }
  attr(wrapper, "dial_boot_wrapper") <- TRUE
  environment(wrapper) <- env
  assign("check_compatibility_simple", wrapper, envir = env)
  invisible(TRUE)
}
.dial_install_boot_wrapper()
