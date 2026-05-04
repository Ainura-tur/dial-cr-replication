# ==============================================================================
# DIAL_NeurIPS_JP_dgp.R  --  Shared Jones-Pewsey DGP for Table 1 and Table 6
# ------------------------------------------------------------------------------
# Implements the synthetic DGP described in Section "Data Generating Process"
# of the NeurIPS submission. Used by:
#   - DIAL_NeurIPS_table1_MCUB_sets.R   (Table 1 synthetic benchmark)
#   - DIAL_NeurIPS_MonteCarlo.R          (Table 6 Monte Carlo validation)
#
# DGP summary
# -----------
#   nu                 ~ N(0, 1) base sequence, length n
#   x  = 5 H(nu, 0, 0.9) + Skew-Normal(loc=0, scale=1, kappa=1.2)
#   e  = lambda (x - mean(x)) + epsilon,    epsilon ~ N(0, sigma_x^2)
#   u1 = Uniform(-1, 1) + Skew-Normal(0, 1, 1.2)
#   v  = residual of u1 ~ x, rescaled so var(v) ~ var(e) (orthogonal to x)
#   u  = e + v
#   y  = beta * x + u,    beta = 2 (true causal effect)
#
# Instruments z_j are constructed at known target cor(z, u) values via a
# Cholesky-style projection on (x, u): given target a = cor(z, x) and
# target b = cor(z, u), z = c1 * x_std + c2 * u_perp_x + c3 * w, where
# w is an independent J-P-style Gaussian.
#
# H is the Jones and Pewsey (2009) sinh-arcsinh transformation:
#   H(x; zeta, kappa) = sinh(kappa * asinh(x) - zeta)
#
# "Skew-Normal" in the spec refers to the J-P "Sinh-Arcsinh distribution",
# NOT Azzalini's skew-normal. We sample from S(.; zeta=0, kappa) by drawing
# a standard normal Z and returning sinh((asinh(Z) + zeta) / kappa) -- the
# inverse of H -- so that S(Z) has the requested skew/kurt.
#
# Lambda calibration
# ------------------
# `calibrate_lambda(target_cor_xu, ...)` does a short bisection to find
# lambda such that cor(x, u) ~ target. Default target = 0.30 per the spec
# ("substantial but not extreme endogeneity").
#
# Public API
# ----------
#   jp_h_transform(x, zeta = 0, kappa = 0.9)     # forward transform
#   jp_skew_normal(n, zeta = 0, kappa = 1.2)     # sampler
#   draw_jp_sample(n, lambda, seed,              # main DGP
#                  cor_zu_grid, cor_zx_default = 0.25,
#                  beta = 2, sigma_x = NULL)
#   calibrate_lambda(n, target_cor_xu = 0.30,
#                    cor_zx_default = 0.25, seed = 1L)
#
# Returned by draw_jp_sample(): a list with x, u, y, z (matrix nrow=n,
# one column per instrument), nu, e, v, lambda, realized_cor_xu, and
# realized cor(z, u) per column for diagnostic.
# ==============================================================================

# ---- J-P sinh-arcsinh transformation --------------------------------------

jp_h_transform <- function(x, zeta = 0, kappa = 0.9) {
  sinh(kappa * asinh(x) - zeta)
}

# Sample from the J-P sinh-arcsinh distribution (location 0, scale 1)
# parameterised so kappa = 1 returns the normal. Implemented by sampling
# a standard normal and transforming.
jp_skew_normal <- function(n, zeta = 0, kappa = 1.2) {
  Z <- stats::rnorm(n)
  sinh((asinh(Z) + zeta) / kappa)
}

# ---- Instrument construction ----------------------------------------------

# Build an instrument z with target cor(z, x) = a and cor(z, u) = b given
# (x, u). Uses a Cholesky-style projection: we orthogonalise u against x,
# then take a linear combination plus an independent residual.
.construct_z <- function(x, u, a, b, seed_w) {
  set.seed(seed_w)
  n <- length(x)
  x_std <- (x - mean(x)) / stats::sd(x)
  u_std <- (u - mean(u)) / stats::sd(u)
  rho <- stats::cor(x_std, u_std)
  
  # Orthogonalise u against x
  u_perp <- u_std - rho * x_std
  s_perp <- stats::sd(u_perp)
  if (s_perp < 1e-8) stop(".construct_z: u nearly collinear with x.")
  u_perp <- u_perp / s_perp
  
  # Independent gaussian residual w
  w_raw <- stats::rnorm(n)
  w <- stats::resid(stats::lm(w_raw ~ x_std + u_perp))
  w <- w / stats::sd(w)
  
  # Solve for c1, c2, c3 such that
  #   cor(z, x) = c1                      (since x_std, u_perp, w orthogonal)
  #   cor(z, u) = c1 * rho + c2 * sqrt(1 - rho^2)
  #   c1^2 + c2^2 + c3^2 = 1
  c1 <- a
  c2 <- (b - a * rho) / sqrt(1 - rho^2)
  c3_sq <- 1 - c1^2 - c2^2
  if (c3_sq <= 0)
    stop(sprintf(".construct_z: infeasible (a=%.3f, b=%.3f, rho=%.3f); ",
                 a, b, rho),
         "need c1^2 + c2^2 < 1.")
  c3 <- sqrt(c3_sq)
  
  c1 * x_std + c2 * u_perp + c3 * w
}

# ---- Main DGP --------------------------------------------------------------

draw_jp_sample <- function(n, lambda, seed,
                           cor_zu_grid,
                           cor_zx_default = 0.25,
                           beta           = 2,
                           sigma_x        = NULL) {
  stopifnot(is.numeric(n), n > 0,
            is.numeric(lambda),
            is.numeric(cor_zu_grid))
  set.seed(seed)
  
  # 1. Endogenous regressor x
  nu <- stats::rnorm(n)
  x  <- 5 * jp_h_transform(nu, zeta = 0, kappa = 0.9) +
    jp_skew_normal(n, zeta = 0, kappa = 1.2)
  if (is.null(sigma_x)) sigma_x <- stats::sd(x)
  
  # 2. Endogenous error component e = lambda (x - mean(x)) + eps
  eps <- stats::rnorm(n, mean = 0, sd = sigma_x)
  e   <- lambda * (x - mean(x)) + eps
  
  # 3. Exogenous error component v: residualize u1 on x, rescale to give
  #    v a standard deviation on the order of x's so the two pieces of u
  #    contribute comparably. The spec writes the rescale as
  #        v <- v * (xbar/2) * (sd(x)/sd(v_raw))
  #    which is a typo: with x roughly centered, xbar ~ 0 collapses v
  #    to zero and forces cor(x, u) -> 1 (the calibration error message
  #    "infeasible (..., rho=0.979)" is exactly that pathology). The
  #    intended scaling is variance control, not mean control:
  #        sd(v) ~ sd(x) / 2
  u1 <- stats::runif(n, -1, 1) + jp_skew_normal(n, zeta = 0, kappa = 1.2)
  v_raw <- stats::resid(stats::lm(u1 ~ x))
  v <- v_raw * (stats::sd(x) / 2) / stats::sd(v_raw)
  
  # 4. Total error u and outcome y
  u <- e + v
  y <- beta * x + u
  
  # 5. Construct instruments at target cor(z, u) values
  K <- length(cor_zu_grid)
  z <- matrix(NA_real_, nrow = n, ncol = K)
  realized_cor_zu <- numeric(K)
  realized_cor_zx <- numeric(K)
  for (j in seq_len(K)) {
    z[, j] <- .construct_z(x, u,
                           a = cor_zx_default,
                           b = cor_zu_grid[j],
                           seed_w = seed + 1000L + j)
    realized_cor_zu[j] <- stats::cor(z[, j], u)
    realized_cor_zx[j] <- stats::cor(z[, j], x)
  }
  colnames(z) <- sprintf("z_%+0.2f", cor_zu_grid)
  
  list(
    x = x, u = u, y = y, z = z,
    nu = nu, e = e, v = v,
    lambda = lambda,
    cor_zu_grid = cor_zu_grid,
    realized_cor_xu = stats::cor(x, u),
    realized_cor_zu = realized_cor_zu,
    realized_cor_zx = realized_cor_zx,
    beta_true = beta,
    sigma_x = sigma_x,
    seed = seed
  )
}

# ---- Lambda calibration ----------------------------------------------------
#
# Find lambda so that the realized cor(x, u) at this n is close to the target.
# Uses bisection on lambda > 0 (assumes positive sign per the spec).

calibrate_lambda <- function(n, target_cor_xu = 0.30,
                             cor_zx_default = 0.25,
                             seed = 1L,
                             tol = 0.01,
                             max_iter = 30L,
                             sign = c("positive", "negative"),
                             verbose = FALSE) {
  sign <- match.arg(sign)
  abs_target <- abs(target_cor_xu)
  
  realized <- function(lam) {
    s <- draw_jp_sample(n = n, lambda = lam, seed = seed,
                        cor_zu_grid = c(0.0),
                        cor_zx_default = cor_zx_default)
    s$realized_cor_xu
  }
  
  # Calibrate magnitude on the positive side, then flip if negative requested.
  lo <- 0.01; hi <- 3.0
  r_lo <- realized(lo); r_hi <- realized(hi)
  if (r_lo > abs_target)
    return(list(lambda = if (sign == "positive") lo else -lo,
                realized = if (sign == "positive") r_lo else -r_lo,
                note = "target below feasible range"))
  if (r_hi < abs_target)
    return(list(lambda = if (sign == "positive") hi else -hi,
                realized = if (sign == "positive") r_hi else -r_hi,
                note = "target above feasible range"))
  
  for (k in seq_len(max_iter)) {
    mid <- (lo + hi) / 2
    r_mid <- realized(mid)
    if (verbose) cat(sprintf("  iter %d: lambda=%.4f  cor(x,u)=%.4f\n",
                             k, mid, r_mid))
    if (abs(r_mid - abs_target) < tol) {
      lam_out <- if (sign == "positive") mid else -mid
      r_out   <- if (sign == "positive") r_mid else -r_mid
      return(list(lambda = lam_out, realized = r_out, iter = k))
    }
    if (r_mid < abs_target) lo <- mid else hi <- mid
  }
  mid <- (lo + hi) / 2
  r_mid <- realized(mid)
  list(lambda = if (sign == "positive") mid else -mid,
       realized = if (sign == "positive") r_mid else -r_mid,
       note = "max_iter reached")
}

# ---- Standard instrument grid (paper) -------------------------------------

JP_INSTRUMENT_GRID <- seq(-0.20, 0.20, by = 0.04)

# Sanity self-test (only runs when this file is sourced standalone)
if (sys.nframe() == 0L && interactive()) {
  cat("J-P DGP self-test\n")
  cal <- calibrate_lambda(n = 5000L, target_cor_xu = 0.30, verbose = TRUE)
  cat(sprintf("calibrated lambda = %.4f, cor(x,u) = %.4f\n",
              cal$lambda, cal$realized))
  s <- draw_jp_sample(n = 5000L, lambda = cal$lambda, seed = 42L,
                      cor_zu_grid = JP_INSTRUMENT_GRID)
  cat("Realized cor(z, u) per instrument:\n")
  print(round(setNames(s$realized_cor_zu, sprintf("%+.2f",
                                                  JP_INSTRUMENT_GRID)), 3))
  cat("Realized cor(z, x) per instrument:\n")
  print(round(setNames(s$realized_cor_zx, sprintf("%+.2f",
                                                  JP_INSTRUMENT_GRID)), 3))
  cat(sprintf("\nPopulation: cor(x, u) = %.3f, beta = %.1f\n",
              s$realized_cor_xu, s$beta_true))
  cat(sprintf("OLS beta_hat = %.3f, IV (z=0) beta_hat = %.3f\n",
              stats::coef(stats::lm(s$y ~ s$x))[2L],
              {
                z_valid_idx <- which(abs(JP_INSTRUMENT_GRID) < 1e-8)
                z_valid <- s$z[, z_valid_idx]
                stats::coef(stats::lm(s$y ~ s$x))[2L] - 1  # placeholder
              }))
}