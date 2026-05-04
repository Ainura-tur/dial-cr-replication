#================================================================
# Fin_Empirical5_clean.R  —  cleaned for DIAL / NeurIPS CR test project
#
# CHANGES FROM ORIGINAL
# ─────────────────────────────────────────────────────────────────
# [CLEAN-1]  AUTO-RUN GUARD added.
#            All top-level execution blocks are now wrapped in:
#              if (.dial_autorun) { ... }
#            Set DIAL_SKIP_AUTORUN <- TRUE before sourcing to suppress
#            auto-execution when using DIAL_NEURIPS_CR_wrapper.R.
#
# [CLEAN-2]  DUPLICATE run_all_applications() REMOVED.
#
# [CLEAN-3]  STALE GALIANI EXAMPLE CALL REMOVED.
#
# [CLEAN-4]  check_compatibility() -> check_compatibility_simple().
#
# [CLEAN-5]  DEAD z <- z0 LINE REMOVED in run_galiani_analysis_fixed.
#
# PHASE 1 CORRECTIONS (DIAL_corrections_todo.md)
# ─────────────────────────────────────────────────────────────────
# [P1-BANERJEE]  rxu_range corrected from c(-0.8, -0.6) to c(-0.6, -0.4).
#                NOTE: PDF page-15 comment says canonical was (-0.6, -0.5).
#                Run DIAL_Phase1_runner.R to compare both and pick the
#                upper bound before updating Table 6.
#
# [P1-BURDE]     run_burde_sensitivity_fixed() ADDED — new function that
#                runs the Burde & Linden analysis with rxu_range = c(-0.4, 0)
#                to test whether the identified set crosses zero when the
#                negative half of the domain is opened.  The existing
#                run_burde_analysis_fixed() (canonical c(0.0, 0.4)) is
#                unchanged.
#
# [P1-GALIANI]   No code change. rxu_range = c(-0.2, 0) confirmed correct.
#
# [P1-DITELLA]   No code change. rxu_range = c(0, 0.25) confirmed correct.
#
# SESSION STARTUP (recommended):
#   DIAL_SKIP_AUTORUN <- TRUE
#   source("Fin_sim3_clean.R")
#   source("Fin_Empirical5_clean.R")     # loads functions only
#   source("DIAL_NEURIPS_CR_wrapper.R")
#   res <- run_all_scenarios()
#================================================================

# Load required packages at the top level
library(MASS)
library(ivreg)
library(lmtest)
library(sandwich)
library(haven)
library(dplyr)
library(Matrix)

#================================================================
# UTILITY FUNCTIONS
#================================================================

g_fun <- function(r, rho_xy, rho_xz, rho_yz) {
  rho_xz * r - (rho_xy * rho_xz - rho_yz) *
    sqrt((1 - r^2) / (1 - rho_xy^2))
}

g_grad <- function(r, rho_xy, rho_xz, rho_yz) {
  S <- sqrt((1 - r^2) / (1 - rho_xy^2))
  c(
    -rho_xz * S + (rho_xy * rho_xz - rho_yz) * S * rho_xy / (1 - rho_xy^2),
    r - rho_xy * S,
    S
  )
}

g_xu_safe <- function(r_xu, delta, tol = 1e-10) {
  if (is.null(dim(delta))) {
    delta <- matrix(delta, nrow = 1)
  } else if (is.matrix(delta) && ncol(delta) != 3) {
    if (nrow(delta) == 3) {
      delta <- t(delta)
    } else {
      stop("delta must be a matrix with 3 columns: [rho_xy, rho_xz, rho_yz]")
    }
  }
  if (ncol(delta) != 3) stop("delta must have 3 columns: [rho_xy, rho_xz, rho_yz]")
  
  rho_xy <- delta[, 1]
  rho_xz <- delta[, 2]
  rho_yz <- delta[, 3]
  
  B <- nrow(delta)
  K <- length(r_xu)
  
  rho_xy_mat <- matrix(rho_xy, nrow = B, ncol = K)
  rho_xz_mat <- matrix(rho_xz, nrow = B, ncol = K)
  rho_yz_mat <- matrix(rho_yz, nrow = B, ncol = K)
  r_xu_mat   <- matrix(r_xu,   nrow = B, ncol = K, byrow = TRUE)
  
  rho_xy_mat <- pmin(pmax(rho_xy_mat, -1 + tol), 1 - tol)
  r_xu_mat   <- pmin(pmax(r_xu_mat,   -1 + tol), 1 - tol)
  
  denom     <- sqrt(pmax(1 - rho_xy_mat^2, tol))
  sqrt_term <- sqrt(pmax(1 - r_xu_mat^2, 0))
  
  gval <- rho_xz_mat * r_xu_mat -
    (rho_xy_mat * rho_xz_mat - rho_yz_mat) / denom * sqrt_term
  
  gval[!is.finite(gval)] <- NA_real_
  return(gval)
}

local_compute_gradient_safe <- function(r_xu, rho_xy, rho_xz, rho_yz, tol = 1e-10) {
  rho_xy <- pmin(pmax(rho_xy, -1 + tol), 1 - tol)
  r_xu   <- pmin(pmax(r_xu,   -1 + tol), 1 - tol)
  
  sqrt_term <- sqrt(pmax((1 - r_xu^2) / (1 - rho_xy^2), 0))
  
  dg_drho_xy <- -rho_xz * sqrt_term +
    (rho_xy * rho_xz - rho_yz) * sqrt_term * rho_xy / (1 - rho_xy^2)
  dg_drho_xz <- r_xu - rho_xy * sqrt_term
  dg_drho_yz <- sqrt_term
  
  grad <- rbind(dg_drho_xy, dg_drho_xz, dg_drho_yz)
  grad[!is.finite(grad)] <- NA_real_
  
  return(grad)
}

estimate_cov_corr_boot <- function(x, y, z, B = 800, seed = 123) {
  set.seed(seed)
  dat <- cbind(x, y, z)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  n   <- nrow(dat)
  
  Rb <- matrix(NA_real_, nrow = B, ncol = 3)
  for (b in 1:B) {
    idx <- sample.int(n, n, replace = TRUE)
    xb <- dat[idx, 1]; yb <- dat[idx, 2]; zb <- dat[idx, 3]
    Rb[b, ] <- c(cor(xb, yb), cor(xb, zb), cor(yb, zb))
  }
  
  S <- cov(Rb, use = "complete.obs")
  S <- (S + t(S)) / 2
  
  if (requireNamespace("Matrix", quietly = TRUE)) {
    S <- as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
  }
  S
}

#================================================================
# ci_simple_union
#================================================================

ci_simple_union <- function(x, y, z,
                            rxu_range   = c(0, 0.8),
                            grid_length = 80,
                            alpha       = 0.05,
                            Sigma_rho   = NULL,
                            cov_method  = c("bootstrap", "diag"),
                            B_boot      = 800,
                            seed        = 123,
                            tol         = 1e-10,
                            make_psd    = TRUE) {
  
  cov_method <- match.arg(cov_method)
  
  dat <- data.frame(x = x, y = y, z = z)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  n   <- nrow(dat)
  if (n < 30) stop("ci_simple_union: need at least 30 complete observations.")
  
  rho_xy  <- suppressWarnings(cor(dat$x, dat$y))
  rho_xz  <- suppressWarnings(cor(dat$x, dat$z))
  rho_yz  <- suppressWarnings(cor(dat$y, dat$z))
  deltahat <- c(rho_xy, rho_xz, rho_yz)
  names(deltahat) <- c("rho_xy", "rho_xz", "rho_yz")
  
  if (any(!is.finite(deltahat))) {
    stop("ci_simple_union: non-finite sample correlations.")
  }
  
  clamp    <- function(a) pmin(pmax(a, -1 + tol), 1 - tol)
  rho_xy_c <- clamp(rho_xy)
  
  if (is.null(Sigma_rho)) {
    if (cov_method == "bootstrap") {
      set.seed(seed)
      Rb <- matrix(NA_real_, nrow = B_boot, ncol = 3)
      for (b in seq_len(B_boot)) {
        idx <- sample.int(n, n, replace = TRUE)
        xb  <- dat$x[idx]; yb <- dat$y[idx]; zb <- dat$z[idx]
        Rb[b, ] <- c(suppressWarnings(cor(xb, yb)),
                     suppressWarnings(cor(xb, zb)),
                     suppressWarnings(cor(yb, zb)))
      }
      Rb        <- Rb[complete.cases(Rb), , drop = FALSE]
      if (nrow(Rb) < 50) stop("ci_simple_union: too many NA bootstrap draws.")
      Sigma_rho <- cov(Rb)
      Sigma_rho <- (Sigma_rho + t(Sigma_rho)) / 2
      if (make_psd && requireNamespace("Matrix", quietly = TRUE)) {
        Sigma_rho <- as.matrix(Matrix::nearPD(Sigma_rho, corr = FALSE)$mat)
      }
    } else {
      var_xy    <- (1 - rho_xy^2)^2 / n
      var_xz    <- (1 - rho_xz^2)^2 / n
      var_yz    <- (1 - rho_yz^2)^2 / n
      Sigma_rho <- diag(c(var_xy, var_xz, var_yz))
    }
  } else {
    if (!is.matrix(Sigma_rho) || any(dim(Sigma_rho) != c(3, 3))) {
      stop("ci_simple_union: Sigma_rho must be a 3x3 matrix.")
    }
    Sigma_rho <- (Sigma_rho + t(Sigma_rho)) / 2
  }
  
  g_fun_local <- function(r, rho_xy, rho_xz, rho_yz) {
    rho_xy <- clamp(rho_xy); r <- clamp(r)
    rho_xz * r - (rho_xy * rho_xz - rho_yz) * sqrt((1 - r^2) / (1 - rho_xy^2))
  }
  
  g_grad_local <- function(r, rho_xy, rho_xz, rho_yz) {
    rho_xy <- clamp(rho_xy); r <- clamp(r)
    S <- sqrt((1 - r^2) / (1 - rho_xy^2))
    c(-rho_xz * S + (rho_xy * rho_xz - rho_yz) * S * rho_xy / (1 - rho_xy^2),
      r - rho_xy * S,
      S)
  }
  
  D_grid <- clamp(seq(rxu_range[1], rxu_range[2], length.out = grid_length))
  
  g_vals <- vapply(D_grid, g_fun_local, numeric(1),
                   rho_xy = rho_xy, rho_xz = rho_xz, rho_yz = rho_yz)
  
  idx_l  <- which.min(g_vals); idx_u  <- which.max(g_vals)
  r_l    <- D_grid[idx_l];     r_u    <- D_grid[idx_u]
  plug_in <- c(g_vals[idx_l], g_vals[idx_u])
  
  grad_l <- g_grad_local(r_l, rho_xy, rho_xz, rho_yz)
  grad_u <- g_grad_local(r_u, rho_xy, rho_xz, rho_yz)
  
  var_l  <- max(as.numeric(t(grad_l) %*% Sigma_rho %*% grad_l), 0)
  var_u  <- max(as.numeric(t(grad_u) %*% Sigma_rho %*% grad_u), 0)
  
  crit <- qnorm(1 - alpha / 4)
  
  CI_lower <- plug_in[1] - crit * sqrt(var_l)
  CI_upper <- plug_in[2] + crit * sqrt(var_u)
  
  list(CI = c(CI_lower, CI_upper), plug_in = plug_in,
       r_star = c(r_l, r_u), rho_hat = deltahat,
       Sigma_rho = Sigma_rho, crit = crit)
}

#================================================================
# MCUB FUNCTIONS
#================================================================

CIhybrid <- function(deltahat, deltaSigma, Al, Au,
                     alpha, alphac, eta,
                     B, Blarge,
                     tol, tol_r, index, g,
                     seed = 123) {
  
  set.seed(seed)
  
  k_dim <- length(deltahat)
  deltastar_demean_large <- MASS::mvrnorm(Blarge, mu = rep(0, k_dim), Sigma = deltaSigma)
  deviation <- apply(abs(deltastar_demean_large / sqrt(diag(deltaSigma))), 1, max)
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
    (alpha - CIproj_p(c_check, delta, alphac, c_bd,
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
      (alpha - CIproj_p(c, delta, alphac, c_bd, deltastar_demean,
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
  CI_c  <- CIcon(deltahat, deltaSigma, Al, Au, c_LF, alphac, tol, tol_r, g)
  CI_h  <- c(min(CI_p[1], CI_c[1]), max(CI_p[2], CI_c[2]))
  
  list(CI_h = CI_h, CI_c = CI_c, CI_p = CI_p, delta1 = delta1, c = c)
}

CIproj_p <- function(c, delta, alphac, c_bd,
                     deltastar_demean,
                     Al, Au,
                     sigma_l, sigma_u,
                     deltaSigma,
                     corr_m, corr_l, corr_u,
                     eta, tol_r, g) {
  
  lambda_l <- g(matrix(delta, nrow = 1))
  lambda_u <- g(matrix(delta, nrow = 1))
  lb <- min(lambda_l); ub <- max(lambda_u); mb <- (lb + ub) / 2
  
  deltastar    <- sweep(deltastar_demean, 2, delta, "+")
  lambdastar_l <- g(deltastar)
  lambdastar_u <- g(deltastar)
  
  delta_sigma <- sqrt(diag(deltaSigma))
  maxdev      <- apply(
    abs(deltastar_demean / matrix(delta_sigma, nrow = nrow(deltastar_demean),
                                  ncol = length(delta_sigma), byrow = TRUE)),
    1, max)
  ind_Delta <- as.numeric(maxdev <= c_bd)
  
  make_Tstar <- function(bound) {
    pmax(
      apply((lambdastar_l - bound) /
              matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE),
            1, min),
      apply((bound - lambdastar_u) /
              matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE),
            1, min)
    )
  }
  
  Tstar_l <- make_Tstar(lb); Tstar_m <- make_Tstar(mb); Tstar_u <- make_Tstar(ub)
  
  ind_proj_l <- as.numeric(Tstar_l <= c)
  ind_proj_m <- as.numeric(Tstar_m <= c)
  ind_proj_u <- as.numeric(Tstar_u <= c)
  
  cLF <- qnorm(1 - eta)
  
  cond_block <- function(bound, select) {
    if (sum(select) == 0) return(rep(0, nrow(deltastar)))
    sub_deltastar <- deltastar[select == 1, , drop = FALSE]
    sub_l_l       <- lambdastar_l[select == 1, , drop = FALSE]
    sub_l_u       <- lambdastar_u[select == 1, , drop = FALSE]
    
    Tc1 <- apply((sub_l_l - bound) /
                   matrix(sigma_l, nrow = nrow(sub_l_l), ncol = length(sigma_l), byrow = TRUE),
                 1, min)
    Tc2 <- apply((bound - sub_l_u) /
                   matrix(sigma_u, nrow = nrow(sub_l_u), ncol = length(sigma_u), byrow = TRUE),
                 1, min)
    Tc  <- pmax(Tc1, Tc2)
    
    th_bounds <- CIcon_TNbounds(bound, sub_deltastar, Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u, cLF, cLF, tol_r, g)
    th1 <- th_bounds[[1]]; th2 <- th_bounds[[2]]
    
    phi     <- (pnorm(Tc) - pnorm(th1)) / (pnorm(th2) - pnorm(th1))
    sub_ind <- as.numeric(phi <= 1 - alphac)
    
    ind_c         <- rep(0, nrow(deltastar))
    ind_c[select == 1] <- sub_ind
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
  return(p)
}

CIcon <- function(deltahat, deltaSigma, Al, Au,
                  c_LF, alphac, tol, tol_r, g) {
  
  kk <- nrow(Al)
  
  Lambda      <- rbind(Al, Au)
  lambdaSigma <- Lambda %*% deltaSigma %*% t(Lambda)
  lambda_sigma <- sqrt(diag(lambdaSigma))
  
  corr_all <- diag(1 / lambda_sigma) %*% lambdaSigma %*% diag(1 / lambda_sigma)
  corr_m   <- corr_all[1:kk, (kk + 1):(2 * kk)]
  corr_l   <- corr_all[1:kk, 1:kk]
  corr_u   <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]
  
  sigma_l  <- lambda_sigma[1:kk]
  sigma_u  <- lambda_sigma[(kk + 1):(2 * kk)]
  
  gval  <- as.numeric(g(matrix(deltahat, nrow = 1)))
  
  lb    <- gval - sigma_l * c_LF; ub    <- gval + sigma_u * c_LF
  lb    <- min(lb); ub    <- max(ub); mid   <- (lb + ub) / 2
  lb_pt <- min(gval); ub_pt <- max(gval)
  
  rej   <- 1; theta <- lb
  
  while (rej == 1 && theta <= min(mid, lb_pt)) {
    Tcl <- min((gval - theta) / sigma_l); Tcu <- min((theta - gval) / sigma_u)
    Tc  <- max(Tcl, Tcu)
    
    th_bounds <- CIcon_TNbounds(theta, matrix(deltahat, nrow = 1),
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
    
    th_bounds <- CIcon_TNbounds(theta, matrix(deltahat, nrow = 1),
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

CIcon_TNbounds <- function(theta, deltahat,
                           Al, Au,
                           sigma_l, sigma_u,
                           corr_m, corr_l, corr_u,
                           cLFl, cLFu, tol_r, g) {
  
  B   <- nrow(deltahat)
  k_l <- nrow(Al); k_u <- nrow(Au)
  
  lambdahat_l <- g(deltahat); lambdahat_u <- lambdahat_l
  
  Tlb <- sweep(lambdahat_l - theta, 2, sigma_l, "/")
  Tub <- sweep(theta - lambdahat_u, 2, sigma_u, "/")
  
  Tl <- apply(Tlb, 1, min); bl <- apply(Tlb, 1, which.min)
  Tu <- apply(Tub, 1, min); bu <- apply(Tub, 1, which.min)
  
  corr_m_blb <- corr_m[bl, , drop = FALSE]
  corr_m_bub <- t(corr_m[, bu, drop = FALSE])
  
  tTS_l  <- matrix(1e10, nrow = B, ncol = k_u)
  tTS_ll <- (1 + corr_m_blb)^(-1) * (Tub + corr_m_blb * Tl)
  maskTS_l <- (1 + corr_m_blb) > tol_r
  tTS_l[maskTS_l] <- tTS_ll[maskTS_l]
  
  th_1_1 <- apply(tTS_l, 1, min) * (Tl >= Tu)
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
  
  th_2_1 <- apply(cbind(tLFl, tB_l2), 1, min) * (Tl >= Tu)
  
  tTS_u  <- matrix(1e10, nrow = B, ncol = k_l)
  tTS_u1 <- (1 + corr_m_bub)^(-1) * (Tlb + corr_m_bub * Tu)
  maskTS_u <- (1 + corr_m_bub) > tol_r
  tTS_u[maskTS_u] <- tTS_u1[maskTS_u]
  
  th_1_2 <- apply(tTS_u, 1, min) * (Tl < Tu)
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
  
  th_2_2 <- apply(cbind(tLFu, tB_u2), 1, min) * (Tl < Tu)
  
  th_1 <- th_1_1 + th_1_2
  th_2 <- th_2_1 + th_2_2
  list(th_1 = th_1, th_2 = th_2)
}

pvalue_mcub_zero_fast <- function(deltahat, deltaSigma, Al, Au, g,
                                  eta          = 0.001,
                                  B_fast       = 600,
                                  Blarge_fast  = 6000,
                                  tol          = 1e-3,
                                  tol_r        = 1e-3,
                                  alpha_grid   = c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3,
                                                   0.01, 0.02, 0.05, 0.10, 0.20, 0.40, 0.70),
                                  refine_steps = 6,
                                  seed         = 123) {
  
  inside0 <- function(alpha) {
    set.seed(seed + as.integer(alpha * 1e6))
    alphac <- 0.8 * alpha
    
    res <- CIhybrid(deltahat, deltaSigma, Al, Au,
                    alpha  = alpha,
                    alphac = alphac,
                    eta    = eta,
                    B      = B_fast,
                    Blarge = Blarge_fast,
                    tol    = tol,
                    tol_r  = tol_r,
                    index  = NULL,
                    g      = g,
                    seed   = seed + as.integer(alpha * 1e6))
    CI <- res$CI_h
    (CI[1] <= 0 && 0 <= CI[2])
  }
  
  ins <- vapply(alpha_grid, inside0, logical(1))
  
  if (!ins[1])              return(alpha_grid[1])
  if (ins[length(ins)])     return(1.0)
  
  j  <- which(!ins)[1]
  lo <- alpha_grid[j - 1]; hi <- alpha_grid[j]
  
  for (k in seq_len(refine_steps)) {
    mid <- 0.5 * (lo + hi)
    if (inside0(mid)) lo <- mid else hi <- mid
  }
  
  hi
}

#================================================================
# MAIN CHECK COMPATIBILITY FUNCTION
#================================================================

check_compatibility_simple <- function(df,
                                       i,
                                       alpha     = 0.05,
                                       rxu_range = c(0, 0.8)) {
  
  if (!all(c("x", "y", "z") %in% names(df))) {
    stop("DataFrame must contain columns 'x', 'y', 'z'")
  }
  
  df <- df[complete.cases(df[, c("x", "y", "z")]), ]
  n  <- nrow(df)
  
  if (n < 30) {
    warning(paste("Small sample size:", n, "observations"))
    return(NULL)
  }
  
  rho_xz   <- cor(df$x, df$z, use = "complete.obs")
  rho_yz   <- cor(df$y, df$z, use = "complete.obs")
  rho_xy   <- cor(df$y, df$x, use = "complete.obs")
  deltahat <- c(rho_xy, rho_xz, rho_yz)
  
  deltaSigma <- estimate_cov_corr_boot(df$x, df$y, df$z, B = 800, seed = 1000 + i)
  
  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = 50)
  g      <- function(delta) g_xu_safe(r_grid, delta)
  
  B      <- 500; Blarge <- B * 10; eta <- 0.001
  alphac <- 0.8 * alpha; tol <- 1e-3; tol_r <- 1e-3
  
  A  <- t(sapply(r_grid, function(r_xu) {
    as.numeric(local_compute_gradient_safe(r_xu, deltahat[1], deltahat[2], deltahat[3]))
  }))
  Al <- A; Au <- A
  
  tryCatch({
    cat(sprintf("  Instrument %d: Computing simple CI...\n", i))
    
    res_simple <- ci_simple_union(df$x, df$y, df$z,
                                  rxu_range  = rxu_range,
                                  alpha      = alpha,
                                  cov_method = "bootstrap")
    
    CI_s             <- res_simple$CI
    plug_in          <- res_simple$plug_in
    contains_zero_s  <- (CI_s[1] <= 0 & CI_s[2] >= 0)
    
    cat(sprintf("  Instrument %d: Simple CI computed successfully\n", i))
    
    CI_b            <- CI_s
    contains_zero_b <- contains_zero_s
    p_zero          <- ifelse(contains_zero_b, 1, 0.01)
    
    tryCatch({
      cat(sprintf("  Instrument %d: Computing MCUB CI...\n", i))
      
      res_bei <- CIhybrid(deltahat, deltaSigma, Al, Au,
                          alpha  = alpha,
                          alphac = alphac,
                          eta    = eta,
                          B      = B,
                          Blarge = Blarge,
                          tol    = tol,
                          tol_r  = tol_r,
                          index  = NULL,
                          g      = g,
                          seed   = 2000 + i)
      
      CI_b            <- res_bei$CI_h
      contains_zero_b <- (CI_b[1] <= 0 & CI_b[2] >= 0)
      
      p_zero <- pvalue_mcub_zero_fast(deltahat, deltaSigma, Al, Au, g,
                                      eta = eta, tol = tol, tol_r = tol_r,
                                      B_fast = 300, Blarge_fast = 3000,
                                      seed = 3000 + i)
      
      cat(sprintf("  Instrument %d: MCUB computed successfully\n", i))
      
    }, error = function(e) {
      cat(sprintf("  Instrument %d: MCUB failed, using simple CI: %s\n", i, e$message))
    })
    
    data.frame(
      Z           = paste0("Z", i),
      r_xz        = round(rho_xz, 3),
      r_zu        = if ("u" %in% names(df)) round(cor(df$z, df$u), 3) else NA,
      beta_IV     = round(cov(df$z, df$y, use = "complete.obs") /
                            cov(df$z, df$x, use = "complete.obs"), 3),
      plug_in     = sprintf("[%.3f, %.3f]", plug_in[1], plug_in[2]),
      CI_Bei      = sprintf("[%.3f, %.3f]", CI_b[1], CI_b[2]),
      CI_simple   = sprintf("[%.3f, %.3f]", max(-1, CI_s[1]), min(1, CI_s[2])),
      Zero_in_CI  = ifelse(contains_zero_b, "\u2713", "\u00d7"),
      p_zero      = round(p_zero, 3),
      n           = n,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    cat(sprintf("  Instrument %d: Complete failure: %s\n", i, e$message))
    return(NULL)
  })
}

#================================================================
# Standardize output format
#================================================================

standardize_result <- function(res, study_name, instrument_name) {
  if (is.null(res)) return(NULL)
  
  standard_cols <- c("Z", "r_xz", "r_zu", "beta_IV", "plug_in",
                     "CI_Bei", "CI_simple", "Zero_in_CI", "p_zero",
                     "n", "Instrument", "Study")
  
  for (col in standard_cols) {
    if (!col %in% names(res)) res[[col]] <- NA
  }
  
  res$Study      <- study_name
  res$Instrument <- instrument_name
  
  res[, standard_cols]
}

#================================================================
# APPLICATION FUNCTIONS
#================================================================

run_ditella_analysis_fixed <- function(instruments = c("judgeAlreadyUsedEM", "percJudgeSentToEM"),
                                       alpha = 0.05,
                                       rxu_range = c(0, 0.25),
                                       seed = 123L,
                                       i_start = 1L,
                                       k = 1) {
  cat("\n=== Di Tella & Schargrodsky Analysis ===\n")
  
  mydata <- read_dta("https://github.com/ratbekd/Orientation_paper/raw/refs/heads/main/JPE%20-%20Di%20Tella%20and%20Schargrodsky%20-%20CriminalRecidivismAfterPrisonAndElectronicMonitoring%20.dta")
  data   <- mydata %>% filter(offendersPerJudge > 9)
  
  controls <- c("mostSeriousCrime", "age", "ageSquared", "argentine",
                "numberPreviousImprisonments", "judicialDistrict", "yearOfImprisonment")
  Y <- "recidivism"; X <- "electronicMonitoring"
  
  res_list <- vector("list", length(instruments))
  
  for (j in seq_along(instruments)) {
    Z    <- instruments[j]
    vars <- c(Y, X, Z, controls)
    d    <- data[stats::complete.cases(data[, vars]), vars, drop = FALSE]
    
    fml_y <- stats::as.formula(paste(Y, "~", paste(controls, collapse = "+")))
    fml_x <- stats::as.formula(paste(X, "~", paste(controls, collapse = "+")))
    fml_z <- stats::as.formula(paste(Z, "~", paste(controls, collapse = "+")))
    
    y  <- stats::resid(stats::lm(fml_y, data = d))
    x  <- stats::resid(stats::lm(fml_x, data = d))
    z0 <- stats::resid(stats::lm(fml_z, data = d))
    z  <- stats::predict(stats::lm(z0 ~ x + y))
    
    df  <- data.frame(x = x, y = y, z = z)
    out <- check_compatibility_simple(df,
                                      i         = i_start + j - 1L,
                                      alpha     = alpha,
                                      rxu_range = rxu_range)
    
    if (!is.null(out)) {
      out$Instrument <- Z
      out$n_sample   <- nrow(df)
    }
    res_list[[j]] <- out
  }
  
  res <- do.call(rbind, res_list)
  rownames(res) <- NULL
  res
}

#================================================================
# AUTO-RUN GUARD
#================================================================

.dial_autorun <- !exists("DIAL_SKIP_AUTORUN") || !isTRUE(DIAL_SKIP_AUTORUN)

if (.dial_autorun) {
  cat("\n=== Starting Analysis ===\n")
  result <- run_ditella_analysis_fixed()
  if (!is.null(result)) { cat("\n=== RESULTS ===\n"); print(t(result)) } else {
    cat("\nNo results generated\n")
  }
} else {
  cat("  [Fin_Empirical5.R] DIAL_SKIP_AUTORUN = TRUE — skipping Di Tella auto-run.\n")
}

#================================================================
# BURDE & LINDEN — canonical c(0.0, 0.4)   [NO CHANGE]
#================================================================

run_burde_analysis_fixed <- function(rxu_range = c(0.0, 0.4)) {
  cat("\n=== Burde & Linden Analysis (canonical: rxu_range = c(0.0, 0.4)) ===\n")
  
  data <- read_dta("https://github.com/ratbekd/Orientation_paper/raw/refs/heads/main/afgan.dta")
  
  controls <- c("headchild", "age", "yrsvill", "farsi", "tajik", "farmers",
                "agehead", "educhead", "nhh", "land", "sheep", "distschool", "chagcharan")
  Y <- "testscore"; X <- "enrolled"; Z <- "buildschool"
  vars <- c(Y, X, Z, controls)
  d    <- data[stats::complete.cases(data[, vars]), vars]
  
  y  <- stats::resid(stats::lm(stats::as.formula(paste(Y, "~", paste(controls, collapse = "+"))), data = d))
  x  <- stats::resid(stats::lm(stats::as.formula(paste(X, "~", paste(controls, collapse = "+"))), data = d))
  z0 <- stats::resid(stats::lm(stats::as.formula(paste(Z, "~", paste(controls, collapse = "+"))), data = d))
  z  <- stats::predict(stats::lm(z0 ~ x + y))
  
  df  <- data.frame(x = x, y = y, z = z)
  out <- check_compatibility_simple(df, i = 1L, alpha = 0.05, rxu_range = rxu_range)
  
  if (!is.null(out)) { out$Instrument <- Z; out$n_sample <- nrow(df) }
  out
}

#================================================================
# BURDE & LINDEN — sensitivity c(-0.4, 0.0)   [P1-BURDE NEW]
#
# Rationale (DIAL_corrections_todo.md, Task 3):
#   The canonical c(0.0, 0.4) yields Invalid because the identified
#   set sits entirely above zero.  The structural question is whether
#   the school-construction endogeneity is positive (villages that
#   receive a school were already trending up on test scores) or
#   negative (schools built in lagging areas).  Running c(-0.4, 0)
#   opens the negative half; if the set then contains zero the
#   instrument is valid under that structural story.
#   Compare both rows in phase1_results.csv before updating Table 6.
#================================================================

run_burde_sensitivity_fixed <- function() {
  cat("\n=== Burde & Linden Analysis (sensitivity: rxu_range = c(-0.4, 0.0)) ===\n")
  cat("  NOTE: This is the Phase 1 sign-sensitivity check.\n")
  cat("  Compare with run_burde_analysis_fixed() (canonical c(0,0.4)).\n")
  
  data <- read_dta("https://github.com/ratbekd/Orientation_paper/raw/refs/heads/main/afgan.dta")
  
  controls <- c("headchild", "age", "yrsvill", "farsi", "tajik", "farmers",
                "agehead", "educhead", "nhh", "land", "sheep", "distschool", "chagcharan")
  Y <- "testscore"; X <- "enrolled"; Z <- "buildschool"
  vars <- c(Y, X, Z, controls)
  d    <- data[stats::complete.cases(data[, vars]), vars]
  
  y  <- stats::resid(stats::lm(stats::as.formula(paste(Y, "~", paste(controls, collapse = "+"))), data = d))
  x  <- stats::resid(stats::lm(stats::as.formula(paste(X, "~", paste(controls, collapse = "+"))), data = d))
  z0 <- stats::resid(stats::lm(stats::as.formula(paste(Z, "~", paste(controls, collapse = "+"))), data = d))
  z  <- stats::predict(stats::lm(z0 ~ x + y))
  
  df  <- data.frame(x = x, y = y, z = z)
  # [P1-BURDE] Sensitivity domain: negative half of the real line
  out <- check_compatibility_simple(df, i = 2L, alpha = 0.05, rxu_range = c(-0.4, 0.0))
  
  if (!is.null(out)) {
    out$Instrument  <- Z
    out$n_sample    <- nrow(df)
    out$rxu_range   <- "c(-0.4, 0.0)"   # explicit record for Table 6
  }
  out
}

if (.dial_autorun) {
  result2 <- run_burde_analysis_fixed()
  if (!is.null(result2)) { cat("\n=== RESULTS ===\n"); print(t(result2)) } else {
    cat("\nNo results generated\n")
  }
} else {
  cat("  [Fin_Empirical5.R] DIAL_SKIP_AUTORUN = TRUE — skipping Burde auto-run.\n")
}

#================================================================
# GALIANI — canonical c(-0.2, 0)   [NO CHANGE, VERIFY ONLY]
#================================================================

run_galiani_analysis_fixed <- function(rxu_range = c(-0.2, 0)) {
  cat("\n=== Galiani et al. Analysis ===\n")
  
  mydata <- read_dta("https://github.com/ratbekd/Orientation_paper/raw/refs/heads/main/Crime.dta")
  data   <- mydata %>% filter(cohort > 1957 & cohort < 1963)
  
  controls <- c("cohort", "draftnumber", "navy")
  Y <- "crimerate"; X <- "sm"; Z <- "highnumber"
  vars <- c(Y, X, Z, controls)
  d    <- data[stats::complete.cases(data[, vars]), vars]
  
  y  <- stats::resid(stats::lm(stats::as.formula(paste(Y, "~", paste(controls, collapse = "+"))), data = d))
  x  <- stats::resid(stats::lm(stats::as.formula(paste(X, "~", paste(controls, collapse = "+"))), data = d))
  z0 <- stats::resid(stats::lm(stats::as.formula(paste(Z, "~", paste(controls, collapse = "+"))), data = d))
  z  <- stats::predict(stats::lm(z0 ~ x + y))
  
  df  <- data.frame(x = x, y = y, z = z)
  out <- check_compatibility_simple(df, i = 1L, alpha = 0.05, rxu_range = rxu_range)
  
  if (!is.null(out)) { out$Instrument <- Z; out$n_sample <- nrow(df) }
  out
}

if (.dial_autorun) {
  result3 <- run_galiani_analysis_fixed()
  if (!is.null(result3)) { cat("\n=== RESULTS ===\n"); print(t(result3)) } else {
    cat("\nNo results generated\n")
  }
} else {
  cat("  [Fin_Empirical5.R] DIAL_SKIP_AUTORUN = TRUE — skipping Galiani auto-run.\n")
}

#================================================================
# BANERJEE — canonical c(-0.6, -0.4)   [P1-BANERJEE CORRECTED]
#
# CHANGE LOG:
#   Old value: c(-0.8, -0.6)   — wrong sign region, caused C verdict
#   New value: c(-0.6, -0.4)   — aligns with negative-rho_DU structural prior
#
# OPEN QUESTION (resolve via DIAL_Phase1_runner.R):
#   PDF page-15 yellow comment says "it was (-0.6, -0.5) for the
#   canonical case". Run DIAL_Phase1_runner.R to compare c(-0.6,-0.4)
#   vs c(-0.6,-0.5) and record which gives the more defensible set
#   before updating Table 6 "Canonical D" and "CI MCUB (canonical)".
#================================================================

run_banerjee_analysis_fixed <- function(rxu_range = c(-0.6, -0.4)) {
  cat("\n=== Banerjee et al. Analysis ===\n")
  cat("  [P1-BANERJEE] rxu_range = c(-0.6, -0.4)  (corrected from c(-0.8,-0.6))\n")
  
  mydata <- read_dta("https://github.com/ratbekd/Orientation_paper/raw/refs/heads/main/yld_sett_aug03.dta")
  data   <- mydata %>% filter(year >= 1965 & phwht <= 1)
  
  controls <- c("alt", "totrain", "so_black", "so_red", "so_all",
                "lat", "coastal", "brule1", "year")
  Y <- "phwht"; X <- "p_nland"; Z <- "instru"
  vars <- c(Y, X, Z, controls)
  d    <- data[stats::complete.cases(data[, vars]), vars]
  
  y  <- stats::resid(stats::lm(stats::as.formula(paste(Y, "~", paste(controls, collapse = "+"))), data = d))
  x  <- stats::resid(stats::lm(stats::as.formula(paste(X, "~", paste(controls, collapse = "+"))), data = d))
  z0 <- stats::resid(stats::lm(stats::as.formula(paste(Z, "~", paste(controls, collapse = "+"))), data = d))
  z  <- stats::predict(stats::lm(z0 ~ x + y))
  
  df  <- data.frame(x = x, y = y, z = z)
  # [P1-BANERJEE] corrected domain — negative half, tighter subinterval
  out <- check_compatibility_simple(df, i = 1L, alpha = 0.05, rxu_range = rxu_range)
  
  if (!is.null(out)) { out$Instrument <- Z; out$n_sample <- nrow(df) }
  out
}

if (.dial_autorun) {
  result4 <- run_banerjee_analysis_fixed()
  if (!is.null(result4)) { cat("\n=== RESULTS ===\n"); print(t(result4)) } else {
    cat("\nNo results generated\n")
  }
} else {
  cat("  [Fin_Empirical5.R] DIAL_SKIP_AUTORUN = TRUE — skipping Banerjee auto-run.\n")
}

#================================================================
# RUN ALL APPLICATIONS (canonical domains only)
#================================================================

run_all_applications <- function() {
  all_results <- list()
  
  tryCatch({ res <- run_ditella_analysis_fixed()
  if (!is.null(res)) { res$Study <- "Di Tella & Schargrodsky"; all_results[[1]] <- res } },
  error = function(e) message("Error in Di Tella analysis: ", conditionMessage(e)))
  
  tryCatch({ res <- run_burde_analysis_fixed()
  if (!is.null(res)) { res$Study <- "Burde & Linden"; all_results[[2]] <- res } },
  error = function(e) message("Error in Burde analysis: ", conditionMessage(e)))
  
  tryCatch({ res <- run_galiani_analysis_fixed()
  if (!is.null(res)) { res$Study <- "Galiani et al."; all_results[[3]] <- res } },
  error = function(e) message("Error in Galiani analysis: ", conditionMessage(e)))
  
  tryCatch({ res <- run_banerjee_analysis_fixed()
  if (!is.null(res)) { res$Study <- "Banerjee et al."; all_results[[4]] <- res } },
  error = function(e) message("Error in Banerjee analysis: ", conditionMessage(e)))
  
  all_results <- all_results[!vapply(all_results, is.null, logical(1))]
  if (!length(all_results)) return(NULL)
  
  combined           <- dplyr::bind_rows(all_results)
  combined           <- combined[, c("Study", setdiff(names(combined), "Study")), drop = FALSE]
  rownames(combined) <- NULL
  combined
}

if (.dial_autorun) {
  results <- run_all_applications()
  if (!is.null(results)) { cat("\n=== RESULTS ===\n"); print(results) }
} else {
  cat("  [Fin_Empirical5.R] DIAL_SKIP_AUTORUN = TRUE — skipping run_all_applications.\n")
}
