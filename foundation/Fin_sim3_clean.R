require("nlme")
require("MASS")
# require("ggplot2")
# library("ggplot2")
require("ivreg")
require("lmtest")
library(stargazer)
library(e1071)
library(dplyr)
library(fGarch)
library(lmtest)
library(zoo)
library(sandwich)
library("gridExtra")
library(gridExtra)
library(grid)
library(knitr)
library('MASS')

fs = function(x,epsilon,delta) dnorm(sinh(delta*asinh(x)-epsilon))*delta*cosh(delta*asinh(x)-epsilon)/sqrt(1+x^2)
seed=123
N=800
set.seed(seed)
b1=15
ds <- 2*b1/N
vec = seq(-b1,b1,ds)
vec <- vec[0:N]
x <- fs(vec, 0, 0.9) * 5 + rsnorm(N, mean = 1, sd = 0.25, xi = 1.5)

u1 <- runif(N, min = -1, max = 1) + rsnorm(N, mean = 1, sd = 1, xi = 1.2)
x1 <- (x - mean(x)) / sd(x)
v <- resid(lm(u1 ~ x1))
v <- v * mean(x) / 2 * sd(x) / sd(v) * 2

e <- (x - mean(x)) * 2 + rnorm(N, sd = sd(x)) * 3
### Computing the total error term  depending on
u<-e+v# for the cov(x,e)>0 case

s=1
z=0
# population regression
beta <-2
y<-beta*x+u
j=0
R=0



mydata <- data.frame(x,u,y,e,v,R,z)  #

y0<-y
x0<-x
V=0
### Generating a vector orthogonal to x
fity<-lm(y0~(x0), data=mydata)
V<-resid(fity)
V<-(V-mean(V))/sd(V)
V<-V*sd(x0)
R<-V
mydata$R<-V
k=1


d=0.616
set.seed(1006)

proxy_u <- resid(lm(y ~ x, data=mydata))  # Use residuals as proxy for u
A0_raw <- rnorm(N)*0.75
mydata$A0 <- resid(lm(A0_raw ~ proxy_u))  # Orthogonalize against proxy
mydata$z <- (mydata$x - k*d*mydata$R) + mydata$A0


cor(mydata$z,mydata$x)
cor(mydata$z,mydata$u)
cor(mydata$x,mydata$u)
cor(mydata$y,mydata$z)
cor(mydata$y,mydata$x)

ivreg(mydata$y~mydata$x|mydata$z,data=mydata)

rho_xz <- cor(mydata$z,mydata$x)
rho_xu <-cor(mydata$x,mydata$u)
rho_xy <-cor(mydata$y,mydata$x)
rho_yz <- cor(mydata$y,mydata$z)
rho_zu <- cor(mydata$z,mydata$u)



df <- mydata

#================================================================
# GENERALIZED SIMULATION FRAMEWORK
#================================================================




CIhybrid <- function(deltahat, deltaSigma, Al, Au,
                     alpha, alphac, eta,
                     B, Blarge,
                     tol, tol_r, index, g,
                     seed = 123) {   # [FIX #11] explicit seed parameter
  
  set.seed(seed)   # [FIX #11] use argument
  
  # ----- 1. Random draws and basic quantities -----
  library(MASS)
  library(stats)
  
  k_dim <- length(deltahat)
  deltastar_demean_large <- MASS::mvrnorm(Blarge, mu = rep(0, k_dim), Sigma = deltaSigma)
  deviation <- apply(abs(deltastar_demean_large / sqrt(diag(deltaSigma))), 1, max)
  c_bd <- quantile(deviation, 1 - eta)
  
  sigma_delta <- sqrt(diag(deltaSigma))
  
  # ----- 2. Construct lambdaSigma and related components -----
  Lambda <- rbind(Al, Au)
  lambdaSigma <- Lambda %*% deltaSigma %*% t(Lambda)
  lambdasigma <- sqrt(diag(lambdaSigma))
  
  kk <- nrow(Al)
  sigma_l <- lambdasigma[1:kk]
  sigma_u <- lambdasigma[(kk + 1):(2 * kk)]
  
  corr_all <- diag(1 / lambdasigma) %*% lambdaSigma %*% diag(1 / lambdasigma)
  corr_l <- corr_all[1:kk, 1:kk]
  corr_u <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]
  corr_m <- corr_all[1:kk, (kk + 1):(2 * kk)]
  
  # ----- 3. Feasible bounds -----
  lb <- deltahat - sigma_delta * c_bd
  ub <- deltahat + sigma_delta * c_bd
  delta1 <- deltahat
  
  lb[index] <- 0
  ub[index] <- 0
  delta1[index] <- 0
  
  cl <- 0
  cu <- qnorm(1 - alpha / 2)
  c <- (cl + cu) / 2
  
  delta_fea <- list()
  c_fea <- numeric()
  
  # ----- 4. Objective wrapper -----
  obj_large <- function(delta, c_check) {
    (alpha - CIproj_p(c_check, delta, alphac, c_bd,
                      deltastar_demean_large,
                      Al, Au, sigma_l, sigma_u,
                      deltaSigma, corr_m, corr_l, corr_u,
                      eta, tol_r, g)) * 100
  }
  
  # ----- 5. Iterative bisection loop -----
  k <- 1
  while ((cu - cl) > tol) {
    set.seed(k)
    k <- k + 1
    
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
      # Optimizer equivalent of fmincon + GlobalSearch
      f_optim <- tryCatch({
        optim(par = delta1, fn = obj, method = "L-BFGS-B",
              lower = lb, upper = ub, control = list(maxit = 1000))
      }, error = function(e) NULL)
      
      if (!is.null(f_optim)) {
        delta1 <- f_optim$par
      }
      p1 <- obj_large(delta1, c)
      
      if (p1 >= 0) {
        cu <- c
      } else {
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
  
  # ----- 6. Final refinement -----
  cl <- c
  cu <- qnorm(1 - alpha / 2)
  
  while ((cu - cl) > tol) {
    c <- (cl + cu) / 2
    p <- obj_large(delta1, c)
    if (p >= 0) {
      cu <- c
    } else {
      cl <- c
    }
  }
  
  # ----- 7. Compute confidence intervals -----
  lambdahat_l <- g(deltahat)
  lambdahat_u <- lambdahat_l
  
  c_LF <- qnorm(1 - eta)
  CI_p <- c(min(lambdahat_l - c * sigma_l), max(lambdahat_u + c * sigma_u))
  CI_c <- CIcon(deltahat, deltaSigma, Al, Au, c_LF, alphac, tol, tol_r, g)
  
  CI_h <- c(min(CI_p[1], CI_c[1]), max(CI_p[2], CI_c[2]))
  
  list(CI_h = CI_h,
       CI_c = CI_c,
       CI_p = CI_p,
       delta1 = delta1,
       c = c)
}

#### helper functions
CIproj_p <- function(c, delta, alphac, c_bd,
                     deltastar_demean,
                     Al, Au,
                     sigma_l, sigma_u,
                     deltaSigma,
                     corr_m, corr_l, corr_u,
                     eta, tol_r, g) {
  # -------------------------------
  # Equivalent to MATLAB CIproj_p.m
  # -------------------------------
  # g : function(delta_matrix) -> lambda_matrix
  # delta, deltastar_demean: row or column vectors/matrices
  # Al, Au: constraint matrices (lower/upper)
  
  # ----- Lambda at point delta -----
  lambda_l <- g(matrix(delta, nrow = 1))
  lambda_u <- g(matrix(delta, nrow = 1))
  lb <- min(lambda_l)
  ub <- max(lambda_u)
  mb <- (lb + ub) / 2
  
  # ----- Randomized samples -----
  deltastar <- sweep(deltastar_demean, 2, delta, "+")
  lambdastar_l <- g(deltastar)
  lambdastar_u <- g(deltastar)
  
  # ----- Condition Δ -----
  delta_sigma <- sqrt(diag(deltaSigma))
  maxdev <- apply(abs(deltastar_demean / matrix(delta_sigma, nrow = nrow(deltastar_demean), ncol = length(delta_sigma), byrow = TRUE)), 1, max)
  ind_Delta <- as.numeric(maxdev <= c_bd)
  
  # ----- Projection statistics -----
  Tstar_l <- pmax(
    apply((lambdastar_l - lb) / matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE), 1, min),
    apply((lb - lambdastar_u) / matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
  )
  Tstar_m <- pmax(
    apply((lambdastar_l - mb) / matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE), 1, min),
    apply((mb - lambdastar_u) / matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
  )
  Tstar_u <- pmax(
    apply((lambdastar_l - ub) / matrix(sigma_l, nrow = nrow(lambdastar_l), ncol = length(sigma_l), byrow = TRUE), 1, min),
    apply((ub - lambdastar_u) / matrix(sigma_u, nrow = nrow(lambdastar_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
  )
  
  ind_proj_l <- as.numeric(Tstar_l <= c)
  ind_proj_m <- as.numeric(Tstar_m <= c)
  ind_proj_u <- as.numeric(Tstar_u <= c)
  
  # ----- Condition cLF -----
  cLF <- qnorm(1 - eta)
  
  # === Helper function for conditional parts ===
  cond_block <- function(bound, select, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, label) {
    if (sum(select) == 0) {
      return(rep(0, nrow(deltastar)))
    }
    sub_deltastar <- deltastar[select == 1, , drop = FALSE]
    sub_l_l <- lambdastar_l[select == 1, , drop = FALSE]
    sub_l_u <- lambdastar_u[select == 1, , drop = FALSE]
    
    Tc1 <- apply((sub_l_l - bound) / matrix(sigma_l, nrow = nrow(sub_l_l), ncol = length(sigma_l), byrow = TRUE), 1, min)
    Tc2 <- apply((bound - sub_l_u) / matrix(sigma_u, nrow = nrow(sub_l_u), ncol = length(sigma_u), byrow = TRUE), 1, min)
    Tc <- pmax(Tc1, Tc2)
    
    # Call to the truncated-normal bounds function
    th_bounds <- CIcon_TNbounds(bound, sub_deltastar, Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u, cLF, cLF, tol_r, g)
    th1 <- th_bounds[[1]]; th2 <- th_bounds[[2]]
    
    phi <- (pnorm(Tc) - pnorm(th1)) / (pnorm(th2) - pnorm(th1))
    sub_ind <- as.numeric(phi <= 1 - alphac)
    
    ind_c <- rep(0, nrow(deltastar))
    ind_c[select == 1] <- sub_ind
    ind_c
  }
  
  # ----- l / m / u conditions -----
  selectl <- (1 - ind_proj_l) * ind_Delta
  ind_c_l <- cond_block(lb, selectl, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, "l")
  
  selectm <- (1 - ind_proj_m) * ind_Delta
  ind_c_m <- cond_block(mb, selectm, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, "m")
  
  selectu <- (1 - ind_proj_u) * ind_Delta
  ind_c_u <- cond_block(ub, selectu, lambdastar_l, lambdastar_u, deltastar, sigma_l, sigma_u, "u")
  
  # ----- Combine indicators -----
  ind_l <- 1 - (1 - ind_c_l) * (1 - ind_proj_l)
  ind_m <- 1 - (1 - ind_c_m) * (1 - ind_proj_m)
  ind_u <- 1 - (1 - ind_c_u) * (1 - ind_proj_u)
  
  # p1 <- mean((1 - ind_l * ind_m) * ind_Delta)
  # p2 <- mean((1 - ind_u * ind_m) * ind_Delta)
  # p <- max(p1, p2) + eta
  p1  <- mean((1 - ind_l * ind_m) * ind_Delta, na.rm = TRUE)
  p2  <- mean((1 - ind_u * ind_m) * ind_Delta, na.rm = TRUE)
  p   <- max(p1, p2, na.rm = TRUE) + eta
  
  if (!is.finite(p)) p <- NA_real_
  return(p)
  
  return(p)
}


CIcon <- function(deltahat, deltaSigma, Al, Au,
                  c_LF, alphac, tol, tol_r, g) {
  # -------------------------------
  # Equivalent to MATLAB CIcon.m
  # -------------------------------
  
  kk <- nrow(Al)
  
  # ----- 1. Lambda covariance and correlations -----
  Lambda <- rbind(Al, Au)
  lambdaSigma <- Lambda %*% deltaSigma %*% t(Lambda)
  lambda_sigma <- sqrt(diag(lambdaSigma))
  
  corr_all <- diag(1 / lambda_sigma) %*% lambdaSigma %*% diag(1 / lambda_sigma)
  corr_m <- corr_all[1:kk, (kk + 1):(2 * kk)]
  corr_l <- corr_all[1:kk, 1:kk]
  corr_u <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]
  
  sigma_l <- lambda_sigma[1:kk]
  sigma_u <- lambda_sigma[(kk + 1):(2 * kk)]
  
  # ----- 2. Compute lambda bounds -----
  gval <- as.numeric(g(matrix(deltahat, nrow = 1)))
  
  lb <- gval - sigma_l * c_LF
  ub <- gval + sigma_u * c_LF
  lb <- min(lb)
  ub <- max(ub)
  
  mid <- (lb + ub) / 2
  lb_pt <- min(gval)
  ub_pt <- max(gval)
  
  # ----- 3. CI lower bound -----
  rej <- 1
  theta <- lb
  
  while (rej == 1 && theta <= min(mid, lb_pt)) {
    Tcl <- min((gval - theta) / sigma_l)
    Tcu <- min((theta - gval) / sigma_u)
    Tc <- max(Tcl, Tcu)
    
    th_bounds <- CIcon_TNbounds(theta, matrix(deltahat, nrow = 1),
                                Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u,
                                c_LF, c_LF, tol_r, g)
    th_1 <- th_bounds[[1]]
    th_2 <- th_bounds[[2]]
    
    t <- qnorm((1 - alphac) * pnorm(th_2) + alphac * pnorm(th_1))
    rej <- as.numeric(Tc > t)
    theta <- theta + tol
  }
  CI_lower <- theta
  
  # ----- 4. CI upper bound -----
  rej <- 1
  theta <- ub
  
  while (rej == 1 && theta >= max(mid, ub_pt)) {
    Tcl <- min((gval - theta) / sigma_l)
    Tcu <- min((theta - gval) / sigma_u)
    Tc <- max(Tcl, Tcu)
    
    th_bounds <- CIcon_TNbounds(theta, matrix(deltahat, nrow = 1),
                                Al, Au, sigma_l, sigma_u,
                                corr_m, corr_l, corr_u,
                                c_LF, c_LF, tol_r, g)
    th_1 <- th_bounds[[1]]
    th_2 <- th_bounds[[2]]
    
    t <- qnorm((1 - alphac) * pnorm(th_2) + alphac * pnorm(th_1))
    rej <- as.numeric(Tc > t)
    theta <- theta - tol
  }
  CI_upper <- theta
  
  # ----- 5. Return final CI -----
  CIcon <- c(CI_lower, CI_upper)
  return(CIcon)
}

CIcon_TNbounds <- function(theta, deltahat,
                           Al, Au,
                           sigma_l, sigma_u,
                           corr_m, corr_l, corr_u,
                           cLFl, cLFu, tol_r, g) {
  B   <- nrow(deltahat)
  k_l <- nrow(Al)
  k_u <- nrow(Au)
  
  lambdahat_l <- g(deltahat)
  lambdahat_u <- lambdahat_l
  
  Tlb <- sweep(lambdahat_l - theta, 2, sigma_l, "/")   # B x k_l
  Tub <- sweep(theta - lambdahat_u, 2, sigma_u, "/")   # B x k_u
  
  Tl <- apply(Tlb, 1, min)         # length B
  bl <- apply(Tlb, 1, which.min)   # argmin per row (1..k_l)
  Tu <- apply(Tub, 1, min)
  bu <- apply(Tub, 1, which.min)   # argmin per row (1..k_u)
  
  # corr_m rows/cols at those argmins
  corr_m_blb <- corr_m[bl, , drop = FALSE]           # B x k_u
  corr_m_bub <- t(corr_m[, bu, drop = FALSE])        # B x k_l
  
  ## ----- Case Tl >= Tu : “lower” side -----
  # TS part
  tTS_l  <- matrix(1e10, nrow = B, ncol = k_u)
  # Broadcast Tl and Tub row-wise
  tTS_ll <- (1 + corr_m_blb)^(-1) * (Tub + corr_m_blb * Tl)
  maskTS_l <- (1 + corr_m_blb) > tol_r
  tTS_l[maskTS_l] <- tTS_ll[maskTS_l]
  
  th_1_1 <- apply(tTS_l, 1, min) * (Tl >= Tu)
  th_1_1[th_1_1 == 1e10] <- -1e10
  th_1_1 <- th_1_1 * (Tl >= Tu)
  
  # LF part
  tLFl <- rep(cLFl, B)
  
  # B part  (*** fixed row-wise broadcasting ***)
  # For each i: taux_1[i,j] = (1 - corr_l[bl[i], j])^{-1} * ( Tlb[i, bl[i]] - corr_l[bl[i], j] * Tl[i] )
  CL_bl   <- corr_l[bl, , drop = FALSE]                                 # B x k_l
  Tlb_min <- matrix(Tlb[cbind(1:B, bl)], nrow = B, ncol = k_l)          # B x k_l (rep each row)
  Tl_rep  <- matrix(Tl, nrow = B, ncol = k_l)
  taux_1  <- (1 - CL_bl)^(-1) * (Tlb_min - CL_bl * Tl_rep)              # B x k_l
  
  tB_l2   <- matrix(1e10, nrow = B, ncol = k_l)
  maskB_l <- (1 > CL_bl + tol_r)
  tB_l2[maskB_l] <- taux_1[maskB_l]
  
  th_2_1 <- apply(cbind(tLFl, tB_l2), 1, min) * (Tl >= Tu)
  
  ## ----- Case Tl < Tu : “upper” side -----
  # TS part
  tTS_u  <- matrix(1e10, nrow = B, ncol = k_l)
  tTS_u1 <- (1 + corr_m_bub)^(-1) * (Tlb + corr_m_bub * Tu)
  maskTS_u <- (1 + corr_m_bub) > tol_r
  tTS_u[maskTS_u] <- tTS_u1[maskTS_u]
  
  th_1_2 <- apply(tTS_u, 1, min) * (Tl < Tu)
  th_1_2[th_1_2 == 1e10] <- -1e10
  th_1_2 <- th_1_2 * (Tl < Tu)
  
  # LF part
  tLFu <- rep(cLFu, B)
  
  # B part  (*** fixed row-wise broadcasting ***)
  CU_bu   <- corr_u[bu, , drop = FALSE]                                 # B x k_u
  Tub_min <- matrix(Tub[cbind(1:B, bu)], nrow = B, ncol = k_u)          # B x k_u
  Tu_rep  <- matrix(Tu, nrow = B, ncol = k_u)
  taux_2  <- (1 - CU_bu)^(-1) * (Tub_min - CU_bu * Tu_rep)              # B x k_u
  
  tB_u2   <- matrix(1e10, nrow = B, ncol = k_u)
  maskB_u <- (1 > CU_bu + tol_r)
  tB_u2[maskB_u] <- taux_2[maskB_u]
  
  th_2_2 <- apply(cbind(tLFu, tB_u2), 1, min) * (Tl < Tu)
  
  th_1 <- th_1_1 + th_1_2
  th_2 <- th_2_1 + th_2_2
  list(th_1 = th_1, th_2 = th_2)
}


##vectorized g
g_xu_safe <- function(r_xu, delta, tol = 1e-10) {
  # ---------------------------------------------------------
  # Safe version of g_xu that accepts both vector and matrix delta
  # ---------------------------------------------------------
  if (is.null(dim(delta))) {
    # plain numeric vector of length 3
    delta <- matrix(delta, nrow = 1)
  } else if (ncol(delta) != 3 && nrow(delta) == 3) {
    # 3×1 column -> convert to 1×3 row
    delta <- t(delta)
  } else if (ncol(delta) != 3) {
    stop("delta must have 3 elements: [rho_xy, rho_xz, rho_yz]")
  }
  
  rho_xy <- delta[, 1]
  rho_xz <- delta[, 2]
  rho_yz <- delta[, 3]
  
  B <- nrow(delta)
  K <- length(r_xu)
  
  rho_xy_mat <- matrix(rho_xy, nrow = B, ncol = K)
  rho_xz_mat <- matrix(rho_xz, nrow = B, ncol = K)
  rho_yz_mat <- matrix(rho_yz, nrow = B, ncol = K)
  r_xu_mat   <- matrix(r_xu,   nrow = B, ncol = K, byrow = TRUE)
  
  # Clip safely
  rho_xy_mat <- pmin(pmax(rho_xy_mat, -1 + tol), 1 - tol)
  r_xu_mat   <- pmin(pmax(r_xu_mat,   -1 + tol), 1 - tol)
  
  denom <- sqrt(pmax(1 - rho_xy_mat^2, tol))
  sqrt_term <- sqrt(pmax(1 - r_xu_mat^2, 0))
  
  gval <- rho_xz_mat * r_xu_mat -
    (rho_xy_mat * rho_xz_mat - rho_yz_mat) / denom * sqrt_term
  
  gval[!is.finite(gval)] <- NA_real_
  return(gval)
}


local_compute_gradient_safe <- function(r_xu, rho_xy, rho_xz, rho_yz, tol = 1e-10) {
  # ---------------------------------------------------------
  # Computes gradient of g_xu(r_xu, rho_xy, rho_xz, rho_yz)
  # with NA/Inf safety.
  #
  # Arguments:
  #   r_xu, rho_xy, rho_xz, rho_yz : numeric scalars or vectors
  #   tol : tolerance for numerical stability
  #
  # Returns:
  #   grad : 3 x length(r_xu) matrix (each column is a gradient vector)
  # ---------------------------------------------------------
  
  # Clamp correlations safely inside (-1, 1)
  rho_xy <- pmin(pmax(rho_xy, -1 + tol), 1 - tol)
  r_xu   <- pmin(pmax(r_xu,   -1 + tol), 1 - tol)
  
  # Common sqrt term
  sqrt_term <- sqrt(pmax((1 - r_xu^2) / (1 - rho_xy^2), 0))
  
  # Partial derivatives
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
  n <- nrow(dat)
  
  Rb <- matrix(NA_real_, nrow = B, ncol = 3)
  for (b in 1:B) {
    idx <- sample.int(n, n, replace = TRUE)
    xb <- dat[idx, 1]; yb <- dat[idx, 2]; zb <- dat[idx, 3]
    Rb[b, ] <- c(cor(xb, yb), cor(xb, zb), cor(yb, zb))
  }
  
  S <- cov(Rb, use = "complete.obs")
  S <- (S + t(S)) / 2  # symmetrize
  
  # optional: make PSD in case bootstrap noise makes it slightly indefinite
  if (requireNamespace("Matrix", quietly = TRUE)) {
    S <- as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
  }
  S
}
# --- Helper for simple union (Bonferroni) CI using delta method ---
# ci_simple_union <- function(x, y, z,
#                             D_grid = seq(rxu_range[1], rxu_range[2], length.out = 80),
#                             alpha = 0.05) {
#   n <- length(x)
#   rho_xy <- cor(x, y); rho_xz <- cor(x, z); rho_yz <- cor(y, z)
#
#   # Variances (simplified)
#   var_xy <- (1 - rho_xy^2)^2 / n
#   var_xz <- (1 - rho_xz^2)^2 / n
#   var_yz <- (1 - rho_yz^2)^2 / n
#   Sigma_rho <- diag(c(var_xy, var_xz, var_yz))
#
#   g_vals <- sapply(D_grid, function(r) g_fun(r, rho_xy, rho_xz, rho_yz))
#   idx_l <- which.min(g_vals); idx_u <- which.max(g_vals)
#   plug_in <- c(g_vals[idx_l], g_vals[idx_u])
#
#   # gradient at those r
#   grad_l <- g_grad(D_grid[idx_l], rho_xy, rho_xz, rho_yz)
#   grad_u <- g_grad(D_grid[idx_u], rho_xy, rho_xz, rho_yz)
#
#   var_l <- t(grad_l) %*% Sigma_rho %*% grad_l
#   var_u <- t(grad_u) %*% Sigma_rho %*% grad_u
#
#   K <- length(D_grid)
#   crit <- qnorm(1 - alpha / (2 * K * 2))  # conservative two-sided
#   CI_lower <- plug_in[1] - crit * sqrt(var_l)
#   CI_upper <- plug_in[2] + crit * sqrt(var_u)
#
#   list(CI = c(CI_lower, CI_upper), plug_in = plug_in)
# }
ci_simple_union <- function(x, y, z,
                            rxu_range   = c(0, 0.8),
                            grid_length = 80,
                            alpha       = 0.05,
                            # If you already computed a full 3x3 cov of (rho_xy,rho_xz,rho_yz), pass it here:
                            Sigma_rho   = NULL,
                            cov_method  = c("bootstrap", "diag"),
                            B_boot      = 800,
                            seed        = 123,
                            tol         = 1e-10,
                            make_psd    = TRUE) {
  
  cov_method <- match.arg(cov_method)
  D_grid = seq(rxu_range[1], rxu_range[2], length.out = 80)
  # -----------------------------
  # 0) Clean data
  # -----------------------------
  dat <- data.frame(x = x, y = y, z = z)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  n <- nrow(dat)
  if (n < 30) stop("ci_simple_union: need at least 30 complete observations.")
  
  # -----------------------------
  # 1) Sample correlations (rho_hat)
  # Order: (rho_xy, rho_xz, rho_yz)
  # -----------------------------
  rho_xy <- suppressWarnings(cor(dat$x, dat$y))
  rho_xz <- suppressWarnings(cor(dat$x, dat$z))
  rho_yz <- suppressWarnings(cor(dat$y, dat$z))
  deltahat <- c(rho_xy, rho_xz, rho_yz)
  names(deltahat) <- c("rho_xy", "rho_xz", "rho_yz")
  
  if (any(!is.finite(deltahat))) {
    stop("ci_simple_union: non-finite sample correlations (check for constant variables).")
  }
  
  # Clamp helper (avoid division by zero in sqrt(1 - rho^2))
  clamp <- function(a) pmin(pmax(a, -1 + tol), 1 - tol)
  rho_xy_c <- clamp(rho_xy)
  
  # -----------------------------
  # 2) Full covariance of correlation vector
  # -----------------------------
  if (is.null(Sigma_rho)) {
    if (cov_method == "bootstrap") {
      set.seed(seed)
      Rb <- matrix(NA_real_, nrow = B_boot, ncol = 3)
      for (b in seq_len(B_boot)) {
        idx <- sample.int(n, n, replace = TRUE)
        xb <- dat$x[idx]; yb <- dat$y[idx]; zb <- dat$z[idx]
        Rb[b, ] <- c(
          suppressWarnings(cor(xb, yb)),
          suppressWarnings(cor(xb, zb)),
          suppressWarnings(cor(yb, zb))
        )
      }
      Rb <- Rb[complete.cases(Rb), , drop = FALSE]
      if (nrow(Rb) < 50) stop("ci_simple_union: too many NA bootstrap draws; check data variation.")
      Sigma_rho <- cov(Rb)
      Sigma_rho <- (Sigma_rho + t(Sigma_rho)) / 2
      
      if (make_psd && requireNamespace("Matrix", quietly = TRUE)) {
        Sigma_rho <- as.matrix(Matrix::nearPD(Sigma_rho, corr = FALSE)$mat)
      }
    } else {
      # Fallback: diagonal-only (NOT consistent with "full covariance" claim)
      var_xy <- (1 - rho_xy^2)^2 / n
      var_xz <- (1 - rho_xz^2)^2 / n
      var_yz <- (1 - rho_yz^2)^2 / n
      Sigma_rho <- diag(c(var_xy, var_xz, var_yz))
    }
  } else {
    # sanity checks
    if (!is.matrix(Sigma_rho) || any(dim(Sigma_rho) != c(3, 3))) {
      stop("ci_simple_union: Sigma_rho must be a 3x3 matrix for (rho_xy, rho_xz, rho_yz).")
    }
    Sigma_rho <- (Sigma_rho + t(Sigma_rho)) / 2
  }
  
  # -----------------------------
  # 3) g(r) and gradient wrt (rho_xy, rho_xz, rho_yz)
  # -----------------------------
  g_fun_local <- function(r, rho_xy, rho_xz, rho_yz) {
    rho_xy <- clamp(rho_xy)
    r      <- clamp(r)
    rho_xz * r - (rho_xy * rho_xz - rho_yz) * sqrt((1 - r^2) / (1 - rho_xy^2))
  }
  
  g_grad_local <- function(r, rho_xy, rho_xz, rho_yz) {
    rho_xy <- clamp(rho_xy)
    r      <- clamp(r)
    S <- sqrt((1 - r^2) / (1 - rho_xy^2))
    c(
      # d g / d rho_xy
      -rho_xz * S + (rho_xy * rho_xz - rho_yz) * S * rho_xy / (1 - rho_xy^2),
      # d g / d rho_xz
      r - rho_xy * S,
      # d g / d rho_yz
      S
    )
  }
  
  # -----------------------------
  # 4) Grid search for min/max over r_xu in D
  # -----------------------------
  D_grid <- seq(rxu_range[1], rxu_range[2], length.out = grid_length)
  D_grid <- clamp(D_grid)
  
  g_vals <- vapply(D_grid, g_fun_local, numeric(1),
                   rho_xy = rho_xy, rho_xz = rho_xz, rho_yz = rho_yz)
  
  idx_l <- which.min(g_vals)
  idx_u <- which.max(g_vals)
  
  r_l <- D_grid[idx_l]
  r_u <- D_grid[idx_u]
  
  plug_in <- c(g_vals[idx_l], g_vals[idx_u])   # [min, max]
  
  # -----------------------------
  # 5) Delta-method SE at argmin/argmax + Bonferroni/union-bound critical value
  # -----------------------------
  grad_l <- g_grad_local(r_l, rho_xy, rho_xz, rho_yz)
  grad_u <- g_grad_local(r_u, rho_xy, rho_xz, rho_yz)
  
  var_l <- as.numeric(t(grad_l) %*% Sigma_rho %*% grad_l)
  var_u <- as.numeric(t(grad_u) %*% Sigma_rho %*% grad_u)
  
  var_l <- max(var_l, 0)
  var_u <- max(var_u, 0)
  
  K <- length(D_grid)
  # union bound across r-grid and both endpoints, two-sided:
  crit <- qnorm(1 - alpha / (4 * K))
  
  CI_lower <- plug_in[1] - crit * sqrt(var_l)
  CI_upper <- plug_in[2] + crit * sqrt(var_u)
  
  list(
    CI        = c(CI_lower, CI_upper),
    plug_in   = plug_in,
    r_star    = c(r_l, r_u),
    rho_hat   = deltahat,
    Sigma_rho = Sigma_rho,
    crit      = crit
  )
}

## --- g(r; rho) and its gradient wrt (rho_xy, rho_xz, rho_yz) ----------------

g_fun <- function(r, rho_xy, rho_xz, rho_yz) {
  rho_xz * r - (rho_xy * rho_xz - rho_yz) * sqrt((1 - r^2) / (1 - rho_xy^2))
}
g_grad <- function(r, rho_xy, rho_xz, rho_yz) {
  S <- sqrt((1 - r^2) / (1 - rho_xy^2))
  c(
    -rho_xz * S + (rho_xy * rho_xz - rho_yz) * S * rho_xy / (1 - rho_xy^2), # ∂/∂rho_xy
    r - rho_xy * S,                                                         # ∂/∂rho_xz
    S                                                                       # ∂/∂rho_yz
  )
}
pvalue_mcub_zero_fast <- function(deltahat, deltaSigma, Al, Au, g,
                                  eta = 0.001,
                                  B_fast = 600, Blarge_fast = 6000,
                                  tol = 1e-3, tol_r = 1e-3,
                                  alpha_grid = c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3,
                                                 0.01, 0.02, 0.05, 0.10, 0.20, 0.40, 0.70),
                                  refine_steps = 6,
                                  seed = 123) {
  
  inside0 <- function(alpha) {
    set.seed(seed + as.integer(alpha * 1e6))  # stabilize across calls
    alpha_C <- alpha / 2
    alphac  <- 0.8 * alpha_C
    
    res <- CIhybrid(deltahat, deltaSigma, Al, Au,
                    alpha = alpha_C, alphac = alphac,
                    eta = eta, B = B_fast, Blarge = Blarge_fast,
                    tol = tol, tol_r = tol_r, index = NULL, g = g)
    CI <- res$CI_h
    (CI[1] <= 0 && 0 <= CI[2])
  }
  
  # 1) Grid search
  ins <- vapply(alpha_grid, inside0, logical(1))
  
  # If 0 excluded even at tiny alpha -> p ~ 0
  if (!ins[1]) return(alpha_grid[1])
  
  # If 0 included even at huge alpha -> p ~ 1
  if (ins[length(ins)]) return(1.0)
  
  # Find first alpha where 0 becomes excluded
  j <- which(!ins)[1]
  lo <- alpha_grid[j - 1]  # included
  hi <- alpha_grid[j]      # excluded
  
  # 2) Short refinement
  for (k in seq_len(refine_steps)) {
    mid <- 0.5 * (lo + hi)
    if (inside0(mid)) lo <- mid else hi <- mid
  }
  
  hi
}
#================================================================
# 1. WRAPPER FUNCTION FOR SINGLE SAMPLE ANALYSIS
#================================================================

analyze_single_sample <- function(sample_data, d_vec, rxu_range = rxu_range,
                                  alpha = 0.05, seed = NULL, k = 1) {
  
  if (!is.null(seed)) set.seed(seed)
  
  df <- sample_data
  n <- nrow(df)
  
  # Storage for results
  res_list <- list()
  
  # Grid for r_xu values
  D_grid <- seq(rxu_range[1], rxu_range[2], length.out = 80)
  
  # Loop over instruments
  for (i in seq_along(d_vec)) {
    d <- d_vec[i]
    df$z <- df$x - k*d*df$R + df$A0
    rho_xz<-cor(df$x,df$z,use = "complete.obs",method=c("pearson"))
    rho_yz<-cor(df$y,df$z,use = "complete.obs",method=c("pearson"))
    rho_xy<-cor(df$y,df$x,use = "complete.obs",method=c("pearson"))
    # --- Asymptotic covariance matrix under bivariate normality ---
    var_rho_xy <- (1 - rho_xy^2)^2 / n
    var_rho_xz <- (1 - rho_xz^2)^2 / n
    var_rho_yz <- (1 - rho_yz^2)^2 / n
    #deltaSigma <- diag(c(var_rho_xy, var_rho_xz, var_rho_yz))
    deltaSigma <- estimate_cov_corr_boot(df$x, df$y, df$z, B = 800, seed = 1000 + i)
    # --- Point estimates ---
    deltahat <- c(rho_xy, rho_xz, rho_yz)
    
    # --- Grid of possible r_xu values ---
    r_grid <- seq(0.4, 0.6, length.out = 50)
    
    # --- g(delta) mapping ---
    #g <- function(delta) g_xu_safe(r_grid, matrix(delta, nrow = 1))
    g <- function(delta) g_xu_safe(r_grid, delta)
    
    # --- Tuning parameters for Union Bound ---
    B      <- 500
    Blarge <- B * 10
    eta    <- 0.001
    alphac <- 0.8 * alpha
    tol    <- 1e-3
    tol_r  <- 1e-3
    
    # --- Compute A matrix (Jacobian of g wrt delta) ---
    A <- t(sapply(r_grid, function(r_xu) {
      grad <- local_compute_gradient_safe(r_xu, deltahat[1], deltahat[2], deltahat[3])
      as.numeric(grad)
    }))
    Al <- A
    Au <- A
    
    # res_simple <- ci_simple_union(df$x, df$y, df$z, D_grid = seq(rxu_range[1], rxu_range[2], length.out = 80),
    #                               alpha = 0.05)
    
    Sigma_rho <- estimate_cov_corr_boot(df$x, df$y, df$z, B = 800, seed = 1000 + i)
    res_simple <- ci_simple_union(df$x, df$y, df$z,
                                  rxu_range  = rxu_range,
                                  alpha      = alpha,
                                  cov_method = "bootstrap")
    
    res_bei <-CIhybrid(deltahat, deltaSigma, Al, Au, alpha, alphac,
                       eta, B, Blarge, tol, tol_r, index = NULL, g = g)
    
    CI_b <- res_bei$CI_h
    CI_s <- res_simple$CI
    plug <- res_bei$CI_c
    plug_in <- res_bei$CI_c
    contains_zero_b <- (CI_b[1] <= 0 & CI_b[2] >= 0)
    contains_zero_s <- (CI_s[1] <= 0 & CI_s[2] >= 0)
    
    # coverage test: does CI cover plug-in interval
    cover_b <- (plug[1] >= CI_b[1] & plug[2] <= CI_b[2])
    cover_s <- (plug[1] >= CI_s[1] & plug[2] <= CI_s[2])
    ### === 9. CI for r_zu = 0 ===
    grad_rzu0 <- g_grad(0, rho_xy, rho_xz, rho_yz)
    
    # rho_xu function
    r_xy <- cor(df$x, df$y)
    r_xz <- cor(df$x, df$z)
    r_yz <- cor(df$y, df$z)
    r_xz<-cor(df$x,df$z,use = "complete.obs",method=c("pearson"))
    r_yz<-cor(df$y,df$z,use = "complete.obs",method=c("pearson"))
    r_xy<-cor(df$y,df$x,use = "complete.obs",method=c("pearson"))
    rho_xu <- function(r_xz, r_xy, r_yz) {
      mult <- 1 - r_xy^2
      if (mult <= 1e-10) return(NA)
      denom <- (r_xy * r_xz - r_yz)^2
      if (denom <= 1e-10) return(NA)
      num <- pmax(0, r_xz^2)
      res <- sqrt(1 / (1 + mult * num / denom))
      if (!is.finite(res)) return(NA)
      res
    }
    
    rxu_point <- rho_xu(r_xz, r_xy, r_yz)
    data <- df
    
    # Delta method for rxu_point
    delta_method_rho_xu <- function(data, r_xy, r_xz, r_yz, n,
                                    grad_rzu0, rxu_point, alpha = 0.05) {
      if (is.null(data) || is.na(rxu_point)) {
        return(list(
          point_estimate = rxu_point,
          point_estimate_bias_corrected = rxu_point,
          ci_bias_corrected = c(L = NA, U = NA),
          bias_correction = 0
        ))
      }
      
      if (!requireNamespace("numDeriv", quietly = TRUE)) {
        se_simple <- 1/sqrt(n)
        z <- qnorm(1 - alpha/2)
        return(list(
          point_estimate = rxu_point,
          point_estimate_bias_corrected = rxu_point,
          ci_bias_corrected = c(L = rxu_point - z*se_simple,
                                U = rxu_point + z*se_simple),
          bias_correction = 0
        ))
      }
      
      theta_hat <- c(r_xy, r_xz, r_yz)
      
      g_fun_local <- function(theta) {
        theta[2] * rxu_point - (theta[1] * theta[2] - theta[3]) * sqrt((1 - rxu_point^2) / (1 - theta[1]^2))
      }
      
      grad_g0 <- grad_rzu0
      
      boot_corrs <- replicate(200, {
        idx <- sample(n, replace = TRUE)
        d <- data[idx, ]
        c(cor(d$x, d$y), cor(d$x, d$z), cor(d$y, d$z))
      })
      Sigma_hat <- cov(t(boot_corrs))
      
      var_g <- as.numeric(t(grad_g0) %*% Sigma_hat %*% grad_g0)
      se_g <- sqrt(var_g / n)
      
      hess_g <- numDeriv::hessian(g_fun_local, theta_hat)
      bias <- sum(0.5 * hess_g * Sigma_hat) / n
      rxu_bc <- rxu_point - bias
      
      z <- qnorm(1 - alpha/2)
      ci <- c(L = rxu_bc - z*se_g, U = rxu_bc + z*se_g)
      
      list(
        point_estimate = rxu_point,
        point_estimate_bias_corrected = rxu_bc,
        ci_bias_corrected = ci,
        bias_correction = bias
      )
    }
    
    rxuf <- delta_method_rho_xu(data, r_xy, r_xz, r_yz, n,
                                grad_rzu0, rxu_point, alpha = 0.05)
    
    rxu_point_corrected <- rxuf$point_estimate_bias_corrected
    ci_rxu <- rxuf$ci_bias_corrected
    
    # Bias correction for target value
    bias_z <- 0
    range_interval <- NULL
    bias_mc <- TRUE
    mc_B <- 200
    
    if (bias_mc && !is.null(data)) {
      boot_corrs <- replicate(mc_B, {
        idx <- sample(nrow(data), replace = TRUE)
        d <- data[idx, ]
        c(xz = cor(d$x, d$z, use = "complete.obs"),
          xy = cor(d$x, d$y, use = "complete.obs"),
          yz = cor(d$y, d$z, use = "complete.obs"))
      })
      valid <- complete.cases(t(boot_corrs))
      boot_corrs <- boot_corrs[, valid]
      
      if (ncol(boot_corrs) > 100) {
        sim_bias_z <- apply(boot_corrs, 2, function(corrs) {
          g_fun(corrs["xz"], corrs["xy"], corrs["yz"], rxu_point_corrected) -
            g_fun(r_xz, r_xy, r_yz, rxu_point_corrected)
        })
        bias_z <- mean(sim_bias_z, na.rm = TRUE)
        range_interval <- quantile(sim_bias_z, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
      }
    }
    
    bias_corrected_g <- 0#g_fun(rxu_point_corrected, r_xy, r_xz, r_yz) - bias_z
    
    
    #############################Experiment with coverage
    # Compute Sigma_g for coverage probability calculation
    # Need to compute covariance matrix of (g_min, g_max)
    n <- nrow(df)
    D_grid <- seq(rxu_range[1], rxu_range[2], length.out = 50)
    
    # Find r values that achieve min/max
    g_vals <- sapply(D_grid, function(r) g_fun(r, rho_xy, rho_xz, rho_yz))
    idx_min <- which.min(g_vals)
    idx_max <- which.max(g_vals)
    r_min <- D_grid[idx_min]
    r_max <- D_grid[idx_max]
    
    # Compute gradients at min/max achieving r values
    grad_min <- g_grad(r_min, rho_xy, rho_xz, rho_yz)
    grad_max <- g_grad(r_max, rho_xy, rho_xz, rho_yz)
    
    # Variance estimates
    var_xy <- (1 - rho_xy^2)^2 / n
    var_xz <- (1 - rho_xz^2)^2 / n
    var_yz <- (1 - rho_yz^2)^2 / n
    Sigma_rho <- diag(c(var_xy, var_xz, var_yz))
    
    
    # Check if zero is in CI
    contains_zero_b <- (CI_b[1] <= 0 & CI_b[2] >= 0)
    
    # Check if CI covers plug-in interval
    cover_b <- (plug_in[1] >= CI_b[1] & plug_in[2] <= CI_b[2])
    
    
    
    
    
    zero_in_CI_Bei <- (CI_b[1] <= 0 && 0 <= CI_b[2])
    
    ## THIS IS CHEAP WAY
    p_zero <- pvalue_mcub_zero_fast(deltahat, deltaSigma, Al, Au, g,
                                    eta = eta, tol = tol, tol_r = tol_r,
                                    B_fast = 600, Blarge_fast = 6000)
    
    
    
    res_list[[i]] <- data.frame(
      Z = paste0("Z", i),
      r_xz = round(cor(df$x, df$z), 2),
      r_zu = round(cor(df$z, df$u), 2),
      beta_IV = round(cov(df$z, df$y, use = "complete.obs") /
                        cov(df$z, df$x, use = "complete.obs"), 3),
      plug_in = sprintf("[%.2f,%.2f]", plug_in[1],plug_in[2]),
      CI_Bei = sprintf("[%.2f,%.2f]", CI_b[1], CI_b[2]),
      CI_simple = sprintf("[%.2f,%.2f]", max(-1,CI_s[1]), min(1,CI_s[2])),
      Zero_in_CI_Bei = ifelse(contains_zero_b, "✓", "×"),
      p_zero= round(p_zero,3)
      #
      #
      
    )
  }
  
  # Combine results
  summary_df <- do.call(rbind, res_list)
  rownames(summary_df) <- NULL
  
  return(summary_df)
}


#
#================================================================
# DATA GENERATING PROCESS
#================================================================



generate_data <- function(N = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Your data generation code
  b1 <- 15
  ds <- 2*b1/N
  vec <- seq(-b1, b1, ds)[1:N]
  x <- fs(vec, 0, 0.9) * 5 + rsnorm(N, mean = 1, sd = 0.25, xi = 1.5)
  
  u1 <- runif(N, min = -1, max = 1) + rsnorm(N, mean = 1, sd = 1, xi = 1.2)
  
  
  x1 <- (x - mean(x)) / sd(x)
  v <- resid(lm(u1 ~ x1))
  v <- v * mean(x) / 2 * sd(x) / sd(v) * 2
  
  e <- (x - mean(x)) * 2 + rnorm(N, sd = sd(x)) * 3
  u <- e + v
  
  beta <- 2
  y <- beta * x + u
  
  fity <- lm(y ~ x)
  V <- resid(fity)
  V <- (V - mean(V)) / sd(V) * sd(x)
  R <- V
  
  proxy_u <- resid(lm(y ~ x))
  A0_raw <- rnorm(N) * 0.75
  A0 <- resid(lm(A0_raw ~ proxy_u))
  
  data.frame(x = x, y = y, u = u, R = R, A0 = A0)
}


#================================================================
# 7. MONTE CARLO SIMULATION
#================================================================

monte_carlo_simulation <- function(n_sim = 10, N = 500, d_vec,
                                   rxu_range = c(0, 0.6),
                                   data_generating_process, sim_0=sim_0,
                                   alpha = 0.05) {
  
  all_results <- list()
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  
  for (sim in 1:n_sim) {
    sim_data <- data_generating_process(N = N, seed = sim_0 + sim)
    
    result <- analyze_single_sample(sim_data, d_vec, rxu_range = rxu_range,
                                    alpha = alpha, seed = sim_0+ sim)
    result$simulation_id <- sim
    
    all_results[[sim]] <- result
    setTxtProgressBar(pb, sim)
  }
  
  close(pb)
  
  combined <- do.call(rbind, all_results)
  rownames(combined) <- NULL
  
  return(combined)
}

#================================================================
# 8. SUMMARIZE RESULTS
#================================================================

summarize_sim_results <- function(sim_results) {
  library(dplyr)
  
  cat("Available columns:", paste(names(sim_results), collapse = ", "), "\n\n")
  
  summary <- sim_results %>%
    group_by(Z) %>%
    summarise(
      mean_r_xz = mean(r_xz, na.rm = TRUE),
      sd_r_xz = sd(r_xz, na.rm = TRUE),
      mean_r_zu = mean(r_zu, na.rm = TRUE),
      sd_r_zu = sd(r_zu, na.rm = TRUE),
      mean_beta_IV = mean(beta_IV, na.rm = TRUE),
      sd_beta_IV = sd(beta_IV, na.rm = TRUE),
      mean_p_zero = mean(p_zero, na.rm = TRUE),
      n_total = n()
    ) %>%
    ungroup()
  
  return(summary)
}

#================================================================
# 9. MAIN EXECUTION
#
# AUTO-RUN GUARD — same pattern as Fin_Empirical5_clean.R
#
# Set DIAL_SKIP_AUTORUN <- TRUE BEFORE sourcing this file when
# using DIAL_NEURIPS_CR_wrapper.R.  All three monte_carlo_simulation()
# calls below run on every source() unless this flag is set.
# Without the flag the three calls take ~10-15 min each time the
# file is sourced, causing the repeated output seen overnight.
#
# Recommended session startup:
#   DIAL_SKIP_AUTORUN <- TRUE
#   source("Fin_sim3_clean.R")
#   source("Fin_Empirical5_clean.R")
#   source("DIAL_NEURIPS_CR_wrapper.R")
#   res <- run_all_scenarios(n_sim = 50)
#================================================================

.dial_autorun_sim <- !exists("DIAL_SKIP_AUTORUN") || !isTRUE(DIAL_SKIP_AUTORUN)

if (.dial_autorun_sim) {
  
  d_vec <- c(1.25, 1.10, 0.95, 0.83, 0.7, 0.62, 0.50, 0.45, 0.35, 0.2, 0.05, -0.05)
  cat("\nRunning Monte Carlo simulation...\n")
  
  mc_results <- monte_carlo_simulation(
    n_sim = 1,
    N = 800,
    d_vec = d_vec,
    rxu_range = c(0.0, 0.8),
    data_generating_process = generate_data,
    alpha = 0.05,
    sim_0 = 9058
  )
  
  cat("\n=== Monte Carlo Results Structure ===\n")
  print(head(mc_results))
  cat("\nDimensions:", nrow(mc_results), "rows,", ncol(mc_results), "columns\n")
  cat("Unique instruments:", length(unique(mc_results$Z)), "\n")
  cat("Number of simulations:", length(unique(mc_results$simulation_id)), "\n")
  
  mc_summary <- summarize_sim_results(mc_results)
  cat("\n=== Monte Carlo Summary by Instrument ===\n")
  print(mc_summary, n = Inf)
  cat("\nSimulation complete! Results saved to mc_summary_clean.csv\n")
  
  # Single sample case A
  d_vec_A <- c(0.62, 0.50, 0.45, 0.35, 0.2, 0.05)
  cat("\nRunning Monte Carlo simulation (case A)...\n")
  mc_results_A <- monte_carlo_simulation(
    n_sim = 1, N = 2000, d_vec = d_vec_A,
    rxu_range = c(0, 0.8),
    data_generating_process = generate_data,
    alpha = 0.05, sim_0 = 9058
  )
  print(t(mc_results_A))
  
  # Single sample case B
  d_vec_B <- c(1.25, 1.10, 0.95, 0.83, 0.7, 0.62)
  cat("\nRunning Monte Carlo simulation (case B)...\n")
  mc_results_B <- monte_carlo_simulation(
    n_sim = 1, N = 2000, d_vec = d_vec_B,
    rxu_range = c(0, 0.8),
    data_generating_process = generate_data,
    alpha = 0.05, sim_0 = 9058
  )
  print(t(mc_results_B))
  
} else {
  cat("  [Fin_sim3.R] DIAL_SKIP_AUTORUN = TRUE — skipping all MC auto-runs.\n")
  cat("  Use run_all_scenarios() or scenario2_mc_power() to run simulations.\n")
}
