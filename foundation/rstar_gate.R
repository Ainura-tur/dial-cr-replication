# ==============================================================================
# rstar_gate.R  -  the validity gate of DIAL: a confidence set for the
#                  implied endogeneity r*, obtained by inverting the CR mapping
# ==============================================================================
#
# IDEA
#   For residualised (D, Y, Z), the CR mapping gives the instrument-error
#   correlation implied by a treatment-error correlation r:
#     g(rho, r) = rho_DZ r - (rho_DY rho_DZ - rho_YZ) sqrt((1 - r^2) / (1 - rho_DY^2))
#   with rho = (rho_DY, rho_DZ, rho_YZ). The implied endogeneity r* is the
#   unique zero of g(rho, .), and the instrument is compatible with the
#   domain D exactly when r* lies in the closure of D (prop:implied(a)).
#
#   For each candidate r, the statistic t(r) = g(rho_hat, r) / se(r) tests
#   g(rho, r) = 0, that is r* = r. Collecting the r that the test does not
#   reject gives a confidence set for r*:
#     CS = { r : |t(r)| <= z_{1 - alpha0} },
#   which covers r* with probability 1 - 2 alpha0 in large samples. The gate
#   rejects validity when CS misses the closed domain. At an edge of D this
#   test has size alpha0; inside D its size tends to zero. Like the
#   Anderson-Rubin test, it never divides by rho_DZ, so a weak first stage
#   does not make t(r) unstable. It does not make CS uninformative either:
#   when the reduced form is weak too, g is flat and CS is wide, but when
#   only the first stage is weak, CS shrinks towards r = -1 and r = 1 (two
#   short pieces near +/-1, as for the two weak MovieLens instruments).
#   Algorithm 1 therefore applies the relevance gate (F >= 10) first.
#
#   se(r) is the bootstrap standard deviation of g(rho_b, r), where rho_b
#   are the correlations of a resample of the raw rows (the residualisation
#   on the controls is redone in each resample).
#
#   The gate uses the correlations of the residualised instrument z0 itself,
#   not of the prepared instrument (its projection on the residualised D and
#   Y, used by the MCUB code): on the prepared scale the correlation matrix
#   is singular and g(rho_hat, r) = sin(arcsin r - arcsin r*), which
#   carries no information beyond r*.
#
# FUNCTIONS
#   dial_raw_rho(eng, w)        raw residualised correlations (DY, DZ, YZ)
#                               for row counts w, plus m = sd(z_hat)/sd(z0)
#   dial_g(rho, r)              the CR mapping (vectorised in r)
#   dial_rstar(rho)             the zero of g(rho, .)
#   dial_rstar_gate(eng, D, alpha0, B, seed, step)
#       returns a list with r_star, the confidence set (cs_lo, cs_hi, its
#       segments), V (1 = validity not rejected), p_D (the largest one-sided
#       p-value over the closed domain; V = 1 iff p_D >= alpha0), the edge
#       nearest to r*, kappa = |g(rho_hat, edge)|, s_edge = sqrt(n) se(edge),
#       and the Gaussian power prediction for the edge.
#   dial_kappa_star(n, s, alpha0, pi_star)
#       the classifiability frontier of the gate: the smallest |g| at the
#       nearest edge that the gate detects with probability 1 - pi_star.
#
# Needs dial_raw_spec() and dial_prep_engine() from bootstrap_rule.R.
# ==============================================================================

dial_g <- function(rho, r)
  rho[2] * r - (rho[1] * rho[2] - rho[3]) * sqrt((1 - r^2) / (1 - rho[1]^2))

dial_rstar <- function(rho) {
  if (any(!is.finite(rho))) return(NA_real_)
  k <- (rho[1] * rho[2] - rho[3]) / (rho[2] * sqrt(1 - rho[1]^2))
  unname(k / sqrt(1 + k^2))
}

dial_raw_rho <- function(eng, w = rep(1, eng$n)) {
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
  c(DY = sxy / sqrt(sxx * syy), DZ = sxz / sqrt(sxx * szz),
    YZ = syz / sqrt(syy * szz), m = sqrt(q / szz))
}

dial_kappa_star <- function(n, s, alpha0 = 0.05, pi_star = 0.05)
  (stats::qnorm(1 - alpha0) + stats::qnorm(1 - pi_star)) * s / sqrt(n)

dial_rstar_gate <- function(eng, D, alpha0 = 0.05, B = 1000L, seed = 1L,
                            step = 0.001) {
  D   <- sort(D)
  rho <- unname(dial_raw_rho(eng)[1:3])
  n   <- eng$n
  set.seed(seed)
  RB <- t(vapply(seq_len(B), function(b) {
    w <- tabulate(sample.int(n, n, replace = TRUE), nbins = n)
    dial_raw_rho(eng, w)[1:3]
  }, numeric(3)))
  RB <- RB[stats::complete.cases(RB) & abs(RB[, 1]) < 1, , drop = FALSE]
  # grid of candidate r, with the edges of D added exactly
  grid <- sort(unique(c(seq(-1 + step, 1 - step, by = step), D)))
  grid <- grid[grid > -1 & grid < 1]
  g_hat <- dial_g(rho, grid)
  a_b <- RB[, 1] * RB[, 2] - RB[, 3]
  gb <- vapply(grid, function(r)
    RB[, 2] * r - a_b * sqrt((1 - r^2) / (1 - RB[, 1]^2)), numeric(nrow(RB)))
  se <- apply(gb, 2, stats::sd)
  tt <- g_hat / se
  zc <- stats::qnorm(1 - alpha0)
  acc <- abs(tt) <= zc
  inD <- grid >= D[1] & grid <= D[2]
  V   <- as.integer(any(acc & inD))
  p_D <- max(1 - stats::pnorm(abs(tt[inD])))
  # segments of the confidence set on the grid
  seg <- NULL
  if (any(acc)) {
    idx <- which(acc); br <- c(0, which(diff(idx) > 1), length(idx))
    seg <- t(vapply(seq_len(length(br) - 1), function(k)
      grid[c(idx[br[k] + 1], idx[br[k + 1]])], numeric(2)))
    colnames(seg) <- c("lo", "hi")
  }
  rs <- dial_rstar(rho)
  edge <- D[which.min(abs(D - rs))]
  if (rs >= D[1] && rs <= D[2]) edge <- NA_real_
  ie <- if (is.na(edge)) NA_integer_ else which(grid == edge)[1]
  kappa  <- if (is.na(edge)) 0 else abs(g_hat[ie])
  s_edge <- if (is.na(edge)) NA_real_ else sqrt(n) * se[ie]
  list(r_star = rs, rho = rho, n = n, D = D, alpha0 = alpha0, B = nrow(RB),
       cs_lo = if (is.null(seg)) NA_real_ else min(seg[, "lo"]),
       cs_hi = if (is.null(seg)) NA_real_ else max(seg[, "hi"]),
       segments = seg, n_segments = if (is.null(seg)) 0L else nrow(seg),
       V = V, p_D = p_D, edge = edge, kappa = kappa, s_edge = s_edge,
       power_gauss = if (is.na(edge)) NA_real_ else
         stats::pnorm(sqrt(n) * kappa / s_edge - zc),
       dist = if (is.na(edge)) 0 else abs(rs - edge))
}

# Convenience wrapper: data frame and column names in, gate out.
dial_rstar_gate_data <- function(data, Y, X, Z, controls = character(0), D,
                                 alpha0 = 0.05, B = 1000L, seed = 1L,
                                 step = 0.001) {
  spec <- dial_raw_spec(data, Y, X, Z, controls)
  dial_rstar_gate(dial_prep_engine(spec), D, alpha0 = alpha0, B = B,
                  seed = seed, step = step)
}
