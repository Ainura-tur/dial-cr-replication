# ==============================================================================
# 14_hsic_rejection_mc.R  -  repeated-sample rejection rates of HSIC and the
#                            r* gate in the design of tab:synthetic
# ==============================================================================
#
# PURPOSE
#   tab:synthetic compares the gate with HSIC on one sample per instrument
#   (script 04). Its generator fixes the sample values of cor(Z, X) and
#   cor(Z, U) at their targets, so it cannot measure size or power (app:mc).
#   This script draws R_REPS samples per instrument from the population
#   version of the same Jones-Pewsey design and reports the share of samples
#   in which each test rejects.
#
#   SAMPLING DESIGN. As in parts A and B of 13, the instrument uses fixed
#   population coefficients and fresh noise in every sample, so the sample
#   correlations vary across samples. In units of sd(X) (app:dgp_jp):
#     X      = (5 H(nu; 0, 0.9) + S(.; 0, 1.2)) / SD_JP, standardised by its
#              population sd
#     U      = lambda X + eps + v / 2, eps ~ N(0, 1), and v the
#              Uniform(-1, 1) + S(.; 0, 1.2) draw standardised by its
#              population sd
#     lambda = rho sqrt(1.25 / (1 - rho^2)), so cor(X, U) = rho = +/-0.30
#     U_perp = (eps + v / 2) / sqrt(1.25), the part of U orthogonal to X
#     Z      = a X + c2 U_perp + c3 W, W ~ N(0, 1), as in .construct_z()
#              of jp_dgp.R but with population coefficients
#     Y      = 2 X + U
#   Then cor(Z, X) = a = 0.25 and cor(Z, U) = rho_ZU hold in the population.
#   z_het follows draw_het_sample() of 04, with z and u0 standardised by
#   their population moments and the loading r set so that the population
#   cor(X, U) = 0.30.
#
#   As in tab:synthetic, rho_ZU < 0 pairs with lambda < 0 and D = (-0.8, 0);
#   the other rows and z_het use lambda > 0 and D = (0, 0.8). With
#   OPPOSITE_SIGN = TRUE, the script adds rho_ZU < 0 under lambda > 0 and
#   D = (0, 0.8), the direction in which the gate has little power.
#
#   TESTS. The gate is dial_rstar_gate() of foundation/rstar_gate.R at
#   alpha0 = 0.05 with B = 500 and grid step 0.002, as in part A of 13.
#   HSIC is the call of 04, dHSIC::dhsic.test(list(z, e_hat),
#   method = "gamma"), at the 5% level, with e_hat = y - beta_IV x. beta_IV
#   is the IV slope with an intercept, the slope that ivreg() returns in 04.
#
# OUTPUT
#   results/hsic_rejection_mc.csv
#   Console: one line per instrument as it finishes, and LaTeX rows for a
#   table in app:mc.
#
# RUNTIME  about 12 minutes on 2 cores at R_REPS = 500 (12 instruments), and
#   less on more cores. The samples run in parallel on Linux and macOS only,
#   as in 13; Windows runs them on one core, in about 17 minutes. The first
#   console line reports the cores in use. DIAL_HSIC_REPS=20 gives a quick
#   test run.
#
# RUN (from the scripts folder)
#   Rscript 14_hsic_rejection_mc.R
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({ library(dplyr); library(dHSIC) })

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

ALPHA0        <- 0.05     # gate level
ALPHA_HSIC    <- 0.05     # HSIC level, as in 04
N_OBS         <- 800L     # sample size of tab:synthetic
R_REPS        <- as.integer(Sys.getenv("DIAL_HSIC_REPS", "500"))
B_GATE        <- 500L     # bootstrap resamples of the gate, as in part A of 13
STEP_MC       <- 0.002    # grid step of the gate, as in part A of 13
COR_ZX        <- 0.25
TARGET_COR_XU <- 0.30
OPPOSITE_SIGN <- FALSE
MASTER_SEED   <- 5151L
N_CORES       <- max(1L, parallel::detectCores() - 0L)
# mclapply() forks, so the samples run in parallel on Linux and macOS only
N_USED        <- if (.Platform$OS.type == "unix") N_CORES else 1L

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
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# The gate needs only the preparation step of bootstrap_rule.R. That file
# also wraps the MCUB code of Fin_Empirical5_clean.R, which this script does
# not use, and warns that it is missing; the warning is suppressed.
suppressWarnings(invisible(capture.output(
  source(file.path(FOUNDATION_DIR, "bootstrap_rule.R")))))
source(file.path(FOUNDATION_DIR, "jp_dgp.R"))
source(file.path(FOUNDATION_DIR, "rstar_gate.R"))

gate_xyz <- function(y, x, z, D, B, seed) {
  eng <- dial_prep_engine(dial_raw_spec(data.frame(Y = y, X = x, Z = z), "Y", "X", "Z"))
  dial_rstar_gate(eng, D, alpha0 = ALPHA0, B = B, seed = seed, step = STEP_MC)
}

par_lapply <- function(X, FUN) {
  if (N_USED > 1L)
    parallel::mclapply(X, FUN, mc.cores = N_USED, mc.preschedule = FALSE) else lapply(X, FUN)
}

# ==============================================================================
# POPULATION DESIGN
# ==============================================================================

# Population sds by numerical integration; every term has mean zero
.e_norm <- function(f) stats::integrate(function(t) f(t) * stats::dnorm(t), -Inf, Inf)$value
VAR_H <- .e_norm(function(t) jp_h_transform(t, zeta = 0, kappa = 0.9)^2)
VAR_S <- .e_norm(function(t) sinh(asinh(t) / 1.2)^2)   # law of jp_skew_normal(kappa = 1.2)
SD_JP <- sqrt(25 * VAR_H + VAR_S)   # sd of 5 H(nu; 0, 0.9) + S(.; 0, 1.2)
SD_U1 <- sqrt(1 / 3 + VAR_S)        # sd of Uniform(-1, 1) + S(.; 0, 1.2)

draw_jp_std <- function(n)
  (5 * jp_h_transform(stats::rnorm(n), zeta = 0, kappa = 0.9) +
     jp_skew_normal(n, zeta = 0, kappa = 1.2)) / SD_JP
draw_u1_std <- function(n)
  (stats::runif(n, -1, 1) + jp_skew_normal(n, zeta = 0, kappa = 1.2)) / SD_U1

# z_het: with Var(z) = Var(u0) = 1, u0 independent of z and
# s(z) = sqrt((1 + z^2) / 2): Cov(x, u) = r E[s(z)], Var(u) = 1 and
# Var(x) = pi^2 + 1, so r = 0.30 sqrt(pi^2 + 1) / E[s(z)].
PI_HET <- COR_ZX / sqrt(1 - COR_ZX^2)
set.seed(MASTER_SEED)
E_S    <- mean(sqrt((1 + draw_jp_std(2e6)^2) / 2))
R_HET  <- TARGET_COR_XU * sqrt(PI_HET^2 + 1) / E_S

sim_grid <- function(n, rho, cor_zu, seed) {
  set.seed(seed)
  lam   <- rho * sqrt(1.25 / (1 - rho^2))
  x     <- draw_jp_std(n)
  eps   <- stats::rnorm(n)
  v     <- draw_u1_std(n)
  u     <- lam * x + eps + 0.5 * v
  uperp <- (eps + 0.5 * v) / sqrt(1.25)
  c2    <- (cor_zu - COR_ZX * rho) / sqrt(1 - rho^2)
  c3    <- sqrt(1 - COR_ZX^2 - c2^2)
  z     <- COR_ZX * x + c2 * uperp + c3 * stats::rnorm(n)
  data.frame(x = x, y = 2 * x + u, z = z, u = u)
}

sim_het <- function(n, seed) {
  set.seed(seed)
  z  <- draw_jp_std(n)
  u0 <- draw_u1_std(n)
  u  <- sqrt((1 + z^2) / 2) * u0
  x  <- PI_HET * z + R_HET * u0 + sqrt(1 - R_HET^2) * stats::rnorm(n)
  data.frame(x = x, y = 2 * x + u, z = z, u = u)
}

# Population r* = cor(X, Y - beta_IV X) of a grid instrument; it depends
# only on rho = cor(X, U), a = cor(Z, X) and delta = cor(Z, U)
rstar_pop <- function(rho, a, delta) {
  k <- delta / a
  (rho - k) / sqrt(1 - 2 * rho * k + k^2)
}

# ==============================================================================
# ONE SAMPLE
# ==============================================================================

one_sample <- function(rw, seed) {
  df <- if (rw$het) sim_het(N_OBS, seed) else sim_grid(N_OBS, rw$rho_xu, rw$cor_zu, seed)
  g  <- gate_xyz(df$y, df$x, df$z, c(rw$d_lo, rw$d_hi), B_GATE, seed + 1000000L)
  beta_iv <- stats::cov(df$z, df$y) / stats::cov(df$z, df$x)
  e_hat   <- df$y - beta_iv * df$x
  hsic_p  <- tryCatch(
    dHSIC::dhsic.test(list(df$z, e_hat), method = "gamma", B = 100)$p.value,
    error = function(e) NA_real_)
  c(V = g$V, hsic_p = hsic_p, rstar = g$r_star, cor_zu = stats::cor(df$z, df$u))
}

# ==============================================================================
# RUN
# ==============================================================================

GRID <- round(JP_INSTRUMENT_GRID, 2)   # -0.20, -0.16, ..., 0.20; exact zero
rows <- lapply(GRID, function(cz)
  list(label = sprintf("z_%+0.2f", cz), cor_zu = cz,
       rho_xu = if (cz < 0) -TARGET_COR_XU else TARGET_COR_XU, het = FALSE))
if (OPPOSITE_SIGN)
  rows <- c(rows, lapply(GRID[GRID < 0], function(cz)
    list(label = sprintf("z_%+0.2f, lambda > 0", cz), cor_zu = cz,
         rho_xu = TARGET_COR_XU, het = FALSE)))
rows <- c(rows, list(list(label = "z_het", cor_zu = 0, rho_xu = TARGET_COR_XU, het = TRUE)))
for (k in seq_along(rows)) {
  rows[[k]]$d_lo <- if (rows[[k]]$rho_xu < 0) -0.8 else 0
  rows[[k]]$d_hi <- if (rows[[k]]$rho_xu < 0) 0 else 0.8
  rows[[k]]$rstar_pop <- if (rows[[k]]$het) TARGET_COR_XU else
    rstar_pop(rows[[k]]$rho_xu, COR_ZX, rows[[k]]$cor_zu)
}

cat("================================================================================\n")
cat(sprintf("HSIC and r* gate over repeated samples: n = %d, %d samples per instrument, %s\n",
            N_OBS, R_REPS,
            if (N_USED > 1L) sprintf("%d cores", N_USED) else
              if (N_CORES > 1L) "1 core (Windows runs the samples on one core)" else "1 core"))
cat("================================================================================\n")
cat(sprintf("  population constants: sd(x) = %.4f, sd(u1) = %.4f, r_het = %.4f\n\n",
            SD_JP, SD_U1, R_HET))

out <- list()
for (k in seq_along(rows)) {
  rw    <- rows[[k]]
  seeds <- MASTER_SEED + 10000L * k + seq_len(R_REPS)
  t0    <- Sys.time()
  M     <- do.call(rbind, par_lapply(seeds, function(sd) one_sample(rw, sd)))
  ok    <- !is.na(M[, "hsic_p"])
  o <- data.frame(instrument  = rw$label, cor_zu = rw$cor_zu, rho_xu = rw$rho_xu,
                  D_lo = rw$d_lo, D_hi = rw$d_hi, rstar_pop = rw$rstar_pop,
                  R = nrow(M),
                  gate_reject = mean(M[, "V"] == 0),
                  hsic_reject = mean(M[ok, "hsic_p"] < ALPHA_HSIC),
                  hsic_ok     = sum(ok),
                  mean_rstar  = mean(M[, "rstar"]),
                  mean_cor_zu = mean(M[, "cor_zu"]))
  out[[k]] <- o
  cat(sprintf("  %-22s r*_pop %+.3f  gate reject %.3f  HSIC reject %.3f  mean r* %+.3f  mean cor(z,u) %+.3f  (%.1f min)\n",
              o$instrument, o$rstar_pop, o$gate_reject, o$hsic_reject,
              o$mean_rstar, o$mean_cor_zu,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, out)
utils::write.csv(res, file.path(RESULTS_DIR, "hsic_rejection_mc.csv"), row.names = FALSE)

if (any(res$hsic_ok < res$R))
  cat("\nNOTE: dhsic.test() failed in some samples; HSIC rates use the samples where it ran.\n")
cat(sprintf("\nMonte Carlo standard error: about %.3f at a rate of 0.05, at most %.3f.\n",
            sqrt(0.05 * 0.95 / R_REPS), sqrt(0.25 / R_REPS)))

# LaTeX rows: Instrument & rho_ZU & D & population r* & gate & HSIC
cat("\nLaTeX rows (Instrument & $\\rho_{ZU}$ & $\\mathcal{D}$ & population $r^{*}$ & gate & HSIC):\n\n")
for (i in seq_len(nrow(res))) {
  r  <- res[i, ]
  nm <- if (r$instrument == "z_het") "$Z_{\\text{het}}$" else
    sprintf("$Z_{%s}$", sub("^z_", "", r$instrument))
  if (r$instrument == "z_het") cat("\\midrule\n")
  rho_txt <- if (r$cor_zu == 0) "$\\phantom{+}0.00$" else sprintf("$%+.2f$", r$cor_zu)
  cat(sprintf("%s & %s & $(%g, %g)$ & $%+.3f$ & $%.3f$ & $%.3f$ \\\\\n",
              nm, rho_txt, r$D_lo, r$D_hi, r$rstar_pop, r$gate_reject, r$hsic_reject))
}

cat(sprintf("\nWrote %s\n", file.path(RESULTS_DIR, "hsic_rejection_mc.csv")))
cat("DONE.\n")
