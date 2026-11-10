# ==============================================================================
# 13_rstar_gate_mc.R  -  size and power of the r* gate
#                        (tab:gate_mc, fig:power_curves, the gate columns of
#                         tab:synthetic and tab:sim)
# ==============================================================================
#
# PURPOSE
#   Monte Carlo checks of the validity gate of Algorithm 1 (the confidence
#   set for r* obtained by inverting the CR mapping; foundation/rstar_gate.R).
#
#   A  Size and power at n = 800: rho_DU = 0.30, cor(Z, D) = 0.25 and
#      delta = cor(Z, U) in {0, 0.04, ..., 0.20} under D = [0, 0.8]; plus
#      the edge cell, a valid instrument with no endogeneity (rho_DU = 0,
#      so r* = 0) under D = [0, 0.8]. Rejection rate of the gate and mean
#      width of the confidence set over R_A samples.
#   B  Power curves with rho_DU = 0.5 and cor(Z, D) = 0.25 (the design of
#      07) for n in N_B and delta on DELTA_B.
#
#   SAMPLING DESIGN. Parts A and B draw (D, U) from a bivariate normal with
#   correlation rho_DU and set Z = a D + c2 V + c3 W, with V the part of U
#   orthogonal to D in the population and W ~ N(0, 1) independent. Then
#   cor(Z, D) = a and cor(Z, U) = delta hold in the population, and the
#   sample correlations vary across samples. The generators of 04, 05, 07
#   and 10 instead residualise W on (D, U) within each sample, which makes
#   the sample values of cor(Z, D) and cor(Z, U) equal their targets in every
#   sample and removes most of the sampling variation in r*; a size or power
#   study needs the population version. With Y = beta D + U, the population
#   r* equals zero exactly when delta = rho_DU a, so the boundary is
#   delta*_inf = 0.075 in part A and 0.125 in part B.
#   C  The gate on the twelve samples of tab:synthetic (the samples of 04).
#   D  The gate on the seven samples of tab:sim (the samples of 06).
#
#   Parts C and D rebuild the samples by evaluating the data-generating
#   part of 04 and 06 (their top-level code up to the first CR call), so
#   the samples are exactly those of the tables.
#
# OUTPUT
#   results/rstar_gate_mc_cells.csv     part A
#   results/rstar_gate_power.csv        part B
#   results/rstar_gate_synthetic.csv    part C
#   results/rstar_gate_sim.csv          part D
#
# RUNTIME  about 40 minutes on 2 cores (part B at n = 50,000 dominates).
#   The cells run in parallel on Linux and macOS only; on Windows the script
#   uses one core and takes about 70 minutes. DIAL_GATE_PARTS=A,C,D skips
#   part B.
#
# RUN (from the scripts folder)
#   Rscript 13_rstar_gate_mc.R
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({ library(dplyr) })

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

ALPHA0  <- 0.05
N_A     <- 800L
R_A     <- 500L
B_A     <- 500L
CELLS_A <- c(0.00, 0.04, 0.08, 0.12, 0.16, 0.20)
N_B     <- c(1107L, 5000L, 50000L)
R_B     <- 200L
B_B     <- 300L
DELTA_B <- c(0, 0.04, 0.08, 0.10, 0.11, 0.12, 0.125, 0.13, 0.14, 0.15, 0.17, 0.20, 0.25)
B_CD    <- 1000L
STEP_MC <- 0.002
MASTER_SEED <- 4242L
N_CORES <- max(1L, parallel::detectCores() - 0L)
PARTS   <- toupper(trimws(strsplit(Sys.getenv("DIAL_GATE_PARTS", "A,B,C,D"), ",")[[1]]))

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

assign("DIAL_SKIP_AUTORUN", TRUE, envir = .GlobalEnv)
invisible(capture.output(source(file.path(FOUNDATION_DIR, "Fin_Empirical5_clean.R"))))
invisible(capture.output(source(file.path(FOUNDATION_DIR, "bootstrap_rule.R"))))
source(file.path(FOUNDATION_DIR, "jp_dgp.R"))
source(file.path(FOUNDATION_DIR, "rstar_gate.R"))

gate_xyz <- function(y, x, z, D, B, seed, controls_df = NULL) {
  d <- data.frame(Y = y, X = x, Z = z)
  ctr <- character(0)
  if (!is.null(controls_df)) { d <- cbind(d, controls_df); ctr <- names(controls_df) }
  eng <- dial_prep_engine(dial_raw_spec(d, "Y", "X", "Z", ctr))
  dial_rstar_gate(eng, D, alpha0 = ALPHA0, B = B, seed = seed, step = STEP_MC)
}

# Evaluate the named top-level function definitions of a script
.defs_from <- function(file, names, envir) {
  ex <- parse(file)
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        is.name(e[[2]]) && as.character(e[[2]]) %in% names)
      eval(e, envir)
  }
}
# Evaluate the top-level code of a script up to the first expression whose
# text starts with `marker`
.eval_until <- function(file, marker, envir) {
  ex <- parse(file, keep.source = TRUE)
  src <- vapply(attr(ex, "srcref"), function(s) paste(as.character(s), collapse = "\n"), "")
  stop_at <- which(startsWith(trimws(src), marker))[1]
  if (is.na(stop_at)) stop("marker not found in ", file, ": ", marker)
  for (i in seq_len(stop_at - 1L)) eval(ex[[i]], envir)
  invisible(envir)
}

par_lapply <- function(X, FUN) {
  if (.Platform$OS.type == "unix" && N_CORES > 1L)
    parallel::mclapply(X, FUN, mc.cores = N_CORES, mc.preschedule = FALSE) else lapply(X, FUN)
}

cat("================================================================================\n")
cat(sprintf("r* gate Monte Carlo: level %.2f, parts %s, %d cores\n",
            ALPHA0, paste(PARTS, collapse = ","), N_CORES))
cat("================================================================================\n\n")

# Population design: (D, U) bivariate normal, Z = a D + c2 V + c3 W
sim_pop <- function(n, rho, a, delta, beta, seed) {
  set.seed(seed)
  Dt <- stats::rnorm(n)
  Vt <- stats::rnorm(n)
  U  <- rho * Dt + sqrt(1 - rho^2) * Vt
  c2 <- (delta - a * rho) / sqrt(1 - rho^2)
  c3 <- sqrt(1 - a^2 - c2^2)
  Z  <- a * Dt + c2 * Vt + c3 * stats::rnorm(n)
  data.frame(x = Dt, y = beta * Dt + U, z = Z)
}
# Population r* of this design
rstar_pop <- function(rho, a, delta, beta) {
  sY <- sqrt(beta^2 + 2 * beta * rho + 1)
  dial_rstar(c((beta + rho) / sY, a, (beta * a + delta) / sY))
}

# ==============================================================================
# A: size and power at n = 800, and the edge cell
# ==============================================================================

if ("A" %in% PARTS) {
  cat("Part A: cells at n = 800 (population design)\n")
  cells <- c(lapply(CELLS_A, function(cz) list(cz = cz, edge = FALSE)),
             list(list(cz = 0, edge = TRUE)))
  outA <- list()
  for (k in seq_along(cells)) {
    cl <- cells[[k]]
    rho <- if (cl$edge) 0 else 0.30
    r0  <- rstar_pop(rho, 0.25, cl$cz, 2)
    seeds <- MASTER_SEED + 1000L * k + seq_len(R_A)
    t0 <- Sys.time()
    rr <- par_lapply(seeds, function(sd) {
      df <- sim_pop(N_A, rho, 0.25, cl$cz, 2, sd)
      g <- gate_xyz(df$y, df$x, df$z, c(0, 0.8), B_A, sd)
      cov <- !is.null(g$segments) &&
        any(g$segments[, "lo"] <= r0 + STEP_MC & g$segments[, "hi"] >= r0 - STEP_MC)
      c(V = g$V, width = g$cs_hi - g$cs_lo, rstar = g$r_star, cover = cov)
    })
    M <- do.call(rbind, rr)
    o <- data.frame(cell = if (cl$edge) "edge" else "grid", cor_zu = cl$cz,
                    rho_DU = rho, rstar_pop = r0, R = nrow(M),
                    reject = mean(M[, "V"] == 0), coverage = mean(M[, "cover"]),
                    mean_width = mean(M[, "width"]), mean_rstar = mean(M[, "rstar"]),
                    sd_rstar = stats::sd(M[, "rstar"]))
    outA[[length(outA) + 1L]] <- o
    cat(sprintf("  %-4s rho_ZU = %+.2f  r*_pop %+.3f  reject %.3f  coverage %.3f  mean width %.3f  r*: mean %+.3f sd %.3f  (%.1f min)\n",
                o$cell, o$cor_zu, o$rstar_pop, o$reject, o$coverage, o$mean_width,
                o$mean_rstar, o$sd_rstar,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  utils::write.csv(do.call(rbind, outA), file.path(RESULTS_DIR, "rstar_gate_mc_cells.csv"),
                   row.names = FALSE)
}

# ==============================================================================
# C: the samples of tab:synthetic (04)
# ==============================================================================

if ("C" %in% PARTS) {
  cat("\nPart C: tab:synthetic samples (04)\n")
  e4 <- new.env()
  .defs_from(file.path(.script_dir, "04_synthetic_benchmark.R"), "draw_het_sample", e4)
  N4 <- 800L; MS4 <- 1234L; COR_ZX4 <- 0.25
  e4$COR_ZX <- COR_ZX4; e4$TARGET_COR_XU <- 0.30
  environment(e4$draw_het_sample) <- e4
  c_pos <- calibrate_lambda(n = N4, target_cor_xu = 0.30, cor_zx_default = COR_ZX4,
                            seed = MS4, sign = "positive", verbose = FALSE)
  c_neg <- calibrate_lambda(n = N4, target_cor_xu = 0.30, cor_zx_default = COR_ZX4,
                            seed = MS4, sign = "negative", verbose = FALSE)
  s_pos <- draw_jp_sample(n = N4, lambda = c_pos$lambda, seed = MS4 + 1L,
                          cor_zu_grid = JP_INSTRUMENT_GRID, cor_zx_default = COR_ZX4)
  s_neg <- draw_jp_sample(n = N4, lambda = c_neg$lambda, seed = MS4 + 2L,
                          cor_zu_grid = JP_INSTRUMENT_GRID, cor_zx_default = COR_ZX4)
  outC <- list()
  for (j in seq_along(JP_INSTRUMENT_GRID)) {
    cz <- JP_INSTRUMENT_GRID[j]
    sm <- if (cz < 0) s_neg else s_pos
    D  <- if (cz < 0) c(-0.8, 0) else c(0, 0.8)
    g  <- gate_xyz(sm$y, sm$x, sm$z[, j], D, B_CD, 7000L + j)
    outC[[j]] <- data.frame(instrument = sprintf("z_%+0.2f", cz), cor_zu = cz,
                            D_lo = D[1], D_hi = D[2], r_star = g$r_star,
                            cs_lo = g$cs_lo, cs_hi = g$cs_hi, V = g$V, p_D = g$p_D)
  }
  sh <- e4$draw_het_sample(N4, seed = MS4 + 3L)
  g  <- gate_xyz(sh$y, sh$x, sh$z, c(0, 0.8), B_CD, 7100L)
  outC[[length(outC) + 1L]] <- data.frame(instrument = "z_het", cor_zu = 0, D_lo = 0, D_hi = 0.8,
                                          r_star = g$r_star, cs_lo = g$cs_lo, cs_hi = g$cs_hi,
                                          V = g$V, p_D = g$p_D)
  resC <- do.call(rbind, outC)
  print(resC, digits = 3, row.names = FALSE)
  utils::write.csv(resC, file.path(RESULTS_DIR, "rstar_gate_synthetic.csv"), row.names = FALSE)
}

# ==============================================================================
# D: the samples of tab:sim (06)
# ==============================================================================

if ("D" %in% PARTS) {
  cat("\nPart D: tab:sim samples (06)\n")
  e6 <- new.env()
  invisible(capture.output(
    .eval_until(file.path(.script_dir, "06_structural_simulation.R"), "ROWS <- list(", e6)))
  rows6 <- list(
    list(lab = "A: homogeneous",          df = e6$df_hom,  ctr = e6$H_AND_C),
    list(lab = "B: heterogeneous effects", df = e6$df_het,  ctr = e6$H_AND_C),
    list(lab = "B: C unobserved",          df = e6$df_hom,  ctr = e6$H_ONLY),
    list(lab = "C: exclusion violation",   df = e6$df_inv,  ctr = e6$H_AND_C),
    list(lab = "C: weak instrument",       df = e6$df_weak, ctr = e6$H_AND_C),
    list(lab = "het. rev.",                df = e6$df_rev,  ctr = e6$H_AND_C),
    list(lab = "het. rev., C hidden",      df = e6$df_rev,  ctr = e6$H_ONLY))
  outD <- list()
  for (k in seq_along(rows6)) {
    rw <- rows6[[k]]
    eng <- dial_prep_engine(dial_raw_spec(rw$df, "Y", "D", "Z", rw$ctr))
    g <- dial_rstar_gate(eng, c(0, 0.8), alpha0 = ALPHA0, B = B_CD, seed = 7200L + k,
                         step = 0.001)
    outD[[k]] <- data.frame(analysis = rw$lab, r_star = g$r_star, cs_lo = g$cs_lo,
                            cs_hi = g$cs_hi, V = g$V, p_D = g$p_D, dist = g$dist)
    cat(sprintf("  %-26s r* = %+.4f  CS = [%+.4f, %+.4f]  V = %d  p_D = %.3f\n",
                rw$lab, g$r_star, g$cs_lo, g$cs_hi, g$V, g$p_D))
  }
  utils::write.csv(do.call(rbind, outD), file.path(RESULTS_DIR, "rstar_gate_sim.csv"),
                   row.names = FALSE)
}

# ==============================================================================
# B: power curves on the DGP of 07
# ==============================================================================

if ("B" %in% PARTS) {
  cat("\nPart B: power curves (design of 07, population sampling)\n")
  outB <- list()
  for (n_val in N_B) for (di in seq_along(DELTA_B)) {
    delta <- DELTA_B[di]
    seeds <- 30000L * (n_val %% 9973L) + 100L * di + seq_len(R_B)
    t0 <- Sys.time()
    rr <- par_lapply(seeds, function(sd) {
      df <- sim_pop(n_val, 0.50, 0.25, delta, 0.30, sd)
      g  <- gate_xyz(df$y, df$x, df$z, c(0, 0.8), B_B, sd)
      c(V = g$V, rstar = g$r_star)
    })
    M <- do.call(rbind, rr)
    outB[[length(outB) + 1L]] <- data.frame(n = n_val, delta = delta, R = nrow(M),
                                            reject = mean(M[, "V"] == 0),
                                            mean_rstar = mean(M[, "rstar"]))
    cat(sprintf("  n = %5d  delta = %.3f  reject %.3f  mean r* %+.3f  (%.1f min)\n",
                n_val, delta, mean(M[, "V"] == 0), mean(M[, "rstar"]),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    utils::write.csv(do.call(rbind, outB), file.path(RESULTS_DIR, "rstar_gate_power.csv"),
                     row.names = FALSE)
  }
}

cat("\nDONE.\n")
