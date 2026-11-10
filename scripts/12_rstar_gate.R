# ==============================================================================
# 12_rstar_gate.R  -  the validity gate for the 13 empirical instruments
#                     (tab:empirical_full, tab:locality, app:frontier)
# ==============================================================================
#
# PURPOSE
#   Applies the validity gate of Algorithm 1 to each instrument under its
#   canonical domain D: the confidence set for the implied endogeneity r*
#   obtained by inverting the CR mapping (foundation/rstar_gate.R), and the
#   verdict V = 1 when the set meets the closed domain. The test has level
#   ALPHA0 = 0.05, so the set has coverage 1 - 2 ALPHA0 = 0.90.
#
#   For each instrument the script reports r*, the confidence set, V, the
#   p-value p_D (the largest one-sided p-value over the closed domain), the
#   breakdown distance from r* to D, and, for an instrument whose r* lies
#   outside D, the distance kappa = |g(rho, e)| from zero to the identified
#   set at the nearest edge e, its scale s = sqrt(n) se(e), the gate's own
#   frontier kappa*_N = (z_{1-alpha0} + z_{1-pi*}) s / sqrt(n), and the
#   Gaussian power prediction. It also reports the verdict of the
#   MCUB-share rule (p_hat_B >= 0.95, from 02 and 03; rem:mcub_margin) for
#   comparison.
#
# INPUTS
#   the data of dial_specs.R; scripts/phase1_results.csv and
#   scripts/phase2_results.csv (02, 03) for the MCUB comparison
#
# OUTPUT
#   results/rstar_gate.csv and LaTeX rows on the console
#
# RUNTIME  about half a minute (B = 1,000 resamples per instrument).
#
# RUN (from the scripts folder)
#   Rscript 12_rstar_gate.R
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({ library(dplyr); library(haven) })

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

ML_DIR     <- if (dir.exists("ml-100k")) "ml-100k" else NULL
NHANES_CSV <- "nhanes_mercury.csv"
ALPHA0  <- 0.05      # level of the gate; the confidence set has coverage 0.90
B_GATE  <- 1000L     # bootstrap resamples for se(r)
STEP    <- 0.001     # grid step for r
PI_STAR <- 0.05      # frontier: detect with probability at least 0.95
F_MIN   <- 10        # relevance gate

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
source(file.path(FOUNDATION_DIR, "rstar_gate.R"))
.f03 <- file.path(.script_dir, "03_phase2_empirical.R")
if (file.exists(.f03)) {
  .wd <- setwd(.script_dir)
  tryCatch(invisible(capture.output(source(.f03))),
           error = function(e) cat("NOTE: could not source 03 (MovieLens reader):",
                                   conditionMessage(e), "\n"))
  setwd(.wd)
}
source(file.path(.script_dir, "dial_specs.R"))

# First-stage F of the residualised regression (Stage 1 of Algorithm 1)
first_stage_F <- function(d, X, Z, controls) {
  f1 <- stats::lm(stats::reformulate(c(Z, controls), X), data = d)
  f0 <- if (length(controls)) stats::lm(stats::reformulate(controls, X), data = d) else
    stats::lm(stats::as.formula(paste(X, "~ 1")), data = d)
  stats::anova(f0, f1)$F[2]
}

# MCUB bootstrap-rule shares of the canonical domain (02, 03)
p1 <- utils::read.csv(file.path(.script_dir, "phase1_results.csv"), stringsAsFactors = FALSE)
p2 <- utils::read.csv(file.path(.script_dir, "phase2_results.csv"), stringsAsFactors = FALSE)
mcub_phat <- function(s) {
  r <- if (s$phase == 1L) p1[p1$instrument == s$instrument & p1$domain_label == "canonical", ] else
    p2[p2$Instrument == s$instrument & p2$domain_label == "canonical", ]
  if (!nrow(r)) NA_real_ else as.numeric(r$p_hat_B[1])
}

# ==============================================================================
# MAIN LOOP
# ==============================================================================

cat("================================================================================\n")
cat(sprintf("r* gate: level %.2f (confidence set %.0f%%), B = %d\n",
            ALPHA0, 100 * (1 - 2 * ALPHA0), B_GATE))
cat("================================================================================\n\n")

rows <- list()
for (s in SPECS) {
  cat(sprintf("--- %s %s (D = [%.2f, %.2f]) ---\n", s$study, s$instrument, s$D[1], s$D[2]))
  d_all <- dial_load_spec_data(s)
  out <- data.frame(study = s$study, instrument = s$instrument, D_lo = s$D[1], D_hi = s$D[2],
                    n = NA_integer_, F_first = NA_real_, r_star = NA_real_,
                    cs_lo = NA_real_, cs_hi = NA_real_, n_segments = NA_integer_,
                    V = NA_integer_, p_D = NA_real_, dist = NA_real_, margin = NA_real_,
                    kappa = NA_real_, s_edge = NA_real_, kappa_star = NA_real_,
                    power_gauss = NA_real_, m = NA_real_, verdict = NA_character_,
                    p_hat_B_mcub = mcub_phat(s), stringsAsFactors = FALSE)
  if (!is.null(d_all)) {
    spec <- dial_raw_spec(d_all, s$Y, s$X, s$Z, s$controls)
    eng  <- dial_prep_engine(spec)
    g    <- dial_rstar_gate(eng, s$D, alpha0 = ALPHA0, B = B_GATE, seed = s$seed, step = STEP)
    out$n <- g$n
    out$F_first <- first_stage_F(spec$data, s$X, s$Z, s$controls)
    out$r_star <- g$r_star; out$cs_lo <- g$cs_lo; out$cs_hi <- g$cs_hi
    out$n_segments <- g$n_segments; out$V <- g$V; out$p_D <- g$p_D
    out$dist <- g$dist
    out$margin <- if (g$dist > 0) -g$dist else min(g$r_star - s$D[1], s$D[2] - g$r_star)
    out$kappa <- g$kappa; out$s_edge <- g$s_edge
    out$kappa_star <- if (is.na(g$s_edge)) NA_real_ else
      dial_kappa_star(g$n, g$s_edge, ALPHA0, PI_STAR)
    out$power_gauss <- g$power_gauss
    out$m <- dial_raw_rho(eng)[["m"]]
    out$verdict <- if (out$F_first < F_MIN) "C (weak)" else if (g$V == 1L) "pass" else "C (invalid)"
    cat(sprintf("    n = %d  F = %.1f  r* = %+.3f  CS = [%+.3f, %+.3f] (%d segment%s)  V = %d  p_D = %.3f  -> %s\n",
                g$n, out$F_first, g$r_star, g$cs_lo, g$cs_hi, g$n_segments,
                if (g$n_segments == 1L) "" else "s", g$V, g$p_D, out$verdict))
    if (g$dist > 0)
      cat(sprintf("    breakdown %.3f  kappa %.3f  s %.2f  kappa*_N %.3f  kappa/kappa* %.1f  Gaussian power %.3f\n",
                  g$dist, g$kappa, g$s_edge, out$kappa_star, g$kappa / out$kappa_star, g$power_gauss))
    cat(sprintf("    MCUB bootstrap rule (02/03): p_hat_B = %.3f -> %s\n", out$p_hat_B_mcub,
                if (isTRUE(out$p_hat_B_mcub >= 0.95)) "pass" else "fail"))
  }
  rows[[length(rows) + 1L]] <- out
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(RESULTS_DIR, "rstar_gate.csv"), row.names = FALSE)

cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")
print(res[, c("study", "instrument", "n", "F_first", "r_star", "cs_lo", "cs_hi", "V", "p_D",
              "dist", "verdict", "p_hat_B_mcub")], digits = 3, row.names = FALSE)
ok <- !is.na(res$V)
cat(sprintf("\nVerdicts: %d pass, %d fail (invalid), %d weak; %d instruments without data.\n",
            sum(res$verdict == "pass", na.rm = TRUE), sum(res$verdict == "C (invalid)", na.rm = TRUE),
            sum(res$verdict == "C (weak)", na.rm = TRUE), sum(!ok)))
chg <- ok & ((res$V == 1L) != (res$p_hat_B_mcub >= 0.95)) & res$F_first >= F_MIN
if (any(chg)) cat("Gate and MCUB bootstrap rule disagree for:",
                  paste(res$instrument[chg], collapse = ", "), "\n")

cat("\nLaTeX rows (instrument & D & r* & 90% set & p_D & distance & verdict):\n\n")
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  if (is.na(r$V)) { cat(sprintf("%% %s %s: data not available\n", r$study, r$instrument)); next }
  cat(sprintf("%s \\texttt{%s} & $(%g,\\,%g)$ & $%+.3f$ & $[%+.3f,\\,%+.3f]$ & $%.3f$ & $%.3f$ & %s \\\\\n",
              r$study, gsub("_", "\\\\_", r$instrument), r$D_lo, r$D_hi, r$r_star,
              r$cs_lo, r$cs_hi, r$p_D, r$dist, r$verdict))
}
cat(sprintf("\nWrote %s\n", file.path(RESULTS_DIR, "rstar_gate.csv")))
cat("\nDONE.\n")
