# =============================================================================
# 08_figures.R
#
# Produces all paper figures and saves them to paper/figures inside the
# scripts folder (the folder that holds this script):
#   fig:architecture     DIAL architecture diagram -> fig1_dial_architecture.png
#   fig:power_curves     Power curves (3 panels)   -> fig2_combined.png
#   fig:domain_sens      Domain sensitivity        -> fig3_domain_sensitivity.png
#   fig:diagnostic_space DIAL diagnostic space     -> fig4_diagnostic_space.png
#
# Figures are named by LaTeX label, not by number, so a float move in the
# paper cannot make these comments stale.
#
# fig:power_curves requires results/cp_grid_summary_lean.csv produced by
# scripts/07_power_curve_grid.R; the script looks for it in results/ inside
# or beside the scripts folder (see Paths below). All other figures are
# self-contained.
# The script prints the 50% crossings delta*_n that Section 4.3 and the
# fig:power_curves caption quote.
#
# NOTE ON THE DOMAIN-SENSITIVITY AND DIAGNOSTIC-SPACE DATA
#   Both figures use hand-coded data frames. Every value must match the
#   empirical tables (tab:empirical_full, tab:domain_sens):
#
#     fig3 rows  -> canonical and within-sign wide MCUB sets.
#                   Phase 1 rows come from 02_phase1_empirical.R
#                   (domain_label "canonical" / "wide_canonical");
#                   Phase 2 rows (Card, Papke) from 03_phase2_empirical.R.
#                   Only instruments with a robustness pair appear.
#                   MovieLens and NHANES have no pair (their canonical D is
#                   already the widest admissible subinterval).
#
#     fig4 rows  -> canonical MCUB endpoints for all 13 instruments plus
#                   R = |beta_IV| / |beta_DML| for the five that pass the
#                   validity gate. Algorithm 1 stops before Stage 3 for
#                   Scenario C (relevance or validity gate), so R is NOT
#                   computed for the eight Scenario C instruments; they
#                   are placed by signed distance only.
#
#   If any value in Phase 1 or Phase 2 shifts, update the data frames
#   below and re-run this script.
#
# Usage: Rscript scripts/08_figures.R, or source() it. To run it line by
# line (for example with Ctrl+Enter in RStudio), first set the working
# directory to the scripts folder.
#
# Runtime: under 2 minutes (excluding fig:power_curves if the CSV is missing).
# =============================================================================

# ---- Paths ------------------------------------------------------------------
# The scripts folder is the folder that holds this script. When the script
# runs line by line, its own path is unknown, and the working directory
# stands in for the scripts folder.
.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) && nchar(f)) dirname(normalizePath(f))
    else getwd()
  }
)
if (!length(.script_dir) || is.na(.script_dir[1])) .script_dir <- getwd()
.repo_root <- normalizePath(file.path(.script_dir, ".."), mustWork = FALSE)

# Figures go to paper/figures inside the scripts folder.
FIG_DIR <- file.path(.script_dir, "paper", "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
FIG_DIR <- normalizePath(FIG_DIR)
cat("Saving figures to:", FIG_DIR, "\n")

# ---- Dependencies -----------------------------------------------------------
install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org/")
}
install_if_missing(c("ggplot2", "ggrepel", "patchwork", "dplyr", "scales"))
suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel)
  library(patchwork); library(dplyr)
})

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;    filter    <- dplyr::filter
mutate    <- dplyr::mutate;    slice     <- dplyr::slice
recode    <- dplyr::recode;    rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange   <- dplyr::arrange;   count     <- dplyr::count
lag       <- dplyr::lag

# =============================================================================
# fig:architecture -- DIAL Architecture Diagram
# =============================================================================

cat("Generating fig:architecture (DIAL architecture diagram)...\n")

# 2800 x 1715 px at 300 dpi: the paper prints the figure at 0.9 textwidth,
# so the labels stay legible (no title; the caption carries it)
png(file.path(FIG_DIR, "fig1_dial_architecture.png"),
    width = 2800, height = 1715, res = 300)

rrect <- function(x, y, w, h, fill, bdr, lwd = 1.5)
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = fill, border = bdr, lwd = lwd)

diam <- function(cx, cy, w, h, fill, bdr, lwd = 1.8)
  polygon(c(cx, cx + w/2, cx, cx - w/2),
          c(cy + h/2, cy, cy - h/2, cy),
          col = fill, border = bdr, lwd = lwd)

arr <- function(x1, y1, x2, y2, col = "gray35", lwd = 1.4)
  arrows(x1, y1, x2, y2, length = 0.08, angle = 22, col = col, lwd = lwd)

arr_L_HV <- function(x1, y1, x2, y2, col = "gray35", lwd = 1.4) {
  segments(x1, y1, x2, y1, col = col, lwd = lwd)
  arrows(x2, y1, x2, y2, length = 0.08, angle = 22, col = col, lwd = lwd)
}

arr_L_VH <- function(x1, y1, x2, y2, col = "gray35", lwd = 1.4) {
  segments(x1, y1, x1, y2, col = col, lwd = lwd)
  arrows(x1, y2, x2, y2, length = 0.08, angle = 22, col = col, lwd = lwd)
}

cA  <- "#7FC97F"; cAb <- "#3D8A3D"
cB  <- "#7FAEDC"; cBb <- "#2A6FB7"
cC  <- "#F4A09A"; cCb <- "#C9433D"
cBg <- "#EFEFE6"; cBgb <- "#8C8C82"
cE  <- "#D6E7F8"; cEb <- "#3F77B5"; cEt <- "#1B4C82"
cR  <- "#E1DDF5"; cRb <- "#6E64C2"; cRt <- "#3C3489"
cT  <- "#D8EFE5"; cTb <- "#3FA083"; cTt <- "#0F6E56"
cF  <- "#FCE3BC"; cFb <- "#E08826"; cFt <- "#9C5A0E"

par(mar = c(0.1, 0.1, 0.1, 0.1), family = "sans")
plot.new()
plot.window(xlim = c(0, 88), ylim = c(26, 96.5), xaxs = "i", yaxs = "i")

xData <- 50;  yData <- 94
xDML  <- 22;  xFS <- 70;  yEst <- 85.5
xF    <- 70;  yF   <- 76.5
xCR   <- 50;  yCR  <- 67
xIV2  <- 50;  yIV2 <- 58.5
xRat  <- 22;  yRat <- 49.5
xR    <- 22;  yR   <- 41.5
xA <- 9;  xB <- 35;  xC <- 78;  yOut <- 31

rrect(xData, yData, 28, 3.2, cBg, cBgb, lwd = 1.3)
text(xData, yData, "Observed data  (Y, D, Z, X)", col = "gray15", cex = 1.0, font = 2)

arr_L_HV(xData - 14, yData - 1.8, xDML, yEst + 2.5)
arr_L_HV(xData + 14, yData - 1.8, xFS,  yEst + 2.5)

rrect(xDML, yEst, 15, 4.5, cE, cEb, lwd = 1.3)
text(xDML, yEst + 1,   "DML",   col = cEt, cex = 1.0, font = 2)
text(xDML, yEst - 1.2, "(partially linear)", col = cEt, cex = 0.72)

arr(xDML, yEst - 2.5, xDML, yRat + 2.75)

rrect(xFS, yEst, 15, 4.5, cE, cEb, lwd = 1.3)
text(xFS, yEst + 1,   "First stage", col = cEt, cex = 0.92, font = 2)
text(xFS, yEst - 1.2, "D ~ Z, X",   col = cEt, cex = 0.72)

arr(xFS, yEst - 2.5, xF, yF + 3.25)

diam(xF, yF, 15, 6, cF, cFb)
text(xF, yF + 0.7, "Weak IV test", col = cFt, cex = 0.78, font = 2)
text(xF, yF - 1.3, "F < 10 ?",    col = cFt, cex = 0.80, font = 2)

arr_L_HV(xF - 7.5, yF, xCR, yCR + 3.75)
text((xF - 7.5 + xCR) / 2, yF + 1.3, "No", col = "gray30", cex = 0.72, font = 2)

arr_L_HV(xF + 7.5, yF, xC + 5, yOut + 4)
text(xF + 10, yF + 1.3, "Yes", col = cCb, cex = 0.72, font = 2)

# Validity gate: the bootstrap rule of Algorithm 1, Stage 2 (share of
# resampled MCUB sets that contain zero at least tau_V), not a single set
diam(xCR, yCR, 20, 7, cR, cRb)
text(xCR, yCR + 1.0, "CR test (MCUB)", col = cRt, cex = 0.86, font = 2)
text(xCR, yCR - 1.4, expression(hat(p)[B] >= tau[V] * "?"), col = cRt, cex = 0.86)

arr_L_HV(xCR + 10, yCR, xC - 6, yOut + 4)
text((xCR + 10 + xC - 6) / 2, yCR + 1.3, "Fail", col = cCb, cex = 0.72, font = 2)

arr(xCR, yCR - 3.75, xIV2, yIV2 + 2.5)
text(xCR + 4, (yCR - 3.75 + yIV2 + 2.5) / 2, "Pass", col = cAb, cex = 0.72, font = 2)

rrect(xIV2, yIV2, 20, 4.5, cE, cEb, lwd = 1.3)
text(xIV2, yIV2 + 1,
     expression(hat(beta)[IV] ~ "(2SLS)"), col = cEt, cex = 0.86, font = 2)
text(xIV2, yIV2 - 1.3, "(LATE)", col = cEt, cex = 0.72)

arr_L_VH(xIV2, yIV2 - 2.5, xRat + 10, yRat)

rrect(xRat, yRat, 20, 5, cT, cTb, lwd = 1.3)
text(xRat, yRat + 1.1, "Effect ratio (R)", col = cTt, cex = 0.86, font = 2)
text(xRat, yRat - 1.3,
     expression(group("|", hat(beta)[IV] / hat(beta)[DML], "|")),
     col = cTt, cex = 0.78)

arr(xRat, yRat - 2.75, xR, yR + 4)

# Ratio gate: Algorithm 1, Stage 3, symmetric (Yes -> A, No -> B). The
# caption adds the sign gate: estimates of opposite sign also go to B.
diam(xR, yR, 18, 7, cF, cFb)
text(xR, yR + 0.2, expression(group("|", log ~ R, "|") < log ~ tau[R] * "?"),
     col = cFt, cex = 0.78, font = 2)

arr_L_HV(xR - 8, yR, xA, yOut + 4)
text(xR - 10.5, yR + 1.5, "Yes", col = cAb, cex = 0.72, font = 2)

arr_L_HV(xR + 8, yR, xB, yOut + 4)
text(xR + 10.5, yR + 1.5, "No",  col = cBb, cex = 0.72, font = 2)

rrect(xA, yOut, 14, 7, cA, cAb)
text(xA, yOut + 2.1, "A",               col = "white",   cex = 1.3, font = 2)
text(xA, yOut - 0.3, "Valid IV",  col = "#0F4F0F", cex = 0.72, font = 2)
text(xA, yOut - 2.3, expression("IV" %~~% "DML"), col = "#0F4F0F", cex = 0.72)

rrect(xB, yOut, 13, 8, cB, cBb)
text(xB, yOut + 2.1, "B",     col = "white",   cex = 1.3, font = 2)
text(xB, yOut - 0.3, "Valid IV",  col = "#0E2E55", cex = 0.72, font = 2)
text(xB, yOut - 2.3, expression("IV" != "DML"), col = "#0E2E55", cex = 0.72)

rrect(xC, yOut, 15, 7, cC, cCb)
text(xC, yOut + 2.1, "C",             col = "white",   cex = 1.3, font = 2)
text(xC, yOut - 0.3, "IV unreliable",  col = "#5C1612", cex = 0.72, font = 2)
text(xC, yOut - 2.3, "Weak or invalid", col = "#5C1612", cex = 0.72)

dev.off()
cat("Saved: fig1_dial_architecture.png\n")

# =============================================================================
# fig:power_curves -- Power curves (3 panels)
# Requires results/cp_grid_summary_lean.csv from 07_power_curve_grid.R
# =============================================================================

cat("Generating fig:power_curves...\n")

# Look for the CSV in results/ one level above the scripts folder (where 07
# writes it), in results/ inside the scripts folder, and in or near the
# working directory. If several copies exist, read the most recently written
# one and list them all, so an old copy is not read by mistake.
csv_name  <- "cp_grid_summary_lean.csv"
csv_tried <- unique(normalizePath(c(
  file.path(.repo_root,  "results", csv_name),
  file.path(.script_dir, "results", csv_name),
  file.path(getwd(), "results", csv_name),
  file.path(getwd(), "..", "results", csv_name),
  file.path(getwd(), csv_name)
), mustWork = FALSE))
csv_found <- csv_tried[file.exists(csv_tried)]
csv_found <- csv_found[order(file.info(csv_found)$mtime, decreasing = TRUE)]
csv_path  <- if (length(csv_found)) csv_found[1] else NA_character_
.stamp    <- function(f) format(file.info(f)$mtime, "%Y-%m-%d %H:%M")
if (length(csv_found) > 1) {
  cat("  Found", length(csv_found), "copies of", csv_name, "(newest first):\n")
  for (f in csv_found) cat("    ", .stamp(f), " ", f, "\n", sep = "")
}

if (is.na(csv_path)) {
  warning("cp_grid_summary_lean.csv not found; run 07_power_curve_grid.R first. ",
          "Skipping fig:power_curves. Looked in:\n  ",
          paste(csv_tried, collapse = "\n  "))
} else {
  cat("  Reading ", csv_path, " (written ", .stamp(csv_path), ")\n", sep = "")
  dat <- read.csv(csv_path, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d cells from cp_grid_summary_lean.csv\n", nrow(dat)))

  ns     <- sort(unique(dat$n))
  n_labs <- paste0("n=", formatC(ns, format = "d", big.mark = ","))

  dat$n_label <- factor(
    paste0("n=", formatC(dat$n, format = "d", big.mark = ",")),
    levels = n_labs
  )

  NCOLORS <- c(
    "n=1,107"  = "#1F77B4",
    "n=1,534"  = "#2CA02C",
    "n=5,000"  = "#D5A03A",
    "n=10,000" = "#D62728",
    "n=50,000" = "#9467BD"
  )

  # ---- Power curves for panels a and b -------------------------------------
  # Straight lines between grid points look jagged where the delta grid is
  # coarse: some sample sizes have only one or two grid points inside the
  # drop. The figure therefore joins each sample size's Monte Carlo
  # estimates with a monotone cubic curve (Fritsch-Carlson; method
  # "monoH.FC" of stats::splinefun). The curve passes through every
  # estimate, stays flat where the estimates are flat, and does not
  # overshoot them. delta*_n is where this curve first crosses 0.5.
  # Set CURVE_FIT <- "linear" for straight lines and linear interpolation,
  # as before.
  CURVE_FIT <- "spline"

  dat  <- dat[order(dat$n, dat$delta), ]
  by_n <- split(dat, dat$n_label, drop = TRUE)

  # Coverage as a function of delta, one function per sample size
  funs <- lapply(by_n, function(d) {
    f <- if (CURVE_FIT == "spline")
      splinefun(d$delta, d$cp_raw, method = "monoH.FC")
    else approxfun(d$delta, d$cp_raw)
    function(x) pmin(1, pmax(0, f(x)))
  })

  # delta*_n: where the curve first crosses 0.5
  cross_half <- function(d, f) {
    cp  <- d$cp_raw
    idx <- which(cp[-length(cp)] > 0.5 & cp[-1] <= 0.5)
    if (!length(idx)) return(NA_real_)
    i <- idx[1]
    if (cp[i + 1] == 0.5) return(d$delta[i + 1])
    uniroot(function(x) f(x) - 0.5, c(d$delta[i], d$delta[i + 1]),
            tol = 1e-10)$root
  }

  delta_star <- do.call(rbind, lapply(names(by_n), function(lab) {
    d <- by_n[[lab]]
    data.frame(n = d$n[1], n_label = d$n_label[1],
               delta_star = cross_half(d, funs[[lab]]))
  }))

  # Curve coordinates for a panel that spans [lo, hi]
  curve_data <- function(lo, hi) {
    do.call(rbind, lapply(names(by_n), function(lab) {
      d <- by_n[[lab]]
      x <- seq(max(lo, min(d$delta)), min(hi, max(d$delta)), length.out = 3001)
      data.frame(n_label = d$n_label[1], delta = x, cp = funs[[lab]](x))
    }))
  }

  # Print the crossings that Section 4.3 and the fig:power_curves caption
  # quote. In the Gaussian design of 07, the identified set contains zero
  # whenever delta <= rho_DU * cor(Z, D) = 0.5 * 0.25 = delta*_inf.
  DELTA_INF <- 0.5 * 0.25
  cat(sprintf("  50%% crossings, %s curves (copy into Section 4.3 and the fig:power_curves caption):\n",
              CURVE_FIT))
  for (k in seq_len(nrow(delta_star))) {
    cat(sprintf("    n = %6d   delta*_n = %.3f   excess over delta*_inf = %.3f\n",
                as.integer(delta_star$n[k]), delta_star$delta_star[k],
                delta_star$delta_star[k] - DELTA_INF))
  }

  theme_fig2 <- function() {
    theme_minimal(base_size = 11) +
      theme(
        plot.title       = element_text(size = 12, face = "plain"),
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(5, 5, 5, 5)
      )
  }

  # Panel a. PANEL_A <- "scaled" plots coverage against the scaled distance
  # sqrt(n) * (delta - delta*_inf). Under the Gaussian approximation of
  # Remark 4, CP_n depends on n and delta near the threshold only through
  # this product, so the five curves should nearly coincide if
  # delta*_n - delta*_inf shrinks like 1/sqrt(n). PANEL_A <- "delta" plots
  # the curves against delta from 0.05 to 0.30 instead.
  PANEL_A <- "scaled"

  if (PANEL_A == "scaled") {
    S_LO <- -0.5; S_HI <- 3
    scaled_curves <- do.call(rbind, lapply(names(by_n), function(lab) {
      d  <- by_n[[lab]]; rn <- sqrt(d$n[1])
      x  <- seq(max(min(d$delta), DELTA_INF + S_LO / rn),
                min(max(d$delta), DELTA_INF + S_HI / rn), length.out = 3001)
      data.frame(n_label = d$n_label[1], s = rn * (x - DELTA_INF),
                 cp = funs[[lab]](x))
    }))
    scaled_points   <- dat
    scaled_points$s <- sqrt(scaled_points$n) * (scaled_points$delta - DELTA_INF)
    scaled_points   <- scaled_points[scaled_points$s >= S_LO &
                                       scaled_points$s <= S_HI, ]

    p_a <- ggplot(scaled_curves,
                  aes(x = s, y = cp, colour = n_label, group = n_label)) +
      geom_hline(yintercept = 0.95, linetype = "dotted",
                 colour = "grey50", linewidth = 0.3) +
      geom_hline(yintercept = 0.5, linetype = "dotted",
                 colour = "grey30", linewidth = 0.3) +
      # delta = delta*_inf, where the identified set starts to exclude zero
      geom_vline(xintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.3) +
      annotate("text", x = 0.06, y = 0.06, hjust = 0, vjust = 0,
               label = "delta[infinity]^'*'", parse = TRUE,
               colour = "grey35", size = 3.2) +
      geom_line(linewidth = 0.8) +
      geom_point(data = scaled_points, aes(x = s, y = cp_raw), size = 1.3) +
      scale_x_continuous(
        limits = c(S_LO, S_HI),
        breaks = seq(-0.5, 3, by = 0.5),
        expand = c(0.005, 0)
      ) +
      scale_y_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1.0),
        limits = c(-0.02, 1.02), expand = c(0, 0)
      ) +
      scale_colour_manual(values = NCOLORS, name = "Sample size") +
      labs(
        title = expression("(a) Coverage on the " * sqrt(n) * " scale"),
        x     = expression(sqrt(n) * (delta - delta[infinity]^"*")),
        y     = expression(CP[n](delta))
      ) +
      theme_fig2() +
      guides(colour = "none")  # the figure shows one legend, panel c's
  } else {
    p_a <- ggplot(curve_data(0.05, 0.30),
                  aes(x = delta, y = cp, colour = n_label, group = n_label)) +
      geom_hline(yintercept = 0.95, linetype = "dotted",
                 colour = "grey50", linewidth = 0.3) +
      geom_hline(yintercept = 0.5, linetype = "dotted",
                 colour = "grey30", linewidth = 0.3) +
      geom_vline(data = delta_star,
                 aes(xintercept = delta_star, colour = n_label),
                 linetype = "dashed", linewidth = 0.4, alpha = 0.7,
                 show.legend = FALSE) +
      geom_line(linewidth = 0.8) +
      scale_x_continuous(
        limits = c(0.05, 0.30),
        breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30),
        expand = c(0.005, 0)
      ) +
      scale_y_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1.0),
        limits = c(-0.02, 1.02), expand = c(0, 0)
      ) +
      scale_colour_manual(values = NCOLORS, name = "Sample size") +
      labs(
        title = "(a) CR test coverage of zero under exclusion violations",
        x     = expression(delta == "|" * rho[ZU] * "|"),
        y     = expression(CP[n](delta))
      ) +
      theme_fig2() +
      guides(colour = "none")  # the figure shows one legend, panel c's
  }

  # Panel b: transition zone zoom; points are the Monte Carlo estimates
  zoom_lo <- 0.125; zoom_hi <- 0.200
  dat_zoom <- dat %>% filter(delta >= zoom_lo & delta <= zoom_hi)

  p_b <- ggplot(curve_data(zoom_lo, zoom_hi),
                aes(x = delta, y = cp, colour = n_label, group = n_label)) +
    geom_hline(yintercept = 0.5, linetype = "dotted",
               colour = "grey30", linewidth = 0.3) +
    geom_vline(
      data = delta_star %>%
        filter(!is.na(delta_star) &
                 delta_star >= zoom_lo & delta_star <= zoom_hi),
      aes(xintercept = delta_star, colour = n_label),
      linetype = "dashed", linewidth = 0.4, alpha = 0.7,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(data = dat_zoom, aes(x = delta, y = cp_raw), size = 1.6) +
    scale_x_continuous(
      limits = c(zoom_lo, zoom_hi),
      breaks = seq(0.13, 0.20, by = 0.01),
      expand = c(0.005, 0)
    ) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      limits = c(-0.02, 1.02), expand = c(0, 0)
    ) +
    scale_colour_manual(values = NCOLORS, name = "Sample size") +
    labs(
      title = expression("(b) Transition zone: " * delta * " near " * delta[n]^"*"),
      x     = expression(delta * " (transition zone)"),
      y     = expression(CP[n](delta))
    ) +
    theme_fig2() +
    guides(colour = "none")  # the figure shows one legend, panel c's

  # Panel c: CI width
  p_c <- ggplot(dat, aes(x = delta, y = ci_w_mean,
                         colour = n_label, group = n_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.0) +
    scale_x_log10(
      breaks = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3),
      labels = c("0.001", "0.003", "0.01", "0.03", "0.1", "0.3")
    ) +
    scale_y_continuous(
      limits = c(0.10, 0.40),
      breaks = seq(0.10, 0.40, by = 0.05)
    ) +
    scale_colour_manual(values = NCOLORS, name = "Sample size") +
    labs(
      title = expression("(c) MCUB CI width: " * 1/sqrt(n) * " behaviour"),
      x     = expression(delta == "|" * rho[ZU] * "|"),
      y     = expression("Mean width of " * widehat(C)[MCUB])
    ) +
    theme_fig2()

  # Panels a and b hide their legends, so patchwork collects a single
  # legend (panel c's) and centres it under the whole figure. Three
  # separate legends do not merge, because the panels use different point
  # sizes, and the third one is cut off at the right edge.
  combined <- (p_a | p_b | p_c) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  # bg = "white": theme_minimal() leaves the plot background blank, so
  # ggsave() would otherwise write a transparent PNG
  ggsave(file.path(FIG_DIR, "fig2_combined.png"), combined,
         width = 16, height = 4.8, units = "in", dpi = 300, bg = "white")
  cat("Saved: fig2_combined.png\n")
}

# =============================================================================
# fig:domain_sens -- Domain sensitivity (within-sign protocol)
#
# Each row's (lower, upper) is a canonical or within-sign wide MCUB set, and
# p_hat is its bootstrap share (columns CI_MCUB / CI_Bei and p_hat_B of the
# CSVs from 02 and 03). All bootstraps redraw the raw rows and redo the
# residualise-and-project step (foundation/bootstrap_rule.R, 2026-10-28).
# Values from the run of 02 and 03 on 2 November 2026.
# The verdict is the bootstrap rule p_hat >= tau_V = 0.95, not zero-in-set:
# Card's sets contain zero but fail the rule.
#   Di Tella judgeAlreadyUsedEM  canonical (0.15, 0.25)  [-0.544, +0.454]  02
#                                wide      (0,    0.8)   [-0.614, +0.928]  02
#   Di Tella percJudgeSentToEM   canonical (0,    0.25)  [-0.204, +0.365]  02
#                                wide      (0,    0.8)   [-0.202, +0.873]  02
#   Galiani highnumber           canonical (-0.2,  0)    [-0.194, +0.016]  02
#                                wide      (-0.8,  0)    [-0.795, +0.016]  02
#   Banerjee instru              canonical (-0.6, -0.4)  [-0.162, +0.224]  02
#                                wide      (-0.8,  0)    [-0.433, +0.596]  02
#   Burde buildschool            canonical (0,    0.4)   [+0.146, +0.736]  02
#                                wide      (0,    0.8)   [+0.146, +0.982]  02
#   Card fatheduc                canonical (0,    0.4)   [-0.071, +0.584]  03
#                                wide      (0,    0.8)   [-0.051, +0.904]  03
#   Papke sole_plan              canonical (0,    0.3)   [+0.101, +0.964]  03
#                                wide      (0,    0.8)   [+0.101, +1.000]  03
#   Papke plan_age               canonical (0,    0.3)   [+0.663, +1.000]  03
#                                wide      (0,    0.8)   [+0.685, +1.000]  03
#   (upper endpoints 1.082, 1.021 and 1.023 before clipping)
# Rows with domain_label "cross_sign_sensitivity" in phase1_results.csv
# are NOT a robustness check and are deliberately excluded.
# The row labels call the wide domain "Robustness", as the paper does.
# =============================================================================

cat("Generating fig:domain_sens (domain sensitivity)...\n")

cValid   <- "#3D8A3D"
cInvalid <- "#C9433D"

# Bootstrap shares p_hat_B (B = 500; column p_hat_B of the 02 and 03
# outputs). _c = canonical domain, _w = within-sign wide domain.
TAU2 <- 0.95
PHAT <- list(
  jaue_c = 0.958, jaue_w = 0.998, pjse_c = 0.996, pjse_w = 0.988,
  gal_c  = 1.000, gal_w  = 1.000, ban_c  = 1.000, ban_w  = 1.000,
  bur_c  = 0.006, bur_w  = 0.004, card_c = 0.814, card_w = 0.798,
  sole_c = 0.284, sole_w = 0.258, age_c  = 0.000, age_w  = 0.000,
  inc    = 0.000, edu    = 0.062, genres = 1.000, year   = 0.016,
  sequel = 0.000)

theme_dial <- function() {
  theme_minimal(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle    = element_text(color = "gray30", size = 9),
      axis.title.x     = element_text(size = 9),
      axis.title.y     = element_text(size = 9),
      axis.text.y      = element_text(size = 7.5, lineheight = 0.95),
      axis.text.x      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 9),
      plot.margin      = margin(10, 14, 8, 8)
    )
}

fig3 <- data.frame(
  row_label = c(
    "Di Tella judgeAlreadyUsedEM\nRobustness (0, 0.8)",
    "Di Tella judgeAlreadyUsedEM\nCanonical  (0.15, 0.25)",
    "Di Tella percJudgeSentToEM\nRobustness (0, 0.8)",
    "Di Tella percJudgeSentToEM\nCanonical  (0, 0.25)",
    "Galiani highnumber\nRobustness (-0.8, 0)",
    "Galiani highnumber\nCanonical  (-0.2, 0)",
    "Banerjee instru\nRobustness (-0.8, 0)",
    "Banerjee instru\nCanonical  (-0.6, -0.4)",
    "Burde buildschool\nRobustness (0, 0.8)",
    "Burde buildschool\nCanonical  (0, 0.4)",
    "Card fatheduc\nRobustness (0, 0.8)",
    "Card fatheduc\nCanonical  (0, 0.4)",
    "Papke sole_plan\nRobustness (0, 0.8)",
    "Papke sole_plan\nCanonical  (0, 0.3)",
    "Papke plan_age\nRobustness (0, 0.8)",
    "Papke plan_age\nCanonical  (0, 0.3)"
  ),
  lower = c(-0.614, -0.544, -0.202, -0.204, -0.795, -0.194,
            -0.433, -0.162,  0.146,  0.146, -0.051, -0.071,
             0.101,  0.101,  0.685,  0.663),
  upper = c( 0.928,  0.454,  0.873,  0.365,  0.016,  0.016,
             0.596,  0.224,  0.982,  0.736,  0.904,  0.584,
             1.000,  0.964,  1.000,  1.000),
  # bootstrap shares p_hat_B (B = 500), same row order
  p_hat = c(PHAT$jaue_w, PHAT$jaue_c, PHAT$pjse_w, PHAT$pjse_c,
            PHAT$gal_w,  PHAT$gal_c,  PHAT$ban_w,  PHAT$ban_c,
            PHAT$bur_w,  PHAT$bur_c,  PHAT$card_w, PHAT$card_c,
            PHAT$sole_w, PHAT$sole_c, PHAT$age_w,  PHAT$age_c),
  stringsAsFactors = FALSE
)
fig3$verdict <- ifelse(fig3$p_hat >= TAU2, "Valid", "Invalid")
fig3$id <- seq_len(nrow(fig3))

# Guard: every pair must share a verdict (the within-sign protocol claim).
stopifnot(all(fig3$verdict[c(TRUE, FALSE)] == fig3$verdict[c(FALSE, TRUE)]))
# Guard: every Valid set contains zero (the rule never passes a set that
# excludes zero in these data); Invalid sets may contain zero (Card).
stopifnot(all(fig3$lower[fig3$verdict == "Valid"] <= 0 &
              fig3$upper[fig3$verdict == "Valid"] >= 0))

fig3$row_factor <- factor(fig3$row_label, levels = rev(fig3$row_label))

p3 <- ggplot(fig3, aes(y = row_factor)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray40", linewidth = 0.4) +
  geom_segment(
    aes(x = lower, xend = upper, yend = row_factor, color = verdict),
    linewidth = 3.2, lineend = "round"
  ) +
  geom_point(aes(x = lower, color = verdict), size = 1.8) +
  geom_point(aes(x = upper, color = verdict), size = 1.8) +
  annotate("text", x = 0.02, y = 0.55,
           label = expression(rho[ZU] == 0),
           color = "gray30", size = 2.9, hjust = 0, vjust = 0) +
  scale_color_manual(
    values = c("Valid" = cValid, "Invalid" = cInvalid),
    breaks = c("Valid", "Invalid")
  ) +
  scale_x_continuous(
    limits = c(-0.85, 1.00),
    breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0)
  ) +
  labs(
    x        = expression("MCUB set for " * rho[ZU]),
    y        = NULL,
    title    = "Domain Sensitivity: MCUB Sets (Within-Sign Protocol)",
    subtitle = "Each instrument shown under two within-sign domains; no verdict flips within the elicited half-line",
    color    = "Verdict"
  ) +
  theme_dial() +
  guides(color = guide_legend(override.aes = list(linewidth = 5, size = 0)))

ggsave(file.path(FIG_DIR, "fig3_domain_sensitivity.png"), p3,
       width = 8.2, height = 8.6, dpi = 300, bg = "white")
cat("Saved: fig3_domain_sensitivity.png\n")

# =============================================================================
# fig:diagnostic_space -- DIAL diagnostic space (13 instruments, two panels)
#
# Both panels use the signed distance d from zero to the nearest endpoint of
# the canonical MCUB set: d = min(-lower, upper), which is the depth of zero
# inside the set when lower <= 0 <= upper (d >= 0) and minus the gap to the
# set otherwise (d < 0).
#
#   Left panel  -- the eight Scenario C instruments, placed by d alone.
#                  Algorithm 1 stops at the relevance or validity gate, so
#                  R is not computed for them. Open triangles mark the two
#                  weak instruments (first-stage F < 10), whose MCUB sets
#                  also exclude zero. Card's set contains zero (d > 0), but
#                  it fails the bootstrap rule (p_hat_B < tau_V); an open
#                  circle marks it.
#   Right panel -- the five instruments that pass both gates, placed by d
#                  and R = |beta_IV| / |beta_DML|. Scenario A is the band
#                  1/tau_R < R < tau_R with estimates of the same sign
#                  (dashed lines at 1/1.5 and 1.5); Scenario B lies outside
#                  it, or has estimates of opposite sign (open circle:
#                  MovieLens num_genres).
#
# Canonical MCUB endpoints (tab:empirical_full):
#   MovieLens num_genres      (0, 0.8)      [-0.400, +0.695]  R = 0.62  B   (item level; signs differ)
#   Di Tella judgeAlreadyUsed (0.15, 0.25)  [-0.544, +0.454]  R = 3.13  B
#   Di Tella percJudgeSent    (0, 0.25)     [-0.204, +0.365]  R = 1.37  A
#   Galiani highnumber        (-0.2, 0)     [-0.194, +0.016]  R = 2.44  B
#   Banerjee instru           (-0.6, -0.4)  [-0.162, +0.224]  R = 1.84  B
#   Papke sole_plan           (0, 0.3)      [+0.101, +0.964]            C
#   Papke plan_age            (0, 0.3)      [+0.663, +1.000]            C
#   NHANES income             (0, 0.8)      [+0.511, +1.000]            C
#   NHANES education          (0, 0.8)      [+0.157, +1.000]            C
#   Card fatheduc             (0, 0.4)      [-0.071, +0.584]            C   (fails the rule)
#   Burde buildschool         (0, 0.4)      [+0.146, +0.736]            C
#   MovieLens release_year    (-0.8, 0)     [-1.000, -0.244]            C   (weak: F = 1.7)
#   MovieLens is_sequel       (-0.8, 0)     [-1.000, -0.524]            C   (weak: F = 1.3)
#   MovieLens rows come from run_movielens_item_level() in 03. R comes from
#   the seeded DML fits of 09_implied_endogeneity.R (run of 4 November
#   2026), which also gives its bootstrap intervals (tab:locality).
# =============================================================================

cat("Generating fig:diagnostic_space (DIAL diagnostic space)...\n")

cA <- "#3D8A3D"; cB <- "#2E5EAA"; cC <- "#C9433D"
TAU_R <- 1.5

fig4 <- data.frame(
  label = c(
    "MovieLens num_genres", "Di Tella judgeAlreadyUsedEM",
    "Di Tella percJudgeSentToEM", "Galiani highnumber",
    "Banerjee instru",
    "Papke sole_plan", "Papke plan_age",
    "NHANES income", "NHANES education",
    "Card fatheduc", "Burde buildschool",
    "MovieLens release_year", "MovieLens is_sequel"
  ),
  lower = c(-0.400, -0.544, -0.204, -0.194, -0.162,
             0.101,  0.663,  0.511,  0.157, -0.071,  0.146,
            -1.000, -1.000),
  upper = c( 0.695,  0.454,  0.365,  0.016,  0.224,
             0.964,  1.000,  1.000,  1.000,  0.584,  0.736,
            -0.244, -0.524),
  p_hat = c(PHAT$genres, PHAT$jaue_c, PHAT$pjse_c, PHAT$gal_c, PHAT$ban_c,
            PHAT$sole_c, PHAT$age_c, PHAT$inc, PHAT$edu, PHAT$card_c,
            PHAT$bur_c, PHAT$year, PHAT$sequel),
  ratio = c(0.62, 3.13, 1.37, 2.44, 1.84, rep(NA_real_, 8)),
  # IV and DML estimates of opposite sign (sign gate of Algorithm 1)
  sign_diff = c(TRUE, FALSE, FALSE, FALSE, FALSE, rep(NA, 8)),
  dial  = c("B", "B", "A", "B", "B", rep("C", 8)),
  weak  = c(rep(FALSE, 11), TRUE, TRUE),
  stringsAsFactors = FALSE
)
# Two-line labels for the right panel
fig4$label2 <- sub("^(Di Tella|MovieLens|Galiani|Banerjee) ", "\\1\n", fig4$label)
fig4$d <- pmin(-fig4$lower, fig4$upper)
# Scenario C instruments whose set contains zero but which fail the rule
fig4$rule_only <- fig4$dial == "C" & fig4$d >= 0

# Guards: the verdict follows Algorithm 1 (relevance, then the bootstrap
# rule, then the ratio), and R exists exactly for the instruments that
# pass both gates.
stopifnot(all((fig4$dial != "C") == (!fig4$weak & fig4$p_hat >= TAU2)))
stopifnot(all(fig4$d[fig4$dial != "C"] >= 0))
stopifnot(all(is.na(fig4$ratio) == (fig4$dial == "C")))
v4 <- fig4$dial != "C"
stopifnot(all(fig4$dial[v4] ==
                ifelse(fig4$sign_diff[v4] | abs(log(fig4$ratio[v4])) >= log(TAU_R),
                       "B", "A")))

valid4   <- fig4[fig4$dial != "C", ]
valid4$label2 <- ifelse(valid4$sign_diff, paste0(valid4$label2, "\n(signs differ)"),
                        valid4$label2)
valid4$band <- ifelse(valid4$ratio >= TAU_R, "high",
                      ifelse(valid4$ratio <= 1 / TAU_R, "low", "mid"))
invalid4 <- fig4[fig4$dial == "C", ]
invalid4$name <- ifelse(invalid4$weak,
                        paste0(invalid4$label, " (weak)"), invalid4$label)
invalid4$name <- ifelse(invalid4$rule_only,
                        sprintf("%s (share %.2f)", invalid4$label, invalid4$p_hat),
                        invalid4$name)
invalid4$mark <- ifelse(invalid4$weak, "weak", ifelse(invalid4$rule_only, "rule", "out"))
# Closest to zero at the top
invalid4$name <- factor(invalid4$name,
                        levels = invalid4$name[order(invalid4$d)])

# ---- Left panel: Scenario C, placed by d alone ------------------------------
pC <- ggplot(invalid4, aes(y = name)) +
  # Shade the whole panel as Scenario C, as the right panel shades A and B.
  # The dashed line at d = 0 still shows that Card's set contains zero.
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = cC, alpha = 0.06) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_segment(aes(x = d, xend = 0, yend = name),
               colour = cC, alpha = 0.35, linewidth = 0.7) +
  geom_point(aes(x = d, shape = mark), colour = cC, size = 2.6, stroke = 0.9) +
  scale_shape_manual(values = c(out = 17, weak = 2, rule = 1), guide = "none") +
  scale_x_continuous(limits = c(-0.85, 0.12), breaks = c(-0.75, -0.5, -0.25, 0),
                     labels = c("-0.75", "-0.50", "-0.25", "0"),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title    = "C: unreliable (invalid or weak)",
    subtitle = "F < 10 or share < 0.95; R not computed",
    x        = expression("Signed distance " * d),
    y        = NULL
  ) +
  theme_dial() +
  theme(
    plot.title         = element_text(face = "bold", colour = cC, size = 10),
    plot.subtitle      = element_text(colour = "gray35", size = 7.5, face = "italic"),
    axis.text.y        = element_text(size = 8, colour = "gray15"),
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(4, 8, 4, 4)
  )

# ---- Right panel: Scenarios A and B, placed by d and R ----------------------
pAB <- ggplot(valid4, aes(x = d, y = ratio, colour = dial)) +
  # Scenario A: the band 1/tau_R < R < tau_R; Scenario B: outside it
  annotate("rect", xmin = 0, xmax = Inf, ymin = 1 / TAU_R, ymax = TAU_R,
           fill = cA, alpha = 0.08) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = TAU_R, ymax = Inf,
           fill = cB, alpha = 0.06) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = 1 / TAU_R,
           fill = cB, alpha = 0.06) +
  geom_hline(yintercept = c(1 / TAU_R, TAU_R), linetype = "dashed",
             colour = "gray50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  annotate("text", x = 0.49, y = 1 / TAU_R + 0.10, hjust = 1, vjust = 0,
           label = "A: valid, IV and DML agree",
           colour = cA, fontface = "bold", size = 3) +
  annotate("text", x = 0.49, y = 4.32, hjust = 1, vjust = 1,
           label = "B: valid, IV and DML disagree",
           colour = cB, fontface = "bold", size = 3) +
  annotate("text", x = 0.49, y = TAU_R + 0.08, hjust = 1, vjust = 0,
           label = "tau[R] == 1.5", parse = TRUE,
           colour = "gray35", size = 2.6) +
  annotate("text", x = 0.02, y = 1 / TAU_R - 0.06, hjust = 0, vjust = 1,
           label = "1 / tau[R]", parse = TRUE,
           colour = "gray35", size = 2.6) +
  geom_point(aes(shape = ifelse(sign_diff, "sign", "same")), size = 2.8, stroke = 1) +
  scale_shape_manual(values = c(same = 16, sign = 1), guide = "none") +
  # Labels stay in their own band
  # The A band is narrow: put the label right beside its point
  geom_text(
    data = valid4[valid4$band == "mid", ], aes(x = d + 0.025, label = label2),
    hjust = 0, size = 2.5, fontface = "bold", lineheight = 0.85,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = valid4[valid4$band == "low", ], aes(label = label2),
    size = 2.5, fontface = "bold", lineheight = 0.85,
    box.padding = 0.4, point.padding = 0.3, force = 3,
    min.segment.length = 0.1, segment.color = "gray60", segment.size = 0.3,
    xlim = c(0.03, 0.36), ylim = c(0.02, 1 / TAU_R - 0.05),
    max.overlaps = Inf, seed = 17
  ) +
  geom_text_repel(
    data = valid4[valid4$band == "high", ], aes(label = label2),
    size = 2.5, fontface = "bold", lineheight = 0.85,
    box.padding = 0.5, point.padding = 0.35, force = 3,
    min.segment.length = 0.1, segment.color = "gray60", segment.size = 0.3,
    # Labels end left of x = 0.43, so the long Di Tella label clears its
    # point at d = 0.454
    xlim = c(0.03, 0.43), ylim = c(TAU_R + 0.3, 4.0),
    max.overlaps = Inf, seed = 17
  ) +
  scale_colour_manual(values = c("A" = cA, "B" = cB), guide = "none") +
  scale_x_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
                     expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(limits = c(0, 4.4), breaks = 0:4,
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title    = "A and B: valid",
    # The subtitles state the gates of Algorithm 1 (relevance, then the
    # bootstrap rule); d is only the x-axis. The \u2265 escape keeps this file ASCII.
    subtitle = "F \u2265 10 and share \u2265 0.95; R computed",
    x        = expression("Signed distance " * d),
    y        = expression("|" * hat(beta)[IV] / hat(beta)[DML] * "|")
  ) +
  theme_dial() +
  theme(
    plot.title    = element_text(face = "bold", colour = "gray20", size = 10),
    plot.subtitle = element_text(colour = "gray35", size = 7.5, face = "italic"),
    plot.margin   = margin(4, 6, 4, 8)
  )

p4 <- (pC | pAB) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title    = "DIAL Diagnostic Space: 13 Instruments Across 8 Studies",
    subtitle = "5 of 13 pass both gates under the canonical domain D",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(colour = "gray30", size = 9)
    )
  )

ggsave(file.path(FIG_DIR, "fig4_diagnostic_space.png"), p4,
       width = 6.8, height = 3.8, dpi = 300, bg = "white")
cat("Saved: fig4_diagnostic_space.png\n")

cat("\nAll figures saved to", FIG_DIR, "\n")
