# =============================================================================
# 08_figures.R
#
# Produces all paper figures and saves them to paper/figures inside the
# scripts folder (the folder that holds this script):
#   fig:architecture     DIAL architecture diagram -> fig1_dial_architecture.png
#   fig:power_curves     Power of the gate         -> fig2_gate_power.png
#   fig:rstar            r* against the domain     -> fig3_rstar_domain.png
#   fig:diagnostic_space DIAL diagnostic space     -> fig4_diagnostic_space.png
#
# Figures are named by LaTeX label, not by number, so a float move in the
# paper cannot make these comments stale.
#
# INPUTS
#   fig:power_curves      results/rstar_gate_power.csv   (13_rstar_gate_mc.R)
#   fig:rstar             results/rstar_gate.csv         (12_rstar_gate.R)
#   fig:diagnostic_space  results/rstar_gate.csv and
#                         results/implied_endogeneity.csv (09: R, its
#                         interval and the B share)
#   fig:architecture is self-contained. The script looks for results/ inside
#   or beside the scripts folder (see Paths below).
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
xCR   <- 50;  yCR  <- 67.5
xIV2  <- 50;  yIV2 <- 58.5
xRat  <- 22                   # y of the effect-ratio box: set below
xR    <- 22;  yR   <- 41.5
xA <- 9;  xB <- 35;  xC <- 78;  yOut <- 31

# Two boxes hold plotmath with hats, subscripts and bars, whose height
# depends on the platform's fonts (Windows fonts are larger than Linux
# ones). Each box is sized to its two lines, so nothing touches its border:
# the first line hangs from the top of the box and the second sits on its
# bottom. A plain string's ink runs below its baseline by about 0.3 of its
# height, so a string on the bottom line is raised by that much.
ink_h <- function(lab, cex, font = 1) {
  h <- strheight(lab, cex = cex, font = font)
  if (is.character(lab)) 1.3 * h else h
}
padB <- 0.55; gapB <- 0.3   # padB leaves room for Windows ink beyond the bounding box
iv_l1  <- expression(hat(beta)[IV] ~ "(2SLS)")
rat_l1 <- "Effect ratio (R)"
rat_l2 <- expression(group("|", hat(beta)[IV] / hat(beta)[DML], "|"))
hIV2 <- max(4.5, 2 * padB + gapB + ink_h(iv_l1, 0.86, 2) + ink_h("(LATE)", 0.72))
hRat <- max(5, 2 * padB + gapB + ink_h(rat_l1, 0.86, 2) + ink_h(rat_l2, 0.78))
# The effect-ratio box keeps its bottom at 47 and grows upward, so the arrow
# to the ratio gate below it keeps its length.
yRat <- 47 + hRat / 2
cat(sprintf("  IV box height %.2f, effect-ratio box height %.2f\n", hIV2, hRat))

rrect(xData, yData, 28, 3.2, cBg, cBgb, lwd = 1.3)
text(xData, yData, "Observed data  (Y, D, Z, X)", col = "gray15", cex = 1.0, font = 2)

arr_L_HV(xData - 14, yData - 1.8, xDML, yEst + 2.5)
arr_L_HV(xData + 14, yData - 1.8, xFS,  yEst + 2.5)

rrect(xDML, yEst, 15, 4.5, cE, cEb, lwd = 1.3)
text(xDML, yEst + 1,   "DML",   col = cEt, cex = 1.0, font = 2)
text(xDML, yEst - 1.2, "(partially linear)", col = cEt, cex = 0.72)

arr(xDML, yEst - 2.5, xDML, yRat + hRat / 2 + 0.25)

rrect(xFS, yEst, 15, 4.5, cE, cEb, lwd = 1.3)
text(xFS, yEst + 1,   "First stage", col = cEt, cex = 0.92, font = 2)
text(xFS, yEst - 1.2, "D ~ Z, X",   col = cEt, cex = 0.72)

# Weak-IV diamond: sized to its two lines, like the gate diamond below, so
# the lines stay inside it whatever fonts the platform uses. The first line
# sits on a baseline just above the centre and the second hangs just below
# it (neither has descenders). The width is the smallest that keeps the
# outer corners of both lines inside the diamond, with a margin.
wk_l1 <- "Weak IV test"; wk_l2 <- "F < 10 ?"
hWk <- 6.5
gWk <- 0.25                                                  # half the gap between the lines
dy1 <- gWk + strheight(wk_l1, cex = 0.78, font = 2) + 0.15   # top of line 1
dy2 <- gWk + strheight(wk_l2, cex = 0.80, font = 2) + 0.15   # bottom of line 2
wWk <- max(15,
           2 * (strwidth(wk_l1, cex = 0.78, font = 2) / 2 + 0.4) / (1 - dy1 / (hWk / 2)),
           2 * (strwidth(wk_l2, cex = 0.80, font = 2) / 2 + 0.4) / (1 - dy2 / (hWk / 2)))
cat(sprintf("  weak-IV diamond: width %.1f\n", wWk))

arr(xFS, yEst - 2.5, xF, yF + hWk / 2 + 0.25)

diam(xF, yF, wWk, hWk, cF, cFb)
text(xF, yF + gWk, wk_l1, col = cFt, cex = 0.78, font = 2, adj = c(0.5, 0))
text(xF, yF - gWk, wk_l2, col = cFt, cex = 0.80, font = 2, adj = c(0.5, 1))

# Validity gate: Algorithm 1, Stage 2 (does the confidence set for r*
# meet the domain D?). The diamond is sized to its label, so the label
# stays inside it whatever fonts the platform uses (Windows fonts are
# wider than Linux ones).
gate_lab <- expression(CS(r^"*") ~ meets ~ D * "?")
hCR <- 8
tw  <- max(strwidth(gate_lab, cex = 0.80), strwidth("CR gate", cex = 0.86, font = 2))
th  <- max(strheight(gate_lab, cex = 0.80), 1.2 * strheight("CS(r) D", cex = 0.80))
dCR <- 1.0 + 0.5 * th + 0.15       # depth of the label's lowest point below the centre
wCR <- max(20, 2 * (tw / 2 + 0.5) / (1 - dCR / (hCR / 2)))
cat(sprintf("  gate diamond: width %.1f (label width %.1f)\n", wCR, tw))

arr_L_HV(xF - wWk / 2, yF, xCR, yCR + hCR / 2 + 0.25)
text((xF - wWk / 2 + xCR) / 2, yF + 1.3, "No", col = "gray30", cex = 0.72, font = 2)

arr_L_HV(xF + wWk / 2, yF, xC + 5, yOut + 4)
text(xF + wWk / 2 + 2.5, yF + 1.3, "Yes", col = cCb, cex = 0.72, font = 2)

diam(xCR, yCR, wCR, hCR, cR, cRb)
text(xCR, yCR + 1.3, "CR gate", col = cRt, cex = 0.86, font = 2)
text(xCR, yCR - 1.0, gate_lab, col = cRt, cex = 0.80)

arr_L_HV(xCR + wCR / 2, yCR, xC - 6, yOut + 4)
text((xCR + wCR / 2 + xC - 6) / 2, yCR + 1.3, "Fail", col = cCb, cex = 0.72, font = 2)

arr(xCR, yCR - hCR / 2 - 0.25, xIV2, yIV2 + hIV2 / 2 + 0.25)
text(xCR + 4, (yCR - hCR / 2 - 0.25 + yIV2 + hIV2 / 2 + 0.25) / 2, "Pass",
     col = cAb, cex = 0.72, font = 2)

rrect(xIV2, yIV2, 20, hIV2, cE, cEb, lwd = 1.3)
text(xIV2, yIV2 + hIV2 / 2 - padB, iv_l1, col = cEt, cex = 0.86, font = 2,
     adj = c(0.5, 1))
text(xIV2, yIV2 - hIV2 / 2 + padB + 0.3 * strheight("(LATE)", cex = 0.72), "(LATE)",
     col = cEt, cex = 0.72, adj = c(0.5, 0))

arr_L_VH(xIV2, yIV2 - hIV2 / 2 - 0.25, xRat + 10, yRat)

rrect(xRat, yRat, 20, hRat, cT, cTb, lwd = 1.3)
text(xRat, yRat + hRat / 2 - padB, rat_l1, col = cTt, cex = 0.86, font = 2,
     adj = c(0.5, 1))
text(xRat, yRat - hRat / 2 + padB, rat_l2, col = cTt, cex = 0.78, adj = c(0.5, 0))

arr(xRat, yRat - hRat / 2 - 0.25, xR, yR + 4)

# Ratio gate: Algorithm 1, Stage 3, symmetric (Yes -> A, No -> B). The
# caption adds the sign gate: estimates of opposite sign also go to B.
diam(xR, yR, 18, 7, cF, cFb)
text(xR, yR + 0.2, expression(group("|", log ~ R, "|") < log ~ tau[R] * "?"),
     col = cFt, cex = 0.78, font = 2)

# The connectors start at the diamond's tips (half-width 9), so they do not
# run into it.
arr_L_HV(xR - 9, yR, xA, yOut + 4)
text(xR - 10.5, yR + 1.5, "Yes", col = cAb, cex = 0.72, font = 2)

arr_L_HV(xR + 9, yR, xB, yOut + 4)
text(xR + 10.5, yR + 1.5, "No",  col = cBb, cex = 0.72, font = 2)

rrect(xA, yOut, 14, 7, cA, cAb)
text(xA, yOut + 2.1, "A",               col = "white",   cex = 1.3, font = 2)
text(xA, yOut - 0.3, "Not rejected",  col = "#0F4F0F", cex = 0.72, font = 2)
text(xA, yOut - 2.3, expression("IV" %~~% "DML"), col = "#0F4F0F", cex = 0.72)

rrect(xB, yOut, 13, 8, cB, cBb)
text(xB, yOut + 2.1, "B",     col = "white",   cex = 1.3, font = 2)
text(xB, yOut - 0.3, "Not rejected",  col = "#0E2E55", cex = 0.72, font = 2)
text(xB, yOut - 2.3, expression("IV" != "DML"), col = "#0E2E55", cex = 0.72)

rrect(xC, yOut, 15, 7, cC, cCb)
text(xC, yOut + 2.1, "C",             col = "white",   cex = 1.3, font = 2)
text(xC, yOut - 0.3, "IV unreliable",  col = "#5C1612", cex = 0.72, font = 2)
text(xC, yOut - 2.3, "Weak or invalid", col = "#5C1612", cex = 0.72)

dev.off()
cat("Saved: fig1_dial_architecture.png\n")

# =============================================================================
# fig:power_curves -- rejection rate of the validity gate (2 panels)
# Requires results/rstar_gate_power.csv from 13_rstar_gate_mc.R (part B).
#
#   (a) rejection rate against delta = cor(Z, U) for each n; dashed line at
#       delta*_inf = rho_DU cor(Z, D) = 0.125, dotted line at alpha0 = 0.05
#   (b) against sqrt(n) (delta - delta*_inf), with the Gaussian prediction
#       Phi(sqrt(n) kappa / s - z_{1 - alpha0}) of prop:gate(d). In this
#       design kappa = (delta - delta*) / sqrt(1 - rho_DU^2) at the edge
#       e = 0; s = s(0) is the asymptotic sd of sqrt(n) g(rho_hat, 0) at
#       delta = delta*, computed here by simulation.
# =============================================================================

cat("Generating fig:power_curves (gate)...\n")
.pw_file <- file.path(.repo_root, "results", "rstar_gate_power.csv")
if (!file.exists(.pw_file)) {
  warning("results/rstar_gate_power.csv not found; run 13_rstar_gate_mc.R first. ",
          "Skipping fig:power_curves.")
} else {
  pw <- utils::read.csv(.pw_file)
  RHO_P <- 0.50; A_P <- 0.25; ALPHA0_P <- 0.05
  DSTAR <- RHO_P * A_P
  # s(0) at delta = delta*: sd of sqrt(n) g(rho_hat, 0) over simulated samples
  set.seed(2027)
  s0 <- local({
    n0 <- 5000L
    c2 <- (DSTAR - A_P * RHO_P) / sqrt(1 - RHO_P^2); c3 <- sqrt(1 - A_P^2 - c2^2)
    g0 <- replicate(2000, {
      Dt <- rnorm(n0); Vt <- rnorm(n0); U <- RHO_P * Dt + sqrt(1 - RHO_P^2) * Vt
      Z <- A_P * Dt + c2 * Vt + c3 * rnorm(n0); Y <- 0.30 * Dt + U
      r <- cor(cbind(Dt, Y, Z))
      (r[2, 3] - r[1, 2] * r[1, 3]) / sqrt(1 - r[1, 2]^2)   # = -g(rho_hat, 0)
    })
    sqrt(n0) * stats::sd(g0)
  })
  cat(sprintf("  s(0) at delta* = %.3f\n", s0))
  pw$n_lab <- factor(sprintf("n = %s", format(pw$n, big.mark = ",", trim = TRUE)),
                     levels = sprintf("n = %s", format(sort(unique(pw$n)), big.mark = ",", trim = TRUE)))
  pw$h <- sqrt(pw$n) * (pw$delta - DSTAR)
  pcols <- c("#7FA6D9", "#2E5EAA", "#0B2A5B")
  th2 <- theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 10))
  # Each n has its own colour and marker shape, so the three sizes stay
  # distinguishable in panel (b), where they fall on one curve.
  pshapes <- c(16, 17, 15)
  pa <- ggplot(pw, aes(delta, reject, colour = n_lab, shape = n_lab)) +
    geom_hline(yintercept = ALPHA0_P, linetype = "dotted", colour = "gray45") +
    geom_vline(xintercept = DSTAR, linetype = "dashed", colour = "gray45") +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
    scale_colour_manual(values = pcols, name = NULL) +
    scale_shape_manual(values = pshapes, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = "(a) Rejection rate against the violation",
         x = expression(delta == cor(Z, U)), y = "Rejection rate") + th2
  hh <- seq(0, 6, length.out = 300)
  pred <- data.frame(h = hh,
    p = stats::pnorm(hh / sqrt(1 - RHO_P^2) / s0 - stats::qnorm(1 - ALPHA0_P)))
  pb <- ggplot(pw[pw$delta <= 0.20, ], aes(h, reject, colour = n_lab, shape = n_lab)) +
    geom_hline(yintercept = ALPHA0_P, linetype = "dotted", colour = "gray45") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "gray45") +
    geom_line(data = pred, aes(h, p), inherit.aes = FALSE, colour = "gray20",
              linewidth = 0.7) +
    geom_point(size = 1.9, alpha = 0.85) +
    annotate("text", x = 2.0, y = 0.42, hjust = 0, size = 2.7, colour = "gray20",
             lineheight = 0.9, label = "Gaussian prediction\n(one curve for all n)") +
    scale_colour_manual(values = pcols, name = NULL) +
    scale_shape_manual(values = pshapes, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    coord_cartesian(xlim = c(-3, 6)) +
    labs(title = "(b) Against the local alternative",
         x = expression(sqrt(n) ~ (delta - delta[infinity]^"*")), y = NULL) + th2 +
    guides(colour = "none", shape = "none")
  p2 <- (pa | pb) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, "fig2_gate_power.png"), p2,
         width = 7.4, height = 3.4, dpi = 300, bg = "white")
  cat("Saved: fig2_gate_power.png\n")
}

# =============================================================================
# fig:rstar -- implied endogeneity r* and its 90% confidence set against D
#
# One row per instrument that passes the relevance gate: the canonical
# domain D (shaded band), r* (point) and the 90% confidence set CS_n of
# the gate (bar), from results/rstar_gate.csv (12_rstar_gate.R). Colour:
# does the set meet D (the gate passes)? For a rejected instrument, the
# label is the breakdown distance d(r*, D).
# =============================================================================

cat("Generating fig:rstar (r* and its confidence set against the domain)...\n")

cValid   <- "#3D8A3D"
cInvalid <- "#C9433D"
TAU_R <- 1.5
TAU_L <- 0.90

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

.gate_file <- file.path(.repo_root, "results", "rstar_gate.csv")
.ie_file   <- file.path(.repo_root, "results", "implied_endogeneity.csv")
if (!file.exists(.gate_file)) stop("results/rstar_gate.csv not found; run 12 first")
if (!file.exists(.ie_file))   stop("results/implied_endogeneity.csv not found; run 09 first")
gt <- utils::read.csv(.gate_file, stringsAsFactors = FALSE)
ie <- utils::read.csv(.ie_file, stringsAsFactors = FALSE)
stopifnot(!anyNA(gt$V))
gt$weak <- gt$F_first < 10
gt$pass <- gt$V == 1L & !gt$weak

f3 <- gt[!gt$weak, ]
f3$label <- paste(f3$study, f3$instrument)
f3 <- f3[order(f3$pass, f3$margin), ]
f3$row <- factor(f3$label, levels = f3$label)
f3$verdict <- ifelse(f3$pass, "Passes the gate", "Fails the gate")
f3$txt <- ifelse(f3$pass, "", sprintf("d = %.3f", f3$dist))

f3$y <- as.numeric(f3$row)
p3 <- ggplot(f3, aes(y = y)) +
  geom_rect(aes(xmin = D_lo, xmax = D_hi, ymin = y - 0.38, ymax = y + 0.38),
            fill = "gray80", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray45", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = cs_lo, xmax = cs_hi, colour = verdict),
                 height = 0.25, linewidth = 0.6) +
  geom_point(aes(x = r_star, colour = verdict), size = 2.4) +
  geom_text(aes(x = pmax(cs_hi, D_hi) + 0.03, label = txt), hjust = 0,
            size = 2.7, colour = cInvalid) +
  scale_colour_manual(values = c("Passes the gate" = cValid, "Fails the gate" = cInvalid),
                      breaks = c("Passes the gate", "Fails the gate")) +
  scale_x_continuous(limits = c(-1, 1.05), breaks = seq(-1, 1, 0.25)) +
  scale_y_continuous(breaks = f3$y, labels = as.character(f3$row),
                     expand = expansion(add = 0.6)) +
  labs(x = expression("Implied endogeneity " * r^"*" * " with its 90% set (shaded: domain " * D * ")"),
       y = NULL, colour = NULL) +
  theme_dial() +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(FIG_DIR, "fig3_rstar_domain.png"), p3,
       width = 7.4, height = 4.6, dpi = 300, bg = "white")
cat("Saved: fig3_rstar_domain.png\n")

# =============================================================================
# fig:diagnostic_space -- DIAL diagnostic space (13 instruments, three panels)
#
#   Left and middle panels -- the six instruments that pass the relevance
#                  and validity gates, sorted by the share of bootstrap
#                  resamples that land in Scenario B.
#                  Left: R with its 95% bootstrap interval (log scale);
#                  green band: Scenario A, blue regions: Scenario B.
#                  Middle: the share in B, which sets the label of
#                  Algorithm 1 (A when at most 1 - TAU_L of the resamples
#                  land in B, B when at least TAU_L do, inconclusive
#                  otherwise). Open circle: IV and DML of opposite sign
#                  (num_genres), with the note "signs differ" under its
#                  share in the middle panel.
#   Right panel -- the seven Scenario C instruments (rejected or weak), by
#                  the signed distance of r* to the closed domain D (grey;
#                  for all seven, minus the breakdown distance). For the
#                  five that the validity gate rejects, the bar is the 90%
#                  set CS(r*) in the same units: it ends short of D, the
#                  "Fail" branch of "CS(r*) meets D?" in fig:architecture.
#                  Open circles: the two weak instruments (no gate, no set).
#   The panels run A and B, then C, as the outcomes do in fig:architecture.
# =============================================================================

cat("Generating fig:diagnostic_space (DIAL diagnostic space)...\n")

cA <- "#3D8A3D"; cB <- "#2E5EAA"; cI <- "#8A6D3B"; cC <- "#C9433D"

dd <- merge(gt[, c("study", "instrument", "margin", "dist", "pass", "weak", "V",
                   "cs_lo", "cs_hi", "D_lo", "D_hi")],
            ie[, c("instrument", "b_IV", "b_DML", "R", "R_lo", "R_hi", "share_B")],
            by = "instrument", sort = FALSE)
dd$label3 <- ifelse(!dd$pass, "C",
               ifelse(dd$share_B >= TAU_L, "B",
                 ifelse(dd$share_B <= 1 - TAU_L, "A", "inconclusive")))
dd$sign_diff <- !is.na(dd$b_DML) & sign(dd$b_IV) != sign(dd$b_DML)

valid4 <- dd[dd$label3 != "C", ]
invalid4 <- dd[dd$label3 == "C", ]
stopifnot(nrow(valid4) + nrow(invalid4) == 13L)
invalid4$mark <- ifelse(invalid4$weak, "weak", "out")
# study on the first line, instrument on the second, to leave room for the panels
invalid4$name <- paste0(invalid4$study, "\n", invalid4$instrument,
                        ifelse(invalid4$weak, " (weak)", ""))
invalid4 <- invalid4[order(invalid4$margin), ]
invalid4$y <- seq_len(nrow(invalid4))
# Signed distance to the closed domain [lo, hi]: positive inside, negative
# outside. For the five instruments that the gate rejects, the 90% set
# CS(r*) is drawn in the same units; it must end short of D. The weak
# instruments stop at the relevance gate, so no set is drawn for them.
sdist <- function(x, lo, hi) ifelse(x < lo, x - lo, ifelse(x > hi, hi - x, pmin(x - lo, hi - x)))
s_a <- ifelse(invalid4$weak, NA, sdist(invalid4$cs_lo, invalid4$D_lo, invalid4$D_hi))
s_b <- ifelse(invalid4$weak, NA, sdist(invalid4$cs_hi, invalid4$D_lo, invalid4$D_hi))
invalid4$s_min <- pmin(s_a, s_b); invalid4$s_max <- pmax(s_a, s_b)
stopifnot(all(invalid4$s_max[!invalid4$weak] < 0))   # every set misses D
# The left panels have nrow(valid4) rows and a header band of 0.8 row units
# (labels 0.45 below the top). Scale the band here so it has the same height
# on the page, which puts the region labels of all three panels on one line.
nC <- nrow(invalid4); nV0 <- nrow(valid4)
hC <- 0.8 * (nC + 0.1) / (nV0 + 0.1)
yMaxC <- nC + 0.5 + hC; yHdC <- yMaxC - 0.45 * hC / 0.8
# x range; D is [0, XC_HI]. XC_LO leaves room left of -1 for the whole circle
# of is_sequel (-0.98), also when wider fonts (Windows) narrow the panel.
XC_LO <- -1.1; XC_HI <- 0.2
gridC <- c(-1, -0.75, -0.5, -0.25)

pC <- ggplot(invalid4, aes(y = y)) +
  annotate("rect", xmin = XC_LO, xmax = 0, ymin = 0.4, ymax = yMaxC,
           fill = cC, alpha = 0.06) +
  annotate("rect", xmin = 0, xmax = XC_HI, ymin = 0.4, ymax = yMaxC,
           fill = "gray55", alpha = 0.22) +
  annotate("segment", x = gridC, xend = gridC, y = 0.4, yend = nC + 0.5,
           colour = "gray85", linewidth = 0.3) +
  annotate("segment", x = 0, xend = 0, y = 0.4, yend = nC + 0.5,
           colour = "gray45", linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = c(XC_LO / 2, XC_HI / 2), y = yHdC, label = c("outside D", "D"),
           colour = c(cC, "gray30"), fontface = "bold", size = c(2.8, 3.3)) +
  geom_errorbarh(aes(xmin = s_min, xmax = s_max), colour = cC, height = 0.25,
                 linewidth = 0.55, na.rm = TRUE) +
  geom_point(aes(x = margin, shape = mark), colour = cC, fill = "white", size = 2.4,
             stroke = 1) +
  scale_shape_manual(values = c(out = 16, weak = 21), guide = "none") +
  scale_x_continuous(limits = c(XC_LO, XC_HI), breaks = c(-1, -0.75, -0.5, -0.25, 0),
                     labels = c("-1", "", "-0.5", "", "0"), expand = c(0, 0)) +
  scale_y_continuous(breaks = invalid4$y, labels = invalid4$name,
                     limits = c(0.4, yMaxC), expand = c(0, 0)) +
  # parallel to "A and B: validity not rejected"; the two weak instruments
  # stop at the relevance gate, before the validity gate. The subtitle uses
  # the wording of the gate in fig:architecture, "CS(r*) meets D?"; italic()
  # because plotmath ignores the italic face set in the theme.
  labs(title = "C: validity rejected (or weak instrument)",
       subtitle = expression(italic(CS(r^"*") ~ "misses" ~ D * "," ~ "or" ~ F < 10)),
       x = expression("Signed distance to " * D), y = NULL) +
  theme_dial() +
  theme(plot.title = element_text(face = "bold", colour = cC, size = 10, hjust = 0.5),
        plot.subtitle = element_text(colour = "gray35", size = 7.5, face = "italic",
                                     hjust = 0.5),
        # centred over names and panel: the title is wider than this panel,
        # which now sits at the right edge of the figure
        plot.title.position = "plot",
        axis.text.y = element_text(size = 8, colour = "gray15"),
        panel.grid = element_blank(),
        plot.margin = margin(4, 10, 4, 4))

# Left and middle panels: one row per instrument that passes, sorted by the
# share of the bootstrap resamples that land in B, the quantity that sets the
# label. Left: R with its 95% bootstrap interval; the green band is
# Scenario A (1/TAU_R < R < TAU_R) and the blue regions on either side are
# Scenario B. A resample lands in B when its R falls in a blue region or its
# IV and DML estimates differ in sign. Middle: that share against the label
# rule of Algorithm 1 (A when at most 1 - TAU_L of the resamples land in B,
# B when at least TAU_L do, inconclusive otherwise). The label thus depends
# on the whole bootstrap distribution of R, not on where the point estimate
# lies, so the left panel draws the points in grey. Open circle: IV and DML
# of opposite sign (num_genres); R uses absolute values, so the middle panel
# adds the note "signs differ" under its share. Grid lines stop below the
# header row, so they do not cross its labels.
valid4 <- valid4[order(valid4$share_B, valid4$R), ]   # largest share at the top
valid4$y <- seq_len(nrow(valid4))
valid4$name <- paste0(valid4$study, "\n", valid4$instrument)
valid4$share_pct <- 100 * valid4$share_B
# shares are multiples of 0.5% (200 resamples): show them exactly, e.g. 2.5% and 99.5%
valid4$share_txt <- paste0(formatC(valid4$share_pct, format = "f", digits = 1,
                                   drop0trailing = TRUE), "%")
valid4$lab_f <- factor(valid4$label3, levels = c("A", "B", "inconclusive"))
nV    <- nrow(valid4)
yHd   <- nV + 0.85           # header row with the region labels
yMax  <- nV + 1.3
X_LO  <- 0.04; X_HI <- 30    # x range of R (log scale)
gridR <- c(0.1, 0.5, 1, 2, 5, 20)
cGrid <- "gray85"
th_side <- theme(plot.title = element_text(face = "bold", colour = "gray20", size = 10,
                                           hjust = 0.5),
                 plot.subtitle = element_text(colour = "gray35", size = 7.5, face = "italic",
                                              hjust = 0.5),
                 plot.title.position = "panel",
                 panel.grid = element_blank())

pAB <- ggplot(valid4, aes(y = y)) +
  annotate("rect", xmin = X_LO, xmax = 1 / TAU_R, ymin = 0.4, ymax = yMax,
           fill = cB, alpha = 0.08) +
  annotate("rect", xmin = 1 / TAU_R, xmax = TAU_R, ymin = 0.4, ymax = yMax,
           fill = cA, alpha = 0.16) +
  annotate("rect", xmin = TAU_R, xmax = X_HI, ymin = 0.4, ymax = yMax,
           fill = cB, alpha = 0.08) +
  annotate("segment", x = gridR, xend = gridR, y = 0.4, yend = nV + 0.5,
           colour = cGrid, linewidth = 0.3) +
  annotate("text", x = c(sqrt(X_LO / TAU_R), 1, sqrt(TAU_R * X_HI)), y = yHd,
           label = c("B", "A", "B"), colour = c(cB, cA, cB), fontface = "bold", size = 3.3) +
  geom_errorbarh(aes(xmin = R_lo, xmax = R_hi), height = 0.25, linewidth = 0.55,
                 colour = "gray35") +
  geom_point(aes(x = R, shape = sign_diff), size = 2.4, stroke = 1, colour = "gray15",
             fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
  scale_x_log10(limits = c(X_LO, X_HI), breaks = gridR,
                labels = c("0.1", "0.5", "1", "2", "5", "20"), expand = c(0, 0)) +
  scale_y_continuous(breaks = valid4$y, labels = valid4$name,
                     limits = c(0.4, yMax), expand = c(0, 0)) +
  labs(title = "A and B: validity not rejected",
       subtitle = "R with its 95% bootstrap interval (log scale)",
       x = expression("R = |" * hat(beta)[IV] / hat(beta)[DML] * "|"), y = NULL) +
  theme_dial() + th_side +
  theme(axis.text.y = element_text(size = 8, colour = "gray15"),
        plot.margin = margin(4, 4, 4, 8))

S_LO <- 100 * (1 - TAU_L); S_HI <- 100 * TAU_L   # 10% and 90%
valid4$s_tx <- ifelse(valid4$share_pct > 50, valid4$share_pct - 5, valid4$share_pct + 5)
valid4$s_hj <- ifelse(valid4$share_pct > 50, 1, 0)
pS <- ggplot(valid4, aes(y = y)) +
  annotate("rect", xmin = 0, xmax = S_LO, ymin = 0.4, ymax = yMax,
           fill = cA, alpha = 0.16) +
  annotate("rect", xmin = S_LO, xmax = S_HI, ymin = 0.4, ymax = yMax,
           fill = cI, alpha = 0.12) +
  annotate("rect", xmin = S_HI, xmax = 100, ymin = 0.4, ymax = yMax,
           fill = cB, alpha = 0.10) +
  annotate("segment", x = 50, xend = 50, y = 0.4, yend = nV + 0.5,
           colour = cGrid, linewidth = 0.3) +
  annotate("text", x = c(S_LO / 2, 50, (S_HI + 100) / 2), y = yHd,
           label = c("A", "inconclusive", "B"), colour = c(cA, cI, cB),
           fontface = c("bold", "bold.italic", "bold"), size = c(3.3, 2.5, 3.3)) +
  geom_point(aes(x = share_pct, colour = lab_f, shape = sign_diff), size = 2.4,
             stroke = 1, fill = "white") +
  geom_text(aes(x = s_tx, label = share_txt, hjust = s_hj), size = 2.5,
            colour = "gray20") +
  geom_text(data = valid4[valid4$sign_diff, ],
            aes(x = s_tx, y = y - 0.3, hjust = s_hj), label = "signs differ",
            size = 2.3, fontface = "italic", colour = "gray35") +
  scale_colour_manual(values = c(A = cA, B = cB, inconclusive = cI), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = c(S_LO, 50, S_HI),
                     labels = paste0(c(S_LO, 50, S_HI), "%"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, yMax), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Agreement label", subtitle = "share of resamples in B",
       x = "% of 200 resamples", y = NULL) +
  theme_dial() + th_side +
  theme(axis.text.y = element_blank(), plot.margin = margin(4, 18, 4, 6))

p4 <- (pAB | pS | pC) +
  plot_layout(widths = c(1.45, 0.95, 1)) +
  plot_annotation(
    title = "DIAL diagnostic space: 13 instruments across 8 studies",
    theme = theme(plot.title = element_text(face = "bold", size = 12,
                                            margin = margin(0, 0, 6, 0))))

ggsave(file.path(FIG_DIR, "fig4_diagnostic_space.png"), p4,
       width = 7.4, height = 4.0, dpi = 300, bg = "white")
cat("Saved: fig4_diagnostic_space.png\n")

cat("\nAll figures saved to", FIG_DIR, "\n")
