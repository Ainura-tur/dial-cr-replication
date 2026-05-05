# =============================================================================
# 08_figures.R
#
# Produces all paper figures:
#   Figure 1: DIAL architecture diagram      -> paper/figures/fig1_dial_architecture.png
#   Figure 2: Power curve (3 panels)         -> paper/figures/fig2_combined.png
#   Figure 3: Domain sensitivity             -> paper/figures/fig3_domain_sensitivity.png
#   Figure 4: DIAL diagnostic space          -> paper/figures/fig4_diagnostic_space.png
#
# Figure 2 requires results/cp_grid_summary_lean.csv produced by
# scripts/07_power_curve_grid.R. All other figures are self-contained.
#
# Usage:
#   Rscript scripts/08_figures.R
#
# Runtime: under 2 minutes (excluding figure 2 if CSV is missing).
# =============================================================================

# ---- Repo-root detection ----------------------------------------------------
.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) && nchar(f)) dirname(normalizePath(f))
    else getwd()
  }
)
.repo_root <- normalizePath(file.path(.script_dir, ".."))
FIG_DIR    <- file.path(.repo_root, "paper", "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

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

# =============================================================================
# FIGURE 1: DIAL Architecture Diagram
# =============================================================================

cat("Generating Figure 1: DIAL architecture diagram...\n")

png(file.path(FIG_DIR, "fig1_dial_architecture.png"),
    width = 3400, height = 2800, res = 300)

par(mar = c(0.2, 0.2, 0.6, 0.2), family = "sans")
plot.new()
plot.window(xlim = c(0, 100), ylim = c(0, 100))

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

cA  <- "#7FC97F"; cAb <- "#3D8A3D"
cB  <- "#7FAEDC"; cBb <- "#2A6FB7"
cC  <- "#F4A09A"; cCb <- "#C9433D"
cBg <- "#EFEFE6"; cBgb <- "#8C8C82"
cE  <- "#D6E7F8"; cEb <- "#3F77B5"; cEt <- "#1B4C82"
cR  <- "#E1DDF5"; cRb <- "#6E64C2"; cRt <- "#3C3489"
cT  <- "#D8EFE5"; cTb <- "#3FA083"; cTt <- "#0F6E56"
cF  <- "#FCE3BC"; cFb <- "#E08826"; cFt <- "#9C5A0E"

xData <- 50; yData <- 95
xDML <- 22; xIV <- 70; yEst <- 84
xF <- 70; yF <- 71
xCR <- 50; yCR <- 57
xRat <- 22; yRat <- 46
xR <- 22; yR <- 35
xA <- 9; xB <- 35; xC <- 78; yOut <- 22

text(50, 99, "DIAL diagnostic framework", cex = 1.30, font = 2, col = "gray10")

rrect(xData, yData, 30, 4, cBg, cBgb, lwd = 1.4)
text(xData, yData, "Observed data  (Y, D, Z, X)", col = "gray15", cex = 0.95, font = 2)

rrect(xDML, yEst, 16, 5.5, cE, cEb, lwd = 1.4)
text(xDML, yEst + 1.1, "DML",   col = cEt, cex = 0.95, font = 2)
text(xDML, yEst - 1.3, "(ATE)", col = cEt, cex = 0.65)

rrect(xIV, yEst, 16, 5.5, cE, cEb, lwd = 1.4)
text(xIV, yEst + 1.1, "IV",     col = cEt, cex = 0.95, font = 2)
text(xIV, yEst - 1.3, "(LATE)", col = cEt, cex = 0.65)

arr_L_HV(xData - 15, yData - 1.7, xDML, yEst + 3, col = "gray35")
arr_L_HV(xData + 15, yData - 1.7, xIV,  yEst + 3, col = "gray35")

arr(xIV, yEst - 3, xIV, yF + 4)
diam(xF, yF, 16, 7, cF, cFb)
text(xF, yF + 0.7, "Weak IV test", col = cFt, cex = 0.70, font = 2)
text(xF, yF - 1.3, "F < 10 ?",     col = cFt, cex = 0.72, font = 2)

arr_L_HV(xF - 8, yF, xCR, yCR + 4.3, col = "gray35")
text((xF - 8 + xCR)/2, yF + 1.6, "No", col = "gray30", cex = 0.65, font = 2)

arr_L_HV(xF + 8, yF, xC + 5, yOut + 4.5, col = cCb, lwd = 1.4)
text(xF + 11, yF + 1.6, "Yes", col = cCb, cex = 0.65, font = 2)

diam(xCR, yCR, 22, 8.5, cR, cRb)
text(xCR, yCR + 1.0, "CR test (MCUB)", col = cRt, cex = 0.78, font = 2)
text(xCR, yCR - 1.4, expression(0 %in% hat(C) * "?"), col = cRt, cex = 0.78)

xPassDrop <- xRat + 7
segments(xCR - 11, yCR, xPassDrop, yCR, col = cAb, lwd = 1.4)
arrows(xPassDrop, yCR, xPassDrop, yRat + 3, length = 0.08, angle = 22, col = cAb, lwd = 1.4)
text((xCR - 11 + xPassDrop)/2, yCR + 1.3, "Pass", col = cAb, cex = 0.65, font = 2)

arr_L_HV(xCR + 11, yCR, xC - 6, yOut + 4.5, col = cCb, lwd = 1.4)
text((xCR + 11 + xC - 6)/2, yCR + 1.6, "Fail", col = cCb, cex = 0.65, font = 2)

arr(xDML, yEst - 3, xDML, yRat + 3, col = "gray35")

rrect(xRat, yRat, 22, 6, cT, cTb, lwd = 1.4)
text(xRat, yRat + 1.1, "Effect ratio (R)", col = cTt, cex = 0.78, font = 2)
text(xRat, yRat - 1.4,
     expression(group("|", hat(beta)[IV]/hat(beta)[DML], "|")),
     col = cTt, cex = 0.70)

arr(xRat, yRat - 3.1, xR, yR + 4)

diam(xR, yR, 17, 9, cF, cFb)
text(xR, yR + 1.2, expression("R" %~~% "1 ?"), col = cFt, cex = 0.80, font = 2)

arr_L_HV(xR - 8.6, yR, xA, yOut + 4.5, col = cAb, lwd = 1.4)
text(xR - 11, yR + 1.5, "Yes", col = cAb, cex = 0.65, font = 2)

arr_L_HV(xR + 8.6, yR, xB, yOut + 4.5, col = cBb, lwd = 1.4)
text(xR + 11, yR + 1.5, "No", col = cBb, cex = 0.65, font = 2)

rrect(xA, yOut, 15, 9, cA, cAb, lwd = 1.5)
text(xA, yOut + 2.4, "A",               col = "white",   cex = 1.30, font = 2)
text(xA, yOut - 0.3, "Global valid IV",  col = "#0F4F0F", cex = 0.62, font = 2)
text(xA, yOut - 2.2, expression("ATE" %~~% "LATE"), col = "#0F4F0F", cex = 0.65)

rrect(xB, yOut, 14, 9, cB, cBb, lwd = 1.5)
text(xB, yOut + 2.4, "B",     col = "white",   cex = 1.30, font = 2)
text(xB, yOut - 0.3, "Local", col = "#0E2E55", cex = 0.65, font = 2)
text(xB, yOut - 2.2, expression("LATE" != "ATE"), col = "#0E2E55", cex = 0.65)

rrect(xC, yOut, 16, 9, cC, cCb, lwd = 1.5)
text(xC, yOut + 2.4, "C",               col = "white",   cex = 1.30, font = 2)
text(xC, yOut - 0.3, "IV unreliable",   col = "#5C1612", cex = 0.65, font = 2)
text(xC, yOut - 2.2, "Weak or invalid", col = "#5C1612", cex = 0.62)

dev.off()
cat("Saved: fig1_dial_architecture.png\n")

# =============================================================================
# FIGURE 2: Power curve (3 panels)
# Requires results/cp_grid_summary_lean.csv from 07_power_curve_grid.R
# =============================================================================

cat("Generating Figure 2: power curves...\n")

csv_path <- file.path(.repo_root, "results", "cp_grid_summary_lean.csv")

if (!file.exists(csv_path)) {
  warning("cp_grid_summary_lean.csv not found. Run 07_power_curve_grid.R first. Skipping Figure 2.")
} else {
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
  
  # Interpolated 50%-crossing per n
  delta_star <- dat %>%
    arrange(n, delta) %>%
    group_by(n, n_label) %>%
    group_modify(~{
      d <- .x$delta; cp <- .x$cp_raw
      idx <- which(cp[-length(cp)] > 0.5 & cp[-1] <= 0.5)
      if (!length(idx)) return(data.frame(delta_star = NA_real_))
      i <- idx[1]
      t <- (0.5 - cp[i]) / (cp[i+1] - cp[i])
      data.frame(delta_star = d[i] + t * (d[i+1] - d[i]))
    }) %>%
    ungroup()
  
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
  
  # Panel a: full power curves, linear x-axis
  p_a <- ggplot(dat, aes(x = delta, y = cp_raw,
                         colour = n_label, group = n_label)) +
    geom_hline(yintercept = 0.95, linetype = "dotted",
               colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dotted",
               colour = "grey30", linewidth = 0.3) +
    geom_vline(data = delta_star,
               aes(xintercept = delta_star, colour = n_label),
               linetype = "dashed", linewidth = 0.4, alpha = 0.7,
               show.legend = FALSE) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.4) +
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
    theme_fig2()
  
  # Panel b: transition zone zoom
  zoom_lo <- 0.125; zoom_hi <- 0.200
  dat_zoom <- dat %>% filter(delta >= zoom_lo & delta <= zoom_hi)
  
  p_b <- ggplot(dat_zoom, aes(x = delta, y = cp_raw,
                              colour = n_label, group = n_label)) +
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
    geom_point(size = 1.6) +
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
      title = expression("(b) Transition zone: " * delta * " near " * delta["*,n"]),
      x     = expression(delta * " (transition zone)"),
      y     = expression(CP[n](delta))
    ) +
    theme_fig2()
  
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
  
  combined <- (p_a | p_b | p_c) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(file.path(FIG_DIR, "fig2_combined.png"), combined,
         width = 16, height = 4.8, units = "in", dpi = 300)
  cat("Saved: fig2_combined.png\n")
}

# =============================================================================
# FIGURE 3: Domain Sensitivity
# =============================================================================

cat("Generating Figure 3: domain sensitivity...\n")

cValid   <- "#3D8A3D"
cInvalid <- "#C9433D"

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
  id = seq_len(12),
  row_label = c(
    "Di Tella judgeAlreadyUsedEM\nWide       (0, 0.8)",
    "Di Tella judgeAlreadyUsedEM\nCanonical  (0.15, 0.25)",
    "Di Tella percJudgeSent\nWide       (0, 0.8)",
    "Di Tella percJudgeSent\nCanonical  (0, 0.25)",
    "Galiani highnumber\nWide       (-0.8, 0)",
    "Galiani highnumber\nCanonical  (-0.2, 0)",
    "Banerjee instru\nWide       (-0.8, 0)",
    "Banerjee instru\nCanonical  (-0.6, -0.4)",
    "Burde buildschool\nWide       (0, 0.8)",
    "Burde buildschool\nCanonical  (0, 0.4)",
    "Card fatheduc\nCanonical  (0, 0.8)",
    "Card fatheduc\nNarrow     (0, 0.4)"
  ),
  lower = c(-0.256, -0.109, -0.046, -0.046, -0.793, -0.189,
            -0.375, -0.101,  0.253,  0.238,  0.063,  0.063),
  upper = c( 0.638,  0.019,  0.775,  0.210,  0.012,  0.012,
             0.544,  0.163,  0.943,  0.664,  0.841,  0.465),
  verdict = c("Valid", "Valid", "Valid", "Valid", "Valid", "Valid",
              "Valid", "Valid", "Invalid", "Invalid", "Invalid", "Invalid"),
  stringsAsFactors = FALSE
)

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
    x        = expression("Identified set for " * rho[ZU]),
    y        = NULL,
    title    = "Domain Sensitivity: MCUB Identified Sets (Within-Sign Protocol)",
    subtitle = "Each instrument shown under 2 within-sign domains; every verdict is robust under v2",
    color    = "Verdict"
  ) +
  theme_dial() +
  guides(color = guide_legend(override.aes = list(linewidth = 5, size = 0)))

ggsave(file.path(FIG_DIR, "fig3_domain_sensitivity.png"), p3,
       width = 8.2, height = 6.8, dpi = 300)
cat("Saved: fig3_domain_sensitivity.png\n")

# =============================================================================
# FIGURE 4: DIAL Diagnostic Space (13 instruments, 2-D)
# =============================================================================

cat("Generating Figure 4: DIAL diagnostic space...\n")

cA <- "#3D8A3D"; cB <- "#2E5EAA"; cC <- "#C9433D"

fig4 <- data.frame(
  label = c(
    "MovieLens\nnum_genres", "Di Tella\njudgeAlreadyUsedEM",
    "Di Tella\npercJudgeSent", "Galiani\nhighnumber",
    "Banerjee\ninstru", "401(k)\nsole plan", "401(k)\nplan age",
    "NHANES\nincome", "NHANES\neducation", "Card\nfatheduc",
    "Burde\nbuildschool", "MovieLens\nrelease_year", "MovieLens\nis_sequel"
  ),
  valid = c(rep(TRUE, 5), rep(FALSE, 8)),
  ratio = c(0.61, 3.47, 1.33, 2.41, 2.71,
            3.40, 1.30, 2.20, 2.40, 1.60, 1.59, 1.30, 2.20),
  dial  = c("A", "B", "A", "B", "B", rep("C", 8)),
  stringsAsFactors = FALSE
)

set.seed(42)
n_valid   <- sum( fig4$valid)
n_invalid <- sum(!fig4$valid)
fig4$x <- NA_real_
fig4$x[ fig4$valid] <- 1 + runif(n_valid,   -0.08, 0.08)
fig4$x[!fig4$valid] <- seq(-0.22, 0.22, length.out = n_invalid)
fig4$nudge_x <- 0
fig4$nudge_x[fig4$label == "Burde\nbuildschool"]      <- -0.06
fig4$nudge_x[fig4$label == "MovieLens\nrelease_year"] <- -0.04

p4 <- ggplot(fig4, aes(x = x, y = ratio)) +
  # Use geom_rect with a dummy data frame so coord_cartesian does not clip them:
  geom_rect(data = data.frame(
    xmin = c(0.5, 0.5, -0.35),
    xmax = c(1.30, 1.30, 0.5),
    ymin = c(0.0, 1.5, 0.0),
    ymax = c(1.5, 4.5, 4.5),
    fill = c(cA, cB, cC),
    alpha = c(0.08, 0.06, 0.05)
  ),
  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      fill = I(fill), alpha = I(alpha)),
  inherit.aes = FALSE)  +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
  annotate("segment", x = 0.5, xend = 1.30, y = 1.5, yend = 1.5,
           linetype = "dashed", color = "gray50", linewidth = 0.3) +
  annotate("text", x = 0.8, y = 0.30, label = "A: Valid\nhomogeneous",
           size = 3.0, color = cA, fontface = "bold", lineheight = 0.85) +
  annotate("text", x = 0.8, y = 4.25, label = "B: Valid\nheterogeneous",
           size = 3.0, color = cB, fontface = "bold", lineheight = 0.85) +
  annotate("text", x = -0.10, y = 4.25, label = "C: Invalid",
           size = 3.0, color = cC, fontface = "bold") +
  geom_point(aes(color = dial), size = 3.2, alpha = 0.9) +
  geom_text_repel(
    aes(label = label, color = dial),
    size = 2.3, fontface = "bold", lineheight = 0.85,
    box.padding = 0.8, point.padding = 0.6,
    force = 4, force_pull = 0.2,
    min.segment.length = 0.02, segment.ncp = -1,
    max.iter = 10000, direction = "y",
    segment.color = "gray60", segment.size = 0.3,
    max.overlaps = Inf, seed = 17
  ) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = c("A" = cA, "B" = cB, "C" = cC), guide = "none") +
  scale_x_continuous(
    limits = c(-0.35, 1.30),
    breaks = c(0, 0.5, 1),
    labels = c("Invalid\n(0 \u2209 CI)", "", "Valid\n(0 \u2208 CI)")
  ) +
  scale_y_continuous(limits = c(0, 4.5), breaks = 0:4) +
  labs(
    x        = "CR-test verdict (MCUB)",
    y        = expression("|" * hat(beta)[IV] / hat(beta)[DML] * "|"),
    title    = "DIAL Diagnostic Space: 13 Instruments Across 8 Studies",
    subtitle = "5 of 13 pass CR test under canonical D; all 13 verdicts robust to within-sign D widening"
  ) +
  theme_dial()

ggsave(file.path(FIG_DIR, "fig4_diagnostic_space.png"), p4,
       width = 7.5, height = 5.3, dpi = 300)
cat("Saved: fig4_diagnostic_space.png\n")

cat("\nAll figures saved to", FIG_DIR, "\n")
