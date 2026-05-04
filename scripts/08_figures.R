# =============================================================================
# DIAL Architecture Diagram v20
# Changes from v19:
#   - CR "Pass" route: goes DOWN from CR left vertex, then LEFT into ratio box.
#     This eliminates the crossing with the DML→Ratio vertical at xDML.
#   - R-gate diamond shows both R ~ 1 ? and delta = 0.5 inside.
#   - Inline expansion (0.5 < R < 1.5) shown to the right of the R-gate.
#   - Outcome labels match reference: A "Global valid IV / ATE ~= LATE",
#                                     B "Local / LATE != ATE",
#                                     C "IV unreliable / Weak or invalid".
# ==============================================================================

png("fig0_dial_architecture.png", width = 3400, height = 2800, res = 300)

par(mar = c(0.2, 0.2, 0.6, 0.2), family = "sans")
plot.new()
plot.window(xlim = c(0, 100), ylim = c(0, 100))

# --- Helpers ---
rrect <- function(x, y, w, h, fill, bdr, lwd = 1.5) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = fill, border = bdr, lwd = lwd)
}

diam <- function(cx, cy, w, h, fill, bdr, lwd = 1.8) {
  polygon(c(cx, cx + w/2, cx, cx - w/2),
          c(cy + h/2, cy, cy - h/2, cy),
          col = fill, border = bdr, lwd = lwd)
}

arr <- function(x1, y1, x2, y2, col = "gray35", lwd = 1.4) {
  arrows(x1, y1, x2, y2, length = 0.08, angle = 22, col = col, lwd = lwd)
}

# L-shape: horizontal first (x1->x2 at y1), then vertical (y1->y2 at x2)
arr_L_HV <- function(x1, y1, x2, y2, col = "gray35", lwd = 1.4) {
  segments(x1, y1, x2, y1, col = col, lwd = lwd)
  arrows(x2, y1, x2, y2, length = 0.08, angle = 22, col = col, lwd = lwd)
}

# L-shape: vertical first (y1->y2 at x1), then horizontal (x1->x2 at y2)
arr_L_VH <- function(x1, y1, x2, y2, col = "gray35", lwd = 1.4) {
  segments(x1, y1, x1, y2, col = col, lwd = lwd)
  arrows(x1, y2, x2, y2, length = 0.08, angle = 22, col = col, lwd = lwd)
}

# --- Color palette ---
cA  <- "#7FC97F"; cAb <- "#3D8A3D"
cB  <- "#7FAEDC"; cBb <- "#2A6FB7"
cC  <- "#F4A09A"; cCb <- "#C9433D"
cBg <- "#EFEFE6"; cBgb <- "#8C8C82"
cE  <- "#D6E7F8"; cEb <- "#3F77B5"; cEt <- "#1B4C82"
cR  <- "#E1DDF5"; cRb <- "#6E64C2"; cRt <- "#3C3489"
cT  <- "#D8EFE5"; cTb <- "#3FA083"; cTt <- "#0F6E56"
cF  <- "#FCE3BC"; cFb <- "#E08826"; cFt <- "#9C5A0E"

# --- Layout ---
xData <- 50; yData <- 95
xDML <- 22; xIV <- 70; yEst <- 84
xF <- 70; yF <- 71
xCR <- 50; yCR <- 57
# Ratio/R-gate/Outcomes block pulled UP (old: yRat=43, yR=30, yOut=14)
xRat <- 22; yRat <- 46
xR <- 22; yR <- 35
xA <- 9; xB <- 35; xC <- 78; yOut <- 22

# --- Title ---
text(50, 99, "DIAL diagnostic framework",
     cex = 1.30, font = 2, col = "gray10")

# ==========================================
# Observed data
# ==========================================
rrect(xData, yData, 30, 4, cBg, cBgb, lwd = 1.4)
text(xData, yData, "Observed data  (Y, D, Z, X)",
     col = "gray15", cex = 0.95, font = 2)

# ==========================================
# Estimators
# ==========================================
rrect(xDML, yEst, 16, 5.5, cE, cEb, lwd = 1.4)
text(xDML, yEst + 1.1, "DML",   col = cEt, cex = 0.95, font = 2)
text(xDML, yEst - 1.3, "(ATE)", col = cEt, cex = 0.65)

rrect(xIV, yEst, 16, 5.5, cE, cEb, lwd = 1.4)
text(xIV, yEst + 1.1, "IV",     col = cEt, cex = 0.95, font = 2)
text(xIV, yEst - 1.3, "(LATE)", col = cEt, cex = 0.65)

# Data -> DML, Data -> IV (symmetric L-shapes)
arr_L_HV(xData - 15, yData - 1.7, xDML, yEst + 3, col = "gray35")
arr_L_HV(xData + 15, yData - 1.7, xIV,  yEst + 3, col = "gray35")

# ==========================================
# F-gate (under IV)
# ==========================================
arr(xIV, yEst - 3, xIV, yF + 4)

diam(xF, yF, 16, 7, cF, cFb)
text(xF, yF + 0.7, "Weak IV test", col = cFt, cex = 0.70, font = 2)
text(xF, yF - 1.3, "F < 10 ?",     col = cFt, cex = 0.72, font = 2)

# F-gate "No" -> goes left then down to CR test top
arr_L_HV(xF - 8, yF, xCR, yCR + 4.3, col = "gray35")
text((xF - 8 + xCR)/2, yF + 1.6, "No", col = "gray30",
     cex = 0.65, font = 2)

# F-gate "Yes" -> goes right then down to C
arr_L_HV(xF + 8, yF, xC + 5, yOut + 4.5, col = cCb, lwd = 1.4)
text(xF + 11, yF + 1.6, "Yes", col = cCb, cex = 0.65, font = 2)

# ==========================================
# CR test (center)
# ==========================================
diam(xCR, yCR, 22, 8.5, cR, cRb)
text(xCR, yCR + 1.0, "CR test (MCUB)",
     col = cRt, cex = 0.78, font = 2)
text(xCR, yCR - 1.4,
     expression(0 %in% hat(C) * "?"),
     col = cRt, cex = 0.78)

# ---- CR "Pass" route (KEY FIX: down, left, down — enters ratio box from top)
# From CR left vertex (xCR - 11, yCR):
#   1) go DOWN to a routing level above the ratio box top
#   2) go LEFT horizontally to above the ratio box's right-of-centre
#   3) go DOWN into the ratio box top edge
yPassRoute <- yRat + 6     # well above ratio box top (which is at yRat + 3)
xPassDrop  <- xRat + 7     # enter top of ratio box, offset right of its center
# Horizontal first (left)
segments(xCR - 11, yCR, xPassDrop, yCR, col = cAb, lwd = 1.4)

# Then vertical down into the box
arrows(xPassDrop, yCR, xPassDrop, yRat + 3,
       length = 0.08, angle = 22, col = cAb, lwd = 1.4)

# Label stays centered on horizontal segment
text((xCR - 11 + xPassDrop)/2, yCR + 1.3, "Pass",
     col = cAb, cex = 0.65, font = 2)


# CR "Fail" -> right then down to C
arr_L_HV(xCR + 11, yCR, xC - 6, yOut + 4.5, col = cCb, lwd = 1.4)
text((xCR + 11 + xC - 6)/2, yCR + 1.6, "Fail",
     col = cCb, cex = 0.65, font = 2)

# ==========================================
# DML -> Effect ratio (straight vertical down the left column)
# ==========================================
arr(xDML, yEst - 3, xDML, yRat + 3, col = "gray35")

# ==========================================
# Effect ratio
# ==========================================
rrect(xRat, yRat, 22, 6, cT, cTb, lwd = 1.4)
text(xRat, yRat + 1.1, "Effect ratio (R)",
     col = cTt, cex = 0.78, font = 2)
text(xRat, yRat - 1.4,
     expression(group("|", hat(beta)[IV]/hat(beta)[DML], "|")),
     col = cTt, cex = 0.70)

# Effect ratio -> R-gate
arr(xRat, yRat - 3.1, xR, yR + 4)

# ==========================================
# R-gate — UPDATED label
# ==========================================
diam(xR, yR, 17, 9, cF, cFb)
text(xR, yR + 1.2,
     expression("R" %~~% "1 ?"),
     col = cFt, cex = 0.80, font = 2)
# text(xR, yR - 1.2,                          # was: expression(delta == 0.5)
#      expression("0.5 < R < 1.5"),
#      col = cFt, cex = 0.62)                 # slightly smaller to fit width

# R-gate "Yes" -> A
arr_L_HV(xR - 8.6, yR, xA, yOut + 4.5, col = cAb, lwd = 1.4)
text(xR - 11, yR + 1.5, "Yes", col = cAb, cex = 0.65, font = 2)

# R-gate "No" -> B  (straight down to B)
arr_L_HV(xR + 8.6, yR, xB, yOut + 4.5, col = cBb, lwd = 1.4)
text(xR + 11, yR + 1.5, "No", col = cBb, cex = 0.65, font = 2)

# ==========================================
# Outcomes
# ==========================================
rrect(xA, yOut, 15, 9, cA, cAb, lwd = 1.5)
text(xA, yOut + 2.4, "A",              col = "white",   cex = 1.30, font = 2)
text(xA, yOut - 0.3, "Global valid IV", col = "#0F4F0F", cex = 0.62, font = 2)
text(xA, yOut - 2.2,
     expression("ATE" %~~% "LATE"),
     col = "#0F4F0F", cex = 0.65)

rrect(xB, yOut, 14, 9, cB, cBb, lwd = 1.5)
text(xB, yOut + 2.4, "B",              col = "white",   cex = 1.30, font = 2)
text(xB, yOut - 0.3, "Local",          col = "#0E2E55", cex = 0.65, font = 2)
text(xB, yOut - 2.2,
     expression("LATE" != "ATE"),
     col = "#0E2E55", cex = 0.65)

rrect(xC, yOut, 16, 9, cC, cCb, lwd = 1.5)
text(xC, yOut + 2.4, "C",               col = "white",   cex = 1.30, font = 2)
text(xC, yOut - 0.3, "IV unreliable",   col = "#5C1612", cex = 0.65, font = 2)
text(xC, yOut - 2.2, "Weak or invalid", col = "#5C1612", cex = 0.62)

dev.off()
cat("Saved: fig0_dial_architecture.png\n")
# =============================================================================
# DIAL_fig2_fig3_v2.R
#
# Regenerates the paper's two empirical figures with the finalized v2 protocol
# data (within-sign domain robustness instead of cross-sign symmetric).
#
# Outputs (in working directory):
#   fig1_diagnostic_space.png    -- Figure 2 in paper (13 instruments, 2-D DIAL space)
#   fig3_domain_sensitivity.png  -- Figure 3 in paper (MCUB bars under two within-sign D per instrument)
#
# Source the finalized Table 6 CIs (post per-instrument canonical correction):
#   Di Tella judgeAlreadyUsedEM:  canonical (0.15, 0.25) -> [-0.109,  0.019]   # CORRECTED
#                                 wide     (0.00, 0.80)  -> [-0.256,  0.638]
#   Di Tella percJudgeSentToEM:   canonical (0.00, 0.25) -> [-0.046,  0.210]   # NEW DOMAIN
#                                 wide     (0.00, 0.80)  -> [-0.046,  0.775]
#   Galiani highnumber:           canonical (-0.20, 0)   -> [-0.189,  0.012]
#                                 wide     (-0.80, 0)    -> [-0.793,  0.012]
#   Banerjee instru:              canonical (-0.60, -0.40)-> [-0.101,  0.163]
#                                 wide     (-0.80, 0)    -> [-0.375,  0.544]
#   Burde buildschool:            canonical (0.00, 0.40) -> [ 0.238,  0.664]
#                                 wide     (0.00, 0.80)  -> [ 0.253,  0.943]
#   Card fatheduc:                canonical (0.00, 0.80) -> [ 0.063,  0.841]
#                                 narrow   (0.00, 0.40)  -> [ 0.063,  0.465]
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# ---- Shared palette --------------------------------------------------------
cA       <- "#3D8A3D"   # Valid A / homogeneous   (green)
cB       <- "#2E5EAA"   # Valid B / heterogeneous (blue)
cC       <- "#C9433D"   # Invalid                 (red)
cValid   <- "#3D8A3D"   # Figure 3: any valid
cInvalid <- "#C9433D"   # Figure 3: invalid

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

# =============================================================================
# FIGURE 2 — DIAL Diagnostic Space (13 instruments, 2-D)
# Changes vs v1:
#   * Banerjee: was (valid=F, C) -> now (valid=T, B) under corrected canonical
#   * Galiani:  was (valid=F, C) -> now (valid=T, B) under its canonical D
#   * 5 valids in total (was 3); 8 invalids (was 10)
# =============================================================================

# NOTE: 'ratio' column is |beta_IV / beta_DML|, which is invariant to domain D.
# These numbers come from Phase 1/Phase 2 runs and are NOT affected by the
# domain corrections — only the `valid` / `dial` columns change.

fig1 <- data.frame(
  label = c(
    "MovieLens\nnum_genres",            # A  valid homogeneous
    "Di Tella\njudgeAlreadyUsedEM",     # B  valid heterogeneous
    "Di Tella\npercJudgeSent",          # A  valid homogeneous
    "Galiani\nhighnumber",              # B  valid heterogeneous  [NEW: was C]
    "Banerjee\ninstru",                 # B  valid heterogeneous  [NEW: was C]
    "401(k)\nsole plan",                # C  invalid
    "401(k)\nplan age",                 # C  invalid
    "NHANES\nincome",                   # C  invalid
    "NHANES\neducation",                # C  invalid
    "Card\nfatheduc",                   # C  invalid
    "Burde\nbuildschool",               # C  invalid
    "MovieLens\nrelease_year",          # C  invalid
    "MovieLens\nis_sequel"              # C  invalid
  ),
  valid = c(rep(TRUE, 5), rep(FALSE, 8)),
  ratio = c(
    0.61, 3.47, 1.33, 2.41, 2.71,       # valids: A, B, A, B, B
    3.40, 1.30, 2.20, 2.40, 1.60,       # invalids cluster 1
    1.59, 1.30, 2.20                    # invalids cluster 2
  ),
  dial  = c("A", "B", "A", "B", "B", rep("C", 8)),
  stringsAsFactors = FALSE
)

# Horizontal positions: valids tight around x = 1; invalids spread across x = 0
set.seed(42)
fig1$x <- NA_real_
n_valid   <- sum( fig1$valid)
n_invalid <- sum(!fig1$valid)
fig1$x[ fig1$valid] <- 1 + runif(n_valid, -0.08, 0.08)
fig1$x[!fig1$valid] <- seq(-0.22, 0.22, length.out = n_invalid)

fig1$nudge_x[fig1$label == "Burde\nbuildschool"]      <- -0.06
fig1$nudge_x[fig1$label == "MovieLens\nrelease_year"] <- -0.04

p1 <- ggplot(fig1, aes(x = x, y = ratio)) +
  # Zone shading
  annotate("rect", xmin = 0.5, xmax = 1.30, ymin = 0.0, ymax = 1.5,
           fill = cA, alpha = 0.08) +
  annotate("rect", xmin = 0.5, xmax = 1.30, ymin = 1.5, ymax = 4.5,
           fill = cB, alpha = 0.06) +
  annotate("rect", xmin = -0.35, xmax = 0.5, ymin = 0.0, ymax = 4.5,
           fill = cC, alpha = 0.05) +
  
  # Threshold lines
  #   Vertical: Invalid | Valid (full height)
  #   Horizontal: A | B boundary — drawn ONLY on the valid side (x >= 0.5),
  #   since R = 1.5 is an A/B distinction that has no meaning in C.
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
  annotate("segment",
           x = 0.5, xend = 1.30, y = 1.5, yend = 1.5,
           linetype = "dashed", color = "gray50", linewidth = 0.3) +
  
  # Zone labels
  annotate("text", x = 0.8, y = 0.30, label = "A: Valid\nhomogeneous",
           size = 3.0, color = cA, fontface = "bold", lineheight = 0.85) +
  annotate("text", x = 0.8, y = 4.25, label = "B: Valid\nheterogeneous",
           size = 3.0, color = cB, fontface = "bold", lineheight = 0.85) +
  annotate("text", x = -0.10, y = 4.25, label = "C: Invalid",
           size = 3.0, color = cC, fontface = "bold") +
  
  # Points
  geom_point(aes(color = dial), size = 3.2, alpha = 0.9) +
  
  geom_text_repel(
    aes(label = label, color = dial,
        size = ifelse(valid, 2.3, 2.1)),   # ✅ move inside aes
    size          = 2.3,
    fontface      = "bold",
    lineheight    = 0.85,
    box.padding   = 0.8,     # ↑ from 0.6
    point.padding = 0.6,     # ↑ from 0.4
    force         = 4,       # ↑ stronger push
    force_pull    = 0.2,     # ↓ weaker pull to point
    min.segment.length = 0.02,
    segment.ncp = -1,
    max.iter      = 10000,
    direction     = "y",
    segment.color = "gray60",
    segment.size  = 0.3,
    max.overlaps  = Inf,
    seed          = 17
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

ggsave("fig1_diagnostic_space.png", p1, width = 7.5, height = 5.3, dpi = 300)
cat("Saved: fig1_diagnostic_space.png\n")


# =============================================================================
# FIGURE 3 — Domain Sensitivity under v2 Within-Sign Protocol
# Changes vs v1:
#   * Cross-sign symmetric (-0.5, 0.5) bars REMOVED
#   * Added wide-canonical bar for Di Tella (both instruments), Banerjee, Burde
#   * Added narrow bar for Card (previously had only canonical)
#   * All same-instrument pairs share a colour (every instrument is robust under v2)
#   * 12 rows (6 instruments x 2 domains) replaces 9 rows from v1
# =============================================================================

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
  lower = c(
    -0.256, -0.109,       # Di Tella judgeAlreadyUsedEM (canonical CORRECTED -0.256 -> -0.109)
    -0.046, -0.046,       # Di Tella percJudgeSent       wide, canonical (now (0, 0.25))
    -0.793, -0.189,       # Galiani                      wide, canonical
    -0.375, -0.101,       # Banerjee                     wide, canonical
    0.253,  0.238,       # Burde                        wide, canonical
    0.063,  0.063        # Card                         canonical, narrow
  ),
  upper = c(
    0.638,  0.019,       # Di Tella judgeAlreadyUsedEM (canonical 0.018 -> 0.019)
    0.775,  0.210,       # Di Tella percJudgeSent
    0.012,  0.012,       # Galiani
    0.544,  0.163,       # Banerjee
    0.943,  0.664,       # Burde
    0.841,  0.465        # Card
  ),
  verdict = c(
    "Valid",   "Valid",   # Di Tella judgeAlreadyUsedEM
    "Valid",   "Valid",   # Di Tella percJudgeSent
    "Valid",   "Valid",   # Galiani
    "Valid",   "Valid",   # Banerjee
    "Invalid", "Invalid", # Burde
    "Invalid", "Invalid"  # Card
  ),
  stringsAsFactors = FALSE
)

# Reverse factor levels so row 1 of the data appears at the top
fig3$row_factor <- factor(fig3$row_label, levels = rev(fig3$row_label))

p3 <- ggplot(fig3, aes(y = row_factor)) +
  # Zero line
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray40", linewidth = 0.4) +
  
  # Identified-set bars
  geom_segment(
    aes(x = lower, xend = upper, yend = row_factor, color = verdict),
    linewidth = 3.2, lineend = "round"
  ) +
  
  # Endpoint ticks
  geom_point(aes(x = lower, color = verdict), size = 1.8) +
  geom_point(aes(x = upper, color = verdict), size = 1.8) +
  
  # rho_ZU = 0 annotation anchored to the dashed vertical
  annotate("text",
           x = 0.02, y = 0.55,
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

ggsave("fig3_domain_sensitivity.png", p3,
       width = 8.2, height = 6.8, dpi = 300)
cat("Saved: fig3_domain_sensitivity.png\n")