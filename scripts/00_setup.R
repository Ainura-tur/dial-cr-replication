# =============================================================================
# 00_setup.R
#
# Run this script once before any other script in the repository.
# It installs all required R packages, sets DIAL_SKIP_AUTORUN to prevent
# long-running auto-runs when foundation files are sourced, and runs a
# smoke test confirming the foundation loads without errors.
#
# Usage:
#   Rscript scripts/00_setup.R
#
# Expected output on success:
#   [OK] All packages installed.
#   [OK] foundation/Fin_sim3_clean.R loaded.
#   [OK] foundation/Fin_Empirical5_clean.R loaded.
#   [OK] Setup complete. Proceed to scripts/02_phase1_empirical.R.
# =============================================================================

cat("=== DIAL replication setup ===\n\n")

# ---- 1. CRAN packages -------------------------------------------------------

cran_pkgs <- c(
  # Core statistics and econometrics
  "MASS", "ivreg", "lmtest", "sandwich", "nlme",
  # Data manipulation
  "dplyr", "tidyr", "haven", "zoo",
  # Distributions and simulation
  "fGarch", "e1071",
  # Double machine learning
  "DoubleML", "mlr3", "mlr3learners", "ranger",
  # Bootstrap and parallel
  "future.apply",
  # Independence testing
  "dHSIC",
  # Tables and output
  "stargazer", "knitr",
  # Graphics
  "ggplot2", "ggrepel", "gridExtra", "grid",
  # Wooldridge datasets (Papke 1995, Card 1995)
  "wooldridge",
  # NHANES data
  "CrossScreening",
  # Numerical differentiation (used in delta method)
  "numDeriv",
  # Matrix utilities
  "Matrix",
  # Remote package installation
  "remotes"
)

install_if_missing <- function(pkgs) {
  new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new)) {
    cat("Installing:", paste(new, collapse = ", "), "\n")
    install.packages(new, repos = "https://cloud.r-project.org/")
  } else {
    cat("All CRAN packages already installed.\n")
  }
}

install_if_missing(cran_pkgs)
cat("[OK] All packages installed.\n\n")

# ---- 2. ivcrtest (GitHub) ---------------------------------------------------

if (!requireNamespace("ivcrtest", quietly = TRUE)) {
  cat("Installing ivcrtest from GitHub...\n")
  remotes::install_git("https://github.com/ratbekd/ivcrtest.git")
} else {
  cat("ivcrtest already installed.\n")
}

# ---- 3. DIAL_SKIP_AUTORUN sentinel ------------------------------------------
# All downstream scripts must set this before sourcing the foundation files.
# Setting it here confirms the pattern works and documents the requirement.

DIAL_SKIP_AUTORUN <- TRUE

# ---- 4. Smoke test: source foundation files ---------------------------------

cat("\nSmoke test: sourcing foundation files...\n")

tryCatch({
  source("foundation/Fin_sim3_clean.R")
  cat("[OK] foundation/Fin_sim3_clean.R loaded.\n")
}, error = function(e) {
  stop("[FAIL] foundation/Fin_sim3_clean.R failed to load:\n  ", conditionMessage(e))
})

tryCatch({
  source("foundation/Fin_Empirical5_clean.R")
  cat("[OK] foundation/Fin_Empirical5_clean.R loaded.\n")
}, error = function(e) {
  stop("[FAIL] foundation/Fin_Empirical5_clean.R failed to load:\n  ", conditionMessage(e))
})

# ---- 5. Verify key functions are available ----------------------------------

required_fns <- c(
  "check_compatibility_simple",
  "CIhybrid",
  "ci_simple_union",
  "pvalue_mcub_zero_fast",
  "g_xu_safe",
  "estimate_cov_corr_boot",
  "generate_data",
  "monte_carlo_simulation"
)

missing_fns <- required_fns[!vapply(required_fns, exists, logical(1))]
if (length(missing_fns)) {
  stop("[FAIL] Missing functions after sourcing foundation:\n  ",
       paste(missing_fns, collapse = ", "))
} else {
  cat("[OK] All required functions available.\n")
}

# ---- 6. Write sentinel file -------------------------------------------------

writeLines(as.character(Sys.time()), ".setup_complete")
cat("\n[OK] Setup complete. Proceed to scripts/02_phase1_empirical.R.\n")
cat("     Recommended run order: 01 -> 02 -> 03 -> 04 -> 05 -> 06 -> 07 -> 08\n")
cat("     See docs/REPRODUCIBILITY.md for details.\n\n")
