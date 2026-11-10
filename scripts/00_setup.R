# =============================================================================
# 00_setup.R
#
# Run this script once before the other scripts. It installs the required
# R packages, checks that the foundation files and the bundled data are in
# place, and sources the two main foundation files as a smoke test.
#
# Usage (from the repository root or from scripts/):
#   Rscript scripts/00_setup.R
#
# Expected output on success:
#   [OK] All packages installed.
#   [OK] foundation/Fin_sim3_clean.R loaded.
#   [OK] foundation/Fin_Empirical5_clean.R loaded.
#   [OK] foundation/jp_dgp.R found.
#   [OK] foundation/bootstrap_rule.R found.
#   [OK] foundation/rstar_gate.R found.
#   [OK] foundation/dml_fallback.R found.
#   [OK] All four data files found.
#   [OK] All required functions available.
#   [OK] Setup complete.
# =============================================================================

cat("=== DIAL replication setup ===\n\n")

# ---- 0. Repo-root detection (works under Rscript and interactively) ----------

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
FOUNDATION_DIR <- file.path(.repo_root, "foundation")
if (!dir.exists(FOUNDATION_DIR))
  stop("foundation/ not found under ", .repo_root,
       ". Run this script as Rscript scripts/00_setup.R from the repository.")

# ---- 1. CRAN packages -------------------------------------------------------

cran_pkgs <- c(
  # Core statistics and econometrics
  "MASS", "ivreg", "lmtest", "sandwich", "nlme", "numDeriv", "Matrix",
  # Data manipulation and input
  "dplyr", "data.table", "haven", "zoo",
  # Distributions and simulation
  "fGarch", "e1071",
  # Double machine learning
  "DoubleML", "mlr3", "mlr3learners", "ranger",
  # Parallel computation
  "future", "future.apply",
  # HSIC independence test (synthetic benchmark baseline)
  "dHSIC",
  # Tables and output
  "stargazer", "knitr",
  # Graphics
  "ggplot2", "ggrepel", "patchwork", "gridExtra",
  # Data packages: Papke (1995) and Card (1995); NHANES extract
  "wooldridge", "CrossScreening"
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
.still_missing <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(.still_missing))
  stop("[FAIL] Could not install: ", paste(.still_missing, collapse = ", "))
cat("[OK] All packages installed.\n\n")

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;    filter    <- dplyr::filter
mutate    <- dplyr::mutate;    slice     <- dplyr::slice
recode    <- dplyr::recode;    rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange   <- dplyr::arrange;   count     <- dplyr::count
lag       <- dplyr::lag

# ---- 2. DIAL_SKIP_AUTORUN sentinel ------------------------------------------
# Set this before sourcing the foundation files. Without it, their auto-run
# blocks start long simulations. Every numbered script sets it.

DIAL_SKIP_AUTORUN <- TRUE

# ---- 3. Smoke test: source the two main foundation files --------------------

cat("Smoke test: sourcing foundation files...\n")

for (.f in c("Fin_sim3_clean.R", "Fin_Empirical5_clean.R")) {
  .path <- file.path(FOUNDATION_DIR, .f)
  tryCatch({
    invisible(capture.output(source(.path)))
    cat(sprintf("[OK] foundation/%s loaded.\n", .f))
  }, error = function(e) {
    stop(sprintf("[FAIL] foundation/%s failed to load:\n  %s",
                 .f, conditionMessage(e)))
  })
}
rm(.f, .path)

# ---- 4. The other foundation files ------------------------------------------
# jp_dgp.R: the Jones-Pewsey DGP, used by 04, 05, 10 and 13.
# bootstrap_rule.R: the MCUB-share rule (tau_V) and the preparation step,
#   used by 02 to 06 and 09 to 13.
# rstar_gate.R: the validity gate of Algorithm 1, used by 12 and 13.
# dml_fallback.R: DoubleML's estimator for machines without DoubleML,
#   used by 06, 09 and 11.

for (.f in c("jp_dgp.R", "bootstrap_rule.R", "rstar_gate.R", "dml_fallback.R")) {
  if (file.exists(file.path(FOUNDATION_DIR, .f))) {
    cat(sprintf("[OK] foundation/%s found.\n", .f))
  } else {
    warning(sprintf("foundation/%s not found (file names are case-sensitive).", .f))
  }
}
rm(.f)

# ---- 5. Bundled data ----------------------------------------------------------
# The Phase 1 applications read these files from data/ (see data/README.md).

.data_files <- c("Di_Tella_Schargrodsky_2013.dta", "Crime.dta",
                 "yld_sett_aug03.dta", "afgan.dta")
.data_ok <- vapply(.data_files, function(f)
  !inherits(try(dial_data_path(f), silent = TRUE), "try-error"), logical(1))
if (all(.data_ok)) {
  cat("[OK] All four data files found.\n")
} else {
  warning("Missing data files in data/: ",
          paste(.data_files[!.data_ok], collapse = ", "))
}
rm(.data_files, .data_ok)

# ---- 6. Verify key functions are available ----------------------------------

required_fns <- c(
  "check_compatibility_simple",
  "CIhybrid",
  "ci_simple_union",
  "pvalue_mcub_zero_fast",
  "g_xu_safe",
  "estimate_cov_corr_boot",
  "dial_data_path",
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

# ---- 7. Write sentinel file -------------------------------------------------

writeLines(as.character(Sys.time()), file.path(.repo_root, ".setup_complete"))
cat("\n[OK] Setup complete.\n")
cat("     Next, from the scripts folder:\n")
cat("       Rscript run_all.R paper   every script behind the paper (about 6-7 hours)\n")
cat("       Rscript run_all.R dry     check the setup and list the steps only\n")
cat("     See README.md and docs/REPRODUCIBILITY.md for details.\n\n")
