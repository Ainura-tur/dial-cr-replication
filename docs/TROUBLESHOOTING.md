# TROUBLESHOOTING.md

Common errors encountered when running the DIAL replication scripts,
with causes and fixes.

---

## "Required functions not loaded" or "could not find function"

**Cause:** The foundation files were not sourced in the correct order, or
were sourced without setting `DIAL_SKIP_AUTORUN <- TRUE` first.

**Fix:** In any interactive session, always begin with:

```r
DIAL_SKIP_AUTORUN <- TRUE
source("foundation/Fin_sim3_clean.R")
source("foundation/Fin_Empirical5_clean.R")
```

`Fin_sim3_clean.R` must come first. Functions such as `CIhybrid`,
`ci_simple_union`, `g_xu_safe`, and `check_compatibility_simple` are
defined across both files and the empirical file calls functions from
the simulation file.

---

## Long unexpected auto-runs when sourcing foundation files

**Cause:** `DIAL_SKIP_AUTORUN` was not set before sourcing. Both foundation
files check for this flag and run full Monte Carlo blocks if it is absent.

**Fix:** Set the flag before any `source()` call:

```r
DIAL_SKIP_AUTORUN <- TRUE
```

This must be done in the same session before sourcing; setting it after
sourcing has no effect on auto-runs that have already started.

---
## Note on `ratbekd` and `ivcrtest`

The GitHub account `ratbekd` and the `ivcrtest` package are third-party
infrastructure from the CR test paper (Dzhumashev & Tursunalieva 2026).
This account belongs to the upstream package author and is unrelated to
the authors of this submission. See `foundation/README.md` for full
provenance details.
## `iv_cr_test` argument mismatch or rxu_range error

**Cause:** Different versions of the `ivcrtest` package have different
argument names. Some versions do not accept `rxu_range` directly in
`iv_cr_test()`; others do not expose a bootstrap size parameter `B`.

**Fix:** Use `check_compatibility_simple()` from `Fin_Empirical5_clean.R`
rather than calling `iv_cr_test()` directly. The wrapper handles argument
detection automatically. If you must call `iv_cr_test()` directly, run:

```r
names(formals(ivcrtest::iv_cr_test))
```

to see which arguments your installed version accepts.

---

## `multisession` workers crash or produce NA results (script 07)

**Cause:** Too many parallel workers for the available memory, or a
worker-level error that is silently swallowed.

**Fix:** Reduce `N_WORKERS` at the top of `scripts/07_power_curve_grid.R`.
Start with `N_WORKERS <- 2` to confirm correctness, then increase. The
script caches each cell to `results/cp_grid/` so completed cells are not
rerun after a crash.

---

## NHANES data missing columns

**Cause:** The `CrossScreening` package version on your machine has a
different schema for `nhanes.fish`, or the package is not installed.

**Fix:**

```r
install.packages("CrossScreening")
library(CrossScreening)
data(nhanes.fish)
names(nhanes.fish)
```

The script expects: `o.LBXTHG`, `fish`, `income`, `education`, `age`,
`gender`. If column names differ, edit the mapping section at the top of
`scripts/01_build_nhanes.R`.

---

## MovieLens download fails

**Cause:** Network access to `files.grouplens.org` is blocked, or the URL
has changed.

**Fix:** Download `ml-100k.zip` manually from:

  https://grouplens.org/datasets/movielens/100k/

Then unzip it to `data/movielens100k/ml-100k/` and set `extract_dir` in
`run_movielens_analysis_fixed()` to point to that path.

---

## `.dta` file download fails (Di Tella, Galiani, Burde, Banerjee)

**Cause:** Network access to the GitHub mirror is blocked, or the file
has been moved.

**Fix:** Download the `.dta` files manually from:

  https://github.com/ratbekd/Orientation_paper

Place them in the corresponding `data/<study>/` subdirectory and update
the `loader` function URL in `scripts/02_phase1_empirical.R` or
`scripts/cr_wrapper.R` to use a local `haven::read_dta("data/...")` call.

---

## `Matrix::nearPD` warning or error

**Cause:** The bootstrap covariance matrix of the correlation triple is
not positive definite due to small sample size or near-constant variables.

**Fix:** This is usually a data quality issue. Check for constant or
near-constant variables in your instrument or outcome:

```r
apply(df[, c("x","y","z")], 2, var)
```

If variances are near zero, the instrument has no first-stage variation
and the analysis cannot proceed for that sample.

---

## `CIhybrid` produces non-finite values or very wide intervals

**Cause:** The correlation triple `(rho_xy, rho_xz, rho_yz)` is near the
boundary of the feasible region, or the bootstrap covariance estimate is
poorly conditioned with small n.

**Fix:** Increase `B_LARGE` (default 5000) and `B` (default 500) in the
call to `CIhybrid()`, or increase the sample size. For the empirical
applications this is not adjustable, but for simulations you can increase
`N` in `generate_data()`.

---

## Bracketed paste corruption in Git Bash on Windows

**Cause:** Pasting multi-line commands into Git Bash with Ctrl+V inserts
`^[[200~` before the command, corrupting it.

**Fix:** Right-click in the Git Bash window and choose **Paste** from the
context menu instead of using Ctrl+V.
