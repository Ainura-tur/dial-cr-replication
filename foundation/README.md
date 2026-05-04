# foundation/

This directory contains two R files that form the upstream infrastructure
for the CR test and MCUB inference used throughout the DIAL paper.

## Files

### `Fin_sim3_clean.R`

Implements the core statistical machinery:

- `CIhybrid()`: hybrid confidence interval (projection + conditional)
- `CIproj_p()`, `CIcon()`, `CIcon_TNbounds()`: interval helpers
- `ci_simple_union()`: simple union-bound CI via delta method
- `g_xu_safe()`, `local_compute_gradient_safe()`: the g(delta) mapping and its gradient
- `estimate_cov_corr_boot()`: bootstrap covariance of the correlation triple
- `pvalue_mcub_zero_fast()`: fast p-value for H0: 0 in MCUB identified set
- `generate_data()`, `monte_carlo_simulation()`, `analyze_single_sample()`: simulation framework

### `Fin_Empirical5_clean.R`

Implements the four Phase 1 empirical analysis functions:

- `check_compatibility_simple()`: main compatibility check wrapper
- `run_ditella_analysis_fixed()`: Di Tella and Schargrodsky (2013)
- `run_burde_analysis_fixed()`, `run_burde_sensitivity_fixed()`: Burde and Linden (2013)
- `run_galiani_analysis_fixed()`: Galiani et al. (2011)
- `run_banerjee_analysis_fixed()`: Banerjee et al. (2007)

## Provenance

These files are from the CR test and MCUB inference codebase developed by
Dzhumashev and Tursunalieva (2026) and are bundled here solely to make the
DIAL replication self-contained. They are included under the same license
terms as the original work. If the upstream repository is public at time of
reading, the canonical source is:

  https://github.com/ratbekd/ivcrtest

## Required source order

Every script in `scripts/` that uses these files sets `DIAL_SKIP_AUTORUN <- TRUE`
before sourcing. This flag suppresses the long-running auto-run blocks at the
bottom of both files. **Do not source these files without setting this flag first**
in an interactive session, or all Monte Carlo blocks will execute immediately.

```r
DIAL_SKIP_AUTORUN <- TRUE
source("foundation/Fin_sim3_clean.R")
source("foundation/Fin_Empirical5_clean.R")
```

`Fin_sim3_clean.R` must always be sourced before `Fin_Empirical5_clean.R`
because the empirical file calls functions defined in the simulation file.
