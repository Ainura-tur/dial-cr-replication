# foundation/

Shared code that the numbered scripts source. With
`DIAL_SKIP_AUTORUN <- TRUE`, sourcing these files defines functions and
starts no simulations (`Fin_sim3_clean.R` still fits one small example at
the top).

## Files

### `Fin_Empirical5_clean.R`

The CR test with MCUB inference, and the four Phase 1 applications.

- `g_fun()`, `g_grad()`, `g_xu_safe()`, `local_compute_gradient_safe()`:
  the CR map g(rho, r) and its gradient.
- `estimate_cov_corr_boot()`: bootstrap covariance of the three sample
  correlations.
- `CIhybrid()`, `CIproj_p()`, `CIcon()`, `CIcon_TNbounds()`: the MCUB
  confidence set (hybrid of projection and conditional sets).
- `ci_simple_union()`, `pvalue_mcub_zero_fast()`: a union-bound set and a
  fast p-value for zero in the MCUB set.
- `check_compatibility_simple()`: the CR test on a prepared (x, y, z)
  frame. `bootstrap_rule.R` wraps it to add the bootstrap validity rule.
- `run_ditella_analysis_fixed()`, `run_galiani_analysis_fixed()`,
  `run_banerjee_analysis_fixed()`, `run_burde_analysis_fixed()`,
  `run_burde_sensitivity_fixed()`, `run_all_applications()`: the Phase 1
  applications.
- `dial_data_path()`: finds a bundled data file in `data/`.

### `Fin_sim3_clean.R`

An earlier copy of the MCUB functions, plus a simulation framework
(`generate_data()`, `analyze_single_sample()`, `monte_carlo_simulation()`,
`summarize_sim_results()`). Scripts 02 and 03 source it before
`Fin_Empirical5_clean.R`, whose definitions then take precedence. Scripts
04 to 10 do not use it. Its top-level code draws a small example sample
(it sets `N = 800`), so scripts 04 and 05 avoid sourcing it.

### `bootstrap_rule.R`

The bootstrap validity rule of Algorithm 1 (Proposition
`prop:finite_sample`): B resamples of the raw rows, one MCUB set per
resample, and the verdict p_hat_B >= tau_V. It also provides the raw-data
bootstrap of the correlation covariance (`dial_cov_corr_boot_raw()`) and
the preparation step (`dial_raw_spec()`, `dial_prep_engine()`,
`dial_prep_rho()`). It caches each rule call in `results/boot_rule/`.
The header lists its options (`dial.boot_rule`, `dial.boot_B`,
`dial.tau2`, `dial.zero_tol`, `dial.boot_mcub`, and others).

### `jp_dgp.R`

The Jones-Pewsey data-generating process of Appendix `app:dgp_jp`:
`draw_jp_sample()`, `calibrate_lambda()` and the instrument grid
`JP_INSTRUMENT_GRID`. Used by scripts 04, 05 and 10.

## Source order

Every numbered script sets `DIAL_SKIP_AUTORUN <- TRUE` and sources what it
needs. In an interactive session, set the flag first:

```r
DIAL_SKIP_AUTORUN <- TRUE
source("foundation/Fin_sim3_clean.R")
source("foundation/Fin_Empirical5_clean.R")
source("foundation/bootstrap_rule.R")
```

Without the flag, both `Fin_` files start long example runs when sourced.
