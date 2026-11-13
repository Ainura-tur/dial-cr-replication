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
  frame. `bootstrap_rule.R` wraps it to add the MCUB-share rule.
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
04 to 14 do not use it. Its top-level code draws a small example sample
(it sets `N = 800`), so scripts 04 and 05 avoid sourcing it.

### `bootstrap_rule.R`

The MCUB-share rule of Appendix `app:mc` (Remark `rem:mcub_margin`): B
resamples of the raw rows, one MCUB set per resample, and the verdict
p_hat_B >= tau_V. The paper compares this rule with the validity gate; the
comments of scripts 02 to 06 and 10 call it the bootstrap validity rule.
The file also provides the raw-data bootstrap of the correlation
covariance (`dial_cov_corr_boot_raw()`) and the preparation step
(`dial_raw_spec()`, `dial_prep_engine()`, `dial_prep_rho()`), which the
validity gate and scripts 09 to 14 use. It caches each rule call in
`results/boot_rule/`. The header lists its options (`dial.boot_rule`,
`dial.boot_B`, `dial.tau2`, `dial.zero_tol`, `dial.boot_mcub`, and
others).

### `rstar_gate.R`

The validity gate of Algorithm 1 (Proposition `prop:gate`).
`dial_g()` and `dial_rstar()` give the CR mapping and its zero r*;
`dial_raw_rho()` gives the correlations of the residualised (D, Y, Z) for
a resample, and m. `dial_rstar_gate()` returns r*, the confidence set
CS_n, the verdict V, the p-value p_D, the nearest edge, kappa, s(e) and
the Gaussian power prediction; `dial_rstar_gate_data()` does the same from
a data frame. `dial_kappa_star()` gives the frontier. Needs
`dial_raw_spec()` and `dial_prep_engine()` from `bootstrap_rule.R`. Used
by scripts 12, 13 and 14.

### `dml_fallback.R`

DoubleML's partially linear estimator (partialling out, DML2, ranger
forests, standard errors aggregated across sample splits as DoubleML
1.0.2 does) for machines where DoubleML or mlr3 cannot be installed.
`dial_dml_engine_report()` prints which engine runs and, with
`DIAL_REQUIRE_DOUBLEML=1`, stops when DoubleML does not load. Used by
scripts 06, 09 and 11; the reported estimates come from DoubleML.

### `jp_dgp.R`

The Jones-Pewsey data-generating process of Appendix `app:dgp_jp`:
`draw_jp_sample()`, `calibrate_lambda()` and the instrument grid
`JP_INSTRUMENT_GRID`. Used by scripts 04, 05, 10, 13 and 14.

## Source order

Every numbered script sets `DIAL_SKIP_AUTORUN <- TRUE` and sources what it
needs. In an interactive session, set the flag first:

```r
DIAL_SKIP_AUTORUN <- TRUE
source("foundation/Fin_sim3_clean.R")
source("foundation/Fin_Empirical5_clean.R")
source("foundation/bootstrap_rule.R")
source("foundation/rstar_gate.R")      # the validity gate
```

Without the flag, both `Fin_` files start long example runs when sourced.
