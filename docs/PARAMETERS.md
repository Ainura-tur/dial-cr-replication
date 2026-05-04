# PARAMETERS.md

All hardcoded constants used across the DIAL replication scripts, their
values, and where they appear. Reviewers can cross-check these against
the paper's appendices and Section 5.

---

## DGP parameters (simulation and power curve scripts)

| Constant | Value | Description | Script(s) |
|----------|-------|-------------|-----------|
| `RHO_DU` | 0.50 | Correlation between treatment D and unobservable U | 04, 07 |
| `COR_ZD` | 0.25 | Instrument-treatment correlation (first-stage strength) | 04, 07 |
| `BETA_TRUE` | 0.30 | True causal effect of D on Y | 04, 07 |
| `BETA_HIGH` | 0.80 | High-type causal effect (heterogeneous DGP) | 06 |
| `BETA_LOW` | 0.05 | Low-type causal effect (heterogeneous DGP) | 06 |
| `PI_Z` | 0.25 | First-stage coefficient in structural simulation | 06 |
| `GAMMA_INVALID` | 0.30 | Exclusion violation magnitude (invalid scenario) | 06 |
| `GAMMA_BORDER` | 0.12 | Exclusion violation magnitude (borderline scenario) | 06 |
| `PI_WEAK` | 0.02 | First-stage coefficient producing F < 10 (weak scenario) | 06 |

---

## Inference parameters

| Constant | Value | Description | Script(s) |
|----------|-------|-------------|-----------|
| `ALPHA` | 0.05 | Nominal significance level | all |
| `ETA` | 0.001 | Hybrid CI tuning parameter | all |
| `TAU2` | 0.95 | Decision threshold for CP curve (p >= TAU2 -> valid) | 07 |
| `PI_STAR` | 0.05 | Target size for the CR test | 07 |
| `B_BOOT` | 500 | Bootstrap draws for MCUB CI (standard) | 04, 05 |
| `B_LARGE` | 5000 | Bootstrap draws for large-sample MCUB check | 04, 05 |
| `B_FAST` | 600 | Bootstrap draws for fast p-value computation | all |
| `B_LARGE_FAST` | 6000 | Large-sample draws for fast p-value computation | all |
| `TOL` | 1e-3 | Bisection tolerance for CI computation | all |
| `TOL_R` | 1e-3 | Truncated-normal bound tolerance | all |
| `ALPHA_C` | 0.8 * ALPHA | Conditional coverage level in hybrid CI | all |

---

## Domain (rxu_range) parameters

These define the assumed range of the correlation between the instrument
residual Z and the unobservable U. Domain choice is the key sensitivity
parameter in the DIAL framework.

| Study | Canonical domain | Notes |
|-------|-----------------|-------|
| Di Tella and Schargrodsky | c(0, 0.25) | Positive endogeneity; narrow |
| Burde and Linden | c(0.0, 0.4) | Positive; pending structural justification for negative |
| Galiani et al. | c(-0.2, 0.0) | Negative endogeneity (conscription lottery) |
| Banerjee et al. | c(-0.6, -0.4) | Negative; corrected from c(-0.8, -0.6) in Phase 1 |
| Papke (1995) | c(0.0, 0.3) | Positive; narrower than wide c(0, 0.8) robustness check |
| NHANES | c(0.0, 0.8) | Positive endogeneity (healthy-eating channel) |
| Card (1995) canonical | c(0.0, 0.8) | Contested direction; see symmetric check |
| Card (1995) symmetric | c(0.0, 0.4) | Robustness check |
| MovieLens num_genres | c(0.0, 0.8) | Positive endogeneity |
| MovieLens release_year | c(-0.8, 0.0) | Negative endogeneity (recency bias) |
| MovieLens is_sequel | c(-0.8, 0.0) | Negative endogeneity (franchise brand capital) |
| Simulation (symmetric default) | c(-0.4, 0.4) | Used in cr_wrapper.R scenarios 1-3, 5 |
| Synthetic Table 1 | c(0.01, 0.80) | Used in 04_synthetic_table1.R |

---

## Grid parameters (power curve script 07)

| Constant | Value | Description |
|----------|-------|-------------|
| `N_GRID` | c(1107, 1534, 5000, 10000, 50000) | Sample sizes on the CP curve grid |
| `DELTA_GRID` | 35 points in [0.001, 0.35] | Exclusion violation magnitudes |
| `R_MC` (n=1107, 1534) | 500 | Monte Carlo reps per cell |
| `R_MC` (n=5000) | 300 | Monte Carlo reps per cell |
| `R_MC` (n=10000) | 200 | Monte Carlo reps per cell |
| `R_MC` (n=50000) | 100 | Monte Carlo reps per cell |
| `N_WORKERS` | 7 | Parallel workers (multisession) |
| `B_BOOT_TARGET` | 199 | Bootstrap override if ivcrtest exposes B argument |
| `RXU_RANGE` | c(0.001, 0.80) | Domain for power curve CR test calls |

---

## Random seeds

| Seed value | Location | Purpose |
|------------|----------|---------|
| 2026 | 04_synthetic_table1.R, 06_structural_table2.R | Global set.seed |
| 123 | estimate_cov_corr_boot(), ci_simple_union() | Bootstrap covariance |
| 42 | 04_synthetic_table1.R DGP calls | Per-scenario data generation |
| 9058 | cr_wrapper.R, Fin_sim3_clean.R auto-run | Simulation seed base |
| 1000 + i | analyze_single_sample() deltaSigma bootstrap | Per-instrument seed |
| 1000*n + d_idx | 07_power_curve_grid.R | Per-cell seed base |
| seed + as.integer(alpha * 1e6) | pvalue_mcub_zero_fast() | Per-alpha-level seed |

---

## Effect ratio thresholds (DIAL classification)

| Threshold | Value | Meaning |
|-----------|-------|---------|
| `TAU_R` lower | 0.5 | Below this, ratio is too small for zone A |
| `TAU_R` upper | 1.5 | Above this, classified as zone B (heterogeneous) |
| F-statistic | 10 | Weak instrument gate threshold |
