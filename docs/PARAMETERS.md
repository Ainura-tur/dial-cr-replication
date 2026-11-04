# Parameters

The settings behind every table and figure, with the script that sets
each one. Each script's header repeats its own settings.

## Decision rule (Algorithm 1)

| Setting | Value | Where |
|---|---|---|
| Relevance gate | first-stage F < 10 gives C (weak) | 03, 06, 09 |
| Validity threshold tau_V | 0.95: valid if p_hat_B >= tau_V | `foundation/bootstrap_rule.R`; option `dial.tau2`, column `Valid_tau2` |
| Bootstrap resamples B | 500 (200 in 10) | option `dial.boot_B`; `B_RULE` in 04, 06, 10 |
| Ratio rule | A if the IV and DML estimates share a sign and abs(log R) < log tau_R; otherwise B, labelled "B (sign)" when the signs differ | 03, 06, 09 |
| tau_R | 1.5, so A needs 1/1.5 < R < 1.5 | `TAU_R` in 06, 08, 09 |
| Margin on the set endpoints | 0 (strict) in the empirical tables; 0.05 in the simulations | option `dial.zero_tol`; `ZERO_TOL` in 04, 05, 06, 10 |

In the code, the threshold tau_V is called `tau2` (`TAU2` in 04 and 07).

## MCUB sets

| Setting | Full-sample set | Each bootstrap resample of the rule |
|---|---|---|
| Draws per bisection step | 500 | 500 |
| Coverage draws | 5,000 | 499 |
| Covariance of the correlations | bootstrap, 800 raw-data draws | bootstrap, 200 raw-data draws |
| Nominal level alpha_0 | 0.05 (empirical), 0.10 (simulations) | as the full-sample set |
| Seed | per instrument | resample b of instrument i uses (5000 + i) * 10000 + b |

Every bootstrap redraws the raw rows and redoes the preparation step:
residualise Y, D and Z on the covariates, then replace the residualised
instrument by its fitted value from a regression on the residualised
(D, Y). `CIhybrid()` uses eta = 0.001, alpha_c = 0.8 alpha_0 and
bisection tolerances of 1e-3.

## DML

Partially linear model, random forests (`ranger`), minimum node size 5,
5-fold cross-fitting. The empirical applications use 500 trees and three
repetitions of the sample split; the structural simulation (06) uses 300
trees. The bootstrap refits of 09 use one repetition and assign the folds
to the original rows.

## Maintained domains D for rho_DU

| Instrument | Canonical | Within-sign wide |
|---|---|---|
| Di Tella `judgeAlreadyUsedEM` | (0.15, 0.25) | (0, 0.8) |
| Di Tella `percJudgeSentToEM` | (0, 0.25) | (0, 0.8) |
| Galiani `highnumber` | (-0.2, 0) | (-0.8, 0) |
| Banerjee `instru` | (-0.6, -0.4) | (-0.8, 0) |
| Burde `buildschool` | (0, 0.4) | (0, 0.8) |
| MovieLens `num_genres` | (0, 0.8) | |
| MovieLens `release_year` | (-0.8, 0) | |
| MovieLens `is_sequel` | (-0.8, 0) | |
| Papke `sole_plan`, `plan_age` | (0, 0.3) | (0, 0.8) |
| NHANES `income`, `education` | (0, 0.8) | |
| Card `fatheduc` | (0, 0.4) | (0, 0.8) |

`run_burde_sensitivity_fixed()` also tests Burde under (-0.4, 0), the
opposite sign. It is a diagnostic and does not enter the tables.

## Implied endogeneity and the ratio R (09)

| Setting | Value |
|---|---|
| Resamples for the r* interval | 2,000 (`B_RSTAR`) |
| Resamples for the R interval | 200 (`B_R`), with a DML refit on each |
| DML seed | one per instrument (1101 to 1110; 43 to 45 for MovieLens) |

## Simulations

| Setting | 04 (`tab:synthetic`) | 05 (`tab:mc`) | 10 (`tab:rule_mc`) |
|---|---|---|---|
| DGP | Jones-Pewsey (`foundation/jp_dgp.R`) | same | same |
| n | 800 | 800 | 800 |
| Instruments | 11 at cor(z, u) from -0.20 to 0.20, plus Z_het | the same 11 | cells 0, 0.12, 0.20 |
| Replications | one sample | 1,000 per instrument | 100 per cell |
| target abs(cor(x, u)), cor(z, x) | 0.30, 0.25 | 0.30, 0.25 | 0.30, 0.25 |
| Domain | (0, 0.8) if cor(z, u) >= 0, else (-0.8, 0) | same | same |
| alpha_0 | 0.10 | 0.10 | 0.10 |
| Master seed | 1234 | 4242 | 4242 |

Script 05 calls `CIhybrid()` once per replication with 500 bisection
draws, 499 coverage draws and a 200-draw raw-data covariance. Its truth
labels are A for abs(rho_ZU) <= 0.04, boundary for 0.08 and 0.12, and C
from 0.16. The error bound is the exact one-sided 95% Clopper-Pearson
bound, compared with pi* = 0.05.

### Structural simulation (06, `tab:sim`)

| Setting | Value |
|---|---|
| n | 50,000 per row |
| Controls H | 5 standard normals; coefficients phi and lambda drawn once with seed 2026 |
| Confounder C | loading 0.40 in both V and U; error SD 0.6 |
| Rows A, B (conf.), C, C (weak) | beta = 0.425 |
| Row B (het.) | beta_H = 0.80, beta_L = 0.05; pi_H = 0.05, pi_L = 0.60 |
| First stage | pi = 0.25; weak row pi = sqrt(5 (0.4^2 + 0.6^2) / n), so F is about 7 |
| Exclusion violation (row C) | gamma = 0.30 |
| Domain, alpha_0, margin | (0, 0.8), 0.10, 0.05 |
| Data seeds | 101 (A and B (conf.)), 102 (B (het.)), 201 (C), 401 (C (weak)) |

### Power curves (07, `fig:power_curves`)

| Setting | Value |
|---|---|
| Sample sizes | 1,107; 1,534; 5,000; 10,000; 50,000 |
| delta grid | 35 base points in [0.001, 0.40] plus 8 points in the transition region |
| Replications per cell | 250, 250, 150, 100, 50 for the five sample sizes |
| MCUB | 500 bisection draws, 499 coverage draws, 200-draw covariance |
| DGP | rho_DU = 0.50, cor(Z, D) = 0.25, beta = 0.30 |
| Domain, alpha_0, tau_V | (0, 0.8), 0.05, 0.95 |

## Classifiability frontier (`tab:frontier`)

sigma = 3.25, the median of the values calibrated on the eight invalid
empirical instruments (range 1.69 to 6.75), with tau_V = 0.95, B = 500
and alpha_0 = 0.05. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md),
section 3, for the formula in R.
