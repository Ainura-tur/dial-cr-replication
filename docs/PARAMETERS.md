# Parameters

The settings behind every table and figure, with the script that sets
each one. Each script's header repeats its own settings.

## Decision rule (Algorithm 1)

| Setting | Value | Where |
|---|---|---|
| Relevance gate | first-stage F < 10 gives C (weak) | `F_MIN` in 12; also 03 and 06 |
| Validity gate | C unless the 90% confidence set CS_n for r* meets the closed domain D | `foundation/rstar_gate.R`, run by 12, 13 and 14 |
| Level alpha_0 of the gate | 0.05, so CS_n has coverage 0.90 | `ALPHA0` in 12, 13 and 14 |
| Resamples for se_n(r) | 1,000 for the empirical instruments and Tables 1 and 9; 500 in Tables 7 and 8; 300 for the power curves | `B_GATE` in 12 and 14; `B_CD`, `B_A`, `B_B` in 13 |
| Grid for r | step 0.001 in 12 and for Table 9 (13, part D); 0.002 for Tables 1, 7 and 8 and the power curves; the edges of D added | `STEP` in 12, `STEP_MC` in 13 and 14 |
| Ratio rule | a resample lands in B when abs(log R) >= log tau_R or the IV and DML estimates differ in sign | 06, 09, 11 |
| tau_R | 1.5, so A needs 1/1.5 < R < 1.5 | `TAU_R` in 06, 08, 09, 11 |
| Agreement label | A if at most 10% of the 200 resamples land in B, B if at least 90% do, inconclusive otherwise | `B_R` in 09; `LABEL_LO`, `LABEL_HI` in 11; `TAU_L` in 08 |

## The MCUB-share rule (Appendix `app:mc`)

Scripts 02 to 06 and 10 also apply this rule, which the paper compares with
the gate (Remark `rem:mcub_margin`). Their comments call it the bootstrap
validity rule or the tau_V rule.

| Setting | Value | Where |
|---|---|---|
| Threshold tau_V | 0.95: valid if p_hat_B >= tau_V | `foundation/bootstrap_rule.R`; option `dial.tau2`, column `Valid_tau2` |
| Bootstrap resamples B | 500 (200 in 10) | option `dial.boot_B`; `B_RULE` in 04, 06, 10 |
| Margin on the set endpoints | 0 (strict) in the empirical tables; 0.05 in the simulations | option `dial.zero_tol`; `ZERO_TOL` in 04, 05, 06, 10 |

In the code, the threshold tau_V is called `tau2` (`TAU2` in 04 and 07).

## MCUB sets

| Setting | Full-sample set | Each bootstrap resample of the MCUB-share rule |
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
5-fold cross-fitting, fitted with DoubleML (`foundation/dml_fallback.R`
implements the same estimator for machines without DoubleML). The
empirical applications use 500 trees and three repetitions of the sample
split; the structural simulation (06) uses 300 trees. The bootstrap refits
of 09 use one repetition and assign the folds to the original rows. The
spread of R in 11 uses ten fits with one split each.

## Maintained domains D for rho_DU

The paper uses the canonical domains. The within-sign wide domains enter
only the domain-sensitivity output of 02 and 03, which the paper does not
report.

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
| Resamples for the IV intervals and the Di Tella difference (11) | 2,000 (`B_IV`) |

## Monte Carlo of the gate (13)

| Setting | Part A (Table 7) | Part B (Figure 3) |
|---|---|---|
| Design | rho_DU = 0.30, cor(Z, D) = 0.25 | rho_DU = 0.50, cor(Z, D) = 0.25 |
| delta = cor(Z, U) | 0, 0.04, ..., 0.20, plus an edge cell (rho_DU = 0, delta = 0) | 13 values from 0 to 0.25 |
| n | 800 | 1,107; 5,000; 50,000 |
| Samples per cell | 500 | 200 |
| Domain, alpha_0 | (0, 0.8), 0.05 | (0, 0.8), 0.05 |
| Master seed | 4242 | 4242 |

Parts C and D apply the gate to the samples of Tables 1 and 9, which they
rebuild by running the data-generating code of 04 and 06.

## HSIC and the gate over repeated samples (14, `tab:hsic_mc`)

| Setting | Value |
|---|---|
| Design | population version of the Jones-Pewsey design of 04: cor(Z, X) = 0.25, cor(X, U) = 0.30 or -0.30, Y = 2X + U; fixed population coefficients and fresh noise in every sample |
| Instruments | the 11 grid instruments, cor(Z, U) from -0.20 to 0.20, and Z_het |
| Domain | (-0.8, 0) with cor(X, U) = -0.30 when cor(Z, U) < 0; otherwise (0, 0.8) with cor(X, U) = 0.30 |
| n | 800 |
| Samples per instrument | 500 (`R_REPS`; the environment variable `DIAL_HSIC_REPS` overrides it) |
| Gate | alpha_0 = 0.05, B = 500, grid step 0.002 (`ALPHA0`, `B_GATE`, `STEP_MC`) |
| HSIC | `dHSIC::dhsic.test()` with the gamma approximation at the 5% level, testing z against e_hat = y - beta_IV x, as in 04 |
| Seeds | master seed 5151; sample i of instrument k uses 5151 + 10000 k + i, and its gate uses that seed plus 1,000,000 |

The results do not depend on the number of cores. `OPPOSITE_SIGN <- TRUE`
adds the instruments with cor(Z, U) < 0 under cor(X, U) = 0.30 and
D = (0, 0.8), the direction in which the gate has little power; the paper
does not report these rows.

## Simulations

| Setting | 04 (`tab:synthetic`) | 05 (MCUB-share study) | 10 (MCUB-share study) |
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

### Power grid of single MCUB sets (07, MCUB-share study)

| Setting | Value |
|---|---|
| Sample sizes | 1,107; 1,534; 5,000; 10,000; 50,000 |
| delta grid | 35 base points in [0.001, 0.40] plus 8 points in the transition region |
| Replications per cell | 250, 250, 150, 100, 50 for the five sample sizes |
| MCUB | 500 bisection draws, 499 coverage draws, 200-draw covariance |
| DGP | rho_DU = 0.50, cor(Z, D) = 0.25, beta = 0.30 |
| Domain, alpha_0, tau_V | (0, 0.8), 0.05, 0.95 |

## Classifiability frontier (`tab:frontier`)

The gate's frontier kappa*_N = (z_{1-alpha_0} + z_{1-pi*}) s(e) / sqrt(n),
with alpha_0 = 0.05 and s(e) = 1. The empirical instruments that pass the
relevance gate have s(e) between 0.86 and 1.06 (`s_edge` in
`results/rstar_gate.csv`). See [REPRODUCIBILITY.md](REPRODUCIBILITY.md),
section 3, for the formula in R.
