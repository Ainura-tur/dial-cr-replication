# Reproducibility

This page gives the run order, what each step needs, and how to check
results quickly. Run every script from the `scripts/` folder.

## 1. Setup

```sh
cd scripts
Rscript 00_setup.R
```

`00_setup.R` installs the missing CRAN packages, sources the two main
foundation files with `DIAL_SKIP_AUTORUN <- TRUE`, checks that
`foundation/jp_dgp.R`, `foundation/bootstrap_rule.R` and the four data
files are present, and confirms that the key functions load. It needs R
4.2 or later.

## 2. Run order

| Step | Script | Needs | Writes |
|---|---|---|---|
| 1 | `01_nhanes_build.R` | `CrossScreening` | `scripts/nhanes_mercury.csv` |
| 2 | `02_phase1_empirical.R` | `data/*.dta` | `scripts/phase1_results.csv` |
| 3 | `03_phase2_empirical.R` | step 1, `wooldridge`, MovieLens (download or `ml-100k/`) | `scripts/phase2_results.csv` |
| 4 | `04_synthetic_benchmark.R` | `foundation/jp_dgp.R` | `scripts/table1_results.csv` |
| 5 | `05_montecarlo_validation.R` | `foundation/jp_dgp.R` | `results/table6_results.csv` |
| 6 | `06_structural_simulation.R` | nothing else | `results/sim_dial_v7_results.rds` |
| 7 | `07_power_curve_grid.R` | nothing else | `results/cp_grid_summary_lean.csv` |
| 8 | `08_figures.R` | step 7 | `scripts/paper/figures/*.png` |
| 9 | `09_implied_endogeneity.R` | step 1, `data/*.dta`, `wooldridge`, MovieLens | `results/implied_endogeneity.csv` |
| 10 | `10_rule_montecarlo.R` | step 5 (its cache, for the comparison column) | `results/rule_mc_results.csv` |

Scripts 04 to 10 source the foundation files themselves and run with
`Rscript`. Figures 3 and 4 use values copied into `08_figures.R` from
the runs of 02, 03 and 09, so 08 needs only the output of 07.

### Scripts 02 and 03

These two scripts define functions and run nothing when sourced. Use the
one-liners in the main README, or run them in an R session. Start a fresh
session for each.

```r
setwd("path/to/dial-cr-replication/scripts")
DIAL_SKIP_AUTORUN <- TRUE
source("../foundation/Fin_sim3_clean.R")
source("../foundation/Fin_Empirical5_clean.R")
source("02_phase1_empirical.R")
res1 <- run_phase1_all()                  # writes phase1_results.csv
```

```r
setwd("path/to/dial-cr-replication/scripts")
DIAL_SKIP_AUTORUN <- TRUE
source("../foundation/Fin_sim3_clean.R")
source("../foundation/Fin_Empirical5_clean.R")
source("03_phase2_empirical.R")
res2 <- run_all_phase2_applications()     # writes phase2_results.csv
# with a local MovieLens copy: run_all_phase2_applications(ml_dir = "ml-100k")
```

The flag `DIAL_SKIP_AUTORUN <- TRUE` must come before the first `source()`
call. Without it, the foundation files start their own long simulations.

## 3. Table 6 (`tab:frontier`)

Table 6 comes from the closed-form frontier of Remark
`rem:frontier_approx`, with tau_V = 0.95, B = 500, alpha_0 = 0.05 and
sigma = 3.25 (the median of the values calibrated on the eight invalid
empirical instruments, Appendix `app:frontier`). These lines reproduce it:

```r
kappa_star <- function(n, pi_star, sigma = 3.25, tau_V = 0.95, B = 500, alpha0 = 0.05)
  sigma / sqrt(n) * (qnorm(1 - alpha0 / 2) -
    qnorm(tau_V - qnorm(1 - pi_star) * sqrt(tau_V * (1 - tau_V) / B)))
n <- c(1107, 1534, 5000, 10000, 50000)
round(outer(n, c(0.01, 0.05, 0.10, 0.20), kappa_star), 3)
```

## 4. Caches and resuming

| Folder | Written by | Holds |
|---|---|---|
| `results/boot_rule/` | `foundation/bootstrap_rule.R` (02, 03, 04, 06, 10) | the resampled MCUB set endpoints of every rule call |
| `results/table6_grid_lean/` | 05 | one file per Monte Carlo cell |
| `results/cp_grid_lean/` | 07 | one file per power-curve cell |
| `results/rule_mc/` | 10 | the per-cell results, saved after each replication |

An interrupted run resumes where it stopped. The cache key of a rule call
records the sample size, domain, alpha, number of resamples, seed, MCUB
settings, sample correlations and resampling unit, so a changed setting
starts a new entry. Delete a folder after you edit the MCUB code, or to
force a fresh run.

## 5. Quick checks

- Bootstrap rule off: `options(dial.boot_rule = FALSE)` before sourcing
  02, 03 or 04 gives the single-sample CR test only. For 02 and 03 this
  removes about 70 of the 120 minutes.
- Script 05: set `R_PER_INSTR <- 50L` for a run of about 10 minutes.
- Script 09: set `B_R <- 0L` to skip the DML refits; r* and its intervals
  then take seconds.
- Script 10: lower `R_REPS` (100 per cell) or `B_RULE` (200 per
  replication).

## 6. Expected differences

The MCUB critical values come from simulation. Rounding differences
between computers can therefore move set endpoints in the third decimal
and change the bootstrap share p_hat_B by a few resamples. Random-forest
fits can also differ across `ranger` versions, which moves the DML
estimates slightly. Run `sessionInfo()` to record the versions you use.
