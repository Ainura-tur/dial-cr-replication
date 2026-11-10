# Reproducibility

This page gives the run order, what each step needs, and how to check
results quickly. Run every script from the `scripts/` folder.

## 1. Setup

```sh
cd scripts
Rscript 00_setup.R
```

`00_setup.R` installs the missing CRAN packages (including DoubleML, mlr3,
mlr3learners and ranger), sources the two main foundation files with
`DIAL_SKIP_AUTORUN <- TRUE`, checks that the other four foundation files
and the four data files are present, and confirms that the key functions
load. It needs R 4.2 or later.

## 2. Run order

`run_all.R` runs the steps in this order, each in a fresh R process:
01, 02, 03, 04, 12, 09, 11, 06, 13, then 08. The scripts of the
MCUB-share rule (05, 07, 10) run only when listed.

| Step | Script | Needs | Writes |
|---|---|---|---|
| 01 | `01_nhanes_build.R` | `CrossScreening` | `scripts/nhanes_mercury.csv` |
| 02 | `02_phase1_empirical.R` | `data/*.dta` | `scripts/phase1_results.csv` |
| 03 | `03_phase2_empirical.R` | 01, `wooldridge`, MovieLens (download or `ml-100k/`) | `scripts/phase2_results.csv` |
| 04 | `04_synthetic_benchmark.R` | `foundation/jp_dgp.R` | `scripts/table1_results.csv` |
| 12 | `12_rstar_gate.R` | 01, 02, 03, `data/*.dta`, `wooldridge`, MovieLens | `results/rstar_gate.csv` |
| 09 | `09_implied_endogeneity.R` | 01, `data/*.dta`, `wooldridge`, MovieLens, DoubleML | `results/implied_endogeneity.csv` |
| 11 | `11_diagnostics.R` | 01, 02, 03, 09, `data/*.dta`, `wooldridge`, MovieLens, DoubleML | `results/diagnostics.csv`, `results/ditella_late_diff.csv`, `results/dml_split_variation.csv` |
| 06 | `06_structural_simulation.R` | DoubleML | `results/sim_dial_v7_results.rds`, `results/coef_vectors.rds` |
| 13 | `13_rstar_gate_mc.R` | the code of 04 and 06 (parts C and D rebuild their samples) | `results/rstar_gate_mc_cells.csv`, `results/rstar_gate_power.csv`, `results/rstar_gate_synthetic.csv`, `results/rstar_gate_sim.csv` |
| 08 | `08_figures.R` | 12, 09, 13 (part B) | `scripts/paper/figures/*.png` |
| 05 | `05_montecarlo_validation.R` | `foundation/jp_dgp.R` | `results/table6_results.csv` |
| 07 | `07_power_curve_grid.R` | nothing else | `results/cp_grid_summary_lean.csv` |
| 10 | `10_rule_montecarlo.R` | 05 (its cache, for the comparison column) | `results/rule_mc_results.csv` |

The steps that need DoubleML fall back to `foundation/dml_fallback.R` when
they run on their own without it. `run_all.R` sets
`DIAL_REQUIRE_DOUBLEML=1` instead, so 09 and 11 stop, and it stops itself
if the log of 06 does not report the DoubleML engine.

### Which script gives which result

| Paper output | Scripts |
|---|---|
| Table 1 (`tab:synthetic`) | 04 (MCUB sets, HSIC); 13, part C (gate) |
| Table 3 (`tab:empirical_full`) | 02 and 03 (MCUB sets); 12 (F, r*, CS_n, V, p_D, m) |
| Table 4 (`tab:locality`) | 09 (r*, R, its interval, B share); 12 (distances) |
| Table 5 (`tab:frontier`) | closed form, section 3 below |
| Table 6 (`tab:frontier_emp`) | 12 |
| Table 7 (`tab:gate_mc`) | 13, part A |
| Table 8 (`tab:sim`) | 06 (estimates, MCUB sets); 13, part D (r*, CS_n, V) |
| Figure 1 (`fig:architecture`) | 08, no inputs |
| Figure 2 (`fig:power_curves`) | 08 from 13, part B |
| Figure 3 (`fig:rstar`) | 08 from 12 |
| Figure 4 (`fig:diagnostic_space`) | 08 from 12 and 09 |
| Numbers in the text | 11: the Di Tella IV difference, the spread of R across sample splits, the scale m |

### Scripts 02 and 03

These two scripts define functions and run nothing when sourced.
`run_all.R` loads the foundation files and calls them for you. To run one
by hand, start a fresh R session in `scripts/` for each:

```r
DIAL_SKIP_AUTORUN <- TRUE
source("../foundation/Fin_sim3_clean.R")
source("../foundation/Fin_Empirical5_clean.R")
source("02_phase1_empirical.R")
res1 <- run_phase1_all()                  # writes phase1_results.csv
```

```r
DIAL_SKIP_AUTORUN <- TRUE
source("../foundation/Fin_sim3_clean.R")
source("../foundation/Fin_Empirical5_clean.R")
source("03_phase2_empirical.R")
res2 <- run_all_phase2_applications()     # writes phase2_results.csv
# with a local MovieLens copy: run_all_phase2_applications(ml_dir = "ml-100k")
```

The flag `DIAL_SKIP_AUTORUN <- TRUE` must come before the first `source()`
call. Without it, the foundation files start their own long simulations.

## 3. Table 5 (`tab:frontier`)

Table 5 gives the gate's frontier from equation `eq:frontier_approx`,
kappa*_N = (z_{1-alpha_0} + z_{1-pi*}) s(e) / sqrt(n), with alpha_0 = 0.05
and s(e) = 1. These lines reproduce it:

```r
kappa_star <- function(n, pi_star, s = 1, alpha0 = 0.05)
  (qnorm(1 - alpha0) + qnorm(1 - pi_star)) * s / sqrt(n)
n <- c(1107, 1534, 5000, 10000, 50000)
round(outer(n, c(0.01, 0.05, 0.10, 0.20), kappa_star), 3)
```

The values scale linearly in s(e). `12_rstar_gate.R` estimates s(e) for
each instrument whose r* lies outside D (column `s_edge`
of `results/rstar_gate.csv`: 0.86 to 1.06 for the instruments that pass
the relevance gate) and gives the frontier at that value (`kappa_star`).

## 4. Caches and resuming

| Folder | Written by | Holds |
|---|---|---|
| `results/boot_rule/` | `foundation/bootstrap_rule.R` (02, 03, 04, 06, 10) | the resampled MCUB set endpoints of every rule call |
| `results/table6_grid_lean/` | 05 | one file per Monte Carlo cell |
| `results/cp_grid_lean/` | 07 | one file per power-curve cell |
| `results/rule_mc/` | 10 | the per-cell results, saved after each replication |
| `results/logs/` | `run_all.R` | one log per step |

An interrupted run of 02 to 07 or 10 resumes where it stopped. The cache
key of a rule call records the sample size, domain, alpha, number of
resamples, seed, MCUB settings, sample correlations and resampling unit,
so a changed setting starts a new entry. Delete a folder after you edit
the MCUB code, or to force a fresh run. Scripts 09, 11, 12 and 13 keep no
cache on disk.

## 5. Quick checks

- `Rscript run_all.R dry` checks the packages, foundation files, data and
  inputs, and lists the steps it would run.
- Script 13: `DIAL_GATE_PARTS=A,C,D` skips part B, the power curves at
  n = 50,000, which take most of its time.
- Script 06: `DIAL_SIM_ROWS=A,Rev,RevHid` runs only those rows of Table 8.
- Script 09: set `B_R <- 0L` to skip the DML refits; r* and its intervals
  then take seconds.
- Script 11: `DIAL_N_SPLITS=0` skips the DML fits of the split
  comparison.
- MCUB-share studies: `options(dial.boot_rule = FALSE)` before sourcing
  02, 03 or 04 skips the rule; in 05, set `R_PER_INSTR <- 50L` for a run
  of about 10 minutes; in 10, lower `R_REPS` (100 per cell) or `B_RULE`
  (200 per replication).

## 6. Expected differences

The MCUB critical values and the gate's standard errors come from
simulation. Rounding differences between computers can therefore move set
endpoints in the third decimal and change bootstrap shares by a few
resamples. Random-forest fits can also differ across `ranger` versions,
which moves the DML estimates slightly, and `dml_fallback.R` gives
slightly different DML estimates from DoubleML. Script 13 runs its cells
in parallel only on Linux and macOS, which changes its run time, not its
results. Run `sessionInfo()` to record the versions you use.
