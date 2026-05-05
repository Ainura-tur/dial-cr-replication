# REPRODUCIBILITY.md

Step-by-step instructions for reproducing all tables and figures in the
DIAL paper from a fresh clone of this repository.

---

## Prerequisites

- R >= 4.2.0
- Git Bash or a Unix-like terminal
- Internet access (for package installation and dataset downloads)
- 7-core machine recommended for the power curve grid (script 07)

---

## Step 0: Clone and set up

```r
git clone https://anonymous.4open.science/r/dial-cr-replication.git
cd dial-cr-replication
Rscript scripts/00_setup.R
```

This installs all CRAN and GitHub packages, sources the foundation files as
a smoke test, and writes `.setup_complete`. If any step fails, see
`docs/TROUBLESHOOTING.md`.

---

## Step 1: Build the NHANES dataset

```r
Rscript scripts/01_build_nhanes.R
```

Requires the `CrossScreening` package and `nhanes.fish` data object.
Writes `nhanes_mercury.csv` to the working directory. Must be run before
script 03.

**Output:** `nhanes_mercury.csv`  
**Runtime:** under 1 minute

---

## Step 2: Phase 1 empirical re-runs (Tables 4-5, rows 1-4)

```r
Rscript scripts/02_phase1_empirical.R
```

Runs the four Phase 1 studies: Di Tella and Schargrodsky, Galiani et al.,
Burde and Linden, and Banerjee et al. Downloads `.dta` files at runtime.

**Output:** `results/phase1_results.csv`  
**Runtime:** 20-40 minutes

---

## Step 3: Phase 2 empirical re-runs (Tables 4-5, rows 5-8)

```r
Rscript scripts/03_phase2_empirical.R
```

Requires `nhanes_mercury.csv` from script 01. Runs Papke (1995), NHANES,
Card (1995), and MovieLens 100K. Downloads MovieLens at runtime.

**Output:** printed summary; update Table 5 rows 5-8 manually from console output  
**Runtime:** 20-40 minutes

---

## Step 4: Table 1 (synthetic MCUB identified sets)

```r
Rscript scripts/04_synthetic_table1.R
```

Generates the six MCUB identified sets for the six synthetic instrument
scenarios in Table 1. Prints paste-ready LaTeX snippets.

**Output:** `results/table1_results.csv`  
**Runtime:** 5 minutes

---

## Step 5: Table 6 (Monte Carlo coverage)

```r
Rscript scripts/05_montecarlo_table6.R
```

Runs the full Monte Carlo simulation grid. This is the longest script.
Results are cached per cell in `results/table6_grid/` so interrupted runs
can be resumed without restarting from scratch.

**Output:** `results/table6_results.csv`  
**Runtime:** 4-8 hours

---

## Step 6: Table 2 (structural simulation, DIAL scenarios)

```r
Rscript scripts/06_structural_table7.R
```

Runs the four DIAL structural scenarios using the DoubleML and ivcrtest
pipeline.

**Output:** `paper/figures/sim_dial_scenarios.png`, `paper/figures/sim_dial_estimates.png`  
**Runtime:** approximately 3 minutes per scenario row

---

## Step 7: Power curve grid (Figure 2)

```r
Rscript scripts/07_power_curve_grid.R
```

Parallelised Monte Carlo over a 5 x 35 grid of (n, delta) cells. Uses 7
workers by default; adjust `N_WORKERS` at the top of the script. Results
cached per cell in `results/cp_grid/`.

**Output:** `results/cp_grid_summary.csv`  
**Runtime:** approximately 75 minutes on 7 cores

---

## Step 8: Figures (Figures 1, 3, 4)

```r
Rscript scripts/08_figures.R
```

Produces the architecture diagram, DIAL diagnostic space plot, and domain
sensitivity figure. Requires no upstream script outputs.

**Output:**
- `paper/figures/fig0_dial_architecture.png`
- `paper/figures/fig1_diagnostic_space.png`
- `paper/figures/fig3_domain_sensitivity.png`

**Runtime:** under 2 minutes

---

## Dependency graph

```
00_setup
    |
    +-- 01_build_nhanes --> 03_phase2_empirical
    |
    +-- 02_phase1_empirical
    |
    +-- 04_synthetic_table1
    |
    +-- 05_montecarlo_table6
    |
    +-- 06_structural_table2
    |
    +-- 07_power_curve_grid
    |
    +-- 08_figures
```

Scripts 02, 04, 05, 06, 07, and 08 are independent of each other and can
be run in any order after setup. Script 03 requires script 01 first.

---

## Interactive session

If reproducing results interactively in RStudio rather than via `Rscript`,
always begin with:

```r
DIAL_SKIP_AUTORUN <- TRUE
source("foundation/Fin_sim3_clean.R")
source("foundation/Fin_Empirical5_clean.R")
```

before sourcing or running any numbered script. Failure to do this is the
single most common source of errors; see `docs/TROUBLESHOOTING.md`.
