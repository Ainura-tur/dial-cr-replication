# dial-cr-replication

Replication code for **"DIAL: A Probabilistic Diagnostic for Instrument Validity and Effect Heterogeneity in Causal Machine Learning"** (NeurIPS 2026).

---

## Quick start

```r
# 1. Clone and open a terminal in the repo root, then start R there.
# 2. Run in order:

Rscript scripts/00_setup.R               # installs packages, smoke-tests foundation
Rscript scripts/02_phase1_empirical.R    # Tables 4-5: Phase 1 empirical re-runs
Rscript scripts/08_figures.R             # Figures 1-3
```

For the full pipeline including Monte Carlo and power curves, see
[docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

---

## Source order requirement

**This is the most common source of errors.** The foundation files contain
long-running auto-run blocks that must be suppressed before sourcing. Every
numbered script sets this automatically, but if you source files manually in
an interactive session, run this first:

```r
DIAL_SKIP_AUTORUN <- TRUE
source("foundation/Fin_sim3_clean.R")
source("foundation/Fin_Empirical5_clean.R")
```

Sourcing in the wrong order or without the flag will trigger unintended
Monte Carlo runs and produce the "required functions not loaded" error in
downstream scripts.

---

## Repository structure

```
dial-cr-replication/
├── scripts/
│   ├── 00_setup.R                  install dependencies; smoke test
│   ├── 01_build_nhanes.R           build nhanes_mercury.csv from CrossScreening
│   ├── 02_phase1_empirical.R       Tables 4-5: Di Tella, Galiani, Burde, Banerjee
│   ├── 03_phase2_empirical.R       Tables 4-5: Papke, NHANES, Card, MovieLens
│   ├── 04_synthetic_table1.R       Table 1: MCUB identified sets (synthetic)
│   ├── 05_montecarlo_table6.R      Table 6: Monte Carlo coverage results
│   ├── 06_structural_table7.R      Table 2: structural simulation (DIAL scenarios)
│   ├── 07_power_curve_grid.R       Figure 2 (power curves): parallelised MC grid
│   ├── 08_figures.R                Figures 1, 3, 4: architecture and empirical plots
│   ├── cr_wrapper.R                unified CR test entry point (all scenarios)
│   └── jp_dgp.R                    Jones-Pewsey data-generating process
│
├── foundation/                     upstream CR test code (see foundation/README.md)
│   ├── Fin_sim3_clean.R
│   └── Fin_Empirical5_clean.R
│
├── data/                           see data/README.md for access instructions
│   ├── ditella2013/
│   ├── galiani2011/
│   ├── banerjee2007/
│   ├── burde2013/
│   ├── papke1995/                  loaded via wooldridge::k401k
│   ├── card1995/                   loaded via wooldridge::card
│   ├── nhanes/                     loaded via CrossScreening::nhanes.fish
│   └── movielens100k/              download instructions only; not re-hosted
│
├── paper/
│   └── figures/                    compiled PNG figures for the paper
│
├── results/                        gitignored; recreated by scripts
│
└── docs/
    ├── REPRODUCIBILITY.md          step-by-step run order with timing
    ├── PARAMETERS.md               all hardcoded constants and their provenance
    └── TROUBLESHOOTING.md          common errors and fixes
```

---

## Script-to-output map

| Script | Output |
|--------|--------|
| `01_build_nhanes.R` | `nhanes_mercury.csv` |
| `02_phase1_empirical.R` | `results/phase1_results.csv`, Tables 4-5 rows 1-4 |
| `03_phase2_empirical.R` | Tables 4-5 rows 5-8 |
| `04_synthetic_table1.R` | `results/table1_results.csv`, Table 1 |
| `05_montecarlo_table6.R` | `results/table6_results.csv`, Table 6 |
| `06_structural_table7.R` | Table 2, `paper/figures/sim_dial_scenarios.png` |
| `07_power_curve_grid.R` | `results/cp_grid/`, Figure 2 power curves |
| `08_figures.R` | `paper/figures/fig0_dial_architecture.png`, `fig1_diagnostic_space.png`, `fig3_domain_sensitivity.png` |

---

## Runtime estimates (7-core machine)

| Script | Approximate wall time |
|--------|-----------------------|
| `00_setup.R` | 5-10 minutes (first run; package installs) |
| `02_phase1_empirical.R` | 20-40 minutes |
| `03_phase2_empirical.R` | 20-40 minutes |
| `04_synthetic_table1.R` | 5 minutes |
| `05_montecarlo_table6.R` | 4-8 hours |
| `06_structural_table7.R` | 3 minutes per row |
| `07_power_curve_grid.R` | 75 minutes (parallelised) |
| `08_figures.R` | 2 minutes |

---

## Software requirements

- R >= 4.2.0
- CRAN packages: see `scripts/00_setup.R` for the full list
- GitHub package: `ivcrtest` (installed automatically by `00_setup.R`)
- Package versions are recorded in `renv.lock`

---

## Data access

All datasets are publicly available. Three are loaded directly via R packages
(`wooldridge`, `CrossScreening`). The remainder are downloaded from the
original supplement mirrors listed in `data/README.md`. MovieLens 100K is
downloaded at runtime from the GroupLens server; do not re-host it.

---

## Anonymisation note

This repository is submitted anonymously for review. After acceptance, the
permanent URL and author information will be added here.
