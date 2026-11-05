# dial-cr-replication

Replication code for the paper *DIAL: Telling Invalid Instruments from
Local Effects When IV and DML Disagree*.

The scripts refer to tables and figures by their LaTeX labels (for example
`tab:sim`), so their comments stay correct when floats move. The table
below also gives the numbers in the submitted version.

## Requirements

- R 4.2 or later.
- The CRAN packages listed in `scripts/00_setup.R`. The script installs any
  that are missing.
- A multi-core CPU. No GPU is needed, and peak memory is under 4 GB. The
  run times below are for 7 cores.

## Quick start

Run every command from the `scripts/` folder.

```sh
cd scripts
Rscript 00_setup.R
Rscript 01_nhanes_build.R

# 02 and 03 define functions; these one-liners load the foundation and run them
Rscript -e 'DIAL_SKIP_AUTORUN <- TRUE; source("../foundation/Fin_sim3_clean.R"); source("../foundation/Fin_Empirical5_clean.R"); source("02_phase1_empirical.R"); run_phase1_all()'
Rscript -e 'DIAL_SKIP_AUTORUN <- TRUE; source("../foundation/Fin_sim3_clean.R"); source("../foundation/Fin_Empirical5_clean.R"); source("03_phase2_empirical.R"); run_all_phase2_applications()'

Rscript 04_synthetic_benchmark.R
Rscript 05_montecarlo_validation.R
Rscript 06_structural_simulation.R
Rscript 07_power_curve_grid.R
Rscript 08_figures.R
Rscript 09_implied_endogeneity.R
Rscript 10_rule_montecarlo.R
```

Run 05 before 10: script 10 compares the bootstrap rule with the
single-set verdicts that 05 caches. Run 07 before 08: Figure 2 reads the
grid summary that 07 writes. [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md)
gives the details, including how to run 02 and 03 in an interactive session.

## Scripts and outputs

| Script | Paper output | Files written | Time (7 cores) |
|---|---|---|---|
| `00_setup.R` | none (installs packages, smoke test) | `.setup_complete` | a few minutes on the first run |
| `01_nhanes_build.R` | NHANES input for Table 3 | `scripts/nhanes_mercury.csv` | seconds |
| `02_phase1_empirical.R` | Tables 3 and 4: Di Tella, Galiani, Banerjee, Burde | `scripts/phase1_results.csv` | about 2 hours for 02 and 03 together |
| `03_phase2_empirical.R` | Tables 3 and 4: MovieLens, Papke, NHANES, Card | `scripts/phase2_results.csv` | (see 02) |
| `04_synthetic_benchmark.R` | Table 1 (`tab:synthetic`) | `scripts/table1_results.csv` | about 45 minutes |
| `05_montecarlo_validation.R` | Table 7 (`tab:mc`) | `results/table6_results.csv`, `results/table6_timing.csv` | about 3 hours |
| `06_structural_simulation.R` | Table 9 (`tab:sim`) | `results/sim_dial_v7_results.rds`, `results/coef_vectors.rds` | about 40 minutes |
| `07_power_curve_grid.R` | data for Figure 2 (`fig:power_curves`) | `results/cp_grid_summary_lean.csv` | about 160 minutes |
| `08_figures.R` | Figures 1 to 4 | `scripts/paper/figures/fig1_*.png` to `fig4_*.png` | about a minute |
| `09_implied_endogeneity.R` | Table 5 (`tab:locality`) | `results/implied_endogeneity.csv` | about 20 minutes |
| `10_rule_montecarlo.R` | Table 8 (`tab:rule_mc`) | `results/rule_mc_results.csv` | about 8.5 hours |

Table 3 is `tab:empirical_full` and Table 4 is `tab:domain_sens`. Table 2
(`tab:denotation`) lists notation. Table 6 (`tab:frontier`) comes from a
closed-form expression; [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md)
gives the few lines of R that reproduce it. Scripts 04, 05, 06, 09 and 10
also print their table rows in LaTeX. The output of 05 keeps the table
number of an earlier layout: `table6_results.csv` holds Table 7.

## Folder layout

```
dial-cr-replication/
├── scripts/      numbered scripts 00-10; 08 writes the figures to scripts/paper/figures/
├── foundation/   shared code: CR test and MCUB sets, bootstrap rule, Jones-Pewsey DGP
├── data/         the four .dta files of the Phase 1 applications
├── results/      CSV outputs; caches are created here and ignored by git
└── docs/         run order, parameters and troubleshooting
```

## Data

- Di Tella and Schargrodsky (2013), Galiani, Rossi and Schargrodsky
  (2011), Banerjee and Iyer (2005) and Burde and Linden (2013): the `.dta`
  files in `data/`. See [data/README.md](data/README.md) for sources and
  terms.
- Papke (1995) and Card (1995): the `wooldridge` R package.
- NHANES: the `CrossScreening` R package. `01_nhanes_build.R` builds the
  analytic frame.
- MovieLens 100K: not included, under the GroupLens terms of use. Script
  03 downloads `u.data` and `u.item` from GroupLens at run time. To work
  offline, unzip `ml-100k.zip` into `scripts/ml-100k/` and call
  `run_all_phase2_applications(ml_dir = "ml-100k")`; script 09 uses that
  folder automatically when it exists.

## Reproducibility notes

- Every script fixes its seeds.
- The MCUB critical values come from simulation, so set endpoints can
  differ in the third decimal between computers, and the bootstrap share
  can differ by a few resamples.
- Bootstrap and Monte Carlo results are cached in `results/` (`boot_rule/`,
  `table6_grid_lean/`, `cp_grid_lean/`, `rule_mc/`). A rerun reuses them.
  Delete a folder to recompute it.

## Licence

The code is released under the MIT licence (see [LICENSE](LICENSE)). The
files in `data/` belong to their original authors and publishers and are
not covered by that licence; see [data/README.md](data/README.md).
