# dial-cr-replication

Replication code for the paper *DIAL: Telling Invalid Instruments from
Local Effects When IV and DML Disagree*.

The scripts refer to tables and figures by their LaTeX labels (for example
`tab:sim`), so their comments stay correct when floats move. The tables
below also give the numbers in the paper.

## Requirements

- R 4.2 or later. The reported results come from R 4.5.2 with DoubleML
  1.0.2, mlr3 1.8.0, mlr3learners 0.16.0 and ranger 0.18.0.
- The CRAN packages listed in `scripts/00_setup.R`. The script installs any
  that are missing.
- A multi-core CPU; no GPU is needed. The run times below are for 7 cores
  unless the table says otherwise.

## Quick start

Run every command from the `scripts/` folder.

```sh
cd scripts
Rscript 00_setup.R
Rscript run_all.R paper      # every script behind the paper, about 6-7 hours
```

`run_all.R` checks the setup first, runs each step in a fresh R process,
writes each step's log to `results/logs/`, and stops at the first step that
fails. It can also run part of the pipeline:

```sh
Rscript run_all.R            # steps 12, 09, 11, 06, 13 and 08; needs the outputs of 01-04
Rscript run_all.R 12 13 08   # only these steps, in the order the pipeline needs
Rscript run_all.R dry        # check the setup and list the steps; run nothing
```

In RStudio, open `run_all.R`, set the working directory to the `scripts/`
folder, and source the file. Every script also runs on its own with
`Rscript`, except 02 and 03, which define functions;
[docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) shows how to run them.

## Scripts and outputs

| Script | Paper output | Files written | Time |
|---|---|---|---|
| `00_setup.R` | none (installs packages, smoke test) | `.setup_complete` | a few minutes on the first run |
| `01_nhanes_build.R` | NHANES input | `scripts/nhanes_mercury.csv` | seconds |
| `02_phase1_empirical.R` | Table 3 (`tab:empirical_full`), MCUB sets: Di Tella, Galiani, Banerjee, Burde | `scripts/phase1_results.csv` | about 2 hours for 02 and 03 together |
| `03_phase2_empirical.R` | Table 3, MCUB sets: MovieLens, Papke, NHANES, Card | `scripts/phase2_results.csv` | (see 02) |
| `04_synthetic_benchmark.R` | Table 1 (`tab:synthetic`), MCUB sets and HSIC | `scripts/table1_results.csv` | about 45 minutes |
| `06_structural_simulation.R` | Table 8 (`tab:sim`), estimates and MCUB sets | `results/sim_dial_v7_results.rds`, `results/coef_vectors.rds` | about 40 minutes |
| `08_figures.R` | Figures 1 to 4 | `scripts/paper/figures/fig1_dial_architecture.png`, `fig2_gate_power.png`, `fig3_rstar_domain.png`, `fig4_diagnostic_space.png` | under 2 minutes |
| `09_implied_endogeneity.R` | Table 4 (`tab:locality`): r*, R, its interval and the B share | `results/implied_endogeneity.csv` | about 20 minutes |
| `11_diagnostics.R` | the Di Tella IV difference (Section 5.3), the spread of R across sample splits (Appendix B.2), the scale m (Appendix G) | `results/diagnostics.csv`, `results/ditella_late_diff.csv`, `results/dml_split_variation.csv` | about 20 minutes |
| `12_rstar_gate.R` | Table 3, gate columns (F, r*, CS_n, V, p_D); Table 6 (`tab:frontier_emp`); data for Figures 3 and 4 | `results/rstar_gate.csv` | under a minute |
| `13_rstar_gate_mc.R` | Table 7 (`tab:gate_mc`, part A); data for Figure 2 (part B); gate columns of Tables 1 and 8 (parts C and D) | `results/rstar_gate_mc_cells.csv`, `rstar_gate_power.csv`, `rstar_gate_synthetic.csv`, `rstar_gate_sim.csv` | about 40 minutes on 2 cores; about 70 minutes on Windows, which runs it on one core |

Table 2 (`tab:denotation`) lists notation. Table 5 (`tab:frontier`) comes
from a closed-form expression; [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md)
gives the R lines that reproduce it. Most scripts also print their table
rows in LaTeX.

Two helper files sit next to the numbered scripts. `dial_specs.R` defines
the 13 empirical instruments (data loaders, domains and seeds) for 09, 11
and 12, so the three scripts analyse the same rows. `run_all.R` is the
runner described above.

### Studies of the MCUB-share rule

Appendix E discusses the MCUB-share rule, which passes an instrument when
at least 95% of bootstrap MCUB sets contain zero, and notes that the
replication files contain Monte Carlo studies of it. No table or figure of
the paper uses their output, so `run_all.R` runs them only when asked
(`Rscript run_all.R 05 07 10`).

| Script | What it studies | Files written | Time |
|---|---|---|---|
| `05_montecarlo_validation.R` | single MCUB sets on the Jones-Pewsey DGP | `results/table6_results.csv`, `results/table6_timing.csv` | about 3 hours |
| `07_power_curve_grid.R` | power of a single MCUB set | `results/cp_grid_summary_lean.csv` | about 160 minutes |
| `10_rule_montecarlo.R` | the MCUB-share rule itself | `results/rule_mc_results.csv` | about 8.5 hours |

These three scripts, and the comments in 02 to 06, call this rule the
"bootstrap validity rule" or the "tau_V rule", and their output names
their own tables (`tab:mc`, `tab:rule_mc`). The DIAL column that 04 and 06
print also uses this rule; the gate verdicts of Tables 1 and 8 come from
13. The MCUB sets of 02, 03, 04 and 06 are the ones the paper reports.

## Results in this repository

`results/` holds the outputs behind the paper: the gate (12 and 13), the
implied endogeneity and the ratio R (09) and the diagnostics (11), from a
run with DoubleML on 9 November 2026, plus the outputs of 05, 07 and 10.
`scripts/` holds the MCUB sets of 02 and 03 and the table of 04, and
`scripts/paper/figures/` the four figures. With these files,
`Rscript run_all.R 08` redraws every figure in under two minutes. The
structural simulation's `results/sim_dial_v7_results.rds` is not included;
06 prints the rows of Table 8.

## Folder layout

```
dial-cr-replication/
├── scripts/      numbered scripts 00-13, dial_specs.R and run_all.R; 08 writes the figures to scripts/paper/figures/
├── foundation/   shared code: CR test and MCUB sets, MCUB-share rule, validity gate, DML fallback, Jones-Pewsey DGP
├── data/         the four .dta files of the Phase 1 applications
├── results/      outputs; caches and logs are created here and ignored by git
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
  offline, unzip `ml-100k.zip` into `scripts/ml-100k/`. `run_all.R` passes
  that folder to 03, and 09, 11 and 12 find it through `dial_specs.R`
  (or through the environment variable `DIAL_ML_DIR`).

## Reproducibility notes

- Every script fixes its seeds.
- The reported DML estimates come from DoubleML. Where DoubleML cannot be
  installed, 06, 09 and 11 use `foundation/dml_fallback.R`, which
  implements the same partialling-out estimator; its estimates differ only
  by the noise of the sample splits and forests. `run_all.R` requires
  DoubleML for these steps.
- The MCUB critical values and the gate's standard errors come from
  simulation, so set endpoints can differ in the third decimal between
  computers, and bootstrap shares by a few resamples.
- Bootstrap and Monte Carlo results are cached in `results/` (`boot_rule/`,
  `table6_grid_lean/`, `cp_grid_lean/`, `rule_mc/`). A rerun reuses them.
  Delete a folder to recompute it.

## Licence

The code is released under the MIT licence (see [LICENSE](LICENSE)). The
files in `data/` belong to their original authors and publishers and are
not covered by that licence; see [data/README.md](data/README.md).
