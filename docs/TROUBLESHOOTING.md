# Troubleshooting

Common errors, their causes and fixes.

## "required functions not loaded" from 02 or 03

The foundation files were not sourced first. Scripts 02 and 03 define
functions only. Start a fresh R session in `scripts/` and run:

```r
DIAL_SKIP_AUTORUN <- TRUE
source("../foundation/Fin_sim3_clean.R")
source("../foundation/Fin_Empirical5_clean.R")
source("02_phase1_empirical.R")      # or 03_phase2_empirical.R
```

Then call `run_phase1_all()` or `run_all_phase2_applications()`.

## Long runs start as soon as a foundation file is sourced

`DIAL_SKIP_AUTORUN <- TRUE` was not set before the first `source()` call.
Setting it afterwards does not stop a run that has started. Restart R and
set the flag first.

## "Data file not found"

`dial_data_path()` could not find a `.dta` file. It looks in
`getOption("dial.data_dir")`, then in `data/` at the repository root. Run
the scripts from `scripts/`, or set
`options(dial.data_dir = "path/to/data")`.

## "nhanes_mercury.csv not found"

Run `Rscript 01_nhanes_build.R` from `scripts/` first. Scripts 03 and 09
read the file from the working folder.

## MovieLens download fails

Script 03 reads `u.data` and `u.item` from GroupLens when no local copy is
given. Without internet access, download `ml-100k.zip` from
https://grouplens.org/datasets/movielens/100k/, unzip it into
`scripts/ml-100k/`, and call
`run_all_phase2_applications(ml_dir = "ml-100k")`. Script 09 uses
`scripts/ml-100k/` automatically when the folder exists.

## Parallel workers fail or return NA

The scripts use `future` multisession workers (up to 7, one fewer than the
physical cores). If workers run out of memory, lower `N_WORKERS` at the top
of the script (05, 07, 09, 10). The bootstrap rule takes its worker count
from `dial_boot_plan()`: in 02 and 03, call `dial_boot_plan(2)` after
sourcing the script and before running it; in 04 and 06, change the
script's `dial_boot_plan()` call to `dial_boot_plan(2)`. Completed cells
and resamples are cached, so a rerun resumes where the failed run stopped.

## Results differ slightly from the paper

The MCUB critical values come from simulation, so set endpoints can move
in the third decimal between computers, and the bootstrap share p_hat_B
by a few resamples. Stale caches can also cause differences after you
edit the code: delete the relevant folder in `results/` (see
[REPRODUCIBILITY.md](REPRODUCIBILITY.md), section 4) and rerun.

## Figure 2 is missing

`08_figures.R` draws Figure 2 from `results/cp_grid_summary_lean.csv`,
which `07_power_curve_grid.R` writes. Run 07 first, or keep the CSV that
the repository ships.
