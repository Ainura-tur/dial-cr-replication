# data/

This directory contains per-dataset subdirectories. Datasets loaded via R
packages do not require manual downloading. Datasets loaded from external
`.dta` files are downloaded at runtime by the analysis scripts from the
mirror URLs listed below.

---

## Di Tella and Schargrodsky (2013)

**Study:** Electronic monitoring and criminal recidivism in Argentina.  
**Script:** `scripts/02_phase1_empirical.R` via `run_ditella_analysis_fixed()`  
**Variables:** `recidivism`, `electronicMonitoring`, `judgeAlreadyUsedEM`,
`percJudgeSentToEM`, and controls.  
**Access:** Downloaded at runtime via `haven::read_dta()` from:

  https://github.com/ratbekd/Orientation_paper/raw/refs/heads/main/JPE%20-%20Di%20Tella%20and%20Schargrodsky%20-%20CriminalRecidivismAfterPrisonAndElectronicMonitoring%20.dta

**Citation:** Di Tella, R. and Schargrodsky, E. (2013). Criminal recidivism
after prison and electronic monitoring. *Journal of Political Economy*,
121(1), 28-73.

---

## Burde and Linden (2013)

**Study:** School construction and test scores in Afghanistan.  
**Script:** `scripts/02_phase1_empirical.R` via `run_burde_analysis_fixed()`  
**Variables:** `testscore`, `enrolled`, `buildschool`, and controls.  
**Access:** Downloaded at runtime from:

  https://github.com/ratbekd/Orientation_paper/raw/refs/heads/main/afgan.dta

**Citation:** Burde, D. and Linden, L. L. (2013). Bringing education to
Afghan girls: A randomized controlled trial of village-based schools.
*American Economic Journal: Applied Economics*, 5(3), 27-40.

---

## Galiani, Rossi, and Schargrodsky (2011)

**Study:** Conscription and crime in Argentina.  
**Script:** `scripts/02_phase1_empirical.R` via `run_galiani_analysis_fixed()`  
**Variables:** `crimerate`, `sm`, `highnumber`, and controls.  
**Access:** Downloaded at runtime from:

  https://github.com/ratbekd/Orientation_paper/raw/refs/heads/main/Crime.dta

**Citation:** Galiani, S., Rossi, M. A., and Schargrodsky, E. (2011).
Conscription and crime: Evidence from the Argentine draft lottery.
*American Economic Journal: Applied Economics*, 3(2), 119-136.

---

## Banerjee, Mookherjee, Munshi, and Ray (2007)

**Study:** Land reform and agricultural productivity in India.  
**Script:** `scripts/02_phase1_empirical.R` via `run_banerjee_analysis_fixed()`  
**Variables:** `phwht`, `p_nland`, `instru`, and controls.  
**Access:** Downloaded at runtime from:

  https://github.com/ratbekd/Orientation_paper/raw/refs/heads/main/yld_sett_aug03.dta

**Citation:** Banerjee, A., Mookherjee, D., Munshi, K., and Ray, D. (2001).
Inequality, control rights, and rent seeking: Sugar cooperatives in
Maharashtra. *Journal of Political Economy*, 109(1), 138-190.

---

## Papke (1995) -- 401(k) participation

**Study:** Employer match rates and 401(k) participation.  
**Script:** `scripts/03_phase2_empirical.R` via `run_papke_analysis_fixed()`  
**Variables:** `prate`, `mrate`, `sole`, `age`, `ltotemp`.  
**Access:** Loaded directly via the `wooldridge` R package:

```r
install.packages("wooldridge")
data <- wooldridge::k401k
```

**Citation:** Papke, L. E. (1995). Participation in and contributions to
401(k) pension plans: Evidence from plan data. *Journal of Human Resources*,
30(2), 311-325.

---

## Card (1995) -- returns to education

**Study:** Returns to schooling using father's education as an instrument.  
**Script:** `scripts/03_phase2_empirical.R` via `run_card_analysis_fixed()`  
**Variables:** `lwage`, `educ`, `fatheduc`, and controls.  
**Access:** Loaded directly via the `wooldridge` R package:

```r
install.packages("wooldridge")
data <- wooldridge::card
```

**Citation:** Card, D. (1995). Using geographic variation in college
proximity to estimate the return to schooling. In L. N. Christofides,
E. K. Grant, and R. Swidinsky (Eds.), *Aspects of Labour Market Behaviour:
Essays in Honour of John Vanderkamp*, pp. 201-222. University of Toronto Press.

---

## NHANES -- fish consumption and mercury exposure

**Study:** Fish consumption frequency and blood mercury levels.  
**Script:** `scripts/01_build_nhanes.R` builds `nhanes_mercury.csv`;
`scripts/03_phase2_empirical.R` consumes it.  
**Variables:** `log_mercury`, `fish_freq`, `income`, `education`, `age`,
`female`, `bmi`.  
**Access:** Loaded via the `CrossScreening` R package (Zhao et al.):

```r
install.packages("CrossScreening")
library(CrossScreening)
data(nhanes.fish)
```

Run `scripts/01_build_nhanes.R` to build the analysis-ready CSV before
running `scripts/03_phase2_empirical.R`.

**Citation:** Rosenbaum, P. R. (2014). Weighted M-statistics with superior
design sensitivity in matched observational studies with multiple controls.
*Journal of the American Statistical Association*, 109(507), 1145-1158.

---

## MovieLens 100K

**Study:** Recommendation effects; movie features as instruments.  
**Script:** `scripts/03_phase2_empirical.R` via `run_movielens_analysis_fixed()`  
**Variables:** `mean_rating`, `log_n`, `num_genres`, `release_year`, `is_sequel`.  
**Access:** Downloaded at runtime from the GroupLens server:

  https://files.grouplens.org/datasets/movielens/ml-100k.zip

This dataset is not re-hosted here. See the GroupLens terms of use at:

  https://grouplens.org/datasets/movielens/100k/

**Citation:** Harper, F. M. and Konstan, J. A. (2015). The MovieLens
datasets: History and context. *ACM Transactions on Interactive Intelligent
Systems*, 5(4), 1-19.
