# data/

This folder holds the four `.dta` files of the Phase 1 applications. The
scripts read them through `dial_data_path()` in
`foundation/Fin_Empirical5_clean.R`. To keep the files elsewhere, set
`options(dial.data_dir = "path/to/folder")` before running a script.

The files belong to their original authors and publishers. They are
included so that the replication runs without downloads, and they are the
exact files behind the paper's results. The MIT licence of the code does
not cover them. Please cite the original papers when you use the data,
and follow the terms of their replication archives.

## Files

| File | Rows x columns | MD5 | Source | Read by |
|---|---|---|---|---|
| `Di_Tella_Schargrodsky_2013.dta` | 24,362 x 42 | `8cc68e1264c2387cd362e604ba027bd9` | Di Tella and Schargrodsky (2013) | 02, 09 |
| `Crime.dta` | 50,000 x 52 | `a5968b6bbc584c2c7a78c0569992e341` | Galiani, Rossi and Schargrodsky (2011) | 02, 09 |
| `yld_sett_aug03.dta` | 8,672 x 204 | `c10134075c9be0ea2074b8ae55d7c76c` | Banerjee and Iyer (2005) | 02, 09 |
| `afgan.dta` | 687 x 17 | `378bca1775ab2baf875d2994bdf56f10` | Burde and Linden (2013) | 02, 09 |

`Di_Tella_Schargrodsky_2013.dta` is the archive file
`JPE - Di Tella and Schargrodsky - CriminalRecidivismAfterPrisonAndElectronicMonitoring .dta`
under a shorter name. `afgan.dta` is an extract of the Burde and Linden
data with 17 variables.

## Variables used

| Application | Outcome Y | Treatment D | Instrument Z | Sample |
|---|---|---|---|---|
| Di Tella | `recidivism` | `electronicMonitoring` | `judgeAlreadyUsedEM`, `percJudgeSentToEM` | `offendersPerJudge > 9` |
| Galiani | `crimerate` | `sm` | `highnumber` | cohorts 1958 to 1962 |
| Banerjee | `phwht` (share of wheat area under high-yielding varieties) | `p_nland` (share of non-landlord area) | `instru` (conquered between 1820 and 1856) | `year >= 1965` and `phwht <= 1` |
| Burde | `testscore` | `enrolled` | `buildschool` | all rows |

The controls of each application are listed in
`scripts/02_phase1_empirical.R`.

## Other data (not stored here)

- Papke (1995): `wooldridge::k401k`.
- Card (1995): `wooldridge::card`.
- NHANES: `CrossScreening::nhanes.fish`, the 1,107-row extract documented
  by Rosenbaum (2014) and Zhao, Small and Rosenbaum (2018).
  `scripts/01_nhanes_build.R` writes the analytic frame to
  `scripts/nhanes_mercury.csv`.
- MovieLens 100K (Harper and Konstan 2015): script 03 downloads it from
  GroupLens at run time. It is not redistributed here, under the GroupLens
  terms of use (https://grouplens.org/datasets/movielens/100k/).

## References

Banerjee, A. and Iyer, L. (2005). History, institutions, and economic
performance: The legacy of colonial land tenure systems in India.
*American Economic Review*, 95(4), 1190-1213.

Burde, D. and Linden, L. L. (2013). Bringing education to Afghan girls: A
randomized controlled trial of village-based schools. *American Economic
Journal: Applied Economics*, 5(3), 27-40.

Card, D. (1995). Using geographic variation in college proximity to
estimate the return to schooling. In L. N. Christofides, E. K. Grant and
R. Swidinsky (Eds.), *Aspects of Labour Market Behaviour: Essays in Honour
of John Vanderkamp*, pp. 201-222. University of Toronto Press.

Di Tella, R. and Schargrodsky, E. (2013). Criminal recidivism after prison
and electronic monitoring. *Journal of Political Economy*, 121(1), 28-73.

Galiani, S., Rossi, M. A. and Schargrodsky, E. (2011). Conscription and
crime: Evidence from the Argentine draft lottery. *American Economic
Journal: Applied Economics*, 3(2), 119-136.

Harper, F. M. and Konstan, J. A. (2015). The MovieLens datasets: History
and context. *ACM Transactions on Interactive Intelligent Systems*, 5(4),
1-19.

Papke, L. E. (1995). Participation in and contributions to 401(k) pension
plans: Evidence from plan data. *Journal of Human Resources*, 30(2),
311-325.

Rosenbaum, P. R. (2014). Weighted M-statistics with superior design
sensitivity in matched observational studies with multiple controls.
*Journal of the American Statistical Association*, 109(507), 1145-1158.

Zhao, Q., Small, D. S. and Rosenbaum, P. R. (2018). Cross-screening in
observational studies that test many hypotheses. *Journal of the American
Statistical Association*, 113(523), 1070-1084.
