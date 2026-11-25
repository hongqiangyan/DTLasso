# Debiased Threshold Lasso (DTLasso) Estimator

This repository contains the R codes and data used to run the simulations and applications in the paper:

**Li, Jiatong and Yan, Hongqiang (2026).**
*Inference for High-Dimensional Regression with a Possible Threshold.*
arXiv preprint [arXiv:2404.08105v4](https://arxiv.org/abs/2404.08105v4).

---

## Repository Structure

### [1] **simulation**
Monte Carlo simulation (Section 4.1).

Files:
- `code/1_simulate.R`
- `code/2_make_tables.R`
- `code/myfnc.R`

Usage:
- Run `1_simulate.R` to simulate the seven designs, with 500 replications each. The results are saved in `results/`. One design takes about 30–60 hours on 16 cores.
- Run `2_make_tables.R` to generate Tables 1 and 2 from `results/`.
- `myfnc.R` provides the supporting functions (the estimators).
- The simulations in the paper were run on Arizona State University's Sol cluster (Linux, R 4.5.1). `MASS::mvrnorm()` uses `eigen()`, whose output depends on the platform's LAPACK library, so a rerun on another platform uses the same seeds but is not bit-identical.

### [2] **cmr**
Application to racial tipping (Section 4.2 and the Online Supplement).

Files:
- `data/analysis_1990_2010_tract2000.zip`
- `data/msa_def.dta`
- `data/msa_coordinates.csv`
- `code/1_estimate.R`
- `code/2_make_tables.R`
- `code/3_make_figures.R`
- `code/preflight_card_cities_gt.R`, `code/tipping_DTLasso_percity_anchor_gt.R`, `code/summarize_results_gt.R`
- `code/myfnc_anchor_gt.R`, `code/common_lambda_cv_gt.R`, `code/sample_search_rules_gt.R`

Usage:
- `analysis_1990_2010_tract2000.zip`: census-tract data for 1990–2000 and 2000–2010 on 2000 tract boundaries, built from IPUMS NHGIS. The outcome is `outcome_white_change`, the threshold variable is `forcing_minority_share`, and there are 50 candidate controls.
- `msa_def.dta`: metropolitan-area names. `msa_coordinates.csv`: map locations of the metropolitan areas.
- `1_estimate.R`: estimates all 197 metropolitan area–decade samples, which takes about 4 hours on 4 cores. It runs `preflight_card_cities_gt.R` (selects the samples), `tipping_DTLasso_percity_anchor_gt.R` (estimates each sample) and `summarize_results_gt.R` (collects the results).
- `2_make_tables.R`: generates Tables 3 and 4, and the coefficient tables for Indianapolis, Nashville and Chicago in the Online Supplement.
- `3_make_figures.R`: generates the two maps of tipping points in the Online Supplement.
- `myfnc_anchor_gt.R`, `common_lambda_cv_gt.R` and `sample_search_rules_gt.R` provide the supporting functions.

### [3] **rz**
Application to fiscal multipliers (Online Supplement).

Files:
- `data/data.rds`
- `code/1_threshold.R`
- `code/2_local_projections.R`
- `code/3_make_figures.R`
- `code/4_make_tables.R`
- `code/config.R`, `code/bootstrap_seeds.csv`
- `code/common.R`, `code/design.R`, `code/lp.R`, `code/multiplier_inference.R`, `code/threshold_node_recovery.R`, `code/threshold_plugin.R`, `code/threshold_plugin_anchor.R`, `code/threshold_test.R`

Usage:
- `data.rds`: quarterly data for 1890Q1–2015Q4. It contains GDP, government spending, military spending news, the unemployment rate and 14 control series. It is built from Ramey and Zubairy (2018) and the other sources listed below.
- `1_threshold.R`: estimates the unemployment threshold and runs the bootstrap test for a threshold effect.
- `2_local_projections.R`: estimates the state-dependent local projections and the cumulative multipliers.
- `3_make_figures.R`: generates the two impulse-response figures.
- `4_make_tables.R`: generates the table of cumulative fiscal multipliers.
- `config.R` holds the tuning constants and seeds, and `bootstrap_seeds.csv` the bootstrap seeds. The other files in `code/` provide the supporting functions.

---

## Software

R 4.5.1 with the packages `glmnet` and `Matrix`. In addition:
- `simulation`: `MASS`, `expm`, `flare`.
- `cmr`: `haven`, `sf`, `maps`, `ragg`.
- `rz`: `ggplot2`, `patchwork`, `ragg`.

## Data

The data are redistributed only to reproduce the results of the paper. The NHGIS terms permit publishing a subset of NHGIS data to meet journal requirements. The series in `rz/data/data.rds` that are derived from the Jordà–Schularick–Taylor Macrohistory Database are shared under its CC BY-NC-SA 4.0 licence. The BAA yield is © Moody's, retrieved from FRED. Please cite the original data sources listed below.

## References

- Li, J., and Yan, H. (2026). Inference for high-dimensional regression with a possible threshold. arXiv preprint arXiv:2404.08105v4.
- Schroeder, J., Van Riper, D., Manson, S., Cooper, G., Krause, Z., Kugler, T., Zhu, T., and Ruggles, S. (2026). IPUMS National Historical Geographic Information System: Version 21.0 [dataset]. Minneapolis, MN: IPUMS. http://doi.org/10.18128/D050.V21.0
- Card, D., Mas, A., and Rothstein, J. (2008). Tipping and the dynamics of segregation. Quarterly Journal of Economics, 123(1), 177–218.
- Ramey, V. A., and Zubairy, S. (2018). Government spending multipliers in good times and in bad: evidence from US historical data. Journal of Political Economy, 126(2), 850–901.
- Bernardini, M., and Peersman, G. (2018). Private debt overhang and the government spending multiplier: evidence for the United States. Journal of Applied Econometrics, 33(4), 485–508.
- Jordà, Ò., Schularick, M., and Taylor, A. M. (2017). Macrofinancial history and the new business cycle facts. NBER Macroeconomics Annual 2016, 31. Data: Jordà–Schularick–Taylor Macrohistory Database, Release 6.
- Shiller, R. J. (2015). Irrational Exuberance, 3rd edition. Princeton University Press. Data: https://shillerdata.com/
- Balke, N. S., and Gordon, R. J. (1986). Appendix B: Historical data. In R. J. Gordon (ed.), The American Business Cycle: Continuity and Change. University of Chicago Press.
- Barro, R. J., and Redlick, C. J. (2011). Macroeconomic effects from government purchases and taxes. Quarterly Journal of Economics, 126(1), 51–102.
- Mertens, K. (2013). Marginal tax rates and income: new time series evidence. NBER Working Paper 19171.
- Piketty, T., and Saez, E. (2003). Income inequality in the United States, 1913–1998. Quarterly Journal of Economics, 118(1), 1–41.
- Goemans, P. (2022). Historical evidence for larger government spending multipliers in uncertain times than in slumps. Economic Inquiry, 60(3), 1164–1185.
- Federal Reserve Bank of St. Louis. FRED, Federal Reserve Economic Data. Moody's Seasoned Baa Corporate Bond Yield [BAA] and other series. https://fred.stlouisfed.org/
