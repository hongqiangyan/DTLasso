# Debiased Threshold Lasso (DTLasso) Estimator

This repository contains the R codes used to run the simulations and applications in the paper:

**Li, Jiatong and Yan, Hongqiang (2025).**  
*Uniform Inference in High-Dimensional Threshold Regression Models.*  
arXiv preprint [arXiv:2404.08105](https://arxiv.org/abs/2404.08105).

---

## Repository Structure

### [1] **Simulations**
Files:
- `myfnc.R`
- `SimulationSyn.R`

Usage:
- Run `SimulationSyn.R` to generate the simulation tables and figures reported in **Section 5.1** of the paper.  
- Ensure that `myfnc.R` is located in the same directory, as it provides supporting functions.

---

### [2] **GrowthRegressions**  
Application to cross-country growth regressions (**Section 5.2**).

Files:
- `GrowthRegressions_data.cvs`  
- `GrowthRegressions_data_lr.cvs`  
- `growth_gdp_DTLasso.R`  
- `growth_lr_DTLasso.R`  

Usage:
- `GrowthRegressions_data.cvs`: dataset for **Table 3** (threshold variable: `gdp60`).  
- `GrowthRegressions_data_lr.cvs`: dataset for **Table B.2** (threshold variable: `lr`).  
- Variable definitions are listed in **TABLE B.1** of the paper.  
- `growth_gdp_DTLasso.R`: generates **Table 3**.  
- `growth_lr_DTLasso.R`: generates **Table B.2**.  
- Requires `myfnc.R`.

---

### [3] **LocalProjectionGovSpending**  
Application to fiscal multipliers (**Section 5.3**).

Files:
- `processed_data.xls`  
- `HDLPTreasuryBillRate.R`  
- `HDLPUnemploymentRate.R`  

Usage:
- `processed_data.xls`: dataset from Valerie A. Ramey’s website, used in **Ramey & Zubairy (2018)**.  
- `HDLPTreasuryBillRate.R`: generates **Figure 3** (Treasury Bill Rate responses).  
- `HDLPUnemploymentRate.R`: generates **Figure 4** (Unemployment Rate responses).  

Data are pre-processed as described in the paper.

---

## References
- Li, J., and Yan, H. (2025). *Uniform inference in high-dimensional threshold regression models.* arXiv preprint arXiv:2404.08105.  
- Ramey, V. A., & Zubairy, S. (2018). Government spending multipliers in good times and in bad: evidence from US historical data. *Journal of Political Economy, 126*(2), 850–901.  
