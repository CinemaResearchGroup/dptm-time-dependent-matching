# Dynamic Propensity Trajectory Matching (DPTM): Simulation Code

Simulation code for the paper *Dynamic Propensity Trajectory Modeling and Matching with Time-Dependent Covariates for Estimating Treatment Effect* (Statistics in Medicine, 2026).

## Files

- **`true_value_function.R`**: `true_value()` approximates the true log hazard ratio among the treated. It simulates a large population, generates a counterfactual never-treated outcome for each treated subject, and fits a frailty Cox model to the paired outcomes. This value is the reference for bias and coverage.

- **`naive_function.R`**: `estimated_naive_value()` runs the naive analysis. Treatment is treated as a fixed binary exposure, the propensity score comes from logistic regression, subjects are grouped by full matching, and the effect is estimated with a Cox model with robust variance.

- **`pure_control_function.R`**: `estimated_purecontrol_value()` runs DPT-based pure control matching (PCM). Each treated subject is matched to a subject who is never treated during follow-up, and the effect is estimated with Cox and frailty Cox models.

- **`all_risk_function.R`**: `estimated_allrisk_value()` runs DPT-based all risk-set matching (ARM). Controls can include subjects who start treatment later; they are censored when they switch. The effect is estimated with frailty Cox and IPCW-adjusted frailty Cox models. The helper `ipcw.new()` at the top of the file computes the IPCW weights and must be run first.
