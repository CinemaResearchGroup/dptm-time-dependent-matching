# Dynamic Propensity Trajectory Matching (DPTM): Simulation Code

Simulation code for the paper:

> Wu X\*, Yang Y\*, Xu Y\*, Li L, Zhang L, Rahman M, McCoy RG, Chinchilli VM, Chen C, Wang M. **Dynamic Propensity Trajectory Modeling and Matching with Time-Dependent Covariates for Estimating Treatment Effect.** *Statistics in Medicine* (2026). (\*co-first authors)

The paper introduces the dynamic propensity trajectory (DPT) framework and DPT-based matching (DPTM) for observational studies in which treatment is initiated at different times, confounders change over time, and the outcome is time to an event. Rather than matching only on the propensity score at the moment of treatment initiation, DPTM matches on the entire propensity trajectory up to that moment.

This repository reproduces the simulation studies in Section 3 of the paper. The simulations compare DPTM under two matching schemes, pure control matching (PCM) and all risk-set matching (ARM), against a naive analysis. The target is the treatment effect on the treated, measured as a log hazard ratio for the post-treatment event time.

## Repository contents

| File | Function | What it does | Results in the paper |
|---|---|---|---|
| `true_value_function.R` | `true_value()` | Approximates the true log hazard ratio from a large simulated population in which each treated subject has an observed and a counterfactual (never-treated) outcome | Reference value for bias and coverage |
| `naive_function.R` | `estimated_naive_value()` | Naive analysis: treatment as a fixed binary exposure, logistic-regression propensity score, full matching | Table 1, "Naive Model" |
| `pure_control_function.R` | `estimated_purecontrol_value()` | DPT-based pure control matching, then Cox and frailty Cox models | Table 1, "After Matching"; Table 2, "Pure Control Matching" |
| `all_risk_function.R` | `estimated_allrisk_value()` | DPT-based all risk-set matching, then frailty Cox and IPCW-adjusted frailty Cox models | Table 2, "All Risk Set Matching" |

`pure_control_function.R` and `all_risk_function.R` also define a helper, `ipcw.new()`, a modified version of `ipcwswitch::ipcw()` that computes the inverse probability of censoring weights used in the all risk-set analysis.

### Arguments

| Argument | Used by | Meaning |
|---|---|---|
| `lambda_z` | all four functions | Baseline hazard of treatment initiation, $\lambda^T_0$. Controls the proportion of subjects treated (see Scenarios). |
| `sample_size` | all four functions | Number of subjects $n$ per simulated dataset. For `true_value()`, the size of the large reference population. |
| `simulation` | the three estimation functions | Number of Monte Carlo replicates. |
| `true_phi` | the three estimation functions | True log hazard ratio from `true_value()`, used to compute bias, MSE and coverage. |

### Return values

`true_value()` prints and returns a data frame with `true` (the true log hazard ratio) and `ratio` (the proportion of subjects treated).

The three estimation functions return a list of per-replicate results:

| Function | Estimates | Standard errors | 95% CI coverage indicators | Squared errors | Other |
|---|---|---|---|---|---|
| `estimated_naive_value()` | `naive_est` | `naive_se` | `naive_ci` | `naive_mse` | |
| `estimated_purecontrol_value()` | `frailty_est`, `cox_est` | `frailty_se`, `cox_se` | `frailty_ci`, `cox_ci` | `frailty_mse`, `cox_mse` | `true.treated_after`, `treated_before`, `treated.to.control` |
| `estimated_allrisk_value()` | `frailty_est`, `ipcw_est` | `frailty_se`, `ipcw_se` | `frailty_ci`, `ipcw_ci` | `frailty_mse`, `ipcw_mse` | `true.treated_after`, `treated_before`, `treated.to.control` |

The extra elements record, for each replicate, the number of treated subjects after matching, the number treated before matching, and the number of treated subjects who were used as controls. Result vectors are pre-allocated to length 1000, so unused entries and replicates that fail are `NA`; summarize them with `na.omit()`.

## Simulation design

### Covariates

Each subject $i = 1, \dots, n$ has two time-invariant covariates, $X_{i1} \sim N(1, 1)$ and $X_{i2} \sim \mathrm{Bernoulli}(0.6)$, and one time-varying covariate $X_{i3}$. The time-varying covariate is recorded at $t \in \{0, 0.5\}$ from a linear mixed-effects model with a subject-specific random intercept $b_{1i}$ and random slope $b_{2i}$,

$$
X_{i3}(t) = \alpha_0 + \alpha_1 t + b_{1i} + b_{2i} t + \epsilon_i(t), \qquad t \in \{0, 0.5\},
$$

and is held constant between measurements, so $X_{i3}(v) = X_{i3}(0)$ for $v < 0.5$ and $X_{i3}(v) = X_{i3}(0.5)$ for $v \ge 0.5$. Here $\alpha_0 = \alpha_1 = 2$, the errors $\epsilon_i(t) \sim N(0, 0.05)$ are drawn independently at each measurement, and $(b_{1i}, b_{2i})^\top$ is bivariate normal with mean zero, $\mathrm{var}(b_{1i}) = 0.1$, $\mathrm{var}(b_{2i}) = 0.025$ and $\mathrm{cov}(b_{1i}, b_{2i}) = 0.05$.

### Time to treatment initiation

Treatment initiation times $T_i$ follow a Cox model with the time-varying covariate:

$$
\lambda^T_i\{v \mid X_i(v)\} = \lambda^T_0 \exp\{\beta_{1T} X_{i1} + \beta_{2T} X_{i2} + \beta_{3T} X_{i3}(v)\},
$$

with $(\beta_{1T}, \beta_{2T}, \beta_{3T}) = (0.3, 0.3, 0.2)$ and $\lambda^T_0$ set by `lambda_z`.

### Time to the clinical event

Event times $E_i$ depend on the covariates and on the current treatment status $A_i(v) = I(v \ge T_i)$:

$$
\lambda^E_i\{v \mid X_i(v), T_i\} = \lambda^E_0 \exp\{\beta_{1E} X_{i1} + \beta_{2E} X_{i2} + \beta_{3E} X_{i3}(v) + \theta A_i(v)\},
$$

with $\lambda^E_0 = 1.5$, $(\beta_{1E}, \beta_{2E}, \beta_{3E}) = (0.1, 0.1, 0.1)$ and $\theta = -0.5$. The counterfactual never-treated hazard drops the $\theta A_i(v)$ term.

Both hazards are piecewise constant in $v$: they change only at $v = 0.5$ and, for the event hazard, at $T_i$. Treatment and event times are therefore drawn exactly by inverting the piecewise-linear cumulative hazard at a uniform random draw.

### Censoring and observed data

Censoring times are non-informative, $C_i \sim U(0.3, 1)$. The observed data and the corresponding variable names in the code are:

| Quantity | Definition | Variable in code |
|---|---|---|
| Treatment initiation time | $T_i$ | `T_z` |
| Event time | $E_i$ | `T_e` |
| Censoring time | $C_i$ | `CensorT` |
| Observed event or censoring time | $U_i = \min(E_i, C_i)$ | `Y_e` |
| Event indicator | $\Delta^E_i = I(E_i \le C_i)$ | `observed.event` |
| Observed time to treatment | $Y_i = \min(T_i, U_i)$ | `Y_z` |
| Treatment indicator | $\Delta^T_i = I(T_i \le U_i)$ | `observed.treat` |
| Post-treatment event time | time from the matched treatment initiation time to event or censoring | `Y_pe` |

### Scenarios

Each treatment ratio is combined with sample sizes $n \in \{500, 1000, 1500\}$ and 1,000 Monte Carlo replicates. The true log hazard ratios below were computed with `true_value()`.

| `lambda_z` ($\lambda^T_0$) | Approximate treatment ratio | True log HR (`true_phi`) | Reported in |
|---|---|---|---|
| 0.30 | 20% | -0.5163079 | Table 1 |
| 0.45 | 30% | -0.4967315 | Table 1 |
| 0.70 | 40% | -0.4900718 | Tables 1 and 2 |
| 1.10 | 50% | -0.4937465 | Table 2 |
| 1.55 | 60% | -0.5002375 | Table 2 |

## Methods

### Dynamic propensity trajectory and matching (PCM and ARM)

A time-dependent Cox model for treatment initiation is fit to the observed data in counting-process format (follow-up split at $v = 0.5$), giving coefficient estimates $\hat\beta$. Each subject's time-dependent propensity score is the estimated relative hazard of treatment, $\psi_i(v) = \exp\{\hat\beta^\top X_i(v)\}$.

A treated subject $i$ who initiates at $T_i$ and a candidate control $j$ are compared through the cumulative squared distance between their propensity trajectories,

$$
M_{i,j} = \int_0^{T_i} \{\psi_i(v) - \psi_j(v)\}^2 \, dv,
$$

which is the paper's general distance with $f(z) = z^2$ and equal weights $\omega(v) = 1$. Because the covariates are piecewise constant, the integral is computed exactly as a sum over the intervals $[0, 0.5)$ and $[0.5, T_i]$, each weighted by its length.

Matching is sequential and without replacement. Treatment initiation times are processed in increasing order, and at each time the treated subject is paired (with `optmatch::pairmatch()`) to the candidate with the smallest $M_{i,j}$. Both subjects are then removed from further matching.

### Pure control matching (`estimated_purecontrol_value`)

Under PCM, the candidates for a treated subject initiating at $T_i$ are the not-yet-matched subjects who are still under observation after $T_i$ and are never treated during their observed follow-up (`observed.treat == 0` and `Y_z > T_i`). In each matched pair, both subjects are followed from $T_i$ to their event or censoring time. The log hazard ratio is estimated with a frailty Cox model with a matched-pair random effect (`coxme`) and with a standard Cox model (`coxph`).

PCM chooses controls using future information (that they will never be treated), so it is expected to be biased, and the simulations show the bias growing with the treatment ratio.

### All risk-set matching (`estimated_allrisk_value`)

Under ARM, the candidates are all not-yet-matched subjects still at risk of initiating treatment after $T_i$ (`Y_z > T_i`), including subjects who will start treatment later. Because matched subjects are removed, a later-treated subject who is used as a control is not matched again as a treated subject; the number of such subjects is returned as `treated.to.control`.

For a per-protocol analysis, a matched control who later initiates treatment at $T_j$ is censored at that time, contributing follow-up $T_j - T_i$ with no event. The log hazard ratio is estimated with a frailty Cox model and with an IPCW-adjusted frailty Cox model that corrects for this informative censoring.

The IPCW steps are as follows. The matched data are split at $v = 0.5$ and at every observed event and switching time using `ipcwswitch::cens.ipw()` and `ipcwswitch::replicRows()`. Weights are computed by `ipcw.new()` as the inverse of the probability of remaining uncensored by treatment switching, estimated from arm-specific Cox models for the switching time with $X_1$, $X_2$ and the time-varying $X_3$. Weights truncated at the 0.1st and 99.9th percentiles are also computed (`weights.trunc`), but the reported IPCW estimates use the untruncated weights. If the weights cannot be computed for a replicate, that replicate's IPCW estimate is `NA`.

### Naive model (`estimated_naive_value`)

The naive approach ignores the timing of treatment initiation. Treatment is coded as a fixed binary exposure (treated before the observed event or censoring time), and the propensity score is estimated by logistic regression on $X_1$, $X_2$ and the baseline value $X_3(0)$. Subjects are grouped by full matching on this score (`MatchIt::matchit(method = "full")`). Within each group, follow-up starts at the treated subject's initiation time, controls whose follow-up ends before that time are dropped, and the treated subject is kept with one randomly selected control. The effect is estimated with a Cox model with a robust variance clustered by group.

### True treatment effect (`true_value`)

Because the hazard ratio is non-collapsible, the target effect among the treated is not simply the data-generating $\theta = -0.5$, so it is approximated numerically. A large population is simulated (`sample_size`, for example 100,000 subjects), and for every treated subject a counterfactual never-treated event time is drawn from the same covariates using an independent uniform draw and the same censoring time. Treated subjects whose counterfactual observed time is at least $T_i$ are kept, and their observed and counterfactual post-treatment outcomes are stacked as pairs. A frailty Cox model with a pair-level random effect is then fit, and the coefficient of the treatment indicator is taken as the true log hazard ratio.

## Performance measures

For $R$ replicates with estimates $\hat\theta_r$, estimated standard errors $\widehat{SE}_r$ and true value $\theta_0$ (`true_phi`):

| Measure | Definition |
|---|---|
| Bias | $\frac{1}{R}\sum_r \hat\theta_r - \theta_0$ |
| MAE (mean absolute error) | $\frac{1}{R}\sum_r \lvert \hat\theta_r - \theta_0 \rvert$ |
| ASE (average standard error) | $\frac{1}{R}\sum_r \widehat{SE}_r$ |
| MCSD (Monte Carlo SD) | Standard deviation of $\hat\theta_1, \dots, \hat\theta_R$ |
| MSE (mean squared error) | $\frac{1}{R}\sum_r (\hat\theta_r - \theta_0)^2$ |
| Coverage | Proportion of replicates whose 95% Wald interval $\hat\theta_r \pm 1.96\,\widehat{SE}_r$ contains $\theta_0$ |

## Getting started

### Requirements

The code is written in R. Install the CRAN packages used by the scripts (`plyr` is called as `plyr::join()`, so it must be installed even though it is not loaded):

```r
install.packages(c("survival", "optmatch", "MASS", "coxme", "dplyr", "MatchIt", "plyr"))
```

The `ipcwswitch` package has been archived on CRAN. Install it from the archived source:

```r
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/ipcwswitch/ipcwswitch_1.0.4.tar.gz",
  repos = NULL,
  type = "source"
)
```

Installing with `repos = NULL` does not install dependencies automatically, so if the installation fails because of a missing package, install that package from CRAN first and then rerun the command above.

### Running the simulations

Each script defines its function and ends with a worked example and summary code. Each script also begins with `rm(list = ls())`, so running a whole script clears the R workspace before defining the function and running its example. Before calling `estimated_allrisk_value()`, run the helper `ipcw.new()` defined at the top of `all_risk_function.R`.

A typical workflow for one scenario:

```r
# 1. True log hazard ratio for the scenario (lambda_z = 0.7, about 40% treated)
tv <- true_value(lambda_z = 0.7, sample_size = 100000)
true_phi <- tv$true            # about -0.4900718

# 2. Monte Carlo simulations for each method
res_naive <- estimated_naive_value(true_phi, sample_size = 1000, simulation = 1000, lambda_z = 0.7)
res_pcm   <- estimated_purecontrol_value(true_phi, sample_size = 1000, simulation = 1000, lambda_z = 0.7)
res_arm   <- estimated_allrisk_value(true_phi, sample_size = 1000, simulation = 1000, lambda_z = 0.7)
```

A helper for summarizing any estimator:

```r
summarize_sim <- function(est, se, ci, truth) {
  ok <- !is.na(est)
  est <- est[ok]; se <- se[ok]; ci <- ci[ok]
  c(Bias     = mean(est) - truth,
    MAE      = mean(abs(est - truth)),
    ASE      = mean(se),
    MCSD     = sd(est),
    MSE      = mean((est - truth)^2),
    Coverage = mean(ci))
}

summarize_sim(res_naive$naive_est, res_naive$naive_se, res_naive$naive_ci, true_phi)
summarize_sim(res_pcm$frailty_est, res_pcm$frailty_se, res_pcm$frailty_ci, true_phi)
summarize_sim(res_pcm$cox_est,     res_pcm$cox_se,     res_pcm$cox_ci,     true_phi)
summarize_sim(res_arm$frailty_est, res_arm$frailty_se, res_arm$frailty_ci, true_phi)
summarize_sim(res_arm$ipcw_est,    res_arm$ipcw_se,    res_arm$ipcw_ci,    true_phi)
```

### Reproducibility and runtime

Replicate `sim` uses `set.seed(sim * 2021 + 1207 * 5 + w * 100)`, where `w` starts at 1 and increases by one whenever the IPCW weights fail for a replicate in `estimated_allrisk_value()`. `true_value()` uses `set.seed(123)`. The matching loops print the index of each processed treatment time, and the PCM and ARM simulations fit several Cox and frailty models per replicate, so a full run of 1,000 replicates takes a while.

## Citation

If you use this code, please cite:

<!-- TODO: add volume, pages and DOI once available -->
```bibtex
@article{wu2026dptm,
  title   = {Dynamic Propensity Trajectory Modeling and Matching with Time-Dependent Covariates for Estimating Treatment Effect},
  author  = {Wu, Xue and Yang, Yifan and Xu, Yifei and Li, Liang and Zhang, Lijun and Rahman, Mahboob and McCoy, Rozalina G. and Chinchilli, Vernon M. and Chen, Chixiang and Wang, Ming},
  journal = {Statistics in Medicine},
  year    = {2026}
}
```

## References

Graffeo N, Latouche A, Le Tourneau C, Chevret S. ipcwswitch: An R package for inverse probability of censoring weighting with an application to switches in clinical trials. *Computers in Biology and Medicine* 2019; 111: 103339.

Li YP, Propert KJ, Rosenbaum PR. Balanced risk set matching. *Journal of the American Statistical Association* 2001; 96(455): 870–882.

Lu B. Propensity score matching with time-dependent covariates. *Biometrics* 2005; 61(3): 721–728.
