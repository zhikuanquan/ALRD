# ALRD

`ALRD` implements the longitudinal robust direct learning methods from
*Doubly Robust Estimation of Individualized Treatment Effects in
Longitudinal Data*. The package directly estimates time-specific
individualized treatment effects (ITEs) from repeated outcome
measurements while accounting for within-subject correlation.

The paper writes treatment as `T = 1` for treatment and `T = -1` for
control and uses the decomposition

``` text
Y_ij = m(X_ij) + T_i * delta(X_ij) / 2 + error_ij.
```

This package accepts treatment coded as either `-1/1` or `0/1`;
internally the learning objective is fit on the `-1/1` scale.

## Algorithms

`fit_lrd()` fits LRD-Learning by minimizing the inverse-propensity
weighted longitudinal loss

``` text
sum_i 1 / pi_T(Z_i) * {Y_i - T_i f(X_i) / 2}' V_i^{-1}
  {Y_i - T_i f(X_i) / 2}.
```

`fit_alrd()` fits ALRD-Learning with the paper's two-stage
augmentation:

1. Estimate the main effect `m(X)` from the weighted working covariance
   loss.
2. Estimate the ITE from the augmented outcome `Y - mhat(X)` using
   `T * f(X) / 2`.

For `fit_alrd()`, use `penalty = "none"` for the non-penalized version
and `penalty = "lasso"` for lasso variable selection. The lasso version
uses clustered cross-validation, so all repeated measurements from the
same subject stay in the same fold. By default, it refits the selected
ITE terms without penalty to reduce shrinkage bias, matching the
real-data analysis strategy in the paper.

For exchangeable and AR(1) working covariance structures, `rho` is estimated
iteratively when it is omitted. The fitting functions repeat: fit the model,
estimate `rho` from residuals, and refit until convergence. To fix the working
correlation instead, pass `rho` directly, for example `rho = 0.5`.

Key lasso outputs include:

``` text
fit$lambda              # selected lambda
fit$lambda_1se          # one-SE lambda
fit$lasso_coefficients  # penalized ITE coefficients
fit$coefficients        # final coefficients, refit if refit = TRUE
fit$selected_terms      # selected ITE design terms
fit$selected_covariates # selected treatment-effect covariates
fit$main_model          # first-stage main-effect fit
```

## Installation

``` r
install.packages(c("devtools", "roxygen2", "testthat", "glmnet"))
devtools::install()
```

## Toy Simulation Example

`simulate_alrd_data()` follows the paper's simulation structure: five
repeated measurements, baseline continuous and binary covariates, a
random subject intercept inducing exchangeable correlation, and
scenarios with correct or misspecified main-effect and propensity-score
models.

``` r
library(ALRD)

dat <- simulate_alrd_data(
  n = 120,
  scenario = 1,
  main_effect_size = "large",
  correlated = TRUE,
  seed = 1
)

ite_covariates <- c(
  "time",
  paste0("z", 1:20),
  paste0("time_z", 1:20)
)
```

### Non-Penalized ALRD

``` r
fit_np <- fit_alrd(
  data = dat,
  id = "id",
  time = "time",
  treatment = "treatment",
  outcome = "outcome",
  covariates = ite_covariates,
  ps = 0.5,
  penalty = "none",
  covariance = "exchangeable"
)

head(predict_ite(fit_np, dat))
```

### Lasso ALRD With Unpenalized Refit

``` r
fit_lasso <- fit_alrd(
  data = dat,
  id = "id",
  time = "time",
  treatment = "treatment",
  outcome = "outcome",
  covariates = ite_covariates,
  ps = 0.5,
  penalty = "lasso",
  refit = TRUE,
  nfolds = 5,
  seed = 10,
  covariance = "exchangeable"
)

fit_lasso$lambda
fit_lasso$selected_covariates
head(predict_ite(fit_lasso, dat))
```

### Evaluating ITE Recovery

The simulation data include the true time-specific ITE, so a quick check
can compare the predicted ITE with the truth.

``` r
pred <- predict_ite(fit_lasso, dat)
cor(pred, dat$true_ite, method = "spearman")
mean(abs(pred - dat$true_ite))
```

## Development

``` r
devtools::document()
devtools::test()
devtools::check()
```

Additional documentation:

- [Package usage guide](docs/package_usage_guide.md)
- [Repository contents and GitHub Desktop guide](docs/repository_contents.md)
- [Collaborator onboarding guide](docs/collaborator_onboarding.md)
- [GitHub collaboration workflow](docs/github_collaboration_workflow.md)

## Version tracker

See [NEWS.md](NEWS.md) for the package development log.
