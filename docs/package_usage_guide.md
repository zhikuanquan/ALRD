# ALRD Package Usage Guide

This guide explains the purpose, main functions, inputs, outputs, and common
workflows for the `ALRD` R package.

## 1. Package purpose

`ALRD` implements longitudinal residual difference learning methods for
estimating individualized treatment effects (ITEs) from repeated outcome data.

The package focuses on the model decomposition

```text
Y_ij = m(X_ij) + T_i * delta(X_ij) / 2 + error_ij,
```

where:

- `Y_ij` is the outcome for subject `i` at time `j`.
- `T_i` is treatment, coded internally as `-1` for control and `1` for
  treatment.
- `m(X_ij)` is the nuisance main effect.
- `delta(X_ij)` is the individualized treatment effect.
- The goal is to estimate `delta(X)`.

The package currently supports:

- LRD learning.
- ALRD learning with main-effect augmentation.
- Non-penalized ALRD.
- Lasso ALRD using `glmnet`.
- Optional unpenalized refit after lasso selection.
- Working covariance structures: independence, exchangeable, and AR(1).
- Propensity scores supplied by the user or estimated internally.
- Toy data simulation for local checks and examples.

## 2. Installation

After uploading the package to GitHub, collaborators can install it with:

```r
install.packages("devtools")
devtools::install_github("YOUR_GITHUB_USERNAME/ALRD")
```

For local development from the package folder:

```r
install.packages(c("devtools", "roxygen2", "testthat", "glmnet"))
devtools::install()
```

## 3. Main functions

### `simulate_alrd_data()`

Creates toy longitudinal data following the simulation structure used in the
ALRD project.

```r
dat <- simulate_alrd_data(
  n = 200,
  scenario = 1,
  main_effect_size = "large",
  correlated = TRUE,
  seed = 1
)
```

Important output columns:

- `id`: subject identifier.
- `time`: repeated-measurement time.
- `treatment`: treatment coded as `-1` or `1`.
- `outcome`: longitudinal outcome.
- `true_ite`: true ITE used for simulation validation.
- `ps`: true propensity score.
- `z1` to `z20`: baseline covariates.
- `time_z1` to `time_z20`: time-by-covariate interactions.

### `fit_lrd()`

Fits longitudinal residual difference learning without main-effect
augmentation.

```r
fit_lrd_obj <- fit_lrd(
  data = dat,
  id = "id",
  time = "time",
  treatment = "treatment",
  outcome = "outcome",
  covariates = c("time", paste0("z", 1:20)),
  ps = 0.5,
  covariance = "exchangeable"
)
```

Use this when you want a direct baseline estimator of the ITE function.

### `fit_alrd()`

Fits augmented longitudinal residual difference learning.

```r
ite_covariates <- c(
  "time",
  paste0("z", 1:20),
  paste0("time_z", 1:20)
)

fit_alrd_obj <- fit_alrd(
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
```

Use `penalty = "none"` for non-penalized ALRD.

Use `penalty = "lasso"` for variable selection:

```r
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
```

### `predict_ite()`

Predicts the estimated ITE for new or training data.

```r
pred_ite <- predict_ite(fit_lasso, dat)
head(pred_ite)
```

For simulated data, compare predictions with the truth:

```r
cor(pred_ite, dat$true_ite, method = "spearman")
mean(abs(pred_ite - dat$true_ite))
```

### `estimate_ps()`

Estimates propensity scores using logistic regression.

```r
ps_fit <- estimate_ps(
  data = dat,
  treatment = "treatment",
  covariates = paste0("z", 1:4)
)

fit_obs <- fit_alrd(
  data = dat,
  id = "id",
  time = "time",
  treatment = "treatment",
  outcome = "outcome",
  covariates = ite_covariates,
  ps = ps_fit$scores,
  penalty = "lasso"
)
```

### `working_covariance()`

Constructs subject-level working correlation matrices.

```r
cov_list <- working_covariance(
  data = dat,
  id = "id",
  time = "time",
  covariance = "ar1",
  rho = 0.5
)
cov_list[[1]]
```

## 4. Propensity-score options

The fitting functions accept several propensity-score workflows.

### Randomized trial

If treatment is randomized 1:1:

```r
ps = 0.5
```

### Known propensity score

If the true or design-based propensity score is known:

```r
ps = dat$ps
```

### Estimated propensity score

Estimate using `estimate_ps()` and pass the estimated vector:

```r
ps_fit <- estimate_ps(dat, treatment = "treatment", covariates = paste0("z", 1:4))
ps = ps_fit$scores
```

## 5. Working covariance options

Use:

```r
covariance = "independence"
covariance = "exchangeable"
covariance = "ar1"
```

If `rho` is omitted for exchangeable or AR(1), the package estimates it
iteratively from residuals:

```r
fit <- fit_alrd(
  data = dat,
  id = "id",
  time = "time",
  treatment = "treatment",
  outcome = "outcome",
  covariates = ite_covariates,
  ps = 0.5,
  penalty = "lasso",
  covariance = "ar1",
  rho = NULL,
  estimate_rho = TRUE
)
```

To fix the working correlation:

```r
rho = 0.5
```

Current simulation evidence suggests that residual-based `rho` estimation may
overestimate the true AR(1) correlation in some settings. Treat `rho` output as
a diagnostic rather than the primary scientific estimand.

## 6. Important fitted-object outputs

For `fit_alrd()` with lasso:

```r
fit_lasso$coefficients
fit_lasso$lasso_coefficients
fit_lasso$selected_terms
fit_lasso$selected_covariates
fit_lasso$lambda
fit_lasso$lambda_1se
fit_lasso$rho
fit_lasso$rho_iterations
fit_lasso$rho_converged
fit_lasso$main_model
```

Interpretation:

- `coefficients`: final ITE coefficients, refit if `refit = TRUE`.
- `lasso_coefficients`: penalized coefficients from the lasso stage.
- `selected_terms`: selected model-matrix terms.
- `selected_covariates`: selected covariates interpreted as treatment-effect
  modifiers.
- `lambda`: cross-validated lasso tuning parameter.
- `rho`: final working correlation.
- `rho_converged`: whether iterative `rho` estimation converged.

## 7. Toy validation example

Run a complete package-use example:

```r
library(ALRD)

dat <- simulate_alrd_data(n = 100, scenario = 1, seed = 123)

covars <- c("time", paste0("z", 1:20), paste0("time_z", 1:20))

fit <- fit_alrd(
  data = dat,
  id = "id",
  time = "time",
  treatment = "treatment",
  outcome = "outcome",
  covariates = covars,
  ps = 0.5,
  penalty = "lasso",
  refit = TRUE,
  covariance = "exchangeable",
  seed = 1
)

fit$selected_covariates
pred <- predict_ite(fit, dat)
cor(pred, dat$true_ite, method = "spearman")
mean(abs(pred - dat$true_ite))
```

This example demonstrates:

- toy data simulation;
- non-penalized ALRD;
- lasso ALRD;
- selected covariates;
- ITE prediction;
- simple validation against the true simulated ITE.

## 8. Recommended first example for collaborators

```r
library(ALRD)

dat <- simulate_alrd_data(n = 100, scenario = 1, seed = 123)

covars <- c("time", paste0("z", 1:20), paste0("time_z", 1:20))

fit <- fit_alrd(
  data = dat,
  id = "id",
  time = "time",
  treatment = "treatment",
  outcome = "outcome",
  covariates = covars,
  ps = 0.5,
  penalty = "lasso",
  refit = TRUE,
  covariance = "exchangeable",
  seed = 1
)

fit$selected_covariates
pred <- predict_ite(fit, dat)
cor(pred, dat$true_ite, method = "spearman")
```
