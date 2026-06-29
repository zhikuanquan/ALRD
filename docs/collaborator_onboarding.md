# Collaborator Onboarding Guide

This guide is intended for collaborators joining the ALRD R package project.

## 1. Project goal

The goal is to develop an R package for ALRD/LRD methods for estimating
individualized treatment effects in longitudinal data.

The current package scope is:

- ALRD and LRD learning;
- non-penalized and lasso versions;
- propensity-score weighting;
- working covariance structures;
- ITE prediction;
- validation simulations.

Joint-model extensions are postponed until the core ALRD package is stable.

## 2. Install the package from GitHub

After the repository is uploaded to GitHub:

```r
install.packages("devtools")
devtools::install_github("YOUR_GITHUB_USERNAME/ALRD")
```

Then load:

```r
library(ALRD)
```

## 3. Clone for development

Using GitHub Desktop:

1. Accept the GitHub collaborator invitation.
2. Open GitHub Desktop.
3. Click **File > Clone Repository**.
4. Select the `ALRD` repository.
5. Choose a local folder.
6. Click **Clone**.

Using command line:

```bash
git clone https://github.com/YOUR_GITHUB_USERNAME/ALRD.git
cd ALRD
```

## 4. Install development dependencies

In R:

```r
install.packages(c(
  "devtools",
  "roxygen2",
  "testthat",
  "glmnet",
  "knitr",
  "rmarkdown"
))
```

Install the local package:

```r
devtools::install()
```

## 5. Run the toy example

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

This is the recommended first check for new collaborators.

## 6. Development workflow

Do not work directly on `main`.

Create a branch for each focused task:

```bash
git checkout main
git pull origin main
git checkout -b feature/my-task
```

Examples:

```text
feature/alrd-loss
feature/lasso-refit
feature/covariance-ar1
feature/predict-method
docs/readme-methodology
tests/core-api
```

Make changes, then commit:

```bash
git status
git add R/file_changed.R tests/testthat/test-file.R
git commit -m "Implement focused task"
git push origin feature/my-task
```

Open a pull request on GitHub.

## 7. Pull request checklist

Before opening a pull request:

- [ ] Code runs locally.
- [ ] Relevant tests pass.
- [ ] Documentation is updated if function arguments changed.
- [ ] No generated outputs or large simulation files are committed.
- [ ] The PR description explains statistical impact.

Recommended checks:

```r
devtools::document()
devtools::test()
devtools::check()
```

## 8. Statistical review checklist

For methodology-related pull requests, reviewers should check:

- Does the code match the ALRD/LRD objective?
- Are treatment codings handled correctly?
- Are propensity-score weights applied correctly?
- Is the working covariance used correctly?
- Are lasso selected covariates interpreted correctly?
- Are simulation metrics computed correctly?

## 9. Coding review checklist

For implementation pull requests, reviewers should check:

- Is the code modular?
- Are function names and arguments consistent?
- Are errors informative?
- Are edge cases tested?
- Does the package still install and load?
- Are examples reproducible?

## 10. Recommended first tasks

1. Freeze v0.1 API and output structure.
2. Add `summary.alrd()` and `summary.lrd()` methods.
3. Add `predict.alrd_fit()` if desired in addition to `predict_ite()`.
4. Strengthen unit tests for lasso/refit outputs.
5. Improve rho diagnostics.
6. Expand Scenario 2 validation.
7. Update README and vignette after API freeze.
