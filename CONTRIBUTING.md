# Contributing to ALRD

This repository is being developed as a research R package for ALRD/LRD
methods for individualized treatment-effect learning in longitudinal data.

## Collaboration model

We use a branch and pull-request workflow.

- `main` contains stable package code only.
- New work should be done on a separate branch.
- Pull requests should be reviewed before merging.
- Generated outputs and large simulation results should not be committed to
  GitHub.

## Roles

### Statistical lead

The statistical lead is responsible for:

- Defining the target estimand and assumptions.
- Specifying the LRD/ALRD objective functions.
- Defining augmentation, propensity-score, and working-covariance choices.
- Designing simulation scenarios.
- Interpreting simulation results.
- Reviewing whether implementation matches the methodology.

### Coding lead

The coding lead is responsible for:

- Refactoring package internals.
- Implementing stable user-facing R functions.
- Standardizing returned objects and S3 methods.
- Adding tests, examples, and documentation.
- Ensuring `R CMD check` or `devtools::check()` passes.

## Branch naming

Use short descriptive branch names:

```text
feature/alrd-loss
feature/lasso-refit
feature/covariance-ar1
feature/predict-method
feature/simulation-scenario2
docs/readme-methodology
tests/core-api
fix/rho-diagnostics
```

## Local development workflow

Start from the latest `main`:

```bash
git checkout main
git pull origin main
git checkout -b feature/my-task
```

After making changes:

```bash
git status
git add R/file_changed.R tests/testthat/test-file.R
git commit -m "Implement focused change"
git push origin feature/my-task
```

Then open a pull request on GitHub.

## Pull request expectations

Each pull request should answer:

- What changed?
- Does this affect the statistical method?
- What tests or examples were run?
- Are there open questions for the reviewer?

Keep pull requests focused. Avoid mixing unrelated changes such as method
changes, simulation changes, and formatting updates in one PR.

## Package checks

Before asking for review, run:

```r
devtools::document()
devtools::test()
devtools::check()
```

If a full check is slow, at least run the relevant tests and a toy example.

## Files that should not be committed

Do not commit:

- `.Rhistory`
- `.RData`
- `.DS_Store`
- generated package tarballs
- `*.Rcheck/`
- generated outputs
- large simulation outputs
- generated preview files

Exploratory simulation scripts and raw simulation outputs should be stored
outside this package-only repository unless they are promoted to formal package
examples, tests, or vignettes.

## Version scope

For v0.1, focus on ALRD/LRD only:

- `fit_alrd()`
- `fit_lrd()`
- `predict_ite()`
- `simulate_alrd_data()`
- non-penalized and lasso versions
- selected covariate output
- covariance options
- toy examples
- validation examples

Joint-model extensions should be postponed until the core ALRD package is
stable.
