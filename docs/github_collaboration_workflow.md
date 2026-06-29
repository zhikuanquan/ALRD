# GitHub Collaboration Workflow for the ALRD R Package

This document describes how to upload the ALRD R package to GitHub and how
collaborators should work together.

## 1. Create the GitHub repository

On GitHub:

1. Go to <https://github.com/new>.
2. Repository name: `ALRD`.
3. Recommended visibility during development: `Private`.
4. Do not initialize with README, `.gitignore`, or license because the local
   package already has these files.
5. Click **Create repository**.

After GitHub creates the empty repository, copy the repository URL. It will look
like one of these:

```text
https://github.com/<username>/ALRD.git
git@github.com:<username>/ALRD.git
```

If SSH keys are not configured, use the HTTPS URL.

## 2. Push the local package to GitHub

From the local package directory:

```bash
cd /Users/zhikuan_quan/Documents/GitHub/ALRD
git remote add origin https://github.com/<username>/ALRD.git
git push -u origin main
```

If GitHub asks for authentication, use your GitHub username and a personal
access token, or authenticate through the Git Credential Manager if prompted.

## 3. Invite collaborators

On GitHub:

1. Open the repository page.
2. Go to **Settings**.
3. Click **Collaborators and teams**.
4. Click **Add people**.
5. Enter each collaborator's GitHub username or email.
6. Recommended permission:
   - `Write` for active coding collaborators.
   - `Read` for collaborators who only review.
   - `Admin` only for a trusted co-maintainer.

## 4. Branch workflow

The `main` branch should contain stable code only. Collaborators should not work
directly on `main`.

Before starting new work:

```bash
git checkout main
git pull origin main
git checkout -b feature/descriptive-task-name
```

Example branch names:

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

## 5. Commit and push changes

After editing files:

```bash
git status
git add R/file_changed.R tests/testthat/test-file.R
git commit -m "Implement focused change"
git push origin feature/descriptive-task-name
```

Then open a pull request on GitHub.

## 6. Pull request expectations

Each pull request should include:

- Summary of what changed.
- Whether the statistical method changed.
- Tests or examples that were run.
- Open questions for reviewers.

Good pull requests are focused. Avoid combining unrelated tasks such as
methodology changes, simulation changes, and documentation rewrites in one PR.

## 7. Review responsibilities

Statistical review should check:

- Does the implementation match the ALRD/LRD objective?
- Are the nuisance-model choices documented correctly?
- Are selected covariates and ITE predictions interpreted correctly?
- Are simulation metrics appropriate?

Coding review should check:

- Is the code modular and readable?
- Are edge cases handled?
- Are tests included?
- Does the package still pass checks?

## 8. Package validation before merge

Before merging important changes, run:

```r
devtools::document()
devtools::test()
devtools::check()
```

For a command-line check:

```bash
R CMD build .
R CMD check --no-manual --no-build-vignettes ALRD_0.0.0.9000.tar.gz
```

The current initial package setup has passed built-package check with:

```text
Status: OK
```

## 9. What should not be committed

Do not commit:

- `.Rhistory`
- `.RData`
- `.DS_Store`
- generated package tarballs
- `*.Rcheck/`
- `outputs/`
- `presentations/`
- large raw simulation outputs

Keep exploratory simulation scripts and generated result files outside this
package-only repository unless they are promoted to formal package examples,
tests, or vignettes.

## 10. Suggested first GitHub issues

Create issues for:

1. Freeze v0.1 API and returned object structure.
2. Refactor ALRD loss construction.
3. Add `summary.alrd()` and `predict.alrd()` methods.
4. Improve covariance and rho diagnostics.
5. Expand Scenario 2 simulations.
6. Add tests for lasso refit and selected covariates.
7. Update README and vignette after API freeze.
