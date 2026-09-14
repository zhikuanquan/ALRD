# ALRD code review and collaborator work plan

Follow-up: the focused weekly update in `weekly_method_check.md` verifies the
method equations against the located 23-page manuscript and fixes covariance
domain checks and correlation-update bounds. The findings below preserve the
original pre-fix audit; its covariance failures and source-unavailable statement
are historical, not the current implementation status.

Review date: 2026-09-13. Scope: the six exported functions, their internal
helpers, registered S3 methods, existing tests, vignette, and package build.
Reviewed development version: 0.0.0.9000 in Documents/ALRD. The GitHub checkout
was not used as the source of this review.

## Assessment

The package contains a useful working prototype of longitudinal direct and
augmented learning. It is ready for controlled development examples with a
complete numeric basis, positive subject-constant weights, and admissible
correlation settings. It is not yet ready for broad scientific use or a JSS
submission: prediction contracts, covariance checks, and full-pipeline
validation need work. Passing the existing tests does not establish statistical
validity across the public API.

The reference article, Huling and Yu (2021), *Subgroup Identification Using the
personalized Package*, JSS 98(5), DOI 10.18637/jss.v098.i05, organizes its
software section around the complete user workflow. ALRD currently implements
the fitting core but lacks much of that diagnostic and validation workflow.
Its multi-treatment and multiple-outcome features should not be copied simply
to match the reference paper's breadth.

This is an implementation review, not verification of the ALRD method paper's
theorems or a new simulation study. The referenced ALRD methodology paper was
not supplied in this task. In particular, two-stage weighting alone does not
justify an unconditional claim of double robustness for every fitted basis,
penalty, covariate history, and covariance specification.

## Reproducible evidence

Run from the development repository:

```sh
Rscript -e 'testthat::test_local(reporter="summary")'
Rscript scripts/audit_package_contracts.R
Rscript scripts/build_help.R
R CMD build . --no-manual
```

The audit is diagnostic, not a regression suite designed to pass. Convert each
defect into a regression test when implementing its fix. Audit output and
session information are saved in `output/review/contract-audit.txt`.

Verification completed: `R CMD build . --no-manual` succeeded and
`R CMD check --no-manual --output=output/review ALRD_0.0.0.9000.tar.gz`
finished with **Status: OK**, including installed examples, Rd cross-references,
tests, and rebuilding vignette outputs. Logs are in `output/review/build.log`
and `output/review/check.log`. This was a normal package check, not
`--as-cran`; PDF manual generation was excluded. R 4.4.2 and glmnet 5.0 were used.
The separate audit defects remain present despite a clean package check.

Parsed R expressions in all six existing public-function source files were
compared to their Git HEAD versions and found identical after documentation
edits. The new package overview adds help metadata, not an estimator.

| Check | Observed result | Interpretation |
|---|---|---|
| Existing suite | Six tests, 21 expectations pass | Basic fitting and return-object checks work. |
| Subject-constant weighted LRD vs direct matrix GLS | Maximum coefficient difference 6.7e-16 | The supported scalar-weight whitening calculation agrees numerically. |
| Reverse input order with aligned row weights, fixed rho | Maximum coefficient difference 3.6e-15 | Ordering is consistent in this numeric example. |
| One zero observation weight | 200 training rows, 199 residuals | Fit vectors and metadata can lose alignment. |
| Exchangeable rho = -0.5 with five visits | Minimum eigenvalue approximately -1 | Public covariance helper returns an invalid correlation matrix. |
| AR(1), times 0 and 0.5, rho = -0.5 | Non-finite off-diagonal entries | Fractional time and negative rho need a domain restriction. |
| Reordered factor levels in prediction | 200 non-finite predictions | Training contrasts and levels are not preserved. |
| One missing prediction covariate | 200 rows, 199 predictions | Missing values silently drop rows. |
| Duplicate design column | 200 non-finite predictions | Rank deficiency propagates through coefficients. |
| Two lasso folds | Error from glmnet | The helper's allowance of two folds conflicts with the solver. |
| rho = 0.5 and estimate_rho = TRUE | rho updated to approximately 0.7517 | Supplied rho is fixed only by default, not unconditionally. |
| Varying row weights: implemented vs symmetric weighted quadratic form | 27.33333 vs 31.57107 | The weight placement changes the loss; define the intended extension. |

## Findings requiring implementation work

### P1: Preserve prediction encoding and row count

In `R/predict_ite.R`, both methods reconstruct a model matrix from stored terms
without stored training factor levels and contrasts. They then index the
coefficients by new matrix column names. Reordered levels can change the
reference category and request nonexistent coefficients. Missing covariates
are checked, but missing values are not, allowing rows to disappear.

Store the training levels, contrasts, and design-column schema. Validate new
levels explicitly, align the new design to the training columns, and define an
NA policy that either errors or returns an NA in the original row position.
Add a test where prediction on the same values with reordered factor levels is
unchanged, plus tests for single-level subsets, novel levels, and missing rows.

### P1: Resolve zero weights and covariance domain validation

`.fit_gls()` drops zero-weight rows internally while the fitted-object metadata
and correlation-update inputs retain all rows. This can shorten residuals and
misalign subsequent rho updates. Decide whether zero weights mean explicit
exclusion or retained rows with zero contribution; preserve a row map and
consistent outputs under that decision.

`.default_rho()` checks only a generic interval, while exchangeable correlation
also requires $\rho > -1/(m-1)$ for every subject. The residual updater clamps to
[-0.95, 0.95], which does not ensure that condition. `.correlation_matrix()`
allows duplicate AR(1) times and negative rho with fractional lags. Validate
finite times, unique subject-time pairs, positive definiteness, and admissible
rho before solving; constrain iterative updates to the same valid domain.

### P1: Specify weights before extending the statistical claims

In `.fit_gls()` and `.transform_gls()`, the Cholesky factor $C$, satisfying
$C^{\mathsf{T}}C = V^{-1}$, is applied first, then $W^{1/2}$.
The quadratic form is

$$
\mathbf{e}^{\mathsf{T}} C^{\mathsf{T}} W C \mathbf{e}.
$$

It equals the intended scalar subject-weighted form when $W$ is constant
within subject. It does not generally equal

$$
\mathbf{e}^{\mathsf{T}} W^{1/2} V^{-1} W^{1/2} \mathbf{e}.
$$

Neither arbitrary-row-weight extension
should be called correct without a defined estimand and loss.

For the first release, consider requiring subject-constant propensity and
analysis weights for the baseline-treatment objective. Otherwise write the
intended row-weighted objective and derive a matching solver before changing
weight placement. Keep the direct GLS oracle test for the supported case.

### P1: Validate the complete two-stage pipeline

`fit_alrd()` computes the main-effect predictions and residualized response
before the ITE-stage `.fit_lasso_gls()` call. Each stage's glmnet folds preserve
subjects, but ITE validation outcomes have already influenced the full-sample
main-effect fit. Automatic propensity estimation and rho estimation also use
the full fitting sample. Clustered folds alone do not provide a held-out
evaluation of the entire procedure.

Add outer subject splits and estimate the propensity model, augmentation,
working correlation, selection, and refit from each training split only. Specify
whether tuning and final performance evaluation use distinct nested layers.
Decide with the method lead whether cross-fitting is also required for the
estimator itself; that is a separate methodological choice from honest testing.

### P1: Enforce the baseline-assignment and rank contracts

The public fit functions accept treatment changes within a subject, although
their stated model uses a single T_i. The existing first two tests and original
vignette used row-varying treatments. Replace those examples/tests with the
intended design and add explicit unsupported-design checks.

Automatic ALRD propensity estimation fits the treatment on all visits using
the outcome basis, including time if supplied. This is not automatically a
baseline-confounder model. With unequal visits it also gives greater likelihood
weight to subjects with more rows. Provide separate propensity covariates or a
subject-level score interface. Decide an explicit common default for ps across
the two fit functions rather than silently assuming observational or randomized
assignment.

QR fitting can return NA coefficients for aliased columns; multiplication by
that coefficient vector can make all predictions NA. Decide whether to reject
rank deficiency or support estimable predictions with a documented alias policy.

### P2: Finish solver controls and diagnostics

- Require at least three folds for the current glmnet workflow; validate integer
  nfolds, rho_max_iter, simulator n/scenario, and logical controls explicitly.
- Separate tuning controls for augmentation and ITE; expose lambda.min versus
  lambda.1se as a deliberate policy. Record fold IDs for reproducibility.
- Use explicit intercept identity rather than matching any column name
  containing "Intercept" to decide penalty exemptions.
- Distinguish selected columns, intercept inclusion, and original variables.
  Forced intercept inclusion during refit may differ from selected_terms when
  its penalized coefficient is numerically zero.
- Warn on failed correlation convergence and expose fitting diagnostics through
  dedicated summary methods. A generic list summary is not a model summary.
- Add row/data validation for Inf, empty data, missing treatment, and treatment
  level counts. Existing complete.cases() checks do not reject Inf.
- Consolidate duplicated GLS and prediction code. Preallocate transformed
  matrices rather than growing them with repeated rbind(). Profile before
  selecting performance targets or claiming speed improvements.

## Function-by-function additions

| Function | Existing capability | Next addition |
|---|---|---|
| fit_lrd | Unpenalized direct GLS; three correlations; iterative rho | Baseline data contract; valid covariance; consistent weight/row handling; explicit propensity policy. |
| fit_alrd | Two-stage fitting; lasso; ITE refit | Separate nuisance/basis controls; full-pipeline validation; convergence diagnostics. |
| estimate_ps | Logistic fit with probability clipping | Subject-level workflow; raw scores, clipping counts, overlap/balance diagnostics, convergence checks. |
| working_covariance | Named subject correlation matrices | Domain/positive-definiteness checks; row/time labels; clearly documented time units. |
| predict_ite | Training/new-data contrast evaluation | Stable factor encoding and output row count; tested new-data schema. |
| simulate_alrd_data | Five scenarios, known ITE/main effect/propensity | Strict argument checks; unbalanced visits and explicit error-correlation controls; verified paper mapping. |
| print methods | Brief fit display | summary and plot methods with structured diagnostics, without inventing inferential quantities. |

## Collaborator tasks aligned with the reference paper

Assignments below are proposals for team agreement, not messages sent to
collaborators. The advisor's suggested writing division is preserved: Zhikuan
leads Sections 1-2; Steven leads Sections 3 onward. Software work crosses that
division when methods and implementation must agree.

| ID / priority | Owner / reviewer | Specific deliverable | Completion criterion |
|---|---|---|---|
| M01 / first | Zhikuan / Shuai | Section 2 specification: estimand, timing, confounders, two losses, weight placement, admissible V, and exact ALRD theorem conditions | Every equation mapped to code and method-paper source; limitations and model classes stated. |
| E01 / P1 | Steven / Zhikuan | Input validation, factor prediction, row mapping, and rank policy | Audit failures become regression tests that pass; numeric predictions stay consistent. |
| E02 / P1 | Zhikuan / Steven | Weight and covariance specification plus implementation | Direct GLS oracle agrees; invalid rho/time combinations reject clearly; zero-weight semantics tested. |
| E03 / P1 | Steven / Zhikuan | Baseline propensity workflow and explicit ps policy | Unequal visit counts do not change subject-level propensity fitting; treatment changes reject clearly. |
| V01 / P1 | Steven / Zhikuan and Shuai | Full-pipeline subject resampling | No test outcome enters any training-stage estimate; all nuisance fitting and tuning recorded by split. |
| E04 / P2 | Steven / Zhikuan | Tuning controls, intercept identity, selected-variable mapping | Separate stage settings work; folds and lambda policy documented; intercept edge cases tested. |
| W01 / P2 | Steven / Shuai | Section 3 diagnostics: overlap, convergence, ITE plots, summary methods | One complete example uses exported functions and reports relevant diagnostics. |
| W02 / P2 | Steven / Zhikuan | Out-of-sample subgroup/policy evaluation for longitudinal outcomes | Scientific target, time aggregation, subject weighting and uncertainty method approved before implementation. |
| D01 / now | Steven / Zhikuan | Adopt and maintain new function help and executable vignette | All examples execute against the installed package; arguments and return fields match implementation. |
| S01 / after M01, V01 | Steven / Zhikuan and Shuai | Section 4 simulation protocol and comparator adapters | Same data/splits and fair tuning; replicate-level output, failures, timing and Monte Carlo uncertainty retained. |
| A01 / after W01, W02 | Steven / Shuai | Section 5 application protocol and reproducible data preparation | Dataset/permissions agreed; missingness, baseline confounding and external evaluation addressed. |
| R01 / before release | Steven / all | Release/replication bundle, CI, package metadata and installation docs | Clean package checks plus substantive audit fixes; fresh environment reproduces examples and manuscript results. |
| M02 / after evidence | Zhikuan: Sections 1-2; Steven: 3-6 / Shuai | Final manuscript with supported claims | Narrative, code outputs and frozen release agree; no unsupported superiority or robustness statements. |

Suggested sequence: agree M01 and the application target first; complete E01-E03
and V01 before final simulation runs; then E04/W01/W02 and D01; finally S01,
A01, R01 and M02. Hold joint-model, multi-category-treatment and new outcome-family
extensions until the core baseline-treatment release meets these criteria.

## Documentation delivered in this review

All six exported function help pages now define parameters, return values,
examples, coding conventions, and known limitations. A package overview and
shared print-method topic were added. The introductory vignette was rewritten
with baseline treatment, subject holdout, both fitting paths, lasso, covariance
inspection, ITE trajectory prediction, and baseline propensity mapping.

The usage guide now uses the actual simulator column `propensity` (previously
incorrectly `ps`) and the actual nested `rho_estimation` fields. AR(1) uses
actual time gaps; supplied rho can be estimated when explicitly requested; the
lasso applies to both stages. New-data and validation limitations are explicit.

The browsable manual is generated by `scripts/build_help.R` into `output/help/`.
Help examples are executable, and the rendered vignette contains real generated
output. R source edits in the six existing function files are documentation
comments only; estimator bodies have not been altered by this review.

Build exclusions were added for local agent instructions, manuscript, output,
scratch files, and development scripts. Scripts remain in the repository and
must be included separately in the JSS replication bundle where appropriate.

## Manuscript and repository follow-up

The earlier manuscript needs the same corrections before circulation: AR(1)
uses actual time differences, complete-case data are required rather than
automatically selected, and clustered glmnet folds do not constitute complete
two-stage validation. Its robustness paragraph needs the actual method-paper
citation and conditions. This task did not regenerate that PDF.

Existing onboarding documents still contain placeholder GitHub URLs. Resolve
those in a later repository-wide documentation pass after the approved local
workflow is established. Paper authorship does not automatically establish
package contributor roles; confirm DESCRIPTION credits separately.
