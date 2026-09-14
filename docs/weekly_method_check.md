# Weekly ALRD update: Professional package documentation

Date: September 13, 2026. Work completed in Documents/ALRD.

## This week's main deliverable

Produce one coherent package help system: an overview, complete help for all
six exported functions and the print methods, and one executable introduction
that takes a new user from installation to fitting and interpreting predictions.
The earlier method check is supporting work for documentation accuracy.

**Zhikuan:** verify the methodological explanations, data requirements,
treatment/propensity conventions, covariance definitions, return-value meanings,
and stated limitations. These checks and an initial documentation draft are
complete; review Steven's revisions for consistency with the implementation.

**Steven:** complete an editorial pass on the function help and introductory
vignette. Improve unclear argument descriptions and example explanations, add
missing user-facing details, and run the examples. Deliver the revised roxygen
comments and vignette plus a short list of unresolved questions. Work in
Documents/ALRD and regenerate the Rd pages rather than editing them directly.

**Definition of done:** every public function has a clear purpose, full argument
descriptions/defaults, structured return values, runnable examples, related-help
links and relevant limitations; the tutorial runs from start to finish; installed
help, examples and vignette pass package checks. No new statistical features or
large simulation runs are assigned this week. Shuai may monitor and advise;
she has no assigned deliverable.

## What I checked

The documentation draft now includes all six function references, package and
print-method help, a data dictionary, a function-selection guide, a fitted-object
reference, randomized and observational examples, lasso, prediction, and
troubleshooting. The executable vignette generates the readable tutorial;
roxygen comments generate installed R help and the browsable reference.
Open `output/help/index.html` to review the complete help system.

The following method checks were completed to support accurate help text.

I compared the implementation with the 23-page local manuscript by Quan, Iosif,
and Chen, *Doubly Robust Estimation of Individualized Treatment Effects in
Longitudinal Data* (2026 draft). Source file:
`/Users/zhikuan_quan/Downloads/Doubly_Robust_Estimation_of_Individualized_Treatment_Effects_in_Longitudinal_Data.pdf`.
It has placeholder publication metadata and is not treated as a published paper.
The older 19-page version has different equation numbering. Huling and Yu's
JSS paper remains the reference for software presentation and workflow.

| Checklist item | Manuscript location | Code and conclusion |
|---|---|---|
| Target and treatment coding | Section 2.1, pp. 3-4, Eq. (1) | `fit_lrd()`, `fit_alrd()`, `predict_ite()`: the design uses signed T times the basis divided by 2, and predictions return the full treated-minus-control contrast. Numeric 0/1 and -1/1 give the same fitted coefficients in the new tests. Baseline treatment constancy is assumed but not yet enforced. |
| Direct and augmented losses | pp. 5-6, Eqs. (3), (5), (6) | The unpenalized LRD fit and both ALRD stages match independent matrix GLS calculations for fixed covariance and positive subject-constant weights. Omitting the common 1/n factor does not change these unpenalized minimizers. |
| Main-effect versus ITE basis | Section 2.3, p. 6 | The paper permits separate classes G and F. The package uses one basis for both: a supported restriction, not the full framework. |
| Propensity weighting | Section 2.1, p. 4; Section 2.5, p. 9 | The code correctly uses p for treated subjects and 1-p for controls. `ps = NULL` assumes 0.5 in LRD but estimates scores from long-data outcome covariates in ALRD. Use explicit design probabilities, or baseline subject-level scores mapped to visits. Arbitrary visit-varying weights are not covered by the paper's scalar subject-weighted loss. |
| Working covariance | Section 2.2, p. 5; Section 2.5, p. 9 | Actual numeric time gaps agree with the manuscript. Cholesky whitening reproduces fixed-covariance GLS. Unit-diagonal working correlation is narrower than a full variance model. The package's residual-based rho updates differ from the EM-type approach discussed in the paper. |

## What I updated

I selected one bounded implementation task: covariance validity. The package
now rejects non-finite/non-numeric times, exchangeable rho at or below
-1/(m-1), duplicate AR(1) times, and negative AR(1) rho with fractional gaps.
AR(1) matrices are checked for numerical positive definiteness. Estimated
exchangeable rho respects the largest subject's lower bound, and fractional-gap
AR(1) estimation searches a nonnegative interval. Help pages and the user guide
were updated accordingly.

The regression tests also check valid positive and negative correlations,
unequal visit counts, both numeric treatment codings, inverse received-arm
probabilities, and the LRD/main-effect/ITE matrix objectives. This checks the
implementation; it is not a new simulation study or proof of consistency.

## What remains unresolved

Theorem 2.1 (pp. 6-7) states population identification under the causal assumptions
and a correct propensity or limiting main-effect function. It does not by itself
validate finite-sample lasso, restricted/misspecified ITE bases, or the full CV
pipeline. In particular, estimating the main effect by Eq. (5) with incorrect
propensity weights need not recover m merely because the chosen basis contains
m. This condition needs methodological clarification before a blanket double-
robustness claim. No change to that estimator was made here.

For example, with constant covariates, true p = 0.8, supplied p = 0.5,
m = 0 and delta = 2, the two population least-squares stages give m-hat = 0.6
and delta-hat = 1.28. This algebraic check concerns the nuisance estimator's
limit; it does not contradict the theorem's premise that the limiting m is true.

The lasso uses glmnet's transformed-data scaling/standardization, so numerical
lambda values should not be equated directly to Eq. (8) without a scaling
derivation. The quadratic spline penalty in Eq. (9) is not implemented.
Prediction factor/NA handling, zero weights, and full-pipeline resampling remain
on the existing backlog. They are not additional assignments for this week.

## Steven's documentation review checklist

1. Read each function's help as a first-time user and revise unclear wording in
   the roxygen comments in R/. Check all arguments and returned components.
2. Edit vignettes/alrd-introduction.Rmd so each example explains its inputs,
   choice of options, output and interpretation; verify installation instructions.
3. Run `Rscript scripts/build_help.R`; inspect output/help/index.html and follow
   the tutorial. Record remaining issues with the exact call/message.
4. Return the revised documentation to Zhikuan for a final methods review.

## Validation

The expanded documentation was regenerated and checked successfully. The latest
`R CMD check --no-manual` finished with **Status: OK**, including examples and
vignette rebuilding; local links resolve across all 11 generated HTML pages.
Logs for this documentation pass are in `output/review/documentation-build.log`,
`documentation-package-build.log`, and `documentation-check.log`.

All 10 tests (53 expectations) passed. The executable guide and source build
succeeded. `R CMD check --no-manual` finished with **Status: OK**, including
installed examples, help-page checks, tests, and vignette rebuilding. Logs:
`output/review/weekly-help-build.log`, `weekly-build.log`, and `weekly-check.log`.
This was not an `--as-cran` check and did not rebuild a PDF reference manual.
GitHub/ALRD was not modified. Synchronization requires separate approval.
