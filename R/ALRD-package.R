#' ALRD: longitudinal treatment-effect estimation
#'
#' Tools for continuous repeated outcomes with binary baseline treatment.
#' Begin with [simulate_alrd_data()], fit with [fit_lrd()] or [fit_alrd()],
#' then evaluate contrasts with [predict_ite()]. [estimate_ps()] supplies
#' logistic assignment scores and [working_covariance()] inspects correlations.
#'
#' @details Read `vignette("alrd-introduction", package = "ALRD")` for a
#'   complete workflow. This development version has known limitations in
#'   factor prediction, zero observation weights,
#'   and full-pipeline cross-validation. See the individual help pages before
#'   use. Missing model data are rejected, not automatically imputed or deleted.
#' @section Choosing a function:
#'   Use [fit_lrd()] for an unaugmented, unpenalized contrast model, or
#'   [fit_alrd()] for main-effect augmentation with optional lasso selection.
#'   [predict_ite()] evaluates the full treated-minus-control contrast rather
#'   than the fitted treatment contribution. Pass known design probabilities
#'   explicitly in randomized studies. For baseline observational treatment,
#'   use [estimate_ps()] on one row per subject and map scores to visits.
#' @seealso [fit_alrd()], [fit_lrd()], [predict_ite()], [estimate_ps()],
#'   [working_covariance()], [simulate_alrd_data()]
#' @keywords package
"_PACKAGE"
