#' Estimate treated-arm propensity scores
#'
#' Fit an additive logistic regression and clip its fitted probabilities.
#' The output is P(treated | covariates), not the probability of the arm received.
#'
#' @param data Complete data frame. For baseline treatment assignment, use
#'   one row per subject with baseline confounders.
#' @param treatment Treatment column name or bare column. Accepts numeric 0/1
#'   or -1/1, logical, two-level factor (second level treated), or character
#'   (alphabetically second value treated). Prefer explicit numeric coding.
#' @param covariates Nonempty character vector of existing covariate columns.
#'   An intercept is automatically included.
#' @param trim Numeric length-two vector of lower and upper clipping bounds,
#'   satisfying 0 < lower < upper < 1; default c(0.01, 0.99).
#'
#' @details
#' This function has no id argument and fits one logistic-regression row per
#' input row. Applying it directly to repeated data overweights subjects with
#' more visits. Post-treatment covariates should not be used to model baseline
#' assignment. Internally, [fit_alrd()] calls it on the full long data when ps
#' is NULL, using the same covariates as the outcome basis.
#'
#' Clipping changes probabilities without removing observations and is not a
#' remedy for a lack of overlap. Inspect score distributions by treatment and
#' the fitted model's convergence. model contains the unclipped glm fit;
#' scores contains clipped probabilities. No cross-fitting, balance diagnostic,
#' or dedicated separation handling is provided.
#'
#' @return A list with model (binomial glm), scores (one clipped probability
#'   per input row), treatment (column name), covariates (column names), and
#'   trim (the clipping bounds).
#' @seealso [fit_alrd()], [fit_lrd()]
#' @examples
#' dat <- simulate_alrd_data(n = 100, scenario = 3, seed = 10)
#' baseline <- dat[!duplicated(dat$id), ]
#' ps_fit <- estimate_ps(baseline, treatment,
#'                       covariates = c("z1", "z2", "z11", "z12"))
#' ps_rows <- ps_fit$scores[match(dat$id, baseline$id)]
#' range(ps_rows)
#' fit <- fit_alrd(dat, id, time, treatment, outcome,
#'                 covariates = c("time", "z1", "z2"), ps = ps_rows)
#' head(predict_ite(fit))
#' @export
estimate_ps <- function(data, treatment, covariates, trim = c(0.01, 0.99)) {
  treatment_expr <- substitute(treatment)
  treatment_name <- .column_name(
    treatment_expr,
    data,
    "treatment",
    .column_value(treatment_expr, parent.frame())
  )
  covariates <- .check_covariates(covariates, data)

  if (length(covariates) == 0L) {
    stop("`covariates` must contain at least one column for propensity estimation.",
         call. = FALSE)
  }

  if (!is.numeric(trim) || length(trim) != 2L || any(!is.finite(trim)) ||
      trim[1L] <= 0 || trim[2L] >= 1 || trim[1L] >= trim[2L]) {
    stop("`trim` must be two finite probabilities with 0 < lower < upper < 1.",
         call. = FALSE)
  }

  data[[treatment_name]] <- .binary_treatment(data[[treatment_name]])
  .validate_complete_data(data, c(treatment_name, covariates))

  formula <- stats::reformulate(covariates, response = treatment_name)
  fit <- stats::glm(formula, data = data, family = stats::binomial())
  scores <- as.numeric(stats::predict(fit, type = "response"))
  scores <- pmin(pmax(scores, trim[1L]), trim[2L])

  list(
    model = fit,
    scores = scores,
    treatment = treatment_name,
    covariates = covariates,
    trim = trim
  )
}
