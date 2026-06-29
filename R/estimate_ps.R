#' Estimate propensity scores
#'
#' Fits a logistic propensity-score model for a binary treatment.
#'
#' @details
#' The LRD and ALRD objectives require `pi_1(Z) = P(T = 1 | Z)` and then use
#' `pi_T(Z)` for the treatment actually received. This helper estimates
#' `pi_1(Z)` by logistic regression when the study is observational. In a
#' randomized 1:1 study, pass `ps = 0.5` directly to [fit_lrd()] or
#' [fit_alrd()] instead of estimating a propensity model.
#'
#' @param data A data frame.
#' @param treatment Treatment column name or bare column. Numeric treatments must
#'   be coded as 0/1.
#' @param covariates Character vector of baseline or time-varying covariate
#'   column names.
#' @param trim Numeric length-two vector giving lower and upper clipping bounds
#'   for predicted propensity scores.
#'
#' @return A list with the fitted `glm` object and a numeric vector of propensity
#'   scores.
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
