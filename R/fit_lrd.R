#' Fit longitudinal residual difference learning
#'
#' Fits the paper's LRD working generalized least-squares model by minimizing
#' the inverse-propensity weighted longitudinal loss
#' `{Y - T f(X) / 2}' V^{-1} {Y - T f(X) / 2}`.
#'
#' @details
#' LRD-Learning corresponds to Equation (3) in the paper. It estimates the
#' treatment scoring function `f(X)` directly, without first modeling the
#' nuisance main effect `m(X)`. The design matrix used in the fitted regression
#' is `T * B(X) / 2`, where `B(X)` is built from `covariates`. The working
#' covariance `V_i(theta)` is incorporated by whitening each subject's block
#' with `V_i(theta)^(-1/2)`, and the inverse propensity contribution is handled
#' through the observation weights `1 / pi_T(Z)`.
#'
#' If `rho` is not supplied and `covariance` is `"exchangeable"` or `"ar1"`,
#' the working correlation parameter is estimated iteratively: fit the ITE
#' model, estimate `rho` from residuals, refit with the updated working
#' covariance, and repeat until convergence or `rho_max_iter` is reached.
#'
#' @param data A data frame with one row per subject-time observation.
#' @param id Subject identifier column name or bare column.
#' @param time Time column name or bare column.
#' @param treatment Binary treatment column name or bare column.
#' @param outcome Outcome column name or bare column.
#' @param covariates Character vector of covariate column names used in the
#'   treatment-effect basis.
#' @param ps Optional propensity scores for treatment `T = 1`. If `NULL`, a
#'   randomized 1:1 study is assumed and `0.5` is used for every row.
#' @param covariance Working structure. One of `"independence"`,
#'   `"exchangeable"`, or `"ar1"`.
#' @param rho Optional working correlation parameter. If supplied, it is fixed.
#'   If `NULL`, `rho` is estimated iteratively unless `estimate_rho = FALSE`.
#' @param estimate_rho Whether to estimate `rho` iteratively from residuals.
#'   Defaults to `TRUE` when `rho = NULL` and the covariance is not
#'   independence.
#' @param rho_max_iter Maximum number of residual-correlation updates.
#' @param rho_tol Convergence tolerance for absolute changes in `rho`.
#' @param weights Optional nonnegative observation weights.
#'
#' @return An object of class `lrd_fit`.
#' @export
fit_lrd <- function(data, id, time, treatment, outcome, covariates = NULL,
                    ps = NULL,
                    covariance = c("independence", "exchangeable", "ar1"),
                    rho = NULL, estimate_rho = NULL, rho_max_iter = 5L,
                    rho_tol = 1e-3, weights = NULL) {
  # Resolve columns up front. The rest of the implementation works with column
  # names, which keeps the model frame and prediction metadata explicit.
  id_expr <- substitute(id)
  time_expr <- substitute(time)
  treatment_expr <- substitute(treatment)
  outcome_expr <- substitute(outcome)
  id_name <- .column_name(id_expr, data, "id", .column_value(id_expr, parent.frame()))
  time_name <- .column_name(time_expr, data, "time", .column_value(time_expr, parent.frame()))
  treatment_name <- .column_name(
    treatment_expr,
    data,
    "treatment",
    .column_value(treatment_expr, parent.frame())
  )
  outcome_name <- .column_name(
    outcome_expr,
    data,
    "outcome",
    .column_value(outcome_expr, parent.frame())
  )
  covariates <- .check_covariates(covariates, data)
  covariance <- match.arg(covariance)
  rho_ctl <- .rho_control(covariance, rho, estimate_rho, rho_max_iter, rho_tol)

  cols <- unique(c(id_name, time_name, treatment_name, outcome_name, covariates))
  .validate_complete_data(data, cols)

  # Convert user-supplied 0/1 or -1/1 treatment to the T in {-1, 1} notation
  # used by the paper's direct-learning objective.
  treatment_signed <- .signed_treatment(data[[treatment_name]])
  ps <- .check_propensity(ps, treatment_signed, nrow(data))
  pi_t <- .propensity_for_received(ps, treatment_signed)

  # B(X) is the user-selected ITE basis. Multiplying by T / 2 turns the paper's
  # target {Y - T f(X) / 2} into a standard weighted GLS regression.
  basis_formula <- .basis_formula(covariates)
  basis <- stats::model.matrix(basis_formula, data = data)
  x <- basis * (treatment_signed / 2)
  colnames(x) <- paste0("ite_", colnames(basis))
  y <- data[[outcome_name]]

  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != nrow(data) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("`weights` must be a nonnegative numeric vector with one value per row.",
           call. = FALSE)
    }
    weights <- weights / pi_t
  } else {
    # The empirical loss weights each observed outcome by inverse probability of
    # the treatment actually received, pi_T(Z).
    weights <- 1 / pi_t
  }

  id_vector <- data[[id_name]]
  time_vector <- data[[time_name]]
  current_rho <- rho_ctl$initial
  rho_history <- current_rho
  rho_converged <- !rho_ctl$estimate

  fit_once <- function(rho_value) {
    .fit_gls(
      x = x,
      y = y,
      id = id_vector,
      time = time_vector,
      covariance = covariance,
      rho = rho_value,
      weights = weights
    )
  }

  fit <- NULL
  for (iter in seq_len(if (rho_ctl$estimate) rho_ctl$max_iter else 1L)) {
    fit <- fit_once(current_rho)
    if (!rho_ctl$estimate) {
      break
    }

    updated_rho <- .estimate_working_rho(
      residuals = fit$residuals,
      id = id_vector,
      time = time_vector,
      covariance = covariance,
      fallback = current_rho
    )
    rho_history <- c(rho_history, updated_rho)
    if (abs(updated_rho - current_rho) < rho_ctl$tol) {
      current_rho <- updated_rho
      rho_converged <- TRUE
      break
    }
    current_rho <- updated_rho
  }

  if (rho_ctl$estimate && !isTRUE(all.equal(fit$rho, current_rho))) {
    fit <- fit_once(current_rho)
  }

  out <- c(
    fit,
    list(
      call = match.call(),
      basis_formula = basis_formula,
      basis_terms = stats::terms(basis_formula, data = data),
      model_frame = data[, cols, drop = FALSE],
      propensity_scores = ps,
      propensity_for_received = pi_t,
      rho_estimation = list(
        estimated = rho_ctl$estimate,
        history = rho_history,
        converged = rho_converged,
        iterations = length(rho_history) - 1L,
        tolerance = rho_ctl$tol
      ),
      id = id_name,
      time = time_name,
      treatment = treatment_name,
      outcome = outcome_name,
      covariates = covariates
    )
  )
  class(out) <- "lrd_fit"
  out
}

#' @export
print.lrd_fit <- function(x, ...) {
  cat("Longitudinal residual difference learning fit\n")
  cat("  Working covariance:", x$covariance, "\n")
  cat("  Rho:", x$rho, "\n")
  if (!is.null(x$rho_estimation) && isTRUE(x$rho_estimation$estimated)) {
    cat("  Rho iterations:", x$rho_estimation$iterations, "\n")
    cat("  Rho converged:", x$rho_estimation$converged, "\n")
  }
  cat("  Coefficients:", length(x$coefficients), "\n")
  invisible(x)
}
