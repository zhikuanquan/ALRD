#' Fit longitudinal residual difference learning
#'
#' Estimate a continuous longitudinal treatment contrast by inverse-propensity
#' weighted generalized least squares without a main-effect augmentation stage.
#'
#' @inheritParams fit_alrd
#' @param covariates Character vector of existing covariate column names.
#'   NULL fits an intercept-only treatment contrast. An intercept is always added.
#' @param ps Treated-arm probability, as a scalar or one probability per input
#'   row, strictly between zero and one. NULL uses 0.5 for every row, unlike
#'   [fit_alrd()], which estimates scores when ps is NULL.
#'
#' @details
#' The design is \eqn{T_i B_i/2}{T_i * B_i / 2}, with signed treatment
#' \eqn{T_i \in \{-1,1\}}{T_i in {-1, 1}} and subject-specific basis matrix
#' \eqn{B_i}{B_i}. For weights constant within subject, the criterion is
#' \deqn{\sum_{i=1}^{n} w_i \left(\mathbf{Y}_i - \frac{T_i}{2}B_i\beta\right)^{\mathsf{T}} V_i^{-1} \left(\mathbf{Y}_i - \frac{T_i}{2}B_i\beta\right).}{sum_i w_i * (Y_i - T_i * B_i * beta / 2)' V_i^(-1) * (Y_i - T_i * B_i * beta / 2).}
#' Here \eqn{w_i = a_i/\pi_{T_i}(Z_i)}{w_i = a_i / pi_T_i(Z_i)} combines
#' the optional analysis weight \eqn{a_i}{a_i} (one by default) and the inverse
#' probability of the treatment received; \eqn{V_i}{V_i} is the working correlation.
#' There is no penalized LRD option in this release.
#' Treatment coding, covariance restrictions, weight limitations, and prediction
#' limitations are described in [fit_alrd()] and [working_covariance()].
#'
#' For a baseline observational treatment, estimate propensity scores on one
#' row per subject using baseline confounders, then map the scores to all visits.
#' Automatic validation of constant within-subject treatment is not implemented.
#'
#' @return An lrd_fit list containing coefficients (with ite_ prefixes),
#'   fitted.values (\eqn{T_i\widehat{\delta}(X_{ij})/2}{T_i * estimated ITE / 2}), residuals (Y minus fitted.values),
#'   rank, df.residual, covariance, rho, propensity_scores,
#'   propensity_for_received, and rho_estimation. The latter contains estimated,
#'   history, converged, iterations, and tolerance. Training metadata are stored
#'   in model_frame, basis_formula, basis_terms, id, time, treatment, outcome,
#'   and covariates. No standard errors or confidence intervals are computed.
#' @seealso [fit_alrd()], [predict_ite()], [working_covariance()]
#' @examples
#' dat <- simulate_alrd_data(n = 60, seed = 8)
#' fit <- fit_lrd(dat, "id", "time", "treatment", "outcome",
#'                covariates = c("time", "z1", "time_z1"), ps = 0.5)
#' print(fit)
#' head(predict_ite(fit, dat))
#' constant <- fit_lrd(dat, id, time, treatment, outcome, ps = 0.5)
#' constant$coefficients
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

#' @rdname print.alrd_fit
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
