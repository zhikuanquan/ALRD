#' Fit augmented longitudinal residual difference learning
#'
#' Fits ALRD using the paper's two-stage longitudinal objective. Stage 1
#' estimates the main effect `m(X)` with inverse-propensity weighted working
#' GLS. Stage 2 estimates the ITE by fitting `Y - mhat(X)` to `T f(X) / 2`
#' under the same working covariance. A lasso version is available for
#' treatment-moderator selection, with optional unpenalized refit on selected
#' terms.
#'
#' @details
#' The paper decomposes the longitudinal potential-outcome contrast as
#' `Y_ij = m(X_ij) + T_i * delta(X_ij) / 2 + epsilon_ij`, where `T_i` is coded
#' as `1` for treatment and `-1` for control. ALRD estimates the ITE function
#' `delta(X)` through the treatment scoring function `f(X)`.
#'
#' Stage 1 estimates the augmentation function `m(X)` by minimizing the
#' inverse-propensity weighted working-covariance loss in Equation (5) of the
#' paper. Stage 2 minimizes the augmented loss in Equation (6), using
#' `Y - mhat(X)` as the outcome and `T * B(X) / 2` as the treatment-effect
#' design matrix, where `B(X)` is the user-specified covariate basis.
#'
#' When `penalty = "lasso"`, the same GLS whitening and inverse-propensity
#' weighting are applied before `glmnet` is called. Cross-validation folds are
#' assigned at the subject level, matching the simulation and data-analysis
#' strategy in the paper. If `refit = TRUE`, the selected ITE terms are refit
#' without penalty to reduce lasso shrinkage bias.
#'
#' If `rho` is not supplied and `covariance` is `"exchangeable"` or `"ar1"`,
#' ALRD estimates the working correlation iteratively. Each iteration fits the
#' main-effect augmentation and ITE model at the current `rho`, estimates a new
#' `rho` from the final residuals `Y - mhat(X) - T*fhat(X)/2`, then refits until
#' the change in `rho` is smaller than `rho_tol` or `rho_max_iter` is reached.
#'
#' @param data A data frame with one row per subject-time observation.
#' @param id Subject identifier column name or bare column.
#' @param time Time column name or bare column.
#' @param treatment Binary treatment column name or bare column.
#' @param outcome Outcome column name or bare column.
#' @param covariates Character vector of covariate column names used in both the
#'   main-effect and treatment-effect models.
#' @param ps Optional propensity scores for treatment `T = 1`. If `NULL`, they
#'   are estimated with [estimate_ps()].
#' @param penalty Either `"none"` for non-penalized ALRD or `"lasso"` for
#'   lasso-penalized main-effect and ITE stages.
#' @param refit For lasso fits, whether to refit the ITE stage without penalty
#'   using the selected treatment-effect terms.
#' @param nfolds Number of clustered cross-validation folds for lasso.
#' @param lambda Optional lambda sequence passed to `glmnet::cv.glmnet()`.
#' @param seed Optional seed for clustered cross-validation folds.
#' @param standardize Whether `glmnet` should standardize columns for lasso.
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
#' @return An object of class `alrd_fit`.
#' @export
fit_alrd <- function(data, id, time, treatment, outcome, covariates,
                     ps = NULL,
                     penalty = c("none", "lasso"),
                     refit = TRUE,
                     nfolds = 5L,
                     lambda = NULL,
                     seed = NULL,
                     standardize = TRUE,
                     covariance = c("independence", "exchangeable", "ar1"),
                     rho = NULL, estimate_rho = NULL, rho_max_iter = 5L,
                     rho_tol = 1e-3, weights = NULL) {
  # Resolve the user-facing column specifications once. Public functions accept
  # both bare column names and character strings so examples can be written in
  # either base-R or package-style notation.
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
  penalty <- match.arg(penalty)
  rho_ctl <- .rho_control(covariance, rho, estimate_rho, rho_max_iter, rho_tol)

  if (length(covariates) == 0L) {
    stop("`covariates` must contain at least one column for ALRD.", call. = FALSE)
  }

  cols <- unique(c(id_name, time_name, treatment_name, outcome_name, covariates))
  .validate_complete_data(data, cols)

  # The manuscript uses T in {-1, 1}. Users may supply treatment as 0/1, but the
  # learning objective is always formed on the signed treatment scale.
  treatment_signed <- .signed_treatment(data[[treatment_name]])

  if (is.null(ps)) {
    # Observational-data workflow: estimate pi_1(Z) when it is not supplied.
    # For randomized studies, users can pass ps = 0.5 to use the known design.
    ps_fit <- estimate_ps(data, treatment = treatment_name, covariates = covariates)
    ps <- ps_fit$scores
  } else {
    ps_fit <- NULL
    ps <- .check_propensity(ps, treatment_signed, nrow(data))
  }

  # Equation (6) weights each subject contribution by 1 / pi_T(Z), where pi_T is
  # pi_1 for treated subjects and 1 - pi_1 for controls.
  pi_t <- .propensity_for_received(ps, treatment_signed)

  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != nrow(data) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("`weights` must be a nonnegative numeric vector with one value per row.",
           call. = FALSE)
    }
    analysis_weights <- weights / pi_t
  } else {
    analysis_weights <- 1 / pi_t
  }

  # B(X) is the working function class for both the main effect g(X) and the ITE
  # f(X). Users can include time, baseline covariates, time interactions, spline
  # basis columns, or other precomputed basis terms in `covariates`.
  basis_formula <- .basis_formula(covariates)
  basis <- stats::model.matrix(basis_formula, data = data)

  # Stage 2: estimate f(X) = delta(X). Since the paper's loss contains
  # T_i f(X_i) / 2, the regression design is B(X) multiplied by T_i / 2 and the
  # response is the augmented outcome Y - mhat(X).
  x <- basis * (treatment_signed / 2)
  colnames(x) <- paste0("ite_", colnames(basis))

  id_vector <- data[[id_name]]
  time_vector <- data[[time_name]]
  y <- data[[outcome_name]]

  fit_once <- function(rho_value) {
    # Stage 1: estimate m(X). This is Equation (5) in the paper. The helper
    # functions implement V_i(theta)^(-1) through a GLS whitening transform.
    if (penalty == "lasso") {
      main_fit <- .fit_lasso_gls(
        x = basis,
        y = y,
        id = id_vector,
        time = time_vector,
        covariance = covariance,
        rho = rho_value,
        weights = analysis_weights,
        nfolds = nfolds,
        lambda = lambda,
        seed = seed,
        standardize = standardize,
        refit = FALSE
      )
    } else {
      main_fit <- .fit_gls(
        x = basis,
        y = y,
        id = id_vector,
        time = time_vector,
        covariance = covariance,
        rho = rho_value,
        weights = analysis_weights
      )
      main_fit$penalty <- "none"
    }

    main_effect <- as.numeric(basis %*% main_fit$coefficients)
    residual <- y - main_effect

    # Penalized ALRD performs lasso selection on the treatment-effect basis only.
    # With `refit = TRUE`, the selected moderators are subsequently refit with
    # the non-penalized working GLS objective, matching the paper's real-data
    # analysis where lasso identifies moderators and refitting reduces bias.
    if (penalty == "lasso") {
      fit <- .fit_lasso_gls(
        x = x,
        y = residual,
        id = id_vector,
        time = time_vector,
        covariance = covariance,
        rho = rho_value,
        weights = analysis_weights,
        nfolds = nfolds,
        lambda = lambda,
        seed = seed,
        standardize = standardize,
        refit = refit
      )
    } else {
      fit <- .fit_gls(
        x = x,
        y = residual,
        id = id_vector,
        time = time_vector,
        covariance = covariance,
        rho = rho_value,
        weights = analysis_weights
      )
      fit$penalty <- "none"
      fit$selected_terms <- names(fit$coefficients)
      fit$selected_covariates <- sub("^ite_", "", names(fit$coefficients))
    }

    list(
      fit = fit,
      main_fit = main_fit,
      main_effect = main_effect
    )
  }

  current_rho <- rho_ctl$initial
  rho_history <- current_rho
  rho_converged <- !rho_ctl$estimate
  fit_bundle <- NULL

  for (iter in seq_len(if (rho_ctl$estimate) rho_ctl$max_iter else 1L)) {
    fit_bundle <- fit_once(current_rho)
    if (!rho_ctl$estimate) {
      break
    }

    updated_rho <- .estimate_working_rho(
      residuals = fit_bundle$fit$residuals,
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

  if (rho_ctl$estimate && !isTRUE(all.equal(fit_bundle$fit$rho, current_rho))) {
    fit_bundle <- fit_once(current_rho)
  }

  fit <- fit_bundle$fit
  main_fit <- fit_bundle$main_fit
  main_effect <- fit_bundle$main_effect

  out <- c(
    fit,
    list(
      call = match.call(),
      main_model = main_fit,
      main_effect = main_effect,
      ps_model = ps_fit,
      propensity_scores = ps,
      propensity_for_received = pi_t,
      rho_estimation = list(
        estimated = rho_ctl$estimate,
        history = rho_history,
        converged = rho_converged,
        iterations = length(rho_history) - 1L,
        tolerance = rho_ctl$tol
      ),
      basis_formula = basis_formula,
      basis_terms = stats::terms(basis_formula, data = data),
      model_frame = data[, cols, drop = FALSE],
      id = id_name,
      time = time_name,
      treatment = treatment_name,
      outcome = outcome_name,
      covariates = covariates
    )
  )
  class(out) <- "alrd_fit"
  out
}

#' @export
print.alrd_fit <- function(x, ...) {
  cat("Augmented longitudinal residual difference learning fit\n")
  cat("  Working covariance:", x$covariance, "\n")
  cat("  Rho:", x$rho, "\n")
  if (!is.null(x$rho_estimation) && isTRUE(x$rho_estimation$estimated)) {
    cat("  Rho iterations:", x$rho_estimation$iterations, "\n")
    cat("  Rho converged:", x$rho_estimation$converged, "\n")
  }
  cat("  Penalty:", x$penalty, "\n")
  if (identical(x$penalty, "lasso")) {
    cat("  Lambda:", signif(x$lambda, 4), "\n")
    cat("  Selected ITE terms:", length(x$selected_terms), "\n")
  }
  cat("  Coefficients:", length(x$coefficients), "\n")
  invisible(x)
}
