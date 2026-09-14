#' Fit augmented longitudinal residual difference learning
#'
#' Estimate a continuous longitudinal treatment contrast in two stages:
#' fit a main-effect model, then regress the residualized outcome on the
#' signed treatment multiplied by a covariate basis divided by two.
#'
#' @param data Data frame with one row per subject and visit. Use complete,
#'   finite outcomes and covariates, unique subject-time pairs, and a baseline
#'   treatment that is constant within subject. These requirements are not all
#'   enforced by the current implementation.
#' @param id Subject identifier, as a bare column name or a single string.
#' @param time Numeric measurement-time column, as a bare name or string.
#'   For AR(1), the numeric time difference determines the correlation.
#' @param treatment Binary treatment column, as a bare name or string.
#'   Numeric -1/1 or 0/1 and logical values are accepted. For factors the second
#'   level is treated; for character values the alphabetically second value is
#'   treated. Explicit numeric coding is recommended.
#' @param outcome Continuous numeric outcome column, as a bare name or string.
#' @param covariates Character vector of existing column names, with at least
#'   one entry. An intercept is added. The same basis is used for both stages.
#'   Construct interactions or nonlinear basis columns before fitting.
#' @param ps Probability of the treated arm, either a scalar or a numeric
#'   vector in input-row order, strictly between zero and one. For baseline
#'   assignment repeat the subject-level probability across visits.
#'   NULL estimates scores from all rows and all covariates via [estimate_ps()].
#'   For a randomized study supply its known assignment probability explicitly.
#' @param penalty "none" (default) or "lasso". Lasso is applied to both stages
#'   and requires the suggested package glmnet.
#' @param refit Logical; for lasso, refit selected ITE columns without penalty.
#'   The main-effect stage is not refitted. Ignored for unpenalized fitting.
#' @param nfolds Requested number of subject-level folds; default five.
#'   Use at least three subjects and three folds. The implementation caps the
#'   number at the observed subject count; two folds fail in cv.glmnet().
#' @param lambda Optional lambda sequence passed separately to each stage's
#'   cv.glmnet() call. Each stage selects its own lambda.min; this is not a
#'   fixed-lambda argument.
#' @param seed Optional seed for subject-fold assignment. The fold helper
#'   restores the caller's random-number state.
#' @param standardize Logical passed to glmnet for the transformed design.
#'   Defaults to TRUE. Ignored for unpenalized fitting.
#' @param covariance "independence" (default), "exchangeable", or "ar1".
#'   See [working_covariance()] for mathematical definitions and restrictions.
#' @param rho Optional correlation parameter. Supplying rho fixes it by
#'   default; with estimate_rho = TRUE it is instead the starting value.
#'   If omitted, the initial value is zero for independence and 0.3 otherwise.
#' @param estimate_rho Logical. Default: update rho only when rho is omitted
#'   and covariance is not independence. Independence always disables updating.
#' @param rho_max_iter Maximum number of correlation updates, default five.
#' @param rho_tol Positive tolerance for the absolute change in rho.
#' @param weights Optional numeric vector in input-row order, multiplied by
#'   inverse received-treatment probabilities. Although zeros and within-subject
#'   variation are accepted, their handling has unresolved limitations; use
#'   strictly positive weights constant within subject for the documented loss.
#'
#' @details
#' With signed treatment T, the model is Y = m(X) + T * delta(X) / 2 + error.
#' Both stages use an intercept-containing model matrix. A subject block is
#' whitened with C = chol(solve(V)); square-root row weights are applied
#' after whitening. The scalar subject-weighted loss is recovered when weights
#' are constant within subject. Arbitrary row weights instead give a quadratic
#' form C' W C and require a separate statistical specification.
#'
#' Lasso leaves basis columns whose names contain "Intercept" unpenalized
#' and selects lambda.min. The saved lambda_1se is informational, not the
#' selected tuning rule. Selected terms refer to design columns, which may
#' include an intercept or individual factor contrasts.
#'
#' @section Validation limits:
#' Subject-fold assignment keeps repeated observations together within each
#' cv.glmnet() call. However, main-effect residualization, any internally
#' estimated propensities, and working correlation use the full fitting sample
#' before ITE-stage validation. This is not full-pipeline out-of-sample validation.
#' Use a separate subject-level test sample for performance assessment.
#' There are no dedicated interval or subgroup-validation methods.
#'
#' Zero observation weights can shorten residual vectors without updating
#' stored row metadata. Rank deficiency can yield NA coefficients.
#' Factor levels and contrasts are not fully preserved for new-data prediction.
#' Inspect rho_estimation$converged: nonconvergence does not issue a warning.
#'
#' @return An alrd_fit list. Important components are:
#'   \describe{
#'     \item{coefficients}{Named ITE-basis coefficients, with an ite_ prefix;
#'       after refitting if requested.}
#'     \item{fitted.values}{Fitted treatment contribution T * delta / 2 to the
#'       residualized response, not the full fitted outcome or the ITE itself.}
#'     \item{residuals}{Outcome minus main effect and treatment contribution.}
#'     \item{main_model, main_effect}{First-stage fit and predictions.}
#'     \item{ps_model}{Result of estimate_ps(), or NULL for supplied scores.}
#'     \item{propensity_scores, propensity_for_received}{Treated-arm scores
#'       and received-arm probabilities.}
#'     \item{rho, covariance, rho_estimation}{Working correlation and structure;
#'       estimation flag, history, convergence, iterations, and tolerance.}
#'     \item{selected_terms, selected_covariates}{Selected design-column names,
#'       with or without the ite_ prefix. All columns for unpenalized fitting.}
#'     \item{lambda, lambda_1se, cv_model}{Lasso-only tuning results for the
#'       ITE stage. Main-stage results are in main_model.}
#'     \item{lasso_coefficients, refit_model, refit}{Lasso-only penalized
#'       coefficients, optional refit, and refit flag.}
#'     \item{model_frame, basis_terms, basis_formula}{Training data and basis
#'       information used for prediction.}
#'   }
#' @seealso [fit_lrd()], [predict_ite()], [estimate_ps()], [working_covariance()]
#' @examples
#' dat <- simulate_alrd_data(n = 60, seed = 4)
#' covs <- c("time", "z1", "z2", "time_z1", "time_z2")
#' fit <- fit_alrd(dat, id, time, treatment, outcome, covs,
#'                 ps = 0.5, covariance = "exchangeable", rho = 0.5)
#' print(fit)
#' head(predict_ite(fit))
#' fit$rho_estimation
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   sparse <- fit_alrd(dat, id, time, treatment, outcome, covs,
#'                     ps = 0.5, penalty = "lasso", nfolds = 3, seed = 4)
#'   sparse$selected_terms
#' }
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

#' Print a longitudinal treatment-effect fit
#'
#' Display the working covariance, correlation parameter, convergence status
#' when estimated, and coefficient count. ALRD also displays the penalty and,
#' for lasso, the selected lambda and selected-term count.
#' @param x An object returned by [fit_alrd()] or [fit_lrd()].
#' @param ... Reserved for future extensions; currently unused.
#' @return The fitted object, invisibly.
#' @details This is a compact fit description, not an inferential summary.
#'   It does not report confidence intervals or validate treatment subgroups.
#' @name print.alrd_fit
#' @seealso [predict_ite()]
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
