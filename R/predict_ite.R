#' Predict individualized treatment-effect contrasts
#'
#' Evaluate the fitted ITE basis for an LRD or ALRD object.
#' This S3 generic returns a treated-minus-control mean contrast, not a
#' potential outcome, treatment probability, or confidence interval.
#'
#' @param object Object returned by [fit_lrd()] or [fit_alrd()].
#' @param newdata Data frame with every covariate used for fitting. NULL uses
#'   the stored training data. Treatment, outcome, id and time are not required
#'   unless included among covariates.
#' @param ... Reserved for future extensions; currently unused.
#'
#' @details
#' Prediction computes B(X) * beta. Positive values favor the treated arm only
#' when larger outcomes are preferred. Negative values favor the control arm
#' under that same convention. For baseline treatment, a sign that varies across
#' visits does not define a validated treatment-switching policy.
#'
#' Supply the same precomputed basis columns as in training. Missing covariate
#' columns cause an error. Missing values may cause model.matrix() to drop rows
#' silently, so check complete.cases() and output length. Factor levels and
#' contrasts are not fully retained by fitted objects: preserve their training
#' encoding, or use explicitly prepared numeric dummy columns with a fixed basis.
#' Novel or reordered factor levels are not reliably handled.
#'
#' @return Numeric vector of ITE estimates in newdata order for complete numeric
#'   bases. With NULL, uses training-row order. Missing values and incompatible
#'   factor encoding can violate length or finiteness; see Details.
#' @seealso [fit_alrd()], [fit_lrd()]
#' @examples
#' dat <- simulate_alrd_data(n = 60, seed = 6)
#' fit <- fit_alrd(dat, id, time, treatment, outcome,
#'                 covariates = c("time", "z1", "time_z1"), ps = 0.5)
#' grid <- data.frame(time = 1:5, z1 = 0.5)
#' grid$time_z1 <- grid$time * grid$z1
#' grid$estimated_ite <- predict_ite(fit, grid)
#' grid
#' @export
predict_ite <- function(object, newdata = NULL, ...) {
  UseMethod("predict_ite")
}

#' @export
predict_ite.lrd_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    newdata <- object$model_frame
  }

  missing_cols <- setdiff(object$covariates, names(newdata))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("`newdata` is missing covariates: %s.",
              paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  # Coefficients are stored with an `ite_` prefix to distinguish them from the
  # first-stage main-effect coefficients. Remove the prefix before aligning the
  # coefficient vector with the model-matrix columns.
  basis <- stats::model.matrix(object$basis_terms, data = newdata)
  coef <- object$coefficients
  names(coef) <- sub("^ite_", "", names(coef))
  as.numeric(basis %*% coef[colnames(basis)])
}

#' @export
predict_ite.alrd_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    newdata <- object$model_frame
  }

  missing_cols <- setdiff(object$covariates, names(newdata))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("`newdata` is missing covariates: %s.",
              paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  # ALRD and LRD use the same prediction target: f(X), the estimated ITE. The
  # augmentation affects estimation efficiency and robustness, not the final
  # prediction formula.
  basis <- stats::model.matrix(object$basis_terms, data = newdata)
  coef <- object$coefficients
  names(coef) <- sub("^ite_", "", names(coef))
  as.numeric(basis %*% coef[colnames(basis)])
}
