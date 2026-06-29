#' Predict individualized treatment effects
#'
#' @details
#' Fitted LRD and ALRD objects store coefficients for the treatment scoring
#' function `f(X)`. Prediction therefore evaluates the same basis `B(X)` used in
#' fitting and returns `B(X) %*% beta`, which estimates the time-specific ITE
#' `delta(X) = E{Y(1) - Y(-1) | X}`. The prediction does not depend on a
#' treatment column in `newdata`.
#'
#' @param object A fitted `lrd_fit` or `alrd_fit` object.
#' @param newdata A data frame containing the covariates used in the fitted
#'   treatment-effect model. If `NULL`, the training data are used.
#' @param ... Reserved for future extensions.
#'
#' @return A numeric vector of predicted individualized treatment effects.
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
