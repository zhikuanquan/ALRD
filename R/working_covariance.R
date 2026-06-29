#' Build working covariance matrices
#'
#' Constructs subject-level working correlation matrices for longitudinal data.
#'
#' @details
#' These matrices represent the `V_i(theta)` component in the paper's
#' longitudinal losses. They are returned mainly for inspection and diagnostics;
#' the fitting functions build the same matrices internally and use
#' `V_i(theta)^(-1/2)` to whiten each subject's contribution.
#'
#' @param data A data frame.
#' @param id Subject identifier column name or bare column.
#' @param time Time column name or bare column.
#' @param covariance Working structure. One of `"independence"`,
#'   `"exchangeable"`, or `"ar1"`.
#' @param rho Correlation parameter. If `NULL`, `0` is used for independence and
#'   `0.3` for exchangeable or AR(1).
#'
#' @return A named list of working correlation matrices, one per subject.
#' @export
working_covariance <- function(data, id, time,
                               covariance = c("independence", "exchangeable", "ar1"),
                               rho = NULL) {
  id_expr <- substitute(id)
  time_expr <- substitute(time)
  id_name <- .column_name(id_expr, data, "id", .column_value(id_expr, parent.frame()))
  time_name <- .column_name(time_expr, data, "time", .column_value(time_expr, parent.frame()))
  covariance <- match.arg(covariance)
  rho <- .default_rho(covariance, rho)
  .validate_complete_data(data, c(id_name, time_name))

  out <- lapply(split(seq_len(nrow(data)), data[[id_name]]), function(idx) {
    idx <- idx[order(data[[time_name]][idx])]
    .correlation_matrix(data[[time_name]][idx], covariance, rho)
  })

  names(out) <- names(split(seq_len(nrow(data)), data[[id_name]]))
  out
}
