#' Build subject-specific working correlation matrices
#'
#' Construct the working correlation matrices used by LRD and ALRD fitting.
#' Despite the function name these are unit-diagonal correlations, without a
#' marginal variance model.
#'
#' @param data Data frame with complete identifiers and numeric measurement times.
#' @param id Subject identifier column name or bare column.
#' @param time Numeric time column name or bare column. Each subject's rows
#'   are sorted by this column in the returned matrix.
#' @param covariance "independence" (default), "exchangeable", or "ar1".
#' @param rho Correlation parameter strictly between -0.99 and 0.99.
#'   NULL gives zero for independence and 0.3 otherwise.
#'   This inspection function does not estimate rho.
#'
#' @details
#' Independence returns an identity matrix. Exchangeable correlation has
#' off-diagonal entries \eqn{\rho}{rho}. For a subject with \eqn{m > 1}{m > 1} measurements, positive
#' definiteness requires \eqn{-1/(m-1) < \rho < 1}{-1/(m-1) < rho < 1}.
#' Invalid values are rejected. Iterative exchangeable updates in the fitting
#' functions stay above the lower bound for the largest subject block.
#'
#' The AR(1) entries are
#' \deqn{(V_i)_{jk} = \rho^{|t_{ij}-t_{ik}|}.}{V_i[j, k] = rho^abs(time[j] - time[k]).}
#' This uses actual numeric time
#' gaps rather than visit positions. Use positive rho with fractional times:
#' negative rho with fractional gaps and repeated AR(1) times are rejected.
#' AR(1) matrices are checked for numerical positive definiteness. Numeric
#' measurement times must be finite. Changing time units changes the interpretation
#' of rho.
#'
#' @section Ordering:
#' The list follows split()'s subject-group ordering and is named by subject.
#' Matrix rows and columns correspond to ascending times within each subject;
#' the matrices do not contain row identifiers. Match by subject name, not by
#' assuming that list order equals the input order.
#'
#' @return Named list of square numeric correlation matrices, one per subject.
#'   Invalid correlation domains cause an informative error.
#' @seealso [fit_alrd()], [fit_lrd()]
#' @examples
#' visits <- data.frame(id = c("b", "a", "a", "b"),
#'                      time = c(2, 0, 2, 0))
#' matrices <- working_covariance(visits, id, time, "ar1", rho = 0.5)
#' matrices[["a"]]
#' eigen(matrices[["a"]], symmetric = TRUE)$values
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
