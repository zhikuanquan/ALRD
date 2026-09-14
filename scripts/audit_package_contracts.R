# Reproducible diagnostic checks of the current implementation.
# Run from Documents/ALRD: Rscript scripts/audit_package_contracts.R
# These checks record defects; they are not a passing regression-test suite.
pkgload::load_all(quiet = TRUE)
set.seed(120)
results <- list()
record <- function(name, expr) {
  ans <- tryCatch(as.character(force(expr)), error = function(e) {
    paste("ERROR:", conditionMessage(e))
  })
  results[[name]] <<- ans
  cat(name, ": ", paste(ans, collapse = "; "), "\n", sep = "")
}
dat <- simulate_alrd_data(n = 40, seed = 12)
fit_args <- list(data = dat, id = "id", time = "time",
                 treatment = "treatment", outcome = "outcome",
                 covariates = c("time", "z1", "z2"), ps = 0.5)
record("zero_weight_residual_length", {
  f <- do.call(fit_alrd, c(fit_args, list(
    weights = replace(rep(1, nrow(dat)), 1, 0))))
  paste("rows", nrow(dat), "residuals", length(f$residuals))
})
record("exchangeable_minimum_eigenvalue", {
  v <- working_covariance(dat, id, time, "exchangeable", rho = -0.5)[[1]]
  min(eigen(v, symmetric = TRUE)$values)
})
record("fractional_time_negative_rho_nonfinite", {
  v <- working_covariance(data.frame(id = c(1, 1), time = c(0, 0.5)),
                          id, time, "ar1", rho = -0.5)[[1]]
  any(!is.finite(v))
})
record("two_fold_lasso", {
  do.call(fit_alrd, c(fit_args, list(penalty = "lasso", nfolds = 2)))
  "completed"
})
dat$group <- factor(rep(c("a", "b", "c", "a", "b"), 40))
factor_fit <- fit_lrd(dat, id, time, treatment, outcome, "group", ps = 0.5)
record("reordered_factor_nonfinite_predictions", {
  nd <- dat
  nd$group <- factor(as.character(nd$group), levels = c("c", "b", "a"))
  sum(!is.finite(predict_ite(factor_fit, nd)))
})
record("missing_prediction_row_count", {
  nd <- dat
  nd$group[1] <- NA
  paste("input", nrow(nd), "predictions", length(predict_ite(factor_fit, nd)))
})
record("rank_deficient_nonfinite_predictions", {
  dat$zdup <- dat$z1
  f <- fit_lrd(dat, id, time, treatment, outcome, c("z1", "zdup"), ps = 0.5)
  sum(!is.finite(predict_ite(f, dat)))
})
record("supplied_rho_can_be_updated", {
  f <- do.call(fit_alrd, c(fit_args, list(covariance = "exchangeable",
                                       rho = 0.5, estimate_rho = TRUE)))
  paste(f$rho_estimation$history, collapse = ", ")
})
record("varying_weight_row_order_difference", {
  w <- seq(0.5, 2, length.out = nrow(dat))
  a <- fit_lrd(dat, id, time, treatment, outcome, c("time", "z1"),
               ps = 0.5, weights = w, covariance = "exchangeable", rho = 0.5)
  ix <- rev(seq_len(nrow(dat)))
  b <- fit_lrd(dat[ix, ], id, time, treatment, outcome, c("time", "z1"),
               ps = 0.5, weights = w[ix], covariance = "exchangeable", rho = 0.5)
  max(abs(a$coefficients - b$coefficients))
})
record("varying_weight_quadratic_form_difference", {
  v <- matrix(0.5, 3, 3); diag(v) <- 1
  cmat <- chol(solve(v)); w <- diag(c(1, 2, 4)); e <- c(1, -1, 2)
  actual <- drop(t(e) %*% t(cmat) %*% w %*% cmat %*% e)
  symmetric <- drop(t(e) %*% sqrt(w) %*% solve(v) %*% sqrt(w) %*% e)
  paste("implemented", signif(actual, 7), "symmetric", signif(symmetric, 7))
})
record("constant_subject_weight_direct_gls_difference", {
  x <- cbind(1, dat$z1) * dat$treatment / 2
  cross <- matrix(0, 2, 2); rhs <- numeric(2)
  w <- 1 + dat$id / 40
  for (s in unique(dat$id)) {
    ix <- which(dat$id == s)
    v <- matrix(0.5, length(ix), length(ix)); diag(v) <- 1
    a <- solve(v) * (w[ix[1]] / 0.5)
    cross <- cross + t(x[ix, ]) %*% a %*% x[ix, ]
    rhs <- rhs + drop(t(x[ix, ]) %*% a %*% dat$outcome[ix])
  }
  f <- fit_lrd(dat, id, time, treatment, outcome, "z1", ps = 0.5,
               weights = w, covariance = "exchangeable", rho = 0.5)
  max(abs(drop(solve(cross, rhs)) - f$coefficients))
})
cat("\n")
print(sessionInfo())
