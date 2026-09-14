# Independent matrix oracle for Equations (3), (5), and (6), fixed V.
direct_subject_gls <- function(x, y, dat, weights, covariance, rho) {
  a <- matrix(0, ncol(x), ncol(x))
  b <- numeric(ncol(x))
  for (id in unique(dat$id)) {
    ii <- which(dat$id == id)
    tt <- dat$time[ii]
    v <- switch(covariance,
      independence = diag(length(ii)),
      exchangeable = {z <- matrix(rho, length(ii), length(ii)); diag(z) <- 1; z},
      ar1 = rho ^ abs(outer(tt, tt, "-")))
    precision <- solve(v) * weights[ii[1]]
    a <- a + crossprod(x[ii, , drop = FALSE], precision %*% x[ii, , drop = FALSE])
    b <- b + drop(crossprod(x[ii, , drop = FALSE], precision %*% y[ii]))
  }
  drop(solve(a, b))
}

test_that("LRD and both ALRD stages match subject-weighted matrix objectives", {
  dat <- simulate_alrd_data(n = 40, times = c(0, 0.5, 2), scenario = 3, seed = 19)
  dat <- dat[-c(2, 8, 20), ] # unbalanced visits
  basis <- model.matrix(~ time + z1 + z2, dat)
  design <- basis * dat$treatment / 2
  user_weights <- 1 + dat$id / 40
  received <- ifelse(dat$treatment == 1, dat$propensity, 1 - dat$propensity)
  weights <- user_weights / received
  for (covariance in c("independence", "exchangeable", "ar1")) {
    args <- list(data = dat, id = "id", time = "time", treatment = "treatment",
                 outcome = "outcome", covariates = c("time", "z1", "z2"),
                 ps = dat$propensity, weights = user_weights,
                 covariance = covariance, rho = 0.4)
    lrd <- do.call(fit_lrd, args)
    alrd <- do.call(fit_alrd, args)
    direct_lrd <- direct_subject_gls(design, dat$outcome, dat, weights, covariance, 0.4)
    main <- direct_subject_gls(basis, dat$outcome, dat, weights, covariance, 0.4)
    residual <- dat$outcome - drop(basis %*% main)
    delta <- direct_subject_gls(design, residual, dat, weights, covariance, 0.4)
    expect_equal(unname(lrd$coefficients), unname(direct_lrd), tolerance = 1e-9)
    expect_equal(unname(alrd$main_model$coefficients), unname(main), tolerance = 1e-9)
    expect_equal(unname(alrd$coefficients), unname(delta), tolerance = 1e-9)
    expect_equal(predict_ite(alrd), unname(drop(basis %*% delta)), tolerance = 1e-9)
    args$data$treatment <- (dat$treatment + 1) / 2
    expect_equal(do.call(fit_alrd, args)$coefficients, alrd$coefficients)
    expect_equal(do.call(fit_lrd, args)$coefficients, lrd$coefficients)
  }
})

test_that("invalid correlation domains fail with actionable errors", {
  d <- data.frame(id = rep(1, 5), time = 1:5)
  expect_error(working_covariance(d, id, time, "exchangeable", -0.5), "must exceed")
  expect_error(working_covariance(d, id, time, "exchangeable", -0.25), "must exceed")
  d$time <- c(0, 0.5, 1, 2, 3)
  expect_error(working_covariance(d, id, time, "ar1", -0.5), "integer time gaps")
  d$time[2] <- 0
  expect_error(working_covariance(d, id, time, "ar1", 0.5), "distinct")
  d$time[1] <- Inf
  expect_error(working_covariance(d, id, time, "independence"), "finite numeric")
  dat <- simulate_alrd_data(n = 20, seed = 12)
  for (f in list(fit_lrd, fit_alrd)) {
    expect_error(f(dat, id, time, treatment, outcome, "z1", ps = 0.5,
                   covariance = "exchangeable", rho = -0.5), "must exceed")
  }
})

test_that("valid correlations retain the actual time-gap definition", {
  d <- data.frame(id = rep(1, 3), time = c(0, 0.5, 2))
  expect_equal(working_covariance(d, id, time, "ar1", 0.4)[[1]],
                0.4 ^ abs(outer(d$time, d$time, "-")))
  d$time <- c(1, 2, 4)
  expect_equal(working_covariance(d, id, time, "ar1", -0.4)[[1]],
                (-0.4) ^ abs(outer(d$time, d$time, "-")))
  expect_gt(min(eigen(working_covariance(d, id, time, "exchangeable", -0.4)[[1]],
                      symmetric = TRUE)$values), 0)
})

test_that("residual updates respect subject size and fractional time domains", {
  residual <- rep(c(1, -1, 1, -1, 0), 6)
  id <- rep(1:6, each = 5)
  rho <- ALRD:::.estimate_working_rho(residual, id, rep(1:5, 6),
                                    "exchangeable", 0.3)
  expect_gt(rho, -0.25)
  v <- ALRD:::.correlation_matrix(1:5, "exchangeable", rho)
  expect_gt(min(eigen(v, symmetric = TRUE)$values), 0)
  rho_ar <- ALRD:::.estimate_working_rho(residual, id, rep(c(0, 0.5, 1, 2, 3), 6),
                                       "ar1", 0.3)
  expect_gte(rho_ar, 0)
  expect_true(all(is.finite(ALRD:::.correlation_matrix(c(0, 0.5, 1, 2, 3), "ar1", rho_ar))))
})
