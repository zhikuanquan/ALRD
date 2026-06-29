test_that("fit_alrd returns finite ITE predictions", {
  set.seed(42)
  n_id <- 60
  n_time <- 3
  dat <- data.frame(
    id = rep(seq_len(n_id), each = n_time),
    time = rep(seq_len(n_time), times = n_id)
  )
  dat$x1 <- rnorm(nrow(dat))
  dat$x2 <- rbinom(nrow(dat), 1, 0.4)
  dat$a <- rbinom(nrow(dat), 1, plogis(-0.1 + 0.5 * dat$x1))
  dat$y <- 1 + 0.4 * dat$x1 - 0.3 * dat$x2 +
    dat$a * (0.6 + 0.2 * dat$x1) + rnorm(nrow(dat), sd = 0.5)

  fit <- fit_alrd(
    data = dat,
    id = "id",
    time = "time",
    treatment = "a",
    outcome = "y",
    covariates = c("x1", "x2"),
    ps = 0.5,
    covariance = "exchangeable"
  )

  pred <- predict_ite(fit, dat)
  expect_s3_class(fit, "alrd_fit")
  expect_length(pred, nrow(dat))
  expect_true(all(is.finite(pred)))
})

test_that("fit_lrd returns finite ITE predictions", {
  set.seed(123)
  dat <- data.frame(
    id = rep(1:30, each = 2),
    time = rep(1:2, times = 30),
    x = rnorm(60)
  )
  dat$a <- rbinom(nrow(dat), 1, 0.5)
  dat$y <- 2 + dat$x + dat$a * (0.5 + dat$x) + rnorm(nrow(dat))

  fit <- fit_lrd(
    data = dat,
    id = "id",
    time = "time",
    treatment = "a",
    outcome = "y",
    covariates = "x",
    ps = 0.5
  )

  pred <- predict_ite(fit, dat)
  expect_s3_class(fit, "lrd_fit")
  expect_length(pred, nrow(dat))
  expect_true(all(is.finite(pred)))
})

test_that("simulate_alrd_data creates paper-style longitudinal data", {
  dat <- simulate_alrd_data(n = 25, scenario = 5, seed = 100)

  expect_equal(nrow(dat), 25 * 5)
  expect_true(all(dat$treatment %in% c(-1, 1)))
  expect_true(all(c("true_main", "true_ite", "propensity", "time_z1") %in% names(dat)))
  expect_true(all(dat$propensity > 0 & dat$propensity < 1))
})

test_that("lasso ALRD exposes selected terms and tuning output", {
  dat <- simulate_alrd_data(n = 45, scenario = 1, seed = 25)
  covs <- c("time", paste0("z", 1:4), paste0("z", 11:12),
            paste0("time_z", c(1:4, 11:12)))

  fit <- fit_alrd(
    data = dat,
    id = "id",
    time = "time",
    treatment = "treatment",
    outcome = "outcome",
    covariates = covs,
    ps = 0.5,
    penalty = "lasso",
    nfolds = 3,
    seed = 3,
    covariance = "exchangeable"
  )

  expect_s3_class(fit, "alrd_fit")
  expect_equal(fit$penalty, "lasso")
  expect_true(is.finite(fit$lambda))
  expect_true(is.numeric(fit$lasso_coefficients))
  expect_true(is.character(fit$selected_covariates))
  expect_length(predict_ite(fit, dat), nrow(dat))
})

test_that("ALRD estimates rho when rho is omitted", {
  dat <- simulate_alrd_data(n = 35, scenario = 1, seed = 222)
  covs <- c("time", "z1", "z2", "z11", "z12")

  fit <- fit_alrd(
    data = dat,
    id = "id",
    time = "time",
    treatment = "treatment",
    outcome = "outcome",
    covariates = covs,
    ps = 0.5,
    penalty = "none",
    covariance = "exchangeable",
    rho_max_iter = 3
  )

  expect_true(fit$rho_estimation$estimated)
  expect_true(is.finite(fit$rho))
  expect_true(length(fit$rho_estimation$history) >= 2)
})

test_that("supplied rho is fixed", {
  dat <- simulate_alrd_data(n = 30, scenario = 1, seed = 333)
  covs <- c("time", "z1", "z2")

  fit <- fit_lrd(
    data = dat,
    id = "id",
    time = "time",
    treatment = "treatment",
    outcome = "outcome",
    covariates = covs,
    ps = 0.5,
    covariance = "exchangeable",
    rho = 0.5
  )

  expect_false(fit$rho_estimation$estimated)
  expect_equal(fit$rho, 0.5)
})
