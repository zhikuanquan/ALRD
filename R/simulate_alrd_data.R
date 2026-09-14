#' Simulate continuous longitudinal treatment-effect data
#'
#' Generate repeated continuous outcomes with baseline binary treatment,
#' known propensity scores, main effects, and treatment contrasts.
#' This is a simulation helper; correspondence with the underlying method
#' paper must be verified separately before claiming exact replication.
#'
#' @param n Positive integer number of independent subjects, default 200.
#' @param times Nonempty numeric vector of distinct, finite visit times,
#'   default 1:5. All subjects use the same schedule.
#' @param scenario Integer 1 through 5. See the scenario definitions below.
#' @param main_effect_size "large" (default) or "small"; scales the coefficients
#'   entering the main effect before any squaring.
#' @param correlated Logical; whether baseline covariates arise from correlated
#'   latent normal variables. This does not control outcome error correlation.
#' @param seed Optional seed passed to set.seed(). This changes the caller's
#'   global random-number state.
#'
#' @details
#' Outcomes follow Y = true_main + treatment * true_ite / 2 +
#' subject_random_intercept + independent_error. The random intercept and error
#' each have variance one, giving conditional exchangeable outcome error
#' correlation 0.5 in every scenario, even with correlated = FALSE.
#'
#' There are ten continuous and ten binary baseline covariates. When correlated
#' is TRUE, the underlying 20-dimensional normal vector has correlation
#' 0.5^abs(j-k); binary covariates are obtained by thresholding at zero.
#' The observed binary correlations need not equal the latent correlations.
#'
#' Scenarios 1 and 2 have randomized treatment probability 0.5, with a linear
#' main effect in scenario 1 and a squared main effect in scenario 2.
#' Scenarios 3 and 4 use logistic assignment with predictor
#' z1 - z2 + z11 - z12, and respectively linear and squared main effects.
#' Scenario 5 uses randomized treatment, a squared main effect, and nonlinear
#' ITE terms z + 5*z^2 for z1, z2, z11, and z12, with signs +, -, +, -.
#' For scenarios 1--4 these moderators enter as z + time*z.
#' Each ITE also has intercept 0.5 and time slope 0.3.
#'
#' These are properties of the data-generating mechanism, not automatic labels
#' of fitted-model correctness; correctness depends on the fitted basis.
#' In particular, squared binary covariates equal the original binary covariates.
#'
#' @return Data frame with n * length(times) rows sorted by id and time:
#'   id, time, treatment (-1/1), outcome, true_main, true_ite, propensity,
#'   z1 through z20, and time_z1 through time_z20. The propensity column,
#'   not a column named ps, contains the known treated-arm probability.
#' @seealso [fit_alrd()], [fit_lrd()], [estimate_ps()]
#' @examples
#' dat <- simulate_alrd_data(n = 40, scenario = 3, seed = 9)
#' head(dat[, c("id", "time", "treatment", "outcome", "true_ite", "propensity")])
#' fit <- fit_lrd(dat, id, time, treatment, outcome,
#'                covariates = c("time", "z1", "z2"),
#'                ps = dat$propensity)
#' head(predict_ite(fit))
#' @export
simulate_alrd_data <- function(n = 200, times = 1:5,
                               scenario = 1L,
                               main_effect_size = c("large", "small"),
                               correlated = TRUE,
                               seed = NULL) {
  scenario <- as.integer(scenario)
  if (!scenario %in% 1:5) {
    stop("`scenario` must be an integer from 1 to 5.", call. = FALSE)
  }
  main_effect_size <- match.arg(main_effect_size)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  p_cont <- 10L
  p_bin <- 10L
  p <- p_cont + p_bin

  # Generate ten continuous and ten binary baseline covariates. The correlated
  # option uses an AR(1)-style covariance matrix, matching the paper's contrast
  # between independent and correlated covariate settings.
  if (correlated) {
    sigma <- 0.5 ^ abs(outer(seq_len(p), seq_len(p), "-"))
    latent <- matrix(stats::rnorm(n * p), nrow = n) %*% chol(sigma)
  } else {
    latent <- matrix(stats::rnorm(n * p), nrow = n)
  }

  z <- matrix(NA_real_, nrow = n, ncol = p)
  z[, seq_len(p_cont)] <- latent[, seq_len(p_cont)]
  z[, p_cont + seq_len(p_bin)] <- as.numeric(latent[, p_cont + seq_len(p_bin)] > 0)
  colnames(z) <- paste0("z", seq_len(p))

  ps <- if (scenario %in% c(3L, 4L)) {
    # Scenarios 3 and 4 create observational treatment assignment through a
    # logistic propensity model; the other scenarios mimic a randomized trial.
    .expit(z[, "z1"] - z[, "z2"] + z[, "z11"] - z[, "z12"])
  } else {
    rep(0.5, n)
  }
  treatment <- ifelse(stats::runif(n) < ps, 1, -1)

  beta_idx <- c(1:5, 11:15)
  # Main-effect coefficients follow Table 1's large/small main-effect settings.
  if (main_effect_size == "large") {
    beta0 <- 0.5
    beta_t <- 0.3
    beta <- c(1, -1, 1, -1, 1, -1, 1, -1, 1, -1)
  } else {
    beta0 <- 0.05
    beta_t <- 0.03
    beta <- c(0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1)
  }
  names(beta) <- paste0("z", beta_idx)

  gamma_idx <- c(1, 2, 11, 12)
  # Only z1, z2, z11, and z12 are true treatment moderators in the linear
  # scenarios; the remaining variables are prognostic covariates or distractors.
  gamma <- c(1, -1, 1, -1)
  names(gamma) <- paste0("z", gamma_idx)
  gamma0 <- 0.5
  gamma_t <- 0.3

  rows <- vector("list", length(times))
  # A subject-specific random intercept induces exchangeable within-subject
  # correlation in the simulated longitudinal outcomes.
  alpha <- stats::rnorm(n, sd = 1)

  for (k in seq_along(times)) {
    tt <- times[k]
    z_df <- as.data.frame(z)
    main_linear <- beta0 + beta_t * tt
    for (nm in names(beta)) {
      main_linear <- main_linear + beta[[nm]] * (z_df[[nm]] + tt * z_df[[nm]])
    }
    main <- if (scenario %in% c(2L, 4L, 5L)) main_linear^2 else main_linear

    ite <- gamma0 + gamma_t * tt
    for (nm in names(gamma)) {
      if (scenario == 5L) {
        # Scenario 5 deliberately violates the linear ITE working model by
        # using a nonlinear moderator effect.
        ite <- ite + gamma[[nm]] * (z_df[[nm]] + 5 * z_df[[nm]]^2)
      } else {
        ite <- ite + gamma[[nm]] * (z_df[[nm]] + tt * z_df[[nm]])
      }
    }

    y <- main + treatment * ite / 2 + alpha + stats::rnorm(n, sd = 1)
    one_time <- data.frame(
      id = seq_len(n),
      time = tt,
      treatment = treatment,
      outcome = y,
      true_main = main,
      true_ite = ite,
      propensity = ps
    )
    one_time <- cbind(one_time, z_df)
    for (q in seq_len(p)) {
      one_time[[paste0("time_z", q)]] <- tt * one_time[[paste0("z", q)]]
    }
    rows[[k]] <- one_time
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$id, out$time), , drop = FALSE]
}
