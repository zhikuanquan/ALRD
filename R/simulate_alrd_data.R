#' Simulate longitudinal ITE data from the paper's simulation design
#'
#' Generates toy longitudinal data following the five simulation scenarios in
#' the ALRD manuscript. Outcomes are generated as
#' `Y_ij = m(Z_i, t_ij) + T_i * delta(Z_i, t_ij) / 2 + alpha_i + error_ij`,
#' with treatment coded as `1` for treatment and `-1` for control.
#'
#' @details
#' Scenarios 1--4 mirror the paper's robustness design: linear versus nonlinear
#' main effects crossed with randomized versus logistic propensity scores.
#' Scenario 5 uses a nonlinear ITE. The generated data also include
#' `time_z*` columns so users can fit the linear ITE specification containing
#' baseline covariates, time, and time-by-baseline interactions used in the
#' simulation study.
#'
#' @param n Number of subjects.
#' @param times Numeric vector of observation times.
#' @param scenario Integer from 1 to 5. Scenarios 1--4 vary correct versus
#'   incorrect main-effect and propensity-score specifications; scenario 5 uses
#'   a nonlinear ITE.
#' @param main_effect_size Either `"large"` or `"small"`.
#' @param correlated Whether baseline covariates should be correlated.
#' @param seed Optional random seed.
#'
#' @return A long-format data frame with outcome, treatment, true ITE,
#'   propensity score, baseline covariates `z1`--`z20`, and time interactions
#'   `time_z1`--`time_z20`.
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
