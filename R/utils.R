.column_name <- function(expr, data, arg, value = NULL) {
  if (is.character(expr) && length(expr) == 1L) {
    name <- expr
  } else if (is.symbol(expr)) {
    name <- as.character(expr)
    if (!name %in% names(data) && is.character(value) && length(value) == 1L) {
      name <- value
    }
  } else {
    stop(sprintf("`%s` must be a column name or string.", arg), call. = FALSE)
  }

  if (!name %in% names(data)) {
    stop(sprintf("`%s` must identify a column in `data`.", arg), call. = FALSE)
  }

  name
}

.column_value <- function(expr, env) {
  tryCatch(eval(expr, env), error = function(e) NULL)
}

.expit <- function(x) {
  1 / (1 + exp(-x))
}

.check_covariates <- function(covariates, data) {
  if (is.null(covariates)) {
    return(character())
  }

  if (!is.character(covariates)) {
    stop("`covariates` must be a character vector of column names.", call. = FALSE)
  }

  missing_cols <- setdiff(covariates, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("The following covariates are missing from `data`: %s.",
              paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  covariates
}

.binary_treatment <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.factor(x)) {
    if (nlevels(x) != 2L) {
      stop("Factor treatment variables must have exactly two levels.", call. = FALSE)
    }
    return(as.integer(x == levels(x)[2L]))
  }

  if (is.character(x)) {
    levels_x <- sort(unique(x))
    if (length(levels_x) != 2L) {
      stop("Character treatment variables must have exactly two unique values.", call. = FALSE)
    }
    return(as.integer(x == levels_x[2L]))
  }

  if (is.numeric(x)) {
    unique_x <- sort(unique(stats::na.omit(x)))
    if (all(unique_x %in% c(0, 1))) {
      return(as.numeric(x))
    }
    if (all(unique_x %in% c(-1, 1))) {
      return(as.numeric(x == 1))
    }
    stop("Numeric treatment variables must be coded as 0/1 or -1/1.", call. = FALSE)
  }

  stop("Treatment must be logical, factor, character, or numeric 0/1 or -1/1.",
       call. = FALSE)
}

.signed_treatment <- function(x) {
  # The paper's treatment coding is T = 1 for treatment and T = -1 for control.
  # This helper is the single conversion point for users who provide a more
  # common 0/1 binary treatment variable.
  if (is.numeric(x)) {
    unique_x <- sort(unique(stats::na.omit(x)))
    if (all(unique_x %in% c(-1, 1))) {
      return(as.numeric(x))
    }
  }

  2 * .binary_treatment(x) - 1
}

.check_propensity <- function(ps, treatment, n) {
  if (is.null(ps)) {
    return(rep(0.5, n))
  }

  if (length(ps) == 1L) {
    ps <- rep(ps, n)
  }

  if (!is.numeric(ps) || length(ps) != n || any(!is.finite(ps)) ||
      any(ps <= 0) || any(ps >= 1)) {
    stop("`ps` must be a probability or a numeric vector of probabilities, one per row.",
         call. = FALSE)
  }

  ps
}

.propensity_for_received <- function(ps, treatment_signed) {
  # Convert pi_1(Z) into pi_T(Z): use pi_1 for treated subjects and 1 - pi_1
  # for controls. This is the denominator in Equations (3), (5), and (6).
  ifelse(treatment_signed == 1, ps, 1 - ps)
}

.validate_complete_data <- function(data, cols) {
  incomplete <- !stats::complete.cases(data[, cols, drop = FALSE])
  if (any(incomplete)) {
    stop("Rows used by the model must not contain missing values.", call. = FALSE)
  }
  invisible(TRUE)
}

.basis_formula <- function(covariates) {
  # B(X) defines the working function class for f(X) and g(X). Keeping this as
  # a formula/model.matrix path allows ordinary columns, factors, and pre-built
  # spline or interaction basis columns to be handled consistently.
  if (length(covariates) == 0L) {
    stats::as.formula("~ 1")
  } else {
    stats::reformulate(covariates)
  }
}

.fit_gls <- function(x, y, id, time, covariance, rho, weights = NULL) {
  # Minimize sum_i w_i (y_i - X_i beta)' V_i^(-1) (y_i - X_i beta).
  # Each subject block is premultiplied by chol(V_i^(-1)), converting the
  # working-covariance loss into ordinary least squares on whitened data.
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  keep <- weights > 0
  if (!all(keep)) {
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
    id <- id[keep]
    time <- time[keep]
    weights <- weights[keep]
  }

  x_out <- matrix(NA_real_, nrow = 0L, ncol = ncol(x))
  y_out <- numeric()

  for (subject in unique(id)) {
    idx <- which(id == subject)
    ord <- order(time[idx])
    idx <- idx[ord]

    # V_i(theta) is represented as a working correlation matrix here. A future
    # extension could multiply by a subject/time-specific variance model before
    # whitening, but the current paper examples focus on correlation structure.
    r <- .correlation_matrix(time[idx], covariance, rho)
    rinv_chol <- chol(solve(r))
    sw <- sqrt(weights[idx])

    # The inverse-propensity and optional user weights enter as sqrt(w) after
    # the covariance whitening, producing the same weighted quadratic objective.
    x_block <- rinv_chol %*% x[idx, , drop = FALSE]
    y_block <- as.vector(rinv_chol %*% y[idx])
    x_out <- rbind(x_out, x_block * sw)
    y_out <- c(y_out, y_block * sw)
  }

  qr_fit <- qr(x_out)
  coef <- qr.coef(qr_fit, y_out)
  names(coef) <- colnames(x)
  fitted_values <- as.vector(x %*% coef)
  residuals <- y - fitted_values

  list(
    coefficients = coef,
    fitted.values = fitted_values,
    residuals = residuals,
    rank = qr_fit$rank,
    df.residual = nrow(x_out) - qr_fit$rank,
    covariance = covariance,
    rho = rho
  )
}

.transform_gls <- function(x, y, id, time, covariance, rho, weights = NULL) {
  # Same whitening transform as `.fit_gls()`, but returns the transformed design
  # and response for penalized fitting. glmnet can then solve the lasso version
  # of the paper's weighted working-covariance objective.
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  keep <- weights > 0
  if (!all(keep)) {
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
    id <- id[keep]
    time <- time[keep]
    weights <- weights[keep]
  }

  x_out <- matrix(NA_real_, nrow = 0L, ncol = ncol(x))
  y_out <- numeric()
  id_out <- id[0L]

  for (subject in unique(id)) {
    idx <- which(id == subject)
    ord <- order(time[idx])
    idx <- idx[ord]
    r <- .correlation_matrix(time[idx], covariance, rho)
    rinv_chol <- chol(solve(r))
    sw <- sqrt(weights[idx])

    x_block <- rinv_chol %*% x[idx, , drop = FALSE]
    y_block <- as.vector(rinv_chol %*% y[idx])
    x_out <- rbind(x_out, x_block * sw)
    y_out <- c(y_out, y_block * sw)
    id_out <- c(id_out, id[idx])
  }

  list(x = x_out, y = y_out, id = id_out)
}

.cluster_fold_id <- function(id, nfolds = 5L, seed = NULL) {
  # The simulation study chooses cross-validation folds at the subject level so
  # repeated observations from one individual never appear in both training and
  # validation folds. This avoids leakage across correlated measurements.
  subjects <- unique(id)
  nfolds <- min(as.integer(nfolds), length(subjects))
  if (nfolds < 2L) {
    stop("`nfolds` must be at least 2 and no larger than the number of subjects.",
         call. = FALSE)
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  shuffled <- sample(subjects, length(subjects))
  subject_fold <- rep(seq_len(nfolds), length.out = length(shuffled))
  names(subject_fold) <- as.character(shuffled)
  as.integer(subject_fold[as.character(id)])
}

.fit_lasso_gls <- function(x, y, id, time, covariance, rho, weights = NULL,
                           nfolds = 5L, lambda = NULL, seed = NULL,
                           standardize = TRUE, refit = TRUE,
                           selection_tol = 1e-8) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("The `glmnet` package is required when `penalty = \"lasso\"`.",
         call. = FALSE)
  }

  transformed <- .transform_gls(x, y, id, time, covariance, rho, weights)
  foldid <- .cluster_fold_id(transformed$id, nfolds = nfolds, seed = seed)

  # Do not penalize the intercept in the linear ITE/main-effect basis. All other
  # columns receive the L1 penalty described in Equation (8).
  penalty_factor <- rep(1, ncol(x))
  penalty_factor[grepl("Intercept", colnames(x), fixed = TRUE)] <- 0

  # glmnet sees an ordinary gaussian lasso problem because `.transform_gls()`
  # has already encoded V_i(theta)^(-1) and 1 / pi_T(Z) in the transformed data.
  cv_fit <- glmnet::cv.glmnet(
    x = transformed$x,
    y = transformed$y,
    alpha = 1,
    family = "gaussian",
    intercept = FALSE,
    standardize = standardize,
    foldid = foldid,
    lambda = lambda,
    penalty.factor = penalty_factor
  )

  lambda_min <- cv_fit$lambda.min
  coef_lasso <- as.numeric(stats::coef(cv_fit, s = lambda_min))[-1L]
  names(coef_lasso) <- colnames(x)
  selected <- abs(coef_lasso) > selection_tol
  selected[penalty_factor == 0] <- TRUE

  # The paper's real-data analysis uses lasso for moderator selection and then
  # refits without penalty to reduce shrinkage bias in the reported ITE
  # coefficients. `refit = TRUE` implements that post-selection GLS refit.
  if (refit && any(selected)) {
    refit_fit <- .fit_gls(
      x = x[, selected, drop = FALSE],
      y = y,
      id = id,
      time = time,
      covariance = covariance,
      rho = rho,
      weights = weights
    )
    coef_final <- rep(0, ncol(x))
    names(coef_final) <- colnames(x)
    coef_final[selected] <- refit_fit$coefficients
    fitted_values <- as.vector(x %*% coef_final)
    residuals <- y - fitted_values
  } else {
    refit_fit <- NULL
    coef_final <- coef_lasso
    fitted_values <- as.vector(x %*% coef_final)
    residuals <- y - fitted_values
  }

  selected_terms <- names(coef_lasso)[abs(coef_lasso) > selection_tol]

  list(
    coefficients = coef_final,
    lasso_coefficients = coef_lasso,
    fitted.values = fitted_values,
    residuals = residuals,
    rank = sum(selected),
    df.residual = length(y) - sum(selected),
    covariance = covariance,
    rho = rho,
    penalty = "lasso",
    lambda = lambda_min,
    lambda_1se = cv_fit$lambda.1se,
    alpha = 1,
    cv_model = cv_fit,
    nfolds = length(unique(foldid)),
    selected_terms = selected_terms,
    selected_covariates = sub("^ite_", "", selected_terms),
    refit_model = refit_fit,
    refit = refit
  )
}

.clamp_rho <- function(rho) {
  min(max(rho, -0.95), 0.95)
}

.estimate_exchangeable_rho <- function(residuals, id, fallback) {
  sigma2 <- mean(residuals^2, na.rm = TRUE)
  if (!is.finite(sigma2) || sigma2 <= 0) {
    return(fallback)
  }

  numerator <- 0
  n_pairs <- 0
  for (subject in unique(id)) {
    r <- residuals[id == subject]
    if (length(r) < 2L) {
      next
    }
    pair_index <- utils::combn(seq_along(r), 2L)
    numerator <- numerator + sum(r[pair_index[1L, ]] * r[pair_index[2L, ]])
    n_pairs <- n_pairs + ncol(pair_index)
  }

  if (n_pairs == 0L) {
    return(fallback)
  }

  rho <- numerator / (n_pairs * sigma2)
  if (!is.finite(rho)) {
    fallback
  } else {
    .clamp_rho(rho)
  }
}

.estimate_ar1_rho <- function(residuals, id, time, fallback) {
  sigma2 <- mean(residuals^2, na.rm = TRUE)
  if (!is.finite(sigma2) || sigma2 <= 0) {
    return(fallback)
  }

  observed_corr <- numeric()
  lag <- numeric()

  for (subject in unique(id)) {
    idx <- which(id == subject)
    if (length(idx) < 2L) {
      next
    }
    idx <- idx[order(time[idx])]
    pair_index <- utils::combn(seq_along(idx), 2L)
    pair_lag <- abs(as.numeric(time[idx[pair_index[1L, ]]]) -
                    as.numeric(time[idx[pair_index[2L, ]]]))
    keep <- pair_lag > 0
    if (!any(keep)) {
      next
    }

    observed_corr <- c(
      observed_corr,
      residuals[idx[pair_index[1L, keep]]] *
        residuals[idx[pair_index[2L, keep]]] / sigma2
    )
    lag <- c(lag, pair_lag[keep])
  }

  if (length(observed_corr) == 0L) {
    return(fallback)
  }

  objective <- function(rho) {
    sum((observed_corr - rho^lag)^2, na.rm = TRUE)
  }

  estimated <- tryCatch(
    stats::optimize(objective, interval = c(-0.95, 0.95))$minimum,
    error = function(e) fallback
  )

  if (!is.finite(estimated)) {
    fallback
  } else {
    .clamp_rho(estimated)
  }
}

.estimate_working_rho <- function(residuals, id, time, covariance, fallback) {
  covariance <- match.arg(covariance, c("independence", "exchangeable", "ar1"))
  if (covariance == "independence") {
    return(0)
  }

  if (covariance == "exchangeable") {
    return(.estimate_exchangeable_rho(residuals, id, fallback))
  }

  .estimate_ar1_rho(residuals, id, time, fallback)
}

.rho_control <- function(covariance, rho, estimate_rho, rho_max_iter, rho_tol) {
  covariance <- match.arg(covariance, c("independence", "exchangeable", "ar1"))
  if (is.null(estimate_rho)) {
    estimate_rho <- is.null(rho) && covariance != "independence"
  }
  if (!is.logical(estimate_rho) || length(estimate_rho) != 1L || is.na(estimate_rho)) {
    stop("`estimate_rho` must be TRUE or FALSE.", call. = FALSE)
  }
  if (covariance == "independence") {
    estimate_rho <- FALSE
  }
  if (!is.numeric(rho_max_iter) || length(rho_max_iter) != 1L ||
      !is.finite(rho_max_iter) || rho_max_iter < 1L) {
    stop("`rho_max_iter` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(rho_tol) || length(rho_tol) != 1L ||
      !is.finite(rho_tol) || rho_tol <= 0) {
    stop("`rho_tol` must be a positive finite number.", call. = FALSE)
  }

  list(
    estimate = estimate_rho,
    initial = .default_rho(covariance, rho),
    max_iter = as.integer(rho_max_iter),
    tol = rho_tol
  )
}

.default_rho <- function(covariance, rho) {
  # `rho` is the working correlation parameter theta in the paper notation. It
  # can be supplied by users or used as the starting value for iterative
  # residual-based estimation.
  covariance <- match.arg(covariance, c("independence", "exchangeable", "ar1"))
  if (!is.null(rho)) {
    if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) ||
        rho <= -0.99 || rho >= 0.99) {
      stop("`rho` must be a finite number between -0.99 and 0.99.", call. = FALSE)
    }
    return(rho)
  }

  if (covariance == "independence") {
    0
  } else {
    0.3
  }
}

.correlation_matrix <- function(time, covariance, rho) {
  # Build V_i(theta) up to a marginal variance scale. AR(1) uses the observed
  # time lag, so irregular and unbalanced visit schedules are supported.
  n <- length(time)
  if (n == 1L || covariance == "independence") {
    return(diag(1, n))
  }

  if (covariance == "exchangeable") {
    out <- matrix(rho, nrow = n, ncol = n)
    diag(out) <- 1
    return(out)
  }

  distance <- abs(outer(as.numeric(time), as.numeric(time), "-"))
  rho ^ distance
}
