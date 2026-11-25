
library(glmnet)

DTLASSO_METHOD_VERSION <- "tipping_common_cv_q60_min100_trim15_gt_20260929"
DTLASSO_THRESHOLD_RULE <- "Q > tau"

# Model: beta0 + X beta + (delta0 + X delta) I(Q > tau).
# beta is the low-Q (Q <= tau) slope; delta is high-Q minus low-Q.
# Equality belongs to the baseline group in every fitting/testing path.
threshold_indicator_gt <- function(q, tau) as.numeric(q > tau)

assert_gt_threshold_object <- function(x, label = "Cached threshold object") {
  if (!identical(x$threshold_rule, DTLASSO_THRESHOLD_RULE))
    stop(label, " must be fitted with I(Q > tau); old or unmarked caches cannot be reused.",
         call. = FALSE)
  invisible(TRUE)
}

# Nodewise estimation of the covariance matrix.
# Always returns a list: Theta (the nodewise-approximate inverse,
# unscaled -- callers multiply by n, matching the previous convention),
# ifnode (always TRUE), plus lambda_node and nz_node, the per-row nodewise-Lasso
# penalty and selected nonzero-coefficient count needed for the addlength
# bias-correction term below (each row of the approximate inverse comes from
# its own separate Lasso fit, with its own lambda and its own sparsity, so a
# single matrix-wide number cannot stand in for either).
safe_inverse <- function(X, cond_threshold = 1e8, ic = 'GIC', context = NULL) {
  # Retain the old name/argument for callers; never check rank or call solve().
  nw <- nodewise_rows(X, rows = 1:ncol(X), ic = ic, context = context)
  list(Theta = nw$Theta / nrow(X), ifnode = TRUE,
       lambda_node = nw$lambda, nz_node = nw$nz, kkt_gap = nw$kkt_gap)
}
# Nodewise-Lasso rows of the precision matrix of Y, restricted to the requested
# column indices `rows` (not divided by n yet -- callers divide, matching the
# original estimate_nodewise_covariance convention). Also returns, per row,
# the selected lambda and the number of nonzero coefficients in that row's
# own Lasso fit (lambda and nz below) -- these are row-specific and must not
# be substituted with a different model's lambda or a matrix-wide sparsity
# fraction (see the addlength computation in Threshold_debiased_Lasso).
nodewise_rows <- function(Y, rows = 1:ncol(Y), ic = 'GIC', context = NULL) {
  p <- ncol(Y)
  n <- nrow(Y)

  C <- matrix(0, length(rows), p)
  tau <- NULL
  lambda_row <- NULL
  nz_row <- NULL
  kkt_gap <- NULL
  column_rms <- sqrt(colMeans(Y^2))

  # Numerical audit only: no new penalty, fit precision, or lambda-grid rule.
  audit_path <- getOption("tipping.nodewise.audit_path", "")
  audit_city <- getOption("tipping.nodewise.city", "unspecified city")
  audit_stage <- if (is.null(context$stage)) "unspecified" else context$stage
  audit_tau <- if (is.null(context$tau)) NA_real_ else context$tau
  write_audit <- function(record) {
    if (!nzchar(audit_path)) return(invisible(NULL))
    append <- file.exists(audit_path)
    write.table(record, audit_path, sep = ",", row.names = FALSE,
                col.names = !append, append = append, na = "NA")
  }

  # Loop over the requested rows only
  for (idx in seq_along(rows)) {
    j <- rows[idx]
    # Estimate the Lasso
    path_warnings <- character()
    path_error <- NULL
    jlas <- tryCatch(withCallingHandlers(glmnet(
      x = Y[, -j],
      y = Y[, j],
      family = 'gaussian',
      intercept = FALSE,
      standardize = TRUE,
      maxit = 1000000L
    ), warning = function(w) {
      path_warnings <<- c(path_warnings, conditionMessage(w))
    }), error = function(e) {
      path_error <<- conditionMessage(e)
      NULL
    })
    record <- data.frame(city = audit_city, stage = audit_stage, tau = audit_tau,
      row = j, n = n, p = p, maxit = 1000000L,
      jerr = if (is.null(jlas)) NA_integer_ else jlas$jerr,
      npasses = if (is.null(jlas)) NA_integer_ else jlas$npasses,
      returned_lambdas = if (is.null(jlas)) 0L else length(jlas$lambda),
      last_lambda = if (is.null(jlas)) NA_real_ else tail(jlas$lambda, 1),
      selected_index = NA_integer_, selected_lambda = NA_real_, selected_df = NA_integer_,
      selected_kkt_relative = NA_real_, selected_relative_duality_gap = NA_real_,
      warnings = paste(path_warnings, collapse = " | "),
      error = if (is.null(path_error)) "" else path_error)
    # Normal automatic early stopping (jerr=0) is allowed. An abnormal
    # truncated path must not be used for GIC or subsequent inference.
    if (is.null(jlas) || jlas$jerr != 0L || !length(jlas$lambda)) {
      write_audit(record)
      stop(sprintf("Nodewise path failed: %s; stage=%s; tau=%s; row=%d; jerr=%s; returned=%d; maxit=1000000; %s",
        audit_city, audit_stage, format(audit_tau, digits = 16), j,
        record$jerr, record$returned_lambdas, record$error), call. = FALSE)
    }
    # Get fit
    jfit <- predict(jlas, newx = Y[, -j], type = "response")
    # residuals
    jres <- matrix(Y[, j], n, length(jlas$lambda)) - jfit
    # std err
    jsig <- colSums(jres^2) / n
    # Computing information criterion
    if (ic == 'WIC')
      jbic  <- log(jsig) + jlas$df * log(n) / n * log(log(p)) # BIC (Wang,2010)
    if (ic == 'GIC')
      jbic  <- log(jsig) + jlas$df * log(p) / n * log(log(n)) # GIC
    if (ic == 'BIC')
      jbic  <- log(jsig) + jlas$df * log(n) / n  #BIC
    if (ic == 'MIC')
      jbic  <- jsig + jlas$df * log(n) * log(log(p)) / n # MC's IC
    if (ic == 'AIC')
      jbic  <- log(jsig) + 2 * jlas$df / n # AIC, on the log(RSS/n) scale
    # Index of selected model
    jind  <- which.min(jbic)
    # Get the parameters
    jpar <- jlas$beta[, jind]
    # Diagnostics of the selected solution only; never refit or change it.
    r_selected <- jres[, jind]
    score_standardized <- as.numeric(crossprod(Y[, -j, drop = FALSE], r_selected)) /
      (n * column_rms[-j])
    selected_lambda <- jlas$lambda[jind]
    violation <- ifelse(jpar != 0,
      abs(score_standardized - selected_lambda * sign(jpar)),
      pmax(abs(score_standardized) - selected_lambda, 0))
    primal <- mean(r_selected^2) / 2 + selected_lambda * sum(column_rms[-j] * abs(jpar))
    dual_residual <- min(1, selected_lambda / max(abs(score_standardized))) * r_selected
    dual <- mean(Y[, j] * dual_residual - dual_residual^2 / 2)
    record$selected_index <- jind
    record$selected_lambda <- selected_lambda
    record$selected_df <- jlas$df[jind]
    record$selected_kkt_relative <- max(violation) / selected_lambda
    record$selected_relative_duality_gap <- (primal - dual) / primal
    write_audit(record)
    # Computing tau squared
    # Agreed normalization correction: match diag(Theta * crossprod(Y)/n)=1.
    # The original glmnet path, standardization and row-wise GIC stay unchanged.
    jtau <- mean(Y[, j] * jres[, jind])
    if (!is.finite(jtau) || jtau <= 0) stop("Nonpositive nodewise normalizer at column ", j)
    kkt_gap <- c(kkt_gap, jtau - (sum(jres[, jind]^2) / n +
      jlas$lambda[jind] * sum(column_rms[-j] * abs(jpar))))
    # Storing the parameters
    C[idx, -j] <- -jpar
    C[idx, j] <- 1
    tau <- c(tau, jtau)
    lambda_row <- c(lambda_row, jlas$lambda[jind])
    nz_row <- c(nz_row, sum(jpar != 0))
  }

  # Construct T-squared inverse
  T2inv <- diag(1 / tau, nrow = length(rows))

  # Construct Theta-hat (rows of the full p x p Theta, restricted to `rows`)
  Theta <- T2inv %*% C

  return(list(Theta = Theta, tau = tau, lambda = lambda_row, nz = nz_row, kkt_gap = kkt_gap))
}

estimate_nodewise_covariance <- function(Y, ic = 'GIC') {
  p <- ncol(Y)
  res <- nodewise_rows(Y, rows = 1:p, ic = ic)
  Theta <- res$Theta

  # sparsity
  sp <- sum(Theta == 0) / (p^2)

  return(list(
    Theta = Theta / nrow(Y),
    ####inverse of t(Y) %*% Y
    sparsity = sp
  ))
}


fit_glmnet_ic <- function(X, y, ic, n, p) {
  fit <- glmnet::glmnet(
    X,
    y,
    standardize = FALSE,
    intercept = FALSE
  )

  y_hat <- predict(fit, newx = X)
  rss <- colSums((matrix(y, nrow = n, ncol = ncol(y_hat)) - y_hat)^2)
  jsig <- pmax(rss / n, .Machine$double.eps)
  df <- fit$df

  ic_vec <- switch(
    tolower(ic),
    bic = log(jsig) + log(n) * df / n,
    gic = log(jsig) + log(p) * df / n * log(log(n)),
    wic = log(jsig) + log(n) * df / n * log(log(p)),
    mic = jsig + df * log(n) * log(log(p)) / n,
    aic = log(jsig) + 2 * df / n,
    stop("Unsupported IC type")
  )

  best_index <- which.min(ic_vec)
  lambda_val <- fit$lambda[best_index]

  list(
    coefficients = as.numeric(coef(fit, s = lambda_val)[-1]),
    lambda_val = lambda_val,
    rss = rss[best_index],
    df = df[best_index],
    ic = ic_vec[best_index]
  )
}



# Same demeaning fix as Threshold_Lasso_fixed_tau (see its comment), applied
# at every tau tried in the grid search, not just the winning one -- the
# Objective values being compared across tau must all come from the same
# (correct) fitting procedure. The winning iteration's state (means, raw
# int2/X2, column norms) is snapshotted into best_* variables and used to
# reconstruct the final output on the raw data scale, exactly matching
# Threshold_Lasso_fixed_tau's return interface. This also fixes a latent bug
# in the original code: it reconstructed the returned X1/lambda_val from
# whatever X11/col.norm1/lambda_val happened to be left over from the LAST
# tau tried in the loop, not the winning tau -- invisible whenever
# intercept2 = FALSE (X11 doesn't depend on tau then, which is all
# SimulationSyn.R uses), but a real bug whenever intercept2 = TRUE.
# Prepare the same weighted design as Threshold_Lasso, without retaining an
# n x 2p matrix for every tau. RMS weights correspond to D(tau) in the paper.
threshold_gic_design <- function(X1, y, thdvar, tau, intercept1, intercept2) {
  X2 <- X1
  X2[thdvar <= tau, ] <- 0
  int2 <- if (intercept2) threshold_indicator_gt(thdvar, tau) else NULL
  X1_mean <- if (intercept1) colMeans(X1) else NULL
  X2_mean <- if (intercept1) colMeans(X2) else NULL
  y_mean <- if (intercept1) mean(y) else NULL
  int2_mean <- if (intercept1 && intercept2) mean(int2) else NULL
  X1_fit <- if (intercept1) sweep(X1, 2, X1_mean, "-") else X1
  X2_fit <- if (intercept1) sweep(X2, 2, X2_mean, "-") else X2
  y_fit <- if (intercept1) y - y_mean else y
  if (intercept2) {
    int2_fit <- if (intercept1) int2 - int2_mean else int2
    X1_fit <- cbind(int2_fit, X1_fit)
  }
  norms1 <- sqrt(colMeans(X1_fit^2))
  norms2 <- sqrt(colMeans(X2_fit^2))
  if (any(!is.finite(c(norms1, norms2))) || any(c(norms1, norms2) <= 0))
    stop("Zero/nonfinite column norm in the threshold design at tau = ", tau)
  list(X = cbind(sweep(X1_fit, 2, norms1, "/"),
                 sweep(X2_fit, 2, norms2, "/")), y = as.numeric(y_fit),
       scales = list(col.norm1 = 1 / norms1, col.norm2 = 1 / norms2,
         p1_fit = ncol(X1_fit), X1_mean = X1_mean, X2_mean = X2_mean,
         y_mean = y_mean, int2_mean = int2_mean))
}

# Among exactly equal minima, return the largest candidate tau.
# No tolerance: a nearly flat profile still uses its actual minimum.
threshold_argmin_upper <- function(scores, taugrid) {
  tied <- which(scores == min(scores))
  tied[which.max(taugrid[tied])]
}

# One numeric lambda grid for the WHOLE threshold search. At each lambda,
# first minimize Qn + 2*c_rho*lambda*|tau-tau_star| over tau, THEN evaluate
# GIC on that selected model. No nodewise fits or post-GIC refits occur here.
threshold_gic_profile <- function(X1, y, thdvar, taugrid, intercept1, intercept2,
                                  tau_star, c_rho, lambda_grid = NULL,
                                  nlambda = 100L, lambda_min_ratio = NULL) {
  n <- nrow(X1)
  p <- 2L * ncol(X1) + as.integer(intercept2)
  if (length(y) != n || length(thdvar) != n || any(!is.finite(X1)) ||
      any(!is.finite(y)) || any(!is.finite(thdvar)))
    stop("X1, y and thdvar must be finite and have matching sample sizes.")
  lambda_grid_input <- lambda_grid
  lambda_min_ratio_input <- lambda_min_ratio
  if (is.null(lambda_grid)) {
    if (length(nlambda) != 1L || !is.finite(nlambda) || nlambda < 2 ||
        nlambda != as.integer(nlambda)) stop("nlambda must be an integer >= 2.")
    if (is.null(lambda_min_ratio)) lambda_min_ratio <- if (n < p) 0.01 else 0.0001
    if (length(lambda_min_ratio) != 1L || !is.finite(lambda_min_ratio) ||
        lambda_min_ratio <= 0 || lambda_min_ratio >= 1)
      stop("lambda_min_ratio must lie strictly between 0 and 1.")
    # KKT entry value: all penalized coefficients are zero at this lambda
    # for EVERY tau. This uses all candidate thresholds, not the anchor alone.
    lambda_max <- 0
    for (tau in taugrid) {
      z <- threshold_gic_design(X1, y, thdvar, tau, intercept1, intercept2)
      lambda_max <- max(lambda_max, max(abs(crossprod(z$X, z$y))) / n)
    }
    if (!is.finite(lambda_max) || lambda_max <= 0)
      stop("Cannot construct a positive common lambda grid for this sample.")
    lambda_grid <- exp(seq(log(lambda_max),
      log(lambda_max * lambda_min_ratio), length.out = nlambda))
    lambda_grid[1] <- lambda_max
    grid_rule <- "global_kkt_log_grid"
  } else {
    if (!is.numeric(lambda_grid) || length(lambda_grid) < 2L ||
        any(!is.finite(lambda_grid)) || any(lambda_grid <= 0) ||
        any(diff(lambda_grid) >= 0))
      stop("lambda_grid must contain at least two positive, strictly decreasing values.")
    lambda_grid <- as.numeric(lambda_grid)
    lambda_max <- lambda_grid[1]
    lambda_min_ratio <- tail(lambda_grid, 1) / lambda_max
    grid_rule <- "supplied_common_grid"
  }

  K <- length(taugrid)
  L <- length(lambda_grid)
  rss_path <- penalty_path <- df_path <- matrix(NA_real_, K, L)
  beta_path <- scales <- vector("list", K)
  for (i in seq_len(K)) {
    z <- threshold_gic_design(X1, y, thdvar, taugrid[i], intercept1, intercept2)
    fit <- glmnet::glmnet(z$X, z$y, lambda = lambda_grid,
                         standardize = FALSE, intercept = FALSE)
    # Never interpolate a missing endpoint or compare different path columns.
    if (!isTRUE(fit$jerr == 0) || length(fit$lambda) != L ||
        !isTRUE(all.equal(as.numeric(fit$lambda), lambda_grid, tolerance = 1e-12)))
      stop("Incomplete common lambda path at tau = ", taugrid[i],
           "; glmnet jerr = ", fit$jerr)
    residuals <- z$y - as.matrix(z$X %*% fit$beta)
    rss_path[i, ] <- colSums(residuals^2)
    penalty_path[i, ] <- as.numeric(Matrix::colSums(abs(fit$beta)))
    df_path[i, ] <- fit$df
    beta_path[[i]] <- fit$beta
    scales[[i]] <- z$scales
  }
  Q_path <- rss_path / n + sweep(penalty_path, 2, 2 * lambda_grid, "*")
  objective_path <- Q_path + outer(abs(taugrid - tau_star), 2 * c_rho * lambda_grid)
  if (any(!is.finite(objective_path))) stop("Nonfinite common-lambda objective.")
  # Exact tau ties: largest tau; GIC ties: first (largest) lambda.
  tau_index <- vapply(seq_len(L), function(j)
    threshold_argmin_upper(objective_path[, j], taugrid), integer(1))
  at <- cbind(tau_index, seq_len(L))
  gic <- log(pmax(rss_path[at] / n, .Machine$double.eps)) +
    log(p) * log(log(n)) * df_path[at] / n
  if (any(!is.finite(gic))) stop("Nonfinite outer GIC.")
  selected_lambda <- which.min(gic)
  chosen_lambda <- lambda_grid[selected_lambda]
  candidates <- lapply(seq_len(K), function(i) c(list(
    coefficients = as.numeric(beta_path[[i]][, selected_lambda]),
    lambda_val = chosen_lambda, rss = rss_path[i, selected_lambda],
    Qn = Q_path[i, selected_lambda], df = df_path[i, selected_lambda],
    penalty_l1 = penalty_path[i, selected_lambda]), scales[[i]]))
  selection <- data.frame(lambda_glmnet = lambda_grid, lambda_paper = 2 * lambda_grid,
    rho_n = 2 * c_rho * lambda_grid, tau = taugrid[tau_index],
    rss = rss_path[at], df = df_path[at], penalty_l1 = penalty_path[at],
    Qn = Q_path[at], anchored_objective = objective_path[at], GIC = gic,
    selected = seq_len(L) == selected_lambda)
  list(X1 = X1, y = y, thdvar = thdvar, taugrid = taugrid, candidates = candidates,
    threshold_rule = DTLASSO_THRESHOLD_RULE,
    lambda = "gic", intercept1 = intercept1, intercept2 = intercept2,
    tuning = "common_lambda_anchored_gic", tau_star = tau_star, c_rho = c_rho,
    tau_tie_rule = "largest_tau_exact",
    lambda_grid_input = lambda_grid_input, nlambda_input = nlambda,
    lambda_min_ratio_input = lambda_min_ratio_input,
    lambda_search = list(table = selection, selected_index = selected_lambda,
      grid_rule = grid_rule, lambda_max = lambda_max, nlambda = L,
      lambda_min_ratio = lambda_min_ratio,
      selected_at_boundary = selected_lambda %in% c(1L, L)),
    table = data.frame(tau = taugrid, Qn = Q_path[, selected_lambda],
      lambda_glmnet = chosen_lambda, lambda_paper = 2 * chosen_lambda,
      rss = rss_path[, selected_lambda], df = df_path[, selected_lambda],
      penalty_l1 = penalty_path[, selected_lambda],
      anchored_objective = objective_path[, selected_lambda]))
}

Threshold_Lasso <- function(X1,
                            y,
                            taugrid = taugrid,
                            lambda = NULL,
                            intercept1 = TRUE,
                            intercept2 = TRUE,
                            nfolds = 10,
                            tau_star = 0.3,
                            c_rho = 1,
                            profile = NULL,
                            lambda_grid = NULL,
                            nlambda = 100L,
                            lambda_min_ratio = NULL) {
  p1 <- p2 <- ncol(X1)

  n <- nrow(X1)

  if (length(c_rho) != 1L || !is.finite(c_rho) || c_rho < 0)
    stop("c_rho must be one finite nonnegative number.")
  if (length(tau_star) != 1L || !is.finite(tau_star) ||
      tau_star < min(taugrid) || tau_star > max(taugrid))
    stop("tau_star must lie within the search interval.")
  taugrid <- sort(unique(c(taugrid, tau_star)))
  if (is.null(profile) && identical(lambda, "gic")) {
    profile <- threshold_gic_profile(X1, y, thdvar, taugrid, intercept1, intercept2,
      tau_star, c_rho, lambda_grid, nlambda, lambda_min_ratio)
  }
  if (!is.null(profile)) assert_gt_threshold_object(profile, "Cached candidate profile")
  if (identical(lambda, "gic") &&
      (!identical(profile$tuning, "common_lambda_anchored_gic") ||
       !identical(profile$tau_tie_rule, "largest_tau_exact") ||
       !identical(profile$tau_star, tau_star) || !identical(profile$c_rho, c_rho) ||
       !identical(profile$lambda_grid_input, lambda_grid) ||
       !identical(profile$nlambda_input, nlambda) ||
       !identical(profile$lambda_min_ratio_input, lambda_min_ratio)))
    stop("Cached GIC profile must match the tau tie rule, anchor, c_rho and common lambda-grid settings.")
  if (is.null(profile)) {
  candidates <- vector("list", length(taugrid))
  # Iterate over each tau in the grid
  for (cnt in 1:length(taugrid)) {
    tau <- taugrid[cnt]
    p1 <- p2 <- ncol(X1)

    int2_raw <- NULL
    if (intercept2 == TRUE) {
      int2_raw <- rep(0, n)
      int2_raw[thdvar > tau] <- 1
    }
    X2_raw <- X1
    X2_raw[thdvar <= tau, ] <- 0

    if (intercept1 == TRUE) {
      X1_mean <- colMeans(X1)
      X2_mean <- colMeans(X2_raw)
      y_mean <- mean(y)
      X1_fit <- sweep(X1, 2, X1_mean, "-")
      X2_fit <- sweep(X2_raw, 2, X2_mean, "-")
      y_fit <- y - y_mean
      if (intercept2 == TRUE) {
        int2_mean <- mean(int2_raw)
        int2_fit <- int2_raw - int2_mean
      }
    } else {
      X1_fit <- X1
      X2_fit <- X2_raw
      y_fit <- y
      if (intercept2 == TRUE) int2_fit <- int2_raw
    }

    X11 <- X1_fit
    if (intercept2 == TRUE) {
      X11 <- cbind(int2_fit, X11)
      p1 <- p1 + 1
    }
    p1_fit <- p1
    p <- p1_fit + p2
    col.norm1 <- 1 / sqrt((1 / n) * diag(t(X11) %*% X11))
    X11 <- X11 %*% diag(col.norm1)
    col.norm2 <- 1 / sqrt((1 / n) * diag(t(X2_fit) %*% X2_fit))

    X2 <- X2_fit %*% diag(col.norm2)


    X <- cbind(X11, X2)

    # Compute the lambda
    if (is.null(lambda)) {
      # ---- Theoretical lambda (default) ----
      lambda_val <- sqrt(qnorm(1 - (0.1 / p)) / n)
      outLas <- glmnet::glmnet(
        X,
        y_fit,
        lambda = lambda_val,
        standardize = FALSE,
        intercept =  FALSE
      )
      coefficients <- coef(outLas)[-1]
      residuals <- y_fit - X %*% coefficients
      rss <- sum(residuals^2)
      penalty <- sum(abs(coefficients))
      # Calculate the Objective Function
      Objective <- rss / n + 2 * lambda_val * penalty
    } else if (lambda %in% c("bic", "gic", "wic", "aic", "mic")) {
      # ---- IC-based lambda selection ----
      ic_fit <- fit_glmnet_ic(X, y_fit, lambda, n, p)
      lambda_val <- ic_fit$lambda_val
      coefficients <- ic_fit$coefficients
      rss <- ic_fit$rss
      penalty <- sum(abs(coefficients))
      Objective <- rss / n + 2 * lambda_val * penalty
    } else if (identical(lambda, "cv")) {
      # ---- CV-based lambda selection ----
      outLas <-  glmnet::cv.glmnet(
        X,
        y_fit,
        nfolds = nfolds,
        standardize = TRUE,
        intercept =  FALSE
      )
      coefficients <- matrix(coef(outLas, s = 'lambda.min'))[-1]  # Adjust depending on the structure of your model
      # Calculate residuals
      predicted <- predict(outLas, X, s = 'lambda.min')
      residuals <- y_fit - predicted
      rss <- sum(residuals^2)
      # Calculate the penalty (L1 norm of the coefficients for Lasso)
      penalty <- sum(abs(coefficients))
      # Calculate the objective function
      Objective <- rss / n + 2 * outLas$lambda.min * penalty
      lambda_val <-   outLas$lambda.min
      #if not given
    } else if (is.numeric(lambda)) {
      # ---- Externally-chosen lambda (e.g. from a caller's own CV loop) --
      # same fit/Objective computation as the theoretical-lambda branch
      # above, just using the given number instead of the sqrt-lasso
      # formula. Needed when a caller wants to score/refit at a SPECIFIC
      # lambda value chosen outside this function (own CV grid, etc.) --
      # 'gic'/'bic'/'cv' above only ever select lambda themselves.
      lambda_val <- lambda
      outLas <- glmnet::glmnet(
        X,
        y_fit,
        lambda = lambda_val,
        standardize = FALSE,
        intercept =  FALSE
      )
      coefficients <- coef(outLas)[-1]
      residuals <- y_fit - X %*% coefficients
      rss <- sum(residuals^2)
      penalty <- sum(abs(coefficients))
      Objective <- rss / n + 2 * lambda_val * penalty
    } else {
      stop("lambda must be one of: NULL, 'bic'/'gic'/'wic'/'aic'/'mic', 'cv', or a single numeric value")
    }

    # Save only coefficients/scales, not a separate n x p matrix per tau.
    # Cache the fitted coefficients/scales for compatible repeated calls.
    candidates[[cnt]] <- list(
      coefficients = as.numeric(coefficients), col.norm1 = col.norm1, col.norm2 = col.norm2,
      p1_fit = p1_fit, lambda_val = lambda_val, rss = rss, Qn = Objective,
      X1_mean = if (intercept1) X1_mean else NULL,
      X2_mean = if (intercept1) X2_mean else NULL,
      y_mean = if (intercept1) y_mean else NULL,
      int2_mean = if (intercept1 && intercept2) int2_mean else NULL)
  }
  profile <- list(
    X1 = X1, y = y, thdvar = thdvar, taugrid = taugrid, candidates = candidates,
    threshold_rule = DTLASSO_THRESHOLD_RULE,
    lambda = lambda, intercept1 = intercept1, intercept2 = intercept2,
    table = data.frame(tau = taugrid,
      Qn = vapply(candidates, function(f) f$Qn, numeric(1)),
      lambda_glmnet = vapply(candidates, function(f) f$lambda_val, numeric(1)),
      lambda_paper = 2 * vapply(candidates, function(f) f$lambda_val, numeric(1))))
  } else {
    if (!identical(profile$X1, X1) || !identical(profile$y, y) ||
        !identical(profile$thdvar, thdvar) || !identical(profile$taugrid, taugrid) ||
        !identical(profile$lambda, lambda) || !identical(profile$intercept1, intercept1) ||
        !identical(profile$intercept2, intercept2)) stop("Cached candidate fits do not match this call.")
  }

  # Keep the previous output-field name. Under common-lambda GIC this is
  # the shared lambda_paper, NOT a lambda tuned separately at the anchor.
  lambda_paper_anchor <- profile$table$lambda_paper[match(tau_star, taugrid)]
  rho_n <- c_rho * lambda_paper_anchor
  scores <- profile$table$Qn + rho_n * abs(taugrid - tau_star)
  selected <- threshold_argmin_upper(scores, taugrid)
  best <- profile$candidates[[selected]]
  best_tau <- taugrid[selected]
  best_objective <- scores[selected]
  best_coefficients <- best$coefficients
  best_col.norm1 <- best$col.norm1
  best_col.norm2 <- best$col.norm2
  best_p1_fit <- best$p1_fit
  best_lambda_val <- best$lambda_val
  best_int2_raw <- if (intercept2) threshold_indicator_gt(thdvar, best_tau) else NULL
  best_X2_raw <- X1
  best_X2_raw[thdvar <= best_tau, ] <- 0
  if (intercept1) {
    best_X1_mean <- best$X1_mean
    best_X2_mean <- best$X2_mean
    best_y_mean <- best$y_mean
    if (intercept2) best_int2_mean <- best$int2_mean
  }

  col.norm <- c(best_col.norm1, best_col.norm2)
  best_coefficients <- best_coefficients %*% diag(col.norm)   # on the (demeaned-)X scale

  off <- as.integer(intercept2)
  deltaint_hat <- if (intercept2 == TRUE) best_coefficients[1] else NULL
  beta_hat  <- best_coefficients[(off + 1):best_p1_fit]
  delta_hat <- best_coefficients[(best_p1_fit + 1):(best_p1_fit + p2)]

  if (intercept1 == TRUE) {
    betaint_hat <- best_y_mean -
      (if (intercept2 == TRUE) deltaint_hat * best_int2_mean else 0) -
      as.numeric(best_X1_mean %*% beta_hat) - as.numeric(best_X2_mean %*% delta_hat)
  }

  # Reassemble to match the ORIGINAL interface exactly: X1 = [1?, int2?, X1],
  # coefficients = [betaint?, deltaint?, beta, delta], all on the raw data scale.
  X11_out <- X1
  if (intercept2 == TRUE) X11_out <- cbind(best_int2_raw, X11_out)
  if (intercept1 == TRUE) X11_out <- cbind(1, X11_out)

  coef_out <- beta_hat
  if (intercept2 == TRUE) coef_out <- c(deltaint_hat, coef_out)
  if (intercept1 == TRUE) coef_out <- c(betaint_hat, coef_out)

  # Output depending on whether intercept is requested
  return(
    list(
      coefficients = matrix(c(coef_out, delta_hat), nrow = 1),
      tau = best_tau,
      objective = best_objective,
      X1 = X11_out,
      X_tau = best_X2_raw,
      lambda_val = best_lambda_val,
      # Means used to reconstruct betaint_hat above -- exposed so
      # Threshold_debiased_Lasso can reconstruct betaint from the DEBIASED
      # slopes with the same formula, instead of feeding a literal constant
      # column into the nodewise/debiasing machinery (see
      # INTERCEPT_DEMEAN_STANDARD_ERROR_NOTES.md sec 4).
      X1_mean = if (intercept1 == TRUE) best_X1_mean else NULL,
      X2_mean = if (intercept1 == TRUE) best_X2_mean else NULL,
      int2_mean = if (intercept1 == TRUE && intercept2 == TRUE) best_int2_mean else NULL,
      y_mean = if (intercept1 == TRUE) best_y_mean else NULL,
      Qn = best$Qn, lambda_paper = 2 * best_lambda_val,
      lambda_paper_anchor = lambda_paper_anchor, rho_n = rho_n, c_rho = c_rho,
      tau_star = tau_star, profile = profile,
      threshold_rule = DTLASSO_THRESHOLD_RULE
    )
  )
}

# Single-tau analog of Threshold_Lasso: fits alpha-hat(tau) at one fixed
# threshold (no grid search / objective comparison). Used by
# threshold_effect_test, which needs the fixed-tau Lasso fit at every tau in
# the testing grid rather than only at the tau that minimizes the objective.
# thdvar is passed explicitly (unlike Threshold_Lasso, which reads it from the
# caller's environment) since this function is called from an lapply/mclapply
# loop over tau.
# See Linear_Lasso's comment for why a literal constant column never gets a
# nonzero coefficient from glmnet(..., intercept=FALSE). When intercept1 is
# TRUE, this function demeans X1 (beta block), X2(tau) (delta block), and
# int2 (if intercept2 is also TRUE) -- EACH BY ITS OWN actual mean, never by
# a different column's mean -- fits with no constant column at all, and
# reconstructs betaint afterward as ybar - (weighted means)'(fitted slopes).
# Demeaning int2 too (not just X1/X2) turned out to matter: leaving int2 raw
# while X1/X2 are demeaned means int2 (nonzero mean) is regressed against
# mean-zero columns during nodewise/debiasing, which is itself a structural
# mismatch that degrades both betaint and deltaint; demeaning int2 as well fixed
# it (verified: betaint/deltaint debiased accuracy improved substantially, delta
# improved too). Verified this demeaning does not distort delta itself
# (exact match to a with-intercept OLS fit in the unpenalized/low-dim case).
# When intercept1 is FALSE, nothing is demeaned -- matches prior behavior
# exactly, including intercept2-alone (which already works fine as a raw
# column, no demeaning needed there).
Threshold_Lasso_fixed_tau <- function(X1,
                                       y,
                                       tau,
                                       thdvar,
                                       lambda = NULL,
                                       intercept1 = TRUE,
                                       intercept2 = TRUE,
                                       nfolds = 10) {
  p1 <- p2 <- ncol(X1)

  n <- nrow(X1)

  int2_raw <- NULL
  if (intercept2 == TRUE) {
    int2_raw <- rep(0, n)
    int2_raw[thdvar > tau] <- 1
  }
  X2_raw <- X1
  X2_raw[thdvar <= tau, ] <- 0

  if (intercept1 == TRUE) {
    X1_mean <- colMeans(X1)
    X2_mean <- colMeans(X2_raw)
    y_mean <- mean(y)
    X1_fit <- sweep(X1, 2, X1_mean, "-")
    X2_fit <- sweep(X2_raw, 2, X2_mean, "-")
    y <- y - y_mean
    if (intercept2 == TRUE) {
      int2_mean <- mean(int2_raw)
      int2_fit <- int2_raw - int2_mean
    }
  } else {
    X1_fit <- X1
    X2_fit <- X2_raw
    if (intercept2 == TRUE) int2_fit <- int2_raw
  }

  X11 <- X1_fit
  if (intercept2 == TRUE) {
    X11 <- cbind(int2_fit, X11)
    p1 <- p1 + 1
  }
  p1_fit <- p1
  p <- p1_fit + p2
  col.norm1 <- 1 / sqrt((1 / n) * diag(t(X11) %*% X11))
  X11 <- X11 %*% diag(col.norm1)
  col.norm2 <- 1 / sqrt((1 / n) * diag(t(X2_fit) %*% X2_fit))

  X2 <- X2_fit %*% diag(col.norm2)

  X <- cbind(X11, X2)

  # Compute the lambda
  if (is.null(lambda)) {
    # ---- Theoretical lambda (default) ----
    lambda_val <- sqrt(qnorm(1 - (0.1 / p)) / n)
    outLas <- glmnet::glmnet(
      X,
      y,
      lambda = lambda_val,
      standardize = FALSE,
      intercept =  FALSE
    )
    coefficients <- coef(outLas)[-1]
  } else if (is.numeric(lambda)) {
    if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
      stop("A supplied fixed-tau lambda must be one finite nonnegative number.")
    lambda_val <- lambda
    outLas <- glmnet::glmnet(X, y, lambda = lambda_val,
                            standardize = FALSE, intercept = FALSE)
    coefficients <- as.numeric(outLas$beta[, 1])
  } else if (lambda %in% c("bic", "gic", "wic", "aic", "mic")) {
    # ---- IC-based lambda selection ----
    ic_fit <- fit_glmnet_ic(X, y, lambda, n, p)
    lambda_val <- ic_fit$lambda_val
    coefficients <- ic_fit$coefficients
  } else if (identical(lambda, "cv")) {
    # ---- CV-based lambda selection ----
    outLas <-  glmnet::cv.glmnet(
      X,
      y,
      nfolds = nfolds,
      standardize = TRUE,
      intercept =  FALSE
    )
    coefficients <- matrix(coef(outLas, s = 'lambda.min'))[-1]
    lambda_val <-   outLas$lambda.min
  } else {
    stop("lambda must be NULL, an IC/CV method, or a single numeric value")
  }

  col.norm <- c(col.norm1, col.norm2)
  coefficients <- coefficients %*% diag(col.norm)   # on the (demeaned-)X scale

  off <- as.integer(intercept2)
  deltaint_hat <- if (intercept2 == TRUE) coefficients[1] else NULL
  beta_hat  <- coefficients[(off + 1):p1_fit]
  delta_hat <- coefficients[(p1_fit + 1):(p1_fit + p2)]

  if (intercept1 == TRUE) {
    betaint_hat <- y_mean -
      (if (intercept2 == TRUE) deltaint_hat * int2_mean else 0) -
      as.numeric(X1_mean %*% beta_hat) - as.numeric(X2_mean %*% delta_hat)
  }

  # Reassemble to match the ORIGINAL interface exactly: X1 = [1?, int2?, X1],
  # coefficients = [betaint?, deltaint?, beta, delta], all on the raw data scale.
  X11_out <- X1
  if (intercept2 == TRUE) X11_out <- cbind(int2_raw, X11_out)
  if (intercept1 == TRUE) X11_out <- cbind(1, X11_out)

  coef_out <- beta_hat
  if (intercept2 == TRUE) coef_out <- c(deltaint_hat, coef_out)
  if (intercept1 == TRUE) coef_out <- c(betaint_hat, coef_out)
  p1_out <- length(coef_out)

  return(
    list(
      coefficients = matrix(c(coef_out, delta_hat), nrow = 1),
      tau = tau,
      X1 = X11_out,
      X_tau = X2_raw,
      lambda_val = lambda_val,
      p1 = p1_out,
      p2 = p2,
      threshold_rule = DTLASSO_THRESHOLD_RULE
    )
  )
}

Threshold_debiased_Lasso <- function (X1,
                                      y,
                                      alpha = 0.05,
                                      taugrid = taugrid,
                                      lambda = NULL,
                                      mu = NULL,
                                      intercept1 = TRUE,
                                      intercept2 = TRUE,
                                      nfolds = 10,
                                      resol = 1.3,
                                      maxiter = 50,
                                      threshold = 1e-2,
                                      verbose = TRUE,
                                      tau_star = 0.3,
                                      c_rho = 1,
                                      lasso_result = NULL) {
  # Compute confidence intervals and p-values.
  #
  # Args:
  #   X1     :  design matrix
  #   y     :  response
  #   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
  #
  # Returns:
  #   noise.sd: Estimate of the noise standard deviation
  #   coef    : Lasso estimated coefficients
  #   unb.coef: Unbiased coefficient estimates
  #   low.lim : Lower limits of confidence intervals
  #   up.lim  : upper limit of confidence intervals
  #   pvals   : p-values for the coefficients
  #
  
  #X1 <- X[,1:50]
  #X2 <- X[,101:200]
  
  #col.norm1 <- 1/sqrt((1/n)*diag(t(X1)%*%X1));
  #X1 <- X1 %*% diag(col.norm1);
  
  #col.norm2 <- 1/sqrt((1/n)*diag(t(X2)%*%X2));
  #X2 <- X2 %*% diag(col.norm2);
  
  p1 <- p2 <- ncol(X1)
  
  n <- nrow(X1)

  if (is.null(lasso_result)) lasso_result <- Threshold_Lasso(
    X1,
    y,
    taugrid = taugrid,
    lambda = lambda,
    intercept1 = intercept1,
    intercept2 = intercept2,
    nfolds = nfolds,
    tau_star = tau_star,
    c_rho = c_rho
  )
  assert_gt_threshold_object(lasso_result, "Initial threshold fit")
  htheta <- t(lasso_result[[1]])

  best_tau <- lasso_result[[2]]
  best_lambda_val <- lasso_result$lambda_val

  X1_out <- lasso_result$X1       # [1?, int2?, X1], raw scale
  X2 <- lasso_result$X_tau        # raw scale

  # betaint (intercept1) has no column of its own in the demeaned fit -- it is
  # a post-hoc reconstruction from means, not a genuine high-dimensional
  # covariate with a population precision-matrix row. Feeding its literal
  # "1" column into the nodewise/debiasing machinery below is a structural
  # mismatch (see INTERCEPT_DEMEAN_STANDARD_ERROR_NOTES.md sec 4): drop it
  # here, debias deltaint/beta/delta normally (unaffected), then reconstruct
  # betaint's point estimate afterward from the DEBIASED slopes via the same
  # mean-weighted formula Threshold_Lasso uses for the raw-Lasso betaint_hat.
  # No valid SE for it -> NA in noise.sd/low.lim/up.lim/pvals.
  Xb1 <- if (intercept1 == TRUE) X1_out[, -1, drop = FALSE] else X1_out

  col.norm1 <- 1 / sqrt((1 / n) * diag(t(Xb1) %*% Xb1))
  col.norm2 <- 1 / sqrt((1 / n) * diag(t(X2) %*% X2))
  col.norm <- c(col.norm1, col.norm2)

  p1 <- ncol(Xb1)   # int2 (if any) + X1, WITHOUT betaint
  p2 <- ncol(X2)
  pp <- p1 + p2

  Xb <- cbind(Xb1, X2)
  # Application with an unpenalized baseline intercept: use the same
  # centered design as the initial Lasso, also in nodewise regressions.
  if (intercept1) Xb <- sweep(Xb, 2, colMeans(Xb), "-")
  sigma.hat <- (1 / n) * (t(Xb) %*% Xb)
  invgram <- safe_inverse(Xb, context = list(stage = "estimation", tau = best_tau))
  M <- invgram$Theta * n
  ifnode <- invgram$ifnode

  htheta_fit <- if (intercept1 == TRUE) htheta[-1] else htheta   # drop betaint, align with Xb

  # Frisch--Waugh residualization of the unpenalized baseline intercept.
  # y_centered - X_centered %*% slopes equals the original fitted residual.
  y_net <- if (intercept1 == TRUE) y - lasso_result$y_mean else y

  unbiased.Lasso <- as.numeric(htheta_fit + (M %*% t(Xb) %*% (y_net - Xb %*% (htheta_fit))) /
                                 n)

  # noise <- NoiseSd( unbiased.Lasso,M, n);
  # s.hat <- noise$sd;

  resids <- y_net - Xb %*% (htheta_fit)
  u2 <- as.numeric(resids)^2
  meat <- crossprod(Xb * sqrt(u2)) / n
  varmx <-   M %*% meat %*% t(M) / n
  stder <- sqrt(diag(varmx))

  interval.sizes <- qnorm(1 - (alpha / 2)) * stder
  low_wald <- unbiased.Lasso - interval.sizes
  up_wald <- unbiased.Lasso + interval.sizes
  addlength <- rep(0, pp)
  if (ifnode) {
    # Each row i's bias-correction uses THAT row's own nodewise-Lasso lambda
    # and nonzero count (invgram$lambda_node[i], invgram$nz_node[i]) -- not
    # the main model's lambda and not a matrix-wide sparsity fraction, since
    # each row of the approximate inverse comes from a separate Lasso fit.
    MM <- M %*% sigma.hat - diag(pp)
    for (i in 1:pp) {
      # nz_node already counts regressors excluding the response column.
      k <- invgram$nz_node[i]
      if (k > 0) {
        effectivemuvec <- sort(abs(MM[i, ]), decreasing = TRUE)[1:k]
        addlength[i] <- sqrt(sum(effectivemuvec^2)) * invgram$lambda_node[i]
      }
    }
  }
  
  # nz_hat <- sum(htheta != 0)
  # MM <- M%*%sigma.hat - diag(pp);
  # for (i in 1:pp){
  #   effectivemuvec <- sort(abs(MM[i,]),decreasing=TRUE);
  #   effectivemuvec <- effectivemuvec[1:nz_hat];
  #   addlength[i] <- sqrt(sum(effectivemuvec*effectivemuvec))*lambda;
  # }
  
  p.vals <- 2 * pnorm(-abs(unbiased.Lasso) / stder)

  # ---- Reconstruct betaint (intercept1) from the DEBIASED slopes, point
  # estimate only -- prepend to match Threshold_Lasso's own
  # [betaint?, deltaint?, beta, delta] column order. No SE/CI/p-value: the
  # quadratic-form SE for this reconstruction is unresolved (see
  # INTERCEPT_DEMEAN_STANDARD_ERROR_NOTES.md sec 4), so those slots are NA. ----
  if (intercept1 == TRUE) {
    off <- as.integer(intercept2)
    deltaint_deb <- if (intercept2 == TRUE) unbiased.Lasso[1] else NULL
    beta_deb  <- unbiased.Lasso[(off + 1):p1]
    delta_deb <- unbiased.Lasso[(p1 + 1):(p1 + p2)]
    betaint_deb <- lasso_result$y_mean -
      (if (intercept2 == TRUE) deltaint_deb * lasso_result$int2_mean else 0) -
      as.numeric(lasso_result$X1_mean %*% beta_deb) - as.numeric(lasso_result$X2_mean %*% delta_deb)

    unb.coef_out <- c(betaint_deb, unbiased.Lasso)
    stder_out    <- c(NA, stder)
    low_out      <- c(NA, unbiased.Lasso - interval.sizes - addlength)
    up_out       <- c(NA, unbiased.Lasso + interval.sizes + addlength)
    pvals_out    <- c(NA, p.vals)
  } else {
    unb.coef_out <- unbiased.Lasso
    stder_out    <- stder
    low_out      <- unbiased.Lasso - interval.sizes - addlength
    up_out       <- unbiased.Lasso + interval.sizes + addlength
    pvals_out    <- p.vals
  }

  returnList <- list(
    "noise.sd" = stder_out,
    #"norm0" =  nz_hat,
    "coef" = htheta,
    "unb.coef" = unb.coef_out,
    "low.lim" = low_out,
    "up.lim" = up_out,
    "pvals" = pvals_out,
    "best_tau" = best_tau,
    "lambda_val" = best_lambda_val,
    "lambda_paper" = 2 * best_lambda_val,
    "low.lim.wald" = c(if (intercept1) NA_real_, low_wald),
    "up.lim.wald" = c(if (intercept1) NA_real_, up_wald),
    "low.lim.widened" = low_out, "up.lim.widened" = up_out,
    "addlength" = c(if (intercept1) NA_real_, addlength),
    "used_node" = TRUE, "inverse_method" = "full_nodewise",
    "node_kkt_gap" = max(abs(invgram$kkt_gap)),
    "gram_diag_error" = max(abs(diag(M %*% sigma.hat) - 1)),
    "threshold_rule" = DTLASSO_THRESHOLD_RULE
  )
  return(returnList)
}



# A literal constant column fed to glmnet(..., intercept=FALSE) can never
# receive a nonzero coefficient (glmnet's zero-variance-column handling makes
# it structurally impossible, independent of lambda/standardize -- verified
# empirically). intercept1 = TRUE therefore demeans y and X1 by their own
# sample means, fits with NO constant column at all, and reconstructs the
# intercept afterward as ybar - xbar'*slopes -- the standard, exact way to
# fit an unpenalized intercept with Lasso (this is what glmnet's own
# intercept=TRUE does internally). Verified exactly equivalent to a
# with-intercept fit in the unpenalized/low-dim case (max abs diff ~1e-15),
# and accurate in high-dim Lasso+debiasing (recovers a true intercept of 7
# as ~7.1). Returned coefficients/X1 still have the same shape as before
# (intercept slot first, X1 unchanged) so callers are unaffected.
Linear_Lasso <- function(X1,
                         y,
                         lambda = NULL,
                         intercept1 = TRUE,
                         col.norm1 = col.norm1,
                         nfolds = 10) {
  p1 <- ncol(X1)

  n <- nrow(X1)

  best_objective <- Inf
  best_tau <- NULL
  best_coefficients <- NULL

  x_mean <- colMeans(X1)
  y_mean <- mean(y)
  X11 <- X1
  if (intercept1 == TRUE) {
    X11 <- sweep(X11, 2, x_mean, "-")
    y <- y - y_mean
  }
  p <- p1
  col.norm1 <- 1 / sqrt((1 / n) * diag(t(X11) %*% X11))
  X11 <- X11 %*% diag(col.norm1)
  X <- X11
  
  # Compute the lambda
  if (is.null(lambda)) {
    # ---- Theoretical lambda (default) ----
    lambda <- sqrt(qnorm(1 - (0.1 / p)) / n)
    lambda_val <- lambda
    outLas <- glmnet::glmnet(
      X,
      y,
      lambda = lambda_val,
      standardize = FALSE,
      intercept =  FALSE
    )
    coefficients <- coef(outLas)[-1]
    residuals <- y - X %*% coefficients
    rss <- sum(residuals^2)
    penalty <- sum(abs(coefficients))
    # Calculate the Objective Function
    Objective <- rss / n + 2 * lambda_val * penalty
    lambda <- NULL
  } else if (lambda %in% c("bic", "gic", "wic", "aic", "mic")) {
    # ---- IC-based lambda selection ----
    ic_fit <- fit_glmnet_ic(X, y, lambda, n, p)
    lambda_val <- ic_fit$lambda_val
    coefficients <- ic_fit$coefficients
    rss <- ic_fit$rss
    penalty <- sum(abs(coefficients))
    Objective <- rss / n + 2 * lambda_val * penalty
  } else if (identical(lambda, "cv")) {
    # ---- CV-based lambda selection ----
    outLas <-  glmnet::cv.glmnet(
      X,
      y,
      nfolds = nfolds,
      standardize = FALSE,
      intercept =  FALSE
    )
    coefficients <- matrix(coef(outLas, s = 'lambda.min'))[-1]  # Adjust depending on the structure of your model
    # Calculate residuals
    predicted <- predict(outLas, X, s = 'lambda.min')
    residuals <- y - predicted
    rss <- sum(residuals^2)
    n <- length(y)
    # Calculate the penalty (L1 norm of the coefficients for Lasso)
    penalty <- sum(abs(coefficients))
    # Calculate the objective function
    Objective <- rss / n + 2 * outLas$lambda.min * penalty
    lambda_val <-   outLas$lambda.min
    #if not given
  } else {
    stop("lambda must be one of: NULL, 'bic', or 'cv'")
  }
  
  
  best_objective <- Objective
  best_coefficients <- coefficients %*% diag(col.norm1)   # slopes, on X1's original scale

  if (intercept1 == TRUE) {
    betaint_hat <- y_mean - as.numeric(x_mean %*% t(best_coefficients))
    best_coefficients <- cbind(betaint_hat, best_coefficients)
    X1_out <- cbind(1, X1)
  } else {
    X1_out <- X1
  }

  # Output depending on whether intercept is requested
  return(
    list(
      coefficients = best_coefficients,
      objective = best_objective,
      X1 = X1_out,
      lambda_val = lambda_val,
      # Exposed for Linear_debiased_Lasso to reconstruct betaint from the
      # debiased slopes -- see the parallel note in Threshold_Lasso's return.
      x_mean = if (intercept1 == TRUE) x_mean else NULL,
      y_mean = if (intercept1 == TRUE) y_mean else NULL
    )
  )
}



Linear_debiased_Lasso <- function (X1,
                                   y,
                                   alpha = 0.05,
                                   lambda = NULL,
                                   mu = NULL,
                                   intercept1 = TRUE,
                                   nfolds = 10,
                                   resol = 1.3,
                                   maxiter = 50,
                                   threshold = 1e-2,
                                   verbose = TRUE) {
  # Compute confidence intervals and p-values.
  #
  # Args:
  #   X1     :  design matrix
  #   y     :  response
  #   alpha :  significance level
  #   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
  #   mu    :  Linfty constraint on M (if null, searches)
  #   resol :  step parameter for the function that computes M
  #   maxiter: iteration parameter for computing M
  #   threshold : tolerance criterion for computing M
  #   verbose : verbose?
  #
  # Returns:
  #   noise.sd: Estimate of the noise standard deviation
  #   norm0   : Estimate of the number of 'significant' coefficients
  #   coef    : Lasso estimated coefficients
  #   unb.coef: Unbiased coefficient estimates
  #   low.lim : Lower limits of confidence intervals
  #   up.lim  : upper limit of confidence intervals
  #   pvals   : p-values for the coefficients
  #
  
  #X1 <- X[,1:50]
  #X2 <- X[,101:200]
  
  #col.norm1 <- 1/sqrt((1/n)*diag(t(X1)%*%X1));
  #X1 <- X1 %*% diag(col.norm1);
  
  #col.norm2 <- 1/sqrt((1/n)*diag(t(X2)%*%X2));
  #X2 <- X2 %*% diag(col.norm2);
  
  p1 <- ncol(X1)
  
  n <- nrow(X1)
  
  col.norm1 <- 1 / sqrt((1 / n) * diag(t(X1) %*% X1))
  
  #print(Lasso (X,y,taugrid=taugrid,lambda=lambda,intercept=intercept,col.norm = col.norm1)[[1]]);
  #print(Lasso (X,y,taugrid=taugrid,lambda=lambda,intercept=intercept,col.norm = col.norm1)[[2]]);
  lasso_result <- Linear_Lasso (
    X1,
    y,
    lambda = lambda,
    intercept1 = intercept1,
    col.norm1 = col.norm1,
    nfolds = nfolds
  )
  htheta <- t(lasso_result[[1]])
  X1_out <- lasso_result[[3]]   # [1?, X1], raw scale

  # betaint has no column of its own in the demeaned fit -- drop it from the
  # debiasing design (same reasoning as Threshold_debiased_Lasso), debias
  # the real slopes normally, then reconstruct betaint's point estimate
  # afterward from the debiased slopes. No valid SE for it -> NA below.
  Xb <- if (intercept1 == TRUE) X1_out[, -1, drop = FALSE] else X1_out
  # Partial out the unpenalized intercept, as in the initial Linear_Lasso fit.
  # Center both the nodewise design and response; this residual equals
  # y - betaint_hat - X1 %*% beta_hat on the original data scale.
  if (intercept1 == TRUE) Xb <- sweep(Xb, 2, lasso_result$x_mean, "-")
  y_net <- if (intercept1 == TRUE) y - lasso_result$y_mean else y
  col.norm1 <- 1 / sqrt((1 / n) * diag(t(Xb) %*% Xb))
  pp <- ncol(Xb)

  sigma.hat <- (1 / n) * (t(Xb) %*% Xb)

  invgram <- safe_inverse(Xb)
  M <- invgram$Theta * n

  htheta_fit <- if (intercept1 == TRUE) htheta[-1] else htheta

  unbiased.Lasso <- as.numeric(htheta_fit + (M %*% t(Xb) %*% (y_net - Xb %*% (htheta_fit))) /
                                 n)
  # noise <- NoiseSd( unbiased.Lasso,M, n);
  # s.hat <- noise$sd;

  resids <- y_net - Xb %*% (htheta_fit)
  u2 <- as.numeric(resids)^2
  meat <- crossprod(Xb * sqrt(u2)) / n
  varmx <-   M %*% meat %*% t(M) / n
  stder <- sqrt(diag(varmx))

  interval.sizes <- qnorm(1 - (alpha / 2)) * stder

  low_wald <- unbiased.Lasso - interval.sizes
  up_wald <- unbiased.Lasso + interval.sizes
  addlength <- rep(0, pp)
  # Use the same widening rule as Threshold_debiased_Lasso, reusing the
  # linear model's already-computed nodewise inverse and row-specific fits.
  MM <- M %*% sigma.hat - diag(pp)
  for (i in seq_len(pp)) {
    k <- invgram$nz_node[i]
    if (k > 0) {
      effectivemuvec <- sort(abs(MM[i, ]), decreasing = TRUE)[seq_len(k)]
      addlength[i] <- sqrt(sum(effectivemuvec^2)) * invgram$lambda_node[i]
    }
  }

  p.vals <- 2 * pnorm(-abs(unbiased.Lasso) / stder)

  if (intercept1 == TRUE) {
    betaint_deb <- lasso_result$y_mean - as.numeric(lasso_result$x_mean %*% unbiased.Lasso)
    unb.coef_out <- c(betaint_deb, unbiased.Lasso)
    stder_out    <- c(NA, stder)
    low_out      <- c(NA, unbiased.Lasso - interval.sizes - addlength)
    up_out       <- c(NA, unbiased.Lasso + interval.sizes + addlength)
    pvals_out    <- c(NA, p.vals)
  } else {
    unb.coef_out <- unbiased.Lasso
    stder_out    <- stder
    low_out      <- unbiased.Lasso - interval.sizes - addlength
    up_out       <- unbiased.Lasso + interval.sizes + addlength
    pvals_out    <- p.vals
  }

  returnList <- list(
    "noise.sd" = stder_out,
    #"norm0" =  nz_hat,
    "coef" = htheta,
    "unb.coef" = unb.coef_out,
    "low.lim" = low_out,
    "up.lim" = up_out,
    "pvals" = pvals_out,
    "lambda_val" = lasso_result$lambda_val,
    "lambda_paper" = 2 * lasso_result$lambda_val,
    "low.lim.wald" = c(if (intercept1) NA_real_, low_wald),
    "up.lim.wald" = c(if (intercept1) NA_real_, up_wald),
    "low.lim.widened" = low_out, "up.lim.widened" = up_out,
    "addlength" = c(if (intercept1) NA_real_, addlength),
    "used_node" = TRUE, "inverse_method" = "full_nodewise",
    "node_kkt_gap" = max(abs(invgram$kkt_gap))
  )
  return(returnList)
}



# Max-type test of H0: delta_0 = 0 (no threshold effect), with a Gaussian
# multiplier bootstrap p-value -- implements T_n, T_n^*, and p_hat_B from the
# "Testing for a Threshold Effect" section of the paper.
#
# taugrid requirement (NOT checked at runtime -- caller's responsibility):
# every tau in taugrid must leave BOTH regimes {thdvar > tau} and
# {thdvar <= tau} non-empty, and no column of X1 may be identically zero on
# the in-regime rows {thdvar > tau}. Otherwise col.norm2 inside
# Threshold_Lasso_fixed_tau divides by zero (Inf), 0*Inf -> NaN feeds into
# glmnet, and that tau crashes the whole call. The tipping application uses
# q60_candidate_grid: nonempty groups and nonconstant full-design columns,
# without percentile trimming. Deliberately no extra defensive
# checks inside the tau loop (no in-regime-count check, no sigma2 floor):
# this function is meant to run inside Monte Carlo loops, where a per-tau
# check adds up over many replications. Get taugrid right at the call site.
#
# intercept1 (plain constant, tau-invariant) is safe to set TRUE or FALSE --
# it appears symmetrically in both the null fit (Linear_Lasso) and the
# alternative-hypothesis fit (Threshold_Lasso_fixed_tau).
#
# intercept2 (deltaint, the regime-jump coefficient on 1{Q_i>tau}) defaults
# to TRUE, matching every other intercept1/intercept2 default in this file
# (callers who want the paper's original no-intercept-jump model must pass
# intercept2 = FALSE explicitly -- never rely on the default here). H0 is
# "every threshold-effect coefficient is 0", and deltaint -- being
# multiplied by 1{Q_i>tau} exactly
# like delta -- IS a threshold-effect coefficient, so it belongs in that
# same null alongside delta, not treated as a free nuisance. Concretely:
# under H0 (deltaint_0 = 0 for every candidate tau, not just one), the model
# collapses to the SAME tau-free Y_i = X_i'beta_0 + U_i for every tau in the
# grid, so Linear_Lasso's null fit and the single shared U0 below are still
# exactly correct -- no tau-varying null model is needed. (This is different
# from treating deltaint as a nuisance ALLOWED to be nonzero under H0, which
# WOULD require a tau-varying null and is NOT what setting intercept2=TRUE
# does here -- deltaint is tested, not exempted.) When intercept2=TRUE,
# deltaint's index is prepended to test_idx below and nodewise-debiased
# exactly like every delta component -- it is a genuine, non-degenerate 0/1
# column (unlike betaint's literal constant column), so this needs no
# special-casing.
#
# Note on Theta scaling: nodewise_rows() returns the *undivided* T2inv %*% C
# (see estimate_nodewise_covariance's comment "inverse of t(Y) %*% Y" for
# Theta/n), so the raw nodewise_rows() output is already on the same scale as
# solve(t(X) %*% X / n) -- i.e. it IS Theta-hat(tau) as defined in the paper,
# with no further rescaling needed.
#
# thdvar must be supplied explicitly (Threshold_Lasso instead reads it from
# the caller's environment; that convention doesn't work here since the tau
# loop below may run under parallel::mclapply).
#
# Does NOT call set.seed() internally -- the bootstrap draws would otherwise
# be silently pinned on every call, which defeats Monte Carlo usage. Callers
# should set.seed() once, exactly as SimulationSyn.R already does.
threshold_effect_test <- function(X1,
                                   y,
                                   taugrid,
                                   thdvar,
                                   lambda = NULL,
                                   lambda0 = NULL,
                                   ic = 'GIC',
                                   B = 500,
                                   alpha = 0.05,
                                   intercept1 = TRUE,
                                   intercept2 = TRUE,
                                   nfolds = 10,
                                   ncores = 1,
                                   verbose = TRUE,
                                   candidate_profile = NULL) {
  n <- nrow(X1)
  y <- as.numeric(y)
  K_n <- length(taugrid)
  if (!is.null(candidate_profile)) {
    assert_gt_threshold_object(candidate_profile, "Global-test candidate profile")
    stopifnot(identical(candidate_profile$X1, X1),
      identical(candidate_profile$y, y), identical(candidate_profile$thdvar, thdvar),
      identical(candidate_profile$taugrid, taugrid),
      identical(candidate_profile$lambda, lambda),
      identical(candidate_profile$intercept1, intercept1),
      identical(candidate_profile$intercept2, intercept2))
  }

  # ---- Null-restricted fit (imposes H0) and its residuals ----
  null_fit <- Linear_Lasso(X1, y, lambda = lambda0, intercept1 = intercept1, col.norm1 = 1, nfolds = nfolds)
  null_coef <- t(null_fit$coefficients)   # full null-model vector [betaint?, beta], not just the intercept
  X11_null <- null_fit$X1
  U0 <- as.numeric(y - X11_null %*% null_coef)

  # ---- Fixed-tau debiased delta estimate + null score, for one tau ----
  fit_one_tau <- function(cnt) {
    tau <- taugrid[cnt]
    if (verbose && ncores == 1) {
      cat(format(Sys.time()), "- threshold_effect_test: tau", cnt, "/", K_n,
          "(", tau, ")\n")
    }

    fit_tau <- if (is.null(candidate_profile)) Threshold_Lasso_fixed_tau(
      X1, y, tau = tau, thdvar = thdvar, lambda = lambda,
      intercept1 = intercept1, intercept2 = intercept2, nfolds = nfolds
    ) else anchor_cached_candidate(candidate_profile, cnt)

    alpha_hat <- as.numeric(fit_tau$coefficients)
    Xtau <- cbind(fit_tau$X1, fit_tau$X_tau)
    p1 <- fit_tau$p1
    p2 <- fit_tau$p2
    delta_idx <- (p1 + 1):(p1 + p2)
    # deltaint (gamma) sits right after betaint (if present) in alpha_hat's
    # [betaint?, deltaint?, beta, delta] order -- see Threshold_Lasso_fixed_tau.
    # Prepend it to the tested set when intercept2=TRUE: it's a
    # threshold-effect coefficient (multiplied by 1{Q>tau}) exactly like
    # delta, tested under the SAME H0, not a nuisance -- see the comment
    # above this function for why this doesn't need a tau-varying null.
    deltaint_idx <- if (intercept2 == TRUE) as.integer(intercept1) + 1 else NULL
    test_idx <- c(deltaint_idx, delta_idx)
    p_test <- length(test_idx)

    # Remove and partial out the unpenalized baseline intercept. A literal
    # constant is not a nodewise predictor. Adjust indices after its removal.
    Xscore <- if (intercept1) Xtau[, -1, drop = FALSE] else Xtau
    if (intercept1) Xscore <- sweep(Xscore, 2, colMeans(Xscore), "-")
    score_idx <- test_idx - as.integer(intercept1)
    nw <- nodewise_rows(Xscore, rows = score_idx, ic = ic,
                        context = list(stage = "global_test", tau = tau))
    Theta_test <- nw$Theta

    resid_tau <- as.numeric(y - Xtau %*% alpha_hat)
    test_hat <- alpha_hat[test_idx]
    xtu <- as.numeric(crossprod(Xscore, resid_tau))
    a_test_hat <- as.numeric(test_hat + (Theta_test %*% xtu) / n)

    Mmat <- (Xscore %*% t(Theta_test)) * U0
    mbar <- colMeans(Mmat)
    sigma2 <- colMeans((Mmat - matrix(mbar, n, p_test, byrow = TRUE))^2)
    if (any(!is.finite(sigma2)) || any(sigma2 <= 0))
      stop("Nonpositive/nonfinite threshold-test score variance at tau = ", tau)
    sigma <- sqrt(sigma2)

    Psi <- sweep(sweep(Mmat, 2, mbar, "-"), 2, sigma, "/")
    Z <- sqrt(n) * a_test_hat / sigma

    list(Z = Z, Psi = Psi)
  }

  if (ncores > 1 && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(1:K_n, fit_one_tau, mc.cores = ncores)
  } else {
    results <- lapply(1:K_n, fit_one_tau)
  }

  var_names <- if (!is.null(colnames(X1))) colnames(X1) else paste0("X", 1:ncol(X1))
  if (intercept2 == TRUE) var_names <- c("deltaint", var_names)   # extra tested column
  Z_mat <- do.call(rbind, lapply(results, `[[`, "Z"))
  rownames(Z_mat) <- taugrid
  colnames(Z_mat) <- var_names
  Tn <- max(abs(Z_mat))

  Psi_all <- do.call(cbind, lapply(results, `[[`, "Psi"))  # n x (p_test * K_n)

  # ---- Gaussian multiplier bootstrap ----
  Xi <- matrix(rnorm(n * B), n, B)
  Tstar_mat <- crossprod(Psi_all, Xi) / sqrt(n)
  Tstar <- apply(abs(Tstar_mat), 2, max)
  pval <- (1 + sum(Tstar >= Tn)) / (B + 1)

  return(list(
    Tn = Tn,
    pval = pval,
    reject = pval <= alpha,
    Z = Z_mat,
    taugrid = taugrid,
    boot_stats = Tstar,
    alpha = alpha,
    threshold_rule = DTLASSO_THRESHOLD_RULE
  ))
}

# Build a candidate taugrid from the ACTUAL OBSERVED VALUES of the threshold
# variable, instead of an arbitrary evenly-spaced seq() grid -- same
# convention growth_gdp_DTLasso.R already hand-codes inline
# (sorted.thres/trim.factor), pulled out here as a reusable function.
#
# Two reasons this is what you want:
#  1. Standard practice in the threshold-regression literature (Hansen-style
#     sup-tests): if a true threshold exists, it coincides with SOME observed
#     value of thdvar, so candidate taus are conventionally the (trimmed)
#     order statistics of thdvar itself, not an arbitrary grid.
#  2. More efficient, not just more standard: any two candidate taus that
#     fall strictly between the SAME pair of adjacent observed thdvar values
#     induce the IDENTICAL sample split (1{thdvar>tau} doesn't change), so
#     every seq()-grid point that isn't itself an observed value is wasted
#     computation -- an identical Threshold_Lasso_fixed_tau/glmnet fit gets
#     redone for nothing. Using the observed values (deduplicated) means
#     every candidate tau in the grid induces a genuinely different split.
#
# trim: fraction trimmed off EACH end (default 0.1 = 10%, matching
# growth_gdp_DTLasso.R and the guidance in threshold_effect_test's header
# comment) so every candidate split leaves both regimes with enough
# observations -- never use the raw min/max as grid endpoints (the extreme
# tau's leave one regime with too few points, or none).
#
# IMPORTANT: the trim boundary is computed as a quantile of thdvar WITH
# duplicates counted (quantile(thdvar, ...) on the raw vector) -- NOT as a
# fraction of the number of UNIQUE values. If thdvar has heavy ties (e.g.
# 90% of observations sit at the same value), those two disagree sharply:
# trimming by unique-value count can throw away a value that represents
# most of the sample, while trimming by observation-weighted quantile
# correctly reflects how much of the actual data sits where. Deduplication
# (for the efficiency reason above) happens AFTER the quantile boundary is
# fixed, only to thin out which candidate taus get tested inside that
# boundary -- it never affects WHERE the boundary itself is set.
#
# digits (default NULL = off): round thdvar to this many decimal places
# BEFORE dedup, to deliberately coarsen resolution for continuous thdvar
# with no natural ties (e.g. runif draws) -- rounding manufactures ties, so
# a bare unique() collapses many nearby values into one candidate. This is
# still data-driven (candidates come from the data's own values, just at
# lower precision), unlike an arbitrary seq() grid picked independently of
# thdvar's actual range/density. The trim boundary (lo/hi) is still computed
# on the UNROUNDED thdvar first, for the same reason given above.
build_taugrid <- function(thdvar, trim = 0.1, digits = NULL) {
  lo <- stats::quantile(thdvar, trim, names = FALSE, type = 7)
  hi <- stats::quantile(thdvar, 1 - trim, names = FALSE, type = 7)
  vals <- if (is.null(digits)) thdvar else round(thdvar, digits)
  vals <- sort(unique(vals))
  out <- vals[vals >= lo & vals <= hi]
  if (length(out) == 0) stop("trim too large: leaves no candidate taus")
  out
}

getBiCop <- function(n,
                     rho,
                     mar.fun = rt,
                     x = NULL,
                     ...) {
  if (!is.null(x) &
      length(x) != n)
    warning("Variable x does not have the same length as n!")
  C <- matrix(rho, nrow = 2, ncol = 2)
  diag(C) <- 1
  
  C <- chol(C)
  
  X1 <- runif(n, 0, 1)
  X2 <- mar.fun(n, 10)
  X <- cbind(X1, X2)
  
  # induce correlation (does not change X1)
  df <- X %*% C
  
  ## if desired: check results
  #all.equal(X1,X[,1])
  #cor(X)
  
  return(df)
}
