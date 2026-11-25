# Anchor-residual feedback for the existing time-series plug-in calibration.
# Conditional on each positive c_rho, update lambda, rho=c_rho*lambda, tau,
# coefficients and residuals. This does NOT estimate or rank c_rho.
# Unknown-threshold calibration and anchor feedback are experimental extensions
# of the HDLP/ASW plug-in idea, not an established joint-tuning theorem.
# Requires lp.R, data_threshold.R, threshold_test.R and threshold_plugin.R.

tap_engine <- function(sample, cfg, seed) {
  pr <- tp_prepare(sample, cfg)
  fit_cache <- new.env(parent = emptyenv())
  score_cache <- new.env(parent = emptyenv())
  number_key <- function(x) sprintf("%a", x)
  stopifnot(all(vapply(pr$designs, function(d) identical(d$unpenalized, 1:3), logical(1))))
  profile <- function(lambda, c_rho) {
    if (length(c_rho) != 1L || !is.finite(c_rho) || c_rho <= 0)
      stop("c_rho must be strictly positive")
    if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0)
      stop("lambda must be strictly positive")
    key <- number_key(lambda)
    if (!exists(key, fit_cache, inherits = FALSE)) {
      fits <- lapply(pr$solvers, function(solve_at) solve_at(lambda))
      table <- data.frame(tau = pr$taus,
        rss_half = vapply(fits, function(f) mean(f$residual^2) / 2, numeric(1)),
        lasso_penalty = vapply(fits, function(f)
          lambda * sum(abs(f$coefficients[f$penalized])), numeric(1)),
        kkt = vapply(fits, `[[`, numeric(1), "kkt"))
      stopifnot(max(table$kkt) < cfg$kkt_tolerance)
      assign(key, list(fits = fits, table = table), fit_cache)
    }
    cached <- get(key, fit_cache, inherits = FALSE)
    table <- cached$table
    table$anchor_penalty <- c_rho * lambda * abs(table$tau - sample$tau_star)
    table$objective_half <- table$rss_half + table$lasso_penalty + table$anchor_penalty
    k <- anchor_argmin_upper(table$objective_half, table$tau)
    list(index = k, tau = table$tau[k], lambda = lambda, c_rho = c_rho,
      rho = c_rho * lambda, fit = cached$fits[[k]], profile = table,
      # Conditional on tau and lambda, c_rho does not affect the coefficient fit.
      state_key = paste(k, key, sep = ":"))
  }
  score <- function(residual, state_key, use_cache = TRUE) {
    if (use_cache && exists(state_key, score_cache, inherits = FALSE))
      return(get(state_key, score_cache, inherits = FALSE))
    M <- tp_scores(pr$Z, residual)
    ell <- tp_score_bandwidth(M)
    maxima <- tp_draw_maxima(M, ell, cfg$plugin, seed)
    cutoff <- tp_quantile(maxima, cfg$plugin)
    value <- list(lambda = cfg$plugin$constant * cutoff, cutoff = cutoff,
      bandwidth = ell, maxima = maxima)
    if (use_cache) assign(state_key, value, score_cache)
    value
  }
  list(prepared = pr, profile = profile, score = score,
    initial_lambda = max(abs(crossprod(pr$Z, pr$y))) / length(pr$y),
    cache_sizes = function() c(profiles = length(ls(fit_cache)), scores = length(ls(score_cache))))
}

tap_relative <- function(x, previous) abs(x - previous) / max(abs(previous), .Machine$double.eps)
tap_residual_change <- function(x, previous)
  sqrt(sum((x - previous)^2) / max(sum(previous^2), .Machine$double.eps))

tap_calibrate <- function(engine, c_rho, cfg, logfile = NULL) {
  if (length(c_rho) != 1L || !is.finite(c_rho) || c_rho <= 0)
    stop("c_rho must be strictly positive")
  st <- cfg$plugin
  stopifnot(st$constant > 0, st$alpha > 0, st$alpha < 1,
    st$maxit >= 1L, st$tolerance > 0, cfg$residual_tolerance > 0,
    cfg$stable_iterations >= 1L)
  residual <- engine$prepared$y
  state_key <- "initial_standardized_y"
  lambda_old <- engine$initial_lambda
  previous <- NULL; previous_ell <- NA_integer_; stable <- 0L
  trace <- list(); status <- "max_iterations"; lookahead <- NULL; cycle_period <- NA_integer_
  for (iteration in seq_len(st$maxit)) {
    calibration <- engine$score(residual, state_key)
    current <- engine$profile(calibration$lambda, c_rho)
    dl <- tap_relative(current$lambda, lambda_old)
    du <- tap_residual_change(current$fit$residual, residual)
    same_tau <- !is.null(previous) && identical(current$tau, previous$tau)
    same_bandwidth <- identical(calibration$bandwidth, previous_ell)
    stable_step <- dl < st$tolerance && du < cfg$residual_tolerance && same_tau && same_bandwidth
    stable <- if (stable_step) stable + 1L else 0L
    values <- sort(current$profile$objective_half)
    trace[[iteration]] <- data.frame(c_rho = c_rho, iteration = iteration,
      lambda = current$lambda, rho = current$rho, tau = current$tau,
      input_tau = if (is.null(previous)) NA_real_ else previous$tau,
      input_lambda = if (is.null(previous)) NA_real_ else previous$lambda,
      cutoff = calibration$cutoff, bandwidth = calibration$bandwidth,
      lambda_relative_change = dl, residual_relative_change = du,
      same_tau = same_tau, same_bandwidth = same_bandwidth, stable_steps = stable,
      objective_half = values[1L], objective_gap = values[2L] - values[1L],
      max_KKT = max(current$profile$kkt))
    anchor_log("Anchor plug-in c_rho", c_rho, "iteration", iteration,
      "lambda", signif(current$lambda, 8), "rho", signif(current$rho, 7),
      "tau", current$tau, "bandwidth", calibration$bandwidth,
      "lambda change", signif(dl, 3), "residual change", signif(du, 3), file = logfile)
    if (stable >= cfg$stable_iterations) {
      next_score <- engine$score(current$fit$residual, current$state_key)
      next_fit <- engine$profile(next_score$lambda, c_rho)
      lookahead <- list(lambda = next_score$lambda, tau = next_fit$tau,
        bandwidth = next_score$bandwidth,
        lambda_relative_change = tap_relative(next_score$lambda, current$lambda),
        residual_relative_change = tap_residual_change(next_fit$fit$residual, current$fit$residual))
      if (lookahead$lambda_relative_change < st$tolerance &&
          lookahead$residual_relative_change < cfg$residual_tolerance &&
          identical(lookahead$tau, current$tau) &&
          identical(lookahead$bandwidth, calibration$bandwidth)) {
        status <- "converged"; break
      }
      stable <- 0L
    }
    # Detect a repeated nonconstant numerical cycle. Never freeze the bandwidth
    # or average alternating thresholds and label the result converged.
    history <- do.call(rbind, trace)
    if (iteration >= 4L) for (period in seq.int(2L, min(8L, floor(iteration / 2)))) {
      a <- seq.int(iteration - period + 1L, iteration); b <- a - period
      repeated <- all(history$tau[a] == history$tau[b]) &&
        all(history$bandwidth[a] == history$bandwidth[b]) &&
        max(abs(history$lambda[a] - history$lambda[b]) / history$lambda[b]) < 1e-10
      varying <- length(unique(history$tau[a])) > 1L ||
        length(unique(history$bandwidth[a])) > 1L ||
        diff(range(history$lambda[a])) / mean(history$lambda[a]) > st$tolerance
      if (repeated && varying) { cycle_period <- period; status <- "cycle"; break }
    }
    if (status == "cycle" || iteration == st$maxit) break
    residual <- current$fit$residual; state_key <- current$state_key
    previous <- current; previous_ell <- calibration$bandwidth; lambda_old <- current$lambda
  }
  anchor_log("Anchor plug-in c_rho", c_rho, status, "after", iteration, "iterations", file = logfile)
  list(c_rho = c_rho, status = status, converged = identical(status, "converged"),
    cycle_period = cycle_period, selected = current, calibration = calibration,
    score_residual = residual,
    score_input = if (is.null(previous)) NULL else
      list(tau = previous$tau, lambda = previous$lambda, coefficients = previous$fit$coefficients),
    trace = do.call(rbind, trace), lookahead = lookahead)
}

tap_verify <- function(sample, result, engine, cfg, replay = TRUE) {
  selected <- result$selected; tab <- selected$profile
  stopifnot(result$c_rho > 0, selected$rho == result$c_rho * selected$lambda,
    selected$lambda == cfg$plugin$constant * tp_quantile(result$calibration$maxima, cfg$plugin),
    identical(tab$tau, engine$prepared$taus), max(result$trace$max_KKT) < cfg$kkt_tolerance)
  Q <- tab$rss_half + tab$lasso_penalty
  objective <- Q + selected$rho * abs(tab$tau - sample$tau_star)
  decimal_objective <- Q + (100 * result$c_rho) * selected$lambda *
    abs(tab$tau / 100 - sample$tau_star / 100)
  stopifnot(max(abs(objective - tab$objective_half)) < 1e-12,
    max(abs(decimal_objective - objective)) < 1e-12,
    anchor_argmin_upper(objective, tab$tau) == selected$index,
    anchor_argmin_upper(2 * objective, tab$tau) == selected$index,
    anchor_argmin_upper(decimal_objective, tab$tau) == selected$index)
  # Fresh design and solver, independent of the profile cache.
  d <- anchor_design(sample$X1, sample$y, sample$q, selected$tau)
  f <- own_lasso_solver(d$X, d$y, d$unpenalized,
    own_lp_control(glmnet_tolerance = cfg$glmnet_tolerance, glmnet_maxit = cfg$glmnet_maxit))(selected$lambda)
  refit_error <- max(abs(f$residual - selected$fit$residual))
  b <- selected$fit$coefficients * d$y_scale / d$weights
  intercept <- d$y_mean - sum(d$means * b)
  raw_error <- max(abs(sample$y - intercept - as.numeric(d$raw %*% b) -
    selected$fit$residual * d$y_scale))
  if (is.null(result$score_input)) {
    reconstructed <- engine$prepared$y
  } else {
    di <- anchor_design(sample$X1, sample$y, sample$q, result$score_input$tau)
    reconstructed <- di$y - as.numeric(di$X %*% result$score_input$coefficients)
  }
  input_error <- max(abs(reconstructed - result$score_residual))
  stopifnot(refit_error < 1e-9, raw_error < 1e-8, input_error < 1e-9,
    identical(f$unpenalized, 1:3), f$kkt < cfg$kkt_tolerance)
  replay_error <- NA_real_
  if (replay) {
    check <- engine$score(result$score_residual, "independent_replay", use_cache = FALSE)
    replay_error <- max(abs(check$maxima - result$calibration$maxima))
    stopifnot(replay_error < 1e-12, check$bandwidth == result$calibration$bandwidth)
  }
  if (result$converged) {
    last <- tail(result$trace, 1L)
    stopifnot(last$stable_steps >= cfg$stable_iterations,
      last$lambda_relative_change < cfg$plugin$tolerance,
      last$residual_relative_change < cfg$residual_tolerance,
      result$lookahead$lambda_relative_change < cfg$plugin$tolerance,
      result$lookahead$residual_relative_change < cfg$residual_tolerance,
      result$lookahead$tau == selected$tau,
      result$lookahead$bandwidth == result$calibration$bandwidth)
  }
  list(numerical_checks_passed = TRUE, converged = result$converged,
    input_residual_error = input_error, raw_residual_error = raw_error,
    independent_refit_error = refit_error, multiplier_replay_error = replay_error,
    max_KKT = max(result$trace$max_KKT))
}
