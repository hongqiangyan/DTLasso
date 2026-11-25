# Numerical recovery only after the existing nodewise fixed-bandwidth retry
# fails. Same residual-dependent lambda equation, normal draws, eigenvector
# convention, bandwidth, and tolerance. Damping halves the update; it does
# NOT halve the fixed-point discrepancy used to accept convergence.

tn_damped_plugin <- function(X, y, unpenalized=integer(), cfg=own_lp_control(),
                             seed=cfg$seed) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X); p <- ncol(X)
  solve_at <- own_lasso_solver(X,y,unpenalized,cfg)
  if (cfg$pi_constant == 0) {
    fit <- solve_at(0); fit$trace <- data.frame(iteration=1L,lambda=0,
      relative_change=0,Q=NA_integer_,cutoff=0,negative_eigenvalues=0L)
    fit$converged <- TRUE; return(fit)
  }
  normals <- own_normal_draws(p,cfg$pi_draws,seed)
  residual <- y-mean(y)
  old_lambda <- max(abs(crossprod(X,y)))/n
  trace <- list(); converged <- FALSE
  for (k in seq_len(cfg$pi_maxit)) {
    lrv <- own_hac(X*residual,cfg$hac_Q)
    eig <- eigen(lrv$covariance,symmetric=TRUE)
    # An eigenvector's sign is arbitrary. Canonicalize it so roundoff from
    # changing variable units does not reflect individual Monte Carlo draws.
    pivots <- max.col(t(abs(eig$vectors)),ties.method="first")
    signs <- sign(eig$vectors[cbind(pivots,seq_len(p))]); signs[signs==0] <- 1
    eig$vectors <- sweep(eig$vectors,2L,signs,"*")
    # Source uses eigenvalue truncation to draw from a PSD approximation.
    root <- sweep(eig$vectors,2L,sqrt(pmax(eig$values,0)),"*")
    draws <- root%*%normals
    maxima <- sort(apply(abs(draws),2L,max))
    # Match the source's zero-based order statistic (951st of 1000 at alpha=.05).
    order <- min(cfg$pi_draws,floor(cfg$pi_draws*(1-cfg$pi_alpha))+1L)
    cutoff <- maxima[order]
    target_lambda <- cfg$pi_constant*cutoff/sqrt(n)
    # Acceptance uses the ORIGINAL fixed-point equation, not the halved step.
    change <- abs(target_lambda-old_lambda)/max(old_lambda,.Machine$double.eps)
    lambda <- if (change < cfg$pi_tolerance) old_lambda else
      0.5 * (old_lambda + target_lambda)
    fit <- solve_at(lambda)
    trace[[k]] <- data.frame(iteration=k,lambda=lambda,relative_change=change, target_lambda=target_lambda,
      Q=lrv$Q,cutoff=cutoff,negative_eigenvalues=sum(eig$values < -1e-10))
    residual <- fit$residual
    if (change < cfg$pi_tolerance) { converged <- TRUE; break }
    old_lambda <- lambda
  }
  fit$trace <- do.call(rbind,trace); fit$converged <- converged
  fit$seed <- seed; fit$numerical_recovery <- "half-step iteration, original fixed-point criterion"; fit
}


# Same finite-draw fixed-point equation; direct bracketing is used only if
# both ordinary and half-step iterations fail. No tolerance is relaxed.
tn_bisect_plugin <- function(X, y, unpenalized, cfg, seed) {
  n <- nrow(X); p <- ncol(X)
  solve_at <- own_lasso_solver(X, y, unpenalized, cfg)
  normals <- own_normal_draws(p, cfg$pi_draws, seed)
  evaluate <- function(lambda) {
    fit <- solve_at(lambda)
    eig <- eigen(own_hac(X*fit$residual, cfg$hac_Q)$covariance, symmetric = TRUE)
    pivots <- max.col(t(abs(eig$vectors)), ties.method = "first")
    signs <- sign(eig$vectors[cbind(pivots,seq_len(p))]); signs[signs == 0] <- 1
    root <- sweep(sweep(eig$vectors,2L,signs,"*"),2L,sqrt(pmax(eig$values,0)),"*")
    maxima <- sort(apply(abs(root %*% normals),2L,max))
    cutoff <- maxima[min(cfg$pi_draws,floor(cfg$pi_draws*(1-cfg$pi_alpha))+1L)]
    target <- cfg$pi_constant*cutoff/sqrt(n)
    list(fit = fit, lambda = lambda, target = target, cutoff = cutoff,
      value = target-lambda, error = abs(target-lambda)/max(lambda,.Machine$double.eps),
      negative_eigenvalues = sum(eig$values < -1e-10))
  }
  hi <- max(abs(crossprod(X,y)))/n; lo <- hi*1e-5
  lower <- evaluate(lo); upper <- evaluate(hi)
  if (!(lower$value > 0 && upper$value < 0)) stop("No positive plug-in root bracket")
  trace <- list(); converged <- FALSE
  for (iteration in seq_len(80L)) {
    mid <- (lo+hi)/2; value <- evaluate(mid)
    trace[[iteration]] <- data.frame(iteration = iteration, lambda = mid,
      relative_change = value$error, target_lambda = value$target, Q = cfg$hac_Q,
      cutoff = value$cutoff, negative_eigenvalues = value$negative_eigenvalues)
    converged <- value$error < cfg$pi_tolerance
    if (converged) break
    if (value$value > 0) lo <- mid else hi <- mid
    if (abs(hi-lo) < 1e-14) break
  }
  fit <- value$fit; fit$converged <- converged
  fit$trace <- do.call(rbind, trace); fit$seed <- seed
  fit$numerical_recovery <- "bisection, original fixed-point criterion"
  fit
}

tn_install_recovery <- function() {
  original <- own_plugin_lasso
  recovered <- new.env(parent = emptyenv())
  wrapper <- function(X, y, unpenalized = integer(), cfg = own_lp_control(), seed = cfg$seed) {
    fit <- original(X, y, unpenalized, cfg, seed)
    if (!fit$converged && !is.null(cfg$hac_Q)) {
      retry_cfg <- cfg; retry_cfg$pi_maxit <- 300L
      retry <- tn_damped_plugin(X, y, unpenalized, retry_cfg, seed)
      recovery_method <- "half_step_after_fixed_bandwidth_failure"
      retry$trace$recovery <- recovery_method
      if (!retry$converged) {
        direct <- tn_bisect_plugin(X, y, unpenalized, retry_cfg, seed)
        recovery_method <- "bisection_after_half_step_failure"
        direct$trace$recovery <- recovery_method
        direct$trace$iteration <- direct$trace$iteration + nrow(retry$trace)
        direct$trace <- rbind(retry$trace, direct$trace)
        retry <- direct
      }
      attr(retry$trace, "original_iterations") <- nrow(fit$trace)
      # te_nodes stores trace row counts and the last Q, so expose the complete
      # numerical history while its final fit is the converged recovery result.
      fit$trace$target_lambda <- fit$trace$lambda
      fit$trace$recovery <- "original_fixed_bandwidth_iterations"
      retry$trace$iteration <- retry$trace$iteration + nrow(fit$trace)
      retry$trace <- rbind(fit$trace[, names(retry$trace)], retry$trace)
      fit <- retry
      assign(as.character(seed), list(error = tail(fit$trace$relative_change, 1L),
        method = recovery_method), envir = recovered)
      message("Nodewise numerical recovery: seed=", seed, ", Q=", cfg$hac_Q,
        ", converged=", fit$converged, ", total fixed-Q iterations=", nrow(fit$trace))
    }
    fit
  }
  assign("own_plugin_lasso", wrapper, envir = .GlobalEnv)
  original_nodes <- te_nodes
  node_wrapper <- function(X, H, cfg, retry_draws = 10000L) {
    result <- original_nodes(X, H, cfg, retry_draws)
    keys <- as.character(cfg$seed + H)
    used <- vapply(keys, exists, logical(1), envir = recovered, inherits = FALSE)
    result$diagnostics$numerical_recovery <- ifelse(used,
      vapply(keys, function(key) if (exists(key, envir = recovered, inherits = FALSE))
        get(key, envir = recovered)$method else "none", character(1)), "none")
    result$diagnostics$fixed_point_error_after_recovery <- vapply(keys, function(key)
      if (exists(key, envir = recovered, inherits = FALSE)) get(key, envir = recovered)$error else NA_real_, numeric(1))
    result
  }
  assign("te_nodes", node_wrapper, envir = .GlobalEnv)
  invisible(NULL)
}
