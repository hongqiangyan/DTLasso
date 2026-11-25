te_integer <- function(x, lower, upper, name) {
  if (length(x) != 1L || !is.finite(x) || x != as.integer(x) ||
      x < lower || x > upper) stop("Invalid ", name)
  as.integer(x)
}

te_multipliers <- function(n, B, bandwidth, seed) {
  n <- te_integer(n, 2L, .Machine$integer.max, "sample size")
  B <- te_integer(B, 1L, .Machine$integer.max, "draw count")
  ell <- te_integer(bandwidth, 1L, n, "bandwidth")
  # Every moving window has ell innovations, including at both boundaries.
  # No wrapping, zero-padding, sample demeaning, or sample rescaling of xi.
  z <- own_normal_draws(n + ell - 1L, B, seed)
  cs <- rbind(0, matrix(apply(z, 2L, cumsum), nrow(z), B))
  (cs[seq_len(n) + ell, , drop = FALSE] -
     cs[seq_len(n), , drop = FALSE]) / sqrt(ell)
}

te_bandwidth <- function(centered_scores) {
  scale <- sqrt(colMeans(centered_scores^2))
  if (any(!is.finite(scale) | scale <= 1e-14)) stop("Degenerate score column")
  # Common, unit-invariant AR(1) plug-in bandwidth across ALL tau/coefficient
  # scores. Empirical tuning rule, not an optimality or uniform-validity claim.
  max(2L, own_andrews_Q(sweep(centered_scores, 2L, scale, "/")))
}

te_null_fit <- function(sample, lambda, cfg) {
  s <- own_standardize(sample$X1, sample$y)
  f <- own_lasso_solver(s$X, s$y, 1L, cfg)(lambda)
  b <- f$coefficients * s$y_scale / s$x_scale
  list(coefficients = b, intercept = s$y_mean - sum(s$x_mean * b),
       residual = f$residual * s$y_scale, lambda = lambda, kkt = f$kkt)
}

te_nodes <- function(X, H, cfg, retry_draws = 10000L) {
  s <- own_standardize(X, X[, 1L]); z <- s$X
  p <- ncol(X); n <- nrow(X)
  Theta <- matrix(0, length(H), p); diagnostics <- vector("list", length(H))
  for (i in seq_along(H)) {
    j <- H[i]; other <- setdiff(seq_len(p), j)
    local_cfg <- cfg
    fit <- own_plugin_lasso(z[, other, drop = FALSE], z[, j], integer(),
                            local_cfg, cfg$seed + j)
    tuning_mode <- "adaptive_bandwidth"
    if (!fit$converged) {
      local_cfg$pi_draws <- retry_draws
      fit <- own_plugin_lasso(z[, other, drop = FALSE], z[, j], integer(),
                              local_cfg, cfg$seed + j)
    }
    # Integer bandwidth updates can enter a deterministic two-cycle. Freeze
    # Q at the largest of its last 20 iterates, then solve the SAME Gaussian
    # plug-in lambda equation conditional on that Q. Record this fallback.
    if (!fit$converged && is.null(cfg$hac_Q)) {
      local_cfg$hac_Q <- max(tail(fit$trace$Q, 20L))
      tuning_mode <- "fixed_bandwidth_after_nonconvergence"
      fit <- own_plugin_lasso(z[, other, drop = FALSE], z[, j], integer(),
                              local_cfg, cfg$seed + j)
    }
    if (!fit$converged) stop("Nodewise tuning failed at ", colnames(X)[j])
    g <- fit$coefficients
    denom <- mean(fit$residual^2) + fit$lambda * sum(abs(g))
    if (!is.finite(denom) || denom <= 1e-12) stop("Degenerate nodewise fit")
    row <- numeric(p); row[j] <- 1; row[other] <- -g
    # Transform the precision row from standardized columns to raw units.
    Theta[i, ] <- row / denom / s$x_scale[j] / s$x_scale
    norm_error <- abs(mean(fit$residual * z[, j]) - denom)
    diagnostics[[i]] <- data.frame(coefficient = colnames(X)[j],
      lambda = fit$lambda, iterations = nrow(fit$trace),
      draws = local_cfg$pi_draws, converged = fit$converged,
      tuning_mode = tuning_mode,
      kkt = fit$kkt, normalization_error = norm_error,
      nonzero = sum(abs(g) > 1e-9),
      Q = tail(fit$trace$Q, 1L))
  }
  list(Theta = Theta, diagnostics = do.call(rbind, diagnostics))
}

te_fit_tau <- function(sample, tau, lambda, null_residual, cfg,
                       retry_draws = 10000L) {
  design <- anchor_design(sample$X1, sample$y, sample$q, tau)
  X <- design$raw; n <- nrow(X)
  # Test state intercept, shock difference, and EVERY lag difference. Baseline
  # slopes and the global intercept are nuisance parameters under H0.
  H <- c(3L, 2L, which(grepl("_difference$", colnames(X)) &
                        colnames(X) != "shock_difference"))
  f <- own_lasso_solver(design$X, design$y, design$unpenalized, cfg)(lambda)
  b <- f$coefficients * design$y_scale / design$weights
  xc <- sweep(X, 2L, design$means, "-")
  residual <- f$residual * design$y_scale
  nw <- te_nodes(X, H, cfg, retry_draws)
  influence_design <- xc %*% t(nw$Theta)
  estimate <- b[H] + as.numeric(crossprod(influence_design, residual)) / n
  scores <- influence_design * null_residual
  colnames(scores) <- colnames(X)[H]
  list(tau = tau, H = H, names = colnames(X)[H], estimate = estimate,
       initial = b[H], scores = scores, coefficients = b,
       intercept = design$y_mean - sum(design$means * b),
       residual = residual, Theta = nw$Theta,
       main_lambda = lambda, main_kkt = f$kkt,
       node_diagnostics = nw$diagnostics,
       n_low = sum(sample$q < tau), n_high = sum(sample$q >= tau))
}

te_atomic_save <- function(object, path) {
  temporary <- paste0(path, ".tmp")
  saveRDS(object, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) stop("Could not save ", path)
}

te_block_sums <- function(scores, block_length) {
  M <- as.matrix(scores); n <- nrow(M)
  ell <- te_integer(block_length, 1L, n - 1L, "block length")
  if (any(!is.finite(M))) stop("Nonfinite score")
  M <- sweep(M, 2L, colMeans(M), "-")
  block <- (seq_len(n) - 1L) %/% ell + 1L
  A <- rowsum(M, block, reorder = FALSE)
  # Retain a final short block. Both the estimates and scores then use all n
  # observations, rather than pairing full-sample estimates with truncated scores.
  omega <- sqrt(colSums(A^2) / n)
  if (any(!is.finite(omega) | omega <= 1e-10)) stop("Degenerate block variance")
  list(A = A, omega = omega, block = block, sizes = tabulate(block), n = n,
       block_length = ell)
}

te_block_bootstrap <- function(scores, estimate, block_length, B = 4999L,
                               alpha = .05, seed = 20261001L,
                               chunk_size = min(250L, B)) {
  a <- te_block_sums(scores, block_length)
  if (length(estimate) != ncol(a$A) || any(!is.finite(estimate)) ||
      length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("Invalid test inputs")
  B <- te_integer(B, 99L, .Machine$integer.max, "draw count")
  chunk_size <- te_integer(chunk_size, 1L, B, "chunk size")
  xi <- own_normal_draws(nrow(a$A), B, seed)
  standardized <- sweep(a$A, 2L, a$omega, "/")
  Z <- sqrt(a$n) * estimate / a$omega
  boot <- numeric(B)
  for (first in seq.int(1L, B, by = chunk_size)) {
    ix <- seq.int(first, min(B, first + chunk_size - 1L))
    z <- crossprod(standardized, xi[, ix, drop = FALSE]) / sqrt(a$n)
    boot[ix] <- apply(abs(z), 2L, max)
  }
  Tn <- max(abs(Z)); exceedances <- sum(boot >= Tn)
  pval <- (1 + exceedances) / (B + 1)
  list(Tn = Tn, pval = pval, reject = pval <= alpha, alpha = alpha,
    critical_value = unname(quantile(boot, 1 - alpha, type = 1)),
    B = B, exceedances = exceedances,
    p_mc_se = sqrt(pval * (1 - pval) / (B + 1)),
    block_length = a$block_length, block_sizes = a$sizes,
    block_count = nrow(a$A), n = a$n, seed = seed,
    omega = a$omega, Z = Z, boot_stats = boot)
}
