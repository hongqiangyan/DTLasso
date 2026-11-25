# Shared-score plug-in calibration for the accepted anchor criterion.
tp_settings <- function() list(constant = 0.4, alpha = 0.05, draws = 1000L,
  maxit = 60L, tolerance = 0.01, seed = 20261002L, workers = 4L)

tp_prepare <- function(sample, cfg) {
  taus <- as.numeric(anchor_tau_grid(sample$q, sample$tau_star))
  designs <- lapply(taus, function(tau) anchor_design(sample$X1, sample$y, sample$q, tau))
  solvers <- lapply(designs, function(d) own_lasso_solver(d$X, d$y, d$unpenalized,
    own_lp_control(glmnet_tolerance = cfg$glmnet_tolerance, glmnet_maxit = cfg$glmnet_maxit)))
  # Exactly the projected, standardized penalized columns used by the solver;
  # do not standardize a second time after projecting out unpenalized columns.
  Z <- do.call(cbind, lapply(designs, function(d)
    qr.resid(qr(d$X[, d$unpenalized, drop = FALSE]), d$X[, d$penalized, drop = FALSE])))
  list(taus = taus, designs = designs, solvers = solvers, Z = Z, y = designs[[1L]]$y)
}

tp_scores <- function(Z, residual) {
  M <- Z * residual
  sweep(M, 2L, colMeans(M), "-")
}

tp_score_bandwidth <- function(M) {
  scales <- sqrt(colMeans(M^2))
  good <- scales > 1e-14
  if (!any(good)) stop("All calibration scores degenerate")
  max(2L, own_andrews_Q(sweep(M[, good, drop = FALSE], 2L, scales[good], "/")))
}

tp_draw_maxima <- function(M, ell, settings, seed) {
  # xi is shared across all coefficients and threshold candidates.
  xi <- te_multipliers(nrow(M), settings$draws, ell, seed)
  groups <- split(seq_len(ncol(M)), rep(seq_len(min(settings$workers, ncol(M))),
                                        length.out = ncol(M)))
  one <- function(ix) apply(abs(crossprod(M[, ix, drop = FALSE], xi)), 2L, max) / nrow(M)
  parts <- if (settings$workers > 1L) parallel::mclapply(groups, one,
    mc.cores = settings$workers, mc.preschedule = TRUE, mc.set.seed = FALSE) else lapply(groups, one)
  if (any(vapply(parts, inherits, logical(1), "try-error"))) stop("Score simulation worker failed")
  Reduce(pmax, parts)
}

tp_quantile <- function(maxima, settings) {
  # Same order-statistic convention as existing HDLP plug-in code.
  sort(maxima)[min(length(maxima), floor(length(maxima) * (1 - settings$alpha)) + 1L)]
}
