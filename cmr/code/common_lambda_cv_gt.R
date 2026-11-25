# Common-lambda 10-fold CV for the anchored threshold estimator.
# All penalties below use glmnet's loss scale; rho_n = 2*c_rho*lambda.
# This file is sourced after myfnc_anchor_gt.R.
# High-Q interactions use strict Q > tau, matching the fitting kernel.

anchor_raw_design <- function(X, q, tau, intercept2 = TRUE) {
  X2 <- X
  X2[q <= tau, ] <- 0
  if (intercept2) cbind(threshold_indicator_gt(q, tau), X, X2) else cbind(X, X2)
}

anchor_cv_design <- function(X, y, q, tau, intercept1, intercept2,
                             allow_fold_constant = FALSE) {
  W <- anchor_raw_design(X, q, tau, intercept2)
  means <- if (intercept1) colMeans(W) else rep(0, ncol(W))
  ymean <- if (intercept1) mean(y) else 0
  Z <- sweep(W, 2, means, "-")
  rms <- sqrt(colMeans(Z^2))
  if (any(!is.finite(rms))) stop("Nonfinite design scale at tau = ", tau)
  zero <- rms == 0
  if (any(zero) && !allow_fold_constant)
    stop("Constant full-sample design column at tau = ", tau)
  # A column constant only in a training fold has coefficient zero in that
  # fold. The baseline intercept absorbs its constant value. Keep positions
  # unchanged, so held-out predictions use exactly the same coefficient map.
  rms[zero] <- 1
  Z <- sweep(Z, 2, rms, "/")
  list(X = Z, y = y - ymean, means = means, ymean = ymean,
       rms = rms, zero = which(zero))
}

anchor_path <- function(z, lambda_grid) {
  fit <- glmnet::glmnet(z$X, z$y, lambda = lambda_grid,
                       standardize = FALSE, intercept = FALSE, maxit = 1000000L)
  if (fit$jerr != 0 || length(fit$lambda) != length(lambda_grid) ||
      !isTRUE(all.equal(as.numeric(fit$lambda), lambda_grid, tolerance = 1e-12)))
    stop("Incomplete common lambda path; glmnet jerr = ", fit$jerr)
  fit
}

anchor_common_cv <- function(X1, y, thdvar, taugrid, tau_star,
                             c_rho = 1, nfolds = 10L,
                             intercept1 = TRUE, intercept2 = TRUE,
                             nlambda = 100L, lambda_min_ratio = NULL,
                             lambda_grid = NULL, foldid = NULL) {
  n <- nrow(X1)
  y <- as.numeric(y)
  d <- 2L * ncol(X1) + as.integer(intercept2)
  stopifnot(nfolds >= 4L, nfolds <= n, length(y) == n,
            length(thdvar) == n, all(is.finite(X1)), all(is.finite(y)),
            all(is.finite(thdvar)), c_rho >= 0, is.finite(tau_star))
  taugrid <- sort(unique(c(taugrid, tau_star)))
  if (is.null(lambda_min_ratio)) lambda_min_ratio <- if (n < d) .01 else .0001
  if (is.null(lambda_grid)) {
    lambda_max <- 0
    for (tau in taugrid) {
      z <- anchor_cv_design(X1, y, thdvar, tau, intercept1, intercept2)
      lambda_max <- max(lambda_max, max(abs(crossprod(z$X, z$y))) / n)
    }
    if (!is.finite(lambda_max) || lambda_max <= 0) stop("Nonpositive lambda_max")
    lambda_grid <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
                           length.out = nlambda))
    lambda_grid[1] <- lambda_max
  }
  lambda_grid <- as.numeric(lambda_grid)
  stopifnot(length(lambda_grid) >= 2L, all(is.finite(lambda_grid)),
            all(lambda_grid > 0), all(diff(lambda_grid) < 0))
  L <- length(lambda_grid)
  if (is.null(foldid)) foldid <- sample(rep(seq_len(nfolds), length.out = n))
  stopifnot(length(foldid) == n, identical(sort(unique(as.integer(foldid))),seq_len(nfolds)))

  # Retain the requested grid endpoints even for a nearly saturated training
  # fit. This disables path truncation, without tightening fit tolerance or
  # refitting any selected penalty. Restore glmnet's process-local settings.
  old_control <- glmnet::glmnet.control()
  on.exit(do.call(glmnet::glmnet.control, old_control), add = TRUE)
  glmnet::glmnet.control(fdev = 0, devmax = 1)

  fold_sse <- fold_tau <- matrix(NA_real_, nfolds, L)
  fold_n <- tabulate(foldid, nbins = nfolds)
  fold_constant_count <- integer(nfolds)
  for (f in seq_len(nfolds)) {
    tr <- foldid != f; va <- !tr
    best_obj <- rep(Inf, L)
    for (tau in taugrid) {
      z <- anchor_cv_design(X1[tr,,drop=FALSE], y[tr], thdvar[tr], tau,
                            intercept1, intercept2, allow_fold_constant = TRUE)
      fold_constant_count[f] <- fold_constant_count[f] + length(z$zero)
      fit <- anchor_path(z, lambda_grid)
      beta <- as.matrix(fit$beta)
      residuals <- z$y - as.matrix(z$X %*% beta)
      obj <- colSums(residuals^2) / sum(tr) +
        2 * lambda_grid * colSums(abs(beta)) +
        2 * c_rho * lambda_grid * abs(tau - tau_star)
      # Ascending tau iteration: <= implements largest-tau exact ties.
      improve <- obj <= best_obj
      if (any(improve)) {
        Wva <- anchor_raw_design(X1[va,,drop=FALSE], thdvar[va], tau, intercept2)
        Zva <- sweep(sweep(Wva, 2, z$means, "-"), 2, z$rms, "/")
        prediction <- z$ymean + Zva %*% beta[,improve,drop=FALSE]
        fold_sse[f,improve] <- colSums((y[va] - prediction)^2)
        fold_tau[f,improve] <- tau
        best_obj[improve] <- obj[improve]
      }
    }
  }
  cv_mse <- colSums(fold_sse) / n
  stopifnot(all(is.finite(cv_mse)))
  selected <- which.min(cv_mse)  # exact lambda ties: larger penalty
  chosen <- lambda_grid[selected]

  # Final full-data fits, at the chosen shared lambda. Retain their path
  # coefficients so estimation and the global test use identical candidates.
  candidates <- vector("list", length(taugrid))
  for (i in seq_along(taugrid)) {
    z <- threshold_gic_design(X1, y, thdvar, taugrid[i], intercept1, intercept2)
    fit <- anchor_path(z, lambda_grid)
    b <- as.numeric(fit$beta[,selected])
    rss <- sum((z$y - as.numeric(z$X %*% b))^2)
    candidates[[i]] <- c(list(coefficients = b, lambda_val = chosen,
      rss = rss, Qn = rss/n + 2*chosen*sum(abs(b)), df = sum(b != 0),
      penalty_l1 = sum(abs(b))), z$scales)
  }
  table <- data.frame(tau = taugrid,
    Qn = vapply(candidates, `[[`, numeric(1), "Qn"),
    lambda_glmnet = chosen, lambda_paper = 2 * chosen,
    rss = vapply(candidates, `[[`, numeric(1), "rss"),
    df = vapply(candidates, `[[`, numeric(1), "df"),
    penalty_l1 = vapply(candidates, `[[`, numeric(1), "penalty_l1"))
  table$anchored_objective <- table$Qn + 2*c_rho*chosen*abs(taugrid-tau_star)
  cv <- list(table = data.frame(lambda_glmnet=lambda_grid,
      lambda_paper=2*lambda_grid, rho_n=2*c_rho*lambda_grid, cv_mse=cv_mse,
      selected=seq_len(L)==selected), selected_index=selected,
    foldid=foldid, fold_sse=fold_sse, fold_tau=fold_tau, fold_n=fold_n,
    training_constant_column_events=fold_constant_count,
    selected_at_boundary=selected %in% c(1L,L), nfolds=nfolds,
    grid_rule="global_kkt_log_grid", nlambda=L,
    lambda_min_ratio=tail(lambda_grid,1)/lambda_grid[1])
  profile <- list(X1=X1,y=y,thdvar=thdvar,taugrid=taugrid,
    threshold_rule=DTLASSO_THRESHOLD_RULE,
    candidates=candidates,lambda=chosen,intercept1=intercept1,intercept2=intercept2,
    tuning="common_lambda_anchored_cv",tau_star=tau_star,c_rho=c_rho,
    tau_tie_rule="largest_tau_exact",table=table,lambda_search=cv)
  list(profile=profile,lambda=chosen,cv=cv)
}

# Reconstruct one cached fixed-tau fit in the existing kernel's output order:
# [baseline intercept, threshold intercept, beta, delta], with optional intercepts.
anchor_cached_candidate <- function(profile, index) {
  assert_gt_threshold_object(profile, "Cached fixed-tau candidate")
  X <- profile$X1; q <- profile$thdvar; tau <- profile$taugrid[index]
  f <- profile$candidates[[index]]
  W <- anchor_raw_design(X, q, tau, profile$intercept2)
  b <- f$coefficients * c(f$col.norm1, f$col.norm2)
  b0 <- if (profile$intercept1) f$y_mean - sum(colMeans(W) * b) else NULL
  X11 <- if (profile$intercept2) cbind(threshold_indicator_gt(q, tau),X) else X
  if (profile$intercept1) X11 <- cbind(1,X11)
  X2 <- X; X2[q<=tau,] <- 0
  list(coefficients=matrix(c(b0,b),nrow=1),tau=tau,X1=X11,X_tau=X2,
       lambda_val=profile$lambda,p1=ncol(X11),p2=ncol(X),
       threshold_rule=DTLASSO_THRESHOLD_RULE)
}
