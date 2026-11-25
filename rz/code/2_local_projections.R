# Ratio of sums of own debiased high-dimensional LP responses, on one common
# origin sample at every h. This is not a direct cumulative-outcome Lasso/IV fit.
# Usage: Rscript code/2_local_projections.R [--verify] [--case=news_anchor|rz_fixed]
entry <- grep("^--file=", commandArgs(FALSE), value = TRUE)
ext <- dirname(normalizePath(sub("^--file=", "", entry)))
source(file.path(ext, "common.R"))
source(file.path(ext, "multiplier_inference.R"))
ctx <- extension_setup(ext); extension_install_lp_recovery()
sample <- ctx$sample; dc <- ctx$bundle$data
news_path <- file.path(dirname(ext), "results/threshold/news/selection.rds")
news <- readRDS(news_path)
stopifnot(unname(news$key$data_md5) == ctx$data_md5, news$result$converged)
cases <- data.frame(case=c("news_anchor","rz_fixed"), shock="news",
 tau=c(100*news$result$selected$tau,6.5), stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly = TRUE)
verify_only <- "--verify" %in% args
filter <- sub("^--case=", "", grep("^--case=", args, value = TRUE))
if (length(filter)) {
  stopifnot(length(filter) == 1L, filter %in% cases$case)
  cases <- cases[cases$case == filter, , drop = FALSE]
}
stopifnot(all(is.finite(cases$tau)))
out <- file.path(dirname(ext), "results/lp")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
hmax <- settings$hmax
rows <- seq_len(nrow(sample$X1) - hmax)
origins <- sample$rows[rows]
n <- length(origins)
stopifnot(n == 480L, all(diff(dc$id[origins]) == 1L), max(origins + hmax) == nrow(dc))
cfg <- do.call(own_lp_control, c(lp_settings, list(seed = 100333L)))
hash <- unname(tools::md5sum(c(ctx$sources, file.path(ext, c("multiplier_inference.R", "2_local_projections.R")))))
irfs <- ratios <- differences <- checks <- inference_checks <- support <- node_checks <- contrasts <- list()

for (cc in seq_len(nrow(cases))) {
  case <- cases$case[cc]; shock_name <- cases$shock[cc]; tau <- cases$tau[cc]
  high <- as.numeric(sample$q[rows] >= tau)
  S <- cbind("Low unemployment" = 1 - high, "High unemployment" = high)
  shock <- dc$newsy[origins]
  Z <- sample$X1[rows, -1L, drop = FALSE]
  X <- cbind(shock_low = shock * (1 - high), shock_high = shock * high,
    state_low = 1 - high, Z * (1 - high), Z * high)
  colnames(X) <- c("shock_low", "shock_high", "state_low",
    paste0(colnames(Z), "_low"), paste0(colnames(Z), "_high"))
  zero <- which(colSums(abs(X)) == 0)
  stopifnot(!any(1:3 %in% zero))
  dropped <- colnames(X)[zero]
  if (length(zero)) X <- X[, -zero, drop = FALSE]
  H <- 1:3
  stopifnot(qr(cbind(1, X[, H]))$rank == 4L, all(colSums(S) > 4L),
    all(colSums(sweep(X, 2L, colMeans(X), "-")^2) > 0))
  threshold_path <- news_path
  key <- list(source_md5 = hash, data_md5 = ctx$data_md5, case = case, shock = shock_name,
    tau = tau, anchor = 6.5, lags = 4L, hmax = hmax, origins = origins,
    threshold_md5 = unname(tools::md5sum(threshold_path)), config = cfg,
    controls = ctx$bundle$controls, dropped_zero_columns = dropped,
    method = "common-sample sum of debiased horizon-specific LP responses")
  directory <- file.path(out, case)
  dir.create(directory, showWarnings = FALSE)
  nodepath <- file.path(directory, "nodes.rds")
  if (file.exists(nodepath)) {
    cached <- readRDS(nodepath)
    stopifnot(identical(cached$key, key), identical(cached$X, X))
    nodes <- cached$nodes
  } else {
    if (verify_only) stop("Missing common-sample nodes: ", case)
    anchor_log(case, "common-sample nodes", n, "tau", tau)
    nodes <- own_nodewise(X, H, cfg)
    te_atomic_save(list(key = key, X = X, nodes = nodes), nodepath)
  }
  stopifnot(all(vapply(nodes$fits, `[[`, logical(1), "converged")))
  node_checks[[case]] <- do.call(rbind, lapply(nodes$fits, function(f) data.frame(case = case,
    target = f$target, lambda = f$lambda, KKT = f$kkt, converged = f$converged,
    fixed_point_error = tail(f$trace$relative_change, 1))))
  support[[case]] <- data.frame(case = case, shock = shock_name, tau = tau,
    regime = colnames(S), n = colSums(S), nonzero_raw_shock = colSums(S * (shock != 0)),
    first_origin = min(dc$quarter[origins]), last_origin = max(dc$quarter[origins]))
  fitted <- list()
  xc <- sweep(X, 2L, colMeans(X), "-"); V <- xc %*% t(nodes$Gamma)
  for (v in c("y", "g")) {
    path <- file.path(directory, paste0("lp_", v, ".rds"))
    if (file.exists(path)) {
      saved <- readRDS(path); stopifnot(identical(saved$key, key), saved$outcome == v)
    } else saved <- list(key = key, outcome = v, horizons = vector("list", hmax + 1L))
    for (h in 0:hmax) {
      y <- dc[[v]][origins + h]
      if (is.null(saved$horizons[[h + 1L]])) {
        if (verify_only) stop("Incomplete LP: ", case, v, h)
        f <- own_fit_horizon(X, y, nodes, H, cfg,
          cfg$seed + 1000L + h + if (v == "y") 10000L else 0L)
        f$origin_rows <- origins; f$target_rows <- origins + h
        saved$horizons[[h + 1L]] <- f
        te_atomic_save(saved, path)
        anchor_log(case, v, "h", h, "lambda", signif(f$lasso$lambda, 5), file = file.path(directory, "lp.log"))
      }
      f <- saved$horizons[[h + 1L]]
      correction <- as.numeric(crossprod(V, f$residual)) / n / nodes$tau2
      covariance <- extension_hac_check(V * f$residual, f$HAC_Q) / outer(nodes$tau2, nodes$tau2) / n
      error <- max(abs(c(y - f$intercept - as.numeric(X %*% f$coefficients) - f$residual,
        f$coefficients[H] + correction - f$estimate, correction - f$correction,
        covariance - f$covariance, f$se - sqrt(pmax(diag(covariance), 0)))))
      stopifnot(error < 1e-8, f$lasso$converged, f$lasso$kkt < 1e-4,
        identical(f$origin_rows, origins), identical(f$target_rows, origins + h),
        all(dc$id[origins + h] - dc$id[origins] == h))
      row <- paste(case, v, h)
      irfs[[row]] <- data.frame(case = case, shock = shock_name, outcome = v, horizon = h,
        tau = tau, n = n, regime = colnames(S), estimate = f$estimate[1:2], se = f$se[1:2],
        lower = f$lower[1:2], upper = f$upper[1:2])
      contrast <- c(-1,1,0)
      difference <- sum(contrast*f$estimate)
      variance <- as.numeric(t(contrast) %*% covariance %*% contrast)
      stopifnot(variance > -1e-10)
      sdifference <- sqrt(max(variance,0))
      contrasts[[row]] <- data.frame(case=case,outcome=v,horizon=h,tau=tau,n=n,
        estimate=difference,se=sdifference,lower=difference-qnorm(.975)*sdifference,
        upper=difference+qnorm(.975)*sdifference,
        p_value=if(sdifference>0) 2*pnorm(-abs(difference)/sdifference) else NA_real_)
      checks[[row]] <- data.frame(case = case, outcome = v, horizon = h, n = n,
        max_error = error, KKT = f$lasso$kkt, mode = f$lasso$tuning_mode,
        converged = f$lasso$converged, fixed_point_error = tail(f$lasso$trace$relative_change, 1))
    }
    fitted[[v]] <- saved$horizons
  }
  # A single set of origins and nodes makes the score of a sum equal the sum of
  # the scores. Its HAC contains cross-horizon, cross-outcome and cross-state terms.
  score_h <- lapply(0:hmax, function(h) {
    w <- sweep(V[, 1:2, drop = FALSE], 2L, nodes$tau2[1:2], "/")
    cbind(w * fitted$y[[h + 1L]]$residual, w * fitted$g[[h + 1L]]$residual)
  })
  joint <- vector("list", hmax + 1L)
  for (h in 0:hmax) {
    ix <- seq_len(h + 1L)
    fy <- fitted$y[ix]; fg <- fitted$g[ix]
    beta <- c(Reduce(`+`, lapply(fy, function(f) f$estimate[1:2])),
              Reduce(`+`, lapply(fg, function(f) f$estimate[1:2])))
    W <- Reduce(`+`, score_h[ix])
    uy <- Reduce(`+`, lapply(fy, `[[`, "residual"))
    ug <- Reduce(`+`, lapply(fg, `[[`, "residual"))
    cy <- rowSums(vapply(0:h, function(j) dc$y[origins + j], numeric(n)))
    cg <- rowSums(vapply(0:h, function(j) dc$g[origins + j], numeric(n)))
    # Scale invariance in the common multivariate bandwidth rule.
    standardized_v <- sweep(V[, 1:2, drop = FALSE], 2L, nodes$x_scale[1:2], "/")
    score_standard <- cbind(standardized_v * (uy / sqrt(mean((cy - mean(cy))^2))),
                            standardized_v * (ug / sqrt(mean((cg - mean(cg))^2))))
    Q <- own_andrews_Q(score_standard)
    covariance <- own_hac(W, Q)$covariance / n
    independent <- extension_hac_check(W, Q) / n
    # A second route: directly debias the sum of horizon-specific coefficient vectors.
    by <- Reduce(`+`, lapply(fy, `[[`, "coefficients"))
    bg <- Reduce(`+`, lapply(fg, `[[`, "coefficients"))
    ay <- sum(vapply(fy, `[[`, numeric(1), "intercept"))
    ag <- sum(vapply(fg, `[[`, numeric(1), "intercept"))
    beta_check <- c(by[1:2] + colMeans(V[, 1:2, drop = FALSE] * uy) / nodes$tau2[1:2],
                    bg[1:2] + colMeans(V[, 1:2, drop = FALSE] * ug) / nodes$tau2[1:2])
    stacked_error <- NA_real_
    if (h %in% c(0L, 7L, 15L, 20L)) {
      all_scores <- do.call(cbind, score_h[ix])
      A <- kronecker(matrix(1, h + 1L, 1L), diag(4))
      stacked <- t(A) %*% extension_hac_check(all_scores, Q) %*% A / n
      stacked_error <- max(abs(stacked - covariance))
      stopifnot(stacked_error < 1e-7)
    }
    inf <- mi_inference(beta, covariance)
    difference_score <- as.numeric(W %*% inf$difference_gradient)
    difference_variance <- extension_hac_check(matrix(difference_score, ncol = 1L), Q)[1, 1] / n
    error <- max(abs(c(independent - covariance, beta_check - beta,
      cy - ay - as.numeric(X %*% by) - uy, cg - ag - as.numeric(X %*% bg) - ug,
      inf$difference_se^2 - difference_variance)))
    stopifnot(error < 1e-7)
    for (s in 1:2) {
      ratios[[paste(case, h, s)]] <- data.frame(case = case, shock = shock_name, tau = tau,
        horizon = h, quarters = h + 1L, n = n, regime = colnames(S)[s],
        cumulative_y = beta[s], cumulative_g = beta[s + 2L], estimate = inf$ratio[s],
        se = inf$ratio_se[s], lower = inf$ratio_lower[s], upper = inf$ratio_upper[s],
        p_value = inf$ratio_p[s], HAC_Q = Q, method = "joint HAC delta")
    }
    differences[[paste(case, h)]] <- data.frame(case = case, shock = shock_name, tau = tau,
      horizon = h, quarters = h + 1L, n = n, estimate = inf$difference, se = inf$difference_se,
      lower = inf$difference_lower, upper = inf$difference_upper, p_value = inf$difference_p,
      HAC_Q = Q, method = "joint HAC delta")
    inference_checks[[paste(case, h)]] <- data.frame(case = case, horizon = h,
      max_error = error, stacked_hac_error = stacked_error,
      minimum_covariance_eigenvalue = min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values))
    joint[[h + 1L]] <- list(beta = beta, covariance = covariance, scores = W, Q = Q, inference = inf)
  }
  if (!verify_only) te_atomic_save(list(key = key, joint = joint), file.path(directory, "joint_inference.rds"))
  anchor_log("Completed common-sample multipliers:", case, "tau", tau)
}
suffix <- if (length(filter)) paste0("_", filter) else ""
export <- function(x, name) write.csv(do.call(rbind, x), file.path(out, paste0(name, suffix, ".csv")), row.names = FALSE)
export(contrasts,"irf_differences");
export(irfs, "common_sample_irf"); export(ratios, "multipliers"); export(differences, "differences")
export(checks, "lp_verification"); export(inference_checks, "inference_verification")
export(support, "sample_counts"); export(node_checks, "node_verification")
write.csv(cases, file.path(out, paste0("cases", suffix, ".csv")), row.names = FALSE)
table <- do.call(rbind, ratios); diff <- do.call(rbind, differences)
table <- table[table$horizon %in% c(7L, 15L), ]
table$years <- (table$horizon + 1L) / 4
names(diff)[names(diff) == "estimate"] <- "high_minus_low"
names(diff)[names(diff) %in% c("p_value", "se", "lower", "upper")] <-
  paste0("difference_", names(diff)[names(diff) %in% c("p_value", "se", "lower", "upper")])
table <- merge(table, diff[, c("case", "horizon", "high_minus_low", "difference_p_value", "difference_se", "difference_lower", "difference_upper")],
  by = c("case", "horizon"), sort = FALSE)
write.csv(table, file.path(out, paste0("two_four_year_table", suffix, ".csv")), row.names = FALSE)
cat("PASS:", length(checks), "LP estimates;", length(inference_checks),
  "joint cumulative inference checks;", n, "common origins at every horizon.\n")
