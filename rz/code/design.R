anchor_log <- function(..., file = NULL) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste(..., collapse = " "))
  cat(line, "\n")
  if (!is.null(file)) cat(line, "\n", file = file, append = TRUE)
  invisible(line)
}

anchor_argmin_upper <- function(score, tau) {
  if (length(score) != length(tau) || any(!is.finite(score))) stop("Invalid threshold scores.")
  tied <- which(score == min(score))
  tied[which.max(tau[tied])] # Exact ties only, no flatness tolerance.
}

anchor_design <- function(X1, y, q, tau) {
  X1 <- as.matrix(X1)
  if (is.null(colnames(X1))) colnames(X1) <- paste0("x",seq_len(ncol(X1)))
  low <- as.numeric(q < tau)
  Z <- X1[, -1L, drop = FALSE]
  # Preserve the original baseline + threshold-difference slope parametrization.
  # Both shock coefficients and the state intercept are unpenalized.
  # The state intercept aligns the conditional-mean model with desla::HDLP.
  raw <- cbind(shock = X1[, 1L], shock_difference = X1[, 1L] * low,
               state_low = low, Z, Z * low)
  colnames(raw) <- c("shock", "shock_difference", "state_low", colnames(Z),
                    paste0(colnames(Z), "_difference"))
  means <- colMeans(raw)
  centered <- sweep(raw, 2L, means, "-")
  weights <- sqrt(colMeans(centered^2))
  if (any(!is.finite(weights) | weights <= 0)) stop("Zero design weight at tau = ", tau)
  X <- sweep(centered, 2L, weights, "/")
  U <- qr(X[, 1:3, drop = FALSE])
  if (U$rank != 3L) stop("Rank-deficient unpenalized block at tau = ", tau)
  y_scale <- sqrt(mean((y - mean(y))^2))
  if (!is.finite(y_scale) || y_scale <= 0) stop("Zero/nonfinite response scale.")
  # Preserve desla's original scale=TRUE convention for the response as well.
  # One h=0 response scale is shared by every tau and lambda.
  list(X = X, y = (y - mean(y)) / y_scale, means = means, weights = weights,
       y_mean = mean(y), y_scale = y_scale,
       unpenalized = 1:3, penalized = seq.int(4L, ncol(X)),
       unpenalized_qr = U, raw = raw)
}

ec_search_sample <- function(dc,controls,lags=4L,outcome="y",anchor=6.5) {
  stopifnot(all(diff(dc$id)==1L),lags>=1L)
  lagged <- do.call(cbind,lapply(controls,function(v) {
    z <- vapply(seq_len(lags),function(k) c(rep(NA_real_,k),head(dc[[v]],-k)),numeric(nrow(dc)))
    colnames(z)<-paste0(v,"_lag_",seq_len(lags)); z
  }))
  qlag <- c(NA_real_,head(dc$unemp,-1L))
  rows <- seq.int(lags+1L,nrow(dc))
  list(X1=cbind(newsy=dc$newsy,lagged)[rows,,drop=FALSE],y=dc[[outcome]][rows],
    q=qlag[rows],rows=rows,quarter=dc$quarter[rows],qlag=qlag,tau_star=anchor)
}

ec_verify_alignment <- function(s,dc,controls,lags,outcome) {
  stopifnot(identical(s$y,dc[[outcome]][s$rows]),
    identical(s$q,dc$unemp[s$rows-1L]),
    identical(as.numeric(s$X1[,1L]),dc$newsy[s$rows]))
  expected <- unlist(lapply(controls,function(v)paste0(v,"_lag_",seq_len(lags))))
  expected <- c("newsy",expected)
  stopifnot(identical(colnames(s$X1),expected))
  for(v in controls) for(k in seq_len(lags)) {
    column <- paste0(v,"_lag_",k)
    stopifnot(identical(as.numeric(s$X1[,column]),dc[[v]][s$rows-k]))
  }
  invisible(TRUE)
}
