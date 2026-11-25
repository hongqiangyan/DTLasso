# Standalone high-dimensional LP building blocks. No desla/HDLP calls.
# Loss convention throughout: RSS/(2*n) + lambda * sum(abs(penalized slopes)).
# glmnet is only the numerical Lasso solver; tuning and inference are below.

own_lp_control <- function(...) {
  defaults <- list(pi_constant = 0.4, pi_alpha = 0.05, pi_draws = 1000L,
    pi_maxit = 15L, pi_tolerance = 0.01, glmnet_tolerance = 1e-12,
    glmnet_maxit = 1000000L, hac_Q = NULL, normalization = "paper",
    alpha = 0.05, seed = 333L, verbose = TRUE)
  x <- list(...)
  if (length(setdiff(names(x), names(defaults)))) stop("Unknown LP setting.")
  cfg <- utils::modifyList(defaults, x)
  if (cfg$normalization != "paper") stop("This version uses the paper normalization.")
  if (!is.finite(cfg$pi_constant) || cfg$pi_constant < 0 ||
      cfg$pi_alpha <= 0 || cfg$pi_alpha >= 1 || cfg$alpha <= 0 || cfg$alpha >= 1 ||
      cfg$pi_draws < 20 || cfg$pi_maxit < 1 || cfg$pi_tolerance <= 0)
    stop("Invalid LP tuning settings.")
  cfg
}

own_standardize <- function(X, y) {
  X <- as.matrix(X); y <- as.numeric(y)
  if (nrow(X) != length(y) || any(!is.finite(X)) || any(!is.finite(y)))
    stop("Nonfinite or misaligned design/response.")
  xm <- colMeans(X); ym <- mean(y)
  xc <- sweep(X, 2L, xm, "-"); yc <- y - ym
  xs <- sqrt(colMeans(xc^2)); ys <- sqrt(mean(yc^2))
  if (any(xs <= 1e-14) || ys <= 1e-14) stop("Constant design column or outcome.")
  list(X = sweep(xc, 2L, xs, "/"), y = yc / ys,
       x_mean = xm, x_scale = xs, y_mean = ym, y_scale = ys)
}

# AR(1) plug-in bandwidth for the Bartlett kernel, with the source's T/2 cap.
own_andrews_moments <- function(scores) {
  W <- as.matrix(scores); n <- nrow(W)
  if (n < 4L || any(!is.finite(W))) stop("Invalid HAC scores.")
  x <- W[-n, , drop = FALSE]; y <- W[-1L, , drop = FALSE]
  x <- sweep(x, 2L, colMeans(x), "-")
  y <- sweep(y, 2L, colMeans(y), "-")
  xx <- colSums(x^2)
  good <- xx > .Machine$double.eps
  rho <- numeric(ncol(W)); rho[good] <- colSums(x*y)[good] / xx[good]
  # Avoid explosive numerical bandwidths for fitted AR coefficients near +/-1.
  rho <- pmax(-0.99, pmin(0.99, rho))
  res <- y - sweep(x, 2L, rho, "*")
  variance <- colSums(res^2) / (n - 2L)
  den <- sum(variance^2 / (1-rho)^4)
  list(n=n, numerator=sum(4*rho^2*variance^2 / ((1-rho)^6*(1+rho)^2)),
       denominator=den)
}

own_andrews_Q <- function(scores) {
  parts <- own_andrews_moments(scores); n <- parts$n
  if (!is.finite(parts$denominator) || parts$denominator <= .Machine$double.eps)
    return(1L)
  a1 <- parts$numerator / parts$denominator
  as.integer(max(1L, min(ceiling(n/2), ceiling(1.1447*(a1*n)^(1/3)))))
}

# Same finite-lag convention as the paper/package: Gamma_l uses n-l.
# No silent centering or prewhitening of the score process.
own_hac <- function(scores, Q = NULL) {
  W <- as.matrix(scores); n <- nrow(W)
  if (any(!is.finite(W))) stop("Invalid HAC score values.")
  if (is.null(Q)) Q <- own_andrews_Q(W)
  if (length(Q) != 1L || !is.finite(Q) || Q < 1 || Q > n || Q != as.integer(Q))
    stop("HAC Q must be an integer between 1 and n.")
  # Form K*W first and use a single W'*(K*W). This equals the lag-by-lag
  # covariance sum but avoids Q separate p-by-p matrix multiplications.
  weighted <- W/n
  if (Q > 1L) for (ell in seq_len(Q-1L)) {
    recent <- seq.int(ell+1L,n); past <- seq_len(n-ell)
    weight <- (1-ell/Q)/(n-ell)
    weighted[recent,] <- weighted[recent,,drop=FALSE]+weight*W[past,,drop=FALSE]
    weighted[past,] <- weighted[past,,drop=FALSE]+weight*W[recent,,drop=FALSE]
  }
  cov <- crossprod(W,weighted)
  list(covariance = (cov+t(cov))/2, Q = Q)
}

own_normal_draws <- function(p, B, seed) {
  had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if (had_seed) old <- get(".Random.seed", envir=.GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old, envir=.GlobalEnv) else
    if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))
      rm(".Random.seed", envir=.GlobalEnv), add=TRUE)
  set.seed(seed)
  matrix(stats::rnorm(p*B), p, B)
}

# Partial out the unpenalized block explicitly. This avoids glmnet's rescaling
# of nonzero penalty.factor entries when some entries are zero.
own_lasso_solver <- function(X, y, unpenalized = integer(), cfg = own_lp_control()) {
  n <- nrow(X); p <- ncol(X)
  U <- sort(unique(as.integer(unpenalized)))
  if (any(!U %in% seq_len(p))) stop("Invalid unpenalized columns.")
  P <- setdiff(seq_len(p), U)
  if (length(U)) {
    qu <- qr(X[,U,drop=FALSE])
    if (qu$rank < length(U)) stop("Rank-deficient unpenalized block.")
    yr <- qr.resid(qu, y)
    Z <- qr.resid(qu, X[,P,drop=FALSE])
  } else { yr <- y; Z <- X[,P,drop=FALSE] }
  function(lambda) {
    if (!is.finite(lambda) || lambda < 0) stop("Invalid lambda.")
    b <- numeric(p)
    if (!length(P)) {
      b[U] <- qr.coef(qu, y)
    } else if (lambda == 0) {
      qx <- qr(X)
      if (qx$rank < p) stop("Zero-penalty regression is not full rank.")
      b <- as.numeric(qr.coef(qx,y))
    } else if (sum(yr^2) > 1e-24) {
      keep <- which(colSums(Z^2) > 1e-20)
      if (length(keep) == 1L) {
        j <- keep; z <- Z[,j]
        zty <- sum(z*yr)/n
        b[P[j]] <- sign(zty)*max(abs(zty)-lambda,0)/(sum(z*z)/n)
      } else if (length(keep) > 1L) {
        fit <- glmnet::glmnet(Z[,keep,drop=FALSE], yr, alpha=1,
          lambda=lambda, standardize=FALSE, intercept=FALSE,
          thresh=cfg$glmnet_tolerance, maxit=cfg$glmnet_maxit)
        if (fit$jerr != 0) stop("Lasso solver did not converge: ", fit$jerr)
        b[P[keep]] <- as.numeric(fit$beta[,1L])
      }
      if (length(U)) b[U] <- qr.coef(qu, y-as.numeric(X[,P,drop=FALSE]%*%b[P]))
    } else if (length(U)) b[U] <- qr.coef(qu,y)
    residual <- y-as.numeric(X%*%b)
    score <- as.numeric(crossprod(X,residual))/n
    active <- P[abs(b[P]) > 1e-9]; inactive <- setdiff(P,active)
    kkt <- max(c(0, abs(score[U]), abs(score[active]-lambda*sign(b[active])),
                 pmax(abs(score[inactive])-lambda,0)))
    if (!is.finite(kkt) || kkt > 1e-4) stop("Lasso KKT check failed: ", kkt)
    list(coefficients=b, residual=residual, lambda=lambda, kkt=kkt,
         unpenalized=U, penalized=P)
  }
}

# Iterative Gaussian plug-in tuning (ASW Section 5.1 / Appendix C.5).
# Common random draws across iterations, not across unrelated regressions.
own_plugin_lasso <- function(X, y, unpenalized=integer(), cfg=own_lp_control(),
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
    lambda <- cfg$pi_constant*cutoff/sqrt(n)
    change <- abs(lambda-old_lambda)/max(old_lambda,.Machine$double.eps)
    fit <- solve_at(lambda)
    trace[[k]] <- data.frame(iteration=k,lambda=lambda,relative_change=change,
      Q=lrv$Q,cutoff=cutoff,negative_eigenvalues=sum(eig$values < -1e-10))
    residual <- fit$residual
    if (change < cfg$pi_tolerance) { converged <- TRUE; break }
    old_lambda <- lambda
  }
  fit$trace <- do.call(rbind,trace); fit$converged <- converged
  fit$seed <- seed; fit
}

# Nodes are fitted once using h=0. Store them in ORIGINAL units so that later
# horizon-specific centering/scaling cannot mix incompatible normalizations.
own_nodewise <- function(X, H, cfg=own_lp_control()) {
  X <- as.matrix(X); p <- ncol(X); n <- nrow(X)
  xm <- colMeans(X); xc <- sweep(X,2L,xm,"-")
  xs <- sqrt(colMeans(xc^2)); z <- sweep(xc,2L,xs,"/")
  if (any(!is.finite(z))) stop("Invalid nodewise design.")
  Gamma <- matrix(0,length(H),p); tau2 <- numeric(length(H)); fits <- list()
  for (i in seq_along(H)) {
    j <- H[i]; other <- setdiff(seq_len(p),j)
    fit <- own_plugin_lasso(z[,other,drop=FALSE],z[,j],integer(),cfg,cfg$seed+i)
    gamma <- fit$coefficients
    denominator <- mean(fit$residual^2)+fit$lambda*sum(abs(gamma))
    if (!is.finite(denominator) || denominator <= 1e-12) stop("Degenerate nodewise residual.")
    Gamma[i,j] <- 1
    Gamma[i,other] <- -gamma*xs[j]/xs[other]
    tau2[i] <- denominator*xs[j]^2
    fit$target <- colnames(X)[j]; fit$target_index <- j
    fit$kkt_normalization_error <- abs(mean(fit$residual*z[,j]) -
      (mean(fit$residual^2)+fit$lambda*sum(abs(gamma))))
    fits[[i]] <- fit
    if (cfg$verbose) cat("nodewise",i,"/",length(H),"lambda",fit$lambda,"\n")
  }
  list(Gamma=Gamma,tau2=tau2,Theta=sweep(Gamma,1L,tau2,"/"),
       H=H,fits=fits,n=n,x_mean=xm,x_scale=xs)
}

own_fit_horizon <- function(X,y,nodes,unpenalized,cfg,seed) {
  s <- own_standardize(X,y); n <- nrow(X)
  fit <- own_plugin_lasso(s$X,s$y,unpenalized,cfg,seed)
  b <- fit$coefficients*s$y_scale/s$x_scale
  xc <- sweep(X,2L,s$x_mean,"-")
  residual <- (y-s$y_mean)-as.numeric(xc%*%b)
  v <- xc%*%t(nodes$Gamma)
  correction <- as.numeric(crossprod(v,residual))/n/nodes$tau2
  estimate <- b[nodes$H]+correction
  # Choose a multivariate bandwidth in standardized score units, as in the
  # reference method. Computing it from raw scores would weight variables
  # differently merely because dollars were changed to thousands of dollars.
  Q <- cfg$hac_Q
  if (is.null(Q)) Q <- own_andrews_Q(
    sweep(v,2L,nodes$x_scale[nodes$H],"/")*(residual/s$y_scale))
  hac <- own_hac(v*residual,Q)
  covariance <- hac$covariance/outer(nodes$tau2,nodes$tau2)/n
  if (any(diag(covariance) < -1e-12)) stop("Negative HAC variance; inspect bandwidth.")
  se <- sqrt(pmax(diag(covariance),0))
  z <- stats::qnorm(1-cfg$alpha/2)
  list(estimate=estimate,se=se,lower=estimate-z*se,upper=estimate+z*se,
       initial=b[nodes$H],correction=correction,coefficients=b,
       intercept=s$y_mean-sum(s$x_mean*b),residual=residual,
       covariance=covariance,HAC_Q=hac$Q,lasso=fit,n=n,
       x_mean=s$x_mean,x_scale=s$x_scale,y_scale=s$y_scale)
}

