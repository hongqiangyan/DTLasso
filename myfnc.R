
library(Matrix)

library(glmnet)

library(expm)

library(flare)

# Nodewise estimation of the covariance matrix
safe_inverse <- function(X, cond_threshold = 1e8) {
  XtX <- t(X) %*% X
  if (kappa(XtX) < cond_threshold) {
    solve(XtX)
  } else {
    cat("solve() failed (ill-conditioned), fallback to nodewise...\n")
    estimate_nodewise_covariance(X)$Theta
  }
}
estimate_nodewise_covariance <- function(Y, ic = 'GIC') {
  # initialization
  p <- ncol(Y)
  n <- nrow(Y)
  #Y <- Y- t((apply(Y,2,mean))%*%matrix(1,1,n)) # Y is de-meaned
  
  C <- matrix(0, p, p)
  diag(C) <- 1
  tau <- NULL
  
  # Loop over the assets
  for (j in 1:p) {
    # Estimate the Lasso
    #print(colSums(is.na(Y))  )
    jlas <- glmnet(
      x = Y[, -j],
      y = Y[, j],
      family = 'gaussian',
      intercept = FALSE,
      standardize = TRUE
    )    # Get fit
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
      jbic  <- log(jsig) + 2 * jlas$df # AIC
    # Index of selected model
    jind  <- which.min(jbic)
    # Get the parameters
    jpar <- jlas$beta[, jind]
    # Computing tau squared
    jtau <- sum(jres[, jind]^2) / n + (1 / 2) * jlas$lambda[jind] * sum(abs(jpar)) # using (10)
    # Storing the parameters
    C[j, -j] <- -jpar
    tau <- c(tau, jtau)
  }
  
  # Construct T-squared inverse
  T2inv <- diag(1 / tau)
  
  # Construct Theta-hat
  Theta <- T2inv %*% C
  
  # sparsity
  sp <- sum(Theta == 0) / (p^2)
  
  return(list(
    Theta = Theta / n,
    ####inverse of t(Y) %*% Y
    sparsity = sp
  ))
}



Threshold_Lasso <- function(X1,
                            y,
                            taugrid = taugrid,
                            lambda = NULL,
                            intercept1 = TRUE,
                            intercept2 = TRUE) {
  p1 <- p2 <- ncol(X1)
  
  n <- nrow(X1)
  
  best_objective <- Inf
  best_tau <- NULL
  best_coefficients <- NULL
  # Iterate over each tau in the grid
  for (cnt in 1:length(taugrid)) {
    tau <- taugrid[cnt]
    X11 <- X1
    p1 <- p2 <- ncol(X1)
    if (intercept2 == TRUE) {
      int2 <- rep(0, n)
      int2[thdvar < tau] <- 1
      X11 <- cbind(int2, X11)
      p1 <- p1 + 1
    }
    if (intercept1 == TRUE) {
      X11 <- cbind(1, X11)
      p1 <- p1 + 1
    }
    p <- p1 + p2
    col.norm1 <- 1 / sqrt((1 / n) * diag(t(X11) %*% X11))
    X11 <- X11 %*% diag(col.norm1)
    X2 <- X1
    X2[thdvar >= tau, ] <- 0
    col.norm2 <- 1 / sqrt((1 / n) * diag(t(X2) %*% X2))
    
    X2 <- X2 %*% diag(col.norm2)
    
    
    X <- cbind(X11, X2)
    
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
      penalty <- sum(abs(coefficients)^2)
      # Calculate the Objective Function
      Objective <- sqrt(rss / n) + lambda_val * penalty
      lambda <- NULL
    } else if (lambda %in% c("bic", "gic", "wic", "aic", "mic")) {
      # ---- IC-based lambda selection ----
      lambda_seq <- exp(seq(log(1e-4), log(1), length.out = 100))
      ic_vec <- numeric(length(lambda_seq))
      
      for (i in seq_along(lambda_seq)) {
        fit <- glmnet::glmnet(
          X,
          y,
          lambda = lambda_seq[i],
          standardize = FALSE,
          intercept = FALSE
        )
        y_hat <- predict(fit, newx = X)
        rss <- sum((y - y_hat)^2)
        df <- sum(coef(fit) != 0)
        jsig <- rss / n
        
        # Compute IC value based on lambda type
        ic_vec[i] <- switch(
          lambda,
          bic = log(jsig) + log(n) * df / n,
          gic = log(jsig) + log(p) * df / n * log(log(n)),
          wic = log(jsig) + log(n) * df / n * log(log(p)),
          # BIC (Wang,2010)
          mic = jsig + df * log(n) * log(log(p)) / n,
          aic = log(jsig) + 2 * df / n,
          stop("Unsupported IC type")
        )
      }
      best_index <- which.min(ic_vec)
      lambda_val <- lambda_seq[best_index]
      final_fit <- glmnet::glmnet(
        X,
        y,
        lambda = lambda_val,
        standardize = FALSE,
        intercept =  FALSE
      )
      coefficients <- coef(final_fit)[-1]
      penalty <- sum(abs(coefficients)^2)
      Objective <-  sqrt(rss / n) + lambda_val * penalty
    } else if (identical(lambda, "cv")) {
      # ---- CV-based lambda selection ----
      outLas <-  glmnet::cv.glmnet(
        X,
        y,
        nfolds = nrow(X),
        standardize = TRUE,
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
      Objective <- sqrt(rss / n) + outLas$lambda.min * penalty
      lambda_val <-   outLas$lambda.min
      #if not given
    } else {
      stop("lambda must be one of: NULL, 'bic', or 'cv'")
    }
    
    # Check if this Objective is the smallest we've encountered
    if (Objective < best_objective) {
      best_objective <- Objective
      best_tau <- tau
      best_coefficients <- coefficients
      X_tau <- X2
      col.norm_tau <- col.norm2
    }
  }
  X_tau <-  X_tau %*% diag(1 / col.norm_tau)
  
  col.norm <- c(col.norm1, col.norm_tau)
  X11 <-  X11 %*% diag(1 / col.norm1)
  
  #
  # cat("Dimensions of best_coefficients:\n")
  # print(dim(best_coefficients))
  # cat("Length of best_coefficients (if vector): ", length(best_coefficients), "\n")
  #
  # cat("Dimensions of diag(col.norm):\n")
  # print(dim(diag(col.norm)))
  # cat("Length of col.norm: ", length(col.norm), "\n")
  best_coefficients <- best_coefficients %*% diag(col.norm)
  
  # Output depending on whether intercept is requested
  return(
    list(
      coefficients = best_coefficients,
      tau = best_tau,
      objective = best_objective,
      X1 = X11,
      X_tau = X_tau,
      lambda_val = lambda_val
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
                                      resol = 1.3,
                                      maxiter = 50,
                                      threshold = 1e-2,
                                      verbose = TRUE) {
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
  
  col.norm1 <- 1 / sqrt((1 / n) * diag(t(X1) %*% X1))
  lasso_result <- Threshold_Lasso(
    X1,
    y,
    taugrid = taugrid,
    lambda = lambda,
    intercept1 = intercept1,
    intercept2 = intercept2
  )
  htheta <- t(lasso_result[[1]])
  
  best_tau <- lasso_result[[2]]
  
  X1 <- lasso_result[[4]]
  X2 <- lasso_result[[5]]
  col.norm1 <- 1 / sqrt((1 / n) * diag(t(X1) %*% X1))
  
  col.norm2 <- 1 / sqrt((1 / n) * diag(t(X2) %*% X2))
  
  
  col.norm <- c(col.norm1, col.norm2)
  if (intercept2 == TRUE) {
    p1 <- p1 + 1
  }
  if (intercept1 == TRUE) {
    p1 <- p1 + 1
  }
  pp <- p1 + p2
  
  Xb <- cbind(X1, X2)
  sigma.hat <- (1 / n) * (t(Xb) %*% Xb)
  invgram<-safe_inverse(Xb)
  if (is.list(invgram)) {
    M <- invgram[[1]] * n
    nz_node<-invgram[[2]]
    ifnode=TRUE
  } else {
    M <- invgram * n
    ifnode=FALSE
  }
  
  unbiased.Lasso <- as.numeric(htheta + (M %*% t(Xb) %*% (y - Xb %*% (htheta))) /
                                 n)
  
  # noise <- NoiseSd( unbiased.Lasso,M, n);
  # s.hat <- noise$sd;
  
  resids <- y - Xb %*% (htheta)
  u2 <- as.numeric(resids)^2
  meat <- crossprod(Xb * sqrt(u2)) / n
  varmx <-   M %*% meat %*% t(M) / n
  stder <- sqrt(diag(varmx))
  
  interval.sizes <- qnorm(1 - (alpha / 2)) * stder
  lambda <- lasso_result$lambda_val
  addlength <- rep(0, pp)
  if (ifnode == TRUE) {
    MM <- M %*% sigma.hat - diag(pp)
    for (i in 1:pp) {
      effectivemuvec <- sort(abs(MM[i, ]), decreasing = TRUE)
      effectivemuvec <- effectivemuvec[1:(nz_node - 1)]
      addlength[i] <- sqrt(sum(effectivemuvec^2)) * lambda
    }
  }
  
  # nz_hat <- sum(htheta != 0)
  # MM <- M%*%sigma.hat - diag(pp);
  # for (i in 1:pp){
  #   effectivemuvec <- sort(abs(MM[i,]),decreasing=TRUE);
  #   effectivemuvec <- effectivemuvec[1:nz_hat];
  #   addlength[i] <- sqrt(sum(effectivemuvec*effectivemuvec))*lambda;
  # }
  
  p.vals <- 2 * (1 - pnorm(abs(unbiased.Lasso) / stder))
  
  
  returnList <- list(
    "noise.sd" = stder,
    #"norm0" =  nz_hat,
    "coef" = htheta,
    "unb.coef" = unbiased.Lasso,
    "low.lim" = unbiased.Lasso - interval.sizes - addlength,
    "up.lim" = unbiased.Lasso + interval.sizes + addlength,
    "pvals" = p.vals,
    "best_tau" = best_tau
  )
  return(returnList)
}



Linear_Lasso <- function(X1,
                         y,
                         lambda = NULL,
                         intercept1 = TRUE,
                         col.norm1 = col.norm1) {
  p1 <- ncol(X1)
  
  n <- nrow(X1)
  
  best_objective <- Inf
  best_tau <- NULL
  best_coefficients <- NULL
  # Iterate over each tau in the grid
  
  X11 <- X1
  p1 <- ncol(X1)
  if (intercept1 == TRUE) {
    X11 <- cbind(1, X11)
    p1 <- p1 + 1
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
    penalty <- sum(abs(coefficients)^2)
    # Calculate the Objective Function
    Objective <- sqrt(rss / n) + lambda_val * penalty
    lambda <- NULL
  } else if (lambda %in% c("bic", "gic", "wic", "aic", "mic")) {
    # ---- IC-based lambda selection ----
    lambda_seq <- exp(seq(log(1e-4), log(1), length.out = 100))
    ic_vec <- numeric(length(lambda_seq))
    
    for (i in seq_along(lambda_seq)) {
      fit <- glmnet::glmnet(
        X,
        y,
        lambda = lambda_seq[i],
        standardize = FALSE,
        intercept = FALSE
      )
      y_hat <- predict(fit, newx = X)
      rss <- sum((y - y_hat)^2)
      df <- sum(coef(fit) != 0)
      jsig <- rss / n
      
      # Compute IC value based on lambda type
      ic_vec[i] <- switch(
        lambda,
        bic = log(jsig) + log(n) * df / n,
        gic = log(jsig) + log(p) * df / n * log(log(n)),
        wic = log(jsig) + log(n) * df / n * log(log(p)),
        # BIC (Wang,2010)
        mic = jsig + df * log(n) * log(log(p)) / n,
        aic = log(jsig) + 2 * df / n,
        stop("Unsupported IC type")
      )
    }
    best_index <- which.min(ic_vec)
    lambda_val <- lambda_seq[best_index]
    final_fit <- glmnet::glmnet(
      X,
      y,
      lambda = lambda_val,
      standardize = FALSE,
      intercept =  FALSE
    )
    coefficients <- coef(final_fit)[-1]
    penalty <- sum(abs(coefficients)^2)
    Objective <-  sqrt(rss / n) + lambda_val * penalty
  } else if (identical(lambda, "cv")) {
    # ---- CV-based lambda selection ----
    outLas <-  glmnet::cv.glmnet(
      X,
      y,
      nfolds = 10,
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
    Objective <- sqrt(rss / n) + outLas$lambda.min * penalty
    lambda_val <-   outLas$lambda.min
    #if not given
  } else {
    stop("lambda must be one of: NULL, 'bic', or 'cv'")
  }
  
  
  best_objective <- Objective
  best_coefficients <- coefficients
  X11 <-  X11 %*% diag(1 / col.norm1)
  
  #
  # cat("Dimensions of best_coefficients:\n")
  # print(dim(best_coefficients))
  # cat("Length of best_coefficients (if vector): ", length(best_coefficients), "\n")
  #
  # cat("Dimensions of diag(col.norm):\n")
  # print(dim(diag(col.norm)))
  # cat("Length of col.norm: ", length(col.norm), "\n")
  best_coefficients <- best_coefficients %*% diag(col.norm1)
  
  # Output depending on whether intercept is requested
  return(
    list(
      coefficients = best_coefficients,
      objective = best_objective,
      X1 = X11,
      lambda_val = lambda_val
    )
  )
}



Linear_debiased_Lasso <- function (X1,
                                   y,
                                   alpha = 0.05,
                                   lambda = NULL,
                                   mu = NULL,
                                   intercept1 = TRUE,
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
    col.norm1 = col.norm1
  )
  htheta <- t(lasso_result[[1]])
  X1 <- lasso_result[[3]]
  col.norm1 <- 1 / sqrt((1 / n) * diag(t(X1) %*% X1))
  if (intercept1 == TRUE) {
    p1 <- p1 + 1
  }
  pp <- p1
  
  Xb <- X1
  sigma.hat <- (1 / n) * (t(Xb) %*% Xb)
  
  M <- (safe_inverse(Xb)) * n
  
  unbiased.Lasso <- as.numeric(htheta + (M %*% t(Xb) %*% (y - Xb %*% (htheta))) /
                                 n)
  # noise <- NoiseSd( unbiased.Lasso,M, n);
  # s.hat <- noise$sd;
  
  resids <- y - Xb %*% (htheta)
  u2 <- as.numeric(resids)^2
  meat <- crossprod(Xb * sqrt(u2)) / n
  varmx <-   M %*% meat %*% t(M) / n
  stder <- sqrt(diag(varmx))
  
  interval.sizes <- qnorm(1 - (alpha / 2)) * stder
  lambda <- lasso_result$lambda_val
  
  addlength <- rep(0, pp)
  # nz_hat <- sum(htheta != 0)
  # MM <- M%*%sigma.hat - diag(pp);
  # for (i in 1:pp){
  #   effectivemuvec <- sort(abs(MM[i,]),decreasing=TRUE);
  #   effectivemuvec <- effectivemuvec[1:nz_hat];
  #   addlength[i] <- sqrt(sum(effectivemuvec*effectivemuvec))*lambda;
  # }
  
  p.vals <- 2 * (1 - pnorm(abs(unbiased.Lasso) / stder))
  
  
  returnList <- list(
    "noise.sd" = stder,
    #"norm0" =  nz_hat,
    "coef" = htheta,
    "unb.coef" = unbiased.Lasso,
    "low.lim" = unbiased.Lasso - interval.sizes - addlength,
    "up.lim" = unbiased.Lasso + interval.sizes + addlength,
    "pvals" = p.vals
  )
  return(returnList)
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
