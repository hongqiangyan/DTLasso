library(MASS)
set.seed(123)
source("myfnc.R")
B <- 100
results <- data.frame(
  Th_FWER = numeric(B),
  Th_Power = numeric(B),
  Lin_FWER = numeric(B),
  Lin_Power = numeric(B)
)
for (k in 1:B) {
  # --- Parameters ---
  p <- 300; n <- 400; s0 <- 15
  b <- 2; b1 <- 1; b0 <- 0
  tau_0 <- 0.5; rho_1 <- 0.5
  taugrid <- seq(0.15, 0.85, 0.01)
  
  # --- Data generation ---
  cov_XG <- rho_1 ^ abs(row(matrix(1:p,p,p)) - col(matrix(1:p,p,p)))
  X1 <- MASS::mvrnorm(n, rep(0,p), cov_XG)
  thdvar <- runif(n)
  
  X2 <- X1; X2[thdvar >= tau_0,] <- 0
  X <- cbind(1, X1, 1, X2)
  theta0 <- c(b0, rep(b,s0), rep(0,p-s0),
              b0, rep(0,s0), rep(b1,s0), rep(0,p-2*s0))
  y <- X %*% theta0 + rt(n, df=10)
  
  # --- Threshold Model ---
  th <- Threshold_debiased_Lasso(X1, y, taugrid=taugrid, intercept1=FALSE,
                                 lambda='gic', intercept2=FALSE, verbose=FALSE)
  
  # --- Linear Model ---
  lin <- Linear_debiased_Lasso(X1, y, intercept1=FALSE, lambda='gic', verbose=FALSE)
  
  # --- Index sets ---
  active_index <- c(1:s0, (p+s0+1):(p+2*s0))
  null_index   <- c((s0+1):p, (p+1):(p+s0), (p+2*s0+1):(2*p))
  
  # --- Holm adjustment ---
  holm_th <- p.adjust(th$pvals, method="holm")
  holm_lin <- p.adjust(lin$pvals, method="holm")
  
  active_index_lin <- 1:s0
  # Null set in Linear: rest of variables
  null_index_lin <- (s0+1):p
  
  # --- Record Threshold results ---
  results$Th_FWER[k] <- any(holm_th[null_index] <= 0.05)
  results$Th_Power[k] <- mean(holm_th[active_index] <= 0.05)
  
  # --- Record Linear results ---
  results$Lin_FWER[k] <- any(holm_lin[null_index_lin] <= 0.05)
  results$Lin_Power[k] <- mean(holm_lin[active_index_lin] <= 0.05)
}

# --- Summary ---
summary_results <- data.frame(
  Model = c("Threshold", "Linear"),
  Mean_FWER = c(mean(results$Th_FWER), mean(results$Lin_FWER)),
  Mean_Power = c(mean(results$Th_Power), mean(results$Lin_Power))
)
print(summary_results)
#--with a threshold effect
coef_path <- th$coef
unb.coef_path <- th$unb.coef

dev.new()
plot(coef_path,ylim=c(-3*b, 3*b), main='Confidence Intervals for debiased LASSO', ylab='', xlab = 'Coefficients');
points(unb.coef_path,col="blue");
theta0 <- c(rep(b,s0), rep(0,p-s0), rep(0,s0), rep(b1,s0), rep(0,p-2*s0));
points(theta0,col="green");
lines(th$up.lim,col="red");
lines(th$low.lim,col="red");
legend('topright', legend=c('LASSO','de-biased LASSO','Ground-truth','Confidence Intervals'), col=c('black', 'blue','green','red'), pch=c(1,1,1,NA_integer_), lty = c(0,0,0,1))

#--without a threshold effect
coef_path_lin<- lin$coef
unb.coef_path_lin <- lin$unb.coef

dev.new()
plot(coef_path_lin,ylim=c(-3*b, 3*b), main='Confidence Intervals for debiased LASSO', ylab='', xlab = 'Coefficients');
points(unb.coef_path_lin,col="blue");
theta0 <- c( rep(b,s0), rep(0,p-s0),rep(0,s0), rep(b1,s0), rep(0,p-2*s0));
points(theta0,col="green");
lines(th$up.lim,col="red");
lines(th$low.lim,col="red");
legend('topright', legend=c('LASSO','de-biased LASSO','Ground-truth','Confidence Intervals'), col=c('black', 'blue','green','red'), pch=c(1,1,1,NA_integer_), lty = c(0,0,0,1))


