# Delta inference from the joint HAC covariance of cumulative responses.
# Order: (GDP_low, GDP_high, spending_low, spending_high).
mi_inference <- function(beta, covariance, alpha = .05) {
  stopifnot(length(beta) == 4L, identical(dim(covariance), c(4L,4L)),
    all(is.finite(c(beta,covariance))), all(abs(beta[3:4]) > 1e-12),
    alpha > 0, alpha < 1)
  scale <- max(abs(covariance), 1e-20)
  stopifnot(min(eigen(covariance,symmetric=TRUE,only.values=TRUE)$values) >= -1e-9*scale)
  ratio <- beta[1:2]/beta[3:4]
  J <- rbind(c(1/beta[3],0,-beta[1]/beta[3]^2,0),
             c(0,1/beta[4],0,-beta[2]/beta[4]^2))
  ratio_covariance <- J %*% covariance %*% t(J)
  gradient <- J[2,]-J[1,]
  variance <- as.numeric(t(gradient)%*%covariance%*%gradient)
  stopifnot(all(diag(ratio_covariance)>=-1e-12),variance>=-1e-12)
  se <- sqrt(pmax(diag(ratio_covariance),0)); difference_se <- sqrt(max(variance,0))
  difference <- ratio[2]-ratio[1]; z <- qnorm(1-alpha/2)
  list(ratio=ratio,ratio_se=se,ratio_covariance=ratio_covariance,jacobian=J,
    ratio_lower=ratio-z*se,ratio_upper=ratio+z*se,
    ratio_p=2*pnorm(-abs(ratio/se)),difference=difference,
    difference_se=difference_se,difference_gradient=gradient,
    difference_lower=difference-z*difference_se,
    difference_upper=difference+z*difference_se,
    difference_p=2*pnorm(-abs(difference/difference_se)))
}
