rm(list=ls())
start.time <- proc.time()[3]
cat(date(), file='growth_app_lr_out.asc', sep='\n', append=F)
source("myfnc.R")

library('lars')
###########################################################
# Data read and construction
###########################################################
g_data <- read.table('GrowthRegressions_data_lr.cvs')

dep <-g_data[,'gr']
reg <-as.matrix(g_data[,c(4:ncol(g_data))])
reg <- cbind(1,reg)
thres <- g_data[,'lr']

nx <- ncol(reg)
nobs <- nrow(reg)
sorted.thres <- sort(thres)
trim.factor <- round(length(sorted.thres)*.10,0)
grid.tau <- sorted.thres[trim.factor:(nobs-trim.factor)]

colnames(reg) <- c('const',colnames(reg)[2:ncol(reg)])
reg.names <- colnames(reg)

###########################################################
# Threshold LASSO function
###########################################################
tlasso <- function (x, y, q, s, grid.tau) {
  
  ngrid <- length(grid.tau)
  nobs <- length(y)
  nreg <- ncol(x)*2
  delta.hat.grid <- matrix(rep(NA), nrow=ngrid, ncol=nreg)
  obj.v <- matrix(rep(NA), nrow=ngrid, ncol=1)
  norm.x <- matrix(rep(NA), nrow=1, ncol=nreg)
  delta.hat <- matrix(rep(NA), nrow=1, ncol=nreg)
  tau.hat <- NA
  for(i in 1:ngrid) {
    ind <- ( q < grid.tau[i] )
    x.reg <- cbind(x,x*ind)
    m <- lars(x.reg,y)
    delta.hat.grid[i,] <- coef(m, s=s, mode="lambda")
    yhat <- predict(m, s=s, mode="lambda",x.reg)$fit
    uhat <- y - yhat
    ssr <- t(uhat)%*%uhat / nobs
    for(j in 1:nreg){
      norm.x[1,j] <- sqrt( t((x.reg[,j]))%*% (x.reg[,j]) / (nobs) )
    }
    p <- norm.x %*% abs(delta.hat.grid[i,])
    obj.v[i,1] <- ssr + s*p    	
  }
  opt<-which.min(obj.v)
  delta.hat <- delta.hat.grid[opt,]
  tau.hat <- grid.tau[opt]
  
  
  
  rval <- list(model=m, est=c(s, tau.hat, delta.hat))
  
  return(rval)
  
}




###########################################################
# Calculate regularization parameter lambda
###########################################################
t_0 <- grid.tau[1]
ind_t0<-matrix(rep(thres < t_0),nobs,nx)
x_tau <- reg * ind_t0
r_n <- min( diag( t(x_tau) %*% x_tau ) / diag(t(reg)%*%reg) )
A <- 2*sqrt(2)	
sigma <- sd(dep)	
lambda_max <- A*sigma*sqrt(log(nx*3)/(nobs*r_n))




####################################################
# Final model with lambda from cross validation
####################################################



lambda <- 0.0044
f.result <- tlasso(x=reg, y=dep, q=thres, s=lambda, grid.tau)
f.model <- f.result$model
f.coef<-f.result$est

initial_lasso <- f.coef[c(3, 50, 4:49, 51:96)]
n<-nrow(reg)
tauest<-f.coef[2]
f.coef[1]
x11 <- reg[,-1]
int1 <- rep(1, n)
int2 <- rep(0, n)
int2[thres< tauest] <- 1
x22<-x11
x22[thres>=tauest] <-0
xx<-cbind(int1,int2,x11,x22)

sigma.hat <- (1 / n) * (t(xx) %*%  xx)
invgram<-safe_inverse(sigma.hat)
if (is.list(invgram)) {
  M <- invgram[[1]] * n
  nz_node<-invgram[[2]]
  ifnode=TRUE
} else {
  M <- invgram * n
  ifnode=FALSE
}
resids <- dep - xx %*% (initial_lasso)
u2 <- as.numeric(resids)^2
meat <- crossprod(xx* sqrt(u2)) / n
varmx <-   M %*% meat %*% t(M) / n
stder <- sqrt(diag(varmx))
unbiased.Lasso <- as.numeric(initial_lasso + (M %*% t(xx) %*% resids) /n)
p.vals <- 2 * (1 - pnorm(abs(unbiased.Lasso) / stder))


base_names <- reg.names[-1]  
threshold_names <- paste0("threshold_", base_names)

colnames_new <- c("intercept",                    
                  "threshold_intercept",           
                  base_names,                     
                  threshold_names)                     
p.vals <- matrix(p.vals, nrow = 1)
colnames(p.vals) <- colnames_new

unbiased.Lasso <- matrix(unbiased.Lasso, nrow= 1) 
colnames(unbiased.Lasso) <- colnames_new

stder <- matrix(stder, nrow= 1) 
colnames(stder) <- colnames_new

significant_indices_01 <- which(p.vals < 0.1)
significant_coefs01 <- unbiased.Lasso[significant_indices_01]

# name
significant_vars01 <- colnames_new[significant_indices_01]

p.vals01 <- p.vals[significant_indices_01]

significant_sd01 <- stder[significant_indices_01]
# combine results
significant_results01 <- data.frame(
  Variable    = significant_vars01,
  p.vals = p.vals01,
  Coefficient = significant_coefs01,
  sd          = significant_sd01,
  row.names   = NULL
)







