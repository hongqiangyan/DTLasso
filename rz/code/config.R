# The single accepted specification. Unemployment is a decimal in selection.
settings <- list(lags=4L,hmax=20L,trim=.10,anchor=.065,c_rho=.75,
 selection_seed=20271002L,node_seed=333L,lp_seed=100333L,
 bootstrap_seed=20261004L,bootstrap_draws=9999L)
selection_cfg <- list(glmnet_tolerance=1e-12,glmnet_maxit=1000000L,
 kkt_tolerance=1e-4,plugin=list(constant=.4,alpha=.05,draws=1000L,
 maxit=60L,tolerance=.01,seed=20261002L,workers=4L),
 residual_tolerance=.01,stable_iterations=2L)
lp_settings <- list(pi_constant=.4,pi_alpha=.05,pi_draws=1000L,
 pi_maxit=60L,pi_tolerance=.01,glmnet_tolerance=1e-12,
 normalization="paper",verbose=FALSE)
