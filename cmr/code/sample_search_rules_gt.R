# Sample and candidate checks only: no outcome-based threshold screening.
q60_city_pool <- function(raw, dec, variables, q_name) {
  take <- raw$decade == dec & raw$pre_city_size_flag == 'True'
  pool <- raw[take, variables, drop=FALSE]
  pool <- pool[complete.cases(pool), , drop=FALSE]
  pool[pool[[q_name]] >= 0 & pool[[q_name]] <= .60, , drop=FALSE]
}

q60_candidate_grid <- function(X, q, min_regime_share = .15) {
  stopifnot(nrow(X)==length(q),all(is.finite(X)),all(is.finite(q)),
            all(q>=0 & q<=.60),ncol(X)>0,
            length(min_regime_share)==1L,is.finite(min_regime_share),
            min_regime_share>=0,min_regime_share<.5)
  # Count actual observations, including ties, on the unchanged raw q.
  # This restricts candidate thresholds; no tract is removed from regression.
  min_regime_obs <- ceiling(min_regime_share * length(q))
  grid <- sort(unique(c(0,round(q[q<=.50],3),.50)))
  check_one <- function(tau) {
    above <- q > tau
    nr <- sum(above); nl <- length(q)-nr  # left includes Q == tau
    W <- cbind(as.numeric(above),X,X*as.numeric(above))
    constant <- apply(W,2,function(v)all(v==v[1]))
    reason <- if(!nl || !nr)'empty_regime' else
      if(min(nl,nr)<min_regime_obs)'below_min_regime_share' else
      if(any(constant))'constant_design_column' else 'admissible'
    data.frame(tau=tau,n_left=nl,n_right=nr,
               min_regime_obs=min_regime_obs,
               constant_columns=sum(constant),reason=reason)
  }
  audit <- do.call(rbind,lapply(grid,check_one))
  valid <- audit$tau[audit$reason=='admissible']
  if(length(valid)<2L)stop('Fewer than two admissible candidates in [0,0.50] after regime-size and design checks.')
  # Data candidates retain three decimals; their midpoint needs at most four.
  # Canonicalize thresholds before deduplication; never round the observed q.
  lo <- min(valid); hi <- max(valid); anchor <- round((lo+hi)/2,4)
  anchor_check <- check_one(anchor)
  if(anchor_check$reason!='admissible')stop('Feasible-interval midpoint fails regime-size or design checks.')
  list(grid=sort(unique(round(c(valid,anchor),4))),tau_star=anchor,tau_lo=lo,tau_hi=hi,
       min_regime_share=min_regime_share,min_regime_obs=min_regime_obs,
       audit=audit,anchor_check=anchor_check,
       excluded=sum(audit$reason!='admissible'))
}
