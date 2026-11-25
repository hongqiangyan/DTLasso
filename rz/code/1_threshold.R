# Consistent trimming at EVERY threshold search: plug-in scores, anchor fit, test.
# Run: Rscript code/1_threshold.R [--verify]. Only the accepted 10--90 range.
# The diagnostic bootstrap does not by itself validate the debiased approximation.
entry <- grep("^--file=",commandArgs(FALSE),value=TRUE)
ext <- dirname(normalizePath(sub("^--file=","",entry)))
source(file.path(ext,"common.R")); ctx <- extension_setup(ext)
source(file.path(ctx$root,"threshold_plugin_anchor.R")); tn_install_recovery()
args <- commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% "--verify"))
trim <- settings$trim; verify <- "--verify" %in% args
label <- "accepted_10_90"
out <- file.path(dirname(ext),"results/threshold")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
logfile <- file.path(out,"run.log")
seed <- settings$selection_seed
# One immutable grid shared by calibration, estimation, and the threshold test.
percent <- ctx$sample
bounds <- as.numeric(quantile(percent$q,c(trim,1-trim),type=7))
observed <- sort(unique(percent$q))
taus <- observed[observed>=bounds[1] & observed<=bounds[2]]
if (percent$tau_star>=bounds[1] && percent$tau_star<=bounds[2])
 taus <- sort(unique(c(taus,percent$tau_star)))
working_data <- ctx$bundle$data; working_data$unemp <- working_data$unemp/100
s <- ec_search_sample(working_data,ctx$bundle$controls,4L,"y",.065)
ec_verify_alignment(s,working_data,ctx$bundle$controls,4L,"y")
anchor_tau_grid <- function(q,tau_star) {
 stopifnot(identical(q,s$q),identical(tau_star,s$tau_star))
 taus/100
}
stopifnot(identical(s$X1,percent$X1),identical(s$y,percent$y),
 all(vapply(taus,function(t) identical(percent$q<t,s$q<t/100),logical(1))))
sources <- c(ctx$sources,file.path(ctx$root,"threshold_plugin_anchor.R"),
 file.path(ext,"1_threshold.R"))
source_hash <- unname(tools::md5sum(sources))
common_key <- list(data_md5=ctx$data_md5,source_md5=source_hash,trim=trim,
 bounds_percent=bounds,taus_percent=taus,anchor=.065,c_rho=.75,
 selection_cfg=selection_cfg,selection_seed=seed)
selections <- summaries <- list()
for (shock in "news") {
 sample <- s
 key <- c(common_key,list(shock=shock,outcome="y",horizon=0L))
 directory <- file.path(out,shock); dir.create(directory,showWarnings=FALSE)
 path <- file.path(directory,"selection.rds")
 if(file.exists(path)) {
  a <- readRDS(path); stopifnot(identical(a$key,key),identical(a$sample,sample))
  if(verify) {
   engine <- tap_engine(sample,selection_cfg,seed)
   a$verification <- tap_verify(sample,a$result,engine,selection_cfg)
  }
 } else {
  if(verify) stop("Missing selection: ",shock)
  anchor_log(label,shock,"recalibrating lambda on",length(taus),"candidates",file=logfile)
  engine <- tap_engine(sample,selection_cfg,seed)
  result <- tap_calibrate(engine,.75,selection_cfg,logfile)
  check <- tap_verify(sample,result,engine,selection_cfg)
  a <- list(key=key,sample=sample,result=result,verification=check)
  te_atomic_save(a,path)
 }
 selections[[shock]] <- a
 z <- a$result$selected; high <- sample$q>=z$tau
 final <- a$result$converged
 summaries[[shock]] <- data.frame(range=label,shock=shock,trim=trim,
  lower_quantile=bounds[1],upper_quantile=bounds[2],candidate_min=min(taus),
  candidate_max=max(taus),candidates=length(taus),c_rho=.75,status=a$result$status,
  iterations=nrow(a$result$trace),lambda=if(final) z$lambda else NA_real_,
  tau_percent=if(final) 100*z$tau else NA_real_,
  last_lambda=z$lambda,last_tau_percent=100*z$tau,
  at_search_boundary=if(final) z$index %in% c(1L,length(taus)) else NA,
  n_high=if(final) sum(high) else NA_integer_,n_low=if(final) sum(!high) else NA_integer_,
  nonzero_shock_high=if(final) sum(high & sample$X1[,1]!=0) else NA_integer_,
  nonzero_shock_low=if(final) sum(!high & sample$X1[,1]!=0) else NA_integer_,
  plugin_bandwidth=a$result$calibration$bandwidth,
  max_KKT=a$verification$max_KKT)
 write.csv(a$result$trace,file.path(directory,"iteration_trace.csv"),row.names=FALSE)
 write.csv(z$profile,file.path(directory,"objective_profile.csv"),row.names=FALSE)
 capture.output(a$verification,file=file.path(directory,"selection_verification.txt"))
 if(exists("engine", inherits=FALSE)) rm(engine)
 gc(FALSE)
}
selection_table <- do.call(rbind,summaries)
write.csv(selection_table,file.path(out,"selection_summary.csv"),row.names=FALSE)
print(selection_table,row.names=FALSE)
if(!selections$news$result$converged) stop("News selection did not converge; no bootstrap at an unconverged lambda")

# GDP/news equation test. Refit every candidate and the null at this range's
# lambda. Saved fits are tied to exact data, settings, seeds, and source hashes.
cfg <- do.call(own_lp_control,c(lp_settings,list(seed=settings$node_seed)))
seed_table <- read.csv(file.path(ext,"bootstrap_seeds.csv"))
stopifnot(identical(seed_table$tau_percent,taus))
seed_ids <- as.integer(seed_table$seed_id)
lambda <- selections$news$result$selected$lambda
null <- te_null_fit(percent,lambda,cfg)
key <- c(common_key,list(lambda=lambda,nodewise=cfg,node_seed_ids=seed_ids,
 B=9999L,bootstrap_seed=20261004L))
directory <- file.path(out,"news"); path <- file.path(directory,"bootstrap.rds")
if(file.exists(path)) {
 saved <- readRDS(path); stopifnot(identical(saved$key,key))
 fits <- saved$fits; result <- saved$result
} else {
 if(verify) stop("Missing bootstrap")
 fit_one <- function(k) {
  local <- cfg;local$seed <- cfg$seed+10000L*seed_ids[k]
  f <- te_fit_tau(percent,taus[k],lambda,null$residual,local)
  f$node_seed <- local$seed; f
 }
 fits <- parallel::mclapply(seq_along(taus),fit_one,
   mc.cores=selection_cfg$plugin$workers,mc.preschedule=TRUE,mc.set.seed=FALSE)
 if(any(vapply(fits,inherits,logical(1),"try-error"))) stop("Precision estimation failed")
 M <- do.call(cbind,lapply(fits,`[[`,"scores"))
 estimate <- unlist(lapply(fits,`[[`,"estimate"),use.names=FALSE)
 ell <- te_bandwidth(sweep(M,2L,colMeans(M),"-"))
 anchor_log(label,"bootstrap B",key$B,"block length",ell,file=logfile)
 result <- te_block_bootstrap(M,estimate,ell,B=key$B,seed=key$bootstrap_seed)
 saved <- list(key=key,sample=percent,selection=selections$news$result$selected,
  null=null,fits=fits,result=result)
 te_atomic_save(saved,path)
}
# Independent reconstruction and multiplier algebra, plus inverse-row diagnostics.
n <- length(percent$y)
null_error <- max(abs(percent$y-null$intercept-as.numeric(percent$X1 %*% null$coefficients)-null$residual))
errors <- numeric(length(fits)); diagnostics <- vector("list",length(fits))
for(k in seq_along(fits)) {
 f <- fits[[k]]; d <- anchor_design(percent$X1,percent$y,percent$q,f$tau)
 xc <- sweep(d$raw,2L,d$means,"-"); V <- xc %*% t(f$Theta)
 errors[k] <- max(abs(c(percent$y-f$intercept-as.numeric(d$raw %*% f$coefficients)-f$residual,
  f$estimate-f$coefficients[f$H]-as.numeric(crossprod(V,f$residual))/n,
  f$scores-V*null$residual)))
 theta_standard <- sweep(sweep(f$Theta,1L,d$weights[f$H],"*"),2L,d$weights,"*")
 error <- theta_standard %*% (crossprod(d$X)/n)-diag(ncol(d$X))[f$H,,drop=FALSE]
 high <- percent$q>=f$tau
 diagnostics[[k]] <- data.frame(tau=f$tau,coefficient=f$names,estimate=f$estimate,
  initial=f$initial,n_high=sum(high),n_low=sum(!high),
  nonzero_news_high=sum(high & percent$X1[,1]!=0),
  nonzero_news_low=sum(!high & percent$X1[,1]!=0),
  inverse_row_max_error=apply(abs(error),1L,max),
  shock_interaction_correlation=cor(d$X[,1],d$X[,2]),
  main_kkt=f$main_kkt,node_kkt=f$node_diagnostics$kkt)
}
M <- do.call(cbind,lapply(fits,`[[`,"scores")); centered <- sweep(M,2L,colMeans(M),"-")
estimate <- unlist(lapply(fits,`[[`,"estimate"),use.names=FALSE)
ell <- te_bandwidth(centered)
starts <- seq.int(1L,n,by=ell)
A <- t(vapply(starts,function(first) colSums(centered[first:min(n,first+ell-1L),,drop=FALSE]),numeric(ncol(M))))
omega <- sqrt(colSums(A^2)/n); Z <- sqrt(n)*estimate/omega
xi <- own_normal_draws(nrow(A),key$B,key$bootstrap_seed)
draw_ids <- unique(as.integer(round(seq(1,key$B,length.out=25L))))
recomputed <- vapply(draw_ids,function(b) max(abs(colSums(A*xi[,b])/(sqrt(n)*omega))),numeric(1))
replay_error <- max(abs(recomputed-result$boot_stats[draw_ids]))
if(verify) {
 replay <- te_block_bootstrap(M,estimate,ell,B=key$B,seed=key$bootstrap_seed,chunk_size=137L)
 replay_error <- max(replay_error,abs(replay$boot_stats-result$boot_stats))
}
diag <- do.call(rbind,diagnostics); diag$Z <- Z; diag$long_run_sd <- omega
stopifnot(null_error<1e-8,max(errors)<1e-8,null$kkt<1e-4,
 max(diag$main_kkt)<1e-4,max(diag$node_kkt)<1e-4,
 ell==result$block_length,max(abs(omega-result$omega))<1e-10,
 max(abs(Z-result$Z))<1e-10,replay_error<1e-9,
 result$pval==(1+sum(result$boot_stats>=result$Tn))/(key$B+1),
 result$critical_value==unname(quantile(result$boot_stats,.95,type=1)))
peak <- diag[which.max(abs(Z)),,drop=FALSE]
bootstrap_table <- data.frame(range=label,candidates=length(taus),coordinates=ncol(M),
 selected_tau=100*selections$news$result$selected$tau,lambda=lambda,
 Tn=result$Tn,critical_95=result$critical_value,
 critical_99=unname(quantile(result$boot_stats,.99,type=1)),
 p_value=result$pval,B=result$B,block_length=ell,block_count=nrow(A),
 max_at_tau=peak$tau,max_coefficient=peak$coefficient,
 peak_high_n=peak$n_high,peak_high_nonzero_news=peak$nonzero_news_high,
 peak_inverse_row_error=peak$inverse_row_max_error,
 all_coordinates_max_inverse_error=max(diag$inverse_row_max_error))
write.csv(diag,file.path(directory,"bootstrap_diagnostics.csv"),row.names=FALSE)
write.csv(peak,file.path(directory,"maximum_diagnostic.csv"),row.names=FALSE)
write.csv(bootstrap_table,file.path(out,"bootstrap_summary.csv"),row.names=FALSE)
check <- list(passed=TRUE,full_draw_replay=verify,null_error=null_error,
 max_fit_error=max(errors),bootstrap_replay_error=replay_error,
 max_main_KKT=max(diag$main_kkt),max_node_KKT=max(diag$node_kkt),
 identical_calibration_selection_test_grid=identical(selections$news$result$selected$profile$tau,taus/100))
# Floating-point unit conversion is verified by states and a tight numerical tolerance.
stopifnot(max(abs(100*selections$news$result$selected$profile$tau-taus))<1e-12)
capture.output(check,file=file.path(out,if(verify) "replay_verification.txt" else "verification.txt"))
capture.output(sessionInfo(),file=file.path(out,"session_info.txt"))
print(bootstrap_table,row.names=FALSE)
anchor_log("COMPLETE",label,"threshold selection and calibrated block test",file=logfile)
