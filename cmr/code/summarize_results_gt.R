args<-commandArgs(FALSE)
root<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
source(file.path(root,'code/myfnc_anchor_gt.R'))
source(file.path(root,'code/common_lambda_cv_gt.R'))
source(file.path(root,'code/sample_search_rules_gt.R'))
audit<-read.csv(file.path(Sys.getenv('TIPPING_CHECK_DIR',file.path(root,'results')),'sample_audit.csv'))
msa_filter<-Sys.getenv('TIPPING_MSAS','')
decade_filter<-Sys.getenv('TIPPING_DECADES','')
if(nzchar(msa_filter))audit<-audit[audit$msa %in% as.numeric(strsplit(msa_filter,',',fixed=TRUE)[[1]]),]
if(nzchar(decade_filter))audit<-audit[audit$decade %in% strsplit(decade_filter,',',fixed=TRUE)[[1]],]
result_dir<-Sys.getenv('TIPPING_OUTPUT_DIR',file.path(root,'results'))
files<-list.files(file.path(result_dir,'rds'),pattern='[.]rds$',recursive=TRUE,full.names=TRUE)
summaries<-coefs<-statuses<-list()
near<-function(a,b,tol=1e-8)stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
for(path in files){
 r<-readRDS(path)
 dec<-gsub('_','-',substr(basename(dirname(path)),1,9))
 m<-as.integer(sub('msa_([0-9]+)[.]rds','\\1',basename(path)))
 if(!paste(dec,m) %in% paste(audit$decade,audit$msa))next
 statuses[[length(statuses)+1L]]<-data.frame(decade=dec,msa=m,
   status=if(is.null(r$summary_row))'not_estimated' else 'success',message=r$log_line)
 if(is.null(r$summary_row))next
 s<-r$summary_row;d<-r$details;z<-d$profile$candidates
 assert_gt_threshold_object(d$profile)
 assert_gt_threshold_object(d$initial)
 assert_gt_threshold_object(d$fit_details)
 assert_gt_threshold_object(d$test_details)
 s$decade<-dec
 stopifnot(s$n_tracts==audit$n[match(paste(dec,m),paste(audit$decade,audit$msa))],
   identical(r$coefs$variable[1:2],c('intercept','threshold_intercept')),
   nrow(r$coefs)==2*(50-s$n_dropped_controls)+2,
   all(is.finite(r$coefs$coefficient)),all(is.finite(r$coefs$sd[-1])),
   all(r$coefs$sd[-1]>0),is.na(r$coefs$p_value[1]),
   all(is.finite(r$coefs$p_value[-1])),
   all(d$profile$table$lambda_glmnet==s$lambda),
   d$fit_details$inverse_method=='full_nodewise')
 near(s$tau_star,(s$tau_lo+s$tau_hi)/2)
 stopifnot(s$n_tracts>=100,all(d$profile$thdvar<=.60),
           all(d$profile$taugrid>=0 & d$profile$taugrid<=.50))
 # Rebuild the grid from the saved raw q and design, then verify that CV,
 # final selection and the bootstrap test all used this same trimmed grid.
 stopifnot(identical(r$run_signature$min_regime_share,.15),
           s$min_regime_share==.15,s$min_regime_obs==ceiling(.15*s$n_tracts))
 search<-q60_candidate_grid(d$profile$X1,d$profile$thdvar,min_regime_share=.15)
 near(d$profile$taugrid,search$grid)
 near(s$tau_star,search$tau_star)
 near(d$test_details$taugrid,d$profile$taugrid)
 stopifnot(all(d$cv$fold_tau %in% d$profile$taugrid))
 nl<-vapply(d$profile$taugrid,function(tau)sum(d$profile$thdvar<=tau),integer(1))
 stopifnot(all(nl>=s$min_regime_obs),all(s$n_tracts-nl>=s$min_regime_obs))
 near(s$rho_n,2*s$c_rho*s$lambda)
 near(d$cv$table$cv_mse,colSums(d$cv$fold_sse)/s$n_tracts)
 stopifnot(d$cv$selected_index==which.min(d$cv$table$cv_mse))
 near(s$lambda,d$cv$table$lambda_glmnet[d$cv$selected_index])
 tab<-d$profile$table;selected<-max(which(tab$anchored_objective==min(tab$anchored_objective)))
 near(s$tau_hat,tab$tau[selected]);near(tab$Qn,tab$rss/s$n_tracts+2*s$lambda*tab$penalty_l1)
 near(tab$anchored_objective,tab$Qn+s$rho_n*abs(tab$tau-s$tau_star))
 fit<-anchor_cached_candidate(d$profile,selected)
 near(fit$coefficients,d$initial$coefficients)
 near(sum((d$profile$y-cbind(fit$X1,fit$X_tau)%*%as.numeric(fit$coefficients))^2),tab$rss[selected])
 near(r$coefs$p_value[-1],2*pnorm(-abs(r$coefs$coefficient[-1]/r$coefs$sd[-1])))
 near(r$coefs$low_wald[-1],r$coefs$coefficient[-1]-qnorm(.975)*r$coefs$sd[-1])
 near(r$coefs$up_widened[-1]-r$coefs$up_wald[-1],r$coefs$addlength[-1])
 stopifnot(ncol(d$test_details$Z)==51-s$n_dropped_controls,
           length(d$test_details$boot_stats)==499L)
 near(s$lintest_Tn,max(abs(d$test_details$Z)))
 near(s$lintest_pval,(1+sum(d$test_details$boot_stats>=s$lintest_Tn))/500)
 s$tau_no_anchor_same_lambda<-max(tab$tau[tab$Qn==min(tab$Qn)])
 s$anchor_selected<-abs(s$tau_hat-s$tau_star)<1e-12
 s$distance_to_anchor<-abs(s$tau_hat-s$tau_star)
 s$anchor_changes_tau_same_lambda<-abs(s$tau_hat-s$tau_no_anchor_same_lambda)>1e-12
 s$node_kkt_gap<-d$fit_details$node_kkt_gap
 s$gram_diag_error<-d$fit_details$gram_diag_error
 s$n_left<-sum(d$profile$thdvar<=s$tau_hat)
 s$n_right<-sum(d$profile$thdvar>s$tau_hat)
 s$selected_at_grid_edge<-s$tau_hat %in% range(d$profile$taugrid)
 s$n_candidates_excluded<-d$n_candidates_excluded
 c<-r$coefs;c$decade<-dec
 # One family contains all finite coefficient tests in this city-period.
 # The reconstructed baseline intercept has no p-value and is excluded.
 finite_p<-is.finite(c$p_value)
 c$p_bonferroni_city<-rep(NA_real_,nrow(c))
 c$p_bonferroni_city[finite_p]<-p.adjust(c$p_value[finite_p],method='bonferroni')
 s$n_bonferroni_5pct<-sum(c$p_bonferroni_city<=.05,na.rm=TRUE)
 is_delta<-grepl('^threshold_',c$variable)
 s$n_delta_bonferroni_5pct<-sum(c$p_bonferroni_city[is_delta]<=.05,na.rm=TRUE)
 summaries[[length(summaries)+1L]]<-s
 coefs[[length(coefs)+1L]]<-c
}
observed<-if(length(statuses))do.call(rbind,statuses) else
 data.frame(decade=character(),msa=integer(),status=character(),message=character())
status<-merge(audit[,c('decade','msa','n_grid')],observed,
 by=c('decade','msa'),all.x=TRUE,sort=TRUE)
status$status[is.na(status$status)]<-'pending'
 status$status[status$status=='pending' & status$n_grid==0]<-'excluded_by_sample_or_design'
status$status[status$status=='not_estimated' & status$n_grid<2]<-'infeasible_original_grid'
status$status[status$status=='not_estimated']<-'estimation_error'
write.csv(status,file.path(result_dir,'run_status.csv'),row.names=FALSE)
if(!length(summaries))stop('No successfully estimated cities yet.')
s<-do.call(rbind,summaries);co<-do.call(rbind,coefs)
s<-s[order(s$decade,s$msa),]
s$global_p_bonferroni_within_decade<-ave(s$lintest_pval,s$decade,FUN=function(p)p.adjust(p,'bonferroni'))
write.csv(s,file.path(result_dir,'all_city_summary_audited.csv'),row.names=FALSE)
write.csv(co,file.path(result_dir,'all_city_coefficients_audited.csv'),row.names=FALSE)
overview<-do.call(rbind,lapply(split(s,s$decade),function(a){
 expected<-audit[audit$decade==a$decade[1],]
 data.frame(decade=a$decade[1],total_cities=nrow(expected),feasible_cities=sum(expected$n_grid>=2),
   successful_cities=nrow(a),tracts=sum(a$n_tracts),
   tau_median=median(a$tau_hat),tau_q25=quantile(a$tau_hat,.25),tau_q75=quantile(a$tau_hat,.75),
   anchor_median=median(a$tau_star),exact_anchor=sum(a$anchor_selected),
   anchor_changes_tau_same_lambda=sum(a$anchor_changes_tau_same_lambda),
   global_rejections_5pct=sum(a$lintest_pval<=.05),
   global_bonferroni_rejections_5pct=sum(a$global_p_bonferroni_within_decade<=.05),
   median_significant_10pct=median(a$n_significant_10pct),
   median_bonferroni_5pct=median(a$n_bonferroni_5pct),median_delta_bonferroni_5pct=median(a$n_delta_bonferroni_5pct),
   cv_boundary_count=sum(a$cv_boundary),median_elapsed=median(a$elapsed_seconds),
   max_gram_diag_error=max(a$gram_diag_error))
}))
cat('All saved successful cities passed the coefficient, CV, anchor, interval, and bootstrap audits.\n')
print(overview,row.names=FALSE)
if(nrow(status)&&any(status$status!='success')){
 cat('\nRun status (including pending and infeasible cities):\n');print(table(status$status))
 problems<-status$status %in% c('infeasible_original_grid','estimation_error')
 if(any(problems))print(status[problems,],row.names=FALSE)
}
cat('\nSuccessful city records:',nrow(s),'of',sum(audit$n_grid>=2),'with feasible grids.\n')
