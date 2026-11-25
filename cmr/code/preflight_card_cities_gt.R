args<-commandArgs(FALSE)
root<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
source(file.path(root,'code/sample_search_rules_gt.R'))
check_dir<-Sys.getenv('TIPPING_CHECK_DIR',file.path(root,'results'))
dir.create(check_dir,recursive=TRUE,showWarnings=FALSE)
e<-new.env(parent=baseenv())
keys<-c('core_controls','extra_controls','control_vars','dep_var','thres_var','msa_var','MIN_NONZERO_OBS','MIN_CITY_OBS','MIN_REGIME_SHARE')
for(x in parse(file.path(root,'code/tipping_DTLasso_percity_anchor_gt.R')))
 if(is.call(x)&&identical(x[[1]],as.name('<-'))&&is.symbol(x[[2]])&&as.character(x[[2]]) %in% keys)eval(x,e)
msa_filter<-Sys.getenv('TIPPING_MSAS','')
raw<-read.csv(unz(file.path(root,'data/analysis_1990_2010_tract2000.zip'),'analysis_1990_2010_tract2000.csv'))
names<-as.data.frame(haven::read_dta(file.path(root,'data/msa_def.dta')))
out<-list();grids<-list()
for(dec in c('1990-2000','2000-2010')){
 pool<-q60_city_pool(raw,dec,c(e$dep_var,e$thres_var,e$msa_var,e$control_vars),e$thres_var)
 ids<-sort(unique(pool[[e$msa_var]]))
 if(nzchar(msa_filter))ids<-intersect(ids,as.numeric(strsplit(msa_filter,',',fixed=TRUE)[[1]]))
 for(m in ids){
  d<-pool[pool[[e$msa_var]]==m,,drop=FALSE]
  X<-as.matrix(d[,e$control_vars,drop=FALSE]);q<-d[[e$thres_var]]
  X<-X[,colSums(abs(X)>1e-8)>=e$MIN_NONZERO_OBS,drop=FALSE]
  search<-if(nrow(d)>=e$MIN_CITY_OBS)tryCatch(q60_candidate_grid(X,q,min_regime_share=e$MIN_REGIME_SHARE),error=function(e)list(error=conditionMessage(e))) else list(error='fewer_than_100_after_Q_filter')
  ok<-is.null(search$error)
  out[[length(out)+1L]]<-data.frame(decade=dec,msa=m,city=as.character(names$place[match(m,names$msar)]),
   n=nrow(d),p=ncol(X),n_grid=if(ok)length(search$grid) else 0,
   tau_lo=if(ok)search$tau_lo else NA,tau_hi=if(ok)search$tau_hi else NA,
   tau_star=if(ok)search$tau_star else NA,
   min_regime_share=e$MIN_REGIME_SHARE,min_regime_obs=ceiling(e$MIN_REGIME_SHARE*nrow(d)),
   n_excluded_candidates=if(ok)search$excluded else NA,
   sample_eligible=nrow(d)>=e$MIN_CITY_OBS,status=if(ok)'ready' else search$error)
  if(ok){z<-search$audit;z$decade<-dec;z$msa<-m;grids[[length(grids)+1L]]<-z}
 }
}
z<-do.call(rbind,out)
write.csv(z,file.path(check_dir,'sample_audit.csv'),row.names=FALSE)
print(z,row.names=FALSE)
cat('\nReady:',sum(z$status=='ready'),'of',nrow(z),'city-periods. This checks sample sizes, the 15% minimum regime share and nonconstant columns, not validity of inference.\n')
