# Portable entry point; all inputs are inside this folder.
extension_setup <- function(directory) {
 root <- normalizePath(directory)
 files <- c("config.R","lp.R","design.R","threshold_test.R",
   "threshold_plugin.R","threshold_plugin_anchor.R","threshold_node_recovery.R")
 for(f in files) sys.source(file.path(root,f),envir=.GlobalEnv)
 bundle <- readRDS(file.path(dirname(root),"data","data.rds"))
 sample <- ec_search_sample(bundle$data,bundle$controls,settings$lags,"y",100*settings$anchor)
 ec_verify_alignment(sample,bundle$data,bundle$controls,settings$lags,"y")
 stopifnot(nrow(sample$X1)==500L,ncol(sample$X1)==57L,length(bundle$controls)==14L)
 tf <- tempfile(); saveRDS(list(data=bundle$data,controls=bundle$controls),tf,version=2)
 data_md5 <- unname(tools::md5sum(tf)); unlink(tf)
 list(root=root,ext=root,bundle=bundle,sample=sample,data_md5=data_md5,
   sources=c(file.path(root,files),file.path(root,"common.R")))
}
anchor_tau_grid <- function(q,tau_star) {
 limits <- as.numeric(quantile(q,c(settings$trim,1-settings$trim),type=7))
 grid <- sort(unique(q[q>=limits[1]&q<=limits[2]]))
 if(tau_star>=limits[1]&&tau_star<=limits[2]) grid<-sort(unique(c(grid,tau_star)))
 grid
}
extension_install_lp_recovery <- function() {
  tn_install_recovery()
  original <- own_plugin_lasso
  assign("own_plugin_lasso", function(X, y, unpenalized = integer(),
      cfg = own_lp_control(), seed = cfg$seed) {
    fit <- original(X, y, unpenalized, cfg, seed)
    mode <- "adaptive"
    if (!fit$converged) {
      precise <- cfg; precise$pi_draws <- 10000L
      fit <- original(X, y, unpenalized, precise, seed)
      mode <- "larger_simulation"
      if (!fit$converged) {
        precise$hac_Q <- max(tail(fit$trace$Q, 20L))
        fit <- original(X, y, unpenalized, precise, seed)
        mode <- "fixed_bandwidth_with_existing_recovery"
      }
    }
    if (!fit$converged) stop("LP tuning did not converge: ", seed)
    fit$tuning_mode <- mode
    fit
  }, envir = .GlobalEnv)
}

# Independent lag-sum implementation for verifying the production HAC routine.
extension_hac_check <- function(scores, Q) {
  n <- nrow(scores)
  out <- crossprod(scores) / n
  if (Q > 1L) for (ell in seq_len(Q - 1L)) {
    cross <- crossprod(scores[seq.int(ell + 1L, n), , drop = FALSE],
                       scores[seq_len(n - ell), , drop = FALSE]) / (n - ell)
    out <- out + (1 - ell / Q) * (cross + t(cross))
  }
  out
}
