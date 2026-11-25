# Anchored re-estimation of the original per-city tipping application.
# GT version: beta is the slope for Q <= tau; delta is high minus low.
# Original tract-quality rules, controls, city seeds, two intercepts,
# mclapply, city RDS checkpoints and combined CSV outputs are retained.
# New sample: Q <= .60 and complete cases, then at least 100 tracts/city-period.
# Search [0,.50] with at least 15% of final city observations on each side;
# also exclude empty regimes and constant full-sample design columns.
# Anchor is the midpoint of the resulting admissible interval.
# Initial penalty: one common lambda, selected by 10-fold prediction CV.
# Anchor: midpoint of the trimmed per-city feasible search interval.
# rho_n = 2*c_rho*lambda_glmnet; c_rho defaults to 1.
# The copied kernel always uses nodewise GIC; intercept residualization is
# made consistent between fitting, debiasing and the threshold-effect test.
# See README.md for the current settings.

# 0. Setup ----
rm(list = ls())
library(glmnet)
library(parallel)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this script with Rscript.")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
PROJECT_DIR <- dirname(SCRIPT_DIR)
source(file.path(SCRIPT_DIR, "myfnc_anchor_gt.R"))
source(file.path(SCRIPT_DIR, "common_lambda_cv_gt.R"))
source(file.path(SCRIPT_DIR, "sample_search_rules_gt.R"))

BASE_SEED <- 42
set.seed(BASE_SEED)

N_CORES <- as.integer(Sys.getenv("TIPPING_N_CORES", "4"))
stopifnot(is.finite(N_CORES), N_CORES >= 1)

# Deterministic per-(decade,MSA) seed (2026-08-26) -- a city's result must
# depend only on its own identity, not on N_CORES, mclapply's scheduling, or
# which other cities happen to run in the same batch. A global set.seed()
# plus mclapply's automatic L'Ecuyer-CMRG substreams is NOT enough for that:
# substreams are handed out by WORKER SLOT as tasks are scheduled, so the
# same city can land on a different substream depending on batch composition
# -- confirmed empirically (Albuquerque run alone vs. alongside 3 other
# cities gave different tau_hat with that approach). Calling set.seed() with
# this city-specific value at the top of process_one_city instead pins every
# stochastic step for that city (cv.glmnet's fold assignment inside
# Threshold_Lasso/Threshold_debiased_Lasso/Linear_Lasso, and the Gaussian
# multiplier bootstrap in threshold_effect_test) to a value that depends only
# on (dec, m) -- reproducible regardless of cores, scheduling, or which other
# cities are run alongside it.
city_seed <- function(dec, m) {
  dec_id <- if (dec == "1990-2000") 1L else 2L
  BASE_SEED * 1000000L + dec_id * 100000L + as.integer(m)
}

data_dir   <- paste0(file.path(PROJECT_DIR, "data"), "/")
output_dir <- paste0(Sys.getenv("TIPPING_OUTPUT_DIR", file.path(PROJECT_DIR, "results")), "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- paste0(output_dir, "run_log_percity_", Sys.Date(), ".log")
cat(format(Sys.time()), "- starting tipping_DTLasso_percity_anchor_gt.R\n", file = log_file, append = TRUE)
log_msg <- function(...) {
  msg <- paste(format(Sys.time()), "-", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# 1. Data ----
raw <- read.csv(unz(paste0(data_dir, "analysis_1990_2010_tract2000.zip"), "analysis_1990_2010_tract2000.csv"))
raw$card_sample_flag <- raw$card_sample_flag == "True"

msa_names <- read.stata_names <- tryCatch({
  library(haven)
  msa_def <- as.data.frame(read_dta(file.path(data_dir, "msa_def.dta")))
  setNames(as.character(msa_def$place), as.character(msa_def$msar))
}, error = function(e) { log_msg("no msa name lookup available:", conditionMessage(e)); NULL })

core_controls <- c(
  "unemployment_rate", "log_mean_family_income", "public_transit_share",
  "vacancy_share", "renter_share", "single_unit_share"
)
extra_controls <- c(
  "earnings_receipt_share", "wage_receipt_share", "social_security_receipt_share",
  "public_assistance_receipt_share", "retirement_receipt_share",
  "log1p_mean_wage_income_among_recipients_1999usd",
  "log1p_mean_social_security_income_among_recipients_1999usd",
  "log1p_mean_public_assistance_income_among_recipients_1999usd",
  "log1p_mean_retirement_income_among_recipients_1999usd",
  "lacking_kitchen_share", "lacking_plumbing_share", "pays_extra_utilities_share",
  "english_only_share", "commute_car_share", "commute_motorcycle_share",
  "commute_bicycle_share", "commute_walked_share", "commute_worked_home_share",
  "log_tract_land_area_sqkm", "poverty_low_under_100_share", "poverty_high_175_plus_share",
  "travel_low_under_30_share", "travel_high_60_plus_share",
  "departure_early_before_0700_share", "departure_late_0900_plus_share",
  "householder_age_low_15_34_share", "householder_age_high_65_plus_share",
  "structure_low_density_share", "structure_middle_density_share", "structure_high_density_share",
  "housing_age_low_0_20_years_share", "housing_age_high_41_plus_years_share",
  "move_duration_low_0_5_years_share", "move_duration_high_21_plus_years_share",
  "vehicles_low_0_1_share", "vehicles_high_4_plus_share",
  "heating_utility_gas_share", "heating_bottled_gas_share", "heating_electricity_share",
  "heating_fuel_oil_share", "heating_coal_coke_share", "heating_wood_share",
  "heating_solar_share", "heating_no_fuel_share"
)
control_vars <- c(core_controls, extra_controls)  # 50 Lasso candidates

dep_var   <- "outcome_white_change"
thres_var <- "forcing_minority_share"
msa_var   <- "metro_1999_msa_pmsa"

# 2. Common-CV, anchor and Q-restricted sample settings ----
N_FOLDS <- as.integer(Sys.getenv("TIPPING_N_FOLDS", "10"))
NLAMBDA <- as.integer(Sys.getenv("TIPPING_NLAMBDA", "100"))
BOOTSTRAP_B <- as.integer(Sys.getenv("TIPPING_BOOTSTRAP_B", "499"))
C_RHO <- as.numeric(Sys.getenv("TIPPING_C_RHO", "1"))
MIN_NONZERO_OBS <- 30
TAU_UPPER <- 0.50
MIN_CITY_OBS <- 100L
MIN_REGIME_SHARE <- .15
OUTPUT_SUFFIX <- "_q60_min100_trim15_anchor_gt"
RESUME_RUN <- tolower(Sys.getenv("TIPPING_RESUME", "true")) == "true"
MSA_FILTER <- Sys.getenv("TIPPING_MSAS", "")
DECADE_FILTER <- Sys.getenv("TIPPING_DECADES", "")
stopifnot(N_FOLDS>=4, NLAMBDA>=2, BOOTSTRAP_B>=1, is.finite(C_RHO), C_RHO>=0)
run_signature <- list(method=DTLASSO_METHOD_VERSION,
  threshold_rule=DTLASSO_THRESHOLD_RULE,
  source_md5=unname(tools::md5sum(file.path(SCRIPT_DIR,
    c("tipping_DTLasso_percity_anchor_gt.R","myfnc_anchor_gt.R","common_lambda_cv_gt.R","sample_search_rules_gt.R")))),
  data_md5=unname(tools::md5sum(file.path(data_dir,"analysis_1990_2010_tract2000.zip"))),
  nfolds=N_FOLDS,nlambda=NLAMBDA,bootstrap_B=BOOTSTRAP_B,c_rho=C_RHO,
  anchor="admissible_interval_midpoint",tau_lower=0,tau_upper=TAU_UPPER,
  sample_q_upper=.60,min_city_obs=MIN_CITY_OBS,
  min_regime_share=MIN_REGIME_SHARE,
  grid_filter="each_regime_at_least_15pct_and_no_constant_full_design_columns",
  min_nonzero=MIN_NONZERO_OBS,base_seed=BASE_SEED)
log_msg("Common-lambda CV:",N_FOLDS,"folds;",NLAMBDA,"penalties; c_rho =",C_RHO,
        "; anchor = feasible interval midpoint; bootstrap B =",BOOTSTRAP_B,
        "; cores =",N_CORES,"; resume =",RESUME_RUN)
log_msg("Sample: original tract-quality flags, Q <= 0.60, complete cases, then n >= 100; search [0,0.50]; each regime >= ceiling(0.15*n) tracts.")
log_msg("Regimes: low Q <= tau, high Q > tau; beta = low-group slope; delta = high minus low.")

# 3. Fit one city ----
fit_one_city <- function(d, control_vars) {
  dep   <- d[[dep_var]]
  thres <- d[[thres_var]]

  # drop controls that are effectively constant in THIS city (see header note)
  nonzero_count <- apply(as.matrix(d[, control_vars]), 2, function(x) sum(abs(x) > 1e-8))
  city_control_vars <- control_vars[nonzero_count >= MIN_NONZERO_OBS]
  dropped <- setdiff(control_vars, city_control_vars)

  reg <- as.matrix(d[, city_control_vars])   # NO manual constant column --
                                              # intercept1=TRUE handles beta0
                                              # internally via demeaning
  colnames(reg) <- city_control_vars
  nobs <- nrow(reg)

  if(nobs < MIN_CITY_OBS)stop("Fewer than 100 complete tracts after Q <= 0.60 filtering.")
  search <- q60_candidate_grid(reg,thres,min_regime_share=MIN_REGIME_SHARE)
  tau_lo <- search$tau_lo; tau_hi <- search$tau_hi
  tau_star <- search$tau_star; grid.tau <- search$grid
  thdvar <<- thres
  started <- proc.time()[["elapsed"]]
  cvfit <- anchor_common_cv(reg, dep, thres, grid.tau, tau_star,
    c_rho=C_RHO,nfolds=N_FOLDS,nlambda=NLAMBDA,
    intercept1=TRUE,intercept2=TRUE)
  initial <- Threshold_Lasso(reg, dep, taugrid=grid.tau,
    lambda=cvfit$lambda,intercept1=TRUE,intercept2=TRUE,
    tau_star=tau_star,c_rho=C_RHO,profile=cvfit$profile)
  desla_fit <- Threshold_debiased_Lasso(reg,dep,taugrid=grid.tau,
    lambda=cvfit$lambda,intercept1=TRUE,intercept2=TRUE,
    alpha=.05,verbose=FALSE,tau_star=tau_star,c_rho=C_RHO,lasso_result=initial)
  estimate_seconds <- proc.time()[["elapsed"]]-started
  # Same cached full-sample candidates; null linear model retains 10-fold CV.
  # The test includes the threshold intercept and all delta coefficients.
  lintest <- threshold_effect_test(reg,dep,taugrid=grid.tau,thdvar=thres,
    lambda=cvfit$lambda,lambda0="cv",nfolds=N_FOLDS,B=BOOTSTRAP_B,
    alpha=.05,intercept1=TRUE,intercept2=TRUE,ncores=1,verbose=FALSE,
    candidate_profile=cvfit$profile)
  elapsed_seconds <- proc.time()[["elapsed"]]-started

  colnames_new <- c("intercept", "threshold_intercept", city_control_vars,
                     paste0("threshold_", city_control_vars))
  # Existing column names/order retained: plain variables are low-group
  # slopes; threshold_* are HIGH-minus-LOW slope/intercept differences.
  # save both the initial (pre-debiasing) LASSO point estimate and the
  # debiased coefficient + its SE/p-value, so the two can be compared.
  # betaint's (intercept's) noise.sd/pvals are NA -- point estimate only, no
  # valid SE for the reasons in INTERCEPT_DEMEAN_STANDARD_ERROR_NOTES.md.
  desla <- data.frame(
    variable = colnames_new,
    lasso_coefficient = as.numeric(desla_fit$coef),
    coefficient = desla_fit$unb.coef,
    sd = desla_fit$noise.sd,
    p_value = desla_fit$pvals,
    low_wald = desla_fit$low.lim.wald,
    up_wald = desla_fit$up.lim.wald,
    low_widened = desla_fit$low.lim.widened,
    up_widened = desla_fit$up.lim.widened,
    addlength = desla_fit$addlength
  )

  list(tau_hat=desla_fit$best_tau,lambda=desla_fit$lambda_val,n=nobs,desla=desla,
       dropped_controls=dropped,n_grid_tau=length(grid.tau),
       lintest_Tn=lintest$Tn,lintest_pval=lintest$pval,lintest_reject=lintest$reject,
       tau_star=tau_star,tau_lo=tau_lo,tau_hi=tau_hi,c_rho=C_RHO,
       lambda_paper=2*cvfit$lambda,rho_n=2*C_RHO*cvfit$lambda,
       cv_boundary=cvfit$cv$selected_at_boundary,cv_mse=min(cvfit$cv$table$cv_mse),
       estimate_seconds=estimate_seconds,elapsed_seconds=elapsed_seconds,
       search_audit=search$audit,anchor_check=search$anchor_check,
       min_regime_share=search$min_regime_share,min_regime_obs=search$min_regime_obs,
       n_candidates_excluded=search$excluded,
       profile=cvfit$profile,cv=cvfit$cv,initial=initial,
       test_details=lintest,fit_details=desla_fit)
}

# 5b. Per-city worker ----
# Runs under mclapply -- must not write to the SHARED output CSVs or cat()
# directly (forked children's writes/console output interleave
# unpredictably). Each city's own .rds (unique filename per city) has no
# such race, though, so this saves one immediately after that city finishes
# -- if the run crashes partway through a decade, every city that had
# already finished is still on disk, not just whatever mclapply had
# returned to the parent before the crash. Also returns coefs/summary_row
# (NULL on failure) plus a pre-formatted log line; the main process
# writes/logs these after mclapply() returns, in MSA order, same as before.
process_one_city <- function(m, pool, dec, rds_dir, resume = FALSE) {
  rds_path <- file.path(rds_dir, sprintf("msa_%s.rds", m))
  if (resume && file.exists(rds_path)) {
    cached <- tryCatch(readRDS(rds_path),error=function(e)NULL)
    if (!is.null(cached$summary_row) && identical(cached$run_signature,run_signature)) {
      cached$log_line <- paste("RESUMED:",cached$log_line)
      return(cached)
    }
  }

  set.seed(city_seed(dec, m))   # see the city_seed() comment in section 0

  # Each worker writes only its own diagnostic file; the original city
  # logs, combined outputs, error handling and RDS resume remain intact.
  audit_dir <- file.path(Sys.getenv("TIPPING_CHECK_DIR", file.path(PROJECT_DIR, "checks")),
                         "nodewise_paths")
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  audit_path <- file.path(audit_dir, sprintf("%s_%s.csv", dec, m))
  if (file.exists(audit_path)) file.remove(audit_path)
  old_audit_options <- options(tipping.nodewise.audit_path = audit_path,
    tipping.nodewise.city = paste(dec, "MSA", m))
  on.exit(options(old_audit_options), add = TRUE)

  d_city <- pool[pool[[msa_var]] == m, ]
  n_city <- nrow(d_city)
  name <- if (!is.null(msa_names)) msa_names[as.character(as.integer(m))] else NA
  if (is.na(name)) name <- as.character(m)

  res <- tryCatch(fit_one_city(d_city, control_vars), error = function(e) {
    list(error = sprintf("%s: MSA %s (%s, n=%d) FAILED: %s",
                          dec, m, name, n_city, conditionMessage(e)))
  })
  if (!is.null(res$error)) {
    out <- list(coefs = NULL, summary_row = NULL, log_line = res$error, run_signature=run_signature)
    saveRDS(out, rds_path)
    return(out)
  }

  coefs <- res$desla
  coefs$msa <- m
  coefs$msa_name <- name
  coefs$n_tracts <- res$n
  coefs$n_tracts_original <- n_city
  coefs$tau_hat <- res$tau_hat
  coefs$lambda <- res$lambda
  coefs$tau_star <- res$tau_star
  coefs$c_rho <- res$c_rho
  coefs$rho_n <- res$rho_n
  coefs$lambda_paper <- res$lambda_paper
  coefs$min_regime_share <- res$min_regime_share
  coefs$min_regime_obs <- res$min_regime_obs
  coefs$method_version <- DTLASSO_METHOD_VERSION
  coefs$threshold_rule <- DTLASSO_THRESHOLD_RULE

  n_sig <- sum(coefs$p_value < 0.1, na.rm = TRUE)
  summary_row <- data.frame(msa = m, msa_name = name, n_tracts = res$n,
                             n_tracts_original = n_city,
                             tau_hat = res$tau_hat, lambda = res$lambda,
                             n_significant_10pct = n_sig,
                             n_grid_tau = res$n_grid_tau,
                             n_dropped_controls = length(res$dropped_controls),
                             lintest_Tn = res$lintest_Tn,
                             lintest_pval = res$lintest_pval,
                             lintest_reject = res$lintest_reject,
                             tau_star=res$tau_star,tau_lo=res$tau_lo,tau_hi=res$tau_hi,
                             min_regime_share=res$min_regime_share,min_regime_obs=res$min_regime_obs,
                             c_rho=res$c_rho,rho_n=res$rho_n,lambda_paper=res$lambda_paper,
                             cv_mse=res$cv_mse,cv_boundary=res$cv_boundary,
                             estimate_seconds=res$estimate_seconds,elapsed_seconds=res$elapsed_seconds,
                             method_version=DTLASSO_METHOD_VERSION,
                             threshold_rule=DTLASSO_THRESHOLD_RULE)

  log_line <- sprintf("%s: MSA %s (%s, n=%d, grid=%d): tau_hat=%.4f lambda=%.5f, %d/%d sig, lintest reject=%s pval=%.4f%s",
                       dec, m, name, res$n, res$n_grid_tau, res$tau_hat, res$lambda,
                       n_sig, nrow(coefs), res$lintest_reject, res$lintest_pval,
                       if (length(res$dropped_controls) > 0)
                         paste0(" [dropped: ", paste(res$dropped_controls, collapse = ", "), "]") else "")

  log_line <- paste0(log_line,sprintf(" | anchor=%.5f rho_n=%.6g c_rho=%.3g CV=%.6g time=%.1fs",
    res$tau_star,res$rho_n,res$c_rho,res$cv_mse,res$elapsed_seconds))
  out <- list(coefs=coefs,summary_row=summary_row,log_line=log_line,
              run_signature=run_signature,details=res)
  saveRDS(out, rds_path)
  out
}

# 6. Main loop: decade x MSA ----
decades <- c("1990-2000", "2000-2010")
if (nzchar(DECADE_FILTER)) decades <- intersect(decades,strsplit(DECADE_FILTER,",",fixed=TRUE)[[1]])
if (!length(decades)) stop("No selected decade matches the data.")

for (dec in decades) {
  dec_tag <- gsub("-", "_", dec)
  pool <- q60_city_pool(raw,dec,c(dep_var,thres_var,msa_var,control_vars),thres_var)
  msas <- sort(unique(pool[[msa_var]]))
  if (nzchar(MSA_FILTER)) msas <- intersect(msas,as.numeric(strsplit(MSA_FILTER,",",fixed=TRUE)[[1]]))
  if (!length(msas)) stop("No selected MSA matches the data.")
  city_counts <- table(pool[[msa_var]])
  excluded_msas <- msas[city_counts[as.character(msas)] < MIN_CITY_OBS]
  for(m in excluded_msas)log_msg(dec,": MSA",m,"excluded after Q <= 0.60 and complete cases: n =",city_counts[as.character(m)])
  msas <- setdiff(msas,excluded_msas)
  if(!length(msas)){log_msg(dec,": no selected city has at least 100 usable tracts.");next}
  pool <- pool[pool[[msa_var]] %in% msas,,drop=FALSE]
  log_msg(sprintf("%s: %d MSAs, %d tracts total", dec, length(msas), nrow(pool)))

  coefs_path   <- paste0(output_dir, dec_tag, OUTPUT_SUFFIX, "_percity_highdim_coefs.csv")
  summary_path <- paste0(output_dir, dec_tag, OUTPUT_SUFFIX, "_percity_highdim_summary.csv")
  # fresh files for this decade -- written once after all cities in the
  # decade finish (see mclapply note below)
  if (file.exists(coefs_path)) file.remove(coefs_path)
  if (file.exists(summary_path)) file.remove(summary_path)

  # per-city .rds crash-safety net -- each forked worker writes its own city
  # here the moment it finishes (see process_one_city). Fresh per decade
  # unless RESUME_RUN=TRUE, in which case existing .rds files are kept and
  # process_one_city() loads from them instead of refitting (see section 2).
  rds_dir <- paste0(output_dir, "rds/", dec_tag, OUTPUT_SUFFIX, "/")
  if (RESUME_RUN) {
    dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
    n_cached <- length(list.files(rds_dir, pattern = "^msa_.*\\.rds$"))
    log_msg(sprintf("%s: RESUME_RUN=TRUE, %d cities already cached in %s", dec, n_cached, rds_dir))
  } else {
    unlink(rds_dir, recursive = TRUE)
    dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # cities are fully independent (no pooling/FE), so fit them in parallel
  # across N_CORES forked workers; results come back in msas order regardless
  # of completion order. The two combined CSVs below are still written once,
  # for the whole decade, after all cities finish -- but every completed
  # city's result is already safe in rds_dir/ as it finishes (see
  # process_one_city), so a crash mid-decade no longer loses that decade's
  # progress: rerun and skip cities whose .rds already exists, or reload
  # from rds_dir/ directly with readRDS().
  results <- mclapply(msas, process_one_city, pool = pool, dec = dec, rds_dir = rds_dir,
                       resume = RESUME_RUN, mc.cores = N_CORES)

  if (any(vapply(results,inherits,logical(1),"try-error")))
    stop("A parallel worker failed outside city error handling; completed city RDS files are retained.")
  for (r in results) log_msg(r$log_line)

  coefs_all   <- do.call(rbind, lapply(results, `[[`, "coefs"))
  summary_all <- do.call(rbind, lapply(results, `[[`, "summary_row"))

  write.table(coefs_all, coefs_path, sep = ",", row.names = FALSE, col.names = TRUE)
  write.table(summary_all, summary_path, sep = ",", row.names = FALSE, col.names = TRUE)

  log_msg(sprintf("%s: done (%d/%d MSAs succeeded)", dec, nrow(summary_all), length(msas)))
}

log_msg("done")
