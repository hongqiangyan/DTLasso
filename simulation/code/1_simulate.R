# ============================================================
# Title: Multi-design Monte Carlo simulation for the threshold-regression
#        paper (Sec. simmo) -- covers the configurations reported in
#        Table tbl:confidence and Table tbl:FWER, extended with the new
#        linearity/threshold-effect test (threshold_effect_test). Runs
#        Threshold_debiased_Lasso, Linear_debiased_Lasso, and
#        threshold_effect_test on the SAME simulated dataset per
#        replication (synchronized, not run separately -- see notes below).
# Author: Hongqiang Yan
# Date: 2026-08-25
# Input: myfnc.R
# Output: results/<design_id>_coef.rds, <design_id>_scalar.rds --
#         ONE pair of files per design (long/tidy data.frames, one row per
#         (rep, coefficient) in _coef.rds, one row per rep in _scalar.rds),
#         not one .rds per replication -- same format Test_Designs.R now
#         uses and Summarize_Designs.R now reads. .rds chosen over CSV:
#         same data, ~3.4x smaller on disk (measured on the smoke test).
#
# Resume behavior (TRUE per-rep resume):
#   - Each rep k within a design is seeded independently:
#     set.seed(BASE_SEED + design_seed_offset(design$id)*1e3 + k), where
#     design_seed_offset() hashes the design's stable id (the original id
#     for the nine renamed designs; see the mapping below) -- NOT its
#     position in the `designs` list. Rep k's simulated data is therefore
#     reproducible no matter what order designs/reps run in, which other
#     reps were already done in a prior run, or whether entries were later
#     added/removed/reordered in `designs` -- unlike the original single
#     sequential set.seed(123) stream (which also meant that running only a
#     subset of designs in a session silently shifted the RNG state seen by
#     later designs relative to a full run) and unlike a naive
#     position-based index (which would reseed every later design the
#     moment you delete/reorder an earlier entry in the list).
#   - On start, an existing _scalar.rds/_coef.rds for a design (if any) is
#     read; already-completed rep numbers are skipped and only the
#     remaining reps run -- safe to invoke this script repeatedly across
#     separate sessions, for any subset of designs, in any order.
#   - New results are checkpointed every CHECKPOINT_EVERY newly-completed
#     reps, merged onto the existing on-disk data (not overwritten from
#     scratch), plus once more at the end of the run. A crash loses at most
#     ~CHECKPOINT_EVERY reps of unsaved work, PLUS whatever's in the batch
#     that was in flight when it died (up to `OUTER_NCORES` reps -- see the
#     two-level core-split note below: checkpointing can only happen at a
#     batch boundary, since mclapply blocks until every rep in the batch
#     returns).
#   - A design is skipped entirely (not even entered) if its _scalar.rds
#     already has >= B_mc distinct rep numbers on disk.
#
# NOTES / OPEN ITEMS (flagged, not silently decided):
#   - rho_QX > 0 (Q correlated with X) is NOT something the original
#     SimulationSyn.R supports (it always drew thdvar independently via
#     runif(n)). The correlation mechanism below (Q built from a shared
#     latent factor with X1's first column, mapped to (0,1) via pnorm) is
#     my own construction -- the paper text doesn't specify the exact DGP
#     for the rho_QX = 0.5 row, so this may not reproduce the exact
#     published numbers for that config. Flag this if you need an exact
#     match to Table tbl:confidence row 5 (ell=0.3433 etc.).
#   - taugrid: fixed seq(0.15,0.85,by=0.01), not build_taugrid() -- true tau_0
#     in {0.4, 0.5} across designs, both covered; build_taugrid() is
#     reserved for real-data application scripts.
#   - Seeding: per-replication -- set.seed(BASE_SEED + design_seed_offset(id)*1e3
#     + k) is drawn fresh inside run_one_rep() (departs from original
#     SimulationSyn.R's single sequential-draw convention). Reproducible
#     for any individual rep in isolation regardless of which batch it
#     lands in, how many reps run alongside it, or how `designs` gets
#     edited later, which is what true per-rep resume requires (see resume
#     note above).
#   - ncores / OUTER_NCORES / INNER_NCORES: two-level split of the total
#     `ncores` core budget (added 2026-08-28, after profiling on Sol showed
#     threshold_effect_test's 71-point tau loop -- not the outer rep loop --
#     is the actual bottleneck within a single rep: each tau needs its own
#     nodewise-Lasso debiasing pass over ~p2 columns). OUTER_NCORES controls
#     how many reps run CONCURRENTLY per batch (parallel::mclapply,
#     fork-based); INNER_NCORES is handed to threshold_effect_test's OWN
#     ncores so run_one_rep() parallelizes ITS internal taugrid loop too,
#     instead of running it on a single core. OUTER_NCORES * INNER_NCORES
#     never exceeds `ncores`, so this nests forks inside forks WITHOUT
#     oversubscribing. Total core-seconds for a full design run is roughly
#     unchanged vs. pure outer-only parallelism (same total work, same total
#     cores) -- the point is cutting wall time for a SINGLE rep, and
#     therefore time-to-first-checkpoint, which pure outer parallelism can't
#     do (you always wait at least one full serial rep before the first
#     batch -- and thus the first checkpoint -- can return, no matter how
#     many cores are available for OTHER reps to run alongside it).
# ============================================================

# 0. Setup ----
library(MASS)           # this script's own dependency (generate_data() uses
                         # mvrnorm) -- myfnc.R declares only what ITS OWN
                         # functions use (Matrix/glmnet/expm/flare), not what
                         # its callers happen to need, since it's shared by
                         # other scripts (Test_Designs.R, Summarize_Designs.R, ...)
source("code/myfnc.R")  # run from the simulation/ folder: Rscript code/1_simulate.R
if (!exists("DTLASSO_METHOD_VERSION", inherits = FALSE) ||
    !identical(DTLASSO_METHOD_VERSION, "common_lambda_anchored_gic_v4_linear_widening"))
  stop("This driver requires myfnc.R from DTLasso_v4_anchor03_linear_widening. Replace BOTH R files together.")

# ---- User-configurable knobs ----
BASE_SEED <- 123         # per-rep seed = BASE_SEED + design_seed_offset(id)*1e3 + k
                          # (1e3 multiplier, not 1e5 -- design_seed_offset()
                          # returns up to ~1e6, and R's set.seed() requires a
                          # 32-bit int (< ~2.1e9); 1e6*1e5 overflows that,
                          # 1e6*1e3 doesn't, and 1e3 is still >> B_mc so
                          # different designs' seed ranges never overlap)

# CLI args (all optional, for cluster submission):
#   arg1 = design id to run alone ("ddd1" through "ddd9");
#          omitted/empty = run all
#   arg2 = B_mc override (e.g. small value for a timing/smoke-test run)
#   arg3 = INNER_NCORES override (cores given to EACH rep's own
#          threshold_effect_test tau-loop; see two-level split below)
#   arg4 = bootstrap_B override (threshold_effect_test's own B; default 499 --
#          use a small value like 20 for a correctness-only smoke test, since
#          this bootstrap x the 71-point taugrid is the actual per-rep
#          bottleneck, not B_mc)
#   arg5 = taugrid step ("by") override (default 0.01 -> 71 points; use e.g.
#          0.1 -> 8 points for a fast correctness check). Affects BOTH
#          taugrid_search and taugrid_test.
#   arg6 = retained for compatibility; only 1 (or empty) is accepted.
#          Alpha and tau now share lambda_paper; c_rho is fixed at 1.
cli_args      <- commandArgs(trailingOnly = TRUE)
design_filter <- if (length(cli_args) >= 1 && nzchar(cli_args[1])) cli_args[1] else NULL

ncores <- {
  slurm_cores <- Sys.getenv("SLURM_CPUS_PER_TASK", "")
  if (nzchar(slurm_cores)) as.integer(slurm_cores) else 8
}                         # <<< TOTAL cores available to this job (from
                          # SLURM's -c/--cpus-per-task when run as a job,
                          # else falls back to 8 for local runs). Split
                          # two ways below -- see OUTER_NCORES/INNER_NCORES.

# ---- Two-level core split ----
# threshold_effect_test()'s own 71-point tau loop (fit_one_tau, called once
# per tau) is the actual bottleneck within a single rep -- each tau needs a
# fresh nodewise-Lasso debiasing pass over ~p2 columns. Previously that loop
# ran on 1 core (ncores=1, fixed) while ALL of `ncores` went to running reps
# concurrently (OUTER only) -- correct, but means you wait for one FULL
# serial rep (the slowest part of the whole script) before seeing even the
# first checkpoint, no matter how many cores you throw at it.
# Splitting cores two ways instead -- INNER_NCORES parallelizes the tau loop
# WITHIN each rep, OUTER_NCORES controls how many reps run concurrently --
# keeps total forked processes at OUTER_NCORES * INNER_NCORES <= ncores (no
# oversubscription) while making the first rep (and thus first checkpoint)
# finish roughly INNER_NCORES times faster. Total core-seconds for a full
# design run is ~unchanged (same total work, same total cores) -- this is
# about latency to first result / ability to monitor progress, not raw
# throughput.
INNER_NCORES <- if (length(cli_args) >= 3 && nzchar(cli_args[3])) {
  as.integer(cli_args[3])
} else {
  min(8L, ncores)         # <<< tune this: higher = faster single-rep wall
                          # time / faster first checkpoint, at the cost of
                          # fewer reps running concurrently
}
OUTER_NCORES <- max(1L, ncores %/% INNER_NCORES)

B_mc   <- 500             # Monte Carlo replications per design
bootstrap_B <- 499        # threshold_effect_test's own bootstrap reps (the
                          # real per-rep cost driver -- see arg4 above)
CHECKPOINT_EVERY <- 10    # merge+save a design's .rds pair every this many NEW reps
taugrid_step   <- 0.01    # see arg5 above
taugrid_search <- seq(0.15, 0.85, by = taugrid_step)
taugrid_test   <- seq(0.15, 0.85, by = taugrid_step)
alpha  <- 0.05
LAMBDA_N <- 100L         # one common decreasing lambda grid per simulated sample
LAMBDA_MIN_RATIO <- NULL # glmnet convention: 0.01 if n < 2p, otherwise 0.0001

output_dir <- "results/"
log_dir    <- "results/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(log_dir, paste0("SimulationSyn_designs_", format(Sys.Date()), ".log"))

log_msg <- function(...) {
  msg <- paste(format(Sys.time()), "-", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# 1. Design configurations ----
# (n, p, s0, b, b1, rho_QX, tau_0) -- "p" here is the per-block dimension
# (paper reports "2p" = total design width = p1+p2, so 2p=600 -> p=300).
# Rows tagged in `source` note which paper table(s) this config feeds.
# ddd1--ddd6 = former dd1--dd6; ddd7 = d8_conf8; ddd8 = d7_conf7;
# ddd9 = d9_fwer4. The original parameters and per-design seeds are retained.
# The paper reports the seven designs below. Designs ddd4 (s0 = 30) and ddd7
# (n = 400) are not reported and were removed; seeds depend only on each
# design's own id, so removing them leaves the other designs unchanged.
designs <- list(
  list(id = "ddd1",           n = 800,  p = 500, s0 = 15, b = 2, b1 = 1,   rho_QX = 0,   tau_0 = 0.5,
       source = "New baseline (n=800, 2p=1000)"),
  list(id = "ddd2",           n = 800,  p = 500, s0 = 15, b = 2, b1 = 0.5, rho_QX = 0,   tau_0 = 0.5,
       source = "New design: weaker threshold effect (b1=0.5)"),
  list(id = "ddd3",           n = 800,  p = 500, s0 = 15, b = 2, b1 = 0,   rho_QX = 0,   tau_0 = 0.5,
       source = "New design: no threshold effect"),
  list(id = "ddd5",           n = 800,  p = 500, s0 = 15, b = 2, b1 = 1,   rho_QX = 0.5, tau_0 = 0.5,
       source = "New design: Q correlated with X (rho_QX=0.5)"),
  list(id = "ddd6",           n = 800,  p = 500, s0 = 15, b = 2, b1 = 1,   rho_QX = 0,   tau_0 = 0.4,
       source = "New design: tau_0=0.4"),
  list(id = "ddd8",           n = 800,  p = 300, s0 = 15, b = 2, b1 = 1,   rho_QX = 0,   tau_0 = 0.5,
       source = "Table confidence row7 (n=800)"),
  list(id = "ddd9",           n = 1000, p = 600, s0 = 15, b = 2, b1 = 1,   rho_QX = 0,   tau_0 = 0.5,
       source = "Table FWER row4 (n=1000, 2p=1200)")
)

if (!is.null(design_filter)) {
  designs <- Filter(function(d) d$id == design_filter, designs)
  if (length(designs) == 0) stop("No design matches id: ", design_filter)
}
if (length(cli_args) >= 2 && nzchar(cli_args[2])) B_mc <- as.integer(cli_args[2])
if (length(cli_args) >= 4 && nzchar(cli_args[4])) bootstrap_B <- as.integer(cli_args[4])
if (length(cli_args) >= 5 && nzchar(cli_args[5])) {
  taugrid_step   <- as.numeric(cli_args[5])
  taugrid_search <- seq(0.15, 0.85, by = taugrid_step)
  taugrid_test   <- seq(0.15, 0.85, by = taugrid_step)
}

# Anchoring parameters; all original CLI positions and runtime knobs stay.
tau_star <- 0.3
c_rho_primary <- if (length(cli_args) >= 6 && nzchar(cli_args[6])) {
  suppressWarnings(as.numeric(cli_args[6]))
} else 1
if (length(c_rho_primary) != 1L || !is.finite(c_rho_primary) || c_rho_primary != 1)
  stop("This version fixes c_rho=1 (rho_n=lambda_paper). Leave argument 6 empty or use 1.")
c_rho_grid <- 1          # retain the existing result/loop structure, with one strength
taugrid_search <- sort(unique(c(taugrid_search, tau_star)))
METHOD_VERSION <- "common_lambda_anchored_gic_v4_linear_widening"

# Statistical settings that must match across replications of one experiment.
# B_mc, core counts, paths and descriptive labels do not change a rep's method.
simulation_settings <- function(design) {
  design_fields <- c("n", "p", "s0", "b", "b1", "rho_QX", "tau_0")
  c(list(method_version = METHOD_VERSION, design_id = design$id),
    lapply(design[design_fields], as.numeric),
    list(base_seed = as.numeric(BASE_SEED), rng_kind = RNGkind(),
      taugrid_search = as.numeric(taugrid_search), taugrid_test = as.numeric(taugrid_test),
      bootstrap_B = as.numeric(bootstrap_B), alpha = as.numeric(alpha),
      tau_star = as.numeric(tau_star), c_rho_primary = as.numeric(c_rho_primary),
      c_rho_grid = as.numeric(c_rho_grid), intercept1 = FALSE, intercept2 = FALSE,
      lambda = "gic", lambda0 = "gic", nodewise_ic = "GIC",
      lambda_tuning = "common_lambda_anchored_gic", lambda_nlambda = LAMBDA_N,
      lambda_grid_rule = "global_kkt_log_grid",
      lambda_min_ratio = if (!is.null(LAMBDA_MIN_RATIO)) LAMBDA_MIN_RATIO else
        if (design$n < 2 * design$p) 0.01 else 0.0001,
      test_lambda = "selected_common_lambda",
      widening_rule = "topk_nodewise_nz"))
}

# Check settings before skipping or appending any saved replications.
check_saved_settings <- function(path, expected_settings) {
  if (!file.exists(path)) return(invisible(NULL))
  s <- readRDS(path)
  if (!nrow(s)) return(invisible(NULL))
  ok <- all(c("method_version", "c_rho", "anchor_results", "stat_settings") %in% names(s)) &&
    all(s$method_version == METHOD_VERSION) && all(s$c_rho == c_rho_primary) &&
    all(vapply(s$anchor_results, function(a) {
      setequal(a$c_rho, c_rho_grid) &&
        identical(as.numeric(a$c_rho[a$model == "threshold"]), as.numeric(c_rho_primary))
    }, logical(1))) &&
    all(vapply(s$stat_settings, function(saved) identical(saved, expected_settings), logical(1)))
  if (!isTRUE(ok)) {
    differences <- if ("stat_settings" %in% names(s)) unique(unlist(lapply(s$stat_settings, function(saved) {
      fields <- union(names(saved), names(expected_settings))
      fields[!vapply(fields, function(field) identical(saved[[field]], expected_settings[[field]]), logical(1))]
    }))) else "missing stat_settings (older result format)"
    if (!length(differences)) differences <- "method/c_rho metadata"
    msg <- paste("Cannot resume", basename(path),
      "because statistical settings differ:", paste(differences, collapse = ", "),
      ". Use a separate run directory with an empty output_designs/. Existing results are unchanged.")
    log_msg(msg)
    stop(msg, call. = FALSE)
  }
  invisible(NULL)
}

# Deterministic per-design seed offset, derived from a stable id (NOT
# from the design's position in the `designs` list above). This is what
# makes it safe to add, remove, comment out, or reorder entries in
# `designs` later -- an existing design's seed (and therefore every rep's
# simulated data) depends only on its stable id string, never on which index
# it happens to sit at. For the ids above, even adjacent offsets are
# separated by 1; multiplying by 1e3 keeps their seed ranges disjoint
# because k only ranges from 1 to B_mc=500. Multiplier is 1e3, not
# 1e5: the offset itself can be up to ~1e6, and set.seed() requires a
# 32-bit int (< ~2.1e9) -- 1e6*1e5 silently overflows that ("supplied seed
# is not a valid integer"), 1e6*1e3 does not, and 1e3 is still far larger
# than B_mc so no two designs' seed ranges can ever collide.
# `%% 1e9 + 7` keeps the running hash from ever exceeding double's exact-
# integer range (2^53) even though R's 32-bit integers would overflow.
design_seed_offset <- function(id) {
  # Preserve the original samples when existing designs are renamed.
  original_ids <- c(ddd1 = "dd1", ddd2 = "dd2", ddd3 = "dd3", ddd4 = "dd4",
                    ddd5 = "dd5", ddd6 = "dd6", ddd7 = "d8_conf8",
                    ddd8 = "d7_conf7", ddd9 = "d9_fwer4")
  if (id %in% names(original_ids)) id <- unname(original_ids[[id]])
  h <- Reduce(function(acc, ch) (acc * 131 + ch) %% 1000000007, utf8ToInt(id), 0)
  as.integer(h %% 1000000)
}

# 2. Data-generating process ----
generate_data <- function(n, p, s0, b, b1, rho_QX, tau_0) {
  cov_XG <- 0.5 ^ abs(row(matrix(1:p, p, p)) - col(matrix(1:p, p, p)))
  X1 <- MASS::mvrnorm(n, rep(0, p), cov_XG)

  if (rho_QX == 0) {
    thdvar <- runif(n)
  } else {
    z <- rho_QX * as.numeric(scale(X1[, 1])) + sqrt(1 - rho_QX^2) * rnorm(n)
    thdvar <- pnorm(z)
  }

  X2 <- X1
  X2[thdvar >= tau_0, ] <- 0

  beta0  <- c(rep(b, s0), rep(0, p - s0))
  delta0 <- c(rep(0, s0), rep(b1, s0), rep(0, p - 2 * s0))
  y <- X1 %*% beta0 + X2 %*% delta0 + rt(n, df = 10)

  list(X1 = X1, thdvar = thdvar, y = as.numeric(y), beta0 = beta0, delta0 = delta0, tau_0 = tau_0)
}

# Unpack one model's [betaint?, deltaint?, beta_1..p, (delta_1..p)] output
# vectors into one row per coefficient. true_delta = NULL for the linear
# model (no delta block at all). true value is NA for betaint/deltaint
# since these designs run intercept1=FALSE/intercept2=FALSE (no DGP
# ground-truth concept for either) -- kept generic so it still works if a
# future design turns them on.
extract_coef_rows <- function(design_id, rep_id, model, coef_raw, coef_deb, se, low, up, pval,
                               intercept1, intercept2, p, true_beta, true_delta = NULL) {
  idx <- 1
  rows <- list()

  if (intercept1) {
    rows[[length(rows)+1]] <- data.frame(
      design_id=design_id, rep=rep_id, model=model, param_type="betaint", param_index=NA_integer_,
      true_value=NA_real_, coef_raw=coef_raw[idx], coef_debiased=coef_deb[idx],
      se=se[idx], low_lim=low[idx], up_lim=up[idx], pval=pval[idx]
    )
    idx <- idx + 1
  }
  if (intercept2) {
    rows[[length(rows)+1]] <- data.frame(
      design_id=design_id, rep=rep_id, model=model, param_type="deltaint", param_index=NA_integer_,
      true_value=NA_real_, coef_raw=coef_raw[idx], coef_debiased=coef_deb[idx],
      se=se[idx], low_lim=low[idx], up_lim=up[idx], pval=pval[idx]
    )
    idx <- idx + 1
  }

  beta_idx <- idx:(idx+p-1)
  rows[[length(rows)+1]] <- data.frame(
    design_id=design_id, rep=rep_id, model=model, param_type="beta", param_index=1:p,
    true_value=true_beta, coef_raw=coef_raw[beta_idx], coef_debiased=coef_deb[beta_idx],
    se=se[beta_idx], low_lim=low[beta_idx], up_lim=up[beta_idx], pval=pval[beta_idx]
  )
  idx <- idx + p

  if (!is.null(true_delta)) {
    delta_idx <- idx:(idx+p-1)
    rows[[length(rows)+1]] <- data.frame(
      design_id=design_id, rep=rep_id, model=model, param_type="delta", param_index=1:p,
      true_value=true_delta, coef_raw=coef_raw[delta_idx], coef_debiased=coef_deb[delta_idx],
      se=se[delta_idx], low_lim=low[delta_idx], up_lim=up[delta_idx], pval=pval[delta_idx]
    )
  }

  do.call(rbind, rows)
}

# Keep the original leading columns and widened low_lim/up_lim aliases.
extract_fit_rows <- function(design, k, model, fit, dat, c_rho = NA_real_, rho_n = NA_real_) {
  linear <- identical(model, "linear")
  truth <- if (linear) dat$linear_target else c(dat$beta0, dat$delta0)
  d <- extract_coef_rows(design$id, k, model, fit$coef, fit$unb.coef, fit$noise.sd,
    fit$low.lim, fit$up.lim, fit$pvals, intercept1 = FALSE, intercept2 = FALSE,
    p = design$p, true_beta = if (linear) dat$linear_target else dat$beta0,
    true_delta = if (linear) NULL else dat$delta0)
  d$c_rho <- c_rho
  d$rho_n <- rho_n
  d$structural_truth <- if (linear) dat$beta0 else c(dat$beta0, dat$delta0)
  d$target <- if (linear) "linear_projection" else "structural_alpha"
  d$param_group <- ifelse(truth == 0, "zero", "nonzero")
  d$pvalue_method <- "wald"
  d$studentized_error <- (fit$unb.coef - truth) / fit$noise.sd
  d$low_lim_wald <- fit$low.lim.wald
  d$up_lim_wald <- fit$up.lim.wald
  d$covered_wald <- fit$low.lim.wald <= truth & fit$up.lim.wald >= truth
  d$interval_length_wald <- fit$up.lim.wald - fit$low.lim.wald
  d$low_lim_widened <- fit$low.lim.widened
  d$up_lim_widened <- fit$up.lim.widened
  d$covered_widened <- fit$low.lim.widened <= truth & fit$up.lim.widened >= truth
  d$interval_length_widened <- fit$up.lim.widened - fit$low.lim.widened
  d$addlength <- fit$addlength
  d$widening_rule <- "topk_nodewise_nz"
  d$lambda_glmnet <- fit$lambda_val
  d$lambda_paper <- fit$lambda_paper
  d$inverse_method <- fit$inverse_method
  d$used_node <- fit$used_node
  d
}

# Merge newly-computed rows onto whatever is already on disk (read fresh
# each call, in case another process/session touched it) and save. Used
# both for periodic checkpoints and the final save of a design's remaining
# reps.
merge_and_save <- function(design_id, new_coef_rows, new_scalar_rows) {
  coef_file   <- file.path(output_dir, paste0(design_id, "_coef.rds"))
  scalar_file <- file.path(output_dir, paste0(design_id, "_scalar.rds"))

  old_coef   <- if (file.exists(coef_file))   readRDS(coef_file)   else NULL
  old_scalar <- if (file.exists(scalar_file)) readRDS(scalar_file) else NULL

  design_coef   <- rbind(old_coef,   do.call(rbind, new_coef_rows))
  design_scalar <- rbind(old_scalar, do.call(rbind, new_scalar_rows))
  rownames(design_coef)   <- NULL
  rownames(design_scalar) <- NULL

  saveRDS(design_coef,   coef_file)
  saveRDS(design_scalar, scalar_file)
  invisible(nrow(design_scalar))
}

# Run ONE replication k of ONE design, start to finish (data draw + all three
# model fits). Self-contained (everything it needs -- design, taugrid_search/
# test, alpha -- is captured from the enclosing design loop's environment via
# closure, which parallel::mclapply's fork inherits automatically). Returns a
# list tagged with k and ok/err so the caller (the main process, never a
# forked worker) can log/accumulate/checkpoint -- workers never touch the log
# file or output files themselves, avoiding any concurrent-write races.
#
# threshold_effect_test's OWN ncores is set to INNER_NCORES here, not the
# outer OUTER_NCORES/ncores -- this call already runs inside an outer-parallel
# worker (one process per rep, up to OUTER_NCORES of them at once);
# INNER_NCORES * OUTER_NCORES <= ncores (see the two-level split note where
# INNER_NCORES/OUTER_NCORES are defined) so this nests forks inside forks
# WITHOUT oversubscribing -- the point being to cut wall time for a single
# rep (and thus time-to-first-checkpoint), since the 71-tau loop inside
# threshold_effect_test, not the outer rep loop, is the actual bottleneck.
run_one_rep <- function(design, k) {
  set.seed(BASE_SEED + design_seed_offset(design$id) * 1000L + k)  # per-rep
                                                     # seed: makes rep k
                                                     # reproducible regardless of
                                                     # run order, batch composition,
                                                     # which other reps were skipped,
                                                     # or later edits to the designs
                                                     # list (id-based, not position-based)
  dat <- generate_data(design$n, design$p, design$s0, design$b, design$b1,
                        design$rho_QX, design$tau_0)
  X1 <- dat$X1
  # Threshold_Lasso resolves `thdvar` by LEXICAL scope (it's a free variable
  # in myfnc.R, sourced at top level -- so it looks in .GlobalEnv), NOT from
  # this function's caller or local frame. A plain `<-` here would make it
  # local to run_one_rep and invisible to Threshold_Lasso ("object 'thdvar'
  # not found" -- caught by the smoke test below). `<<-` writes it into
  # .GlobalEnv instead. Safe under mclapply (fork-based): each forked worker
  # has its own private copy of the whole process, including .GlobalEnv, so
  # concurrent reps assigning thdvar<<- don't race or clobber each other.
  thdvar <<- dat$thdvar
  y <- dat$y
  # All current designs have independent Q or tau_0=0.5 with a symmetric DGP.
  # This is a reporting target, not a change to generate_data() or its RNG.
  dat$linear_target <- dat$beta0 + design$tau_0 * dat$delta0

  t0 <- proc.time()[3]

  tryCatch({
    # One lambda path per tau; select tau WITH anchoring at each lambda,
    # then use outer GIC to choose one common lambda. c_rho is fixed at 1.
    lasso_primary <- Threshold_Lasso(X1, y, taugrid = taugrid_search,
      lambda = 'gic', intercept1 = FALSE, intercept2 = FALSE,
      tau_star = tau_star, c_rho = c_rho_primary,
      nlambda = LAMBDA_N, lambda_min_ratio = LAMBDA_MIN_RATIO)
    lasso_all <- lapply(c_rho_grid, function(c) Threshold_Lasso(X1, y,
      taugrid = taugrid_search, lambda = 'gic', intercept1 = FALSE, intercept2 = FALSE,
      tau_star = tau_star, c_rho = c, profile = lasso_primary$profile,
      nlambda = LAMBDA_N, lambda_min_ratio = LAMBDA_MIN_RATIO))
    selected_taus <- vapply(lasso_all, function(f) f$tau, numeric(1))
    first_at_tau <- which(!duplicated(selected_taus))
    # Complete beta/delta inference only at selected thresholds, once per distinct tau.
    th_unique <- lapply(first_at_tau, function(i) Threshold_debiased_Lasso(X1, y,
      taugrid = taugrid_search, intercept1 = FALSE, lambda = 'gic', intercept2 = FALSE,
      alpha = alpha, verbose = FALSE, lasso_result = lasso_all[[i]]))
    th_all <- th_unique[match(selected_taus, selected_taus[first_at_tau])]
    th <- th_all[[match(c_rho_primary, c_rho_grid)]]
    lin <- Linear_debiased_Lasso(X1, y, intercept1 = FALSE, lambda = 'gic', alpha = alpha, verbose = FALSE)
    lintest <- threshold_effect_test(X1, y, taugrid = taugrid_test, thdvar = thdvar,
                                      lambda = lasso_primary$lambda_val, lambda0 = 'gic', B = bootstrap_B, alpha = alpha,
                                      intercept1 = FALSE, intercept2 = FALSE,
                                      ncores = INNER_NCORES, verbose = FALSE)
    anchor_results <- do.call(rbind, lapply(seq_along(c_rho_grid), function(i) {
      fit <- lasso_all[[i]]; deb <- th_all[[i]]
      data.frame(design_id = design$id, rep = k,
        seed = BASE_SEED + design_seed_offset(design$id) * 1000L + k,
        model = if (c_rho_grid[i] == c_rho_primary) "threshold" else paste0("threshold_c", c_rho_grid[i]),
        c_rho = c_rho_grid[i], n = design$n, p = design$p,
        tau_0 = if (design$b1 == 0) NA_real_ else design$tau_0,
        tau_star = tau_star, th_best_tau = fit$tau, rho_n = fit$rho_n,
        lambda_glmnet = fit$lambda_val, lambda_paper = fit$lambda_paper,
        lambda_paper_anchor = fit$lambda_paper_anchor, tuning = "common_lambda_anchored_gic",
        Qn = fit$Qn, anchored_objective = fit$objective,
        distance_anchor = abs(fit$tau - tau_star),
        distance_truth = if (design$b1 == 0) NA_real_ else abs(fit$tau - design$tau_0),
        inverse_method = deb$inverse_method, used_node = deb$used_node,
        node_kkt_gap = deb$node_kkt_gap, gram_diag_error = deb$gram_diag_error,
        lintest_Tn = lintest$Tn, lintest_pval = lintest$pval, lintest_reject = lintest$reject)
    }))
    list(ok = TRUE, k = k, dat = dat, th = th, lin = lin, lintest = lintest,
         th_all = th_all, anchor_results = anchor_results, profile = lasso_primary$profile$table,
         lambda_search = lasso_primary$profile$lambda_search,
         elapsed = proc.time()[3] - t0)
  }, error = function(e) list(ok = FALSE, k = k, err = conditionMessage(e)))
}

# 3. Run all designs ----
# Safe to re-invoke this script for any subset of designs, in any order,
# across separate sessions -- see the per-rep resume note in the header.
# Reps within a design run in batches of `OUTER_NCORES`, each batch executed
# concurrently via parallel::mclapply (fork-based -- Unix/Mac only, same
# guard convention as threshold_effect_test's own mclapply use). Each of
# those OUTER_NCORES concurrent reps in turn uses INNER_NCORES cores for its
# own threshold_effect_test call -- see the two-level split note above
# run_one_rep(). OUTER_NCORES * INNER_NCORES <= ncores throughout.
log_msg("Method:", METHOD_VERSION, "| tau_star =", tau_star,
        "| c_rho = 1, rho_n = lambda_paper | common lambda candidates =", LAMBDA_N,
        "| outer GIC includes anchored tau selection")
for (design_idx in seq_along(designs)) {
  design <- designs[[design_idx]]
  scalar_file <- file.path(output_dir, paste0(design$id, "_scalar.rds"))
  run_settings <- simulation_settings(design)
  check_saved_settings(scalar_file, run_settings)

  done_ks <- if (file.exists(scalar_file)) unique(readRDS(scalar_file)$rep) else integer(0)
  if (length(done_ks) >= B_mc) {
    log_msg("=== design", design$id, "already complete (", length(done_ks), "reps found), skipping ===")
    next
  }

  remaining_ks <- setdiff(1:B_mc, done_ks)
  log_msg("=== design", design$id, "(", design$source, ") --",
          length(done_ks), "reps already done,", length(remaining_ks), "remaining, batch size", OUTER_NCORES,
          "x", INNER_NCORES, "inner cores/rep ===")

  new_coef_rows      <- list()
  new_scalar_rows    <- list()
  n_since_checkpoint <- 0

  batches <- split(remaining_ks, ceiling(seq_along(remaining_ks) / OUTER_NCORES))

  for (batch_i in seq_along(batches)) {
    batch <- batches[[batch_i]]

    if (OUTER_NCORES > 1 && .Platform$OS.type != "windows") {
      batch_results <- parallel::mclapply(batch, function(k) run_one_rep(design, k),
                                           mc.cores = OUTER_NCORES)
    } else {
      batch_results <- lapply(batch, function(k) run_one_rep(design, k))
    }

    # Collect/log/accumulate sequentially in the main process -- the only
    # process that ever writes the log file or the output .rds files.
    for (res in batch_results) {
      k <- res$k
      if (!res$ok) {
        log_msg(design$id, "rep", k, "ERROR:", res$err)
        next
      }

      dat <- res$dat; th <- res$th; lin <- res$lin; lintest <- res$lintest

      for (i in seq_along(c_rho_grid)) {
        a <- res$anchor_results[i, ]
        new_coef_rows[[length(new_coef_rows)+1]] <- extract_fit_rows(
          design, k, a$model, res$th_all[[i]], dat, a$c_rho, a$rho_n)
      }
      new_coef_rows[[length(new_coef_rows)+1]] <- extract_fit_rows(design, k, "linear", lin, dat)
      scalar <- res$anchor_results[match(c_rho_primary, c_rho_grid), ]
      scalar$model <- NULL
      scalar$elapsed <- res$elapsed
      scalar$tau_0_dgp <- design$tau_0
      scalar$method_version <- METHOD_VERSION
      scalar$stat_settings <- I(list(run_settings))
      old_columns <- c("design_id", "rep", "n", "p", "tau_0", "elapsed", "th_best_tau",
                       "lintest_Tn", "lintest_pval", "lintest_reject")
      scalar <- scalar[c(old_columns, setdiff(names(scalar), old_columns))]
      scalar$anchor_results <- I(list(res$anchor_results))
      scalar$profile <- I(list(res$profile))
      scalar$lambda_search <- I(list(res$lambda_search))
      scalar$test_details <- I(list(list(Tn = lintest$Tn, pval = lintest$pval,
        reject = lintest$reject, taugrid = lintest$taugrid,
        lambda_glmnet = scalar$lambda_glmnet,
        lambda_tuning = "selected_common_lambda", lambda0 = "gic")))
      new_scalar_rows[[length(new_scalar_rows)+1]] <- scalar

      log_msg(design$id, "rep", k, "/", B_mc, "done in", round(res$elapsed, 1), "sec  |",
              "th_best_tau:", th$best_tau, "| c_rho:", scalar$c_rho, "rho_n:", scalar$rho_n,
              "| common lambda_glmnet:", scalar$lambda_glmnet,
              "| lambda grid boundary:", res$lambda_search$selected_at_boundary,
              "| lintest reject:", lintest$reject,
              "pval:", round(lintest$pval, 4))

      n_since_checkpoint <- n_since_checkpoint + 1
    }

    # Can only checkpoint at a batch boundary (mclapply blocks until the
    # WHOLE batch returns) -- a crash loses at most one in-flight batch
    # (<= OUTER_NCORES reps) on top of whatever's short of CHECKPOINT_EVERY.
    is_last_batch <- (batch_i == length(batches))
    if (n_since_checkpoint >= CHECKPOINT_EVERY || is_last_batch) {
      n_total <- merge_and_save(design$id, new_coef_rows, new_scalar_rows)
      log_msg(design$id, "checkpoint saved --", n_total, "/", B_mc, "reps total on disk")
      new_coef_rows      <- list()
      new_scalar_rows    <- list()
      n_since_checkpoint <- 0
    }
  }
}

log_msg("All designs complete.")
