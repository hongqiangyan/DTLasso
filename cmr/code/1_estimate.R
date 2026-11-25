# ============================================================
# Step 1 -- data/ -> results/
# Estimates all eligible metropolitan area-decade samples and collects the
# results in results/all_city_summary_audited.csv and
# results/all_city_coefficients_audited.csv.
# Run from the cmr/ folder:  Rscript code/1_estimate.R
# ============================================================

entry <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", entry))))
setwd(root)

# Settings used in the paper. A value already set in the environment is kept,
# e.g. TIPPING_MSAS=3980 TIPPING_DECADES=1990-2000 Rscript code/1_estimate.R
settings <- c(
  TIPPING_DECADES = "1990-2000,2000-2010", TIPPING_MSAS = "",
  TIPPING_N_FOLDS = "10", TIPPING_NLAMBDA = "100", TIPPING_BOOTSTRAP_B = "499",
  TIPPING_C_RHO = "1", TIPPING_N_CORES = "4", TIPPING_RESUME = "true",
  TIPPING_OUTPUT_DIR = file.path(root, "results"),
  TIPPING_CHECK_DIR = file.path(root, "results"),
  OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)
for (name in names(settings)) {
  if (!nzchar(Sys.getenv(name))) do.call(Sys.setenv, as.list(settings[name]))
}
dir.create(file.path(root, "results"), showWarnings = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script) {
  cat("Running", script, "\n")
  status <- system2(rscript, shQuote(file.path(root, "code", script)))
  if (status != 0) stop(script, " failed.")
}
run("preflight_card_cities_gt.R")           # eligible samples -> results/sample_audit.csv
run("tipping_DTLasso_percity_anchor_gt.R")  # one result per sample -> results/rds/
run("summarize_results_gt.R")               # -> the two audited CSVs
