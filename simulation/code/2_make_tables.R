# ============================================================
# Title: Step 2 -- results/ -> output/: the simulation tables of the paper (Tables tbl:confidence,
#        tbl:FWER)
# Author: Hongqiang Yan
# Date: 2026-11-24
# Input:  results/<design>_coef.rds, <design>_scalar.rds
#         (written by code/1_simulate.R)
# Output: output/tbl_confidence.tex, output/tbl_FWER.tex
# Run from the simulation/ folder:  Rscript code/2_make_tables.R
# No simulation is run here.
# ============================================================

# 0. Setup ----
result_dir <- "results/"
table_dir <- "output/"
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# Rows of both tables, in the paper's order.
designs <- c(
  "Baseline" = "ddd1",
  "$p=300$" = "ddd8",
  "$b_2=0.5$" = "ddd2",
  "Dependent $Q$" = "ddd5",
  "$\\tau_0=0.4$" = "ddd6",
  "$b_2=0$" = "ddd3",
  "$n=1000,\\ p=600$" = "ddd9"
)

# 1. Metrics per design ----
# Confidence intervals: all 2p threshold-model coefficients alpha = (beta, delta),
# reported (widened) intervals low_lim/up_lim, grouped by whether alpha_0 is zero.
# Tests: per replication, FP/TP are rejection rates of individual 5% Wald tests among
# zero/nonzero coefficients; FWER is any Bonferroni rejection among zero coefficients.
design_metrics <- function(id) {
  co <- readRDS(file.path(result_dir, paste0(id, "_coef.rds")))
  sc <- readRDS(file.path(result_dir, paste0(id, "_scalar.rds")))
  stopifnot(nrow(sc) == 500L, !anyDuplicated(sc$rep))
  th <- co[co$model == "threshold" & co$param_type %in% c("beta", "delta"), ]
  len <- th$up_lim - th$low_lim
  covered <- th$low_lim <= th$true_value & th$up_lim >= th$true_value
  nonzero <- th$true_value != 0
  test_rates <- function(d) {
    by_rep <- vapply(split(d, d$rep), function(a) {
      zero <- a$structural_truth == 0
      c(FP = mean(a$pval[zero] <= 0.05), TP = mean(a$pval[!zero] <= 0.05),
        FWER = as.numeric(any(a$pval[zero] <= 0.05 / nrow(a))))
    }, numeric(3))
    rowMeans(by_rep)
  }
  delta_true <- th$true_value[th$param_type == "delta"]
  no_effect <- all(delta_true == 0)
  data.frame(
    design = id, reps = nrow(sc), no_effect = no_effect,
    mae = if (no_effect) NA_real_ else mean(abs(sc$th_best_tau - sc$tau_0)),
    tau_exact = if (no_effect) NA_real_ else mean(abs(sc$th_best_tau - sc$tau_0) < 1e-12),
    anchor = mean(abs(sc$th_best_tau - sc$tau_star) < 1e-12),
    anchor_count = sum(abs(sc$th_best_tau - sc$tau_star) < 1e-12),
    len_all = mean(len), len_nonzero = mean(len[nonzero]), len_zero = mean(len[!nonzero]),
    cov_all = mean(covered), cov_nonzero = mean(covered[nonzero]),
    cov_zero = mean(covered[!nonzero]),
    t(setNames(test_rates(th), paste0("th_", c("FP", "TP", "FWER")))),
    t(setNames(test_rates(co[co$model == "linear" & co$param_type == "beta", ]),
      paste0("lin_", c("FP", "TP", "FWER")))),
    reject = mean(sc$lintest_reject)
  )
}
m <- do.call(rbind, lapply(designs, design_metrics))
m$label <- names(designs)
f3 <- function(x) sprintf("%.3f", x)

# 2. Table tbl:confidence ----
confidence_rows <- sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", m$label,
  ifelse(is.na(m$mae), "n.a.", sprintf("%.5f", m$mae)), f3(m$anchor),
  f3(m$len_all), f3(m$len_nonzero), f3(m$len_zero),
  f3(m$cov_all), f3(m$cov_nonzero), f3(m$cov_zero))
writeLines(c(
  r"(\begin{table*}[t])",
  r"(\caption{Threshold Estimation and Confidence Intervals for the Slope Coefficients})",
  r"(\label{tbl:confidence})",
  r"(\tabcolsep=0pt)",
  r"(\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrrrrrrr@{\extracolsep{\fill}}})",
  r"(\toprule)",
  r"(& \multicolumn{2}{c}{Threshold estimation})",
  r"(& \multicolumn{3}{c}{Average length})",
  r"(& \multicolumn{3}{c}{Coverage} \\)",
  r"(\cmidrule(lr){2-3}\cmidrule(lr){4-6}\cmidrule(lr){7-9})",
  r"(Changed parameters & MAE & \shortstack{Anchor\\selected})",
  r"(& All & Nonzero & Zero & All & Nonzero & Zero \\)",
  r"(\midrule)",
  confidence_rows,
  r"(\botrule)",
  r"(\end{tabular*})",
  "",
  r"(\begin{tablenotes})",
  r"(\item \textit{Notes.} The anchor is $\tau_\ast=0.3$. When $b_2=0$,)",
  r"($\tau_0$ is unidentified and all $p$ entries of $\bm{\delta}_0$ belong)",
  r"(to the zero group.)",
  r"(\end{tablenotes})",
  r"(\end{table*})"
), file.path(table_dir, "tbl_confidence.tex"))

# 3. Table tbl:FWER ----
fwer_rows <- unlist(lapply(seq_len(nrow(m)), function(i) c(
  m$label[i],
  paste0("& ", f3(m$th_FP[i]), " & ", f3(m$th_TP[i]), " & ", f3(m$th_FWER[i])),
  paste0("& ", f3(m$lin_FP[i]), " & ", f3(m$lin_TP[i]), " & ", f3(m$lin_FWER[i])),
  paste0("& ", f3(m$reject[i]), " \\\\")
)))
writeLines(c(
  r"(\begin{table*}[t])",
  r"(\caption{Coefficient Tests and the Threshold-Effect Test})",
  r"(\label{tbl:FWER})",
  r"(\tabcolsep=0pt)",
  r"(\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrrrrrr@{\extracolsep{\fill}}})",
  r"(\toprule)",
  r"(& \multicolumn{3}{c}{Debiased threshold Lasso})",
  r"(& \multicolumn{3}{c}{Linear benchmark})",
  r"(& \shortstack{Threshold-effect\\test} \\)",
  r"(\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-8})",
  r"(Changed parameters)",
  r"(& FP & TP & FWER)",
  r"(& FP & TP & FWER)",
  r"(& Rej. \\)",
  r"(\midrule)",
  fwer_rows,
  r"(\botrule)",
  r"(\end{tabular*})",
  "",
  r"(\begin{tablenotes})",
  r"(\item \textit{Notes.} The designs are the same as in Table~\ref{tbl:confidence}.)",
  r"(The Bonferroni cutoffs are $0.05/(2p)$ for the threshold model)",
  r"(and $0.05/p$ for the linear benchmark. Rej.\ is the rejection)",
  r"(rate of the threshold-effect test at the 5\% level.)",
  r"(\end{tablenotes})",
  r"(\end{table*})"
), file.path(table_dir, "tbl_FWER.tex"))

cat(format(Sys.time()), "- wrote tbl_confidence.tex and tbl_FWER.tex to", table_dir, "\n")
