# ============================================================
# Title: Step 4 -- results/ -> output/fiscal_multipliers.tex
#        (supplement Table tab:fiscal-multipliers)
# Author: Hongqiang Yan
# Input:  results/lp/multipliers.csv, results/lp/differences.csv,
#         results/threshold/selection_summary.csv
# Output: output/fiscal_multipliers.tex
# Run from the rz/ folder:  Rscript code/4_make_tables.R
# ============================================================

entry <- grep("^--file=", commandArgs(FALSE), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", entry))))
read <- function(f) read.csv(file.path(base, "results", f), stringsAsFactors = FALSE)
multipliers <- read("lp/multipliers.csv")
differences <- read("lp/differences.csv")
threshold <- read("threshold/selection_summary.csv")
out <- file.path(base, "output")
dir.create(out, showWarnings = FALSE)

fmt <- function(x) formatC(x, format = "f", digits = 3)
cell <- function(r) paste0("\\shortstack{", fmt(r$estimate), "\\\\$[", fmt(r$lower), ",",
  fmt(r$upper), "]$}")

rows <- character()
for (case in c("news_anchor", "rz_fixed")) for (h in c(7L, 15L)) {
  z <- multipliers[multipliers$case == case & multipliers$horizon == h, ]
  high <- z[z$regime == "High unemployment", ]
  low <- z[z$regime == "Low unemployment", ]
  diff <- differences[differences$case == case & differences$horizon == h, ]
  stopifnot(nrow(high) == 1L, nrow(low) == 1L, nrow(diff) == 1L)
  label <- if (case == "news_anchor") {
    paste0("Estimated $", fmt(threshold$tau_percent), "\\%$")
  } else {
    "Fixed $6.5\\%$"
  }
  rows <- c(rows,
    paste0(label, " & ", (h + 1) / 4, " years"),
    paste0("& ", cell(high)),
    paste0("& ", cell(low)),
    paste0("& ", cell(diff), " \\\\[4pt]"))
}

writeLines(c(
  r"(\begin{table}[!tbp])",
  r"(\centering)",
  r"(\caption{Cumulative fiscal multipliers under alternative)",
  r"(unemployment thresholds.})",
  r"(\label{tab:fiscal-multipliers})",
  r"(\small)",
  r"(\setlength{\tabcolsep}{4pt})",
  r"(\begin{tabular}{@{}llccc@{}})",
  r"(\toprule)",
  r"(Specification & Window & High unemployment)",
  r"(& Low unemployment & High minus low \\)",
  r"(\midrule)",
  rows,
  r"(\bottomrule)",
  r"(\end{tabular})",
  r"(\par\smallskip)",
  r"(\begin{minipage}{0.97\linewidth})",
  r"(\footnotesize)",
  r"(\textit{Notes:} Brackets report 95\% delta-method confidence)",
  r"(intervals. All specifications use the same 480 shock origins.)",
  r"(\end{minipage})",
  r"(\end{table})"
), file.path(out, "fiscal_multipliers.tex"))
cat("Saved output/fiscal_multipliers.tex.\n")
