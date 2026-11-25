# ============================================================
# Title: Step 2 -- results/ -> output/: the tipping tables in the paper
# Author: Hongqiang Yan
# Date: 2026-11-24
# Input:  results/all_city_summary_audited.csv
#         results/all_city_coefficients_audited.csv
#         (both written by code/summarize_results_gt.R)
# Output: output/
#           main paper:  tb_tipping-tests.tex, tb_indianapolis-coefficients.tex
#           supplement:  tb_citycoef-<msa>-<period>.tex for Indianapolis (3480),
#                        Nashville (5360) and Chicago (1600)
# No model is fitted here. Run from any directory:
#   Rscript code/2_make_tables.R
# ============================================================

# 0. Setup ----
args <- commandArgs(FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", args,
  value = TRUE)[1]))))
result_dir <- file.path(root, "results")
table_dir <- file.path(root, "output")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
old <- list.files(table_dir, pattern = "[.](tex|csv)$", full.names = TRUE)
invisible(file.remove(old))  # outputs are fully regenerated on every run

summary_df <- read.csv(file.path(result_dir, "all_city_summary_audited.csv"),
  stringsAsFactors = FALSE)
coef_df <- read.csv(file.path(result_dir, "all_city_coefficients_audited.csv"),
  stringsAsFactors = FALSE)
decades <- c("1990-2000", "2000-2010")
by_decade <- lapply(setNames(decades, decades), function(x) summary_df[summary_df$decade == x, ])

# CMR's six original controls (core_controls in tipping_DTLasso_percity_anchor_gt.R).
cmr_controls <- c("unemployment_rate", "log_mean_family_income", "public_transit_share",
  "vacancy_share", "renter_share", "single_unit_share")

# Cities with full coefficient tables in the supplement, in the supplement's order.
supplement_cities <- c("3480" = "Indianapolis", "5360" = "Nashville", "1600" = "Chicago")

# Row order and labels of the supplement coefficient tables.
control_labels <- c(
  poverty_high_175_plus_share = "Income at least 175\\% of the poverty line",
  poverty_low_under_100_share = "Below poverty line",
  commute_bicycle_share = "Bicycle commute",
  heating_bottled_gas_share = "Bottled-gas heat",
  commute_car_share = "Car commute",
  heating_coal_coke_share = "Coal/coke heat",
  travel_high_60_plus_share = "Commute 60+ min",
  travel_low_under_30_share = "Commute under 30 min",
  departure_late_0900_plus_share = "Depart 9 a.m. or later",
  departure_early_before_0700_share = "Depart before 7 a.m.",
  earnings_receipt_share = "Earnings receipt",
  heating_electricity_share = "Electric heat",
  english_only_share = "English-only persons age 5+",
  vehicles_high_4_plus_share = "Four+ vehicles",
  heating_fuel_oil_share = "Fuel-oil heat",
  structure_high_density_share = "High-density housing",
  householder_age_low_15_34_share = "Householder age 15--34",
  householder_age_high_65_plus_share = "Householder age 65+",
  housing_age_high_41_plus_years_share = "Housing built 41+ yrs ago",
  housing_age_low_0_20_years_share = "Housing built prev. 20 yrs",
  intercept = "Intercept",
  lacking_kitchen_share = "Lacking kitchen",
  lacking_plumbing_share = "Lacking plumbing",
  log_mean_family_income = "Log mean fam. income",
  log_tract_land_area_sqkm = "Log tract land area",
  structure_low_density_share = "Low-density housing",
  log1p_mean_public_assistance_income_among_recipients_1999usd =
    "Log mean public-assist. income",
  log1p_mean_retirement_income_among_recipients_1999usd = "Log mean retirement income",
  log1p_mean_social_security_income_among_recipients_1999usd =
    "Log mean Social Security income",
  log1p_mean_wage_income_among_recipients_1999usd = "Log mean wage income",
  structure_middle_density_share = "Middle-density housing",
  commute_motorcycle_share = "Motorcycle commute",
  move_duration_high_21_plus_years_share = "Moved in 21+ yrs ago",
  move_duration_low_0_5_years_share = "Moved in prev. 5 yrs",
  heating_no_fuel_share = "No-fuel heat",
  pays_extra_utilities_share = "Utilities paid separately from rent",
  public_assistance_receipt_share = "Public assist. receipt",
  public_transit_share = "Public transit share",
  renter_share = "Renter-occupied share",
  retirement_receipt_share = "Retirement receipt",
  single_unit_share = "Single-unit share",
  social_security_receipt_share = "Social Security receipt",
  heating_solar_share = "Solar heat",
  unemployment_rate = "Unemployment rate",
  heating_utility_gas_share = "Utility-gas heat",
  vacancy_share = "Vacancy share",
  wage_receipt_share = "Wage receipt",
  commute_walked_share = "Walked to work",
  heating_wood_share = "Wood heat",
  commute_worked_home_share = "Worked at home",
  vehicles_low_0_1_share = "Zero/one vehicle"
)
stopifnot(length(control_labels) == 51L)

# Rows of the Indianapolis table in the main text: variable, label, group.
indianapolis_rows <- data.frame(
  variable = c("single_unit_share", "move_duration_low_0_5_years_share",
    "housing_age_high_41_plus_years_share", "pays_extra_utilities_share",
    "departure_late_0900_plus_share"),
  label = c("Single-unit housing share", "Recent-mover share",
    "Share of housing at least 41 years old",
    "Separate utility payments (share of rental units)",
    "Share of commuters leaving at or after 9 a.m."),
  stringsAsFactors = FALSE)

# 1. Helpers ----
period_label <- function(decade) gsub("-", "--", decade, fixed = TRUE)
write_table <- function(lines, name) {
  writeLines(lines, file.path(table_dir, paste0(name, ".tex")), useBytes = TRUE)
}
city_row <- function(msa, decade) {
  row <- summary_df[summary_df$msa == as.integer(msa) & summary_df$decade == decade, ]
  stopifnot(nrow(row) == 1L)
  row
}
city_coefs <- function(msa, decade) {
  coef_df[coef_df$msa == as.integer(msa) & coef_df$decade == decade, ]
}
get_coef <- function(cc, variable) {
  r <- cc[cc$variable == variable, ]
  stopifnot(nrow(r) == 1L)
  r
}
rejects <- function(d) d$lintest_pval <= 0.05
count_pct <- function(k, n) sprintf("%d (%.1f\\%%)", k, 100 * k / n)
stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) "^{***}" else if (p < 0.05) "^{**}" else if (p < 0.10) "^{*}" else ""
}
# Two significant digits; scientific values as m\!\times\!10^{e}.
fmt2 <- function(x) {
  if (x == 0) return("0")
  s <- sprintf("%.2g", x)
  if (!grepl("e", s, fixed = TRUE)) return(s)
  paste0(sub("e.*$", "", s), "\\!\\times\\!10^{", as.integer(sub("^.*e", "", s)), "}")
}

# Slope differences significant after Bonferroni (5%), split by control group.
jumps <- coef_df[grepl("^threshold_", coef_df$variable) &
  coef_df$variable != "threshold_intercept", ]
jumps$control <- sub("^threshold_", "", jumps$variable)
jumps$significant <- !is.na(jumps$p_bonferroni_city) & jumps$p_bonferroni_city <= 0.05
cities_with_jump <- function(decade, cmr) {
  z <- jumps[jumps$decade == decade & jumps$significant &
    (jumps$control %in% cmr_controls) == cmr, ]
  length(unique(z$msa))
}

# 2. Main paper: threshold effects across metropolitan areas (tb:tipping-tests) ----
n_msa <- vapply(by_decade, nrow, 1L)
write_table(c(
  r"(\begin{table*}[t])",
  r"(\caption{Threshold Effects Across Metropolitan Areas})",
  r"(\label{tb:tipping-tests})",
  r"(\tabcolsep=0pt)",
  r"(\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrr@{\extracolsep{\fill}}})",
  r"(\toprule)",
  r"(& 1990--2000 & 2000--2010 \\)",
  r"(\midrule)",
  sprintf("Metropolitan areas & %d & %d \\\\", n_msa[[1]], n_msa[[2]]),
  "Threshold-effect test rejections",
  sprintf("& %s & %s \\\\", count_pct(sum(rejects(by_decade[[1]])), n_msa[[1]]),
    count_pct(sum(rejects(by_decade[[2]])), n_msa[[2]])),
  r"(\addlinespace)",
  r"(\multicolumn{3}{l}{\textit{At least one significant slope difference}} \\)",
  sprintf("CMR's original controls & %d & %d \\\\", cities_with_jump(decades[1], TRUE),
    cities_with_jump(decades[2], TRUE)),
  sprintf("Additional controls & %d & %d \\\\", cities_with_jump(decades[1], FALSE),
    cities_with_jump(decades[2], FALSE)),
  r"(\botrule)",
  r"(\end{tabular*})",
  r"(\begin{tablenotes})",
  r"(\item Notes: Entries count metropolitan areas. The threshold-effect)",
  r"(test uses a 5\% level. The final two rows count areas with)",
  r"(at least one significant $\delta_j$ at 5\% after Bonferroni)",
  r"(correction across all coefficient tests within the)",
  r"(metropolitan area and decade.)",
  r"(\end{tablenotes})",
  r"(\end{table*})"
), "tb_tipping-tests")

# 3. Main paper: selected coefficients for Indianapolis ----
indy <- lapply(setNames(decades, decades), function(x) city_coefs("3480", x))
n_tests <- vapply(indy, function(cc) sum(is.finite(cc$p_value)), 1L)
stopifnot(n_tests[[1]] == n_tests[[2]])
indy_block <- function(variable, label) {
  est <- character(); se <- character()
  for (x in decades) {
    b <- get_coef(indy[[x]], variable); d <- get_coef(indy[[x]], paste0("threshold_", variable))
    est <- c(est, sprintf("$%.3f%s$", b$coefficient, stars(b$p_bonferroni_city)),
      sprintf("$%.3f%s$", d$coefficient, stars(d$p_bonferroni_city)))
    se <- c(se, sprintf("$(%.3f)$", b$sd), sprintf("$(%.3f)$", d$sd))
  }
  c(label, paste0("& ", paste(est, collapse = " & "), " \\\\"),
    paste0("& ", paste(se, collapse = " & "), " \\\\"), r"(\addlinespace)")
}
indy_summary <- lapply(decades, function(x) city_row("3480", x))
two_col <- function(a, b) c(sprintf("& \\multicolumn{2}{r}{%s}", a),
  sprintf("& \\multicolumn{2}{r}{%s} \\\\", b))
write_table(c(
  r"(\begin{table*}[t])",
  r"(\caption{Selected Coefficient Estimates for Indianapolis})",
  r"(\label{tb:indianapolis-coefficients})",
  r"(\tabcolsep=0pt)",
  r"(\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrrr@{\extracolsep{\fill}}})",
  r"(\toprule)",
  r"(& \multicolumn{2}{c}{1990--2000})",
  r"(& \multicolumn{2}{c}{2000--2010} \\)",
  r"(\cmidrule(lr){2-3}\cmidrule(lr){4-5})",
  r"(& $\widehat\beta_j$ & $\widehat\delta_j$)",
  r"(& $\widehat\beta_j$ & $\widehat\delta_j$ \\)",
  r"(\midrule)",
  r"(\multicolumn{5}{l}{\textit{CMR's original control}} \\)",
  indy_block(indianapolis_rows$variable[1], indianapolis_rows$label[1]),
  r"(\multicolumn{5}{l}{\textit{Additional controls}} \\)",
  unlist(lapply(2:nrow(indianapolis_rows), function(i)
    indy_block(indianapolis_rows$variable[i], indianapolis_rows$label[i]))),
  r"(\addlinespace)",
  "Tracts",
  two_col(indy_summary[[1]]$n_tracts, indy_summary[[2]]$n_tracts),
  r"(Tipping point (\%))",
  two_col(sprintf("%.1f", 100 * indy_summary[[1]]$tau_hat),
    sprintf("%.1f", 100 * indy_summary[[2]]$tau_hat)),
  r"(Threshold-effect test $p$-value)",
  two_col(sprintf("%.3f", indy_summary[[1]]$lintest_pval),
    sprintf("%.3f", indy_summary[[2]]$lintest_pval)),
  r"(\botrule)",
  r"(\end{tabular*})",
  r"(\begin{tablenotes})",
  r"(\item Notes: Debiased threshold Lasso estimates with)",
  r"(heteroskedasticity-robust standard errors in parentheses.)",
  r"($\beta_j$ is the slope for $Q_i\leq\widehat\tau$, and)",
  r"($\delta_j$ is its change for $Q_i>\widehat\tau$.)",
  r"(Stars indicate significance based on two-sided Wald $p$-values with a Bonferroni)",
  sprintf("correction across all %d coefficient tests in each decade", n_tests[[1]]),
  r"(($^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$).)",
  r"(\end{tablenotes})",
  r"(\end{table*})"
), "tb_indianapolis-coefficients")

# 4. Supplement: all candidate controls for the three cities ----
se_note <- function(x) paste0("{\\fontsize{7}{8}\\selectfont ", x, "}")
coef_row <- function(cc, variable) {
  base <- variable
  jump <- if (variable == "intercept") "threshold_intercept" else paste0("threshold_", variable)
  b <- cc[cc$variable == base, ]; d <- cc[cc$variable == jump, ]
  label <- control_labels[[variable]]
  if (!nrow(b) && !nrow(d)) return(paste0(label, " & -- & -- & -- & -- \\\\"))
  stopifnot(nrow(b) == 1L, nrow(d) == 1L)
  debiased <- function(r) {
    se <- if (is.na(r$sd)) "(n/a)" else paste0("($", fmt2(r$sd), "$)")
    paste0("$", fmt2(r$coefficient), stars(r$p_value), "$ ", se_note(se))
  }
  paste0(label, " & $", fmt2(b$lasso_coefficient), "$ & $", fmt2(d$lasso_coefficient),
    "$ & ", debiased(b), " & ", debiased(d), " \\\\")
}
panel <- function(rows, name) c(
  r"(\begin{minipage}[t]{0.49\textwidth})",
  r"(\centering)",
  paste0(r"(\textit{Panel )", name, r"(}\par\vspace{0pt})"),
  r"(\begin{tabular}{@{}>{\raggedright\arraybackslash}p{0.33\linewidth}rrrr@{}})",
  r"(\toprule)",
  r"(Control & $\beta_L$ & $\delta_L$ & $\beta_D$ & $\delta_D$ \\)",
  r"(\midrule)",
  rows,
  r"(\bottomrule)",
  r"(\end{tabular})"
)
for (msa in names(supplement_cities)) for (x in decades) {
  cc <- city_coefs(msa, x); s <- city_row(msa, x)
  rows <- vapply(names(control_labels), function(v) coef_row(cc, v), "")
  split_at <- ceiling(length(rows) / 2)
  label <- paste0(r"(\label{tb:citycoef-)", msa, "-", x, "}",
    if (x == decades[1]) paste0(r"(\label{tb:citycoef-)", msa, "}") else "")
  write_table(c(
    r"(\noindent\begin{minipage}{\linewidth})",
    r"(\centering)",
    r"(\setlength{\tabcolsep}{1.8pt})",
    r"(\renewcommand{\arraystretch}{1.02})",
    paste0(r"(\captionof{table}{All Candidate Controls: Lasso and Debiased Estimates for )",
      supplement_cities[[msa]], ", ", period_label(x), "}"),
    label,
    r"({\fontsize{8}{9.5}\selectfont)",
    panel(rows[seq_len(split_at)], "A"),
    r"(\end{minipage}\hfill)",
    panel(rows[-seq_len(split_at)], "B"),
    r"(\end{minipage})",
    "}",
    r"(\begin{flushleft})",
    r"(\fontsize{8}{9.5}\selectfont)",
    r"(Notes: Subscripts $L$ and $D$ denote initial Lasso and)",
    r"(debiased threshold Lasso estimates, respectively.)",
    r"(Heteroskedasticity-robust standard errors are in parentheses.)",
    r"(The Intercept row reports $\mu$ in the $\beta$ columns)",
    r"(and $\gamma$ in the $\delta$ columns.)",
    r"(Dashes indicate excluded controls, and n/a indicates)",
    r"(that no standard error is reported.)",
    r"(Stars indicate significance based on unadjusted two-sided Wald $p$-values, with)",
    r"($^{*}p<0.10$, $^{**}p<0.05$, and $^{***}p<0.01$.)",
    sprintf(r"($N=%d$, $\widehat\tau=%.4f$, $\tau_\ast=%.4f$.)", s$n_tracts, s$tau_hat,
      s$tau_star),
    sprintf(r"(Threshold-effect test: $T_n=%.2f$, bootstrap $p=%.3f$.)", s$lintest_Tn,
      s$lintest_pval),
    r"(\end{flushleft})",
    r"(\end{minipage})"
  ), paste0("tb_citycoef-", msa, "-", x))
}

cat(format(Sys.time()), "- wrote", length(list.files(table_dir, pattern = "[.]tex$")),
  "tables to", table_dir, "\n")
