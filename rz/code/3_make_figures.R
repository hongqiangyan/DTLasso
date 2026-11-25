# Step 3 -- results/ -> output/: the two IRF figures in the supplement.
# Run from the rz/ folder:  Rscript code/3_make_figures.R
# Only legends, numerical ticks and panel letters are drawn in R.
entry <- grep("^--file=", commandArgs(FALSE), value=TRUE)
here <- if(length(entry)) dirname(normalizePath(sub("^--file=","",entry[1]))) else normalizePath(getwd())
ext <- here
results <- file.path(dirname(ext),"results/lp")
out <- file.path(dirname(here),"output"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
library(ggplot2); library(patchwork)
a <- readRDS(file.path(dirname(ext),"results/threshold/news/selection.rds"))
irf <- read.csv(file.path(results,"common_sample_irf.csv"))
diff <- read.csv(file.path(results,"irf_differences.csv"))
counts <- read.csv(file.path(results,"sample_counts.csv"))
stopifnot(a$result$converged,a$key$c_rho==.75,
 abs(100*a$result$selected$tau-9.273)<1e-10,all(irf$n==480),
 !anyDuplicated(irf[c("case","outcome","horizon","regime")]))
style <- theme_classic(base_size=10,base_family="serif") +
 theme(axis.title=element_blank(),axis.text=element_text(size=9,colour="grey15"),
 axis.line=element_line(linewidth=.3,colour="grey40"),
 axis.ticks=element_line(linewidth=.3,colour="grey40"),
 panel.grid.major.y=element_line(colour="grey93",linewidth=.2),
 plot.title=element_text(size=11,face="plain",margin=margin(b=7)),
 legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=9),
 plot.margin=margin(7,9,7,7))
palette <- c("High unemployment"="#2463A6","Low unemployment"="#C3423F")
line_types <- c("High unemployment"="solid","Low unemployment"="dashed")
save_plot <- function(name,panels) {
 p <- wrap_plots(panels,ncol=2,guides="collect") & theme(legend.position="bottom")
 # The supplement includes the PNG versions only.
 ggsave(file.path(out,paste0(name,".png")),p,device=ragg::agg_png,
  width=7.2,height=5,dpi=240,bg="white")
}
make_irf <- function(cases,name) {
 panels <- list()
 for(v in c("g","y")) {
  row <- irf[irf$case %in% cases & irf$outcome==v,]
  ylim <- range(0,row$lower,row$upper)
  for(cc in cases) {
   d <- row[row$case==cc,]; d$regime <- factor(d$regime,levels=names(palette))
   stopifnot(nrow(d)==42,all(is.finite(as.matrix(d[c("estimate","lower","upper")]))),
    all(d$lower<=d$estimate),all(d$estimate<=d$upper))
   panels[[length(panels)+1L]] <- ggplot(d,aes(horizon,estimate,colour=regime,fill=regime,linetype=regime)) +
    geom_ribbon(aes(ymin=lower,ymax=upper),colour=NA,linetype=0,alpha=.15,show.legend=FALSE) +
    geom_hline(yintercept=0,colour="grey50",linewidth=.3,linetype="dotted") +
    geom_line(linewidth=.7) + scale_colour_manual(values=palette,limits=names(palette)) +
    scale_fill_manual(values=palette,limits=names(palette)) + scale_linetype_manual(values=line_types) +
    scale_x_continuous(breaks=seq(0,20,5)) + scale_y_continuous(breaks=scales::breaks_pretty(n=5)) +
    coord_cartesian(ylim=ylim) + labs(title=paste0("(",letters[length(panels)+1L],")")) + style
  }
 }
 save_plot(name,panels)
}
make_irf(c("rz_fixed","news_anchor"),"fiscal_irfs")
panels <- list()
for(v in c("g","y")) {
 row <- diff[diff$case %in% c("rz_fixed","news_anchor") & diff$outcome==v,]
 ylim <- range(0,row$lower,row$upper)
 for(cc in c("rz_fixed","news_anchor")) {
  d <- row[row$case==cc,]
  panels[[length(panels)+1L]] <- ggplot(d,aes(horizon,estimate)) +
   geom_ribbon(aes(ymin=lower,ymax=upper),fill="#2463A6",alpha=.15) +
   geom_hline(yintercept=0,colour="grey50",linewidth=.3,linetype="dotted") +
   geom_line(aes(colour="High minus low"),linewidth=.7) +
   scale_colour_manual(values=c("High minus low"="#2463A6")) +
   scale_x_continuous(breaks=seq(0,20,5)) + scale_y_continuous(breaks=scales::breaks_pretty(n=5)) +
   coord_cartesian(ylim=ylim) + labs(title=paste0("(",letters[length(panels)+1L],")")) + style
 }
}
save_plot("fiscal_irf_differences",panels)
cat("Saved output/fiscal_irfs.png and output/fiscal_irf_differences.png.\n")
