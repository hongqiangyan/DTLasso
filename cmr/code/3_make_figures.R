# Step 3 -- results/ -> output/: the two maps in the supplement (no model is fitted).
# Run from the cmr/ folder:  Rscript code/3_make_figures.R
# Geographic points: mean of 2000 Census tract internal points within each
# original MSA/PMSA, using the same point for both periods.
library(sf)
library(maps)
library(ragg)

arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=', '', arg[1]))))
here <- file.path(root, 'output')  # PNG outputs
dir.create(here, showWarnings = FALSE)
saved <- read.csv(file.path(root, 'results/all_city_summary_audited.csv'))
# Cities rejecting at 5%, with their map location (data/msa_coordinates.csv).
coords <- read.csv(file.path(root, 'data', 'msa_coordinates.csv'))
dat <- merge(saved[saved$lintest_pval <= .05, c('msa', 'msa_name', 'decade', 'n_tracts',
  'tau_hat', 'tau_star', 'lintest_pval', 'lintest_reject')], coords, by = 'msa')
dat <- dat[order(dat$decade, dat$msa), ]; rownames(dat) <- NULL
dat$lintest_reject <- as.logical(dat$lintest_reject)
key <- function(d) paste(d$msa, d$decade)
target <- saved[saved$lintest_reject,]
stopifnot(nrow(dat)==92L, setequal(key(dat),key(target)),
          all(dat$lintest_reject), all(dat$lintest_pval<=.05),
          all(dat$tau_hat>=0 & dat$tau_hat<=.4),
          all(table(dat$decade)==c(65L,27L)))
ii <- match(key(dat),key(saved))
stopifnot(all(dat$tau_hat==saved$tau_hat[ii]),
          all(dat$lintest_pval==saved$lintest_pval[ii]))

state <- st_as_sf(maps::map('state',plot=FALSE,fill=TRUE))
state <- st_transform(state,5070)
state <- st_make_valid(state)
state <- st_simplify(state,dTolerance=6000,preserveTopology=TRUE)
state <- st_make_valid(state)
points <- st_transform(st_as_sf(dat,coords=c('longitude','latitude'),crs=4269),5070)
xy <- st_coordinates(points)
dat$x <- xy[,1];dat$y <- xy[,2]
# Both periods and the 1990-2000 Northeast inset share one point-area scale.
main_bbox <- c(xmin=-2400000,ymin=240000,xmax=2400000,ymax=3250000)
ne_corners <- st_transform(st_as_sfc(st_bbox(c(xmin=-80.8,ymin=38.6,xmax=-69.5,ymax=44.5),crs=4269)),5070)
ne_bbox <- as.numeric(st_bbox(ne_corners));names(ne_bbox)<-c('xmin','ymin','xmax','ymax')
ne_state <- suppressWarnings(st_crop(state,st_bbox(ne_corners)))
main_rect <- c(.05,.15,11.85,6.75)
ne_rect <- c(12.15,2.50,16.25,6.25)
W <- 16.5;H <- 6.9
# Area, not radius, is proportional to the final anchored threshold.
radius <- function(t) .145*sqrt(t/.4)
tex <- character()
emit <- function(x) tex <<- c(tex,x)
escape <- function(s) {
  ch<-strsplit(s,'',fixed=TRUE)[[1]]
  paste0(vapply(ch,function(x) if(x %in% c('%','&','_','#'))paste0('\\',x) else x,''),collapse='')
}
tikcolor <- function(x) paste0('{rgb,255:red,',strtoi(substr(x,2,3),16L),';green,',
                              strtoi(substr(x,4,5),16L),';blue,',strtoi(substr(x,6,7),16L),'}')
coords <- function(x,y) paste(sprintf('(%.3f,%.3f)',x,y),collapse=' -- ')
transformer <- function(b,r) {
  sc<-min((r[3]-r[1])/(b[3]-b[1]),(r[4]-r[2])/(b[4]-b[2]))
  mid<-c(mean(r[c(1,3)]),mean(r[c(2,4)]))
  bm<-c(mean(b[c(1,3)]),mean(b[c(2,4)]))
  function(x,y) cbind(mid[1]+(x-bm[1])*sc,mid[2]+(y-bm[2])*sc)
}
trmain <- transformer(main_bbox,main_rect);trne <- transformer(ne_bbox,ne_rect)
stroke <- function(x,y,col='#AAB0B3',width=.4) {
  lines(x,y,col=col,lwd=width)
  emit(sprintf('\\draw[draw=%s,line width=%.2fpt] %s;',tikcolor(col),width*.75,coords(x,y)))
}
rectangle <- function(x0,y0,x1,y1,fill='#FFFFFF',border=NA,width=.4) {
  rect(x0,y0,x1,y1,col=fill,border=border,lwd=width)
  draw<-if(is.na(border))'none' else tikcolor(border)
  emit(sprintf('\\path[fill=%s,draw=%s,line width=%.2fpt] (%.3f,%.3f) rectangle (%.3f,%.3f);',tikcolor(fill),draw,width*.75,x0,y0,x1,y1))
}
label <- function(x,y,s,pt=9,adj=.5,col='#202629',bold=FALSE) {
  text(x,y,s,cex=pt/10,adj=c(adj,.5),col=col,font=if(bold)2 else 1)
  anchor<-if(adj==0)'west' else if(adj==1)'east' else 'center'
  emit(sprintf('\\node[anchor=%s,text=%s,inner sep=0pt,font=\\sffamily\\fontsize{%.1f}{%.1f}\\selectfont%s] at (%.3f,%.3f) {%s};',
    anchor,tikcolor(col),pt,pt*1.15,if(bold)'\\bfseries' else '',x,y,escape(s)))
}
dot <- function(x,y,t) {
  rr<-radius(t);cl<-'#226B8E'
  for(i in order(rr,decreasing=TRUE)){
    symbols(x[i],y[i],circles=rr[i],inches=FALSE,add=TRUE,bg=cl,fg='#FFFFFF',lwd=.35)
    emit(sprintf('\\path[fill=%s,draw=white,line width=.22pt] (%.3f,%.3f) circle[radius=%.3fcm];',tikcolor(cl),x[i],y[i],rr[i]))
  }
}
drawstates <- function(obj,tr) {
  # Promote mixed polygon/multipolygon features first, so every component
  # is retained rather than taking only the first island of a multipolygon.
  ps<-st_cast(st_cast(st_geometry(obj),'MULTIPOLYGON'),'POLYGON')
  for(p in ps){
    rr<-lapply(p,function(v)tr(v[,1],v[,2]))
    rr<-Filter(function(v)nrow(v)>=3L && all(is.finite(v)),rr)
    if(!length(rr))next
    xx<-unlist(lapply(rr,function(v)c(v[,1],NA_real_)))
    yy<-unlist(lapply(rr,function(v)c(v[,2],NA_real_)))
    xx<-head(xx,-1L);yy<-head(yy,-1L)
    polypath(xx,yy,col='#F0F2F3',border='#B1B8BC',lwd=.35,rule='evenodd')
    paths<-paste(vapply(rr,function(v)paste0(coords(v[,1],v[,2]),' -- cycle'),''),collapse=' ')
    emit(paste0('\\path[even odd rule,fill=',tikcolor('#F0F2F3'),',draw=',tikcolor('#B1B8BC'),',line width=.24pt] ',paths,';'))
  }
}
# Place every numeric value beside its marker, with short leader lines.
# In 1990-2000, northeastern values are printed once in the enlarged inset.
# This layout changes annotation positions only, never city coordinates.
number_labels <- function(pos,t,bounds,obstacles=pos,obstacle_t=t,pt=7.1) {
  if(!length(t))return(invisible(NULL))
  # Conventional half-up rounding for display only; saved estimates stay exact.
  txt<-sprintf('%.1f',floor(1000*t+.5+1e-9)/10)
  ww<-strwidth(txt,cex=pt/10)+.06
  hh<-rep(strheight('0',cex=pt/10)+.06,length(t))
  radii<-radius(obstacle_t)+.025
  d<-as.matrix(dist(pos));diag(d)<-Inf
  # Give the most crowded locations the first choice of nearby space.
  ord<-order(apply(d,1,min),seq_along(t))
  placed<-matrix(NA_real_,nrow=length(t),ncol=4)
  centres<-matrix(NA_real_,nrow=length(t),ncol=2)
  for(i in ord){
    angles<-seq(0,2*pi,length.out=25)[-25]
    ca<-expand.grid(angle=angles,extra=seq(0,.95,by=.055))
    cs<-cos(ca$angle);sn<-sin(ca$angle)
    reach<-pmin(ifelse(abs(cs)>1e-6,ww[i]/2/abs(cs),Inf),
                ifelse(abs(sn)>1e-6,hh[i]/2/abs(sn),Inf))+radius(t[i])+.045+ca$extra
    cx<-pos[i,1]+cs*reach;cy<-pos[i,2]+sn*reach
    bb<-cbind(cx-ww[i]/2,cy-hh[i]/2,cx+ww[i]/2,cy+hh[i]/2)
    cost<-reach+0.02*abs(sn)
    outside<-bb[,1]<bounds[1]|bb[,2]<bounds[2]|bb[,3]>bounds[3]|bb[,4]>bounds[4]
    cost[outside]<-cost[outside]+10000
    for(j in seq_len(nrow(obstacles))){
      dx<-pmax(bb[,1]-obstacles[j,1],0,obstacles[j,1]-bb[,3])
      dy<-pmax(bb[,2]-obstacles[j,2],0,obstacles[j,2]-bb[,4])
      cost<-cost+1000*pmax(0,radii[j]-sqrt(dx^2+dy^2))
    }
    for(j in which(is.finite(placed[,1]))){
      ox<-pmax(0,pmin(bb[,3],placed[j,3]+.03)-pmax(bb[,1],placed[j,1]-.03))
      oy<-pmax(0,pmin(bb[,4],placed[j,4]+.03)-pmax(bb[,2],placed[j,2]-.03))
      cost<-cost+10000*ox*oy
    }
    k<-which.min(cost);placed[i,]<-bb[k,];centres[i,]<-c(cx[k],cy[k])
  }
  for(i in seq_along(t)){
    stroke(c(pos[i,1],centres[i,1]),c(pos[i,2],centres[i,2]),col='#73828B',width=.45)
  }
  for(i in seq_along(t)){
    b<-placed[i,];rectangle(b[1],b[2],b[3],b[4])
    label(centres[i,1],centres[i,2],txt[i],pt=pt)
  }
  invisible(centres)
}
panels <- list()
for(dec in c('1990-2000','2000-2010')){
  z<-dat[dat$decade==dec,];tex<-character()
  show_ne<-dec=='1990-2000'
  file<-file.path(here,paste0('tipping_map_',gsub('-','_',dec),'.png'))
  ragg::agg_png(file,width=W,height=H,units='cm',res=360,background='white',pointsize=10)
  par(mar=c(0,0,0,0),xaxs='i',yaxs='i',family='sans',xpd=NA)
  plot.new();plot.window(c(0,W),c(0,H),asp=1)
  emit(sprintf('\\begin{tikzpicture}[x=1cm,y=1cm]\n\\path[use as bounding box] (0,0) rectangle (%.2f,%.2f);',W,H))
  drawstates(state,trmain)
  zz<-trmain(z$x,z$y)
  stopifnot(all(zz[,1]>=main_rect[1] & zz[,1]<=main_rect[3]),
            all(zz[,2]>=main_rect[2] & zz[,2]<=main_rect[4]))
  in_ne<-show_ne & z$x>=ne_bbox[1]&z$x<=ne_bbox[3]&z$y>=ne_bbox[2]&z$y<=ne_bbox[4]
  number_labels(zz[!in_ne,,drop=FALSE],z$tau_hat[!in_ne],main_rect,
                obstacles=zz,obstacle_t=z$tau_hat)
  dot(zz[,1],zz[,2],z$tau_hat)
  if(show_ne){
    rectangle(ne_rect[1],ne_rect[2],ne_rect[3],ne_rect[4],border='#ABB2B6')
    drawstates(ne_state,trne)
    zn<-z[in_ne,]
    xyne<-trne(zn$x,zn$y)
    number_labels(xyne,zn$tau_hat,ne_rect+c(.025,.025,-.025,-.025))
    dot(xyne[,1],xyne[,2],zn$tau_hat)
  }
  # Only the data labels and a compact size legend are printed in the figure.
  vv<-c(5,10,20,30,40);xx<-seq(12.55,15.85,length.out=length(vv))
  legend_y<-if(show_ne)1.28 else 3.30
  dot(xx,rep(legend_y,length(vv)),vv/100)
  for(i in seq_along(vv)){
    label(xx[i],legend_y-.40,as.character(vv[i]),pt=8)
  }
  label(mean(xx),legend_y+.60,'Tipping point (%)',pt=8.5)
  emit('\\end{tikzpicture}')
  dev.off()
  panels[[dec]]<-paste(tex,collapse='\n')
  cat(dec,':',nrow(z),'cities;',sum(in_ne),'repeated in Northeast inset\n')
}

# The supplement includes the two PNG maps (fig:tipping-map-1990, -2000).
cat('Exported two PNG figures:',
  file.path(here,c('tipping_map_1990_2000.png','tipping_map_2000_2010.png')),sep='\n')
