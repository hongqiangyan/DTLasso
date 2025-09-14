library(dplyr)
library(purrr)
library(purrr)
library(glmnet)

library(tsDyn)
library(dplyr)
library(ggplot2)
library(vars)
#dc <- dc %>%

library(desla) #0.3.0
library(ggplot2) #3.4.2
library(ggpubr) #0.6.0
library(ggpattern) #1.0.1
library(reshape2) #1.4.4

#RZ<-readxl::read_xls(system.file("extdata", "processed_data.xls", package="HDLPrepro", mustWork = TRUE)) #This will throw some warnings about dates in the wrong format. We don't use that variable so don't bother fixing it
#RZdf<-data.frame(quarter=RZ$quarter, newsy=RZ$newsy, g=RZ$g, y=RZ$y, taxy=RZ$taxy,
#lag_slack=dplyr::lag(RZ$slack), lag_zlb=dplyr::lag(RZ$zlb_dummy), lag_recession=dplyr::lag(RZ$recession), unemp = RZ$unemp)

#RZdf$interestrate[is.na(RZdf$interestrate)] <- 4
RZ<-readxl::read_xls("processed_data.xls") 
#This will throw some warnings about dates in the wrong format. We don't use that variable so don't bother fixing it
RZdf<-data.frame(quarter=RZ$quarter, newsy=RZ$newsy, g=RZ$g, y=RZ$y, taxy=RZ$taxy,
                 lag_slack=dplyr::lag(RZ$slack), lag_zlb=dplyr::lag(RZ$zlb_dummy), lag_recession=dplyr::lag(RZ$recession), unemp = RZ$unemp, interestrate = RZ$tbill)

RZdf$interestrate[is.na(RZdf$interestrate)] <- 4

# Remove NA values created by lagging
dc<-na.omit(RZdf)

# Create lagged variables for z_t
for (i in 1:40) {
  dc <- dc %>%
    mutate(!!paste0("newsy_t_lag_", i) := lag(newsy, i),
           !!paste0("y_t_lag_", i) := lag(y, i),
           !!paste0("g_t_lag_", i) := lag(g, i),
           !!paste0("taxy_t_lag_", i) := lag(taxy, i))
}

# Create the lead variable for y_t
dc <- dc %>%
  mutate(y_t_plus_1 = lead(g, 0))

# Prepare the regression formula
# Generate lagged variable names
newsy_lag_vars <- paste0("newsy_t_lag_", 1:40)
y_lag_vars <- paste0("y_t_lag_", 1:40)
g_lag_vars <- paste0("g_t_lag_", 1:40)
taxy_lag_vars <- paste0("taxy_t_lag_", 1:40)

# Combine all lagged variable names
all_lag_vars <- c(newsy_lag_vars, y_lag_vars, g_lag_vars,taxy_lag_vars)

#dc <- dc %>%
#mutate(
#low_interest = ifelse(interestrate < 0.36333, 1, 0),
#lag_unemp = lag(unemp)
#)

dc <- dc %>%
  mutate(
    #low_interest = ifelse(interestrate < 0.36333, 1, 0),
    lag_interestrate = lag(interestrate)
  )

X <- dc[, all_lag_vars, drop = FALSE]

#X <- cbind(dc$lag_unemp,dc$newsy,X)
X <- cbind(dc$lag_interestrate,dc$newsy,X)

mydata <- cbind(dc$y_t_plus_1,X)

mydata <- na.omit(mydata)

y <-  mydata[, 1]
X1 <- mydata[, -c(1,2)]

thdvar = mydata[,2]

X2 <- X1

#X2[thdvar>=6,] <- 0
#X2[thdvar>=0.5,] <- 0

quantiles <- quantile(thdvar, probs = c(0.10, 0.90))
quantile_10 <- quantiles[1]
quantile_90 <- quantiles[2]
interval <- c(quantile_10, quantile_90)
#0.170666 7.020640

#taugrid = seq(3,8,0.1)
sorted.thres <- sort(unique(thdvar))

trim.factor <- round(length(sorted.thres) * 0.1, 0)

taugrid <- sorted.thres[trim.factor:(length(sorted.thres) - trim.factor)]
X <- cbind(X1,X2)

X <- as.matrix(X)

y <- as.numeric(y)
####################################
#th <- SSLasso(X,y,taugrid=taugrid,verbose = TRUE)
#######################################
p <- ncol(X);
n <- nrow(X);
pp <- p;

H<-1

best_objective   <- Inf
best_tau         <- NA_real_
best_coefficients <- NULL
best_fit_obj     <- NULL
set.seed(222)
for (cnt in 1:length(taugrid)) {
  tau <- taugrid[cnt]
  X1 <- mydata[, -c(1,2)]
  X2 <- X1
  X2[thdvar>=tau,] <- 0
  xx <- cbind(X1, X2)
  xx <- as.matrix(xx)
  

  d <- desla(X = xx, y = y, H = H,PI_constant = 0.4)
  

  res_init <- as.numeric(d$residuals$init)
  beta     <- as.numeric(d$betahat)       
  lambda   <- as.numeric(d$lambdas$init)
  
 
  rss <- sum(res_init^2)
  Objective_lasso <- sqrt(rss / n) + lambda * sum(abs(beta))
  

  if (Objective_lasso < best_objective) {
    best_objective    <- Objective_lasso
    best_tau          <- tau
    best_coefficients <- beta
    best_fit_obj      <- d
  }
}


RZdf<-data.frame(quarter=RZ$quarter, newsy=RZ$newsy, g=RZ$g, y=RZ$y, taxy=RZ$taxy,
                 lag_slack=dplyr::lag(RZ$slack), lag_zlb=dplyr::lag(RZ$zlb_dummy), lag_recession=dplyr::lag(RZ$recession), unemp = RZ$unemp, interestrate = RZ$tbill)

RZdf$interestrate[is.na(RZdf$interestrate)] <- 4
dc<-na.omit(RZdf)

library(desla) #0.3.0
library(ggplot2) #3.4.2
library(ggpubr) #0.6.0
library(ggpattern) #1.0.1
library(reshape2) #1.4.4

dc <- dc %>%
  mutate(
    low_interest = ifelse(interestrate < best_tau, 1, 0),
    lag_low_interest = lag(low_interest)
  ) %>%
  na.omit()
hmax=20 # maximum horizon - the x axis of the plot will be 0:hmax
lags=40 # number of lags included in the local projection equations
PIconstant=0.4 # this is the plug-in constant used for the data-dependent selection of the lasso penalization. Generally, higher value gives stronger penalization. For details, see Algorithm 1 in the supplementary appendix C.5 of https://doi.org/10.1016/j.jeconom.2022.08.008
threads<-parallel::detectCores()-2 # the number of cores used in parallel computation 
# seed which controls the random number generation for reproducibility. We use set.seed(1) in our application
################################################################################
# estimating these HDLPs can take a few minutes

linear_g<-HDLP(x=dc$newsy, y=dc$g, q=cbind(dc$y,dc$taxy), lags=lags, hmax=hmax)
linear_y<-HDLP(x=dc$newsy, y=dc$y, q=cbind(dc$g,dc$taxy), lags=lags, hmax=hmax)
nl_g<-HDLP(x=dc$newsy, y=dc$g, q=cbind(dc$y,dc$taxy), state_variables=dc$lag_low_interest, lags=lags, hmax=hmax)
nl_y<-HDLP(x=dc$newsy, y=dc$y, q=cbind(dc$g,dc$taxy), state_variables=dc$lag_low_interest, lags=lags, hmax=hmax)

# code for Figure 4
g_df<-data.frame(mean_all=linear_g$intervals[,2,1], lb_all=linear_g$intervals[,1,1], ub_all=linear_g$intervals[,3,1],
                 mean_hu=nl_g$intervals[,2,1], lb_hu=nl_g$intervals[,1,1], ub_hu=nl_g$intervals[,3,1],
                 mean_lu=nl_g$intervals[,2,2], lb_lu=nl_g$intervals[,1,2], ub_lu=nl_g$intervals[,3,2],
                 Horizon=0:hmax)
y_df<-data.frame(mean_all=linear_y$intervals[,2,1], lb_all=linear_y$intervals[,1,1], ub_all=linear_y$intervals[,3,1],
                 mean_hu=nl_y$intervals[,2,1], lb_hu=nl_y$intervals[,1,1], ub_hu=nl_y$intervals[,3,1],
                 mean_lu=nl_y$intervals[,2,2], lb_lu=nl_y$intervals[,1,2], ub_lu=nl_y$intervals[,3,2],
                 Horizon=0:hmax)

col2="blue"
col1="red"


p2<-ggplot(data=g_df, aes(x=Horizon)) +
  geom_hline(yintercept = 0, color="black", linewidth=0.5)+
  geom_line(aes(y=mean_all), color = "black", linewidth = 1)+
  scale_colour_manual(values = c("black",col1, col2), labels=c("Linear model","Normal", "ZLB"), breaks = c("linear","high", "low"),
                      guide = guide_legend(direction = "horizontal", title.position = "top",
                                           label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                           label.theme = element_text(angle = 90)))+
  geom_ribbon(aes(ymin = lb_all, ymax = ub_all), fill = "black", alpha=0.25)+
  ylim(c(min(c(g_df$lb_all,g_df$lb_hu,g_df$lb_lu)),max(c(g_df$ub_all,g_df$ub_hu,g_df$ub_lu))))+
  xlab("Horizon") +
  ylab("Government Spending")+theme_bw()+
  labs(title="Linear model")+
  theme(plot.title = element_text(face="bold", size=14, hjust=0.5))

df<-data.frame(Horizon=0:hmax, high=g_df$mean_hu, low=g_df$mean_lu)
df_melt<-melt(df, id="Horizon"); names(df_melt)[2]<-"State"
p3<-ggplot(data=g_df, aes(x=Horizon)) +
  geom_hline(yintercept = 0, color="black", linewidth=0.5)+
  geom_line(data=df_melt, mapping=aes(y=value, color=State), linewidth = 0.75)+
  scale_colour_manual(values = c("black",col1, col2), labels=c("Linear model","Normal", "ZLB"), breaks = c("linear","high", "low"),
                      guide = guide_legend(direction = "horizontal", title.position = "top",
                                           label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                           label.theme = element_text(angle = 90)))+
  geom_ribbon_pattern(aes(ymin = lb_hu, ymax = ub_hu), fill = col1, alpha=0.1,
                      color = col1,
                      pattern = "circle",
                      pattern_color=col1,
                      pattern_fill = col1,
                      pattern_angle = 90,
                      pattern_density = 0.2,
                      pattern_spacing = 0.05) +
  geom_ribbon_pattern(aes(ymin = lb_lu, ymax = ub_lu), fill = col2, alpha=0.1,
                      color = col2,
                      pattern = "none",
                      pattern_color=col2,
                      pattern_fill = col2,
                      pattern_angle = 90,
                      pattern_density = 0.001,
                      pattern_spacing = 0.05) +
  ylim(c(min(c(g_df$lb_all,g_df$lb_hu,g_df$lb_lu)),max(c(g_df$ub_all,g_df$ub_hu,g_df$ub_lu))))+
  xlab("Horizon") +
  ylab("Government Spending")+theme_bw()+
  labs(title="HDLPT")+
  theme(plot.title = element_text(face="bold", size=14, hjust=0.5))

p5<-ggplot(data=y_df, aes(x=Horizon)) +
  geom_hline(yintercept = 0, color="black", linewidth=0.5)+
  geom_line(aes(y=mean_all), color = "black", linewidth = 1)+
  scale_colour_manual(values = c("black",col1, col2), labels=c("Linear model","Normal", "ZLB"), breaks = c("linear","high", "low"),
                      guide = guide_legend(direction = "horizontal", title.position = "top",
                                           label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                           label.theme = element_text(angle = 90)))+
  geom_ribbon(aes(ymin = lb_all, ymax = ub_all), fill = "black", alpha=0.25)+
  ylim(c(min(c(y_df$lb_all,y_df$lb_hu,y_df$lb_lu)),max(c(y_df$ub_all,y_df$ub_hu,y_df$ub_lu))))+
  xlab("Horizon") +
  ylab("GDP")+theme_bw()

df<-data.frame(Horizon=0:hmax, high=y_df$mean_hu, low=y_df$mean_lu)
df_melt<-melt(df, id="Horizon"); names(df_melt)[2]<-"State"
p6<-ggplot(data=y_df, aes(x=Horizon)) +
  geom_hline(yintercept = 0, color="black", linewidth=0.5)+
  geom_line(data=df_melt, mapping=aes(y=value, color=State), linewidth = 1)+
  scale_colour_manual(values = c(col1, col2), labels=c("Normal", "ZLB"), breaks = c("high", "low"),
                      guide = guide_legend(direction = "horizontal", title.position = "top",
                                           label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                           label.theme = element_text(angle = 90)))+
  geom_ribbon_pattern(aes(ymin = lb_hu, ymax = ub_hu), fill = col1, alpha=0.1,
                      color = col1,
                      pattern = "circle",
                      pattern_color=col1,
                      pattern_fill = col1,
                      pattern_angle = 90,
                      pattern_density = 0.2,
                      pattern_spacing = 0.05) +
  geom_ribbon_pattern(aes(ymin = lb_lu, ymax = ub_lu), fill = col2, alpha=0.1,
                      color = col2,
                      pattern = "none",
                      pattern_color=col2,
                      pattern_fill = col2,
                      pattern_angle = 90,
                      pattern_density = 0.001,
                      pattern_spacing = 0.05) +
  ylim(c(min(c(y_df$lb_all,y_df$lb_hu,y_df$lb_lu)),max(c(y_df$ub_all,y_df$ub_hu,y_df$ub_lu))))+
  xlab("Horizon") +
  ylab("GDP")+theme_bw()

combined1<-ggarrange(p2,p3,p5,p6, common.legend = TRUE, legend="right")



