# plot functions
plot_survival_CI <- function(
    df, 
    time_col           ="time",
    surv_treated_col   ="surv1",
    surv_treated_sd_col="surv1.sd",
    surv_control_col   ="surv0",
    surv_control_sd_col="surv0.sd",
    ci_multiplier      =1.96,
    color_treated      ="blue",
    color_control      ="red",
    alpha_fill         =0.2,
    line_size          =0.8,
    xlab               ="Time (day)",
    ylab               ="Survival Probability",
    fig.title          ="CCOD"
){
  ggplot(df, aes_string(x=time_col)) +
    geom_line(aes_string(y=surv_treated_col, color="'Treated'"), size=line_size) +
    geom_ribbon(aes_string(
      ymin=paste0(surv_treated_col, " - ", ci_multiplier, " * ", surv_treated_sd_col),
      ymax=paste0(surv_treated_col, " + ", ci_multiplier, " * ", surv_treated_sd_col),
      fill="'Treated'"), alpha=alpha_fill) +
    geom_line(aes_string(y=surv_control_col, color="'Control'"), size=line_size) +
    geom_ribbon(aes_string(
      ymin=paste0(surv_control_col, " - ", ci_multiplier, " * ", surv_control_sd_col),
      ymax=paste0(surv_control_col, " + ", ci_multiplier, " * ", surv_control_sd_col),
      fill="'Control'"), alpha=alpha_fill) +
    scale_color_manual(name ="Group", values=c("Treated"=color_treated, "Control"=color_control)) +
    scale_fill_manual(name  ="Group", values=c("Treated"=color_treated, "Control"=color_control)) +
    labs(x=xlab, y=ylab, title=fig.title) + ylim(0.88,1) + 
    theme_minimal()
}

plot_fedweights <- function(results.list, site.names=c("SA", "OA", "BP", "US")) {
  df1 <- results.list$weights[,1:4] %>% as.data.frame()
  df0 <- results.list$weights[,5:8] %>% as.data.frame()
  df1$time <- seq_len(nrow(df1)) 
  df0$time <- seq_len(nrow(df0))
  colnames(df1) <- colnames(df0) <- c(site.names, "time")
  
  df_long1 <- df1 %>%
    pivot_longer(cols=-time, names_to="Region", values_to="Weight")
  df_long0 <- df0 %>%
    pivot_longer(cols=-time, names_to="Region", values_to="Weight")
  df_long1$Region <- factor(df_long1$Region, levels=c("SA", "OA", "BP", "US"))
  df_long0$Region <- factor(df_long0$Region, levels=c("SA", "OA", "BP", "US"))
  
  ggplot(df_long1, aes(x=time, y=Weight, color=Region)) +
    geom_line(alpha=0.5) + ylim(0,1) +            
    geom_smooth(se=FALSE, size=1) +
    labs(
      title="Treated group",
      x="Time (day)",
      y="Federated weight"
    ) +
    theme_minimal(base_size=14) +
    theme(legend.position="top") -> p.wt.1   
  
  ggplot(df_long0, aes(x=time, y=Weight, color=Region)) +
    geom_line(alpha=0.5) + ylim(0,1) +                      
    geom_smooth(se=FALSE, size=1) +
    labs(
      title="Control group",
      x="Time (day)",
      y="Federated weight"
    ) +
    theme_minimal(base_size=14) +
    theme(legend.position="top") -> p.wt.0
  
  return(list(p.wt.1=p.wt.1, p.wt.0=p.wt.0))
}

# data pre-processing
dat <- read.csv("amp_survival.csv")
head(dat)
mean(dat$hiv1event)
dat$hiv1event
library(dplyr)
library(ggplot2)
library(gridExtra)
library(SuperLearner)
library(survSuperLearner)
library(CFsurvival)
library(glmnet)
library(caret)
library(tidyr)
library(xtable)
library(survival)
source("TrtSurvCurves.R")
source("EIFestimates.R")

Delta <- dat$hiv1event
mean(Delta)
Y <- dat$hiv1survday
mean(Y<601)
sum(Y==601)
mean(Y==601)
A <- as.numeric(dat$rx_pool=="T1+T2")
sum(Y==601 & Delta==1)
X <- data.frame(bweight=dat$bweight, score=dat$standardized_risk_score, age=dat$bbmi)

unique(dat$country)
site <- dat$country
site[site%in%c("Tanzania, Mozambique, Kenya", "Zimbabwe", "Botswana", "Malawi")] <- "African country other than South Africa"
site[site%in%c("Peru", "Brazil")] <- "Brazil or Peru"
site[site%in%c("United States", "Switzerland")] <- "United States or Switzerland"
site <- factor(site, levels=c("South Africa", "African country other than South Africa", "Brazil or Peru", "United States or Switzerland"))
unique(site)

dat.hiv <- data.frame(cbind(A, Y, Delta, X, site))
head(dat.hiv)
range(dat.hiv$Y)

result.SA <- TrtSurvCurves(data=dat.hiv, 
                           covar.name=c("age","bweight"),
                           tgt.name="South Africa", 
                           n.folds=5, 
                           s=2388)

result.SA.naive <- POOL_IVW(data=dat.hiv, 
                            covar.name=c("age","bweight"),
                            tgt.name="South Africa",
                            prop.SL.library=c("SL.glm"), 
                            event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            n.folds=5, 
                            s=2388)
save(file="result_main_SA_2.Rdata", result.SA, result.SA.naive)

fit_cluster_cox <- coxph(
  Surv(Y, Delta) ~ A + age + bweight + strata(site),
  data = dat.hiv,
  cluster = site,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

dat <- subset(dat.hiv, site == "South Africa")
sf_cluster <- survfit(fit_cluster_cox, newdata = dat, se.fit = TRUE)

surv_list <- split(sf_cluster$surv, rep(seq_along(sf_cluster$strata), sf_cluster$strata))
se_list   <- split(sf_cluster$std.err, rep(seq_along(sf_cluster$strata), sf_cluster$strata))
time_list <- split(sf_cluster$time, rep(seq_along(sf_cluster$strata), sf_cluster$strata))

df.cluster <- data.frame(
  time = time_list[[1]],
  surv1 = surv_list[[2]],
  surv1.sd = se_list[[2]],
  surv0 = surv_list[[1]],
  surv0.sd = se_list[[1]]
)

### Results plotting and summary for South Africa (main results in the paper)
load("result_main_SA_2.Rdata")
plot_survival_CI(df=result.SA$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt
plot_survival_CI(df=result.SA$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod
plot_survival_CI(df=result.SA$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed
plot_survival_CI(df = df.cluster, color_treated = "plum4", color_control = "maroon2", fig.title = "CLCOX") -> p.cluster
plot_survival_CI(df=result.SA.naive$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw
plot_survival_CI(df=result.SA.naive$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool

pdf(file="AMP_SA_2.pdf", width=12, height=6)
grid.arrange(p.fed, p.tgt, p.ccod, p.ivw, p.pool, p.cluster, ncol=3)
dev.off()

p.wt.SA0 <- plot_fedweights(result.SA)$p.wt.0
p.wt.SA1 <- plot_fedweights(result.SA)$p.wt.1

pdf(file="AMP_SAwts_2.pdf", width=10, height=3.4)
grid.arrange(p.wt.SA1, p.wt.SA0, ncol=2)
dev.off()

### Relative efficiency

### RD and RMST
source("Extends.R")
load("result_main_SA_2.Rdata")
extend.results <- FuseSurv_Extend(eval.times=result.SA$eval.times,
                                  site=result.SA$site,
                                  IF.00=result.SA$IF.00,
                                  IF.01=result.SA$IF.01,
                                  S.00=result.SA$S.00, 
                                  S.01=result.SA$S.01,
                                  Aug.00.mean=result.SA$Aug.00.mean, 
                                  Aug.01.mean=result.SA$Aug.01.mean,
                                  Aug.R0.mean=result.SA$Aug.R0.mean, 
                                  Aug.R1.mean=result.SA$Aug.R1.mean,
                                  Aug.R0.mean.sour=result.SA$Aug.R0.mean.sour, 
                                  Aug.R1.mean.sour=result.SA$Aug.R1.mean.sour,
                                  IF.R0=result.SA$IF.R0, 
                                  IF.R1=result.SA$IF.R1, 
                                  IF.CCOD.0=result.SA$IF.CCOD.0, 
                                  IF.CCOD.1=result.SA$IF.CCOD.1,
                                  ind.R1.ccod=result.SA$ind.R1.ccod,
                                  s=1)

## selected time points
time.sel <- c(148, 330, 512)

## subset
df.RD.TGT  <- extend.results$df.RD.TGT [extend.results$df.RD.TGT$time %in% time.sel, ]
df.RD.FED  <- extend.results$df.RD.FED [extend.results$df.RD.FED$time %in% time.sel, ]
df.RD.CCOD <- extend.results$df.RD.CCOD[extend.results$df.RD.CCOD$time %in% time.sel, ]

df.SR.TGT  <- extend.results$df.SR.TGT [extend.results$df.SR.TGT$time %in% time.sel, ]
df.SR.FED  <- extend.results$df.SR.FED [extend.results$df.SR.FED$time %in% time.sel, ]
df.SR.CCOD <- extend.results$df.SR.CCOD[extend.results$df.SR.CCOD$time %in% time.sel, ]

## two-sided Wald p-values
## RD: H0 = 0
df.RD.TGT$pval  <- 2 * pnorm(-abs(df.RD.TGT$RD  / df.RD.TGT$sd))
df.RD.FED$pval  <- 2 * pnorm(-abs(df.RD.FED$RD  / df.RD.FED$sd))
df.RD.CCOD$pval <- 2 * pnorm(-abs(df.RD.CCOD$RD / df.RD.CCOD$sd))

## SR: H0 = 1
df.SR.TGT$pval  <- 2 * pnorm(-abs((df.SR.TGT$SR  - 1) / df.SR.TGT$sd))
df.SR.FED$pval  <- 2 * pnorm(-abs((df.SR.FED$SR  - 1) / df.SR.FED$sd))
df.SR.CCOD$pval <- 2 * pnorm(-abs((df.SR.CCOD$SR - 1) / df.SR.CCOD$sd))

## display
round(df.RD.TGT,  3)
round(df.RD.FED,  3)
round(df.RD.CCOD, 3)

round(df.SR.TGT,  3)
round(df.SR.FED,  3)
round(df.SR.CCOD, 3)

extend.results$df.RMST.diff.TGT
extend.results$df.RMST.diff.FED
extend.results$df.RMST.diff.CCOD

p_rmst_diff_tgt  <- 2 * pnorm(-abs(extend.results$df.RMST.diff.TGT$RMST  / extend.results$df.RMST.diff.TGT$sd))
p_rmst_diff_fed  <- 2 * pnorm(-abs(extend.results$df.RMST.diff.FED$RMST  / extend.results$df.RMST.diff.FED$sd))
p_rmst_diff_ccod <- 2 * pnorm(-abs(extend.results$df.RMST.diff.CCOD$RMST / extend.results$df.RMST.diff.CCOD$sd))

round(c(TGT = p_rmst_diff_tgt,
        FED = p_rmst_diff_fed,
        CCOD = p_rmst_diff_ccod), 3)

extend.results$df.RMST.0.TGT
extend.results$df.RMST.0.FED
extend.results$df.RMST.0.CCOD

extend.results$df.RMST.1.TGT
extend.results$df.RMST.1.FED
extend.results$df.RMST.1.CCOD

## add p-values for RMST(0), RMST(1), RMST diff
## RMST.0 / RMST.1: H0 = 0
extend.results$df.RMST.0.TGT$pval  <- 2 * pnorm(-abs(extend.results$df.RMST.0.TGT$RMST  / extend.results$df.RMST.0.TGT$sd))
extend.results$df.RMST.0.FED$pval  <- 2 * pnorm(-abs(extend.results$df.RMST.0.FED$RMST  / extend.results$df.RMST.0.FED$sd))
extend.results$df.RMST.0.CCOD$pval <- 2 * pnorm(-abs(extend.results$df.RMST.0.CCOD$RMST / extend.results$df.RMST.0.CCOD$sd))

extend.results$df.RMST.1.TGT$pval  <- 2 * pnorm(-abs(extend.results$df.RMST.1.TGT$RMST  / extend.results$df.RMST.1.TGT$sd))
extend.results$df.RMST.1.FED$pval  <- 2 * pnorm(-abs(extend.results$df.RMST.1.FED$RMST  / extend.results$df.RMST.1.FED$sd))
extend.results$df.RMST.1.CCOD$pval <- 2 * pnorm(-abs(extend.results$df.RMST.1.CCOD$RMST / extend.results$df.RMST.1.CCOD$sd))

## RMST diff: H0 = 0
extend.results$df.RMST.diff.TGT$pval  <- 2 * pnorm(-abs(extend.results$df.RMST.diff.TGT$RMST  / extend.results$df.RMST.diff.TGT$sd))
extend.results$df.RMST.diff.FED$pval  <- 2 * pnorm(-abs(extend.results$df.RMST.diff.FED$RMST  / extend.results$df.RMST.diff.FED$sd))
extend.results$df.RMST.diff.CCOD$pval <- 2 * pnorm(-abs(extend.results$df.RMST.diff.CCOD$RMST / extend.results$df.RMST.diff.CCOD$sd))

add_CI <- function(df, est_col, sd_col, digits_ci = 3, digits_p = 3) {
  df$lower <- df[[est_col]] - 1.96 * df[[sd_col]]
  df$upper <- df[[est_col]] + 1.96 * df[[sd_col]]
  df$CI <- sprintf(paste0("%.", digits_ci, "f, %.", digits_ci, "f"), df$lower, df$upper)
  df$CI <- paste0("(", df$CI, ")")
  df$pval <- round(df$pval, digits_p)
  df
}

df.RD <- rbind(
  cbind(Method = "TGT",  add_CI(df.RD.TGT,  "RD", "sd", digits_ci = 3, digits_p = 3)),
  cbind(Method = "FED",  add_CI(df.RD.FED,  "RD", "sd", digits_ci = 3, digits_p = 3)),
  cbind(Method = "CCOD", add_CI(df.RD.CCOD, "RD", "sd", digits_ci = 3, digits_p = 3))
)[, c("Method", "time", "RD", "sd", "CI", "pval")]

xtable(df.RD, digits = c(0, 0, 0, 3, 3, 0, 3))

df.SR <- rbind(
  cbind(Method = "TGT",  add_CI(df.SR.TGT,  "SR", "sd", digits_ci = 3, digits_p = 3)),
  cbind(Method = "FED",  add_CI(df.SR.FED,  "SR", "sd", digits_ci = 3, digits_p = 3)),
  cbind(Method = "CCOD", add_CI(df.SR.CCOD, "SR", "sd", digits_ci = 3, digits_p = 3))
)[, c("Method", "time", "SR", "sd", "CI", "pval")]

xtable(df.SR, digits = c(0, 0, 0, 3, 3, 0, 3))

add_CI <- function(df, est_col, sd_col, digits_ci = 2, digits_p = 3) {
  df$lower <- df[[est_col]] - 1.96 * df[[sd_col]]
  df$upper <- df[[est_col]] + 1.96 * df[[sd_col]]
  df$CI <- sprintf(paste0("%.", digits_ci, "f, %.", digits_ci, "f"), df$lower, df$upper)
  df$CI <- paste0("(", df$CI, ")")
  df$pval <- round(df$pval, digits_p)
  df
}

df.RMST.0 <- rbind(
  data.frame(Method = "TGT",  add_CI(extend.results$df.RMST.0.TGT,   "RMST", "sd", digits_ci = 2, digits_p = 3)),
  data.frame(Method = "FED",  add_CI(extend.results$df.RMST.0.FED,   "RMST", "sd", digits_ci = 2, digits_p = 3)),
  data.frame(Method = "CCOD", add_CI(extend.results$df.RMST.0.CCOD,  "RMST", "sd", digits_ci = 2, digits_p = 3))
)[, c("Method", "time.max", "RMST", "sd", "CI", "pval")]

xtable(df.RMST.0, digits = c(0, 0, 0, 2, 2, 0, 3))
 
df.RMST.1 <- rbind(
  data.frame(Method = "TGT",  add_CI(extend.results$df.RMST.1.TGT,   "RMST", "sd", digits_ci = 2, digits_p = 3)),
  data.frame(Method = "FED",  add_CI(extend.results$df.RMST.1.FED,   "RMST", "sd", digits_ci = 2, digits_p = 3)),
  data.frame(Method = "CCOD", add_CI(extend.results$df.RMST.1.CCOD,  "RMST", "sd", digits_ci = 2, digits_p = 3))
)[, c("Method", "time.max", "RMST", "sd", "CI", "pval")]

xtable(df.RMST.1, digits = c(0, 0, 0, 2, 2, 0, 3)) 

df.RMST.diff <- rbind(
  data.frame(Method = "TGT",  add_CI(extend.results$df.RMST.diff.TGT,  "RMST", "sd", digits_ci = 2, digits_p = 3)),
  data.frame(Method = "FED",  add_CI(extend.results$df.RMST.diff.FED,  "RMST", "sd", digits_ci = 2, digits_p = 3)),
  data.frame(Method = "CCOD", add_CI(extend.results$df.RMST.diff.CCOD, "RMST", "sd", digits_ci = 2, digits_p = 3))
)[, c("Method", "time.max", "RMST", "sd", "CI", "pval")]

xtable(df.RMST.diff, digits = c(0, 0, 0, 2, 2, 0, 3))