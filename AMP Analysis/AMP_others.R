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

# Supplemental analyses
source("TrtSurvCurves.R")
source("EIFestimates.R")

### Other SA countries
result.OA <- TrtSurvCurves(data=dat.hiv, 
                           tgt.name="African country other than South Africa",
                           prop.SL.library=c("SL.glm"), 
                           event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                           cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                           n.folds=5, 
                           s=2222)

result.OA.naive <- POOL_IVW(data=dat.hiv, 
                            tgt.name="South Africa",
                            prop.SL.library=c("SL.glm"), 
                            event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            n.folds=5, 
                            s=2388)

save(file="result_main_OA.Rdata", result.OA, result.OA.naive)

fit_cluster_cox <- coxph(
  Surv(Y, Delta) ~ A + age + score + bweight + strata(site),
  data = dat.hiv,
  cluster = site,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

dat <- subset(dat.hiv, site == "African country other than South Africa")
sf_cluster <- survfit(fit_cluster_cox, newdata = dat, se.fit = TRUE)

surv_list <- split(sf_cluster$surv, rep(seq_along(sf_cluster$strata), sf_cluster$strata))
se_list   <- split(sf_cluster$std.err, rep(seq_along(sf_cluster$strata), sf_cluster$strata))
time_list <- split(sf_cluster$time, rep(seq_along(sf_cluster$strata), sf_cluster$strata))

df.cluster <- data.frame(
  time = time_list[[1]],
  surv1 = surv_list[[1]],
  surv1.sd = se_list[[1]],
  surv0 = surv_list[[2]],
  surv0.sd = se_list[[2]]
)

load("result_main_OA.Rdata")
plot_survival_CI(df=result.OA$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt
plot_survival_CI(df=result.OA$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod
plot_survival_CI(df=result.OA$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed
plot_survival_CI(df = df.cluster, color_treated = "plum4", color_control = "maroon2", fig.title = "CLCOX") -> p.cluster
plot_survival_CI(df=result.OA$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw
plot_survival_CI(df=result.OA$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool

pdf(file="AMP_OA.pdf", width=12, height=6)
grid.arrange(p.fed, p.tgt, p.ccod, p.ivw, p.pool, p.cluster, ncol=3)
dev.off()

p.wt.0 <- plot_fedweights(result.OA, site.names=c("OA", "SA", "BP", "US"))$p.wt.0
p.wt.1 <- plot_fedweights(result.OA, site.names=c("OA", "SA", "BP", "US"))$p.wt.1

pdf(file="AMP_OAwts.pdf", width=10, height=3.4)
grid.arrange(p.wt.1, p.wt.0, ncol=2)
dev.off()


### Brazil or Peru
result.BP <- TrtSurvCurves(data=dat.hiv, 
                           tgt.name="Brazil or Peru",
                           prop.SL.library=c("SL.glm"), 
                           event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                           cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                           n.folds=5, 
                           s=2222)

result.BP.naive <- POOL_IVW(data=dat.hiv, 
                            tgt.name="South Africa",
                            prop.SL.library=c("SL.glm"), 
                            event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            n.folds=5, 
                            s=2388)

save(file="result_main_BP.Rdata", result.BP, result.BP.naive)

dat <- subset(dat.hiv, site == "Brazil or Peru")
sf_cluster <- survfit(fit_cluster_cox, newdata = dat, se.fit = TRUE)

surv_list <- split(sf_cluster$surv, rep(seq_along(sf_cluster$strata), sf_cluster$strata))
se_list   <- split(sf_cluster$std.err, rep(seq_along(sf_cluster$strata), sf_cluster$strata))
time_list <- split(sf_cluster$time, rep(seq_along(sf_cluster$strata), sf_cluster$strata))

df.cluster <- data.frame(
  time = time_list[[1]],
  surv1 = surv_list[[1]],
  surv1.sd = se_list[[1]],
  surv0 = surv_list[[2]],
  surv0.sd = se_list[[2]]
)

load("result_main_BP.Rdata")
plot_survival_CI(df=result.BP$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt
plot_survival_CI(df=result.BP$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod
plot_survival_CI(df=result.BP$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed
plot_survival_CI(df = df.cluster, color_treated = "plum4", color_control = "maroon2", fig.title = "CLCOX") -> p.cluster
plot_survival_CI(df=result.BP$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw
plot_survival_CI(df=result.BP$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool

pdf(file="AMP_BP.pdf", width=12, height=6)
grid.arrange(p.fed, p.tgt, p.ccod, p.ivw, p.pool, p.cluster, ncol=3)
dev.off()

p.wt.0 <- plot_fedweights(result.BP, site.names=c("BP", "OA", "SA", "US"))$p.wt.0
p.wt.1 <- plot_fedweights(result.BP, site.names=c("BP", "OA", "SA", "US"))$p.wt.1

pdf(file="AMP_BPwts.pdf", width=10, height=3.4)
grid.arrange(p.wt.1, p.wt.0, ncol=2)
dev.off()

### United States or Switzerland
result.US <- TrtSurvCurves(data=dat.hiv, 
                           tgt.name="United States or Switzerland",
                           prop.SL.library=c("SL.glm"), 
                           event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                           cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                           n.folds=5, 
                           s=2222)

result.US.naive <- POOL_IVW(data=dat.hiv, 
                            tgt.name="South Africa",
                            prop.SL.library=c("SL.glm"), 
                            event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            n.folds=5, 
                            s=2388)

save(file="result_main_US.Rdata", result.US, result.US.naive)

dat <- subset(dat.hiv, site == "United States or Switzerland")
sf_cluster <- survfit(fit_cluster_cox, newdata = dat, se.fit = TRUE)

surv_list <- split(sf_cluster$surv, rep(seq_along(sf_cluster$strata), sf_cluster$strata))
se_list   <- split(sf_cluster$std.err, rep(seq_along(sf_cluster$strata), sf_cluster$strata))
time_list <- split(sf_cluster$time, rep(seq_along(sf_cluster$strata), sf_cluster$strata))

df.cluster <- data.frame(
  time = time_list[[1]],
  surv1 = surv_list[[1]],
  surv1.sd = se_list[[1]],
  surv0 = surv_list[[2]],
  surv0.sd = se_list[[2]]
)

load("result_main_US.Rdata")
plot_survival_CI(df=result.US$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt
plot_survival_CI(df=result.US$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod
plot_survival_CI(df=result.US$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed
plot_survival_CI(df = df.cluster, color_treated = "plum4", color_control = "maroon2", fig.title = "CLCOX") -> p.cluster
plot_survival_CI(df=result.US$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw
plot_survival_CI(df=result.US$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool

pdf(file="AMP_US.pdf", width=12, height=6)
grid.arrange(p.fed, p.tgt, p.ccod, p.ivw, p.pool, p.cluster, ncol=3)
dev.off()

p.wt.0 <- plot_fedweights(result.US, site.names=c("US", "OA", "SA", "BP"))$p.wt.0
p.wt.1 <- plot_fedweights(result.US, site.names=c("US", "OA", "SA", "BP"))$p.wt.1

pdf(file="AMP_USwts.pdf", width=10, height=3.4)
grid.arrange(p.wt.1, p.wt.0, ncol=2)
dev.off()
