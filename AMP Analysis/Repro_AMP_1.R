### ----------------------------------------------------------------------------
### Reproducibility for AMP Data Analysis: Part I
### --- This script reproduces Table 1--3 in the main text,
### --- Figures 1 and 5 in the main text, and Figures A.1--A.3 and
### --- Tables A.1--A.6 in Online Supplemental Material Appendix A.
### ----------------------------------------------------------------------------

### --- Packages and data pre-processing
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
library(survminer)
library(EValue)

### Read the AMP trials dataset and construct the analysis variables.
dat <- read.csv("amp_survival.csv")
Delta <- dat$hiv1event
Y <- dat$hiv1survday
A <- as.numeric(dat$rx_pool == "T1+T2")
X <- data.frame(
  bweight = dat$bweight,
  score   = dat$standardized_risk_score,
  age     = dat$bbmi
)
site <- dat$country
site[site%in%c("Tanzania, Mozambique, Kenya", "Zimbabwe", "Botswana", "Malawi")] <- "African country other than South Africa"
site[site%in%c("Peru", "Brazil")] <- "Brazil or Peru"
site[site%in%c("United States", "Switzerland")] <- "United States or Switzerland"
site <- factor(site, levels=c("South Africa", "African country other than South Africa", "Brazil or Peru", "United States or Switzerland"))
unique(site)
dat.hiv <- data.frame(cbind(A, Y, Delta, X, site))

### ----------------------------------------------------------------------------
#### Reproducibility 1: Numbers in Table 1 in Section 2
### ----------------------------------------------------------------------------
fmt_mean_sd <- function(mean, sd, digits=2) sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean, sd)
fmt_count_pct <- function(count, rate, digits=2) sprintf(paste0("%d (%.", digits, "f%%)"), count, 100 * rate)
make_amp_table1 <- function(dat.hiv) {
  site_levels <- levels(dat.hiv$site)
  site_labels <- c(
    "South Africa" = "SA (women)",
    "African country other than South Africa" = "OA (women)",
    "Brazil or Peru" = "BP (men, TG)",
    "United States or Switzerland" = "US (men, TG)"
  )
  
  make_block <- function(a) {
    dat.a <- dat.hiv %>% filter(A == a)
    total <- dat.a %>% summarize(
      n=n(), age.mean=mean(age), age.sd=sd(age),
      bweight.mean=mean(bweight), bweight.sd=sd(bweight),
      score.mean=mean(score), score.sd=sd(score),
      event.count=sum(Delta), event.rate=mean(Delta)
    )
    bysite <- dat.a %>% group_by(site) %>% summarize(
      n=n(), age.mean=mean(age), age.sd=sd(age),
      bweight.mean=mean(bweight), bweight.sd=sd(bweight),
      score.mean=mean(score), score.sd=sd(score),
      event.count=sum(Delta), event.rate=mean(Delta), .groups="drop"
    )
    
    vals <- c("Total", unname(site_labels[site_levels]))
    out <- data.frame(Variable=c("Age (years)", "Weight (kg)", "ML risk score", "HIV-1 diagnosis"), check.names=FALSE)
    for (j in seq_along(vals)) out[[vals[j]]] <- NA_character_
    
    fill_col <- function(col, x) {
      out[out$Variable=="Age (years)", col] <<- fmt_mean_sd(x$age.mean, x$age.sd, 2)
      out[out$Variable=="Weight (kg)", col] <<- fmt_mean_sd(x$bweight.mean, x$bweight.sd, 2)
      out[out$Variable=="ML risk score", col] <<- fmt_mean_sd(x$score.mean, x$score.sd, 2)
      out[out$Variable=="HIV-1 diagnosis", col] <<- fmt_count_pct(x$event.count, x$event.rate, 2)
    }
    
    fill_col("Total", total)
    for (s in site_levels) fill_col(site_labels[[s]], bysite %>% filter(site == s))
    out
  }
  
  tab.treated <- make_block(1)
  tab.control <- make_block(0)
  
  cat("\nTreated (bnAb) group\n")
  print(tab.treated, row.names=FALSE)
  cat("\nControl (placebo) group\n")
  print(tab.control, row.names=FALSE)
  
  invisible(list(treated=tab.treated, control=tab.control))
}

### --- Print Table 1 in Section 2
make_amp_table1(dat.hiv)

### ----------------------------------------------------------------------------
### Reproducibility 2: Sensitivity analysis for censoring in Section 5.1
### ----------------------------------------------------------------------------

### Compute the E-value reported in Section 5.1.
### Calibrate the E-value by examining associations of standardized observed
### covariates with the HIV-1 outcome and with censoring.
model <- coxph(Surv(Y, Delta)~A+age+score+bweight, data=dat.hiv)
summary_fit <- summary(model)$conf.int
HR <- as.numeric(summary_fit["A", "exp(coef)"])
lo_CI <- as.numeric(summary_fit["A", "lower .95"])
hi_CI <- as.numeric(summary_fit["A", "upper .95"])
evalues.RR(est = HR, lo = lo_CI, hi = hi_CI, true=1)[2,1] # E-value reported in Section 5.1

dat.cal <- dat.hiv %>%
  mutate(
    age_z = as.numeric(scale(age)),
    score_z = as.numeric(scale(score)),
    bweight_z = as.numeric(scale(bweight)),
    censor_event = 1 - Delta
  )

## 1) Association of observed covariates with outcome
fit.outcome <- coxph(
  Surv(Y, Delta) ~ age_z + score_z + bweight_z,
  data = dat.cal
)
summary(fit.outcome)

## 2) Association of observed covariates with censoring
## Treat censoring as the event
fit.censor <- coxph(
  Surv(Y, censor_event) ~ age_z + score_z + bweight_z,
  data = dat.cal
)

summary(fit.censor)

### Extract hazard ratios and 95% confidence intervals from a Cox model.
get_hr_table <- function(fit, model_name) {
  s <- summary(fit)
  out <- data.frame(
    Variable = rownames(s$conf.int),
    HR = s$conf.int[, "exp(coef)"],
    Lower = s$conf.int[, "lower .95"],
    Upper = s$conf.int[, "upper .95"],
    Model = model_name,
    row.names = NULL
  )
  out
}

tab.outcome <- get_hr_table(fit.outcome, "Outcome")
tab.censor  <- get_hr_table(fit.censor, "Censoring")
tab.calib <- rbind(tab.outcome, tab.censor)
print(tab.calib)

### ----------------------------------------------------------------------------
### Reproducibility 3: Figure 5 in Section 5.2
### ----------------------------------------------------------------------------

### Estimate treatment-specific survival curves with SA as the target region.
### TGT, FED, and CCOD are obtained from TrtSurvCurves(), while POOL and IVW
### are obtained from POOL_IVW().
source("TrtSurvCurves.R")
source("EIFestimates.R")

result.SA <- TrtSurvCurves(data=dat.hiv,
                           tgt.name="South Africa",
                           prop.SL.library=c("SL.glm"),
                           event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                           cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                           n.folds=5,
                           s=2388)

result.SA.naive <- POOL_IVW(data=dat.hiv,
                            tgt.name="South Africa",
                            prop.SL.library=c("SL.glm"),
                            event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                            n.folds=5,
                            s=2388)
save(file="result_main_SA.Rdata", result.SA, result.SA.naive)

### Fit the stratified clustered Cox model used as an additional comparator.
fit_cluster_cox <- coxph(
  Surv(Y, Delta) ~ A + age + score + bweight + strata(site),
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

### Plot treatment-specific survival curves and pointwise 95% confidence intervals.
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

### Plot the estimated federated weights over time for each analytical region.
plot_fedweights <- function(results.list, site.names=c("SA", "OA", "BP", "US")) {
  df0 <- results.list$weights[,1:4] %>% as.data.frame()
  df1 <- results.list$weights[,5:8] %>% as.data.frame()
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

load("result_main_SA.Rdata")
plot_survival_CI(df=result.SA$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt
plot_survival_CI(df=result.SA$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod
plot_survival_CI(df=result.SA$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed
plot_survival_CI(df = df.cluster, color_treated = "plum4", color_control = "maroon2", fig.title = "CLCOX") -> p.cluster
plot_survival_CI(df=result.SA.naive$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw
plot_survival_CI(df=result.SA.naive$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool

### --- Below code prints panels (A) and (B) of Figure 5
pdf(file="AMP_SA.pdf", width=12, height=6)
grid.arrange(p.fed, p.tgt, p.ccod, p.ivw, p.pool, p.cluster, ncol=3)
dev.off()

p.wt.SA0 <- plot_fedweights(result.SA)$p.wt.0
p.wt.SA1 <- plot_fedweights(result.SA)$p.wt.1

pdf(file="AMP_SAwts.pdf", width=10, height=3.4)
grid.arrange(p.wt.SA1, p.wt.SA0, ncol=2)
dev.off()

### --- Below code shows numbers in panel (A) of Figure 5
### --- These numbers are manually added in Figure 5 separately using Microsoft PowerPoint
load("result_main_SA.Rdata")
### --- Below are the 6 numbers shown in the FED panel of Figure 5(A)
result.SA$df.FED[148,3] / result.SA$df.TGT[148,3]
result.SA$df.FED[148,5] / result.SA$df.TGT[148,5]
result.SA$df.FED[330,3] / result.SA$df.TGT[330,3]
result.SA$df.FED[330,5] / result.SA$df.TGT[330,5]
result.SA$df.FED[512,3] / result.SA$df.TGT[512,3]
result.SA$df.FED[512,5] / result.SA$df.TGT[512,5]

### --- Below are the 6 numbers shown in the CCOD panel of Figure 5(A)
result.SA$df.CCOD[148,3] / result.SA$df.TGT[148,3]
result.SA$df.CCOD[148,5] / result.SA$df.TGT[148,5]
result.SA$df.CCOD[330,3] / result.SA$df.TGT[330,3]
result.SA$df.CCOD[330,5] / result.SA$df.TGT[330,5]
result.SA$df.CCOD[512,3] / result.SA$df.TGT[512,3]
result.SA$df.CCOD[512,5] / result.SA$df.TGT[512,5]

### --- Below are the 6 numbers shown in the IVW panel of Figure 5(A)
result.SA.naive$df.IVW[148,3] / result.SA$df.TGT[148,3]
result.SA.naive$df.IVW[148,5] / result.SA$df.TGT[148,5]
result.SA.naive$df.IVW[330,3] / result.SA$df.TGT[330,3]
result.SA.naive$df.IVW[330,5] / result.SA$df.TGT[330,5]
result.SA.naive$df.IVW[512,3] / result.SA$df.TGT[512,3]
result.SA.naive$df.IVW[512,5] / result.SA$df.TGT[512,5]

### --- Below are the 6 numbers shown in the POOL panel of Figure 5(A)
result.SA.naive$df.POOL[148,3] / result.SA$df.TGT[148,3]
result.SA.naive$df.POOL[148,5] / result.SA$df.TGT[148,5]
result.SA.naive$df.POOL[330,3] / result.SA$df.TGT[330,3]
result.SA.naive$df.POOL[330,5] / result.SA$df.TGT[330,5]
result.SA.naive$df.POOL[512,3] / result.SA$df.TGT[512,3]
result.SA.naive$df.POOL[512,5] / result.SA$df.TGT[512,5]

### --- Below are the 6 numbers shown in the CLCOX panel of Figure 5(A)
df.cluster[df.cluster$time==148,3] / result.SA$df.TGT[148,3]
df.cluster[df.cluster$time==148,5] / result.SA$df.TGT[148,5]
df.cluster[df.cluster$time==332,3] / result.SA$df.TGT[332,3]
df.cluster[df.cluster$time==332,5] / result.SA$df.TGT[332,5]
df.cluster[df.cluster$time==513,3] / result.SA$df.TGT[513,3]
df.cluster[df.cluster$time==513,5] / result.SA$df.TGT[513,5]

### ----------------------------------------------------------------------------
### Reproducibility 4: Numbers in Tables 2 and 3 in Section 5.2
### ----------------------------------------------------------------------------
source("Extends.R")
time.sel <- c(148, 330, 512)

fmt_num <- function(x, d) sprintf(paste0("%.", d, "f"), x)
fmt_ci  <- function(est, se, d) paste0("(", fmt_num(est-1.96*se,d), ", ", fmt_num(est+1.96*se,d), ")")
fmt_p   <- function(p) sprintf("%.3f", p)

### Compute survival difference (SD), survival ratio (SR), and RMST summaries
### from a fitted TrtSurvCurves object.
### ----------------------------------------------------------------------------
### Input:
###   res      : fitted object returned by TrtSurvCurves()
###   time.sel : evaluation days for SD and SR
###
### Output:
###   formatted SD/SR and RMST tables, together with the underlying estimates.
### ----------------------------------------------------------------------------
get_contrasts <- function(res, time.sel=c(148,330,512), s=1222) {
  
  zz <- file(tempfile(), open="wt"); sink(zz, type="message")
  
  ### Obtain SD, SR, and RMST estimates using FuseSurv_Extend().
  ext <- suppressWarnings(suppressMessages(FuseSurv_Extend(
    eval.times=res$eval.times, site=res$site, IF.00=res$IF.00, IF.01=res$IF.01,
    S.00=res$S.00, S.01=res$S.01, Aug.00.mean=res$Aug.00.mean, Aug.01.mean=res$Aug.01.mean,
    Aug.R0.mean=res$Aug.R0.mean, Aug.R1.mean=res$Aug.R1.mean,
    Aug.R0.mean.sour=res$Aug.R0.mean.sour, Aug.R1.mean.sour=res$Aug.R1.mean.sour,
    IF.R0=res$IF.R0, IF.R1=res$IF.R1, IF.CCOD.0=res$IF.CCOD.0, IF.CCOD.1=res$IF.CCOD.1,
    ind.R1.ccod=res$ind.R1.ccod, s=s)))
  sink(type="message"); close(zz)
  
  methods <- c("TGT","FED","CCOD")
  
  SD <- lapply(methods, function(m) ext[[paste0("df.RD.",m)]][ext[[paste0("df.RD.",m)]]$time %in% time.sel,])
  SR <- lapply(methods, function(m) ext[[paste0("df.SR.",m)]][ext[[paste0("df.SR.",m)]]$time %in% time.sel,])
  names(SD) <- names(SR) <- methods
  
  ### Compute two-sided Wald p-values for the reported causal contrasts.
  for(m in methods) {
    SD[[m]]$pval <- 2*pnorm(-abs(SD[[m]]$RD/SD[[m]]$sd))
    SR[[m]]$pval <- 2*pnorm(-abs((SR[[m]]$SR-1)/SR[[m]]$sd))
    ext[[paste0("df.RMST.diff.",m)]]$pval <-
      2*pnorm(-abs(ext[[paste0("df.RMST.diff.",m)]]$RMST/ext[[paste0("df.RMST.diff.",m)]]$sd))
  }
  
  ### Format the estimates, standard errors, confidence intervals, and p-values
  ### for presentation in Tables 2 and 3.
  tab.SD.SR <- do.call(rbind, lapply(time.sel, function(day) do.call(rbind, lapply(methods, function(m) {
    sd <- SD[[m]][SD[[m]]$time==day,]; sr <- SR[[m]][SR[[m]]$time==day,]
    data.frame(
      Day=day, Method=m,
      `SD Est. (95% CI)`=paste0(fmt_num(sd$RD,3)," ",fmt_ci(sd$RD,sd$sd,3)),
      `SE(SD)`=fmt_num(sd$sd,3), `p-value_SD`=fmt_p(sd$pval),
      `SR Est. (95% CI)`=paste0(fmt_num(sr$SR,3)," ",fmt_ci(sr$SR,sr$sd,3)),
      `SE(SR)`=fmt_num(sr$sd,3), `p-value_SR`=fmt_p(sr$pval), check.names=FALSE)
  }))))
  
  rmst_block <- function(label, prefix, show.p=FALSE) do.call(rbind, lapply(methods, function(m) {
    df <- ext[[paste0("df.RMST.",prefix,".",m)]]
    data.frame(Group=label, Method=m, `RMST Est.`=fmt_num(df$RMST,2), SE=fmt_num(df$sd,2),
               `95% CI`=fmt_ci(df$RMST,df$sd,2),
               `p-value`=if(show.p) fmt_p(df$pval) else "--", check.names=FALSE)
  }))
  
  tab.RMST <- rbind(rmst_block("Control group","0"), rmst_block("Treated group","1"),
                    rmst_block("RMST difference","diff",TRUE))
  
  list(ext=ext, SD.SR=tab.SD.SR, RMST=tab.RMST)
}

### --- SA target: Tables 2 and 3
load("result_main_SA.Rdata")
contrast.SA <- get_contrasts(result.SA)

tab.SD.SR <- contrast.SA$SD.SR
tab.RMST  <- contrast.SA$RMST

print(tab.SD.SR, row.names=FALSE)
print(tab.RMST, row.names=FALSE)

### ----------------------------------------------------------------------------
### Reproducibility 5: Figures A.1, A.2 and A.3 in Supplement A
### ----------------------------------------------------------------------------

### Repeat the treatment-specific survival analysis using OA, BP, and US
### successively as the target region.

### --- OA as the target region
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
plot_survival_CI(df=result.OA.naive$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw
plot_survival_CI(df=result.OA.naive$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool

### --- Below code prints panels (A) and (B) of Figure A.1
pdf(file="AMP_OA.pdf", width=12, height=6)
grid.arrange(p.fed, p.tgt, p.ccod, p.ivw, p.pool, p.cluster, ncol=3)
dev.off()

p.wt.0 <- plot_fedweights(result.OA, site.names=c("OA", "SA", "BP", "US"))$p.wt.0
p.wt.1 <- plot_fedweights(result.OA, site.names=c("OA", "SA", "BP", "US"))$p.wt.1

pdf(file="AMP_OAwts.pdf", width=10, height=3.4)
grid.arrange(p.wt.1, p.wt.0, ncol=2)
dev.off()

### --- BP as the target region
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
plot_survival_CI(df=result.BP.naive$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw
plot_survival_CI(df=result.BP.naive$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool

### --- Below code prints panels (A) and (B) of Figure A.2
pdf(file="AMP_BP.pdf", width=12, height=6)
grid.arrange(p.fed, p.tgt, p.ccod, p.ivw, p.pool, p.cluster, ncol=3)
dev.off()

p.wt.0 <- plot_fedweights(result.BP, site.names=c("BP", "OA", "SA", "US"))$p.wt.0
p.wt.1 <- plot_fedweights(result.BP, site.names=c("BP", "OA", "SA", "US"))$p.wt.1

pdf(file="AMP_BPwts.pdf", width=10, height=3.4)
grid.arrange(p.wt.1, p.wt.0, ncol=2)
dev.off()


### --- US as the target region
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
plot_survival_CI(df=result.US.naive$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw
plot_survival_CI(df=result.US.naive$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool

### --- Below code prints panels (A) and (B) of Figure A.3
pdf(file="AMP_US.pdf", width=12, height=6)
grid.arrange(p.fed, p.tgt, p.ccod, p.ivw, p.pool, p.cluster, ncol=3)
dev.off()

p.wt.0 <- plot_fedweights(result.US, site.names=c("US", "OA", "SA", "BP"))$p.wt.0
p.wt.1 <- plot_fedweights(result.US, site.names=c("US", "OA", "SA", "BP"))$p.wt.1

pdf(file="AMP_USwts.pdf", width=10, height=3.4)
grid.arrange(p.wt.1, p.wt.0, ncol=2)
dev.off()

### ----------------------------------------------------------------------------
### Reproducibility 6: Figure 1 in Section 2
### ----------------------------------------------------------------------------

### Combine the target-only survival curves from all four analytical regions
### to reproduce Figure 1.

load("result_main_SA.Rdata")
plot_survival_CI(df=result.SA$df.TGT, 
                 color_treated="springgreen4", 
                 color_control="springgreen", 
                 fig.title="SA (women)") -> p.tgt.SA

load("result_main_OA.Rdata")
plot_survival_CI(df=result.OA$df.TGT, 
                 color_treated="springgreen4", 
                 color_control="springgreen", 
                 fig.title="OA (women)") -> p.tgt.OA

load("result_main_BP.Rdata")
plot_survival_CI(df=result.BP$df.TGT, 
                 color_treated="springgreen4", 
                 color_control="springgreen", 
                 fig.title="BP (men, TG)") -> p.tgt.BP

load("result_main_US.Rdata")
plot_survival_CI(df=result.US$df.TGT, 
                 color_treated="springgreen4", 
                 color_control="springgreen", 
                 fig.title="US (men, TG)") -> p.tgt.US

pdf("AMP_siteSurvs.pdf", width=9, height=4.5)
grid.arrange(p.tgt.SA, p.tgt.OA, p.tgt.BP, p.tgt.US, ncol=2)
dev.off()


### ----------------------------------------------------------------------------
### Reproducibility 7: Causal contrasts for OA, BP, and US target regions
### --- Tables A.1--A.6 in Web Appendix A
### ----------------------------------------------------------------------------

### Apply the same contrast function used for SA to each additional target region
### to reproduce Tables A.1--A.6.

load("result_main_OA.Rdata")
contrast.OA <- get_contrasts(result.OA)
load("result_main_BP.Rdata")
contrast.BP <- get_contrasts(result.BP)
load("result_main_US.Rdata")
contrast.US <- get_contrasts(result.US)

### --- Tables A.1 and A.2: OA
print(contrast.OA$SD.SR, row.names=FALSE)
print(contrast.OA$RMST, row.names=FALSE)

### --- Tables A.3 and A.4: BP
print(contrast.BP$SD.SR, row.names=FALSE)
print(contrast.BP$RMST, row.names=FALSE)

### --- Tables A.5 and A.6: US
print(contrast.US$SD.SR, row.names=FALSE)
print(contrast.US$RMST, row.names=FALSE)
