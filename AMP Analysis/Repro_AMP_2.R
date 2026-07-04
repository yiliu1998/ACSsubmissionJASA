### ----------------------------------------------------------------------------
### Reproducibility for AMP Data Analysis: Part II
### --- This file reproduces Figures A.4, Tables A.1 and A.2 in Online
### --- Supplemental Material A
### ----------------------------------------------------------------------------

### Packages and data pre-processing
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
### Reproducibility 1: Figure A.4 in Supplement A
### ----------------------------------------------------------------------------
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

### --- Run and save the results
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

load("result_main_SA_2.Rdata")
plot_survival_CI(df=result.SA$df.TGT, color_treated="springgreen4", color_control="springgreen", fig.title="TGT") -> p.tgt
plot_survival_CI(df=result.SA$df.CCOD, color_treated="darkslateblue", color_control="darkorchid1", fig.title="CCOD") -> p.ccod
plot_survival_CI(df=result.SA$df.FED, color_treated="cyan4", color_control="cyan", fig.title="FED") -> p.fed
plot_survival_CI(df = df.cluster, color_treated = "plum4", color_control = "maroon2", fig.title = "CLCOX") -> p.cluster
plot_survival_CI(df=result.SA.naive$df.IVW, color_treated="darkorange3", color_control="orange", fig.title="IVW") -> p.ivw
plot_survival_CI(df=result.SA.naive$df.POOL, color_treated="firebrick4", color_control="tomato1", fig.title="POOL") -> p.pool

### --- Output panels (A) and (B) of Figure A.4
pdf(file="AMP_SA_2.pdf", width=12, height=6)
grid.arrange(p.fed, p.tgt, p.ccod, p.ivw, p.pool, p.cluster, ncol=3)
dev.off()

p.wt.SA0 <- plot_fedweights(result.SA)$p.wt.0
p.wt.SA1 <- plot_fedweights(result.SA)$p.wt.1

pdf(file="AMP_SAwts_2.pdf", width=10, height=3.4)
grid.arrange(p.wt.SA1, p.wt.SA0, ncol=2)
dev.off()

### ----------------------------------------------------------------------------
### Reproducibility 2: Tables A.1 and A.2 in Supplement A
### ----------------------------------------------------------------------------
### --- Run results
source("Extends.R")
load("result_main_SA_2.Rdata")
extend.results <- NULL
zz <- file(tempfile(), open = "wt")
sink(zz, type = "message")
invisible(capture.output({
  extend.results <- suppressWarnings(suppressMessages(
    FuseSurv_Extend(eval.times=result.SA$eval.times,
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
  ))
}, type = "output"))
sink(type = "message")
close(zz)

time.sel <- c(148, 330, 512)
df.RD.TGT  <- extend.results$df.RD.TGT [extend.results$df.RD.TGT$time %in% time.sel, ]
df.RD.FED  <- extend.results$df.RD.FED [extend.results$df.RD.FED$time %in% time.sel, ]
df.RD.CCOD <- extend.results$df.RD.CCOD[extend.results$df.RD.CCOD$time %in% time.sel, ]

df.SR.TGT  <- extend.results$df.SR.TGT [extend.results$df.SR.TGT$time %in% time.sel, ]
df.SR.FED  <- extend.results$df.SR.FED [extend.results$df.SR.FED$time %in% time.sel, ]
df.SR.CCOD <- extend.results$df.SR.CCOD[extend.results$df.SR.CCOD$time %in% time.sel, ]

df.RD.TGT$pval  <- 2 * pnorm(-abs(df.RD.TGT$RD  / df.RD.TGT$sd))
df.RD.FED$pval  <- 2 * pnorm(-abs(df.RD.FED$RD  / df.RD.FED$sd))
df.RD.CCOD$pval <- 2 * pnorm(-abs(df.RD.CCOD$RD / df.RD.CCOD$sd))

df.SR.TGT$pval  <- 2 * pnorm(-abs((df.SR.TGT$SR  - 1) / df.SR.TGT$sd))
df.SR.FED$pval  <- 2 * pnorm(-abs((df.SR.FED$SR  - 1) / df.SR.FED$sd))
df.SR.CCOD$pval <- 2 * pnorm(-abs((df.SR.CCOD$SR - 1) / df.SR.CCOD$sd))

extend.results$df.RMST.0.TGT$pval  <- NA
extend.results$df.RMST.0.FED$pval  <- NA
extend.results$df.RMST.0.CCOD$pval <- NA

extend.results$df.RMST.1.TGT$pval  <- NA
extend.results$df.RMST.1.FED$pval  <- NA
extend.results$df.RMST.1.CCOD$pval <- NA

extend.results$df.RMST.diff.TGT$pval  <- 2 * pnorm(-abs(extend.results$df.RMST.diff.TGT$RMST  / extend.results$df.RMST.diff.TGT$sd))
extend.results$df.RMST.diff.FED$pval  <- 2 * pnorm(-abs(extend.results$df.RMST.diff.FED$RMST  / extend.results$df.RMST.diff.FED$sd))
extend.results$df.RMST.diff.CCOD$pval <- 2 * pnorm(-abs(extend.results$df.RMST.diff.CCOD$RMST / extend.results$df.RMST.diff.CCOD$sd))

fmt_num <- function(x, digits) sprintf(paste0("%.", digits, "f"), x)
fmt_ci <- function(est, se, digits) paste0("(", fmt_num(est - 1.96 * se, digits), ", ", fmt_num(est + 1.96 * se, digits), ")")
fmt_p <- function(p) ifelse(is.na(p), "--", sprintf("%.3f", p))

make_rd_sr_rows <- function(rd.list, sr.list, time.sel=c(148,330,512)) {
  methods <- c("TGT", "FED", "CCOD")
  do.call(rbind, lapply(time.sel, function(day) {
    do.call(rbind, lapply(methods, function(m) {
      rd <- rd.list[[m]][rd.list[[m]]$time == day, ]
      sr <- sr.list[[m]][sr.list[[m]]$time == day, ]
      data.frame(
        Day = day, Method = m,
        `RD Est. (95% CI)` = paste0(fmt_num(rd$RD, 3), " ", fmt_ci(rd$RD, rd$sd, 3)),
        `SE(RD)` = fmt_num(rd$sd, 3),
        `p-value RD` = fmt_p(rd$pval),
        `SR Est. (95% CI)` = paste0(fmt_num(sr$SR, 3), " ", fmt_ci(sr$SR, sr$sd, 3)),
        `SE(SR)` = fmt_num(sr$sd, 3),
        `p-value SR` = fmt_p(sr$pval),
        check.names = FALSE
      )
    }))
  }))
}

tab.A1 <- make_rd_sr_rows(
  rd.list = list(TGT=df.RD.TGT, FED=df.RD.FED, CCOD=df.RD.CCOD),
  sr.list = list(TGT=df.SR.TGT, FED=df.SR.FED, CCOD=df.SR.CCOD),
  time.sel = time.sel
)

make_rmst_rows <- function(df0.list, df1.list, dfdiff.list) {
  methods <- c("TGT", "FED", "CCOD")
  bind_block <- function(label, dfs, show_p=FALSE) {
    do.call(rbind, lapply(methods, function(m) {
      df <- dfs[[m]]
      data.frame(
        Group = label, Method = m,
        `RMST Est.` = fmt_num(df$RMST, 2),
        SE = fmt_num(df$sd, 2),
        `95% CI` = fmt_ci(df$RMST, df$sd, 2),
        `p-value` = if (show_p) fmt_p(df$pval) else "--",
        check.names = FALSE
      )
    }))
  }
  rbind(
    bind_block("Control group", df0.list, show_p=FALSE),
    bind_block("Treated group", df1.list, show_p=FALSE),
    bind_block("RMST difference", dfdiff.list, show_p=TRUE)
  )
}

tab.A2 <- make_rmst_rows(
  df0.list = list(TGT=extend.results$df.RMST.0.TGT, FED=extend.results$df.RMST.0.FED, CCOD=extend.results$df.RMST.0.CCOD),
  df1.list = list(TGT=extend.results$df.RMST.1.TGT, FED=extend.results$df.RMST.1.FED, CCOD=extend.results$df.RMST.1.CCOD),
  dfdiff.list = list(TGT=extend.results$df.RMST.diff.TGT, FED=extend.results$df.RMST.diff.FED, CCOD=extend.results$df.RMST.diff.CCOD)
)

### --- Print Table A.1
print(tab.A1, row.names = FALSE)
# xtable(tab.A1, align = c("r","r","r","c","c","c","c","c","c"),
#        caption = "Estimated risk RD and SR at days 148, 330, and 512.",
#        label = "tab:RD-AMP-supp")

### --- Print Table A.2
print(tab.A2, row.names = FALSE)
# xtable(tab.A2, align = c("r","l","r","c","c","c","c"),
#        caption = "Estimated RMST by treatment group and RMST difference up to day 601.",
#        label = "tab:RMST-AMP-supp")