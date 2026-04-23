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
library(cobalt)
library(randomForestSRC)
library(cowplot)

## =========================================================
## 0. Read data and construct analysis dataset
## =========================================================
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
site[site == "South Africa"] <- "SA"
site[site %in% c("Tanzania, Mozambique, Kenya", "Zimbabwe", "Botswana", "Malawi")] <- "OA"
site[site %in% c("Peru", "Brazil")] <- "BP"
site[site %in% c("United States", "Switzerland")] <- "US"

site <- factor(
  site,
  levels = c("SA", "OA", "BP", "US")
)

dat.hiv <- data.frame(
  A = A,
  Y = Y,
  Delta = Delta,
  bweight = X$bweight,
  score = X$score,
  age = X$age,
  site = site
)

## basic setup
bal.vars <- c("age", "score", "bweight")
X.ps <- dat.hiv[, bal.vars, drop = FALSE]
tgt.name <- "SA"

site.cols <- c(
  "SA" = "#E41A1C",
  "OA" = "#377EB8",
  "BP" = "#4DAF4A",
  "US" = "#984EA3"
)

## =========================================================
## 1. Helper functions
## =========================================================

## helper: truncate PS and build ATE weights
make_ate_wt <- function(ps, A, eps = 0.01) {
  ps <- pmin(pmax(ps, eps), 1 - eps)
  wt <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
  return(list(ps = ps, wt = wt))
}

## interaction learner for SuperLearner
SL.glm.interaction <- function(Y, X, newX, family, obsWeights, ...) {
  fit <- glm(
    Y ~ (age + score + bweight)^2,
    data = data.frame(Y = Y, X),
    family = family,
    weights = obsWeights
  )
  pred <- predict(fit, newdata = newX, type = "response")
  fit.out <- list(object = fit)
  class(fit.out) <- "SL.glm.interaction"
  out <- list(pred = pred, fit = fit.out)
  return(out)
}

predict.SL.glm.interaction <- function(object, newdata, ...) {
  predict(object$object, newdata = newdata, type = "response")
}

## fit one PS model within each site
fit_sitewise_ps <- function(dat, method = c("GLM", "GLM.interaction", "LASSO", "Ensemble.All", "Ensemble.GLM"),
                            bal.vars = c("age", "score", "bweight"), eps = 0.01) {
  method <- match.arg(method)
  ps <- rep(NA_real_, nrow(dat))
  wt <- rep(NA_real_, nrow(dat))
  
  for (s in levels(dat$site)) {
    idx <- dat$site == s
    dat.s <- dat[idx, , drop = FALSE]
    
    if (method == "GLM") {
      fit <- glm(A ~ age + score + bweight, data = dat.s, family = binomial())
      ps.s <- predict(fit, type = "response")
    }
    
    if (method == "GLM.interaction") {
      fit <- glm(A ~ (age + score + bweight)^2, data = dat.s, family = binomial())
      ps.s <- predict(fit, type = "response")
    }
    
    if (method == "LASSO") {
      x.mat <- model.matrix(
        A ~ age + score + bweight + age:score + age:bweight + score:bweight,
        data = dat.s
      )[, -1, drop = FALSE]
      
      cvfit <- cv.glmnet(
        x = x.mat,
        y = dat.s$A,
        family = "binomial",
        alpha = 1,
        nfolds = min(10, nrow(dat.s))
      )
      ps.s <- as.numeric(
        predict(cvfit, newx = x.mat, s = "lambda.min", type = "response")
      )
    }
    
    if (method == "Ensemble.All") {
      fit <- SuperLearner(
        Y = dat.s$A,
        X = dat.s[, bal.vars, drop = FALSE],
        family = binomial(),
        SL.library = c("SL.glm", "SL.glm.interaction", "SL.glmnet"),
        method = "method.NNLS"
      )
      ps.s <- fit$SL.predict
    }
    
    if (method == "Ensemble.GLM") {
      fit <- SuperLearner(
        Y = dat.s$A,
        X = dat.s[, bal.vars, drop = FALSE],
        family = binomial(),
        SL.library = c("SL.glm", "SL.glm.interaction"),
        method = "method.NNLS"
      )
      ps.s <- fit$SL.predict
    }
    
    tmp <- make_ate_wt(ps.s, dat.s$A, eps = eps)
    ps[idx] <- tmp$ps
    wt[idx] <- tmp$wt
  }
  
  return(list(ps = ps, wt = wt))
}

## =========================================================
## 2. Site-specific PS models
## =========================================================

ps.glm.out      <- fit_sitewise_ps(dat.hiv, method = "GLM",             bal.vars = bal.vars)
ps.glm.int.out  <- fit_sitewise_ps(dat.hiv, method = "GLM.interaction", bal.vars = bal.vars)
ps.lasso.out    <- fit_sitewise_ps(dat.hiv, method = "LASSO",           bal.vars = bal.vars)
ps.ens.out      <- fit_sitewise_ps(dat.hiv, method = "Ensemble.All",    bal.vars = bal.vars)
ps.ens.small.out<- fit_sitewise_ps(dat.hiv, method = "Ensemble.GLM",    bal.vars = bal.vars)

ps.glm       <- ps.glm.out$ps
wt.glm       <- ps.glm.out$wt
ps.glm.int   <- ps.glm.int.out$ps
wt.glm.int   <- ps.glm.int.out$wt
ps.lasso     <- ps.lasso.out$ps
wt.lasso     <- ps.lasso.out$wt
ps.ens       <- ps.ens.out$ps
wt.ens       <- ps.ens.out$wt
ps.ens.small <- ps.ens.small.out$ps
wt.ens.small <- ps.ens.small.out$wt

## =========================================================
## 3. Love plot based on site-specific weights
## =========================================================

wt.list <- data.frame(
  GLM = wt.glm,
  GLM.interaction = wt.glm.int,
  LASSO = wt.lasso,
  Ensemble.All = wt.ens,
  Ensemble.GLM = wt.ens.small
)

bal.all <- bal.tab(
  x = X.ps,
  treat = dat.hiv$A,
  weights = wt.list,
  method = "weighting",
  estimand = "ATE",
  un = TRUE
)

p.all <- love.plot(
  bal.all,
  stats = "mean.diffs",
  abs = TRUE,
  threshold = 0.03,
  sample.names = c("Unadjusted", "GLM", "GLM.int", "LASSO", "Ens.All", "Ens.GLM"),
  title = "(A) Covariate balance love plot",
  colors = c("grey50", "#1B9E77", "#D95F02", "blue", "#E7298A", "#66A61E"),
  shapes = c(15, 17, 15, 17, 16, 17),
  alpha = 0.65,
  line = FALSE,
  size = 3
) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

print(p.all)

## =========================================================
## 4. Site-specific PS histograms
## =========================================================

ps.site.df <- data.frame(
  site = dat.hiv$site,
  A = factor(dat.hiv$A, levels = c(0, 1), labels = c("Control", "Treatment")),
  GLM = ps.glm,
  GLM.interaction = ps.glm.int,
  LASSO = ps.lasso,
  Ensemble.All = ps.ens,
  Ensemble.GLM = ps.ens.small
)

ps.site.long <- ps.site.df %>%
  pivot_longer(
    cols = c(GLM, GLM.interaction, LASSO, Ensemble.All, Ensemble.GLM),
    names_to = "Method",
    values_to = "PS"
  )

ps.site.long$Method <- factor(
  ps.site.long$Method,
  levels = c("GLM", "GLM.interaction", "LASSO", "Ensemble.All", "Ensemble.GLM")
)

p.ps.hist.site <- ggplot(ps.site.long, aes(x = PS, fill = A)) +
  geom_histogram(
    bins = 20,
    position = "identity",
    alpha = 0.45,
    color = "black",
    linewidth = 0.25
  ) +
  facet_grid(site ~ Method) +
  scale_fill_manual(values = c("Control" = "#E41A1C", "Treatment" = "#377EB8")) +
  labs(
    title = "(B) Site-specific propensity score histograms",
    x = "Estimated propensity score",
    y = "Count",
    fill = "Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom",
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

pdf("AMP_pslove.pdf", height = 3, width = 10)
grid.arrange(
  p.all, p.ps.hist.site,
  ncol = 2,
  widths = c(1,2)
)
dev.off()

## =========================================================
## 5. Positivity diagnostics for nuisance functions
##    survival / censoring / density ratio
## =========================================================

## NOTE:
## This part assumes estimate_omega_np() is available in your environment,
## as in your main analysis code.

fit.times <- 1:601
t.check   <- 600
t.ind     <- match(t.check, fit.times)
n.folds   <- 5

dat.use <- dat.hiv %>%
  mutate(id = 1:n())

## ---------- full target-site survival models ----------
dat0 <- dat.use %>% filter(site == tgt.name)

surv.fit.0.tgt <- survSuperLearner(
  time = dat0$Y[dat0$A == 0],
  event = dat0$Delta[dat0$A == 0],
  X = dat0[dat0$A == 0, bal.vars, drop = FALSE],
  new.times = fit.times,
  event.SL.library = c("survSL.km", "survSL.coxph", "survSL.rfsrc"),
  cens.SL.library  = c("survSL.km", "survSL.coxph", "survSL.rfsrc")
)

surv.fit.1.tgt <- survSuperLearner(
  time = dat0$Y[dat0$A == 1],
  event = dat0$Delta[dat0$A == 1],
  X = dat0[dat0$A == 1, bal.vars, drop = FALSE],
  new.times = fit.times,
  event.SL.library = c("survSL.km", "survSL.coxph", "survSL.rfsrc"),
  cens.SL.library  = c("survSL.km", "survSL.coxph", "survSL.rfsrc")
)

## ---------- collect out-of-fold nuisance predictions ----------
set.seed(12345)
nuis.list <- list()

## target site
folds.tgt <- createFolds(1:nrow(dat0), k = n.folds, list = TRUE)

for (i in seq_along(folds.tgt)) {
  pred.ind  <- folds.tgt[[i]]
  train.ind <- setdiff(seq_len(nrow(dat0)), pred.ind)
  A.train   <- dat0$A[train.ind]
  X.pred    <- dat0[pred.ind, bal.vars, drop = FALSE]
  
  fit0 <- survSuperLearner(
    time = dat0$Y[train.ind][A.train == 0],
    event = dat0$Delta[train.ind][A.train == 0],
    X = dat0[train.ind, bal.vars, drop = FALSE][A.train == 0, , drop = FALSE],
    new.times = fit.times,
    event.SL.library = c("survSL.km", "survSL.coxph", "survSL.rfsrc"),
    cens.SL.library  = c("survSL.km", "survSL.coxph", "survSL.rfsrc")
  )
  
  fit1 <- survSuperLearner(
    time = dat0$Y[train.ind][A.train == 1],
    event = dat0$Delta[train.ind][A.train == 1],
    X = dat0[train.ind, bal.vars, drop = FALSE][A.train == 1, , drop = FALSE],
    new.times = fit.times,
    event.SL.library = c("survSL.km", "survSL.coxph", "survSL.rfsrc"),
    cens.SL.library  = c("survSL.km", "survSL.coxph", "survSL.rfsrc")
  )
  
  pred0 <- predict.survSuperLearner(fit0, newdata = X.pred, new.times = fit.times)
  pred1 <- predict.survSuperLearner(fit1, newdata = X.pred, new.times = fit.times)
  
  nuis.list[[length(nuis.list) + 1]] <- data.frame(
    id     = dat0$id[pred.ind],
    site   = dat0$site[pred.ind],
    A      = dat0$A[pred.ind],
    S0_end = pred0$event.SL.predict[, t.ind],
    S1_end = pred1$event.SL.predict[, t.ind],
    G0_end = pred0$cens.SL.predict[, t.ind],
    G1_end = pred1$cens.SL.predict[, t.ind],
    omega  = 1
  )
}

## source sites
source.sites <- setdiff(levels(dat.use$site), tgt.name)
X0 <- as.matrix(dat0[, bal.vars, drop = FALSE])

for (s in source.sites) {
  dat.r <- dat.use %>% filter(site == s)
  folds.r <- createFolds(1:nrow(dat.r), k = n.folds, list = TRUE)
  
  for (i in seq_along(folds.r)) {
    pred.ind  <- folds.r[[i]]
    train.ind <- setdiff(seq_len(nrow(dat.r)), pred.ind)
    A.train   <- dat.r$A[train.ind]
    X.pred    <- dat.r[pred.ind, bal.vars, drop = FALSE]
    
    ## density ratio
    omega.hat <- estimate_omega_np(
      x = as.matrix(dat.r[train.ind, bal.vars, drop = FALSE]),
      x_target = X0,
      x.pred = as.matrix(X.pred),
      method = "logistic"
    )
    
    ## target survival model
    predS0 <- predict.survSuperLearner(surv.fit.0.tgt, newdata = X.pred, new.times = fit.times)
    predS1 <- predict.survSuperLearner(surv.fit.1.tgt, newdata = X.pred, new.times = fit.times)
    
    ## source censoring model
    fit0.src <- survSuperLearner(
      time = dat.r$Y[train.ind][A.train == 0],
      event = dat.r$Delta[train.ind][A.train == 0],
      X = dat.r[train.ind, bal.vars, drop = FALSE][A.train == 0, , drop = FALSE],
      new.times = fit.times,
      event.SL.library = c("survSL.km", "survSL.coxph", "survSL.rfsrc"),
      cens.SL.library  = c("survSL.km", "survSL.coxph", "survSL.rfsrc")
    )
    
    fit1.src <- survSuperLearner(
      time = dat.r$Y[train.ind][A.train == 1],
      event = dat.r$Delta[train.ind][A.train == 1],
      X = dat.r[train.ind, bal.vars, drop = FALSE][A.train == 1, , drop = FALSE],
      new.times = fit.times,
      event.SL.library = c("survSL.km", "survSL.coxph", "survSL.rfsrc"),
      cens.SL.library  = c("survSL.km", "survSL.coxph", "survSL.rfsrc")
    )
    
    predG0 <- predict.survSuperLearner(fit0.src, newdata = X.pred, new.times = fit.times)
    predG1 <- predict.survSuperLearner(fit1.src, newdata = X.pred, new.times = fit.times)
    
    nuis.list[[length(nuis.list) + 1]] <- data.frame(
      id     = dat.r$id[pred.ind],
      site   = dat.r$site[pred.ind],
      A      = dat.r$A[pred.ind],
      S0_end = predS0$event.SL.predict[, t.ind],
      S1_end = predS1$event.SL.predict[, t.ind],
      G0_end = predG0$cens.SL.predict[, t.ind],
      G1_end = predG1$cens.SL.predict[, t.ind],
      omega  = as.numeric(omega.hat)
    )
  }
}

nuis.df <- bind_rows(nuis.list)

## summary table
nuis.summary <- nuis.df %>%
  group_by(site) %>%
  summarise(
    n = n(),
    S0_min = min(S0_end, na.rm = TRUE),
    S0_q01 = quantile(S0_end, 0.01, na.rm = TRUE),
    S1_min = min(S1_end, na.rm = TRUE),
    S1_q01 = quantile(S1_end, 0.01, na.rm = TRUE),
    G0_min = min(G0_end, na.rm = TRUE),
    G0_q01 = quantile(G0_end, 0.01, na.rm = TRUE),
    G1_min = min(G1_end, na.rm = TRUE),
    G1_q01 = quantile(G1_end, 0.01, na.rm = TRUE),
    omega_max = max(omega, na.rm = TRUE),
    omega_q99 = quantile(omega, 0.99, na.rm = TRUE)
  )

print(nuis.summary)

## =========================================================
## 5A-5C. Compact nuisance positivity plots
## Layout: A | B | C
## A and B each contain two vertically stacked subplots
## =========================================================

## ---------- survival ----------
surv.A0 <- ggplot(nuis.df, aes(x = S0_end, fill = site)) +
  geom_histogram(
    bins = 30,
    alpha = 0.55,
    color = "black",
    linewidth = 0.25,
    position = "identity"
  ) +
  scale_fill_manual(values = site.cols) +
  labs(
    title = paste0("Predicted S(t | A=0, X) at day ", t.check),
    x = "Predicted survival probability",
    y = "Count",
    fill = "Site"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

surv.A1 <- ggplot(nuis.df, aes(x = S1_end, fill = site)) +
  geom_histogram(
    bins = 30,
    alpha = 0.55,
    color = "black",
    linewidth = 0.25,
    position = "identity"
  ) +
  scale_fill_manual(values = site.cols) +
  labs(
    title = paste0("Predicted S(t | A=1, X) at day ", t.check),
    x = "Predicted survival probability",
    y = "Count",
    fill = "Site"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

## ---------- censoring ----------
cens.A0 <- ggplot(nuis.df, aes(x = G0_end, fill = site)) +
  geom_histogram(
    bins = 30,
    alpha = 0.55,
    color = "black",
    linewidth = 0.25,
    position = "identity"
  ) +
  scale_fill_manual(values = site.cols) +
  labs(
    title = paste0("Predicted G(t | A=0, X) at day ", t.check),
    x = "Predicted censoring survival probability",
    y = "Count",
    fill = "Site"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    # legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

cens.A1 <- ggplot(nuis.df, aes(x = G1_end, fill = site)) +
  geom_histogram(
    bins = 30,
    alpha = 0.55,
    color = "black",
    linewidth = 0.25,
    position = "identity"
  ) +
  scale_fill_manual(values = site.cols) +
  labs(
    title = paste0("Predicted G(t | A=1, X) at day ", t.check),
    x = "Predicted censoring survival probability",
    y = "Count",
    fill = "Site"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    # legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

## ---------- density ratio ----------
omega.df <- nuis.df %>%
  filter(site != tgt.name)

p.omega <- ggplot(omega.df, aes(x = omega, fill = site)) +
  geom_histogram(
    bins = 30,
    alpha = 0.60,
    color = "black",
    linewidth = 0.25,
    position = "identity"
  ) +
  facet_wrap(~ site, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = site.cols[names(site.cols) != tgt.name]) +
  labs(
    title = "Individual density-ratio estimates",
    x = "Estimated density-ratio from source to target",
    y = "Count",
    fill = "Site"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.position = "bottom",
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

## ---------- shared legend with all regions ----------
legend.df <- data.frame(
  x = 1:4,
  y = 1,
  site = factor(c("SA", "OA", "BP", "US"), levels = c("SA", "OA", "BP", "US"))
)

p.legend <- ggplot(legend.df, aes(x = x, y = y, fill = site)) +
  geom_col() +
  scale_fill_manual(
    values = site.cols,
    breaks = c("SA", "OA", "BP", "US"),
    labels = c("SA", "OA", "BP", "US"),
    name = "Region"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )

legend.shared <- cowplot::get_legend(p.legend)

## remove legends from main panels
p.omega.noleg <- p.omega + theme(legend.position = "none")

## ---------- assemble A, B, C ----------
panel.A <- plot_grid(
  surv.A0, surv.A1,
  ncol = 1,
  labels = c("A", ""),
  label_size = 13,
  rel_heights = c(1, 1)
)

panel.B <- plot_grid(
  cens.A0, cens.A1,
  ncol = 1,
  labels = c("B", ""),
  label_size = 13,
  rel_heights = c(1, 1)
)

panel.C <- plot_grid(
  p.omega.noleg,
  ncol = 1,
  labels = c("C"),
  label_size = 13
)

p.nuis <- plot_grid(
  panel.A, panel.B, panel.C,
  nrow = 1,
  rel_widths = c(0.8, 1, 0.65),
  align = "h"
)

p.nuis.final <- plot_grid(
  p.nuis,
  legend.shared,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

print(p.nuis.final)

pdf("AMP_SComega.pdf", width = 14, height = 4.5)
print(p.nuis.final)
dev.off()