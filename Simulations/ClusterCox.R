library(survival)
library(dplyr)

fit_cluster_cox_one <- function(dat, eval.times, target_site = 0) {
  dat$siteF <- factor(dat$site)
  
  fit <- coxph(
    Surv(Y, Delta) ~ A + X1 + X2 + X3 + strata(siteF),
    data = dat,
    cluster = site,
    ties = "breslow"
  )
  dat.tgt <- dat[dat$site == target_site, , drop = FALSE]
  newdata <- data.frame(
    A = c(1, 0),
    X1 = mean(dat.tgt$X1),
    X2 = mean(dat.tgt$X2),
    X3 = mean(dat.tgt$X3),
    siteF = factor(target_site, levels = levels(dat$siteF))
  )
  sf <- survfit(fit, newdata = newdata, se.fit = TRUE)
  ss <- summary(sf, times = eval.times, extend = TRUE)
  out <- data.frame(
    time = ss$time,
    surv = ss$surv,
    std.err = ss$std.err,
    strata = rep(1:2, each = length(eval.times))
  )
  df.CLCOX <- data.frame(
    time = eval.times,
    surv1 = out$surv[out$strata == 1],
    surv1.sd = out$std.err[out$strata == 1],
    surv0 = out$surv[out$strata == 2],
    surv0.sd = out$std.err[out$strata == 2]
  )
  list(df.CLCOX = df.CLCOX)
}

run_cluster_cox_sim <- function(datlist, eval.times, target_site = 0) {
  M <- length(datlist)
  out <- vector("list", M)
  
  for (j in seq_len(M)) {
    cat("replicate", j, "of", M, "\n")
    out[[j]] <- fit_cluster_cox_one(
      dat = datlist[[j]],
      eval.times = eval.times,
      target_site = target_site
    )
  }
  
  out
}

time.RE <- c(30, 60, 90)

load("obsdata_s.Rdata")
result.CLCOX.homo <- run_cluster_cox_sim(dat.homo, eval.times = time.RE, target_site = 0)
result.CLCOX.diffX <- run_cluster_cox_sim(dat.diffX, eval.times = time.RE, target_site = 0)
result.CLCOX.diffT <- run_cluster_cox_sim(dat.diffT, eval.times = time.RE, target_site = 0)
result.CLCOX.diffC <- run_cluster_cox_sim(dat.diffC, eval.times = time.RE, target_site = 0)
result.CLCOX.diffAll <- run_cluster_cox_sim(dat.diffAll, eval.times = time.RE, target_site = 0)

save(
  result.CLCOX.homo,
  result.CLCOX.diffX,
  result.CLCOX.diffT,
  result.CLCOX.diffC,
  result.CLCOX.diffAll,
  file = "Res_CLCOX_s.Rdata"
)

load("obsdata_l.Rdata")
result.CLCOX.homo <- run_cluster_cox_sim(dat.homo, eval.times = time.RE, target_site = 0)
result.CLCOX.diffX <- run_cluster_cox_sim(dat.diffX, eval.times = time.RE, target_site = 0)
result.CLCOX.diffT <- run_cluster_cox_sim(dat.diffT, eval.times = time.RE, target_site = 0)
result.CLCOX.diffC <- run_cluster_cox_sim(dat.diffC, eval.times = time.RE, target_site = 0)
result.CLCOX.diffAll <- run_cluster_cox_sim(dat.diffAll, eval.times = time.RE, target_site = 0)

save(
  result.CLCOX.homo,
  result.CLCOX.diffX,
  result.CLCOX.diffT,
  result.CLCOX.diffC,
  result.CLCOX.diffAll,
  file = "Res_CLCOX_l.Rdata"
)

load("obsdata_l2.Rdata")
result.CLCOX.homo <- run_cluster_cox_sim(dat.homo, eval.times = time.RE, target_site = 0)
result.CLCOX.diffX <- run_cluster_cox_sim(dat.diffX, eval.times = time.RE, target_site = 0)
result.CLCOX.diffT <- run_cluster_cox_sim(dat.diffT, eval.times = time.RE, target_site = 0)
result.CLCOX.diffC <- run_cluster_cox_sim(dat.diffC, eval.times = time.RE, target_site = 0)
result.CLCOX.diffAll <- run_cluster_cox_sim(dat.diffAll, eval.times = time.RE, target_site = 0)

save(
  result.CLCOX.homo,
  result.CLCOX.diffX,
  result.CLCOX.diffT,
  result.CLCOX.diffC,
  result.CLCOX.diffAll,
  file = "Res_CLCOX_l2.Rdata"
)


load("obsdata_limO.Rdata")
result.CLCOX.homo <- run_cluster_cox_sim(dat.homo, eval.times = time.RE, target_site = 0)
result.CLCOX.diffX <- run_cluster_cox_sim(dat.diffX, eval.times = time.RE, target_site = 0)
result.CLCOX.diffT <- run_cluster_cox_sim(dat.diffT, eval.times = time.RE, target_site = 0)
result.CLCOX.diffC <- run_cluster_cox_sim(dat.diffC, eval.times = time.RE, target_site = 0)
result.CLCOX.diffAll <- run_cluster_cox_sim(dat.diffAll, eval.times = time.RE, target_site = 0)

save(
  result.CLCOX.homo,
  result.CLCOX.diffX,
  result.CLCOX.diffT,
  result.CLCOX.diffC,
  result.CLCOX.diffAll,
  file = "Res_CLCOX_limO.Rdata"
)

load("obsdata_imbal.Rdata")
result.CLCOX.imbal <- run_cluster_cox_sim(dat.imbal, eval.times = time.RE, target_site = 0)
result.CLCOX.moreK <- run_cluster_cox_sim(dat.moreK, eval.times = time.RE, target_site = 0)

save(
  result.CLCOX.imbal,
  result.CLCOX.moreK,
  file = "Res_CLCOX_imbal.Rdata"
)
