# ============================================================
# Clustered Cox benchmark for the simulation study
#
# This script fits a stratified Cox proportional hazards model
# with clustering by site to each simulated dataset and evaluates
# treatment-specific survival probabilities in the target site.
#
# The analysis is repeated across all Monte Carlo replicates and
# across the simulation scenarios and sample-size settings generated
# in the main data-generation script.
# ============================================================

library(survival)
library(dplyr)

# ------------------------------------------------------------
# Fit the clustered Cox model to one simulated dataset
#
# Inputs:
#   dat         : one simulated dataset containing site, treatment,
#                 covariates, observed time, and event indicator
#   eval.times  : vector of time points at which survival is evaluated
#   target_site : site used to define the target covariate profile
#
# Output:
#   a list containing df.CLCOX, which reports estimated survival
#   probabilities and standard errors under treatment and control
#   at the requested evaluation times
# ------------------------------------------------------------
fit_cluster_cox_one <- function(dat, eval.times, target_site = 0) {
  
  # Convert site indicator to a factor for stratification in the Cox model.
  dat$siteF <- factor(dat$site)
  
  # Fit a Cox model with a common treatment/covariate effect across sites,
  # site-specific baseline hazards through strata(), and cluster-robust
  # variance estimation by site.
  fit <- coxph(Surv(Y, Delta) ~ A + X1 + X2 + X3 + strata(siteF), data = dat, cluster = site, ties = "breslow")
  
  # Extract target-site observations to define the target covariate profile.
  dat.tgt <- dat[dat$site == target_site, , drop = FALSE]
  
  # Construct treatment and control prediction profiles at the mean target-site covariate values.
  newdata <- data.frame(A = c(1, 0), X1 = mean(dat.tgt$X1), X2 = mean(dat.tgt$X2), X3 = mean(dat.tgt$X3), siteF = factor(target_site, levels = levels(dat$siteF)))
  
  # Obtain treatment-specific survival curves and associated standard errors.
  sf <- survfit(fit, newdata = newdata, se.fit = TRUE)
  ss <- summary(sf, times = eval.times, extend = TRUE)
  
  # Organize survival estimates from the two treatment profiles.
  out <- data.frame(time = ss$time, surv = ss$surv, std.err = ss$std.err, strata = rep(1:2, each = length(eval.times)))
  
  # Store survival estimates and standard errors for treated and control groups.
  df.CLCOX <- data.frame(time = eval.times, surv1 = out$surv[out$strata == 1], surv1.sd = out$std.err[out$strata == 1], surv0 = out$surv[out$strata == 2], surv0.sd = out$std.err[out$strata == 2])
  
  return(list(df.CLCOX = df.CLCOX))
}


# ------------------------------------------------------------
# Apply the clustered Cox benchmark to all Monte Carlo replicates
#
# Inputs:
#   datlist     : list of simulated datasets
#   eval.times  : vector of time points at which survival is evaluated
#   target_site : target-site index
#
# Output:
#   list of clustered Cox results, one element per Monte Carlo replicate
# ------------------------------------------------------------
run_cluster_cox_sim <- function(datlist, eval.times, target_site = 0) {
  
  # Number of Monte Carlo replicates.
  M <- length(datlist)
  out <- vector("list", M)
  
  # Fit the clustered Cox model separately to each replicate.
  for (j in seq_len(M)) {
    if(j%%100==0) cat("replicate", j, "of", M, "\n")
    out[[j]] <- fit_cluster_cox_one(dat = datlist[[j]], eval.times = eval.times, target_site = target_site)
  }
  
  return(out)
}

# Survival evaluation times used in the simulation summaries.
time.RE <- c(30, 60, 90)

### --- clustered Cox results for source-site sample size n = 300
# Load simulated datasets from the main data-generation script.
load("obsdata_s.Rdata")

# Apply the clustered Cox benchmark to each simulation scenario.
result.CLCOX.homo <- run_cluster_cox_sim(dat.homo, eval.times = time.RE, target_site = 0)
result.CLCOX.diffX <- run_cluster_cox_sim(dat.diffX, eval.times = time.RE, target_site = 0)
result.CLCOX.diffT <- run_cluster_cox_sim(dat.diffT, eval.times = time.RE, target_site = 0)
result.CLCOX.diffC <- run_cluster_cox_sim(dat.diffC, eval.times = time.RE, target_site = 0)
result.CLCOX.diffAll <- run_cluster_cox_sim(dat.diffAll, eval.times = time.RE, target_site = 0)

# Save clustered Cox results for downstream simulation summaries.
save(result.CLCOX.homo, result.CLCOX.diffX, result.CLCOX.diffT, result.CLCOX.diffC, result.CLCOX.diffAll, file = "Res_CLCOX_s.Rdata")


### --- clustered Cox results for source-site sample size n = 600
load("obsdata_l.Rdata")

result.CLCOX.homo <- run_cluster_cox_sim(dat.homo, eval.times = time.RE, target_site = 0)
result.CLCOX.diffX <- run_cluster_cox_sim(dat.diffX, eval.times = time.RE, target_site = 0)
result.CLCOX.diffT <- run_cluster_cox_sim(dat.diffT, eval.times = time.RE, target_site = 0)
result.CLCOX.diffC <- run_cluster_cox_sim(dat.diffC, eval.times = time.RE, target_site = 0)
result.CLCOX.diffAll <- run_cluster_cox_sim(dat.diffAll, eval.times = time.RE, target_site = 0)

save(result.CLCOX.homo, result.CLCOX.diffX, result.CLCOX.diffT, result.CLCOX.diffC, result.CLCOX.diffAll, file = "Res_CLCOX_l.Rdata")

### --- clustered Cox results for source-site sample size n = 1000
load("obsdata_l2.Rdata")

result.CLCOX.homo <- run_cluster_cox_sim(dat.homo, eval.times = time.RE, target_site = 0)
result.CLCOX.diffX <- run_cluster_cox_sim(dat.diffX, eval.times = time.RE, target_site = 0)
result.CLCOX.diffT <- run_cluster_cox_sim(dat.diffT, eval.times = time.RE, target_site = 0)
result.CLCOX.diffC <- run_cluster_cox_sim(dat.diffC, eval.times = time.RE, target_site = 0)
result.CLCOX.diffAll <- run_cluster_cox_sim(dat.diffAll, eval.times = time.RE, target_site = 0)

save(result.CLCOX.homo, result.CLCOX.diffX, result.CLCOX.diffT, result.CLCOX.diffC, result.CLCOX.diffAll, file = "Res_CLCOX_l2.Rdata")

### --- clustered Cox results for the limited-overlap setting
load("obsdata_limO.Rdata")

result.CLCOX.homo <- run_cluster_cox_sim(dat.homo, eval.times = time.RE, target_site = 0)
result.CLCOX.diffX <- run_cluster_cox_sim(dat.diffX, eval.times = time.RE, target_site = 0)
result.CLCOX.diffT <- run_cluster_cox_sim(dat.diffT, eval.times = time.RE, target_site = 0)
result.CLCOX.diffC <- run_cluster_cox_sim(dat.diffC, eval.times = time.RE, target_site = 0)
result.CLCOX.diffAll <- run_cluster_cox_sim(dat.diffAll, eval.times = time.RE, target_site = 0)

save(result.CLCOX.homo, result.CLCOX.diffX, result.CLCOX.diffT, result.CLCOX.diffC, result.CLCOX.diffAll, file = "Res_CLCOX_limO.Rdata")
