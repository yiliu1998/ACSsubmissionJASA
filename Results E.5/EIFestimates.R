# ============================================================
# Utility functions for density-ratio estimation and
# treatment-specific survival estimation
#
# This script defines:
#   (1) estimate_omega_np(): estimates source-to-target density ratios;
#   (2) get.survival(): computes the local/federated survival estimator;
#   (3) get.survival.CCOD(): computes the CCOD survival estimator.
# ============================================================


# ------------------------------------------------------------
# Estimate the source-to-target density ratio
#
# Inputs:
#   x         : covariate matrix from the source site
#   x_target  : covariate matrix from the target site
#   x.pred    : covariate matrix at which the density ratio is evaluated
#   method    : classification method used to distinguish source and
#               target observations; options are "logistic" or "glmnet"
#
# Output:
#   omega     : estimated target-to-source density ratio evaluated at x.pred,
#               truncated to the interval [0.05, 20]
# ------------------------------------------------------------
estimate_omega_np = function(x, x_target, x.pred, method="logistic") {
  
  # Stack source and target covariates and create an indicator for source membership.
  x.all = rbind(x, x_target)
  z.all = c(rep(1, nrow(x)), rep(0, nrow(x_target)))
  
  # Logistic-regression estimator of the source-membership probability.
  if(method=="logistic") {
    colnames(x.all) <- colnames(x.pred) <- paste0("X", 1:ncol(x.all))
    fit <- glm(z.all~.-z.all, data=data.frame(z.all, x.all), family=binomial(link="logit"))
    src.predict = predict(fit, newdata=data.frame(x.pred), type="response")
  }
  
  # Penalized logistic-regression alternative using cross-validated glmnet.
  if(method=="glmnet") {
    fit <- cv.glmnet(x=as.matrix(x.all), y=z.all, nfolds=5, family='binomial')
    src.predict = predict(fit, x.pred, type="response", lambda=fit$lambda.min)
  }
  
  # Convert the estimated source-membership odds into the target-to-source density ratio.
  omega = (1-src.predict)/src.predict * nrow(x)/nrow(x_target)
  omega = pmax(pmin(omega, 20), 0.05)
  return(omega)
}

# ------------------------------------------------------------
# Compute the treatment-specific survival estimator
#
# Inputs:
#   Y          : observed follow-up time
#   Delta      : event indicator
#   A          : treatment indicator
#   R          : site-membership indicator used in the influence-function
#                representation; default R=0
#   fit.times  : time grid at which the survival curve is estimated
#   S.hats     : estimated conditional event-free survival probabilities
#   G.hats     : estimated conditional censoring survival probabilities
#   g.hats     : estimated treatment propensity scores
#   omega.hats : estimated density-ratio weights; default is 1
#
# Output:
#   IF.vals    : estimated influence-function contributions over fit.times
#   AUG.means  : average augmentation term over fit.times
#   surv       : estimated marginal survival probabilities
#   surv.sd    : estimated standard errors
# ------------------------------------------------------------
get.survival <- function(Y, Delta, A, R=0,
                         fit.times, 
                         S.hats, G.hats, g.hats, omega.hats=1) {
  
  # Remove time zero, sort evaluation times, and reorder nuisance predictions accordingly.
  fit.times <- fit.times[fit.times > 0]
  n <- length(Y)
  ord <- order(fit.times)
  fit.times <- fit.times[ord]
  S.hats <- S.hats[, ord]
  G.hats <- G.hats[, ord]
  
  # Compute the discrete-time integral term appearing in the augmentation component.
  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S.hats[i,])* 1/ G.hats[i,-ncol(G.hats)]
    if(any(fit.times[-1] > Y[i])) vals[fit.times[-1] > Y[i]] <- 0
    c(0, cumsum(vals))
  }))
  
  # Evaluate the estimated event and censoring survival functions at each subject's observed time.
  S.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
  G.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))
  
  # Initialize storage for influence-function values, survival estimates, and augmentation means.
  IF.vals <- matrix(NA, nrow=n, ncol=length(fit.times))
  surv <- AUG.means <- rep(NA, length(fit.times))
  
  # Evaluate the estimator separately at each requested time point.
  for(t0 in fit.times) {
    k <- min(which(fit.times>=t0))
    S.hats.t0 <- S.hats[,k]
    
    # Event contribution and cumulative integral component.
    inner.func.1 <- ifelse(Y<=t0 & Delta==1, 1/(S.hats.Y*G.hats.Y), 0 )
    inner.func.2 <- int.vals[,k]
    k1 <- which(fit.times==t0)
    
    # Density-ratio weighted augmentation term for treated subjects.
    augment <- omega.hats*S.hats.t0*as.numeric(A==1)*(inner.func.1 - inner.func.2)/g.hats
    
    # Plug-in plus augmentation representation of the survival estimator.
    if.func <- S.hats.t0 - augment
    surv[k1] <- mean(if.func)
    
    # Store target- and source-specific influence-function contributions.
    IF.vals[,k1] <- if.func*I(R==0) - augment*I(R!=0)
    
    AUG.means[k1] <- mean(augment)
  }
  
  # Restrict estimated survival probabilities to the admissible range.
  surv = pmin(1, pmax(0, surv))
  
  # Estimate pointwise standard errors from the empirical variance of the IF values.
  surv.sd <- sqrt(apply(IF.vals, 2, var, na.rm=T)/n)
  
  return(list(IF.vals=IF.vals, AUG.means=AUG.means, surv=surv, surv.sd=surv.sd))
}


# ------------------------------------------------------------
# Compute the CCOD treatment-specific survival estimator
#
# Inputs:
#   Y         : observed follow-up time
#   Delta     : event indicator
#   A         : treatment indicator
#   R         : indicator for membership in the target population
#   fit.times : time grid at which the survival curve is estimated
#   S.hats    : estimated conditional event-free survival probabilities
#   G.hats    : estimated conditional censoring survival probabilities
#   g.hats    : estimated treatment propensity scores
#   eta0.hats : estimated CCOD weighting/selection component
#
# Output:
#   IF.vals   : estimated influence-function contributions over fit.times
#   surv      : estimated marginal survival probabilities
#   surv.sd   : estimated standard errors
# ------------------------------------------------------------
get.survival.CCOD <- function(Y, Delta, A, R,
                              fit.times, 
                              S.hats, G.hats, g.hats, eta0.hats) {
  
  # Remove time zero, sort evaluation times, and reorder nuisance predictions accordingly.
  fit.times <- fit.times[fit.times > 0]
  n <- length(Y)
  ord <- order(fit.times)
  fit.times <- fit.times[ord]
  S.hats <- S.hats[, ord]
  G.hats <- G.hats[, ord]
  
  # Compute the discrete-time integral term appearing in the augmentation component.
  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S.hats[i,])* 1/ G.hats[i,-ncol(G.hats)]
    if(any(fit.times[-1] > Y[i])) vals[fit.times[-1] > Y[i]] <- 0
    c(0, cumsum(vals))
  }))
  
  # Evaluate the estimated event and censoring survival functions at each subject's observed time.
  S.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
  G.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))
  
  # Initialize storage for influence-function values, survival estimates, and standard errors.
  IF.vals <- matrix(NA, nrow=n, ncol=length(fit.times))
  surv <- surv.sd <- rep(NA, length(fit.times))
  
  # Evaluate the CCOD estimator separately at each requested time point.
  for(t0 in fit.times) {
    k <- min(which(fit.times>=t0))
    S.hats.t0 <- S.hats[,k]
    
    # Event contribution and cumulative integral component.
    inner.func.1 <- ifelse(Y<=t0 & Delta==1, 1/(S.hats.Y*G.hats.Y), 0 )
    inner.func.2 <- int.vals[,k]
    k1 <- which(fit.times==t0)
    
    # Augmentation term for treated subjects.
    augment <- S.hats.t0*as.numeric(A==1)*(inner.func.1 - inner.func.2)/g.hats
    
    # CCOD influence-function representation.
    if.func <- (S.hats.t0*I(R==1) - eta0.hats*augment) / mean(R==1)
    surv[k1] <- mean(if.func)
    IF.vals[,k1] <- if.func 
    
    # Estimate the variance by combining contributions from target and non-target observations.
    surv.sd[k1] <- sqrt(var(if.func[R==1])*mean(R==1) + var(if.func[R!=1])*mean(R!=1)) / sqrt(n)
  }
  
  # Restrict estimated survival probabilities to the admissible range.
  surv = pmin(1, pmax(0, surv))
  
  return(list(IF.vals=IF.vals, surv=surv, surv.sd=surv.sd))
}
