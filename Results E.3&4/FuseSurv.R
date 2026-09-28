# ============================================================
# Main estimation function for the simulation study
#
# FuseSurv() computes five estimators of treatment-specific survival
# probabilities for a target site:
#   TGT  : target-site estimator
#   IVW  : inverse-variance weighted estimator
#   POOL : simple pooled estimator
#   CCOD : estimator under the common conditional outcome distribution
#   FED  : proposed data-adaptive federated estimator
#
# Inputs:
#   data              : combined individual-level dataset containing all sites
#   covar.name        : names of baseline covariates
#   site.var          : site indicator; site 0 is treated as the target site
#   trt.name          : binary treatment variable
#   time.var          : observed follow-up time
#   event             : event indicator
#   fit.times         : time grid used to fit/evaluate nuisance survival functions
#   eval.times        : time points at which final survival estimates are reported
#   prop.SL.library   : Super Learner library for propensity-score models
#   event.SL.library  : survival Super Learner library for event-time models
#   cens.SL.library   : survival Super Learner library for censoring models
#   n.folds           : number of folds used for cross-fitting
#   s                 : random seed
#
# Output:
#   a list containing:
#     df.TGT   : target-site survival estimates and standard errors
#     df.IVW   : inverse-variance weighted estimates and standard errors
#     df.POOL  : simple pooled estimates and standard errors
#     df.CCOD  : CCOD estimates and standard errors
#     df.FED   : federated estimates and standard errors
#     weights  : estimated federated weights for control and treatment
#     chi      : estimated target-source discrepancy measures used in FED
# ============================================================

FuseSurv <- function(data, 
                     covar.name=c("X1","X2","X3"), 
                     site.var="site", 
                     trt.name="A", 
                     time.var="Y", 
                     event="Delta", 
                     fit.times=1:90, 
                     eval.times=c(30,60,90),
                     prop.SL.library=c("SL.mean", "SL.glm"), 
                     event.SL.library=c("survSL.km", "survSL.coxph", "survSL.weibreg"), 
                     cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.weibreg"),
                     n.folds=5,
                     s=1) {
  
  # Basic site/sample-size information and evaluation-time indices.
  site <- data[, site.var]
  K <- length(unique(site))
  n.site <- table(site)
  prop.site <- n.site/sum(n.site)
  fit.times <- fit.times[fit.times>0]
  N.time <- length(eval.times)
  eval.ind <- which(fit.times%in%eval.times)
  
  # Generate reproducible seeds for the cross-fitting and weighting steps.
  set.seed(seed=s)
  seeds <- round(runif(20*K, 0, 20e5))
  
  
  ##################################################################
  ## ~~~~~~~~~~~~~~~ Step 1: target site estimate ~~~~~~~~~~~~~~~ ##
  ##################################################################
  
  # Restrict to the target site (site 0) and extract analysis variables.
  dat0 <- data[site==0, ]
  A <- dat0[, trt.name]
  Y <- dat0[, time.var]
  Delta <- dat0[, event]
  X <- dat0[, covar.name]
  n <- length(Y)
  
  #### data splitting
  
  # Create cross-fitting folds within the target site.
  set.seed(seeds[1])
  pred.folds <- createFolds(1:n, k=n.folds, list=T)
  IF.01 <- IF.00 <- S.00 <- S.01 <- NULL
  theta.01 <- theta.00 <- theta.00.sd <- theta.01.sd <- Aug.00.mean <- Aug.01.mean <- rep(0, N.time)
  
  for(i in 1:n.folds) {
    pred.ind <- pred.folds[[i]]
    train.ind <- (1:n)[-pred.ind]
    A.train <- A[train.ind]
    
    #### fit nuisance functions
    
    # Estimate the treatment propensity score using Super Learner on the training fold.
    ps.fit=SuperLearner(Y=A[train.ind], X=X[train.ind,], 
                        family=binomial(), SL.library=prop.SL.library)
    g.hats=predict(ps.fit, X[pred.ind, ])$pred
    
    # Fit conditional event and censoring survival models among untreated subjects.
    surv.fit.0=survSuperLearner(time=Y[train.ind][A.train==0], 
                                event=Delta[train.ind][A.train==0], 
                                X=X[train.ind,][A.train==0,], 
                                new.times=fit.times,
                                event.SL.library=event.SL.library, 
                                cens.SL.library=cens.SL.library)
    
    # Fit conditional event and censoring survival models among treated subjects.
    surv.fit.1=survSuperLearner(time=Y[train.ind][A.train==1], 
                                event=Delta[train.ind][A.train==1], 
                                X=X[train.ind,][A.train==1,], 
                                new.times=fit.times,
                                event.SL.library=event.SL.library, 
                                cens.SL.library=cens.SL.library)
    
    # Obtain cross-fitted nuisance predictions for the held-out fold.
    surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
    surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
    
    S.hats.0 <- surv.pred.0$event.SL.predict
    S.hats.1 <- surv.pred.1$event.SL.predict
    G.hats.0 <- surv.pred.0$cens.SL.predict
    G.hats.1 <- surv.pred.1$cens.SL.predict
    
    #### calculate counterfactual survivals
    
    # Compute cross-fitted target-site survival estimates under treatment and control.
    S1 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=A[pred.ind],
                       fit.times=fit.times, S.hats=S.hats.1, G.hats=G.hats.1, g.hats=g.hats)
    S0 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=1-A[pred.ind], 
                       fit.times=fit.times, S.hats=S.hats.0, G.hats=G.hats.0, g.hats=1-g.hats)
    
    # Accumulate fold-specific estimates, standard errors, influence-function values,
    # nuisance predictions, and augmentation terms.
    theta.00 <- theta.00 + S0$surv[eval.ind]
    theta.01 <- theta.01 + S1$surv[eval.ind]
    theta.00.sd <- theta.00.sd + S0$surv.sd[eval.ind]
    theta.01.sd <- theta.01.sd + S1$surv.sd[eval.ind]
    IF.01 <- rbind(IF.01, S1$IF.vals[,eval.ind])
    IF.00 <- rbind(IF.00, S0$IF.vals[,eval.ind])
    S.00 <- rbind(S.00, S.hats.0[,eval.ind])
    S.01 <- rbind(S.01, S.hats.1[,eval.ind])
    Aug.01.mean <- Aug.01.mean + S1$AUG.means[eval.ind]
    Aug.00.mean <- Aug.00.mean + S0$AUG.means[eval.ind]
  }
  
  # Average the estimates and augmentation terms across cross-fitting folds.
  Aug.01.mean <- Aug.01.mean / n.folds
  Aug.00.mean <- Aug.00.mean / n.folds
  theta.00 <- theta.00 / n.folds
  theta.01 <- theta.01 / n.folds
  theta.00.sd <- theta.00.sd / (n.folds*sqrt(n.folds))
  theta.01.sd <- theta.01.sd / (n.folds*sqrt(n.folds))
  
  # Store target-site estimates at the requested evaluation times.
  df.TGT <- data.frame(time=eval.times, 
                       surv1=theta.01, surv1.sd=theta.01.sd, 
                       surv0=theta.00, surv0.sd=theta.00.sd )
  
  ### train models from the target site
  
  # Refit the target-site conditional survival models using all target-site observations.
  # These models are subsequently evaluated on source-site covariates in Step 2.
  surv.fit.0.tgt=survSuperLearner(time=Y[A==0], 
                                  event=Delta[A==0], 
                                  X=X[A==0,], 
                                  new.times=fit.times,
                                  event.SL.library=event.SL.library, 
                                  cens.SL.library=cens.SL.library)
  
  surv.fit.1.tgt=survSuperLearner(time=Y[A==1], 
                                  event=Delta[A==1], 
                                  X=X[A==1,], 
                                  new.times=fit.times,
                                  event.SL.library=event.SL.library, 
                                  cens.SL.library=cens.SL.library)
  
  
  ##################################################################
  ## ~~~~~~~~~~~~~~ Step 2: source site estimates ~~~~~~~~~~~~~~~ ##
  ##################################################################
  
  # Target-site covariates are used to estimate source-to-target density ratios.
  X0 <- as.matrix(dat0[, covar.name])
  
  # Initialize storage for source-site augmentation terms, influence functions,
  # and naive source-specific survival estimates.
  Aug.R0.mean <- Aug.R1.mean <- Aug.R0.mean.sour <- Aug.R1.mean.sour <- matrix(0, nrow=N.time, ncol=K-1)
  IF.R0 <- IF.R1 <- df.SOUR <- list()
  
  # Process each source site separately.
  for(r in 1:(K-1)) {
    dat.r <- data[site==r, ]
    A <- dat.r[, trt.name]
    Y <- dat.r[, time.var]
    Delta <- dat.r[, event]
    X <- dat.r[, covar.name]
    n <- length(Y)
    
    #### data splitting
    
    # Create cross-fitting folds within the current source site.
    set.seed(seeds[r+1])
    pred.folds <- createFolds(1:n, k=n.folds, list=T)
    IF.R0[[r]] <- IF.R1[[r]] <- matrix(NA, nrow=1, ncol=N.time)
    theta.R1 <- theta.R0 <- theta.R0.sd <- theta.R1.sd <- rep(0, N.time) 
    
    for(i in 1:n.folds) {
      pred.ind <- pred.folds[[i]]
      train.ind <- (1:n)[-pred.ind]
      A.train <- A[train.ind]
      
      ### fit density ratio and propensity scores
      
      # Estimate the target-to-source covariate density ratio using logistic regression.
      omega.hats <- estimate_omega_np(x=as.matrix(X)[train.ind,], 
                                      x_target=X0, 
                                      x.pred=as.matrix(X)[pred.ind,], 
                                      method="logistic")
      
      # Estimate the source-site treatment propensity score using Super Learner.
      ps.fit <- SuperLearner(Y=A[train.ind], X=X[train.ind,], 
                             family=binomial(), SL.library=prop.SL.library)
      g.hats <- predict(ps.fit, X[pred.ind, ])$pred
      
      ### predict conditional event survival from the target site model
      
      # Evaluate the target-trained event-survival models on source-site covariates.
      surv.pred.0 <- predict.survSuperLearner(surv.fit.0.tgt, newdata=X[pred.ind,], new.times=fit.times)
      surv.pred.1 <- predict.survSuperLearner(surv.fit.1.tgt, newdata=X[pred.ind,], new.times=fit.times)
      S.hats.0 <- surv.pred.0$event.SL.predict
      S.hats.1 <- surv.pred.1$event.SL.predict
      
      ### for censoring, use the source site model
      
      # Fit source-specific censoring nuisance models separately by treatment group.
      surv.fit.0=survSuperLearner(time=Y[train.ind][A.train==0], 
                                  event=Delta[train.ind][A.train==0], 
                                  X=X[train.ind,][A.train==0,], 
                                  new.times=fit.times,
                                  event.SL.library=event.SL.library, 
                                  cens.SL.library=cens.SL.library)
      
      surv.fit.1=survSuperLearner(time=Y[train.ind][A.train==1], 
                                  event=Delta[train.ind][A.train==1], 
                                  X=X[train.ind,][A.train==1,], 
                                  new.times=fit.times,
                                  event.SL.library=event.SL.library, 
                                  cens.SL.library=cens.SL.library)
      
      # Obtain source-specific censoring-survival predictions.
      surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
      surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
      G.hats.0 <- surv.pred.0$cens.SL.predict
      G.hats.1 <- surv.pred.1$cens.SL.predict
      
      # Construct source contributions to the target estimand using density-ratio
      # weighting, target-trained event models, and source-specific censoring models.
      S1 <- get.survival(Y[pred.ind], Delta[pred.ind], A=A[pred.ind], R=r,
                         fit.times=fit.times, S.hats=S.hats.1, G.hats=G.hats.1, 
                         g.hats=g.hats, omega.hats=omega.hats)
      S0 <- get.survival(Y[pred.ind], Delta[pred.ind], A=1-A[pred.ind], R=r,
                         fit.times=fit.times, S.hats=S.hats.0, G.hats=G.hats.0, 
                         g.hats=1-g.hats, omega.hats=omega.hats)
      
      # Store source-site augmentation and influence-function contributions used by FED.
      Aug.R1.mean[,r] <- Aug.R1.mean[,r] + S1$AUG.means[eval.ind]
      Aug.R0.mean[,r] <- Aug.R0.mean[,r] + S0$AUG.means[eval.ind]
      IF.R1[[r]] <- rbind(IF.R1[[r]], S1$IF.vals[,eval.ind])
      IF.R0[[r]] <- rbind(IF.R0[[r]], S0$IF.vals[,eval.ind])
      
      ### naive source site estimates
      ### use the source site model prediction for survival
      
      # Replace target-trained survival predictions by source-trained predictions
      # to obtain the naive source-specific estimates used by IVW.
      S.hats.0 <- surv.pred.0$event.SL.predict
      S.hats.1 <- surv.pred.1$event.SL.predict
      
      S1 <- get.survival(Y[pred.ind], Delta[pred.ind], A=A[pred.ind], 
                         fit.times=fit.times, S.hats=S.hats.1, G.hats=G.hats.1, g.hats=g.hats)
      S0 <- get.survival(Y[pred.ind], Delta[pred.ind], A=1-A[pred.ind],
                         fit.times=fit.times, S.hats=S.hats.0, G.hats=G.hats.0, g.hats=1-g.hats)
      
      # Accumulate naive source-specific estimates and augmentation terms across folds.
      theta.R0 <- theta.R0 + S0$surv[eval.ind]
      theta.R1 <- theta.R1 + S1$surv[eval.ind]
      theta.R0.sd <- theta.R0.sd + S0$surv.sd[eval.ind]
      theta.R1.sd <- theta.R1.sd + S1$surv.sd[eval.ind]
      Aug.R1.mean.sour[,r] <- Aug.R1.mean.sour[,r] + S1$AUG.means[eval.ind]
      Aug.R0.mean.sour[,r] <- Aug.R0.mean.sour[,r] + S0$AUG.means[eval.ind]
    }
    
    # Average source-specific estimates across cross-fitting folds.
    theta.R0 <- theta.R0 / n.folds
    theta.R1 <- theta.R1 / n.folds
    theta.R0.sd <- theta.R0.sd / (n.folds*sqrt(n.folds))
    theta.R1.sd <- theta.R1.sd / (n.folds*sqrt(n.folds))
    
    df.SOUR[[r]] <- data.frame(time=eval.times, 
                               surv1=theta.R1, surv1.sd=theta.R1.sd, 
                               surv0=theta.R0, surv0.sd=theta.R0.sd )
    
    # Remove initialization rows and average augmentation terms across folds.
    IF.R1[[r]] <- IF.R1[[r]][-1,]
    IF.R0[[r]] <- IF.R0[[r]][-1,]
    Aug.R1.mean[,r] <- Aug.R1.mean[,r] / n.folds
    Aug.R0.mean[,r] <- Aug.R0.mean[,r] / n.folds
    Aug.R1.mean.sour[,r] <- Aug.R1.mean.sour[,r] / n.folds
    Aug.R0.mean.sour[,r] <- Aug.R0.mean.sour[,r] / n.folds
  }
  
  
  ##################################################################
  ## ~~~~~~~~~~~ Step 3: inverse variance weighting ~~~~~~~~~~~~~ ##
  ##################################################################
  
  # Combine the target and naive source estimates using inverse-variance weights.
  df.IVW <- data.frame(time=df.TGT$time, surv1=NA, surv1.sd=NA, surv0=NA, surv0.sd=NA)
  
  for (i in seq_along(df.TGT$time)) {
    tgt.surv1 <- df.TGT$surv1[i]
    tgt.var1 <- df.TGT$surv1.sd[i]^2
    tgt.surv0 <- df.TGT$surv0[i]
    tgt.var0 <- df.TGT$surv0.sd[i]^2
    
    w.surv1 <- 0
    w.var1 <- 0
    w.surv0 <- 0
    w.var0 <- 0
    
    # Add inverse-variance weighted contributions from each source site.
    for (r in 1:(K-1)) {
      src.surv1 <- df.SOUR[[r]]$surv1[i]
      src.var1 <- df.SOUR[[r]]$surv1.sd[i]^2
      src.surv0 <- df.SOUR[[r]]$surv0[i]
      src.var0 <- df.SOUR[[r]]$surv0.sd[i]^2
      
      # inverse variance weights
      w.surv1 <- w.surv1 + src.surv1 / src.var1
      w.var1 <- w.var1 + 1 / src.var1
      w.surv0 <- w.surv0 + src.surv0 / src.var0
      w.var0 <- w.var0 + 1 / src.var0
    }
    
    # add the target site into the weighted sums
    w.surv1 <- w.surv1 + tgt.surv1 / tgt.var1
    w.var1 <- w.var1 + 1 / tgt.var1
    w.surv0 <- w.surv0 + tgt.surv0 / tgt.var0
    w.var0 <- w.var0 + 1 / tgt.var0
    
    # compute the IVW estimates and variances
    df.IVW$surv1[i] <- w.surv1 / w.var1
    df.IVW$surv1.sd[i] <- sqrt(1 / w.var1)
    df.IVW$surv0[i] <- w.surv0 / w.var0
    df.IVW$surv0.sd[i] <- sqrt(1 / w.var0)
  }
  
  
  ##################################################################
  ## ~~~~~~~ Step 4: data-adaptive weighting (federated) ~~~~~~~~ ##
  ##################################################################
  
  # Estimate FED source weights separately by treatment level and evaluation time.
  set.seed(seeds[K+5])
  wt1 <- wt0 <- chi0 <- chi1 <- augdiff0 <- augdiff1 <- matrix(NA, nrow=N.time, ncol=K-1)
  
  for(i in 1:N.time) {
    
    # Embed target influence-function contributions in the combined sample.
    IF1.tgt=c(IF.01[,i], rep(0, length(site[site!=0])))
    IF0.tgt=c(IF.00[,i], rep(0, length(site[site!=0])))
    IF0.diff <- IF1.diff <- matrix(0, ncol=K-1, nrow=length(IF0.tgt))
    ind0 <- which(site==0)
    
    for(r in 1:(K-1)) {
      
      # Construct target-minus-source influence-function contrasts for each source.
      IF0.diff[,r][ind0] <- IF0.tgt[ind0]
      IF1.diff[,r][ind0] <- IF1.tgt[ind0]
      
      ind <- which(site==r)
      IF0.diff[,r][ind] <- -IF.R0[[r]][,i] 
      IF1.diff[,r][ind] <- -IF.R1[[r]][,i] 
      
      # chi measures target-source discrepancy and enters the adaptive penalty.
      chi0[i,r] <- Aug.00.mean[i]-Aug.R0.mean[i,r]
      chi1[i,r] <- Aug.01.mean[i]-Aug.R1.mean[i,r]
      
      # Difference between the target augmentation and naive source augmentation,
      # used to construct the final federated estimator.
      augdiff0[i,r] <- Aug.00.mean[i]-Aug.R0.mean.sour[i,r]
      augdiff1[i,r] <- Aug.01.mean[i]-Aug.R1.mean.sour[i,r]
    } 
    
    # Cross-validate the regularization level and estimate nonnegative source
    # weights for the control survival curve. If fitting fails, assign zero weights.
    cvfit0=try(cv.glmnet(x=IF0.diff, y=IF0.tgt))
    if(class(cvfit0)[1]!="try-error") {
      fit0=try(glmnet(x=IF0.diff, y=IF0.tgt, 
                      penalty.factor=chi0[i,]^2,
                      intercept=FALSE,
                      alpha=1,
                      lambda=cvfit0$lambda.1se,
                      lower.limits=0,
                      upper.limits=1))
      if(class(fit0)[1]!="try-error") { 
        wt0[i,]=coef(fit0, s=cvfit0$lambda.1se)[-1] } else { wt0[i,]=rep(0, K-1) }
    } else { wt0[i,]=rep(0, K-1) }
    
    # Repeat the data-adaptive weight estimation for the treated survival curve.
    cvfit1=try(cv.glmnet(x=IF1.diff, y=IF1.tgt))
    if(class(cvfit1)[1]!="try-error") {
      fit1=try(glmnet(x=IF1.diff, y=IF1.tgt, 
                      penalty.factor=chi1[i,]^2,
                      intercept=FALSE,
                      alpha=1,
                      lambda=cvfit1$lambda.1se,
                      lower.limits=0,
                      upper.limits=1))
      if(class(fit1)[1]!="try-error") { 
        wt1[i,]=coef(fit1, s=cvfit1$lambda.1se)[-1] } else { wt1[i,]=rep(0, K-1) }
    } else { wt1[i,]=rep(0, K-1) }
  } 
  
  # The target weight is one minus the sum of source weights.
  wt0.tgt <- 1-apply(wt0,1,sum)
  wt1.tgt <- 1-apply(wt1,1,sum)
  weights <- cbind(wt0.tgt,wt0, wt1.tgt,wt1)
  
  # Construct the final federated survival estimates.
  theta0.fed <- apply(augdiff0*wt0, 1, sum) + theta.00
  theta1.fed <- apply(augdiff1*wt1, 1, sum) + theta.01
  
  # Estimate the variance of the federated estimator using target- and source-specific influence-function contributions.
  all.var0 <- (colMeans(IF.00^2)*((wt0.tgt)^2+2*wt0.tgt*(1-wt0.tgt)) + apply(S.00, 2, var)*(1-wt0.tgt)^2) / n.site[1] 
  all.var1 <- (colMeans(IF.01^2)*((wt1.tgt)^2+2*wt1.tgt*(1-wt1.tgt)) + apply(S.01, 2, var)*(1-wt1.tgt)^2) / n.site[1] 
  
  for(k in 1:(K-1)) {
    all.var0 <- all.var0 + apply(IF.R0[[k]], 2, var)*wt0[,k]^2 / n.site[k+1]
    all.var1 <- all.var1 + apply(IF.R1[[k]], 2, var)*wt1[,k]^2 / n.site[k+1]
  }
  
  theta0.fed.sd <- sqrt(all.var0) 
  theta1.fed.sd <- sqrt(all.var1) 
  
  # Store FED estimates and standard errors.
  df.FED <- data.frame(time=eval.times, 
                       surv1=theta1.fed, surv1.sd=theta1.fed.sd, 
                       surv0=theta0.fed, surv0.sd=theta0.fed.sd )
  
  
  #########################################################
  ## ~~~~ Step 5: simple pooling and CCOD estimates ~~~~ ##
  #########################################################
  
  # Use all observations jointly for the POOL and CCOD estimators.
  A <- data[, trt.name]
  Y <- data[, time.var]
  Delta <- data[, event]
  X <- data[, covar.name]
  R <- as.numeric(data[, site.var]==0)
  n <- length(Y)
  
  #### data splitting
  
  # Create cross-fitting folds over the pooled sample.
  set.seed(seeds[10*K])
  pred.folds <- createFolds(1:n, k=n.folds, list=T)
  theta.pool.1 <- theta.pool.0 <- theta.pool.0.sd <- theta.pool.1.sd <- rep(0, N.time)
  theta.ccod.1 <- theta.ccod.0 <- theta.ccod.0.sd <- theta.ccod.1.sd <- rep(0, N.time)
  
  for(i in 1:n.folds) {
    pred.ind <- pred.folds[[i]]
    train.ind <- (1:n)[-pred.ind]
    A.train <- A[train.ind]
    
    #### fit nuisance functions
    
    # Pooled treatment propensity score.
    ps.fit=SuperLearner(Y=A[train.ind], X=X[train.ind,], 
                        family=binomial(), SL.library=prop.SL.library)
    g.hats=predict(ps.fit, X[pred.ind, ])$pred
    
    # Target-site propensity score. 
    eta0.fit=SuperLearner(Y=R[train.ind], X=X[train.ind,], 
                          family=binomial(), SL.library=prop.SL.library)
    eta0.hats=predict(eta0.fit, X[pred.ind, ])$pred
    
    # Fit pooled conditional event and censoring survival nuisance models.
    surv.fit.0=survSuperLearner(time=Y[train.ind][A.train==0], 
                                event=Delta[train.ind][A.train==0], 
                                X=X[train.ind,][A.train==0,], 
                                new.times=fit.times,
                                event.SL.library=event.SL.library, 
                                cens.SL.library=cens.SL.library)
    
    surv.fit.1=survSuperLearner(time=Y[train.ind][A.train==1], 
                                event=Delta[train.ind][A.train==1], 
                                X=X[train.ind,][A.train==1,], 
                                new.times=fit.times,
                                event.SL.library=event.SL.library, 
                                cens.SL.library=cens.SL.library)
    
    # Obtain cross-fitted event and censoring survival predictions.
    surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
    surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
    
    S.hats.0 <- surv.pred.0$event.SL.predict
    S.hats.1 <- surv.pred.1$event.SL.predict
    G.hats.0 <- surv.pred.0$cens.SL.predict
    G.hats.1 <- surv.pred.1$cens.SL.predict
    
    #### calculate simple pooling's conditional survival nuisances
    S1 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=A[pred.ind],
                       fit.times=fit.times, S.hats=S.hats.1, G.hats=G.hats.1, g.hats=g.hats)
    S0 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=1-A[pred.ind], 
                       fit.times=fit.times, S.hats=S.hats.0, G.hats=G.hats.0, g.hats=1-g.hats)
    
    theta.pool.0 <- theta.pool.0 + S0$surv[eval.ind]
    theta.pool.1 <- theta.pool.1 + S1$surv[eval.ind]
    theta.pool.0.sd <- theta.pool.0.sd + S0$surv.sd[eval.ind]
    theta.pool.1.sd <- theta.pool.1.sd + S1$surv.sd[eval.ind]
    
    #### calculate CCOD's conditional survival nuisances
    S1 <- get.survival.CCOD(Y=Y[pred.ind], Delta=Delta[pred.ind], A=A[pred.ind], R=R[pred.ind],
                            fit.times=fit.times, S.hats=S.hats.1, G.hats=G.hats.1, g.hats=g.hats, eta0.hats=eta0.hats)
    S0 <- get.survival.CCOD(Y=Y[pred.ind], Delta=Delta[pred.ind], A=1-A[pred.ind], R=R[pred.ind],
                            fit.times=fit.times, S.hats=S.hats.0, G.hats=G.hats.0, g.hats=1-g.hats, eta0.hats=eta0.hats)
    
    theta.ccod.0 <- theta.ccod.0 + S0$surv[eval.ind]
    theta.ccod.1 <- theta.ccod.1 + S1$surv[eval.ind]
    theta.ccod.0.sd <- theta.ccod.0.sd + S0$surv.sd[eval.ind]
    theta.ccod.1.sd <- theta.ccod.1.sd + S1$surv.sd[eval.ind]
  }
  
  ### the simple pooling estimates
  theta.pool.0 <- theta.pool.0 / n.folds
  theta.pool.1 <- theta.pool.1 / n.folds
  theta.pool.0.sd <- theta.pool.0.sd / (n.folds*sqrt(n.folds))
  theta.pool.1.sd <- theta.pool.1.sd / (n.folds*sqrt(n.folds))
  
  df.POOL <- data.frame(time=eval.times, 
                        surv1=theta.pool.1, surv1.sd=theta.pool.1.sd, 
                        surv0=theta.pool.0, surv0.sd=theta.pool.0.sd )
  
  ### the CCOD estimates
  theta.ccod.0 <- theta.ccod.0 / n.folds
  theta.ccod.1 <- theta.ccod.1 / n.folds
  theta.ccod.0.sd <- theta.ccod.0.sd / (n.folds*sqrt(n.folds))
  theta.ccod.1.sd <- theta.ccod.1.sd / (n.folds*sqrt(n.folds))
  
  df.CCOD <- data.frame(time=eval.times, 
                        surv1=theta.ccod.1, surv1.sd=theta.ccod.1.sd, 
                        surv0=theta.ccod.0, surv0.sd=theta.ccod.0.sd )
  
  # Return all estimators together with the FED weights and discrepancy measures.
  return(list(df.TGT=df.TGT, df.IVW=df.IVW, df.POOL=df.POOL, 
              df.CCOD=df.CCOD, df.FED=df.FED, 
              weights=weights, chi=cbind(chi0, chi1)))
}
