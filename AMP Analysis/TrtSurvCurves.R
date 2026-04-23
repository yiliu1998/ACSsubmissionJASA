TrtSurvCurves <- function(data, 
                          covar.name=c("age","score","bweight"), 
                          site.var="site",
                          tgt.name="South Africa",
                          trt.name="A", 
                          time.var="Y", 
                          event="Delta", 
                          fit.times=1:601, 
                          eval.times=1:601,
                          prop.SL.library=c("SL.glm"), 
                          event.SL.library=c("survSL.km", "survSL.coxph", "survSL.rfsrc"),
                          cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.rfsrc"),
                          n.folds=2, 
                          s=1333) {
  
  site <- as.character(data[, site.var])
  tgt.name <- as.character(tgt.name)
  site.names <- unique(site)
  non.tgt.site.names <- site.names[site.names!=tgt.name]
  K <- length(site.names)
  site[site==tgt.name] <- "0"
  for(r in 1:(K-1)) {
    site[site==non.tgt.site.names[r]] <- as.character(r)
  }
  site <- as.numeric(site)
  n.site <- table(site)
  prop.site <- n.site/sum(n.site)
  
  fit.times <- fit.times[fit.times>0]
  N.time <- length(eval.times)
  
  if(!all(eval.times%in%fit.times)) stop("eval.times must be a subset of fit.times")
  eval.ind <- which(fit.times%in%eval.times)
  
  data <- data[order(site),]
  site <- site[order(site)]
  
  set.seed(seed=s)
  seeds <- round(runif(20*K, 0, 20e5))
  
  ################################################################
  ## ~~~~~~~~~~~~~~~ Target-site-only estimator ~~~~~~~~~~~~~~~ ##
  ################################################################
  dat0 <- data[site==0, ]
  A <- dat0[, trt.name]
  Y <- dat0[, time.var]
  Delta <- dat0[, event]
  X <- dat0[, covar.name]
  n <- length(Y)
  
  #### data splitting
  set.seed(seeds[1])
  pred.folds <- createFolds(1:n, k=n.folds, list=T)
  IF.01 <- IF.00 <- S.00 <- S.01 <- NULL
  theta.01 <- theta.00 <- theta.00.sd <- theta.01.sd <- Aug.00.mean <- Aug.01.mean <- rep(0, N.time)
  
  for(i in 1:n.folds) {
    pred.ind <- pred.folds[[i]]
    train.ind <- (1:n)[-pred.ind]
    A.train <- A[train.ind]
    
    #### fit nuisance functions
    ps.fit=SuperLearner(Y=A[train.ind], X=X[train.ind,], 
                        family=binomial(), SL.library=prop.SL.library)
    g.hats=predict(ps.fit, X[pred.ind, ])$pred
    
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
    
    surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
    surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
    
    S.hats.0 <- surv.pred.0$event.SL.predict
    S.hats.1 <- surv.pred.1$event.SL.predict
    G.hats.0 <- surv.pred.0$cens.SL.predict
    G.hats.1 <- surv.pred.1$cens.SL.predict
    
    #### calculate counterfactual survivals
    S1 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=A[pred.ind],
                       fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, g.hats=g.hats)
    S0 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=1-A[pred.ind], 
                       fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, g.hats=1-g.hats)
    
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
  Aug.01.mean <- Aug.01.mean / n.folds
  Aug.00.mean <- Aug.00.mean / n.folds
  theta.00 <- theta.00 / n.folds
  theta.01 <- theta.01 / n.folds
  theta.00.sd <- theta.00.sd / (n.folds*sqrt(n.folds))
  theta.01.sd <- theta.01.sd / (n.folds*sqrt(n.folds))
  
  df.TGT <- data.frame(time=eval.times, 
                       surv1=theta.01, surv1.sd=theta.01.sd, 
                       surv0=theta.00, surv0.sd=theta.00.sd )
  
  ### train models from the target site
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
  ## ~~~~~~~~~ Density-ratio adjusted local estimates ~~~~~~~~~~~ ##
  ##################################################################
  X0 <- as.matrix(dat0[, covar.name])
  Aug.R0.mean <- Aug.R1.mean <- Aug.R0.mean.sour <- Aug.R1.mean.sour <- 
    matrix(0, nrow=N.time, ncol=K-1)
  IF.R0 <- IF.R1 <- df.SOUR <- list()
  for(r in 1:(K-1)) {
    dat.r <- data[site==r, ]
    A <- dat.r[, trt.name]
    Y <- dat.r[, time.var]
    Delta <- dat.r[, event]
    X <- dat.r[, covar.name]
    n <- length(Y)
    
    #### data splitting
    set.seed(seeds[r+1])
    pred.folds <- createFolds(1:n, k=n.folds, list=T)
    IF.R0[[r]] <- IF.R1[[r]] <- matrix(NA, nrow=1, ncol=N.time)
    theta.R1 <- theta.R0 <- theta.R0.sd <- theta.R1.sd <- rep(0, N.time) 
    
    for(i in 1:n.folds) {
      pred.ind <- pred.folds[[i]]
      train.ind <- (1:n)[-pred.ind]
      A.train <- A[train.ind]
      
      ### fit density ratio and propensity scores
      omega.hats <- estimate_omega_np(x=as.matrix(X)[train.ind,], 
                                      x_target=X0, 
                                      x.pred=as.matrix(X)[pred.ind,], 
                                      method="logistic")
      
      ps.fit <- SuperLearner(Y=A[train.ind], X=X[train.ind,], 
                             family=binomial(), SL.library=prop.SL.library)
      g.hats <- predict(ps.fit, X[pred.ind, ])$pred
      
      ### predict conditional event survival from the target site model
      surv.pred.0 <- predict.survSuperLearner(surv.fit.0.tgt, newdata=X[pred.ind,], new.times=fit.times)
      surv.pred.1 <- predict.survSuperLearner(surv.fit.1.tgt, newdata=X[pred.ind,], new.times=fit.times)
      S.hats.0 <- surv.pred.0$event.SL.predict
      S.hats.1 <- surv.pred.1$event.SL.predict
      
      ### for censoring, use the source site model
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
      
      surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
      surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
      G.hats.0 <- surv.pred.0$cens.SL.predict
      G.hats.1 <- surv.pred.1$cens.SL.predict
      
      S1 <- get.survival(Y[pred.ind], Delta[pred.ind], A=A[pred.ind], R=r,
                         fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, 
                         g.hats=g.hats, omega.hats=omega.hats)
      S0 <- get.survival(Y[pred.ind], Delta[pred.ind], A=1-A[pred.ind], R=r,
                         fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, 
                         g.hats=1-g.hats, omega.hats=omega.hats)
      
      Aug.R1.mean[,r] <- Aug.R1.mean[,r] + S1$AUG.means[eval.ind]
      Aug.R0.mean[,r] <- Aug.R0.mean[,r] + S0$AUG.means[eval.ind]
      IF.R1[[r]] <- rbind(IF.R1[[r]], S1$IF.vals[,eval.ind])
      IF.R0[[r]] <- rbind(IF.R0[[r]], S0$IF.vals[,eval.ind])
      
      ### naive source site augmented term estimates
      S.hats.0 <- surv.pred.0$event.SL.predict
      S.hats.1 <- surv.pred.1$event.SL.predict
      
      S1 <- get.survival(Y[pred.ind], Delta[pred.ind], A=A[pred.ind], 
                         fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, g.hats=g.hats)
      S0 <- get.survival(Y[pred.ind], Delta[pred.ind], A=1-A[pred.ind],
                         fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, g.hats=1-g.hats)
    
      Aug.R1.mean.sour[,r] <- Aug.R1.mean.sour[,r] + S1$AUG.means[eval.ind]
      Aug.R0.mean.sour[,r] <- Aug.R0.mean.sour[,r] + S0$AUG.means[eval.ind]
    }
    theta.R0 <- theta.R0 / n.folds
    theta.R1 <- theta.R1 / n.folds
    theta.R0.sd <- theta.R0.sd / (n.folds*sqrt(n.folds))
    theta.R1.sd <- theta.R1.sd / (n.folds*sqrt(n.folds))
    
    IF.R1[[r]] <- IF.R1[[r]][-1,]
    IF.R0[[r]] <- IF.R0[[r]][-1,]
    Aug.R1.mean[,r] <- Aug.R1.mean[,r] / n.folds
    Aug.R0.mean[,r] <- Aug.R0.mean[,r] / n.folds
    Aug.R1.mean.sour[,r] <- Aug.R1.mean.sour[,r] / n.folds
    Aug.R0.mean.sour[,r] <- Aug.R0.mean.sour[,r] / n.folds
  }
  
  ##################################################################
  ## ~~~~~~~~~~~~~~~ Federated weighting estimate ~~~~~~~~~~~~~~~ ##
  ##################################################################
  set.seed(seeds[K+5])
  wt1 <- wt0 <- chi0 <- chi1 <- augdiff0 <- augdiff1 <- matrix(NA, nrow=N.time, ncol=K-1)
  for(i in 1:N.time) {
    IF1.tgt=c(IF.01[,i]-mean(IF.01[,i]), rep(0, length(site[site!=0])))
    IF0.tgt=c(IF.00[,i]-mean(IF.00[,i]), rep(0, length(site[site!=0])))
    IF0.diff <- IF1.diff <- matrix(0, ncol=K-1, nrow=length(IF0.tgt))
    ind0 <- which(site==0)
    for(r in 1:(K-1)) {
      IF0.diff[,r][ind0] <- IF0.tgt[ind0]
      IF1.diff[,r][ind0] <- IF1.tgt[ind0]
      
      ind <- which(site==r)
      IF0.diff[,r][ind] <- -IF.R0[[r]][,i] 
      IF1.diff[,r][ind] <- -IF.R1[[r]][,i] 
      
      chi0[i,r] <- Aug.00.mean[i]-Aug.R0.mean[i,r]
      chi1[i,r] <- Aug.01.mean[i]-Aug.R1.mean[i,r]
      
      augdiff0[i,r] <- Aug.00.mean[i]-Aug.R0.mean.sour[i,r]
      augdiff1[i,r] <- Aug.01.mean[i]-Aug.R1.mean.sour[i,r]
    } 
    
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
  wt0.tgt <- 1-apply(wt0,1,sum)
  wt1.tgt <- 1-apply(wt1,1,sum)
  weights <- cbind(wt0.tgt, wt0, wt1.tgt, wt1)
  
  theta0.fed <- apply(augdiff0*wt0, 1, sum) + theta.00
  theta1.fed <- apply(augdiff1*wt1, 1, sum) + theta.01
  
  all.var0 <- (apply(IF.00,2,var)*(wt0.tgt^2+2*wt0.tgt*(1-wt0.tgt)) + apply(S.00,2,var)*(1-wt0.tgt)^2) / n.site[1] 
  all.var1 <- (apply(IF.01,2,var)*(wt1.tgt^2+2*wt1.tgt*(1-wt1.tgt)) + apply(S.01,2,var)*(1-wt1.tgt)^2) / n.site[1] 
  for(k in 1:(K-1)) {
    all.var0 <- all.var0 + apply(IF.R0[[k]],2,var)*wt0[,k]^2 / n.site[k+1]
    all.var1 <- all.var1 + apply(IF.R1[[k]],2,var)*wt1[,k]^2 / n.site[k+1]
  }
  theta0.fed.sd <- sqrt(all.var0) 
  theta1.fed.sd <- sqrt(all.var1) 
  
  df.FED <- data.frame(time=eval.times, 
                       surv1=theta1.fed, surv1.sd=theta1.fed.sd, 
                       surv0=theta0.fed, surv0.sd=theta0.fed.sd )
  
  ##################################################################
  ## ~~~~~~~~~~~~~~~~~~~~~~~ CCOD estimate ~~~~~~~~~~~~~~~~~~~~~~ ##
  ##################################################################
  A <- data[, trt.name]
  Y <- data[, time.var]
  Delta <- data[, event]
  X <- data[, covar.name]
  R <- as.numeric(site==0)
  n <- length(Y)
  
  #### data splitting
  set.seed(seeds[10*K])
  pred.folds <- createFolds(1:n, k=n.folds, list=T)
  IF.CCOD.0 <- IF.CCOD.1 <- NULL
  theta.ccod.1 <- theta.ccod.0 <- theta.ccod.0.sd <- theta.ccod.1.sd <- rep(0, N.time)
  ind.R1.ccod <- list()
  for(i in 1:n.folds) {
    pred.ind <- pred.folds[[i]]
    train.ind <- (1:n)[-pred.ind]
    A.train <- A[train.ind]
    ind.R1.ccod[[i]] <- c(sum(R[pred.ind]==1),length(pred.ind))
    
    #### fit nuisance functions
    ps.fit=SuperLearner(Y=A[train.ind], X=X[train.ind,], 
                        family=binomial(), SL.library=prop.SL.library)
    g.hats=predict(ps.fit, X[pred.ind, ])$pred
    
    # propensity score of the target site R=0
    eta0.fit=SuperLearner(Y=R[train.ind], X=X[train.ind,], 
                          family=binomial(), SL.library=prop.SL.library)
    eta0.hats=predict(eta0.fit, X[pred.ind, ])$pred
    
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
    
    surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
    surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
    
    S.hats.0 <- surv.pred.0$event.SL.predict
    S.hats.1 <- surv.pred.1$event.SL.predict
    G.hats.0 <- surv.pred.0$cens.SL.predict
    G.hats.1 <- surv.pred.1$cens.SL.predict
    
    #### calculate counterfactual survivals (CCOD)
    S1 <- get.survival.CCOD(Y=Y[pred.ind], Delta=Delta[pred.ind], A=A[pred.ind], R=R[pred.ind],
                            fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, g.hats=g.hats, eta0.hats=eta0.hats)
    S0 <- get.survival.CCOD(Y=Y[pred.ind], Delta=Delta[pred.ind], A=1-A[pred.ind], R=R[pred.ind],
                            fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, g.hats=1-g.hats, eta0.hats=eta0.hats)
    
    IF.CCOD.1 <- rbind(IF.CCOD.1, S1$IF.vals[,eval.ind])
    IF.CCOD.0 <- rbind(IF.CCOD.0, S0$IF.vals[,eval.ind])
    
    theta.ccod.0 <- theta.ccod.0 + S0$surv[eval.ind]
    theta.ccod.1 <- theta.ccod.1 + S1$surv[eval.ind]
    theta.ccod.0.sd <- theta.ccod.0.sd + S0$surv.sd[eval.ind]
    theta.ccod.1.sd <- theta.ccod.1.sd + S1$surv.sd[eval.ind]
  }
  theta.ccod.0 <- theta.ccod.0 / n.folds
  theta.ccod.1 <- theta.ccod.1 / n.folds
  theta.ccod.0.sd <- theta.ccod.0.sd / (n.folds*sqrt(n.folds))
  theta.ccod.1.sd <- theta.ccod.1.sd / (n.folds*sqrt(n.folds))
  
  df.CCOD <- data.frame(time=eval.times, 
                        surv1=theta.ccod.1, surv1.sd=theta.ccod.1.sd, 
                        surv0=theta.ccod.0, surv0.sd=theta.ccod.0.sd )
  
  ind.R1.ccod <- unlist(lapply(seq_along(ind.R1.ccod), function(j) {
    k <- ind.R1.ccod[[j]][1]
    n <- ind.R1.ccod[[j]][2]
    start <- (j - 1) * n + 1
    end   <- (j - 1) * n + k
    start:end
  }))
  return(list(df.TGT=df.TGT, df.FED=df.FED, df.CCOD=df.CCOD,
              IF.00=IF.00, IF.01=IF.01, S.00=S.00, S.01=S.01, 
              IF.R0=IF.R0, IF.R1=IF.R1, IF.CCOD.0=IF.CCOD.0, IF.CCOD.1=IF.CCOD.1, 
              Aug.00.mean=Aug.00.mean, Aug.01.mean=Aug.01.mean, 
              Aug.R0.mean=Aug.R0.mean, Aug.R1.mean=Aug.R1.mean,
              Aug.R0.mean.sour=Aug.R0.mean.sour, Aug.R1.mean.sour=Aug.R1.mean.sour,
              site=site, ind.R1.ccod=ind.R1.ccod, 
              eval.times=eval.times, weights=weights, chi=cbind(chi0, chi1)))
}


POOL_IVW <- function(data, 
                     covar.name=c("age","score","bweight"), 
                     site.var="site",
                     tgt.name="South Africa",
                     trt.name="A", 
                     time.var="Y", 
                     event="Delta", 
                     fit.times=1:601, 
                     eval.times=1:601,
                     prop.SL.library=c("SL.glm"), 
                     event.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                     cens.SL.library=c("survSL.km", "survSL.coxph", "survSL.gam"),
                     n.folds=5,
                     s=2388) {
  
  site <- as.character(data[, site.var])
  tgt.name <- as.character(tgt.name)
  site.names <- unique(site)
  non.tgt.site.names <- site.names[site.names!=tgt.name]
  K <- length(site.names)
  site[site==tgt.name] <- "0"
  for(r in 1:(K-1)) {
    site[site==non.tgt.site.names[r]] <- as.character(r)
  }
  site <- as.numeric(site)
  n.site <- table(site)
  prop.site <- n.site/sum(n.site)
  
  fit.times <- fit.times[fit.times>0]
  N.time <- length(eval.times)
  
  if(!all(eval.times%in%fit.times)) stop("eval.times must be a subset of fit.times")
  eval.ind <- which(fit.times%in%eval.times)
  
  set.seed(seed=s)
  seeds <- round(runif(20*K, 0, 20e5))
  
  #### Target-site estimates ####
  dat0 <- data[site==0, ]
  A <- dat0[, trt.name]
  Y <- dat0[, time.var]
  Delta <- dat0[, event]
  X <- dat0[, covar.name]
  n <- length(Y)
  
  set.seed(seeds[1])
  pred.folds <- createFolds(1:n, k=n.folds, list=T)
  theta.01 <- theta.00 <- theta.00.sd <- theta.01.sd <- rep(0, N.time)
  
  for(i in 1:n.folds) {
    pred.ind <- pred.folds[[i]]
    train.ind <- (1:n)[-pred.ind]
    A.train <- A[train.ind]
    
    ps.fit=SuperLearner(Y=A[train.ind], X=X[train.ind,], 
                        family=binomial(), SL.library=prop.SL.library)
    g.hats=predict(ps.fit, X[pred.ind, ])$pred
    
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
    
    surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
    surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
    
    S.hats.0 <- surv.pred.0$event.SL.predict
    S.hats.1 <- surv.pred.1$event.SL.predict
    G.hats.0 <- surv.pred.0$cens.SL.predict
    G.hats.1 <- surv.pred.1$cens.SL.predict
    
    S1 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=A[pred.ind],
                       fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, g.hats=g.hats)
    S0 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=1-A[pred.ind], 
                       fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, g.hats=1-g.hats)
    
    theta.00 <- theta.00 + S0$surv[eval.ind]
    theta.01 <- theta.01 + S1$surv[eval.ind]
    theta.00.sd <- theta.00.sd + S0$surv.sd[eval.ind]
    theta.01.sd <- theta.01.sd + S1$surv.sd[eval.ind]
  }
  theta.00 <- theta.00 / n.folds
  theta.01 <- theta.01 / n.folds
  theta.00.sd <- theta.00.sd / (n.folds*sqrt(n.folds))
  theta.01.sd <- theta.01.sd / (n.folds*sqrt(n.folds))
  
  df.TGT <- data.frame(time=eval.times, 
                       surv1=theta.01, surv1.sd=theta.01.sd, 
                       surv0=theta.00, surv0.sd=theta.00.sd )
  
  #### Source-site estimates ####
  X0 <- as.matrix(dat0[, covar.name])
  df.SOUR <- list()
  
  for(r in 1:(K-1)) {
    dat.r <- data[site==r, ]
    A <- dat.r[, trt.name]
    Y <- dat.r[, time.var]
    Delta <- dat.r[, event]
    X <- dat.r[, covar.name]
    n <- length(Y)
    
    set.seed(seeds[r+1])
    pred.folds <- createFolds(1:n, k=n.folds, list=T)
    theta.R1 <- theta.R0 <- theta.R0.sd <- theta.R1.sd <- rep(0, N.time) 
    
    for(i in 1:n.folds) {
      pred.ind <- pred.folds[[i]]
      train.ind <- (1:n)[-pred.ind]
      A.train <- A[train.ind]
      
      ps.fit <- SuperLearner(Y=A[train.ind], X=X[train.ind,], 
                             family=binomial(), SL.library=prop.SL.library)
      g.hats <- predict(ps.fit, X[pred.ind, ])$pred
      
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
      
      surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
      surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
      
      G.hats.0 <- surv.pred.0$cens.SL.predict
      G.hats.1 <- surv.pred.1$cens.SL.predict
      S.hats.0 <- surv.pred.0$event.SL.predict
      S.hats.1 <- surv.pred.1$event.SL.predict
      
      S1 <- get.survival(Y[pred.ind], Delta[pred.ind], A=A[pred.ind], 
                         fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, g.hats=g.hats)
      S0 <- get.survival(Y[pred.ind], Delta[pred.ind], A=1-A[pred.ind],
                         fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, g.hats=1-g.hats)
      
      theta.R0 <- theta.R0 + S0$surv[eval.ind]
      theta.R1 <- theta.R1 + S1$surv[eval.ind]
      theta.R0.sd <- theta.R0.sd + S0$surv.sd[eval.ind]
      theta.R1.sd <- theta.R1.sd + S1$surv.sd[eval.ind]
    }
    theta.R0 <- theta.R0 / n.folds
    theta.R1 <- theta.R1 / n.folds
    theta.R0.sd <- theta.R0.sd / (n.folds*sqrt(n.folds))
    theta.R1.sd <- theta.R1.sd / (n.folds*sqrt(n.folds))
    
    df.SOUR[[r]] <- data.frame(time=eval.times, 
                               surv1=theta.R1, surv1.sd=theta.R1.sd, 
                               surv0=theta.R0, surv0.sd=theta.R0.sd )
  }
  
  #### IVW ####
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
  
  #### Simple pooling (POOL) ####
  A <- data[, trt.name]
  Y <- data[, time.var]
  Delta <- data[, event]
  X <- data[, covar.name]
  R <- as.numeric(site==0)
  n <- length(Y)
  
  #### data splitting
  set.seed(seeds[10*K])
  pred.folds <- createFolds(1:n, k=n.folds, list=T)
  theta.pool.1 <- theta.pool.0 <- theta.pool.0.sd <- theta.pool.1.sd <- rep(0, N.time)
  for(i in 1:n.folds) {
    pred.ind <- pred.folds[[i]]
    train.ind <- (1:n)[-pred.ind]
    A.train <- A[train.ind]
    
    #### fit nuisance functions
    ps.fit=SuperLearner(Y=A[train.ind], X=X[train.ind,], 
                        family=binomial(), SL.library=prop.SL.library)
    g.hats=predict(ps.fit, X[pred.ind, ])$pred
    
    # propensity score of the target site R=0
    eta0.fit=SuperLearner(Y=R[train.ind], X=X[train.ind,], 
                          family=binomial(), SL.library=prop.SL.library)
    eta0.hats=predict(eta0.fit, X[pred.ind, ])$pred
    
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
    
    surv.pred.0 <- predict.survSuperLearner(surv.fit.0, newdata=X[pred.ind,], new.times=fit.times)
    surv.pred.1 <- predict.survSuperLearner(surv.fit.1, newdata=X[pred.ind,], new.times=fit.times)
    
    S.hats.0 <- surv.pred.0$event.SL.predict
    S.hats.1 <- surv.pred.1$event.SL.predict
    G.hats.0 <- surv.pred.0$cens.SL.predict
    G.hats.1 <- surv.pred.1$cens.SL.predict
    
    #### calculate counterfactual survivals (simple pooling)
    S1 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=A[pred.ind],
                       fit.times=fit.times, S.hats=S.hats.1, G.hat=G.hats.1, g.hats=g.hats)
    S0 <- get.survival(Y=Y[pred.ind], Delta=Delta[pred.ind], A=1-A[pred.ind], 
                       fit.times=fit.times, S.hats=S.hats.0, G.hat=G.hats.0, g.hats=1-g.hats)
    
    theta.pool.0 <- theta.pool.0 + S0$surv[eval.ind]
    theta.pool.1 <- theta.pool.1 + S1$surv[eval.ind]
    theta.pool.0.sd <- theta.pool.0.sd + S0$surv.sd[eval.ind]
    theta.pool.1.sd <- theta.pool.1.sd + S1$surv.sd[eval.ind]
  }
  ### the simple pooling estimates
  theta.pool.0 <- theta.pool.0 / n.folds
  theta.pool.1 <- theta.pool.1 / n.folds
  theta.pool.0.sd <- theta.pool.0.sd / (n.folds*sqrt(n.folds))
  theta.pool.1.sd <- theta.pool.1.sd / (n.folds*sqrt(n.folds))
  
  df.POOL <- data.frame(time=eval.times, 
                        surv1=theta.pool.1, surv1.sd=theta.pool.1.sd, 
                        surv0=theta.pool.0, surv0.sd=theta.pool.0.sd )
  
  return(list(df.IVW=df.IVW, df.POOL=df.POOL))
}
