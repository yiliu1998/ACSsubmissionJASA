FuseSurv_Extend <- function(site,
                            eval.times,
                            IF.00, IF.01,
                            S.00, S.01,
                            Aug.00.mean, Aug.01.mean,
                            Aug.R0.mean, Aug.R1.mean,
                            Aug.R0.mean.sour, Aug.R1.mean.sour,
                            IF.R0, IF.R1, 
                            IF.CCOD.0, IF.CCOD.1, 
                            ind.R1.ccod,
                            s=4399) {
  
  set.seed(seed=s)
  
  ##################################################################
  ## ~~~~~~ Target-only and CCOD estimates for RD and RMST ~~~~~~ ##
  ##################################################################
  
  ### Risk difference
  n0 <- nrow(IF.00)
  n.all <- nrow(IF.CCOD.0)
  prop.R1 <- n0/n.all
  K <- length(unique(site))
  n.site <- table(site)
  
  IF.TGT.RD <- IF.01-IF.00
  RD.TGT <- apply(IF.TGT.RD, 2, mean)
  RD.TGT.sd <- apply(IF.TGT.RD, 2, sd)/sqrt(n0)
  df.RD.TGT <- data.frame(time=eval.times, RD=RD.TGT, sd=RD.TGT.sd)
  
  IF.CCOD.RD <- IF.CCOD.1-IF.CCOD.0
  RD.CCOD <- apply(IF.CCOD.RD, 2, mean)
  RD.CCOD.sd <- sqrt((apply(IF.CCOD.RD[ind.R1.ccod,],2,var)*prop.R1 + 
                        apply(IF.CCOD.RD[-ind.R1.ccod,],2,var)*(1-prop.R1)) /n.all)
  df.RD.CCOD <- data.frame(time=eval.times, RD=RD.CCOD, sd=RD.CCOD.sd)
  
  ### Survival ratio
  IF.00.center <- sweep(IF.00, 2, colMeans(IF.00), `-`)
  IF.01.center <- sweep(IF.01, 2, colMeans(IF.01), `-`)
  S0_hat <- colMeans(IF.00)
  S1_hat <- colMeans(IF.01)
  
  IF.TGT.SR <- sweep(IF.01.center,2,1/S0_hat, `*`)-sweep(IF.00.center,2,(S1_hat/S0_hat^2), `*`)
  SR.TGT <- S1_hat/S0_hat
  SR.TGT.sd <- apply(IF.TGT.SR, 2, sd) / sqrt(n0)
  df.SR.TGT <- data.frame(time=eval.times, SR=SR.TGT, sd=SR.TGT.sd)
  
  IF.CCOD.0.center <- sweep(IF.CCOD.0, 2, colMeans(IF.CCOD.0), `-`)
  IF.CCOD.1.center <- sweep(IF.CCOD.1, 2, colMeans(IF.CCOD.1), `-`)
  S0.ccod_hat <- colMeans(IF.CCOD.0)
  S1.ccod_hat <- colMeans(IF.CCOD.1)
  
  IF.CCOD.SR <- sweep(IF.CCOD.1.center,2,1/S0.ccod_hat, `*`) - sweep(IF.CCOD.0.center,2,(S1.ccod_hat/S0.ccod_hat^2), `*`)
  SR.CCOD <- S1.ccod_hat/S0.ccod_hat
  SR.CCOD.sd <- sqrt((apply(IF.CCOD.SR[ind.R1.ccod,],2,var)*prop.R1 + 
                        apply(IF.CCOD.SR[-ind.R1.ccod,],2,var)*(1-prop.R1)) /n.all)
  df.SR.CCOD <- data.frame(time=eval.times, SR=SR.CCOD, sd=SR.CCOD.sd)
  
  ### RMST (by treatment arm and difference)
  dt <- diff(c(0, eval.times))
  
  IF.TGT.RMST.0 <- IF.00 %*% dt
  RMST.0.TGT <- mean(IF.TGT.RMST.0)
  RMST.0.TGT.sd <- sd(IF.TGT.RMST.0) / sqrt(n0)
  df.RMST.0.TGT <- data.frame(time.max=max(eval.times), RMST=RMST.0.TGT, sd=RMST.0.TGT.sd)
  
  IF.TGT.RMST.1 <- IF.01 %*% dt
  RMST.1.TGT <- mean(IF.TGT.RMST.1)
  RMST.1.TGT.sd <- sd(IF.TGT.RMST.1) / sqrt(n0)
  df.RMST.1.TGT <- data.frame(time.max=max(eval.times), RMST=RMST.1.TGT, sd=RMST.1.TGT.sd)
  
  IF.TGT.RMST.diff <- IF.TGT.RD %*% dt
  RMST.diff.TGT <- mean(IF.TGT.RMST.diff)
  RMST.diff.TGT.sd <- sd(IF.TGT.RMST.diff) / sqrt(n0)
  df.RMST.diff.TGT <- data.frame(time.max=max(eval.times), RMST=RMST.diff.TGT, sd=RMST.diff.TGT.sd)
  
  IF.CCOD.RMST.0 <- IF.CCOD.0 %*% dt
  RMST.0.CCOD <- mean(IF.CCOD.RMST.0)
  RMST.0.CCOD.sd <-  sqrt((var(IF.CCOD.RMST.0[ind.R1.ccod])*prop.R1 +
                             var(IF.CCOD.RMST.0[-ind.R1.ccod])*(1-prop.R1)) /n.all)
  df.RMST.0.CCOD <- data.frame(time.max=max(eval.times), RMST=RMST.0.CCOD, sd=RMST.0.CCOD.sd)
  
  IF.CCOD.RMST.1 <- IF.CCOD.1 %*% dt
  RMST.1.CCOD <- mean(IF.CCOD.RMST.1)
  RMST.1.CCOD.sd <- sqrt((var(IF.CCOD.RMST.1[ind.R1.ccod])*prop.R1 + 
                            var(IF.CCOD.RMST.1[-ind.R1.ccod])*(1-prop.R1)) /n.all)
  df.RMST.1.CCOD <- data.frame(time.max=max(eval.times), RMST=RMST.1.CCOD, sd=RMST.1.CCOD.sd)
  
  IF.CCOD.RMST.diff <- IF.CCOD.RD %*% dt
  RMST.diff.CCOD <- mean(IF.CCOD.RMST.diff)
  RMST.diff.CCOD.sd <- sqrt((var(IF.CCOD.RMST.diff[ind.R1.ccod])*prop.R1 + 
                               var(IF.CCOD.RMST.diff[-ind.R1.ccod])*(1-prop.R1)) /n.all)
  df.RMST.diff.CCOD <- data.frame(time.max=max(eval.times), RMST=RMST.diff.CCOD, sd=RMST.diff.CCOD.sd)
  
  ##################################################################
  ## ~~~~~~ Federated weighting estimate for RD, RR, RMST ~~~~~~~ ##
  ##################################################################
  
  ### Survival difference
  N.time <- length(eval.times)
  Aug.TGT.mean <- Aug.01.mean-Aug.00.mean
  Aug.R.mean <- Aug.R1.mean-Aug.R0.mean
  Aug.R.mean.sour <- Aug.R1.mean.sour-Aug.R0.mean.sour
  wt.RD <- chi.RD <- augdiff.RD <- matrix(NA, nrow=N.time, ncol=K-1)
  for(i in 1:N.time) {
    IF.TGT.RD.center <- c(IF.TGT.RD[,i]-mean(IF.TGT.RD[,i]), rep(0, n.all-n0))
    IF.RD.diff <- matrix(0, ncol=K-1, nrow=n.all)
    ind0 <- which(site==0)
    
    for(r in 1:(K-1)) {
      IF.RD.diff[,r][ind0] <- IF.TGT.RD.center[ind0]
      ind <- which(site==r)
      IF.RD.diff[,r][ind] <- -(IF.R1[[r]][,i]-IF.R0[[r]][,i]) 
      chi.RD[i,r] <- Aug.TGT.mean[i]-Aug.R.mean[i,r]
      augdiff.RD[i,r] <- Aug.TGT.mean[i]-Aug.R.mean.sour[i,r]
    } 
    cv.fit=try(cv.glmnet(x=IF.RD.diff, y=IF.TGT.RD.center))
    if(class(cv.fit)[1]!="try-error") {
      fit0=try(glmnet(x=IF.RD.diff, y=IF.TGT.RD.center, 
                      penalty.factor=chi.RD[i,]^2,
                      intercept=FALSE,
                      alpha=1,
                      lambda=cv.fit$lambda.1se,
                      lower.limits=0,
                      upper.limits=1))
      if(class(fit0)[1]!="try-error") { 
        wt.RD[i,]=coef(fit0, s=cv.fit$lambda.1se)[-1] } else { wt.RD[i,]=rep(0, K-1) }
    } else { wt.RD[i,]=rep(0, K-1) }
  } 
  wt.RD.tgt <- 1-apply(wt.RD,1,sum)
  RD.FED <- apply(augdiff.RD*wt.RD, 1, sum) + RD.TGT
  all.var <- (apply(IF.TGT.RD,2,var)*(wt.RD.tgt^2+2*wt.RD.tgt*(1-wt.RD.tgt)) + 
                apply(S.01-S.00,2,var)*(1-wt.RD.tgt)^2) / n.site[1] 
  for(k in 1:(K-1)) {
    all.var <- all.var + apply(IF.R1[[k]]-IF.R0[[k]], 2, var)*wt.RD[,k]^2 / n.site[k+1]
  }
  RD.FED.sd <- sqrt(all.var) 
  df.RD.FED <- data.frame(time=eval.times, RD=RD.FED, sd=RD.FED.sd)
  
  
  ### Survival ratio
  N.time <- length(eval.times)
  Aug.SR.TGT.mean <- (Aug.01.mean/S0_hat) - Aug.00.mean*(S1_hat/S0_hat^2)
  Aug.SR.R.mean <- sweep(Aug.R1.mean,1,1/S0_hat,`*`)-sweep(Aug.R0.mean,1,S1_hat/(S0_hat^2),`*`)
  Aug.SR.R.mean.sour <-sweep(Aug.R1.mean.sour,1,1/S0_hat,`*`)-sweep(Aug.R0.mean.sour,1,S1_hat/(S0_hat^2),`*`)
  wt.SR <- chi.SR <- augdiff.SR <- matrix(NA, nrow=N.time, ncol=K-1)
  for(i in 1:N.time) {
    IF.TGT.SR.center <- c(IF.TGT.SR[,i]-mean(IF.TGT.SR[,i]), rep(0, n.all-n0))
    IF.SR.diff <- matrix(0, ncol=K-1, nrow=n.all)
    ind0 <- which(site==0)
    
    for(r in 1:(K-1)) {
      IF.SR.diff[,r][ind0] <- IF.TGT.SR.center[ind0]
      ind <- which(site==r)
      IF.SR.diff[,r][ind] <- -(IF.R1[[r]][,i]/S0_hat[i]-IF.R0[[r]][,i]*(S1_hat[i]/S0_hat[i]^2)) 
      chi.SR[i,r] <- Aug.SR.TGT.mean[i]-Aug.SR.R.mean[i,r]
      augdiff.SR[i,r] <- Aug.SR.TGT.mean[i]-Aug.SR.R.mean.sour[i,r]
    } 
    cv.fit=try(cv.glmnet(x=IF.SR.diff, y=IF.TGT.SR.center))
    if(class(cv.fit)[1]!="try-error") {
      fit0=try(glmnet(x=IF.SR.diff, y=IF.TGT.SR.center, 
                      penalty.factor=chi.SR[i,]^2,
                      intercept=FALSE,
                      alpha=1,
                      lambda=cv.fit$lambda.1se,
                      lower.limits=0,
                      upper.limits=1))
      if(class(fit0)[1]!="try-error") { 
        wt.SR[i,]=coef(fit0, s=cv.fit$lambda.1se)[-1] } else { wt.SR[i,]=rep(0, K-1) }
    } else { wt.SR[i,]=rep(0, K-1) }
  } 
  wt.SR.tgt <- 1-apply(wt.SR,1,sum)
  SR.FED <- apply(augdiff.SR*wt.SR, 1, sum) + SR.TGT
  all.var <- (apply(IF.TGT.SR,2,var)*(wt.SR.tgt^2+2*wt.SR.tgt*(1-wt.SR.tgt)) + 
                apply(sweep(Aug.R1.mean,1,1/S0_hat,`*`)-sweep(Aug.R0.mean,1,S1_hat/(S0_hat^2),`*`),2,var)*(1-wt.SR.tgt)^2) / n.site[1] 
  for(k in 1:(K-1)) {
    all.var <- all.var + apply(IF.R1[[k]]-IF.R0[[k]], 2, var)*wt.SR[,k]^2 / n.site[k+1]
  }
  SR.FED.sd <- sqrt(all.var) 
  df.SR.FED <- data.frame(time=eval.times, SR=SR.FED, sd=SR.FED.sd)
  
  
  ### Restrict mean survival times
  Aug.RMST.0.TGT.mean <- sum(Aug.00.mean*dt)
  Aug.RMST.1.TGT.mean <- sum(Aug.01.mean*dt)
  Aug.RMST.0.R.mean <- t(t(Aug.R0.mean) %*% dt)
  Aug.RMST.1.R.mean <- t(t(Aug.R1.mean) %*% dt)
  Aug.RMST.0.R.mean.sour <- t(t(Aug.R0.mean.sour) %*% dt)
  Aug.RMST.1.R.mean.sour <- t(t(Aug.R1.mean.sour) %*% dt)
  wt.RMST.0 <- chi.RMST.0 <- augdiff.RMST.0 <- wt.RMST.1 <- chi.RMST.1 <- augdiff.RMST.1 <- 
    wt.RMST <- chi.RMST <- augdiff.RMST <- matrix(NA, nrow=1, ncol=K-1)
  IF.TGT.RMST.0.center <- c(IF.TGT.RMST.0-mean(IF.TGT.RMST.0), rep(0, n.all-n0))
  IF.TGT.RMST.1.center <- c(IF.TGT.RMST.1-mean(IF.TGT.RMST.1), rep(0, n.all-n0))
  IF.RMST.0.diff <- IF.RMST.1.diff <- matrix(0, ncol=K-1, nrow=n.all)
  
  ind0 <- which(site==0)
  for(r in 1:(K-1)) {
    IF.RMST.0.diff[,r][ind0] <- IF.TGT.RMST.0.center[ind0]
    IF.RMST.1.diff[,r][ind0] <- IF.TGT.RMST.1.center[ind0]
    
    ind <- which(site==r)
    IF.RMST.0.diff[,r][ind] <- -IF.R0[[r]] %*% dt
    IF.RMST.1.diff[,r][ind] <- -IF.R1[[r]] %*% dt
    
    chi.RMST.0[,r] <- Aug.RMST.0.TGT.mean-Aug.RMST.0.R.mean[,r]
    chi.RMST.1[,r] <- Aug.RMST.1.TGT.mean-Aug.RMST.1.R.mean[,r]
    chi.RMST[,r] <- chi.RMST.1[,r] - chi.RMST.0[,r]
    
    augdiff.RMST.0[,r] <- Aug.RMST.0.TGT.mean-Aug.RMST.0.R.mean.sour[,r]
    augdiff.RMST.1[,r] <- Aug.RMST.1.TGT.mean-Aug.RMST.1.R.mean.sour[,r]
    augdiff.RMST[,r] <- augdiff.RMST.1[,r] - augdiff.RMST.0[,r]
  }
  IF.TGT.RMST.center <- IF.TGT.RMST.1.center - IF.TGT.RMST.0.center
  IF.RMST.diff <- IF.RMST.1.diff - IF.RMST.0.diff
  
  cv.fit0=try(cv.glmnet(x=IF.RMST.0.diff, y=IF.TGT.RMST.0.center))
  if(class(cv.fit0)[1]!="try-error") {
    fit0=try(glmnet(x=IF.RMST.0.diff, y=IF.TGT.RMST.0.center, 
                    penalty.factor=chi.RMST.0^2,
                    intercept=FALSE,
                    alpha=1,
                    lambda=cv.fit0$lambda.1se,
                    lower.limits=0,
                    upper.limits=1))
    if(class(fit0)[1]!="try-error") { 
      wt.RMST.0=coef(fit0, s=cv.fit0$lambda.1se)[-1] } else { wt.RMST.0=rep(0, K-1) }
  } else { wt.RMST.0=rep(0, K-1) }
  
  cv.fit1=try(cv.glmnet(x=IF.RMST.1.diff, y=IF.TGT.RMST.1.center))
  if(class(cv.fit1)[1]!="try-error") {
    fit1=try(glmnet(x=IF.RMST.1.diff, y=IF.TGT.RMST.1.center, 
                    penalty.factor=chi.RMST.1^2,
                    intercept=FALSE,
                    alpha=1,
                    lambda=cv.fit1$lambda.1se,
                    lower.limits=0,
                    upper.limits=1))
    if(class(fit1)[1]!="try-error") { 
      wt.RMST.1=coef(fit1, s=cv.fit1$lambda.1se)[-1] } else { wt.RMST.1=rep(0, K-1) }
  } else { wt.RMST.1=rep(0, K-1) }
  
  cv.fit=try(cv.glmnet(x=IF.RMST.diff, y=IF.TGT.RMST.center))
  if(class(cv.fit)[1]!="try-error") {
    fit=try(glmnet(x=IF.RMST.diff, y=IF.TGT.RMST.center, 
                    penalty.factor=chi.RMST^2,
                    intercept=FALSE,
                    alpha=1,
                    lambda=cv.fit$lambda.1se,
                    lower.limits=0,
                    upper.limits=1))
    if(class(fit)[1]!="try-error") { 
      wt.RMST=coef(fit, s=cv.fit$lambda.1se)[-1] } else { wt.RMST=rep(0, K-1) }
  } else { wt.RMST=rep(0, K-1) }
  
  wt.RMST.0.tgt <- 1-sum(wt.RMST.0)
  wt.RMST.1.tgt <- 1-sum(wt.RMST.1)
  wt.RMST.tgt <- 1-sum(wt.RMST)
  
  RMST.0.FED <- sum(augdiff.RMST.0*wt.RMST.0) + RMST.0.TGT
  RMST.1.FED <- sum(augdiff.RMST.1*wt.RMST.1) + RMST.1.TGT
  RMST.diff.FED <- sum(augdiff.RMST*wt.RMST) + RMST.diff.TGT
  
  all.var0 <- (var(IF.TGT.RMST.0)*(wt.RMST.0.tgt^2+2*wt.RMST.0.tgt*(1-wt.RMST.0.tgt)) +
                var(S.00%*%dt)*(1-wt.RMST.0.tgt)^2) / n.site[1] 
  all.var1 <- (var(IF.TGT.RMST.1)*(wt.RMST.1.tgt^2+2*wt.RMST.1.tgt*(1-wt.RMST.1.tgt)) +
                var(S.01%*%dt)*(1-wt.RMST.1.tgt)^2) / n.site[1] 
  all.var <- (var(IF.TGT.RMST.diff)*(wt.RMST.tgt^2+2*wt.RMST.tgt*(1-wt.RMST.tgt)) +
                var((S.01-S.00)%*%dt)*(1-wt.RMST.tgt)^2) / n.site[1] 
  for(k in 1:(K-1)) {
    all.var0 <- all.var0 + var(IF.R0[[k]]%*%dt)*wt.RMST.0[k]^2 / n.site[k+1]
    all.var1 <- all.var1 + var(IF.R1[[k]]%*%dt)*wt.RMST.1[k]^2 / n.site[k+1]
    all.var <- all.var + var((IF.R1[[k]]-IF.R0[[k]])%*%dt)*wt.RMST[k]^2 / n.site[k+1]
  }
  RMST.0.FED.sd <- sqrt(all.var0) 
  RMST.1.FED.sd <- sqrt(all.var1) 
  RMST.diff.FED.sd <- sqrt(all.var) 
  
  df.RMST.0.FED <- data.frame(time.max=max(eval.times), RMST=RMST.0.FED, sd=RMST.0.FED.sd)
  df.RMST.1.FED <- data.frame(time.max=max(eval.times), RMST=RMST.1.FED, sd=RMST.1.FED.sd)
  df.RMST.diff.FED <- data.frame(time.max=max(eval.times), RMST=RMST.diff.FED, sd=RMST.diff.FED.sd)
  
  return(list(df.RD.TGT=df.RD.TGT, 
              df.RD.CCOD=df.RD.CCOD, 
              df.RD.FED=df.RD.FED,
              
              df.SR.TGT=df.SR.TGT, 
              df.SR.CCOD=df.SR.CCOD, 
              df.SR.FED=df.SR.FED,
              
              df.RMST.0.TGT=df.RMST.0.TGT, 
              df.RMST.1.TGT=df.RMST.1.TGT, 
              df.RMST.diff.TGT=df.RMST.diff.TGT, 
              
              df.RMST.0.CCOD=df.RMST.0.CCOD,
              df.RMST.1.CCOD=df.RMST.1.CCOD, 
              df.RMST.diff.CCOD=df.RMST.diff.CCOD, 
              
              df.RMST.0.FED=df.RMST.0.FED,
              df.RMST.1.FED=df.RMST.1.FED,
              df.RMST.diff.FED=df.RMST.diff.FED))
}
