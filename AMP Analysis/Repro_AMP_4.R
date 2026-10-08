### ----------------------------------------------------------------------------
### Reproducibility for additional SA-target diagnostics
### --- This script reproduces Tables A.12, A.13 and A.14 in
### --- Online Supplemental Material Appendix A.
### ----------------------------------------------------------------------------

library(glmnet)
library(dplyr)

load("result_main_SA.Rdata")

time.sel <- c(148, 330, 512)
ind.sel <- match(time.sel, result.SA$eval.times)
src <- c("OA", "BP", "US")


### ----------------------------------------------------------------------------
### 1. Selected federated weights and chi: Table A.12
### ----------------------------------------------------------------------------

### Extract the selected FED weights at the reported evaluation days.
### weights = cbind(control target, control sources,
###                 treated target, treated sources)

tab.weights <- rbind(
  data.frame(Day=time.sel, Arm="Control",
             SA=result.SA$weights[ind.sel,1],
             OA=result.SA$weights[ind.sel,2],
             BP=result.SA$weights[ind.sel,3],
             US=result.SA$weights[ind.sel,4]),
  data.frame(Day=time.sel, Arm="Treated",
             SA=result.SA$weights[ind.sel,5],
             OA=result.SA$weights[ind.sel,6],
             BP=result.SA$weights[ind.sel,7],
             US=result.SA$weights[ind.sel,8])
)

### Extract the corresponding estimated source-target discrepancy terms.
### chi = cbind(chi0, chi1), where each block contains OA, BP, US

tab.chi <- rbind(
  data.frame(Day=time.sel, Arm="Control",
             OA=result.SA$chi[ind.sel,1],
             BP=result.SA$chi[ind.sel,2],
             US=result.SA$chi[ind.sel,3]),
  data.frame(Day=time.sel, Arm="Treated",
             OA=result.SA$chi[ind.sel,4],
             BP=result.SA$chi[ind.sel,5],
             US=result.SA$chi[ind.sel,6])
)

cat("\nSelected FED weights:\n")
print(tab.weights, row.names=FALSE)

cat("\nEstimated chi discrepancies:\n")
print(tab.chi, row.names=FALSE)


### ----------------------------------------------------------------------------
### 2. Sensitivity to lambda: Table A.13
### ----------------------------------------------------------------------------

### Construct the response, source IF contrasts, discrepancy terms, and
### augmentation differences used in the FED weighting problem.
get_input <- function(res, i, arm=0) {
  site <- res$site
  ind0 <- which(site==0)
  n <- length(site)
  
  IF.tgt <- if(arm==0) res$IF.00[,i] else res$IF.01[,i]
  IF.src <- if(arm==0) res$IF.R0 else res$IF.R1
  chi <- if(arm==0) res$chi[i,1:3] else res$chi[i,4:6]
  
  y <- c(IF.tgt-mean(IF.tgt), rep(0,n-length(IF.tgt)))
  X <- matrix(0,nrow=n,ncol=3)
  
  for(k in 1:3) {
    X[ind0,k] <- y[ind0]
    indk <- which(site==k)
    X[indk,k] <- -IF.src[[k]][,i]
  }
  
  aug.tgt <- if(arm==0) res$Aug.00.mean[i] else res$Aug.01.mean[i]
  aug.src <- if(arm==0) res$Aug.R0.mean.sour[i,] else res$Aug.R1.mean.sour[i,]
  theta <- if(arm==0) res$df.TGT$surv0[i] else res$df.TGT$surv1[i]
  
  list(X=X, y=y, chi=chi, augdiff=aug.tgt-aug.src, theta=theta)
}


### Refit the FED weighting step for a specified lambda value.
### If exclude=TRUE, the BP and US source weights are constrained to zero.
fit_fed <- function(res, i, arm, lambda, exclude=FALSE) {
  z <- get_input(res,i,arm)
  
  upper <- if(exclude) c(1,0,0) else c(1,1,1)
  
  fit <- glmnet(z$X,z$y,penalty.factor=z$chi^2,intercept=FALSE,alpha=1,
                lambda=lambda,lower.limits=0,upper.limits=upper)
  
  w <- as.numeric(coef(fit,s=lambda)[-1])
  wtgt <- 1-sum(w)
  est <- z$theta + sum(z$augdiff*w)
  
  ## Same variance calculation as TrtSurvCurves()
  n.site <- table(res$site)
  
  if(arm==0) {
    v <- (var(res$IF.00[,i])*(wtgt^2+2*wtgt*(1-wtgt)) +
            var(res$S.00[,i])*(1-wtgt)^2)/n.site[1]
    for(k in 1:3) v <- v + var(res$IF.R0[[k]][,i])*w[k]^2/n.site[k+1]
  } else {
    v <- (var(res$IF.01[,i])*(wtgt^2+2*wtgt*(1-wtgt)) +
            var(res$S.01[,i])*(1-wtgt)^2)/n.site[1]
    for(k in 1:3) v <- v + var(res$IF.R1[[k]][,i])*w[k]^2/n.site[k+1]
  }
  
  c(Est=unname(est), SE=unname(sqrt(v)),
    SA=unname(wtgt), OA=unname(w[1]), BP=unname(w[2]), US=unname(w[3]))
}


### --- Recover the original CV-selected lambda values used in the SA analysis

K <- length(unique(result.SA$site))

set.seed(2388)
seeds <- round(runif(20*K,0,20e5))
set.seed(seeds[K+5])

lambda0 <- lambda1 <- rep(NA,length(result.SA$eval.times))

for(i in seq_along(result.SA$eval.times)) {
  z0 <- get_input(result.SA,i,0)
  cv0 <- try(cv.glmnet(x=z0$X,y=z0$y),silent=TRUE)
  if(!inherits(cv0,"try-error")) lambda0[i] <- cv0$lambda.1se
  
  z1 <- get_input(result.SA,i,1)
  cv1 <- try(cv.glmnet(x=z1$X,y=z1$y),silent=TRUE)
  if(!inherits(cv1,"try-error")) lambda1[i] <- cv1$lambda.1se
}


### --- Evaluate FED estimates at 0.5x, 1x, and 2x the CV-selected lambda

lambda.mult <- c(0.5,1,2)
out.lambda <- list(); cc <- 1

for(i in ind.sel) {
  for(a in 0:1) {
    lambda.cv <- if(a==0) lambda0[i] else lambda1[i]
    
    for(m in lambda.mult) {
      x <- fit_fed(result.SA,i,a,m*lambda.cv)
      
      out.lambda[[cc]] <- data.frame(
        Day=result.SA$eval.times[i],
        Arm=ifelse(a==0,"Control","Treated"),
        Lambda.mult=m,
        Est=x["Est"], SE=x["SE"],
        SA=x["SA"], OA=x["OA"], BP=x["BP"], US=x["US"]
      )
      cc <- cc+1
    }
  }
}

tab.lambda <- do.call(rbind,out.lambda)
tab.lambda$Lower <- tab.lambda$Est-1.96*tab.lambda$SE
tab.lambda$Upper <- tab.lambda$Est+1.96*tab.lambda$SE

print(round(tab.lambda[,c("Day","Lambda.mult","Est","SE","Lower","Upper",
                          "SA","OA","BP","US")],4),row.names=FALSE)


### ----------------------------------------------------------------------------
### 3. FED sensitivity excluding BP and US: Table A.14
### ----------------------------------------------------------------------------
### --- Compare the original FED estimator with a version borrowing only from OA

out.exclude <- list(); cc <- 1

for(i in ind.sel) {
  for(a in 0:1) {
    lam <- if(a==0) lambda0[i] else lambda1[i]
    
    x.full <- fit_fed(result.SA,i,a,lam,exclude=FALSE)
    x.OA   <- fit_fed(result.SA,i,a,lam,exclude=TRUE)
    
    out.exclude[[cc]] <- data.frame(
      Day=result.SA$eval.times[i],
      Arm=ifelse(a==0,"Control","Treated"),
      FED=x.full["Est"],
      SE.FED=x.full["SE"],
      FED.OAonly=x.OA["Est"],
      SE.OAonly=x.OA["SE"],
      SA.weight=x.OA["SA"],
      OA.weight=x.OA["OA"]
    )
    cc <- cc+1
  }
}

tab.exclude <- do.call(rbind,out.exclude)
print(tab.exclude,row.names=FALSE)
