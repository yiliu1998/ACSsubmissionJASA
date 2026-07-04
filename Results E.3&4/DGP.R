### ---  true survival values at 30, 60 and 90 days
N <- 10^8
r <- 0
X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r
X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)

rho <- 1.2
lambda <- 0.6
delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2)
ph.t.0 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2))
ph.t.1 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) + delta.t)

t <- c(30,60,90)
S0.true <- sapply(t, function(time) mean(exp(-ph.t.0*lambda*time^rho)))
S1.true <- sapply(t, function(time) mean(exp(-ph.t.1*lambda*time^rho)))
S0.true
S1.true

save(file="truth.Rdata", S1.true, S0.true, eval.times=t)

### --- data generation process for the main setting (Section E.3)
DGP <- function(M=600, N=rep(1000, 5), case="homo", s=311) {
  expit <- function(x) 1/(1+exp(-x))
  N.site <- length(N)
  site <- list()
  for(i in 1:N.site) { site[[i]] <- rep(i-1, N[i]) }
  site <- unlist(site)
  R <- max(site)
  dat <- list()
  
  set.seed(s)
  for(i in 1:M) {
    r <- 0
    X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
    X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
    X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
    
    if(case=="homo") {
      D.C <- D.T <- 0
    }
    if(case=="diffX") {
      for(r in 2:R) {
        X1[site==r] <- x1 <- 33*rbeta(N[r+1], shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
        X2[site==r] <- 52*rbeta(N[r+1], shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
        X3[site==r] <- (4+2*r)*rbeta(N[r+1], shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
      }
      D.C <- D.T <- 0
    }
    if(case=="diffT") {
      D.T <- site
      D.C <- 0
    }
    if(case=="diffC") {
      D.T <- 0
      D.C <- site
    }
    if(case=="diffAll") {
      for(r in 2:R) {
        X1[site==r] <- x1 <- 33*rbeta(N[r+1], shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
        X2[site==r] <- 52*rbeta(N[r+1], shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
        X3[site==r] <- (4+2*r)*rbeta(N[r+1], shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
      }
      D.C <- D.T <- site
    }
    
    eta.src <- -1.05+log(1.3+exp(-12+X1/10)+exp(-2+X3/3)+exp(-2+X2/12))
    eta.lim <- -1.05+log(0.3+exp(-120+X1)+exp(-6+X3)+exp(-6+X2/4))
    eta.tgt <- 0.7*eta.src+0.3*eta.lim
    
    g0s <- expit(ifelse(site==0, eta.tgt, eta.src))
    A <- rbinom(sum(N), size=1, prob=g0s)
    
    delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2) + D.T*0.02*(X1+X3-25)
    ph.t <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) - D.T*0.03*(X3-X1+20) + A*delta.t)
    
    delta.c <- -0.4+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X1+X2+X3-50)
    ph.c <- exp(-4.87+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X2-25) + A*delta.c)
    
    u1 <- runif(sum(N))
    u2 <- runif(sum(N))
    rho <- 1.2
    lambda <- 0.6
    event.time <- (-log(u1)/(ph.t * lambda))^(1/rho)
    cens.time <- (-log(u2)/(ph.c * lambda))^(1/rho)
    cens.time[cens.time > 200] <- 200
    obs.time <- pmin(event.time, cens.time)
    obs.event <- as.numeric(event.time <= cens.time)
    
    dat[[i]] <- data.frame(site=site, X1=X1, X2=X2, X3=X3, A=A, Y=obs.time, Delta=obs.event)
  }
  return(dat=dat)
}

dat.homo <- DGP(case="homo", N=c(300, rep(300, 4)))
dat.diffX <- DGP(case="diffX", N=c(300, rep(300, 4)))
dat.diffT <- DGP(case="diffT", N=c(300, rep(300, 4)))
dat.diffC <- DGP(case="diffC", N=c(300, rep(300, 4)))
dat.diffAll <- DGP(case="diffAll", N=c(300, rep(300, 4)))
save(file="obsdata_s.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)

dat.homo <- DGP(case="homo", N=c(300, rep(600, 4)))
dat.diffX <- DGP(case="diffX", N=c(300, rep(600, 4)))
dat.diffT <- DGP(case="diffT", N=c(300, rep(600, 4)))
dat.diffC <- DGP(case="diffC", N=c(300, rep(600, 4)))
dat.diffAll <- DGP(case="diffAll", N=c(300, rep(600, 4)))
save(file="obsdata_l.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)

dat.homo <- DGP(case="homo", N=c(300, rep(1000, 4)))
dat.diffX <- DGP(case="diffX", N=c(300, rep(1000, 4)))
dat.diffT <- DGP(case="diffT", N=c(300, rep(1000, 4)))
dat.diffC <- DGP(case="diffC", N=c(300, rep(1000, 4)))
dat.diffAll <- DGP(case="diffAll", N=c(300, rep(1000, 4)))
save(file="obsdata_l2.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)

### --- data generation process for the limited overlap setting (Section E.4)
DGP_limO <- function(M=600, N=rep(1000, 5), case="homo", s=311) {
  expit <- function(x) 1/(1+exp(-x))
  N.site <- length(N)
  site <- list()
  for(i in 1:N.site) { site[[i]] <- rep(i-1, N[i]) }
  site <- unlist(site)
  R <- max(site)
  dat <- list()
  
  set.seed(s)
  for(i in 1:M) {
    r <- 0
    X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
    X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
    X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
    
    if(case=="homo") {
      D.C <- D.T <- 0
    }
    if(case=="diffX") {
      for(r in 2:R) {
        X1[site==r] <- x1 <- 33*rbeta(N[r+1], shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
        X2[site==r] <- 52*rbeta(N[r+1], shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
        X3[site==r] <- (4+2*r)*rbeta(N[r+1], shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
      }
      D.C <- D.T <- 0
    }
    if(case=="diffT") {
      D.T <- site
      D.C <- 0
    }
    if(case=="diffC") {
      D.T <- 0
      D.C <- site
    }
    if(case=="diffAll") {
      for(r in 2:R) {
        X1[site==r] <- x1 <- 33*rbeta(N[r+1], shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
        X2[site==r] <- 52*rbeta(N[r+1], shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
        X3[site==r] <- (4+2*r)*rbeta(N[r+1], shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
      }
      D.C <- D.T <- site
    }
    
    logit.g0s <- -1.05+log( 0.3+exp(-120+X1)+exp(-6+X3)+exp(-6+X2/4) )
    g0s <- expit(logit.g0s)
    A.tgt <- rbinom(sum(N), size=1, prob=g0s)
    
    logit.g0s <- -1.05+log( 1.3+exp(-12+X1/10)+exp(-2+X3/3)+exp(-2+X2/12) )
    g0s <- expit(logit.g0s)
    A.src <- rbinom(sum(N), size=1, prob=g0s)
    
    A <- A.tgt*I(site==0) + A.src*I(site!=0)
    
    delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2) + D.T*0.02*(X1+X3-25)
    ph.t <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) - D.T*0.03*(X3-X1+20) + A*delta.t)
    
    delta.c <- -0.4+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X1+X2+X3-50)
    ph.c <- exp(-4.87+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X2-25) + A*delta.c)
    
    u1 <- runif(sum(N))
    u2 <- runif(sum(N))
    rho <- 1.2
    lambda <- 0.6
    event.time <- (-log(u1)/(ph.t * lambda))^(1/rho)
    cens.time <- (-log(u2)/(ph.c * lambda))^(1/rho)
    cens.time[cens.time > 200] <- 200
    obs.time <- pmin(event.time, cens.time)
    obs.event <- as.numeric(event.time <= cens.time)
    
    dat[[i]] <- data.frame(site=site, X1=X1, X2=X2, X3=X3, A=A, Y=obs.time, Delta=obs.event)
  }
  return(dat=dat)
}

dat.homo <- DGP_limO(case="homo", N=c(300, rep(300, 4)))
dat.diffX <- DGP_limO(case="diffX", N=c(300, rep(300, 4)))
dat.diffT <- DGP_limO(case="diffT", N=c(300, rep(300, 4)))
dat.diffC <- DGP_limO(case="diffC", N=c(300, rep(300, 4)))
dat.diffAll <- DGP_limO(case="diffAll", N=c(300, rep(300, 4)))
save(file="obsdata_limO.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)

### --- run clustered Cox simulation results (takes only 10 mins or so)
source("ClusterCox.R")
