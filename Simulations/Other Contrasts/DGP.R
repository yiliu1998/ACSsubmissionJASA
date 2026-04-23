DGP <- function(M=520, N=rep(1000, 5), case="homo", s=1245) {
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
    logit.g0s <- -1.05+log( 1.3+exp(-12+X1/10)+exp(-2+X3/3)+exp(-2+X2/12) )
    g0s <- expit(logit.g0s)
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

dat.homo <- DGP(case="homo", N=c(300, rep(600, 4)))
dat.diffX <- DGP(case="diffX", N=c(300, rep(600, 4)))
dat.diffT <- DGP(case="diffT", N=c(300, rep(600, 4)))
dat.diffC <- DGP(case="diffC", N=c(300, rep(600, 4)))
dat.diffAll <- DGP(case="diffAll", N=c(300, rep(600, 4)))
save(file="obsdata_main.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
