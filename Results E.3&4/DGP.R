# ============================================================
# Reproducibility code for simulation data generation
#
# This script:
#   (1) computes Monte Carlo approximations to the true target-site
#       survival probabilities at days 30, 60, and 90;
#   (2) generates the main simulation datasets used in Section E.3;
#   (3) generates the limited-overlap simulation datasets used in
#       Section E.4; and
#   (4) runs the clustered Cox benchmark through ClusterCox.R.
#
# Saved objects:
#   truth.Rdata        : true target-site survival probabilities
#   obsdata_s.Rdata    : main setting, source-site n = 300
#   obsdata_l.Rdata    : main setting, source-site n = 600
#   obsdata_l2.Rdata   : main setting, source-site n = 1000
#   obsdata_limO.Rdata : limited-overlap setting, source n = 300
# ============================================================


### --- true survival values at 30, 60 and 90 days
# Use a large Monte Carlo sample to approximate the true target-site survival probabilities.
N <- 10^8
r <- 0

# Generate baseline covariates under the target-site distribution.
X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r
X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)

# Weibull-type event-time model parameters.
rho <- 1.2
lambda <- 0.6

# Define the heterogeneous treatment effect and treatment-specific event-time hazard factors.
delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2)
ph.t.0 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2))
ph.t.1 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) + delta.t)

# Compute true marginal survival probabilities at the three evaluation times.
t <- c(30,60,90)
S0.true <- sapply(t, function(time) mean(exp(-ph.t.0*lambda*time^rho)))
S1.true <- sapply(t, function(time) mean(exp(-ph.t.1*lambda*time^rho)))
S0.true
S1.true

# Save truth values for downstream bias, RMSE, and coverage calculations.
save(file="truth.Rdata", S1.true, S0.true, eval.times=t)


### --- data generation process for the main setting (Section E.3)
# M: number of Monte Carlo replicates.
# N: vector of site-specific sample sizes; site 0 is the target site.
# case: distribution-shift scenario.
# s: random seed.
DGP <- function(M=600, N=rep(1000, 5), case="homo", s=311) {
  
  # Logistic inverse-link function.
  expit <- function(x) 1/(1+exp(-x))
  
  # Construct the site indicator: site 0 is the target and sites 1,...,K-1 are sources.
  N.site <- length(N)
  site <- list()
  for(i in 1:N.site) { site[[i]] <- rep(i-1, N[i]) }
  site <- unlist(site)
  R <- max(site)
  
  # Store all Monte Carlo replicates.
  dat <- list()
  
  # Fix the random seed for reproducibility.
  set.seed(s)
  
  for(i in 1:M) {
    
    # --------------------------------------------------------
    # Step 1: Generate baseline covariates
    # --------------------------------------------------------
    
    # Start from the target-site covariate distribution.
    r <- 0
    X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
    X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
    X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
    
    # --------------------------------------------------------
    # Step 2: Specify the cross-site distribution-shift scenario
    #
    # "homo"    : homogeneous DGP across sites
    # "diffX"   : covariate shift
    # "diffT"   : outcome/event-time shift
    # "diffC"   : censoring shift
    # "diffAll" : covariate, outcome, and censoring shifts
    # --------------------------------------------------------
    
    if(case=="homo") {
      D.C <- D.T <- 0
    }
    
    if(case=="diffX") {
      
      # Induce covariate-distribution shifts in Source Sites 2--4 while keeping Source Site 1 aligned with the target.
      for(r in 2:R) {
        X1[site==r] <- x1 <- 33*rbeta(N[r+1], shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
        X2[site==r] <- 52*rbeta(N[r+1], shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
        X3[site==r] <- (4+2*r)*rbeta(N[r+1], shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
      }
      
      # Keep the event-time and censoring models otherwise unchanged.
      D.C <- D.T <- 0
    }
    
    if(case=="diffT") {
      # Induce progressively stronger site-specific shifts in the event-time model.
      D.T <- site
      D.C <- 0
    }
    
    if(case=="diffC") {
      # Induce progressively stronger site-specific shifts in the censoring model.
      D.T <- 0
      D.C <- site
    }
    
    if(case=="diffAll") {
      
      # Simultaneously induce covariate-distribution shifts in Source Sites 2--4.
      for(r in 2:R) {
        X1[site==r] <- x1 <- 33*rbeta(N[r+1], shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
        X2[site==r] <- 52*rbeta(N[r+1], shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
        X3[site==r] <- (4+2*r)*rbeta(N[r+1], shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
      }
      
      # Also induce site-specific shifts in both the event-time and censoring models.
      D.C <- D.T <- site
    }
    
    # --------------------------------------------------------
    # Step 3: Generate treatment assignment
    # --------------------------------------------------------
    
    # Propensity-score models for the source and target sites.
    eta.src <- -1.05 + log(1.3 + exp(-12+X1/10) + exp(-2+X3/3) + exp(-2+X2/12))
    eta.lim <- -1.05 + log(0.3 + exp(-120+X1) + exp(-6+X3) + exp(-6+X2/4))
    
    # The target-site propensity score combines the standard and more limited-overlap specifications.
    eta.tgt <- 0.7*eta.src+0.3*eta.lim
    g0s <- expit(ifelse(site==0, eta.tgt, eta.src))
    
    # Generate the binary treatment indicator.
    A <- rbinom(sum(N), size=1, prob=g0s)
    
    # --------------------------------------------------------
    # Step 4: Generate event times
    # --------------------------------------------------------
    
    # Heterogeneous treatment effect with optional site-specific outcome shift controlled by D.T.
    delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2) + D.T*0.02*(X1+X3-25)
    
    # Multiplicative event-time hazard factor.
    ph.t <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) - D.T*0.03*(X3-X1+20) + A*delta.t)
    
    # --------------------------------------------------------
    # Step 5: Generate censoring times
    # --------------------------------------------------------
    
    # Treatment effect in the censoring model with optional site-specific censoring shift controlled by D.C.
    delta.c <- -0.4+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X1+X2+X3-50)
    
    # Multiplicative censoring-time hazard factor.
    ph.c <- exp(-4.87+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X2-25) + A*delta.c)
    
    # Generate event and censoring times from Weibull-type models by inverse-transform sampling.
    u1 <- runif(sum(N))
    u2 <- runif(sum(N))
    rho <- 1.2
    lambda <- 0.6
    event.time <- (-log(u1)/(ph.t * lambda))^(1/rho)
    cens.time <- (-log(u2)/(ph.c * lambda))^(1/rho)
    
    # Apply administrative censoring at day 200.
    cens.time[cens.time > 200] <- 200
    
    # Construct the observed follow-up time and event indicator.
    obs.time <- pmin(event.time, cens.time)
    obs.event <- as.numeric(event.time <= cens.time)
    
    # Store one simulated dataset.
    dat[[i]] <- data.frame(site=site, X1=X1, X2=X2, X3=X3, A=A, Y=obs.time, Delta=obs.event)
  }
  
  return(dat=dat)
}


### --- generate the main simulation datasets under three source-site sample-size configurations
# The target-site sample size is fixed at n=300. The four source-site sample sizes are varied across 300, 600, and 1000.
# Equal target/source sample sizes: n=300 per site.
dat.homo <- DGP(case="homo", N=c(300, rep(300, 4)))
dat.diffX <- DGP(case="diffX", N=c(300, rep(300, 4)))
dat.diffT <- DGP(case="diffT", N=c(300, rep(300, 4)))
dat.diffC <- DGP(case="diffC", N=c(300, rep(300, 4)))
dat.diffAll <- DGP(case="diffAll", N=c(300, rep(300, 4)))
save(file="obsdata_s.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)

# Source-site sample size n=600.
dat.homo <- DGP(case="homo", N=c(300, rep(600, 4)))
dat.diffX <- DGP(case="diffX", N=c(300, rep(600, 4)))
dat.diffT <- DGP(case="diffT", N=c(300, rep(600, 4)))
dat.diffC <- DGP(case="diffC", N=c(300, rep(600, 4)))
dat.diffAll <- DGP(case="diffAll", N=c(300, rep(600, 4)))
save(file="obsdata_l.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)

# Source-site sample size n=1000.
dat.homo <- DGP(case="homo", N=c(300, rep(1000, 4)))
dat.diffX <- DGP(case="diffX", N=c(300, rep(1000, 4)))
dat.diffT <- DGP(case="diffT", N=c(300, rep(1000, 4)))
dat.diffC <- DGP(case="diffC", N=c(300, rep(1000, 4)))
dat.diffAll <- DGP(case="diffAll", N=c(300, rep(1000, 4)))
save(file="obsdata_l2.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)


### --- data generation process for the limited-overlap setting (Section E.4)
# This setting uses a more extreme target-site propensity-score model to create poorer treatment overlap while retaining the same outcome and censoring mechanisms.
DGP_limO <- function(M=600, N=rep(1000, 5), case="homo", s=311) {
  
  # Logistic inverse-link function.
  expit <- function(x) 1/(1+exp(-x))
  
  # Construct the site indicator: site 0 is the target and sites 1,...,K-1 are sources.
  N.site <- length(N)
  site <- list()
  for(i in 1:N.site) { site[[i]] <- rep(i-1, N[i]) }
  site <- unlist(site)
  R <- max(site)
  
  # Store all Monte Carlo replicates.
  dat <- list()
  
  # Fix the random seed for reproducibility.
  set.seed(s)
  
  for(i in 1:M) {
    
    # --------------------------------------------------------
    # Step 1: Generate baseline covariates
    # --------------------------------------------------------
    # Start from the target-site covariate distribution.
    r <- 0
    X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
    X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
    X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
    
    # --------------------------------------------------------
    # Step 2: Specify the cross-site distribution-shift scenario
    # --------------------------------------------------------
    
    if(case=="homo") {
      D.C <- D.T <- 0
    }
    
    if(case=="diffX") {
      
      # Induce covariate-distribution shifts in Source Sites 2--4 while keeping Source Site 1 aligned with the target.
      for(r in 2:R) {
        X1[site==r] <- x1 <- 33*rbeta(N[r+1], shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
        X2[site==r] <- 52*rbeta(N[r+1], shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
        X3[site==r] <- (4+2*r)*rbeta(N[r+1], shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
      }
      
      D.C <- D.T <- 0
    }
    
    if(case=="diffT") {
      # Induce site-specific shifts in the event-time model.
      D.T <- site
      D.C <- 0
    }
    
    if(case=="diffC") {
      # Induce site-specific shifts in the censoring model.
      D.T <- 0
      D.C <- site
    }
    
    if(case=="diffAll") {
      
      # Simultaneously induce covariate-distribution shifts in Source Sites 2--4.
      for(r in 2:R) {
        X1[site==r] <- x1 <- 33*rbeta(N[r+1], shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
        X2[site==r] <- 52*rbeta(N[r+1], shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
        X3[site==r] <- (4+2*r)*rbeta(N[r+1], shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
      }
      
      # Also induce site-specific shifts in the event-time and censoring models.
      D.C <- D.T <- site
    }
    
    # --------------------------------------------------------
    # Step 3: Generate treatment assignment under poorer target-site overlap
    # --------------------------------------------------------
    # More extreme propensity-score model for the target site.
    logit.g0s <- -1.05+log(0.3+exp(-120+X1)+exp(-6+X3)+exp(-6+X2/4))
    g0s <- expit(logit.g0s)
    A.tgt <- rbinom(sum(N), size=1, prob=g0s)
    
    # Standard propensity-score model for the source sites.
    logit.g0s <- -1.05+log(1.3+exp(-12+X1/10)+exp(-2+X3/3)+exp(-2+X2/12))
    g0s <- expit(logit.g0s)
    A.src <- rbinom(sum(N), size=1, prob=g0s)
    
    # Use the limited-overlap treatment assignment in the target site and the standard assignment in source sites.
    A <- A.tgt*I(site==0) + A.src*I(site!=0)
    
    # --------------------------------------------------------
    # Step 4: Generate event times
    # --------------------------------------------------------
    
    # Heterogeneous treatment effect with optional site-specific outcome shift controlled by D.T.
    delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2) + D.T*0.02*(X1+X3-25)
    
    # Multiplicative event-time hazard factor.
    ph.t <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) - D.T*0.03*(X3-X1+20) + A*delta.t)
    
    # --------------------------------------------------------
    # Step 5: Generate censoring times
    # --------------------------------------------------------
    
    # Treatment effect in the censoring model with optional site-specific censoring shift controlled by D.C.
    delta.c <- -0.4+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X1+X2+X3-50)
    
    # Multiplicative censoring-time hazard factor.
    ph.c <- exp(-4.87+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X2-25) + A*delta.c)
    
    # Generate event and censoring times from Weibull-type models by inverse-transform sampling.
    u1 <- runif(sum(N))
    u2 <- runif(sum(N))
    rho <- 1.2
    lambda <- 0.6
    event.time <- (-log(u1)/(ph.t * lambda))^(1/rho)
    cens.time <- (-log(u2)/(ph.c * lambda))^(1/rho)
    
    # Apply administrative censoring at day 200.
    cens.time[cens.time > 200] <- 200
    
    # Construct the observed follow-up time and event indicator.
    obs.time <- pmin(event.time, cens.time)
    obs.event <- as.numeric(event.time <= cens.time)
    
    # Store one simulated dataset.
    dat[[i]] <- data.frame(site=site, X1=X1, X2=X2, X3=X3, A=A, Y=obs.time, Delta=obs.event)
  }
  
  return(dat=dat)
}

### --- generate the limited-overlap simulation datasets
# Use n=300 participants in the target site and each of the four source sites.
dat.homo <- DGP_limO(case="homo", N=c(300, rep(300, 4)))
dat.diffX <- DGP_limO(case="diffX", N=c(300, rep(300, 4)))
dat.diffT <- DGP_limO(case="diffT", N=c(300, rep(300, 4)))
dat.diffC <- DGP_limO(case="diffC", N=c(300, rep(300, 4)))
dat.diffAll <- DGP_limO(case="diffAll", N=c(300, rep(300, 4)))
save(file="obsdata_limO.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)

### --- run clustered Cox simulation results
# Run the clustered Cox benchmark implemented in a separate script.
# Approximate runtime on the authors' computing environment is about 10 minutes.
source("ClusterCox.R")
