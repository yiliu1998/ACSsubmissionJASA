### ----------------------------------------------------------------------------
### Reproducibility: Data-generating process and true values for Section E.5
###
### This script:
###   (1) computes Monte Carlo approximations to the true target-site survival
###       contrasts, including survival difference, survival ratio, and RMST;
###   (2) generates the simulation datasets under the five distribution-shift
###       scenarios considered in Section E.5; and
###   (3) saves the truth values and simulated datasets for downstream analyses.
###
### Output files:
###   truth.Rdata    : true survival curves and contrast values
###   obsdata_s.Rdata: simulated datasets with n = 300 per site
### ----------------------------------------------------------------------------


### --- true values of survival contrasts

# Fix the seed for reproducibility of the large-sample Monte Carlo approximation.
set.seed(20260702)

# Large Monte Carlo sample used to approximate the true target-site quantities.
N <- 10^7
r <- 0

# Generate baseline covariates under the target-site distribution.
X1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
X2 <- 52*rbeta(sum(N), shape1=1.5+(X1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r
X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(X1-50+3*r)/20, shape2=3+0.1*r)

# Weibull-type event-time model parameters.
rho <- 1.2
lambda <- 0.6

# Define the heterogeneous treatment effect and treatment-specific event-time hazard factors.
delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2)
ph.t.0 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2))
ph.t.1 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) + delta.t)

# Use the day 1--60 grid for survival contrasts and RMST up to day 60.
t.RMST <- seq(1,60,1)
t <- 1:60
dt <- diff(c(0,t.RMST))

# Compute marginal survival probabilities under control and treatment.
S0.true <- sapply(t, function(time) mean(exp(-ph.t.0*lambda*time^rho)))
S1.true <- sapply(t, function(time) mean(exp(-ph.t.1*lambda*time^rho)))

# Compute survival probabilities on the RMST integration grid.
S0.true.RMST <- sapply(t.RMST, function(time) mean(exp(-ph.t.0*lambda*time^rho)))
S1.true.RMST <- sapply(t.RMST, function(time) mean(exp(-ph.t.1*lambda*time^rho)))

# Compute true RMST values under treatment and control and their difference.
RMST.1.true <- S1.true.RMST %*% dt
RMST.0.true <- S0.true.RMST %*% dt
RMST.diff.true <- RMST.1.true-RMST.0.true

# Compute the true survival difference and survival ratio over time.
RD.true <- S1.true-S0.true
SR.true <- S1.true/S0.true

# Display the true survival contrasts at days 30 and 60 and RMST values up to day 60.
RD.true[c(30,60)]
SR.true[c(30,60)]
RMST.1.true
RMST.0.true
RMST.diff.true

# Save truth values for downstream bias, RMSE, and coverage calculations.
save(file="truth.Rdata", 
     S0.true, S1.true, RD.true, SR.true,
     RMST.1.true, RMST.0.true, RMST.diff.true)


### --- data-generating process for the Section E.5 simulation study
# ==============================================================================
# Inputs:
#   M    : number of Monte Carlo datasets generated
#   N    : vector of site-specific sample sizes; site 0 is the target site
#   case : distribution-shift scenario
#   s    : random seed
#
# Output:
#   a list containing M independently generated multi-site datasets
# ==============================================================================
DGP <- function(M=600, N=rep(1000, 5), case="homo", s=311) {
  
  # Logistic inverse-link function.
  expit <- function(x) 1/(1+exp(-x))
  
  # Construct the site indicator: site 0 is the target and sites 1,...,K-1 are sources.
  N.site <- length(N)
  site <- list()
  for(i in 1:N.site) { site[[i]] <- rep(i-1, N[i]) }
  site <- unlist(site)
  R <- max(site)
  
  # Store the independently generated Monte Carlo datasets.
  dat <- list()
  
  # Fix the seed for reproducibility.
  set.seed(s)
  
  for(i in 1:M) {
    
    # --------------------------------------------------------
    # Step 1: Generate baseline covariates
    # --------------------------------------------------------
    
    # Begin with the target-site covariate distribution for all observations.
    r <- 0
    X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
    X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
    X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
    
    
    # --------------------------------------------------------
    # Step 2: Specify the cross-site distribution-shift scenario
    #
    # "homo"    : homogeneous DGP across sites
    # "diffX"   : covariate shift
    # "diffT"   : event-time/outcome shift
    # "diffC"   : censoring shift
    # "diffAll" : covariate, event-time, and censoring shifts
    # --------------------------------------------------------
    
    if(case=="homo") {
      D.C <- D.T <- 0
    }
    
    if(case=="diffX") {
      
      # Induce covariate-distribution shifts in Source Sites 2--4 while
      # keeping Source Site 1 aligned with the target site.
      for(r in 2:R) {
        X1[site==r] <- x1 <- 33*rbeta(N[r+1], shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
        X2[site==r] <- 52*rbeta(N[r+1], shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
        X3[site==r] <- (4+2*r)*rbeta(N[r+1], shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
      }
      
      # Keep the event-time and censoring models otherwise unchanged.
      D.C <- D.T <- 0
    }
    
    if(case=="diffT") {
      # Induce progressively stronger source-specific shifts in the event-time model.
      D.T <- site
      D.C <- 0
    }
    
    if(case=="diffC") {
      # Induce progressively stronger source-specific shifts in the censoring model.
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
    # Step 3: Generate treatment assignment
    # --------------------------------------------------------
    
    # Define source-site and limited-overlap propensity-score linear predictors.
    eta.src <- -1.05+log(1.3+exp(-12+X1/10)+exp(-2+X3/3)+exp(-2+X2/12))
    eta.lim <- -1.05+log(0.3+exp(-120+X1)+exp(-6+X3)+exp(-6+X2/4))
    
    # The target-site propensity score is a mixture of the two specifications.
    eta.tgt <- 0.7*eta.src+0.3*eta.lim
    
    # Generate treatment according to the target- or source-specific propensity score.
    g0s <- expit(ifelse(site==0, eta.tgt, eta.src))
    A <- rbinom(sum(N), size=1, prob=g0s)
    
    
    # --------------------------------------------------------
    # Step 4: Generate event times
    # --------------------------------------------------------
    
    # Heterogeneous treatment effect with optional source-specific event-time shift.
    delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2) + D.T*0.02*(X1+X3-25)
    
    # Multiplicative event-time hazard factor.
    ph.t <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) - D.T*0.03*(X3-X1+20) + A*delta.t)
    
    
    # --------------------------------------------------------
    # Step 5: Generate censoring times
    # --------------------------------------------------------
    
    # Treatment effect in the censoring model with optional source-specific censoring shift.
    delta.c <- -0.4+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X1+X2+X3-50)
    
    # Multiplicative censoring-time hazard factor.
    ph.c <- exp(-4.87+0.01*(X1-25)-0.02*(X2-25)+0.01*(X3-2) - D.C*0.025*(X2-25) + A*delta.c)
    
    # Draw event and censoring times using inverse-transform sampling.
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
    
    # Store one complete simulated multi-site dataset.
    dat[[i]] <- data.frame(site=site, X1=X1, X2=X2, X3=X3, A=A, Y=obs.time, Delta=obs.event)
  }
  
  return(dat=dat)
}


### --- generate simulation datasets for Section E.5

# Generate 600 candidate Monte Carlo datasets with n=300 observations per site.
# Downstream code retains the requested number of successful simulation replicates.
dat.homo <- DGP(case="homo", N=c(300, rep(300, 4)))
dat.diffX <- DGP(case="diffX", N=c(300, rep(300, 4)))
dat.diffT <- DGP(case="diffT", N=c(300, rep(300, 4)))
dat.diffC <- DGP(case="diffC", N=c(300, rep(300, 4)))
dat.diffAll <- DGP(case="diffAll", N=c(300, rep(300, 4)))

# Save all five simulation scenarios for the Section E.5 analysis.
save(file="obsdata_s.Rdata", dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
