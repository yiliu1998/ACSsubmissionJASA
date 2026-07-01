library(ggplot2)
library(gridExtra)

#### Label shift
N <- 10^4
r <- 0
X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)

rho <- 1.2
lambda <- 0.6
t <- 0:90
p.survs <- list()
delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2) 
ph.t.0 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2))
ph.t.1 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) + delta.t)

# Baseline calculations
S0_baseline <- sapply(t, function(time) mean(exp(-ph.t.0 * lambda * time^rho)))
S1_baseline <- sapply(t, function(time) mean(exp(-ph.t.1 * lambda * time^rho)))

baseline_data <- data.frame(
  Time = rep(t, 2),
  Survival = c(S0_baseline, S1_baseline),
  Group = rep(c("Control", "Treated"), each = length(t))
)

# Create baseline plot
p.survs[[1]] <- ggplot(baseline_data, aes(x=Time, y=Survival, color=Group)) +
  geom_line(size=1) + ylim(0.3,1) + 
  labs(title="Target Site", x="Time", y="Survival Probability") +
  theme(legend.position="none")

for(D.T in 1:3) {
  delta.t <- -0.36 - 0.1*(X1-25) + 0.05*(X2-25) + 0.05*(X3-2) + D.T*0.02*(X1+X3-25)
  ph.t.0 <- exp(-5.02 + 0.1*(X1-25) - 0.1*(X2-25) + 0.05*(X3-2) - D.T*0.03*(X3 - X1 + 20))
  ph.t.1 <- exp(-5.02 + 0.1*(X1-25) - 0.1*(X2-25) + 0.05*(X3-2) - D.T*0.03*(X3 - X1 + 20) + delta.t)
  
  S0 <- sapply(t, function(time) mean(exp(-ph.t.0 * lambda * time^rho)))
  S1 <- sapply(t, function(time) mean(exp(-ph.t.1 * lambda * time^rho)))
  
  survival_data <- data.frame(
    Time = rep(t, 2),
    Survival = c(S0, S1),
    Group = rep(c("Control", "Treated"), each = length(t))
  )
  
  p.survs[[D.T+1]] <- ggplot() +
    # Source site curves
    geom_line(data=survival_data, aes(x=Time, y=Survival, color=Group), size=1) +
    # Add baseline curves in a neutral color or linetype
    geom_line(data=baseline_data, aes(x=Time, y=Survival, group=Group, color=Group), 
              linetype="dashed", alpha=0.5, size=0.8) +
    ylim(0.3,1) +
    labs(title=paste0("Source Site ", D.T), x="Time", y="") +
    theme_minimal() + theme(legend.position="none")
}

D.T <- 4
delta.t <- -0.36 - 0.1*(X1-25) + 0.05*(X2-25) + 0.05*(X3-2) + D.T*0.02*(X1+X3-25)
ph.t.0 <- exp(-5.02 + 0.1*(X1-25) - 0.1*(X2-25) + 0.05*(X3-2) - D.T*0.03*(X3 - X1 + 20))
ph.t.1 <- exp(-5.02 + 0.1*(X1-25) - 0.1*(X2-25) + 0.05*(X3-2) - D.T*0.03*(X3 - X1 + 20) + delta.t)

S0 <- sapply(t, function(time) mean(exp(-ph.t.0 * lambda * time^rho)))
S1 <- sapply(t, function(time) mean(exp(-ph.t.1 * lambda * time^rho)))

survival_data <- data.frame(
  Time = rep(t, 2),
  Survival = c(S0, S1),
  Group = rep(c("Control", "Treated"), each = length(t))
)

p.survs[[D.T+1]] <- ggplot() +
  geom_line(data=survival_data, aes(x=Time, y=Survival, color=Group), size=1) +
  geom_line(data=baseline_data, aes(x=Time, y=Survival, group=Group, color=Group), 
            linetype="dashed", alpha=0.5, size=0.8) +
  ylim(0.3,1) + 
  labs(title=paste0("Source Site ", D.T),
       x="Time",
       y="",
       color="Group") +
  theme_minimal()

pdf("surv_curves_ls.pdf", height=2.5, width=15)
grid.arrange(p.survs[[1]], p.survs[[2]], p.survs[[3]], p.survs[[4]], p.survs[[5]], 
             nrow=1, widths=c(1,1,1,1,1.4))
dev.off()


#### Covariate shift
N <- 10^4
r <- 0
X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)

rho <- 1.2
lambda <- 0.6
t <- 0:90
p.survs <- list()
delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2) 
ph.t.0 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2))
ph.t.1 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) + delta.t)

S0 <- sapply(t, function(time) mean(exp(-ph.t.0*lambda*time^rho)))
S1 <- sapply(t, function(time) mean(exp(-ph.t.1*lambda*time^rho)))

survival_data <- data.frame(
  Time=rep(t, 2),
  Survival=c(S0, S1),
  Group=rep(c("Control", "Treated"), each=length(t))
)

ggplot(survival_data, aes(x=Time, y=Survival, color=Group)) +
  geom_line(size=1) + ylim(0.3,1) + 
  labs(title=paste0("Target Site"),
       x="Time",
       y="Survival Probability") +
  theme(legend.position="none") -> p.survs[[1]]

for(r in 1:3) {
  X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
  X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
  X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
  delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2)
  ph.t.0 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2))
  ph.t.1 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) + delta.t)
  
  S0 <- sapply(t, function(time) mean(exp(-ph.t.0*lambda*time^rho)))
  S1 <- sapply(t, function(time) mean(exp(-ph.t.1*lambda*time^rho)))
  
  survival_data <- data.frame(
    Time=rep(t, 2),
    Survival=c(S0, S1),
    Group=rep(c("Control", "Treated"), each=length(t))
  )
  
  p.survs[[r+1]] <- ggplot() +
    # Source site curves
    geom_line(data=survival_data, aes(x=Time, y=Survival, color=Group), size=1) +
    # Add baseline curves in a neutral color or linetype
    geom_line(data=baseline_data, aes(x=Time, y=Survival, group=Group, color=Group), 
              linetype="dashed", alpha=0.5, size=0.8) +
    ylim(0.3,1) +
    labs(title=paste0("Source Site ", r), x="Time", y="") +
    theme_minimal() + theme(legend.position="none") 
}

r <- 4
X1 <- x1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
X2 <- 52*rbeta(sum(N), shape1=1.5+(x1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r 
X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(x1-50+3*r)/20, shape2=3+0.1*r)
delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2)
ph.t.0 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2))
ph.t.1 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) + delta.t)

S0 <- sapply(t, function(time) mean(exp(-ph.t.0*lambda*time^rho)))
S1 <- sapply(t, function(time) mean(exp(-ph.t.1*lambda*time^rho)))

survival_data <- data.frame(
  Time=rep(t, 2),
  Survival=c(S0, S1),
  Group=rep(c("Control", "Treated"), each=length(t))
)

p.survs[[r+1]] <- ggplot() +
  # Source site curves
  geom_line(data=survival_data, aes(x=Time, y=Survival, color=Group), size=1) +
  # Add baseline curves in a neutral color or linetype
  geom_line(data=baseline_data, aes(x=Time, y=Survival, group=Group, color=Group), 
            linetype="dashed", alpha=0.5, size=0.8) +
  ylim(0.3,1) +
  labs(title=paste0("Source Site ", 4), x="Time", y="") +
  theme_minimal() 

pdf("surv_curves_covs.pdf", height=2.5, width=15)
grid.arrange(p.survs[[1]], p.survs[[2]], p.survs[[3]], p.survs[[4]], 
             p.survs[[5]], nrow=1, widths=c(1,1,1,1,1.4))
dev.off()


