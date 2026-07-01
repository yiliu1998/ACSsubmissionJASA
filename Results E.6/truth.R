#### true survival values
N <- 10^6
r <- 0
X1 <- 33*rbeta(sum(N), shape1=1.1-0.05*r, shape2=1.1+0.2*r) + 9 + 2*r
X2 <- 52*rbeta(sum(N), shape1=1.5+(X1+0.5*r)/20, shape2=4+2*r) + 7 + 2*r
X3 <- (4+2*r)*rbeta(sum(N), shape1=1.5+abs(X1-50+3*r)/20, shape2=3+0.1*r)

rho <- 1.2
lambda <- 0.6
delta.t <- -0.36-0.1*(X1-25)+0.05*(X2-25)+0.05*(X3-2)
ph.t.0 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2))
ph.t.1 <- exp(-5.02+0.1*(X1-25)-0.1*(X2-25)+0.05*(X3-2) + delta.t)

t.RMST <- seq(1,60,0.1)
t <- 1:60
dt <- diff(c(0,t.RMST))
S0.true <- sapply(t, function(time) mean(exp(-ph.t.0*lambda*time^rho)))
S1.true <- sapply(t, function(time) mean(exp(-ph.t.1*lambda*time^rho)))
S0.true.RMST <- sapply(t.RMST, function(time) mean(exp(-ph.t.0*lambda*time^rho)))
S1.true.RMST <- sapply(t.RMST, function(time) mean(exp(-ph.t.1*lambda*time^rho)))
RMST.1.true <- S1.true.RMST %*% dt
RMST.0.true <- S0.true.RMST %*% dt
RMST.diff.true <- RMST.1.true-RMST.0.true
RD.true <- S1.true-S0.true
SR.true <- S1.true/S0.true

save(file="truth.Rdata", 
     S0.true, S1.true, RD.true, SR.true,
     RMST.1.true, RMST.0.true, RMST.diff.true)
