estimate_omega_np = function(x, x_target, x.pred, method="logistic") {
  x.all = rbind(x, x_target)
  z.all = c(rep(1, nrow(x)), rep(0, nrow(x_target)))
  
  if(method=="logistic") {
    colnames(x.all) <- colnames(x.pred) <- paste0("X", 1:ncol(x.all))
    fit <- glm(z.all~.-z.all, data=data.frame(z.all, x.all), family=binomial(link="logit"))
    src.predict = predict(fit, newdata=data.frame(x.pred), type="response")
  }
  if(method=="glmnet") {
    fit <- cv.glmnet(x=as.matrix(x.all), y=z.all, nfolds=5, family='binomial')
    src.predict = predict(fit, x.pred, type="response", lambda=fit$lambda.min)
  }
  omega = (1-src.predict)/src.predict * nrow(x)/nrow(x_target)
  omega = pmax(pmin(omega, 20), 0.05)
  return(omega)
}

get.survival <- function(Y, Delta, A, R=0,
                         fit.times, 
                         S.hats, G.hats, g.hats, omega.hats=1) {
  
  fit.times <- fit.times[fit.times > 0]
  n <- length(Y)
  ord <- order(fit.times)
  fit.times <- fit.times[ord]
  S.hats <- S.hats[, ord]
  G.hats <- G.hats[, ord]
  
  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S.hats[i,])* 1/ G.hats[i,-ncol(G.hats)]
    if(any(fit.times[-1] > Y[i])) vals[fit.times[-1] > Y[i]] <- 0
    c(0, cumsum(vals))
  }))
  
  S.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
  G.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))
  IF.vals <- matrix(NA, nrow=n, ncol=length(fit.times))
  surv <- AUG.means <- rep(NA, length(fit.times))
  
  for(t0 in fit.times) {
    k <- min(which(fit.times>=t0))
    S.hats.t0 <- S.hats[,k]
    inner.func.1 <- ifelse(Y<=t0 & Delta==1, 1/(S.hats.Y*G.hats.Y), 0 )
    inner.func.2 <- int.vals[,k]
    k1 <- which(fit.times==t0)
    augment <- omega.hats*S.hats.t0*as.numeric(A==1)*(inner.func.1 - inner.func.2)/g.hats
    if.func <- S.hats.t0 - augment
    surv[k1] <- mean(if.func)
    # IF.vals[,k1] <- (if.func - mean(if.func))*I(R==0) - augment*I(R!=0)
    IF.vals[,k1] <- if.func*I(R==0) - augment*I(R!=0)
    AUG.means[k1] <- mean(augment)
  }
  surv = pmin(1, pmax(0, surv))
  # surv.sd <- sqrt(colMeans(IF.vals^2, na.rm=T)/n)
  surv.sd <- sqrt(apply(IF.vals, 2, var, na.rm=T)/n)
  return(list(IF.vals=IF.vals, AUG.means=AUG.means, surv=surv, surv.sd=surv.sd))
}

get.survival.CCOD <- function(Y, Delta, A, R,
                              fit.times, 
                              S.hats, G.hats, g.hats, eta0.hats) {
  
  fit.times <- fit.times[fit.times > 0]
  n <- length(Y)
  ord <- order(fit.times)
  fit.times <- fit.times[ord]
  S.hats <- S.hats[, ord]
  G.hats <- G.hats[, ord]
  
  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S.hats[i,])* 1/ G.hats[i,-ncol(G.hats)]
    if(any(fit.times[-1] > Y[i])) vals[fit.times[-1] > Y[i]] <- 0
    c(0, cumsum(vals))
  }))
  
  S.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
  G.hats.Y <- sapply(1:n, function(i) stepfun(fit.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))
  IF.vals <- matrix(NA, nrow=n, ncol=length(fit.times))
  surv <- surv.sd <- rep(NA, length(fit.times))
  
  for(t0 in fit.times) {
    k <- min(which(fit.times>=t0))
    S.hats.t0 <- S.hats[,k]
    inner.func.1 <- ifelse(Y<=t0 & Delta==1, 1/(S.hats.Y*G.hats.Y), 0 )
    inner.func.2 <- int.vals[,k]
    k1 <- which(fit.times==t0)
    augment <- S.hats.t0*as.numeric(A==1)*(inner.func.1 - inner.func.2)/g.hats
    if.func <- (S.hats.t0*I(R==1) - eta0.hats*augment) / mean(R==1)
    surv[k1] <- mean(if.func)
    IF.vals[,k1] <- if.func # - mean(if.func)
    surv.sd[k1] <- sqrt(var(if.func[R==1])*mean(R==1) + var(if.func[R!=1])*mean(R!=1)) / sqrt(n)
  }
  surv = pmin(1, pmax(0, surv))
  return(list(IF.vals=IF.vals, surv=surv, surv.sd=surv.sd))
}