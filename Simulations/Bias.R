library(tidyverse)
library(RColorBrewer)

M <- 500
all.bias.list <- function(datlist, coxlist, M=500) {
  TGT.1 <- IVW.1 <- POOL.1 <- CCOD.1 <- FED.1 <- CLCOX.1 <- c()
  TGT.0 <- IVW.0 <- POOL.0 <- CCOD.0 <- FED.0 <- CLCOX.0 <- c()
  
  for(j in 1:M) {
    TGT.1 <- c(TGT.1, datlist[[j]]$df.TGT[, "surv1"] - S1.true)
    IVW.1 <- c(IVW.1, datlist[[j]]$df.IVW[, "surv1"] - S1.true)
    POOL.1 <- c(POOL.1, datlist[[j]]$df.POOL[, "surv1"] - S1.true)
    CCOD.1 <- c(CCOD.1, datlist[[j]]$df.CCOD[, "surv1"] - S1.true)
    FED.1 <- c(FED.1, datlist[[j]]$df.FED[, "surv1"] - S1.true)
    CLCOX.1 <- c(CLCOX.1, coxlist[[j]]$df.CLCOX[, "surv1"] - S1.true)
    
    TGT.0 <- c(TGT.0, datlist[[j]]$df.TGT[, "surv0"] - S0.true)
    IVW.0 <- c(IVW.0, datlist[[j]]$df.IVW[, "surv0"] - S0.true)
    POOL.0 <- c(POOL.0, datlist[[j]]$df.POOL[, "surv0"] - S0.true)
    CCOD.0 <- c(CCOD.0, datlist[[j]]$df.CCOD[, "surv0"] - S0.true)
    FED.0 <- c(FED.0, datlist[[j]]$df.FED[, "surv0"] - S0.true)
    CLCOX.0 <- c(CLCOX.0, coxlist[[j]]$df.CLCOX[, "surv0"] - S0.true)
  }
  
  return(list(TGT.1, IVW.1, POOL.1, CCOD.1, FED.1, CLCOX.1,
              TGT.0, IVW.0, POOL.0, CCOD.0, FED.0, CLCOX.0))
}

load("Res_homo_s.Rdata")
load("Res_diffX_s.Rdata")
load("Res_diffT_s.Rdata")
load("Res_diffC_s.Rdata")
load("Res_diffAll_s.Rdata")
load("Res_CLCOX_s.Rdata")
load("truth.Rdata")

homo <- all.bias.list(result.homo, result.CLCOX.homo)
diffX <- all.bias.list(result.diffX, result.CLCOX.diffX)
diffT <- all.bias.list(result.diffT, result.CLCOX.diffT)
diffC <- all.bias.list(result.diffC, result.CLCOX.diffC)
diffAll <- all.bias.list(result.diffAll, result.CLCOX.diffAll)

all.bias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.bias)
t <- c(30, 60, 90)

time <- rep(t, n / length(t))

Case <- factor(
  c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
    rep("Censoring Shift", n/5), rep("All Shift", n/5)),
  levels = c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift")
)

Method <- factor(
  rep(c(rep("TGT", length(t)*M),
        rep("IVW", length(t)*M),
        rep("POOL", length(t)*M),
        rep("CCOD", length(t)*M),
        rep("FED", length(t)*M),
        rep("CLCOX", length(t)*M)), 10),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df <- data.frame(Bias = all.bias, Case = Case, Method = Method, time = time, Treatment = Treatment)

colors <- c(brewer.pal(11, "Paired")[c(1, 3, 7, 8, 10)], "plum3")

ggplot(df, aes(x = factor(time), y = Bias, fill = Method)) + 
  geom_boxplot(lwd = 0.3, outlier.size = 0.6, color = "gray50") + 
  scale_fill_manual(values = colors) + 
  geom_hline(yintercept = 0, col = "hotpink2", lty = 2, linewidth = 0.6) + 
  facet_grid(Treatment ~ Case) + 
  labs(x = "Time (day)", y = "Bias") + 
  theme_bw() -> p.bias

pdf(file = "bias_s.pdf", width = 9, height = 4)
p.bias
dev.off()


load("Res_homo_l.Rdata")
load("Res_diffX_l.Rdata")
load("Res_diffT_l.Rdata")
load("Res_diffC_l.Rdata")
load("Res_diffAll_l.Rdata")
load("Res_CLCOX_l.Rdata")

homo <- all.bias.list(result.homo, result.CLCOX.homo)
diffX <- all.bias.list(result.diffX, result.CLCOX.diffX)
diffT <- all.bias.list(result.diffT, result.CLCOX.diffT)
diffC <- all.bias.list(result.diffC, result.CLCOX.diffC)
diffAll <- all.bias.list(result.diffAll, result.CLCOX.diffAll)

all.bias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.bias)

time <- rep(t, n / length(t))

Case <- factor(
  c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
    rep("Censoring Shift", n/5), rep("All Shift", n/5)),
  levels = c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift")
)

Method <- factor(
  rep(c(rep("TGT", length(t)*M),
        rep("IVW", length(t)*M),
        rep("POOL", length(t)*M),
        rep("CCOD", length(t)*M),
        rep("FED", length(t)*M),
        rep("CLCOX", length(t)*M)), 10),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df <- data.frame(Bias = all.bias, Case = Case, Method = Method, time = time, Treatment = Treatment)

ggplot(df, aes(x = factor(time), y = Bias, fill = Method)) + 
  geom_boxplot(lwd = 0.3, outlier.size = 0.6, color = "gray50") + 
  scale_fill_manual(values = colors) + 
  geom_hline(yintercept = 0, col = "hotpink2", lty = 2, linewidth = 0.6) + 
  facet_grid(Treatment ~ Case) + 
  labs(x = "Time (day)", y = "Bias") + 
  theme_bw() -> p.bias

pdf(file = "bias_l.pdf", width = 9, height = 4)
p.bias
dev.off()

load("Res_homo_l2.Rdata")
load("Res_diffX_l2.Rdata")
load("Res_diffT_l2.Rdata")
load("Res_diffC_l2.Rdata")
load("Res_diffAll_l2.Rdata")
load("Res_CLCOX_l2.Rdata")

homo <- all.bias.list(result.homo, result.CLCOX.homo)
diffX <- all.bias.list(result.diffX, result.CLCOX.diffX)
diffT <- all.bias.list(result.diffT, result.CLCOX.diffT)
diffC <- all.bias.list(result.diffC, result.CLCOX.diffC)
diffAll <- all.bias.list(result.diffAll, result.CLCOX.diffAll)

all.bias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.bias)

time <- rep(t, n / length(t))

Case <- factor(
  c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
    rep("Censoring Shift", n/5), rep("All Shift", n/5)),
  levels = c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift")
)

Method <- factor(
  rep(c(rep("TGT", length(t)*M),
        rep("IVW", length(t)*M),
        rep("POOL", length(t)*M),
        rep("CCOD", length(t)*M),
        rep("FED", length(t)*M),
        rep("CLCOX", length(t)*M)), 10),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df <- data.frame(Bias = all.bias, Case = Case, Method = Method, time = time, Treatment = Treatment)

ggplot(df, aes(x = factor(time), y = Bias, fill = Method)) + 
  geom_boxplot(lwd = 0.3, outlier.size = 0.6, color = "gray50") + 
  scale_fill_manual(values = colors) + 
  geom_hline(yintercept = 0, col = "hotpink2", lty = 2, linewidth = 0.6) +
  facet_grid(Treatment ~ Case) + 
  labs(x = "Time (day)", y = "Bias") + 
  theme_bw() -> p.bias

pdf(file = "bias_l2.pdf", width = 9, height = 4)
p.bias
dev.off()


load("Res_homo_limO.Rdata")
result.homo <- results
load("Res_diffX_limO.Rdata")
result.diffX <- results
load("Res_diffT_limO.Rdata")
result.diffT <- results
load("Res_diffC_limO.Rdata")
result.diffC <- results
load("Res_diffAll_limO.Rdata")
result.diffAll <- results
load("Res_CLCOX_limO.Rdata")

homo <- all.bias.list(result.homo, result.CLCOX.homo)
diffX <- all.bias.list(result.diffX, result.CLCOX.diffX)
diffT <- all.bias.list(result.diffT, result.CLCOX.diffT)
diffC <- all.bias.list(result.diffC, result.CLCOX.diffC)
diffAll <- all.bias.list(result.diffAll, result.CLCOX.diffAll)

all.bias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.bias)

time <- rep(t, n / length(t))

Case <- factor(
  c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
    rep("Censoring Shift", n/5), rep("All Shift", n/5)),
  levels = c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift")
)

Method <- factor(
  rep(c(rep("TGT", length(t)*M),
        rep("IVW", length(t)*M),
        rep("POOL", length(t)*M),
        rep("CCOD", length(t)*M),
        rep("FED", length(t)*M),
        rep("CLCOX", length(t)*M)), 10),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df <- data.frame(Bias = all.bias, Case = Case, Method = Method, time = time, Treatment = Treatment)

ggplot(df, aes(x = factor(time), y = Bias, fill = Method)) + 
  geom_boxplot(lwd = 0.3, outlier.size = 0.6, color = "gray50") + 
  scale_fill_manual(values = colors) + 
  geom_hline(yintercept = 0, col = "hotpink2", lty = 2, linewidth = 0.6) +
  facet_grid(Treatment ~ Case) + 
  labs(x = "Time (day)", y = "Bias") + 
  theme_bw() -> p.bias

pdf(file = "bias_limO.pdf", width = 9, height = 4)
p.bias
dev.off()