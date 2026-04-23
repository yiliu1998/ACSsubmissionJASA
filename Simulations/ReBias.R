M <- 500
RelBias.results <- function(datlist, coxlist, M=500) {
  TGT.1 <- IVW.1 <- POOL.1 <- CCOD.1 <- FED.1 <- CLCOX.1 <- c()
  TGT.0 <- IVW.0 <- POOL.0 <- CCOD.0 <- FED.0 <- CLCOX.0 <- c()
  
  for(j in 1:M) {
    time.use <- datlist[[j]]$df.TGT$time
    df.cox <- coxlist[[j]]$df.CLCOX
    df.cox <- df.cox[match(time.use, df.cox$time), ]
    
    TGT.1 <- c(TGT.1, ((datlist[[j]]$df.TGT[, "surv1"] - S1.true) / S1.true))
    IVW.1 <- c(IVW.1, ((datlist[[j]]$df.IVW[, "surv1"] - S1.true) / S1.true))
    POOL.1 <- c(POOL.1, ((datlist[[j]]$df.POOL[, "surv1"] - S1.true) / S1.true))
    CCOD.1 <- c(CCOD.1, ((datlist[[j]]$df.CCOD[, "surv1"] - S1.true) / S1.true))
    FED.1 <- c(FED.1, ((datlist[[j]]$df.FED[, "surv1"] - S1.true) / S1.true))
    CLCOX.1 <- c(CLCOX.1, ((df.cox[, "surv1"] - S1.true) / S1.true))
    
    TGT.0 <- c(TGT.0, ((datlist[[j]]$df.TGT[, "surv0"] - S0.true) / S0.true))
    IVW.0 <- c(IVW.0, ((datlist[[j]]$df.IVW[, "surv0"] - S0.true) / S0.true))
    POOL.0 <- c(POOL.0, ((datlist[[j]]$df.POOL[, "surv0"] - S0.true) / S0.true))
    CCOD.0 <- c(CCOD.0, ((datlist[[j]]$df.CCOD[, "surv0"] - S0.true) / S0.true))
    FED.0 <- c(FED.0, ((datlist[[j]]$df.FED[, "surv0"] - S0.true) / S0.true))
    CLCOX.0 <- c(CLCOX.0, ((df.cox[, "surv0"] - S0.true) / S0.true))
  }
  
  return(list(absrelbias.TGT.1 = TGT.1,
              absrelbias.IVW.1 = IVW.1,
              absrelbias.POOL.1 = POOL.1,
              absrelbias.CCOD.1 = CCOD.1,
              absrelbias.FED.1 = FED.1,
              absrelbias.CLCOX.1 = CLCOX.1,
              absrelbias.TGT.0 = TGT.0,
              absrelbias.IVW.0 = IVW.0,
              absrelbias.POOL.0 = POOL.0,
              absrelbias.CCOD.0 = CCOD.0,
              absrelbias.FED.0 = FED.0,
              absrelbias.CLCOX.0 = CLCOX.0))
}

library(ggplot2)
library(dplyr)
library(reshape2)
load("truth.Rdata")

load("Res_homo_s.Rdata")
load("Res_diffX_s.Rdata")
load("Res_diffT_s.Rdata")
load("Res_diffC_s.Rdata")
load("Res_diffAll_s.Rdata")
load("Res_CLCOX_s.Rdata")

homo <- RelBias.results(result.homo, result.CLCOX.homo)
diffX <- RelBias.results(result.diffX, result.CLCOX.diffX)
diffT <- RelBias.results(result.diffT, result.CLCOX.diffT)
diffC <- RelBias.results(result.diffC, result.CLCOX.diffC)
diffAll <- RelBias.results(result.diffAll, result.CLCOX.diffAll)

all.relbias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.relbias)

t <- c(30, 60, 90)
time <- rep(rep(t, times = M), times = 5 * 2 * 6)

Case <- factor(
  c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
    rep("Censoring Shift", n/5), rep("All Shift", n/5)),
  levels = c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift")
)

Method <- factor(
  rep(c(rep("TGT", length(t)*M), rep("IVW", length(t)*M), rep("POOL", length(t)*M), 
        rep("CCOD", length(t)*M), rep("FED", length(t)*M), rep("CLCOX", length(t)*M)), 10),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df <- data.frame(RelBias = all.relbias, Case = Case, Method = Method, time = time, Treatment = Treatment)

df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>%
  summarize(RelBias = abs(mean(RelBias * 100)), .groups = 'drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RelBias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RelBias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "white", "green4"),
    values=scales::rescale(c(0, 0.05, 1)),
    name="ARBias%"
  ) +
  facet_grid(Treatment~Case) +
  labs(x="Time (day)", y="") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.relbias

pdf(file = "relbias_s.pdf", width = 9.5, height = 4.5)
p.relbias
dev.off()


load("Res_homo_l.Rdata")
load("Res_diffX_l.Rdata")
load("Res_diffT_l.Rdata")
load("Res_diffC_l.Rdata")
load("Res_diffAll_l.Rdata")
load("Res_CLCOX_l.Rdata")

homo <- RelBias.results(result.homo, result.CLCOX.homo)
diffX <- RelBias.results(result.diffX, result.CLCOX.diffX)
diffT <- RelBias.results(result.diffT, result.CLCOX.diffT)
diffC <- RelBias.results(result.diffC, result.CLCOX.diffC)
diffAll <- RelBias.results(result.diffAll, result.CLCOX.diffAll)

all.relbias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.relbias)

time <- rep(rep(t, times = M), times = 5 * 2 * 6)

Case <- factor(
  c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
    rep("Censoring Shift", n/5), rep("All Shift", n/5)),
  levels = c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift")
)

Method <- factor(
  rep(c(rep("TGT", length(t)*M), rep("IVW", length(t)*M), rep("POOL", length(t)*M), 
        rep("CCOD", length(t)*M), rep("FED", length(t)*M), rep("CLCOX", length(t)*M)), 10),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df <- data.frame(RelBias = all.relbias, Case = Case, Method = Method, time = time, Treatment = Treatment)

df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>%
  summarize(RelBias = abs(mean(RelBias * 100)), .groups = 'drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RelBias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RelBias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "white", "green4"),
    values=scales::rescale(c(0, 0.05, 1)),
    name="ARBias%"
  ) +
  facet_grid(Treatment~Case) +
  labs(x="Time (day)", y="") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.relbias

pdf(file = "relbias_l.pdf", width = 9.5, height = 4.5)
p.relbias
dev.off()


load("Res_homo_l2.Rdata")
load("Res_diffX_l2.Rdata")
load("Res_diffT_l2.Rdata")
load("Res_diffC_l2.Rdata")
load("Res_diffAll_l2.Rdata")
load("Res_CLCOX_l2.Rdata")

homo <- RelBias.results(result.homo, result.CLCOX.homo)
diffX <- RelBias.results(result.diffX, result.CLCOX.diffX)
diffT <- RelBias.results(result.diffT, result.CLCOX.diffT)
diffC <- RelBias.results(result.diffC, result.CLCOX.diffC)
diffAll <- RelBias.results(result.diffAll, result.CLCOX.diffAll)

all.relbias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.relbias)

time <- rep(rep(t, times = M), times = 5 * 2 * 6)

Case <- factor(
  c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
    rep("Censoring Shift", n/5), rep("All Shift", n/5)),
  levels = c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift")
)

Method <- factor(
  rep(c(rep("TGT", length(t)*M), rep("IVW", length(t)*M), rep("POOL", length(t)*M), 
        rep("CCOD", length(t)*M), rep("FED", length(t)*M), rep("CLCOX", length(t)*M)), 10),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df <- data.frame(RelBias = all.relbias, Case = Case, Method = Method, time = time, Treatment = Treatment)

df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>%
  summarize(RelBias = abs(mean(RelBias * 100)), .groups = 'drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RelBias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RelBias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "white", "green4"),
    values=scales::rescale(c(0, 0.05, 1)),
    name="ARBias%"
  ) +
  facet_grid(Treatment~Case) +
  labs(x="Time (day)", y="") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.relbias

pdf(file = "relbias_l2.pdf", width = 9.5, height = 4.5)
p.relbias
dev.off()


### main results in the paper
df.main <- df %>% filter(time == 30)

df_heatmap <- df.main %>%
  group_by(Case, Method, Treatment) %>%
  summarize(RelBias = abs(mean(RelBias * 100)), .groups = 'drop')

ggplot(df_heatmap, aes(x=factor(Case), y=Method, fill=RelBias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RelBias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "white", "green4"),
    values=scales::rescale(c(0, 0.05, 1)),
    name="ARBias%"
  ) +
  facet_grid(.~Treatment) +
  labs(x="", y="") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=35, hjust=1)) -> p.relbias.main

pdf(file = "relbias_main.pdf", width = 5.8, height = 2.8)
p.relbias.main
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

homo <- RelBias.results(result.homo, result.CLCOX.homo)
diffX <- RelBias.results(result.diffX, result.CLCOX.diffX)
diffT <- RelBias.results(result.diffT, result.CLCOX.diffT)
diffC <- RelBias.results(result.diffC, result.CLCOX.diffC)
diffAll <- RelBias.results(result.diffAll, result.CLCOX.diffAll)

all.relbias <- c(unlist(homo), unlist(diffX), unlist(diffT), unlist(diffC), unlist(diffAll))
n <- length(all.relbias)

time <- rep(rep(t, times = M), times = 5 * 2 * 6)

Case <- factor(
  c(rep("Homogeneous", n/5), rep("Covariate Shift", n/5), rep("Outcome Shift", n/5), 
    rep("Censoring Shift", n/5), rep("All Shift", n/5)),
  levels = c("Homogeneous", "Covariate Shift", "Outcome Shift", "Censoring Shift", "All Shift")
)

Method <- factor(
  rep(c(rep("TGT", length(t)*M), rep("IVW", length(t)*M), rep("POOL", length(t)*M), 
        rep("CCOD", length(t)*M), rep("FED", length(t)*M), rep("CLCOX", length(t)*M)), 10),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/10), rep("Control", n/10)), 5)

df <- data.frame(RelBias = all.relbias, Case = Case, Method = Method, time = time, Treatment = Treatment)

df_heatmap <- df %>%
  group_by(Case, Method, Treatment, time) %>%
  summarize(RelBias = abs(mean(RelBias * 100)), .groups = 'drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RelBias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RelBias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "white", "green4"),
    values=scales::rescale(c(0, 0.05, 1)),
    name="ARBias%"
  ) +
  facet_grid(Treatment~Case) +
  labs(x="Time (day)", y="") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.relbias

pdf(file = "relbias_limO.pdf", width = 9.5, height = 4.5)
p.relbias
dev.off()



load("Res_imbal_500.Rdata")
load("Res_CLCOX_imbal.Rdata")
all.bias <- unlist(RelBias.results(results, result.CLCOX.imbal, M=500))
n <- length(all.bias)
time <- rep(rep(t, times = M), times = 2*6)

Method <- factor(
  rep(c(rep("TGT", length(t)*M), rep("IVW", length(t)*M), rep("POOL", length(t)*M), 
        rep("CCOD", length(t)*M), rep("FED", length(t)*M), rep("CLCOX", length(t)*M)), 2),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/2), rep("Control", n/2)), 1)

df <- data.frame(RelBias = all.bias, Method = Method, time = time, Treatment = Treatment)

df_heatmap <- df %>%
  group_by(Method, Treatment, time) %>%
  summarize(RelBias = abs(mean(RelBias * 100, na.rm=T)), .groups = 'drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RelBias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RelBias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "white", "green4"),
    values=scales::rescale(c(0, 0.05, 1)),
    name="ARBias%"
  ) +
  facet_grid(Treatment~.) +
  labs(x="Time (day)", y="") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.rbias

pdf("relbias_imbal.pdf", width=3.5, height=5)
p.rbias
dev.off()



load("Res_moreK_500.Rdata")
load("Res_CLCOX_imbal.Rdata")
all.bias <- unlist(RelBias.results(results, result.CLCOX.moreK, M=500))
n <- length(all.bias)
time <- rep(rep(t, times = M), times = 2*6)

Method <- factor(
  rep(c(rep("TGT", length(t)*M), rep("IVW", length(t)*M), rep("POOL", length(t)*M), 
        rep("CCOD", length(t)*M), rep("FED", length(t)*M), rep("CLCOX", length(t)*M)), 2),
  levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
)

Treatment <- rep(c(rep("Treated", n/2), rep("Control", n/2)), 1)

df <- data.frame(RelBias = all.bias, Method = Method, time = time, Treatment = Treatment)

df_heatmap <- df %>%
  group_by(Method, Treatment, time) %>%
  summarize(RelBias = abs(mean(RelBias * 100, na.rm=T)), .groups = 'drop')

ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RelBias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RelBias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white", "white", "green4"),
    values=scales::rescale(c(0, 0.05, 1)),
    name="ARBias%"
  ) +
  facet_grid(Treatment~.) +
  labs(x="Time (day)", y="") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.rbias

pdf("relbias_moreK.pdf", width=3.5, height=5)
p.rbias
dev.off()
