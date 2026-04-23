load("truth.Rdata")
RD.true[c(30,60)]
SR.true[c(30,60)]
RMST.0.true
RMST.1.true
RMST.diff.true

library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyverse)
library(RColorBrewer)

all.bias.list <- function(results, M=500) {
  TGT.RD <- CCOD.RD <- FED.RD <- c()
  TGT.SR <- CCOD.SR <- FED.SR <- c()
  TGT.RMST.1 <- CCOD.RMST.1 <- FED.RMST.1 <- c()
  TGT.RMST.0 <- CCOD.RMST.0 <- FED.RMST.0 <- c()
  TGT.RMST.diff <- CCOD.RMST.diff <- FED.RMST.diff <- c()
  for(j in 1:M) {
    TGT.RD <- c(TGT.RD, results[[j]]$df.RD.TGT[,"RD"]-RD.true)
    CCOD.RD <- c(CCOD.RD, results[[j]]$df.RD.CCOD[,"RD"]-RD.true)
    FED.RD <- c(FED.RD, results[[j]]$df.RD.FED[,"RD"]-RD.true)
    
    TGT.SR <- c(TGT.SR, results[[j]]$df.SR.TGT[,"SR"]-SR.true)
    CCOD.SR <- c(CCOD.SR, results[[j]]$df.SR.CCOD[,"SR"]-SR.true)
    FED.SR <- c(FED.SR, results[[j]]$df.SR.FED[,"SR"]-SR.true)
    
    TGT.RMST.1 <- c(TGT.RMST.1, results[[j]]$df.RMST.1.TGT[,"RMST"]-RMST.1.true)
    CCOD.RMST.1 <- c(CCOD.RMST.1, results[[j]]$df.RMST.1.CCOD[,"RMST"]-RMST.1.true)
    FED.RMST.1 <- c(FED.RMST.1, results[[j]]$df.RMST.1.FED[,"RMST"]-RMST.1.true)
    
    TGT.RMST.0 <- c(TGT.RMST.0, results[[j]]$df.RMST.0.TGT[,"RMST"]-RMST.0.true)
    CCOD.RMST.0 <- c(CCOD.RMST.0, results[[j]]$df.RMST.0.CCOD[,"RMST"]-RMST.0.true)
    FED.RMST.0 <- c(FED.RMST.0, results[[j]]$df.RMST.0.FED[,"RMST"]-RMST.0.true)
    
    TGT.RMST.diff <- c(TGT.RMST.diff, results[[j]]$df.RMST.diff.TGT[,"RMST"]-RMST.diff.true)
    CCOD.RMST.diff <- c(CCOD.RMST.diff, results[[j]]$df.RMST.diff.CCOD[,"RMST"]-RMST.diff.true)
    FED.RMST.diff <- c(FED.RMST.diff, results[[j]]$df.RMST.diff.FED[,"RMST"]-RMST.diff.true)
  }
  return(list(TGT.RD=TGT.RD, CCOD.RD=CCOD.RD, FED.RD=FED.RD, 
              TGT.SR=TGT.SR, CCOD.SR=CCOD.SR, FED.SR=FED.SR, 
              TGT.RMST.0=TGT.RMST.0, CCOD.RMST.0=CCOD.RMST.0, FED.RMST.0=FED.RMST.0,
              TGT.RMST.1=TGT.RMST.1, CCOD.RMST.1=CCOD.RMST.1, FED.RMST.1=FED.RMST.1,
              TGT.RMST.diff=TGT.RMST.diff, CCOD.RMST.diff=CCOD.RMST.diff, FED.RMST.diff=FED.RMST.diff))
}

RRMSE.results <- function(results, M=500) {
  
  TGT.RD <- CCOD.RD <- FED.RD <- c()
  TGT.SR <- CCOD.SR <- FED.SR <- c()
  for(j in 1:M) {
    TGT.RD <- c(TGT.RD, (results[[j]]$df.RD.TGT[,"RD"]-RD.true)^2)
    CCOD.RD <- c(CCOD.RD, (results[[j]]$df.RD.CCOD[,"RD"]-RD.true)^2)
    FED.RD <- c(FED.RD, (results[[j]]$df.RD.FED[,"RD"]-RD.true)^2)
    
    TGT.SR <- c(TGT.SR, (results[[j]]$df.SR.TGT[,"SR"]-SR.true)^2)
    CCOD.SR <- c(CCOD.SR, (results[[j]]$df.SR.CCOD[,"SR"]-SR.true)^2)
    FED.SR <- c(FED.SR, (results[[j]]$df.SR.FED[,"SR"]-SR.true)^2)
  }
  TGT.RD <- matrix(TGT.RD, nrow=length(results[[1]]$df.RD.TGT$time))
  CCOD.RD <- matrix(CCOD.RD, nrow=length(results[[1]]$df.RD.TGT$time))
  FED.RD <- matrix(FED.RD, nrow=length(results[[1]]$df.RD.TGT$time))
  
  TGT.SR <- matrix(TGT.SR, nrow=length(results[[1]]$df.RD.TGT$time))
  CCOD.SR <- matrix(CCOD.SR, nrow=length(results[[1]]$df.RD.TGT$time))
  FED.SR <- matrix(FED.SR, nrow=length(results[[1]]$df.RD.TGT$time))
  
  TGT.RMST.1 <- CCOD.RMST.1 <- FED.RMST.1 <- c()
  TGT.RMST.0 <- CCOD.RMST.0 <- FED.RMST.0 <- c()
  TGT.RMST.diff <- CCOD.RMST.diff <- FED.RMST.diff <- c()
  for(j in 1:M) {
    TGT.RMST.1 <- c(TGT.RMST.1, (results[[j]]$df.RMST.1.TGT[,"RMST"]-RMST.1.true)^2)
    CCOD.RMST.1 <- c(CCOD.RMST.1, (results[[j]]$df.RMST.1.CCOD[,"RMST"]-RMST.1.true)^2)
    FED.RMST.1 <- c(FED.RMST.1, (results[[j]]$df.RMST.1.FED[,"RMST"]-RMST.1.true)^2)
    
    TGT.RMST.0 <- c(TGT.RMST.0, (results[[j]]$df.RMST.0.TGT[,"RMST"]-RMST.0.true)^2)
    CCOD.RMST.0 <- c(CCOD.RMST.0, (results[[j]]$df.RMST.0.CCOD[,"RMST"]-RMST.0.true)^2)
    FED.RMST.0 <- c(FED.RMST.0, (results[[j]]$df.RMST.0.FED[,"RMST"]-RMST.0.true)^2)
    
    TGT.RMST.diff <- c(TGT.RMST.diff, (results[[j]]$df.RMST.diff.TGT[,"RMST"]-RMST.diff.true)^2)
    CCOD.RMST.diff <- c(CCOD.RMST.diff, (results[[j]]$df.RMST.diff.CCOD[,"RMST"]-RMST.diff.true)^2)
    FED.RMST.diff <- c(FED.RMST.diff, (results[[j]]$df.RMST.diff.FED[,"RMST"]-RMST.diff.true)^2)
  }
  
  rrmse.TGT.RD <- c(TGT.RD/apply(TGT.RD,1,median))
  rrmse.CCOD.RD <- c(CCOD.RD/apply(TGT.RD,1,median))
  rrmse.FED.RD <- c(FED.RD/apply(TGT.RD,1,median))
  
  rrmse.TGT.SR <- c(TGT.SR/apply(TGT.SR,1,median))
  rrmse.CCOD.SR <- c(CCOD.SR/apply(TGT.SR,1,median))
  rrmse.FED.SR <- c(FED.SR/apply(TGT.SR,1,median))
  
  rrmse.TGT.RMST.1 <- TGT.RMST.1/mean(TGT.RMST.1)
  rrmse.CCOD.RMST.1 <- CCOD.RMST.1/mean(TGT.RMST.1)
  rrmse.FED.RMST.1 <- FED.RMST.1/mean(TGT.RMST.1)

  rrmse.TGT.RMST.0 <- TGT.RMST.0/mean(TGT.RMST.0)
  rrmse.CCOD.RMST.0 <- CCOD.RMST.0/mean(TGT.RMST.0)
  rrmse.FED.RMST.0 <- FED.RMST.0/mean(TGT.RMST.0)
  
  rrmse.TGT.RMST.diff <- TGT.RMST.diff/mean(TGT.RMST.diff)
  rrmse.CCOD.RMST.diff <- CCOD.RMST.diff/mean(TGT.RMST.diff)
  rrmse.FED.RMST.diff <- FED.RMST.diff/mean(TGT.RMST.diff)
  
  return(list(rrmse.TGT.RD=rrmse.TGT.RD, rrmse.CCOD.RD=rrmse.CCOD.RD, 
              rrmse.FED.RD=rrmse.FED.RD,
              
              rrmse.TGT.SR=rrmse.TGT.SR, rrmse.CCOD.SR=rrmse.CCOD.SR, 
              rrmse.FED.SR=rrmse.FED.SR,
              
              rrmse.TGT.RMST.0=rrmse.TGT.RMST.0, rrmse.CCOD.RMST.0=rrmse.CCOD.RMST.0, 
              rrmse.FED.RMST.0=rrmse.FED.RMST.0,
              
              rrmse.TGT.RMST.1=rrmse.TGT.RMST.1, rrmse.CCOD.RMST.1=rrmse.CCOD.RMST.1, 
              rrmse.FED.RMST.1=rrmse.FED.RMST.1,
              
              rrmse.TGT.RMST.diff=rrmse.TGT.RMST.diff, rrmse.CCOD.RMST.diff=rrmse.CCOD.RMST.diff, 
              rrmse.FED.RMST.diff=rrmse.FED.RMST.diff))
}

CP.results <- function(results, M=500, quant=1.96) {
  cp.TGT.RD <- cp.CCOD.RD <- cp.FED.RD <- c()
  cp.TGT.SR <- cp.CCOD.SR <- cp.FED.SR <- c()
  cp.TGT.RMST.0 <- cp.CCOD.RMST.0 <- cp.FED.RMST.0 <- c()
  cp.TGT.RMST.1 <- cp.CCOD.RMST.1 <- cp.FED.RMST.1 <- c()
  cp.TGT.RMST.diff <- cp.CCOD.RMST.diff <- cp.FED.RMST.diff <- c()

  for(i in 1:length(t)) {
    TGT.RD <- CCOD.RD <- FED.RD <- matrix(NA, ncol=2, nrow=M)
    TGT.SR <- CCOD.SR <- FED.SR <- matrix(NA, ncol=2, nrow=M)
    for(j in 1:M) {
      TGT.RD[j,] <- results[[j]]$df.RD.TGT[,c("RD","sd")][i,] %>% as.numeric()
      CCOD.RD[j,] <- results[[j]]$df.RD.CCOD[,c("RD","sd")][i,] %>% as.numeric()
      FED.RD[j,] <- results[[j]]$df.RD.FED[,c("RD","sd")][i,] %>% as.numeric()
      
      TGT.SR[j,] <- results[[j]]$df.SR.TGT[,c("SR","sd")][i,] %>% as.numeric()
      CCOD.SR[j,] <- results[[j]]$df.SR.CCOD[,c("SR","sd")][i,] %>% as.numeric()
      FED.SR[j,] <- results[[j]]$df.SR.FED[,c("SR","sd")][i,] %>% as.numeric()
    }
    cp.TGT.RD[i] <- mean((RD.true[i]<TGT.RD[,1] + quant*TGT.RD[,2]) & (RD.true[i]>TGT.RD[,1] - quant*TGT.RD[,2]), na.rm=T)
    cp.CCOD.RD[i] <- mean((RD.true[i]<CCOD.RD[,1] + quant*CCOD.RD[,2]) & (RD.true[i]>CCOD.RD[,1] - quant*CCOD.RD[,2]), na.rm=T)
    cp.FED.RD[i] <- mean((RD.true[i]<FED.RD[,1] + quant*FED.RD[,2]) & (RD.true[i]>FED.RD[,1] - quant*FED.RD[,2]), na.rm=T)
    
    cp.TGT.SR[i] <- mean((SR.true[i]<TGT.SR[,1] + quant*TGT.SR[,2]) & (SR.true[i]>TGT.SR[,1] - quant*TGT.SR[,2]), na.rm=T)
    cp.CCOD.SR[i] <- mean((SR.true[i]<CCOD.SR[,1] + quant*CCOD.SR[,2]) & (SR.true[i]>CCOD.SR[,1] - quant*CCOD.SR[,2]), na.rm=T)
    cp.FED.SR[i] <- mean((SR.true[i]<FED.SR[,1] + quant*FED.SR[,2]) & (SR.true[i]>FED.SR[,1] - quant*FED.SR[,2]), na.rm=T)
  }
  
  TGT.RMST.0 <- CCOD.RMST.0 <- FED.RMST.0 <- matrix(NA, ncol=2, nrow=M)
  TGT.RMST.1 <- CCOD.RMST.1 <- FED.RMST.1 <- matrix(NA, ncol=2, nrow=M)
  TGT.RMST.diff <- CCOD.RMST.diff <- FED.RMST.diff <- matrix(NA, ncol=2, nrow=M)
  for(j in 1:M) {
    TGT.RMST.0[j,] <- results[[j]]$df.RMST.0.TGT[,c("RMST","sd")] %>% as.numeric()
    CCOD.RMST.0[j,] <- results[[j]]$df.RMST.0.CCOD[,c("RMST","sd")] %>% as.numeric()
    FED.RMST.0[j,] <- results[[j]]$df.RMST.0.FED[,c("RMST","sd")] %>% as.numeric()
    
    TGT.RMST.1[j,] <- results[[j]]$df.RMST.1.TGT[,c("RMST","sd")] %>% as.numeric()
    CCOD.RMST.1[j,] <- results[[j]]$df.RMST.1.CCOD[,c("RMST","sd")] %>% as.numeric()
    FED.RMST.1[j,] <- results[[j]]$df.RMST.1.FED[,c("RMST","sd")] %>% as.numeric()
    
    TGT.RMST.diff[j,] <- results[[j]]$df.RMST.diff.TGT[,c("RMST","sd")] %>% as.numeric()
    CCOD.RMST.diff[j,] <- results[[j]]$df.RMST.diff.CCOD[,c("RMST","sd")] %>% as.numeric()
    FED.RMST.diff[j,] <- results[[j]]$df.RMST.diff.FED[,c("RMST","sd")] %>% as.numeric()
  }
  
  RMST.0.true <- as.numeric(RMST.0.true)
  cp.TGT.RMST.0 <- mean((RMST.0.true<TGT.RMST.0[,1] + quant*TGT.RMST.0[,2]) & (RMST.0.true>TGT.RMST.0[,1] - quant*TGT.RMST.0[,2]), na.rm=T)
  cp.CCOD.RMST.0 <- mean((RMST.0.true<CCOD.RMST.0[,1] + quant*CCOD.RMST.0[,2]) & (RMST.0.true>CCOD.RMST.0[,1] - quant*CCOD.RMST.0[,2]), na.rm=T)
  cp.FED.RMST.0 <- mean((RMST.0.true<FED.RMST.0[,1] + quant*FED.RMST.0[,2]) & (RMST.0.true>FED.RMST.0[,1] - quant*FED.RMST.0[,2]), na.rm=T)
  
  RMST.1.true <- as.numeric(RMST.1.true)
  cp.TGT.RMST.1 <- mean((RMST.1.true<TGT.RMST.1[,1] + quant*TGT.RMST.1[,2]) & (RMST.1.true>TGT.RMST.1[,1] - quant*TGT.RMST.1[,2]), na.rm=T)
  cp.CCOD.RMST.1 <- mean((RMST.1.true<CCOD.RMST.1[,1] + quant*CCOD.RMST.1[,2]) & (RMST.1.true>CCOD.RMST.1[,1] - quant*CCOD.RMST.1[,2]), na.rm=T)
  cp.FED.RMST.1 <- mean((RMST.1.true<FED.RMST.1[,1] + quant*FED.RMST.1[,2]) & (RMST.1.true>FED.RMST.1[,1] - quant*FED.RMST.1[,2]), na.rm=T)
  
  RMST.diff.true <- as.numeric(RMST.diff.true)
  cp.TGT.RMST.diff <- mean((RMST.diff.true<TGT.RMST.diff[,1] + quant*TGT.RMST.diff[,2]) & (RMST.diff.true>TGT.RMST.diff[,1] - quant*TGT.RMST.diff[,2]), na.rm=T)
  cp.CCOD.RMST.diff <- mean((RMST.diff.true<CCOD.RMST.diff[,1] + quant*CCOD.RMST.diff[,2]) & (RMST.diff.true>CCOD.RMST.diff[,1] - quant*CCOD.RMST.diff[,2]), na.rm=T)
  cp.FED.RMST.diff <- mean((RMST.diff.true<FED.RMST.diff[,1] + quant*FED.RMST.diff[,2]) & (RMST.diff.true>FED.RMST.diff[,1] - quant*FED.RMST.diff[,2]), na.rm=T)
  
  return(list(cp.TGT.RD=cp.TGT.RD, cp.CCOD.RD=cp.CCOD.RD, cp.FED.RD=cp.FED.RD,
              cp.TGT.SR=cp.TGT.SR, cp.CCOD.SR=cp.CCOD.SR, cp.FED.SR=cp.FED.SR,
              cp.TGT.RMST.0=cp.TGT.RMST.0, cp.CCOD.RMST.0=cp.CCOD.RMST.0, cp.FED.RMST.0=cp.FED.RMST.0,
              cp.TGT.RMST.1=cp.TGT.RMST.1, cp.CCOD.RMST.1=cp.CCOD.RMST.1, cp.FED.RMST.1=cp.FED.RMST.1,
              cp.TGT.RMST.diff=cp.TGT.RMST.diff, cp.CCOD.RMST.diff=cp.CCOD.RMST.diff, cp.FED.RMST.diff=cp.FED.RMST.diff))
}

### Relative bias
M <- 500

load("homo_500.Rdata")
homo <- all.bias.list(results)
load("diffX_500.Rdata")
diffX <- all.bias.list(results)
load("diffT_500.Rdata")
diffT <- all.bias.list(results)
load("diffC_500.Rdata")
diffC <- all.bias.list(results)
load("diffAll_500.Rdata")
diffAll <- all.bias.list(results)

all.bias.RD <- c(homo$TGT.RD, homo$CCOD.RD, homo$FED.RD,
                 diffX$TGT.RD, diffX$CCOD.RD, diffX$FED.RD,
                 diffT$TGT.RD, diffT$CCOD.RD, diffT$FED.RD,
                 diffC$TGT.RD, diffC$CCOD.RD, diffC$FED.RD,
                 diffAll$TGT.RD, diffAll$CCOD.RD, diffAll$FED.RD) / RD.true
n.RD <- length(all.bias.RD)
t <- 1:60
time <- rep(t, n.RD/length(t))
Case <- factor(c(rep("Homogeneous", n.RD/5), rep("Covariate Shift", n.RD/5), 
                 rep("Outcome Shift", n.RD/5), rep("Censoring Shift", n.RD/5), 
                 rep("All Shifts", n.RD/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", 
                        "Censoring Shift", "All Shifts"))
Method <- factor(rep(c(rep("TGT", length(t)*M), rep("CCOD", length(t)*M), rep("FED", length(t)*M)), 5), 
                 levels=c("FED", "TGT", "CCOD"))

df <- data.frame(Bias=all.bias.RD, Case=Case, Method=Method, time=time) %>% 
  filter(time%in%c(30,60))
df_heatmap <- df %>%
  group_by(Case, Method, time) %>% summarize(Bias=abs(mean(Bias*100)), .groups='drop')
ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=Bias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(Bias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "green4"), 
    values=scales::rescale(c(0, 5, 100)),
    name="ARBias%"
  ) +
  facet_grid(~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.bias.RD

pdf(file="bias_RD.pdf", width=7, height=1.8)
p.bias.RD
dev.off()

all.bias.SR <- c(homo$TGT.SR, homo$CCOD.SR, homo$FED.SR,
                 diffX$TGT.SR, diffX$CCOD.SR, diffX$FED.SR,
                 diffT$TGT.SR, diffT$CCOD.SR, diffT$FED.SR,
                 diffC$TGT.SR, diffC$CCOD.SR, diffC$FED.SR,
                 diffAll$TGT.SR, diffAll$CCOD.SR, diffAll$FED.SR) / SR.true
n.SR <- length(all.bias.SR)
t <- 1:60
time <- rep(t, n.SR/length(t))
Case <- factor(c(rep("Homogeneous", n.SR/5), rep("Covariate Shift", n.SR/5), 
                 rep("Outcome Shift", n.SR/5), rep("Censoring Shift", n.SR/5), 
                 rep("All Shifts", n.SR/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", 
                        "Censoring Shift", "All Shifts"))
Method <- factor(rep(c(rep("TGT", length(t)*M), rep("CCOD", length(t)*M), rep("FED", length(t)*M)), 5), 
                 levels=c("FED", "TGT", "CCOD"))

df <- data.frame(Bias=all.bias.SR, Case=Case, Method=Method, time=time) %>% filter(time%in%c(30,60))
df_heatmap <- df %>%
  group_by(Case, Method, time) %>% summarize(Bias=abs(mean(Bias*100)), .groups='drop')
ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=Bias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(Bias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "green4"), 
    values=scales::rescale(c(0, 5, 100)),
    name="ARBias%"
  ) +
  facet_grid(~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.bias.SR

pdf(file="bias_SR.pdf", width=7, height=1.8)
p.bias.SR
dev.off()

all.bias.RMST <- c(homo$TGT.RMST.0 / RMST.0.true, homo$TGT.RMST.1 / RMST.1.true, homo$TGT.RMST.diff / RMST.diff.true,
                   homo$CCOD.RMST.0 / RMST.0.true, homo$CCOD.RMST.1 / RMST.1.true, homo$CCOD.RMST.diff / RMST.diff.true,
                   homo$FED.RMST.0 / RMST.0.true, homo$FED.RMST.1 / RMST.1.true, homo$FED.RMST.diff / RMST.diff.true,
                   
                   diffX$TGT.RMST.0 / RMST.0.true, diffX$TGT.RMST.1 / RMST.1.true, diffX$TGT.RMST.diff / RMST.diff.true,
                   diffX$CCOD.RMST.0 / RMST.0.true, diffX$CCOD.RMST.1 / RMST.1.true, diffX$CCOD.RMST.diff / RMST.diff.true,
                   diffX$FED.RMST.0 / RMST.0.true, diffX$FED.RMST.1 / RMST.1.true, diffX$FED.RMST.diff / RMST.diff.true,
                   
                   diffT$TGT.RMST.0 / RMST.0.true, diffT$TGT.RMST.1 / RMST.1.true, diffT$TGT.RMST.diff / RMST.diff.true,
                   diffT$CCOD.RMST.0 / RMST.0.true, diffT$CCOD.RMST.1 / RMST.1.true, diffT$CCOD.RMST.diff / RMST.diff.true,
                   diffT$FED.RMST.0 / RMST.0.true, diffT$FED.RMST.1 / RMST.1.true, diffT$FED.RMST.diff / RMST.diff.true,
                   
                   diffC$TGT.RMST.0 / RMST.0.true, diffC$TGT.RMST.1 / RMST.1.true, diffC$TGT.RMST.diff / RMST.diff.true,
                   diffC$CCOD.RMST.0 / RMST.0.true, diffC$CCOD.RMST.1 / RMST.1.true, diffC$CCOD.RMST.diff / RMST.diff.true,
                   diffC$FED.RMST.0 / RMST.0.true, diffC$FED.RMST.1 / RMST.1.true, diffC$FED.RMST.diff / RMST.diff.true,
                   
                   diffAll$TGT.RMST.0 / RMST.0.true, diffAll$TGT.RMST.1 / RMST.1.true, diffAll$TGT.RMST.diff / RMST.diff.true,
                   diffAll$CCOD.RMST.0 / RMST.0.true, diffAll$CCOD.RMST.1 / RMST.1.true, diffAll$CCOD.RMST.diff / RMST.diff.true,
                   diffAll$FED.RMST.0 / RMST.0.true, diffAll$FED.RMST.1 / RMST.1.true, diffAll$FED.RMST.diff / RMST.diff.true)

n.RMST <- length(all.bias.RMST)
Case <- factor(c(rep("Homogeneous", n.RMST/5), rep("Covariate Shift", n.RMST/5), 
                 rep("Outcome Shift", n.RMST/5), rep("Censoring Shift", n.RMST/5), 
                 rep("All Shifts", n.RMST/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", 
                        "Censoring Shift", "All Shifts"))
Method <- factor(rep(c(rep("TGT", 3*M), rep("CCOD", 3*M), rep("FED", 3*M)), 5), 
                 levels=c("FED", "TGT", "CCOD"))
Measure <- factor(rep(c(rep("A = 0", M), rep("A = 1", M), rep("Difference", M)), 15),
                 levels=c("A = 0", "A = 1", "Difference"))

df <- data.frame(Bias=all.bias.RMST, Case=Case, Method=Method, Measure=Measure) 
df_heatmap <- df %>%
  group_by(Case, Method, Measure) %>% summarize(Bias=abs(mean(Bias*100)), .groups='drop')
ggplot(df_heatmap, aes(x=Measure, y=Method, fill=Bias)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(Bias, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "green4"), 
    values=scales::rescale(c(0, 5, 100)),
    name="ARBias%"
  ) +
  facet_grid(~Case) + labs(x="", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), axis.text.x=element_blank()) -> p.bias.RMST

pdf(file="bias_RMST.pdf", width=8, height=2)
p.bias.RMST
dev.off()

### RRMSE
load("homo_500.Rdata")
homo <- RRMSE.results(results)
load("diffX_500.Rdata")
diffX <- RRMSE.results(results)
load("diffT_500.Rdata")
diffT <- RRMSE.results(results)
load("diffC_500.Rdata")
diffC <- RRMSE.results(results)
load("diffAll_500.Rdata")
diffAll <- RRMSE.results(results)

all.rrmse.RD <- c(homo$rrmse.TGT.RD, homo$rrmse.CCOD.RD, homo$rrmse.FED.RD,
                  diffX$rrmse.TGT.RD, diffX$rrmse.CCOD.RD, diffX$rrmse.FED.RD,
                  diffT$rrmse.TGT.RD, diffT$rrmse.CCOD.RD, diffT$rrmse.FED.RD,
                  diffC$rrmse.TGT.RD, diffC$rrmse.CCOD.RD, diffC$rrmse.FED.RD,
                  diffAll$rrmse.TGT.RD, diffAll$rrmse.CCOD.RD, diffAll$rrmse.FED.RD)
n.RD <- length(all.rrmse.RD)
t <- 1:60
time <- rep(t, n.RD/length(t))
Case <- factor(c(rep("Homogeneous", n.RD/5), rep("Covariate Shift", n.RD/5), 
                 rep("Outcome Shift", n.RD/5), rep("Censoring Shift", n.RD/5), 
                 rep("All Shifts", n.RD/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", 
                        "Censoring Shift", "All Shifts"))
Method <- factor(rep(c(rep("TGT", length(t)*M), rep("CCOD", length(t)*M), 
                       rep("FED", length(t)*M)), 5), 
                 levels=c("FED", "TGT", "CCOD"))

df <- data.frame(RRMSE=all.rrmse.RD, Case=Case, Method=Method, time=time) %>% 
  filter(time%in%c(30,60))
df_heatmap <- df %>%
  group_by(Case, Method, time) %>% summarize(RRMSE=median(RRMSE), .groups='drop')
ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RRMSE)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RRMSE,2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "orange"), 
    values=scales::rescale(c(0, 1, 2.5)),
    name="RRMSE"
  ) +
  facet_grid(~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.rrmse.RD

pdf(file="rrmse_RD.pdf", width=7, height=1.8)
p.rrmse.RD
dev.off()

all.rrmse.SR <- c(homo$rrmse.TGT.SR, homo$rrmse.CCOD.SR, homo$rrmse.FED.SR,
                  diffX$rrmse.TGT.SR, diffX$rrmse.CCOD.SR, diffX$rrmse.FED.SR,
                  diffT$rrmse.TGT.SR, diffT$rrmse.CCOD.SR, diffT$rrmse.FED.SR,
                  diffC$rrmse.TGT.SR, diffC$rrmse.CCOD.SR, diffC$rrmse.FED.SR,
                  diffAll$rrmse.TGT.SR, diffAll$rrmse.CCOD.SR, diffAll$rrmse.FED.SR)
n.RD <- length(all.rrmse.SR)
t <- 1:60
time <- rep(t, n.SR/length(t))
Case <- factor(c(rep("Homogeneous", n.SR/5), rep("Covariate Shift", n.SR/5), 
                 rep("Outcome Shift", n.SR/5), rep("Censoring Shift", n.SR/5), 
                 rep("All Shifts", n.SR/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", 
                        "Censoring Shift", "All Shifts"))
Method <- factor(rep(c(rep("TGT", length(t)*M), rep("CCOD", length(t)*M), 
                       rep("FED", length(t)*M)), 5), 
                 levels=c("FED", "TGT", "CCOD"))

df <- data.frame(RRMSE=all.rrmse.SR, Case=Case, Method=Method, time=time) %>% 
  filter(time%in%c(30,60))
df_heatmap <- df %>%
  group_by(Case, Method, time) %>% summarize(RRMSE=median(RRMSE), .groups='drop')
ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=RRMSE)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RRMSE,2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "orange"), 
    values=scales::rescale(c(0, 1, 1.35)),
    name="RRMSE"
  ) +
  facet_grid(~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.rrmse.SR

pdf(file="rrmse_SR.pdf", width=7, height=1.8)
p.rrmse.SR
dev.off()

all.rrmse.RMST <- c(homo$rrmse.TGT.RMST.0, homo$rrmse.TGT.RMST.1, homo$rrmse.TGT.RMST.diff,
                    homo$rrmse.CCOD.RMST.0, homo$rrmse.CCOD.RMST.1, homo$rrmse.CCOD.RMST.diff,
                    homo$rrmse.FED.RMST.0, homo$rrmse.FED.RMST.1, homo$rrmse.FED.RMST.diff,
                    
                    diffX$rrmse.TGT.RMST.0, diffX$rrmse.TGT.RMST.1, diffX$rrmse.TGT.RMST.diff,
                    diffX$rrmse.CCOD.RMST.0, diffX$rrmse.CCOD.RMST.1, diffX$rrmse.CCOD.RMST.diff,
                    diffX$rrmse.FED.RMST.0, diffX$rrmse.FED.RMST.1, diffX$rrmse.FED.RMST.diff,
                    
                    diffT$rrmse.TGT.RMST.0, diffT$rrmse.TGT.RMST.1, diffT$rrmse.TGT.RMST.diff,
                    diffT$rrmse.CCOD.RMST.0, diffT$rrmse.CCOD.RMST.1, diffT$rrmse.CCOD.RMST.diff,
                    diffT$rrmse.FED.RMST.0, diffT$rrmse.FED.RMST.1, diffT$rrmse.FED.RMST.diff,
                    
                    diffC$rrmse.TGT.RMST.0, diffC$rrmse.TGT.RMST.1, diffC$rrmse.TGT.RMST.diff,
                    diffC$rrmse.CCOD.RMST.0, diffC$rrmse.CCOD.RMST.1, diffC$rrmse.CCOD.RMST.diff,
                    diffC$rrmse.FED.RMST.0, diffC$rrmse.FED.RMST.1, diffC$rrmse.FED.RMST.diff,
                    
                    diffAll$rrmse.TGT.RMST.0, diffAll$rrmse.TGT.RMST.1, diffAll$rrmse.TGT.RMST.diff,
                    diffAll$rrmse.CCOD.RMST.0, diffAll$rrmse.CCOD.RMST.1, diffAll$rrmse.CCOD.RMST.diff,
                    diffAll$rrmse.FED.RMST.0, diffAll$rrmse.FED.RMST.1, diffAll$rrmse.FED.RMST.diff)

n.RMST <- length(all.rrmse.RMST)
Case <- factor(c(rep("Homogeneous", n.RMST/5), rep("Covariate Shift", n.RMST/5), 
                 rep("Outcome Shift", n.RMST/5), rep("Censoring Shift", n.RMST/5), 
                 rep("All Shifts", n.RMST/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", 
                        "Censoring Shift", "All Shifts"))
Method <- factor(rep(c(rep("TGT", 3*M), rep("CCOD", 3*M), rep("FED", 3*M)), 5), 
                 levels=c("FED", "TGT", "CCOD"))
Measure <- factor(rep(c(rep("A = 0", M), rep("A = 1", M), rep("Difference", M)), 15),
                  levels=c("A = 0", "A = 1", "Difference"))

df <- data.frame(RRMSE=all.rrmse.RMST, Case=Case, Method=Method, Measure=Measure) 
df_heatmap <- df %>%
  group_by(Case, Method, Measure) %>% summarize(RRMSE=mean(RRMSE), .groups='drop')
ggplot(df_heatmap, aes(x=Measure, y=Method, fill=RRMSE)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(RRMSE, 2)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("white","white", "orange"), 
    values=scales::rescale(c(0, 1, 100)),
    name="RRMSE"
  ) +
  facet_grid(~Case) + labs(x="", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), axis.text.x=element_blank()) -> p.rrmse.RMST

pdf(file="rrmse_RMST.pdf", width=8, height=2)
p.rrmse.RMST
dev.off()

### Coverage probability
load("homo_500.Rdata")
homo <- CP.results(results)
load("diffX_500.Rdata")
diffX <- CP.results(results)
load("diffT_500.Rdata")
diffT <- CP.results(results)
load("diffC_500.Rdata")
diffC <- CP.results(results)
load("diffAll_500.Rdata")
diffAll <- CP.results(results)

all.cp.RD <- c(homo$cp.TGT.RD, homo$cp.CCOD.RD, homo$cp.FED.RD,
                 diffX$cp.TGT.RD, diffX$cp.CCOD.RD, diffX$cp.FED.RD,
                 diffT$cp.TGT.RD, diffT$cp.CCOD.RD, diffT$cp.FED.RD,
                 diffC$cp.TGT.RD, diffC$cp.CCOD.RD, diffC$cp.FED.RD,
                 diffAll$cp.TGT.RD, diffAll$cp.CCOD.RD, diffAll$cp.FED.RD)
n.RD <- length(all.cp.RD)
t <- 1:60
time <- rep(t, n.RD/length(t))
Case <- factor(c(rep("Homogeneous", n.RD/5), rep("Covariate Shift", n.RD/5), 
                 rep("Outcome Shift", n.RD/5), rep("Censoring Shift", n.RD/5), 
                 rep("All Shifts", n.RD/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", 
                        "Censoring Shift", "All Shifts"))
Method <- factor(rep(c(rep("TGT", length(t)), rep("CCOD", length(t)), rep("FED", length(t))), 5), 
                 levels=c("FED", "TGT", "CCOD"))

df <- data.frame(CP=all.cp.RD, Case=Case, Method=Method, time=time) %>% 
  filter(time%in%c(30,60))
df_heatmap <- df %>%
  group_by(Case, Method, time) %>% summarize(CP=CP*100, .groups='drop')
ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=CP)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CP, 1)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("mediumvioletred", "white","mediumvioletred"), 
    values=scales::rescale(c(0, 87, 95, 99, 100)), # Soft transition around 95%
    name="CP%",
    limits=c(0, 100)
  ) +
  facet_grid(~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.cp.RD

pdf(file="cp_RD.pdf", width=7, height=1.8)
p.cp.RD
dev.off()

all.cp.SR <- c(homo$cp.TGT.SR, homo$cp.CCOD.SR, homo$cp.FED.SR,
               diffX$cp.TGT.SR, diffX$cp.CCOD.SR, diffX$cp.FED.SR,
               diffT$cp.TGT.SR, diffT$cp.CCOD.SR, diffT$cp.FED.SR,
               diffC$cp.TGT.SR, diffC$cp.CCOD.SR, diffC$cp.FED.SR,
               diffAll$cp.TGT.SR, diffAll$cp.CCOD.SR, diffAll$cp.FED.SR)
n.SR <- length(all.cp.SR)
t <- 1:60
time <- rep(t, n.SR/length(t))
Case <- factor(c(rep("Homogeneous", n.SR/5), rep("Covariate Shift", n.SR/5), 
                 rep("Outcome Shift", n.SR/5), rep("Censoring Shift", n.SR/5), 
                 rep("All Shifts", n.SR/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", 
                        "Censoring Shift", "All Shifts"))
Method <- factor(rep(c(rep("TGT", length(t)), rep("CCOD", length(t)), rep("FED", length(t))), 5), 
                 levels=c("FED", "TGT", "CCOD"))

df <- data.frame(CP=all.cp.SR, Case=Case, Method=Method, time=time) %>% 
  filter(time%in%c(30,60))
df_heatmap <- df %>%
  group_by(Case, Method, time) %>% summarize(CP=CP*100, .groups='drop')
ggplot(df_heatmap, aes(x=factor(time), y=Method, fill=CP)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CP, 1)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("mediumvioletred", "white","mediumvioletred"), 
    values=scales::rescale(c(0, 87, 95, 99, 100)), # Soft transition around 95%
    name="CP%",
    limits=c(0, 100)
  ) +
  facet_grid(~Case) + labs(x="Time (day)", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) -> p.cp.SR

pdf(file="cp_SR.pdf", width=7, height=1.8)
p.cp.SR
dev.off()

all.cp.RMST <- c(homo$cp.TGT.RMST.0, homo$cp.TGT.RMST.1, homo$cp.TGT.RMST.diff,
                   homo$cp.CCOD.RMST.0, homo$cp.CCOD.RMST.1, homo$cp.CCOD.RMST.diff,
                   homo$cp.FED.RMST.0, homo$cp.FED.RMST.1, homo$cp.FED.RMST.diff,
                   
                   diffX$cp.TGT.RMST.0, diffX$cp.TGT.RMST.1, diffX$cp.TGT.RMST.diff,
                   diffX$cp.CCOD.RMST.0, diffX$cp.CCOD.RMST.1, diffX$cp.CCOD.RMST.diff,
                   diffX$cp.FED.RMST.0, diffX$cp.FED.RMST.1, diffX$cp.FED.RMST.diff,
                   
                   diffT$cp.TGT.RMST.0, diffT$cp.TGT.RMST.1, diffT$cp.TGT.RMST.diff,
                   diffT$cp.CCOD.RMST.0, diffT$cp.CCOD.RMST.1, diffT$cp.CCOD.RMST.diff,
                   diffT$cp.FED.RMST.0, diffT$cp.FED.RMST.1, diffT$cp.FED.RMST.diff,
                   
                   diffC$cp.TGT.RMST.0, diffC$cp.TGT.RMST.1, diffC$cp.TGT.RMST.diff,
                   diffC$cp.CCOD.RMST.0, diffC$cp.CCOD.RMST.1, diffC$cp.CCOD.RMST.diff,
                   diffC$cp.FED.RMST.0, diffC$cp.FED.RMST.1, diffC$cp.FED.RMST.diff,
                   
                   diffAll$cp.TGT.RMST.0, diffAll$cp.TGT.RMST.1, diffAll$cp.TGT.RMST.diff,
                   diffAll$cp.CCOD.RMST.0, diffAll$cp.CCOD.RMST.1, diffAll$cp.CCOD.RMST.diff,
                   diffAll$cp.FED.RMST.0, diffAll$cp.FED.RMST.1, diffAll$cp.FED.RMST.diff)

n.RMST <- length(all.cp.RMST)
Case <- factor(c(rep("Homogeneous", n.RMST/5), rep("Covariate Shift", n.RMST/5), 
                 rep("Outcome Shift", n.RMST/5), rep("Censoring Shift", n.RMST/5), 
                 rep("All Shifts", n.RMST/5)), 
               levels=c("Homogeneous", "Covariate Shift", "Outcome Shift", 
                        "Censoring Shift", "All Shifts"))
Method <- factor(rep(c(rep("TGT", 3), rep("CCOD", 3), rep("FED", 3)), 5), 
                 levels=c("FED", "TGT", "CCOD"))
Measure <- factor(rep(c(rep("A = 0", 1), rep("A = 1", 1), rep("Difference", 1)), 15),
                  levels=c("A = 0", "A = 1", "Difference"))

df <- data.frame(CP=all.cp.RMST, Case=Case, Method=Method, Measure=Measure) 
df_heatmap <- df %>%
  group_by(Case, Method, Measure) %>% summarize(CP=CP*100, .groups='drop')
ggplot(df_heatmap, aes(x=Measure, y=Method, fill=CP)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(CP, 1)), size=3, color="black") +
  scale_fill_gradientn(
    colors=c("mediumvioletred", "white","mediumvioletred"), 
    values=scales::rescale(c(0, 87, 95, 99, 100)), # Soft transition around 95%
    name="CP%",
    limits=c(0, 100)
  ) +
  facet_grid(~Case) + labs(x="", y="") + theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1)) -> p.cp.RMST

pdf(file="cp_RMST.pdf", width=8, height=2.5)
p.cp.RMST
dev.off()
