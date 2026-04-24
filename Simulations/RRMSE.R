M <- 500
t <- c(30, 60, 90)

RRMSE.results <- function(datlist, coxlist, M = 500, t = c(30, 60, 90)) {
  
  methods <- c("TGT", "IVW", "POOL", "CCOD", "FED", "CLCOX")
  trts <- c("1", "0")
  
  se.array <- array(
    NA_real_,
    dim = c(M, length(t), length(methods), length(trts)),
    dimnames = list(
      sim = 1:M,
      time = t,
      method = methods,
      treatment = c("Treated", "Control")
    )
  )
  
  for (j in 1:M) {
    time.use <- datlist[[j]]$df.TGT$time
    df.cox <- coxlist[[j]]$df.CLCOX
    df.cox <- df.cox[match(time.use, df.cox$time), ]
    
    se.array[j, , "TGT", "Treated"]  <- (datlist[[j]]$df.TGT[, "surv1"]  - S1.true)^2
    se.array[j, , "IVW", "Treated"]  <- (datlist[[j]]$df.IVW[, "surv1"]  - S1.true)^2
    se.array[j, , "POOL", "Treated"] <- (datlist[[j]]$df.POOL[, "surv1"] - S1.true)^2
    se.array[j, , "CCOD", "Treated"] <- (datlist[[j]]$df.CCOD[, "surv1"] - S1.true)^2
    se.array[j, , "FED", "Treated"]  <- (datlist[[j]]$df.FED[, "surv1"]  - S1.true)^2
    se.array[j, , "CLCOX", "Treated"] <- (df.cox[, "surv1"] - S1.true)^2
    
    se.array[j, , "TGT", "Control"]  <- (datlist[[j]]$df.TGT[, "surv0"]  - S0.true)^2
    se.array[j, , "IVW", "Control"]  <- (datlist[[j]]$df.IVW[, "surv0"]  - S0.true)^2
    se.array[j, , "POOL", "Control"] <- (datlist[[j]]$df.POOL[, "surv0"] - S0.true)^2
    se.array[j, , "CCOD", "Control"] <- (datlist[[j]]$df.CCOD[, "surv0"] - S0.true)^2
    se.array[j, , "FED", "Control"]  <- (datlist[[j]]$df.FED[, "surv0"]  - S0.true)^2
    se.array[j, , "CLCOX", "Control"] <- (df.cox[, "surv0"] - S0.true)^2
  }
  
  rmse.array <- (apply(se.array, c(2, 3, 4), mean, na.rm = TRUE))
  
  out <- expand.grid(
    time = as.numeric(dimnames(rmse.array)$time),
    Method = dimnames(rmse.array)$method,
    Treatment = dimnames(rmse.array)$treatment,
    KEEP.OUT.ATTRS = FALSE
  )
  
  out$RMSE <- mapply(
    function(tt, mm, aa) rmse.array[as.character(tt), mm, aa],
    out$time, out$Method, out$Treatment
  )
  
  out$RMSE.TGT <- mapply(
    function(tt, aa) rmse.array[as.character(tt), "TGT", aa],
    out$time, out$Treatment
  )
  
  out$RRMSE <- out$RMSE / out$RMSE.TGT
  
  out$Method <- factor(
    out$Method,
    levels = c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
  )
  
  out
}

make_rrmse_df <- function(result.list, clcox.list, case.name, M = 500, t = c(30, 60, 90)) {
  df <- RRMSE.results(result.list, clcox.list, M = M, t = t)
  df$Case <- case.name
  df
}

plot_rrmse_heatmap <- function(df_heatmap) {
  ggplot(df_heatmap, aes(x = factor(time), y = Method, fill = RRMSE)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(RRMSE, 2)), size = 3, color = "black") +
    scale_fill_gradientn(
      colors = c("white", "white", "orange"),
      values = scales::rescale(c(0, 1, 100)),
      name = "RRMSE"
    ) +
    facet_grid(Treatment ~ Case) +
    labs(x = "Time (day)", y = "") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
}

library(ggplot2)
library(dplyr)

load("truth.Rdata")

## Small sample setting
load("Res_homo_s.Rdata")
load("Res_diffX_s.Rdata")
load("Res_diffT_s.Rdata")
load("Res_diffC_s.Rdata")
load("Res_diffAll_s.Rdata")
load("Res_CLCOX_s.Rdata")

df_heatmap <- bind_rows(
  make_rrmse_df(result.homo, result.CLCOX.homo, "Homogeneous", M, t),
  make_rrmse_df(result.diffX, result.CLCOX.diffX, "Covariate Shift", M, t),
  make_rrmse_df(result.diffT, result.CLCOX.diffT, "Outcome Shift", M, t),
  make_rrmse_df(result.diffC, result.CLCOX.diffC, "Censoring Shift", M, t),
  make_rrmse_df(result.diffAll, result.CLCOX.diffAll, "All Shift", M, t)
) %>%
  mutate(
    Case = factor(
      Case,
      levels = c("Homogeneous", "Covariate Shift", "Outcome Shift",
                 "Censoring Shift", "All Shift")
    )
  )

p.rrmse <- plot_rrmse_heatmap(df_heatmap)

pdf(file = "rrmse_s.pdf", width = 8.5, height = 4.5)
p.rrmse
dev.off()


## Large sample setting
load("Res_homo_l.Rdata")
load("Res_diffX_l.Rdata")
load("Res_diffT_l.Rdata")
load("Res_diffC_l.Rdata")
load("Res_diffAll_l.Rdata")
load("Res_CLCOX_l.Rdata")

df_heatmap <- bind_rows(
  make_rrmse_df(result.homo, result.CLCOX.homo, "Homogeneous", M, t),
  make_rrmse_df(result.diffX, result.CLCOX.diffX, "Covariate Shift", M, t),
  make_rrmse_df(result.diffT, result.CLCOX.diffT, "Outcome Shift", M, t),
  make_rrmse_df(result.diffC, result.CLCOX.diffC, "Censoring Shift", M, t),
  make_rrmse_df(result.diffAll, result.CLCOX.diffAll, "All Shift", M, t)
) %>%
  mutate(
    Case = factor(
      Case,
      levels = c("Homogeneous", "Covariate Shift", "Outcome Shift",
                 "Censoring Shift", "All Shift")
    )
  )

p.rrmse <- plot_rrmse_heatmap(df_heatmap)

pdf(file = "rrmse_l.pdf", width = 8.5, height = 4.5)
p.rrmse
dev.off()


## Larger sample setting 2
load("Res_homo_l2.Rdata")
load("Res_diffX_l2.Rdata")
load("Res_diffT_l2.Rdata")
load("Res_diffC_l2.Rdata")
load("Res_diffAll_l2.Rdata")
load("Res_CLCOX_l2.Rdata")

df_heatmap <- bind_rows(
  make_rrmse_df(result.homo, result.CLCOX.homo, "Homogeneous", M, t),
  make_rrmse_df(result.diffX, result.CLCOX.diffX, "Covariate Shift", M, t),
  make_rrmse_df(result.diffT, result.CLCOX.diffT, "Outcome Shift", M, t),
  make_rrmse_df(result.diffC, result.CLCOX.diffC, "Censoring Shift", M, t),
  make_rrmse_df(result.diffAll, result.CLCOX.diffAll, "All Shift", M, t)
) %>%
  mutate(
    Case = factor(
      Case,
      levels = c("Homogeneous", "Covariate Shift", "Outcome Shift",
                 "Censoring Shift", "All Shift")
    )
  )

p.rrmse <- plot_rrmse_heatmap(df_heatmap)

pdf(file = "rrmse_l2.pdf", width = 8.5, height = 4.5)
p.rrmse
dev.off()


## Main text, time = 30
df.main <- df_heatmap %>% filter(time == 30)

p.rrmse.main <- ggplot(df.main, aes(x = Case, y = Method, fill = RRMSE)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(RRMSE, 2)), size = 3, color = "black") +
  scale_fill_gradientn(
    colors = c("white", "white", "orange"),
    values = scales::rescale(c(0, 1, 100)),
    name = "RRMSE"
  ) +
  facet_grid(. ~ Treatment) +
  labs(x = "", y = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

pdf(file = "rrmse_main.pdf", width = 5.8, height = 2.8)
p.rrmse.main
dev.off()


## Limited overlap setting
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

df_heatmap <- bind_rows(
  make_rrmse_df(result.homo, result.CLCOX.homo, "Homogeneous", M, t),
  make_rrmse_df(result.diffX, result.CLCOX.diffX, "Covariate Shift", M, t),
  make_rrmse_df(result.diffT, result.CLCOX.diffT, "Outcome Shift", M, t),
  make_rrmse_df(result.diffC, result.CLCOX.diffC, "Censoring Shift", M, t),
  make_rrmse_df(result.diffAll, result.CLCOX.diffAll, "All Shift", M, t)
) %>%
  mutate(
    Case = factor(
      Case,
      levels = c("Homogeneous", "Covariate Shift", "Outcome Shift",
                 "Censoring Shift", "All Shift")
    )
  )

p.rrmse <- plot_rrmse_heatmap(df_heatmap)

pdf(file = "rrmse_limO.pdf", width = 8.5, height = 4.5)
p.rrmse
dev.off()
x