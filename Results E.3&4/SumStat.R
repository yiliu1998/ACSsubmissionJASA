### ----------------------------------------------------------------------------
### Reproducibility: This script plots Figures E.2--E.13 in the Online
### Supplemental Material. It reproduces the simulation summary figures for
### treatment-specific survival curves under the main settings and the
### additional poor-overlap setting.
###
### Required input files: truth.Rdata, Res_*_s.Rdata, Res_*_l.Rdata,
### Res_*_l2.Rdata, Res_*_limO.Rdata, and Res_CLCOX_*.Rdata.
### Output files are the PDF files used in Figures E.2--E.13.
### ----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(RColorBrewer)

M <- 500L
t <- c(30, 60, 90)
quant <- 1.96
methods <- c("TGT", "IVW", "POOL", "CCOD", "FED", "CLCOX")
method_levels <- c("FED", "TGT", "IVW", "POOL", "CCOD", "CLCOX")
case_files <- c("Homogeneous" = "homo", "Covariate Shift" = "diffX", "Outcome Shift" = "diffT", "Censoring Shift" = "diffC", "All Shift" = "diffAll")
case_levels <- names(case_files)
plot_colors <- c(brewer.pal(11, "Paired")[c(1, 3, 7, 8, 10)], "plum3")

load("truth.Rdata")

load_result_set <- function(suffix) {
  out <- list()
  for (case_key in case_files) {
    env <- new.env(); load(paste0("Res_", case_key, "_", suffix, ".Rdata"), envir = env)
    out[[case_key]] <- env$results
  }
  env <- new.env(); load(paste0("Res_CLCOX_", suffix, ".Rdata"), envir = env)
  list(results = out, clcox = list(homo = env$result.CLCOX.homo, diffX = env$result.CLCOX.diffX, diffT = env$result.CLCOX.diffT, diffC = env$result.CLCOX.diffC, diffAll = env$result.CLCOX.diffAll))
}

get_df <- function(datlist, coxlist, j, method) {
  if (method == "CLCOX") {
    time_use <- datlist[[j]]$df.TGT$time
    df <- coxlist[[j]]$df.CLCOX
    df[match(time_use, df$time), ]
  } else {
    datlist[[j]][[paste0("df.", method)]]
  }
}

make_long_estimates <- function(datlist, coxlist, case_name, M = 500L) {
  bind_rows(lapply(methods, function(method) {
    bind_rows(lapply(seq_len(M), function(j) {
      df <- get_df(datlist, coxlist, j, method)
      tibble(
        Case = case_name,
        Method = method,
        time = df$time,
        est1 = df$surv1,
        sd1 = df$surv1.sd,
        est0 = df$surv0,
        sd0 = df$surv0.sd
      )
    }))
  }))
}

make_suffix_df <- function(suffix, M = 500L) {
  obj <- load_result_set(suffix)
  bind_rows(lapply(names(case_files), function(case_name) {
    key <- unname(case_files[[case_name]])
    make_long_estimates(obj$results[[key]], obj$clcox[[key]], case_name, M)
  })) %>%
    mutate(
      Case = factor(Case, levels = case_levels),
      Method = factor(Method, levels = method_levels)
    )
}

make_bias_df <- function(df) {
  bind_rows(
    df %>% mutate(Bias = est1 - rep(S1.true, length.out = n()), Treatment = "Treated") %>% select(Bias, Case, Method, time, Treatment),
    df %>% mutate(Bias = est0 - rep(S0.true, length.out = n()), Treatment = "Control") %>% select(Bias, Case, Method, time, Treatment)
  ) %>% mutate(Treatment = factor(Treatment, levels = c("Treated", "Control")))
}

make_summary_df <- function(df, metric = c("relbias", "rrmse", "cp", "ciw")) {
  metric <- match.arg(metric)
  long <- bind_rows(
    df %>% mutate(Treatment = "Treated", est = est1, sd = sd1, truth = rep(S1.true, length.out = n())),
    df %>% mutate(Treatment = "Control", est = est0, sd = sd0, truth = rep(S0.true, length.out = n()))
  ) %>% mutate(Treatment = factor(Treatment, levels = c("Treated", "Control")))
  
  if (metric == "relbias") {
    long %>% group_by(Case, Method, Treatment, time) %>% summarize(Value = abs(mean((est - truth) / truth * 100, na.rm = TRUE)), .groups = "drop")
  } else if (metric == "rrmse") {
    mse <- long %>% group_by(Case, Method, Treatment, time) %>% summarize(MSE = mean((est - truth)^2, na.rm = TRUE), .groups = "drop")
    tgt <- mse %>% filter(Method == "TGT") %>% select(Case, Treatment, time, MSE.TGT = MSE)
    mse %>% left_join(tgt, by = c("Case", "Treatment", "time")) %>% mutate(Value = MSE / MSE.TGT) %>% select(Case, Method, Treatment, time, Value)
  } else if (metric == "cp") {
    long %>% group_by(Case, Method, Treatment, time) %>% summarize(Value = mean((truth < est + quant * sd) & (truth > est - quant * sd), na.rm = TRUE) * 100, .groups = "drop")
  } else {
    long %>% group_by(Case, Method, Treatment, time) %>% summarize(Value = median(2 * quant * sd, na.rm = TRUE), .groups = "drop")
  }
}

plot_bias <- function(df, file) {
  p <- ggplot(make_bias_df(df), aes(x = factor(time), y = Bias, fill = Method)) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.6, color = "gray50") +
    scale_fill_manual(values = plot_colors) +
    geom_hline(yintercept = 0, col = "hotpink2", lty = 2, linewidth = 0.6) +
    facet_grid(Treatment ~ Case) + labs(x = "Time (day)", y = "Bias") + theme_bw()
  pdf(file = file, width = 9, height = 4); print(p); dev.off(); invisible(p)
}

plot_heatmap <- function(df, metric, file) {
  d <- make_summary_df(df, metric)
  if (metric == "relbias") {
    fill <- scale_fill_gradientn(colors = c("white", "white", "green4"), values = scales::rescale(c(0, 0.05, 1)), name = "ARBias%")
    legend_lab <- NULL
  } else if (metric == "rrmse") {
    fill <- scale_fill_gradientn(colors = c("white", "white", "orange"), values = scales::rescale(c(0, 1, 100)), name = "RRMSE")
    legend_lab <- NULL
  } else if (metric == "cp") {
    fill <- scale_fill_gradientn(colors = c("mediumvioletred", "white", "mediumvioletred"), values = scales::rescale(c(0, 80, 95, 99, 100)), name = "CP%", limits = c(0, 100))
    legend_lab <- NULL
  } else {
    fill <- scale_fill_gradientn(colors = c("white", "slateblue2"), name = "CI width")
    legend_lab <- NULL
  }
  p <- ggplot(d, aes(x = factor(time), y = Method, fill = Value)) +
    geom_tile(color = "white") + geom_text(aes(label = round(Value, ifelse(metric == "cp", 1, 2))), size = 3, color = "black") +
    fill + facet_grid(Treatment ~ Case) + labs(x = "Time (day)", y = legend_lab) +
    theme_minimal(base_size = 12) + theme(panel.grid = element_blank())
  pdf(file = file, width = ifelse(metric == "relbias", 9.5, 8.5), height = 4.5); print(p); dev.off(); invisible(p)
}

run_one_setting <- function(suffix, tag) {
  message("Generating figures for setting: ", suffix)
  df <- make_suffix_df(suffix, M)
  plot_bias(df, paste0("bias_", tag, ".pdf"))
  plot_heatmap(df, "relbias", paste0("relbias_", tag, ".pdf"))
  plot_heatmap(df, "rrmse", paste0("rrmse_", tag, ".pdf"))
  plot_heatmap(df, "cp", paste0("cp_", tag, ".pdf"))
  plot_heatmap(df, "ciw", paste0("ciw_", tag, ".pdf"))
}

### Figure mapping:
### E.2:  bias_s.pdf
### E.3:  relbias_s.pdf + rrmse_s.pdf
### E.4:  cp_s.pdf + ciw_s.pdf
### E.5:  bias_l.pdf
### E.6:  relbias_l.pdf + rrmse_l.pdf
### E.7:  cp_l.pdf + ciw_l.pdf
### E.8:  bias_l2.pdf
### E.9:  relbias_l2.pdf + rrmse_l2.pdf
### E.10: cp_l2.pdf + ciw_l2.pdf
### E.11: bias_limO.pdf
### E.12: relbias_limO.pdf + rrmse_limO.pdf
### E.13: cp_limO.pdf + ciw_limO.pdf

run_one_setting("s", "s")       # n_k = 300, main setting
run_one_setting("l", "l")       # n_k = 600, main setting
run_one_setting("l2", "l2")     # n_k = 1000, main setting
run_one_setting("limO", "limO") # n_k = 300, poor-overlap setting

message("Done. Figures E.2--E.13 have been generated as PDF files.")
