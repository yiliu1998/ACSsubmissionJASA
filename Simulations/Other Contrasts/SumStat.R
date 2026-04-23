## Summary statistics and heatmaps
## Metrics: absolute relative bias, relative RMSE, and coverage probability

load("truth.Rdata")

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

M <- 500
methods <- c("TGT", "CCOD", "FED")
method_levels <- c("FED", "TGT", "CCOD")
time_plot <- c(30, 60)

case_files <- c(
  "Homogeneous" = "homo_500.Rdata",
  "Covariate Shift" = "diffX_500.Rdata",
  "Outcome Shift" = "diffT_500.Rdata",
  "Censoring Shift" = "diffC_500.Rdata",
  "All Shifts" = "diffAll_500.Rdata"
)
case_levels <- names(case_files)

load_results <- function(file) {
  env <- new.env()
  load(file, envir = env)
  env$results
}

extract_time_metric <- function(results, metric, method, M = 500) {
  value_col <- metric
  df_name <- paste0("df.", metric, ".", method)
  time_use <- results[[1]][[df_name]]$time
  est <- sd <- matrix(NA_real_, nrow = M, ncol = length(time_use))
  colnames(est) <- colnames(sd) <- time_use
  
  for (j in seq_len(M)) {
    tmp <- results[[j]][[df_name]]
    est[j, ] <- tmp[, value_col]
    sd[j, ] <- tmp[, "sd"]
  }
  
  list(time = time_use, est = est, sd = sd)
}

extract_rmst_metric <- function(results, measure, method, M = 500) {
  df_name <- paste0("df.RMST.", measure, ".", method)
  est <- sd <- rep(NA_real_, M)
  
  for (j in seq_len(M)) {
    tmp <- results[[j]][[df_name]]
    est[j] <- tmp[, "RMST"]
    sd[j] <- tmp[, "sd"]
  }
  
  list(est = est, sd = sd)
}

summarize_time_metric <- function(results, case_name, metric, true_value, M = 500) {
  extracted <- setNames(
    lapply(methods, function(m) extract_time_metric(results, metric, m, M)),
    methods
  )
  time_use <- extracted[["TGT"]]$time
  
  out <- lapply(methods, function(m) {
    est <- extracted[[m]]$est
    sd <- extracted[[m]]$sd
    err <- sweep(est, 2, true_value, "-")
    tgt_err <- sweep(extracted[["TGT"]]$est, 2, true_value, "-")
    
    tibble(
      Case = case_name,
      Metric = metric,
      Method = m,
      time = time_use,
      ARBias = abs(colMeans(sweep(err, 2, true_value, "/") * 100, na.rm = TRUE)),
      RRMSE = sqrt(colMeans(err^2, na.rm = TRUE)) /
        sqrt(colMeans(tgt_err^2, na.rm = TRUE)),
      CP = colMeans(
        sweep(est + 1.96 * sd, 2, true_value, ">") &
          sweep(est - 1.96 * sd, 2, true_value, "<"),
        na.rm = TRUE
      ) * 100
    )
  })
  
  bind_rows(out)
}

summarize_rmst_metric <- function(results, case_name, M = 500) {
  rmst_info <- tibble(
    Measure = c("A = 0", "A = 1", "Difference"),
    suffix = c("0", "1", "diff"),
    true_value = c(as.numeric(RMST.0.true), as.numeric(RMST.1.true), as.numeric(RMST.diff.true))
  )
  
  bind_rows(lapply(seq_len(nrow(rmst_info)), function(i) {
    measure <- rmst_info$Measure[i]
    suffix <- rmst_info$suffix[i]
    true_value <- rmst_info$true_value[i]
    
    extracted <- setNames(
      lapply(methods, function(m) extract_rmst_metric(results, suffix, m, M)),
      methods
    )
    tgt_err <- extracted[["TGT"]]$est - true_value
    
    bind_rows(lapply(methods, function(m) {
      est <- extracted[[m]]$est
      sd <- extracted[[m]]$sd
      err <- est - true_value
      
      tibble(
        Case = case_name,
        Metric = "RMST",
        Method = m,
        Measure = measure,
        ARBias = abs(mean(err / true_value * 100, na.rm = TRUE)),
        RRMSE = sqrt(mean(err^2, na.rm = TRUE)) / sqrt(mean(tgt_err^2, na.rm = TRUE)),
        CP = mean(
          (true_value < est + 1.96 * sd) &
            (true_value > est - 1.96 * sd),
          na.rm = TRUE
        ) * 100
      )
    }))
  }))
}

summarize_case <- function(case_name, file, M = 500) {
  results <- load_results(file)
  bind_rows(
    summarize_time_metric(results, case_name, "RD", RD.true, M),
    summarize_time_metric(results, case_name, "SR", SR.true, M),
    summarize_rmst_metric(results, case_name, M)
  )
}

sumstat <- bind_rows(
  Map(function(case_name, file) summarize_case(case_name, file, M),
      names(case_files), case_files)
) %>%
  mutate(
    Case = factor(Case, levels = case_levels),
    Method = factor(Method, levels = method_levels),
    Measure = factor(Measure, levels = c("A = 0", "A = 1", "Difference"))
  )

plot_time_heatmap <- function(data, value_col, legend_name, fill_colors, fill_values,
                              file, width = 7, height = 1.8, limits = NULL) {
  p <- ggplot(data, aes(x = factor(time), y = Method, fill = .data[[value_col]])) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(.data[[value_col]], 2)), size = 3, color = "black") +
    scale_fill_gradientn(
      colors = fill_colors,
      values = scales::rescale(fill_values),
      name = legend_name,
      limits = limits
    ) +
    facet_grid(~ Case) +
    labs(x = "Time (day)", y = "") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
  
  pdf(file = file, width = width, height = height)
  print(p)
  dev.off()
  invisible(p)
}

plot_rmst_heatmap <- function(data, value_col, legend_name, fill_colors, fill_values,
                              file, width = 8, height = 2, limits = NULL,
                              show_x_text = FALSE) {
  p <- ggplot(data, aes(x = Measure, y = Method, fill = .data[[value_col]])) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(.data[[value_col]], 2)), size = 3, color = "black") +
    scale_fill_gradientn(
      colors = fill_colors,
      values = scales::rescale(fill_values),
      name = legend_name,
      limits = limits
    ) +
    facet_grid(~ Case) +
    labs(x = "", y = "") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
  
  if (!show_x_text) {
    p <- p + theme(axis.text.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  pdf(file = file, width = width, height = height)
  print(p)
  dev.off()
  invisible(p)
}

## Bias plots
bias_colors <- c("white", "white", "green4")
bias_values <- c(0, 5, 100)

plot_time_heatmap(
  sumstat %>% filter(Metric == "RD", time %in% time_plot),
  value_col = "ARBias",
  legend_name = "ARBias%",
  fill_colors = bias_colors,
  fill_values = bias_values,
  file = "bias_RD.pdf"
)

plot_time_heatmap(
  sumstat %>% filter(Metric == "SR", time %in% time_plot),
  value_col = "ARBias",
  legend_name = "ARBias%",
  fill_colors = bias_colors,
  fill_values = bias_values,
  file = "bias_SR.pdf"
)

plot_rmst_heatmap(
  sumstat %>% filter(Metric == "RMST"),
  value_col = "ARBias",
  legend_name = "ARBias%",
  fill_colors = bias_colors,
  fill_values = bias_values,
  file = "bias_RMST.pdf",
  width = 8,
  height = 2,
  show_x_text = FALSE
)

## RRMSE plots
rrmse_colors <- c("white", "white", "orange")

plot_time_heatmap(
  sumstat %>% filter(Metric == "RD", time %in% time_plot),
  value_col = "RRMSE",
  legend_name = "RRMSE",
  fill_colors = rrmse_colors,
  fill_values = c(0, 1, 2.5),
  file = "rrmse_RD.pdf"
)

plot_time_heatmap(
  sumstat %>% filter(Metric == "SR", time %in% time_plot),
  value_col = "RRMSE",
  legend_name = "RRMSE",
  fill_colors = rrmse_colors,
  fill_values = c(0, 1, 1.35),
  file = "rrmse_SR.pdf"
)

plot_rmst_heatmap(
  sumstat %>% filter(Metric == "RMST"),
  value_col = "RRMSE",
  legend_name = "RRMSE",
  fill_colors = rrmse_colors,
  fill_values = c(0, 1, 100),
  file = "rrmse_RMST.pdf",
  width = 8,
  height = 2,
  show_x_text = FALSE
)

## Coverage probability plots
cp_colors <- c("mediumvioletred", "white", "mediumvioletred")
cp_values <- c(0, 87, 95, 99, 100)

plot_time_heatmap(
  sumstat %>% filter(Metric == "RD", time %in% time_plot),
  value_col = "CP",
  legend_name = "CP%",
  fill_colors = cp_colors,
  fill_values = cp_values,
  file = "cp_RD.pdf",
  limits = c(0, 100)
)

plot_time_heatmap(
  sumstat %>% filter(Metric == "SR", time %in% time_plot),
  value_col = "CP",
  legend_name = "CP%",
  fill_colors = cp_colors,
  fill_values = cp_values,
  file = "cp_SR.pdf",
  limits = c(0, 100)
)

plot_rmst_heatmap(
  sumstat %>% filter(Metric == "RMST"),
  value_col = "CP",
  legend_name = "CP%",
  fill_colors = cp_colors,
  fill_values = cp_values,
  file = "cp_RMST.pdf",
  width = 8,
  height = 2.5,
  limits = c(0, 100),
  show_x_text = TRUE
)

## Optional: save the numerical summary used by all plots
saveRDS(sumstat, file = "sumstat_corrected.rds")
