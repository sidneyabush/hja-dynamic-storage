# -----------------------------------------------------------------------------
# Plots: Thermal and Low-Flow Prediction Models
# -----------------------------------------------------------------------------
# Inputs:
#   - Storage_Thermal_LowFlow_ModelSummary.csv
#   - HJA_StorageMetrics_Annual_All.csv
#
# Outputs:
#   - QA_Storage_Thermal_Scatterplots.png
#   - QA_Model_Diagnostics.png
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

rm(list = ls())

script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    getwd()
  }
})
if (is.null(script_dir) || script_dir == "" || script_dir == ".") {
  script_dir <- getwd()
}

config_path <- file.path(script_dir, "config.R")
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(script_dir), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- file.path(getwd(), "config.R")
}
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found. Please ensure config.R exists in the repo root.")
}

output_dir <- OUTPUT_DIR
plot_dir <- file.path(FIGURES_DIR, "Statistical")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

summary_file <- file.path(output_dir, "Storage_Thermal_LowFlow_ModelSummary.csv")
if (!file.exists(summary_file)) {
  stop("Missing file: Storage_Thermal_LowFlow_ModelSummary.csv")
}

model_summary <- read_csv(summary_file, show_col_types = FALSE)
annual_file <- file.path(output_dir, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(BASE_DATA_DIR, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
}

analysis_data <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

if (!("temp_during_min_Q_7d_C" %in% names(analysis_data)) && ("temp_at_min_Q_7d_C" %in% names(analysis_data))) {
  analysis_data <- analysis_data %>%
    mutate(temp_during_min_Q_7d_C = temp_at_min_Q_7d_C)
}

fit_from_summary <- function(df, response, predictor_string) {
  preds <- trimws(unlist(strsplit(predictor_string, ";")))
  preds <- preds[preds != "" & preds %in% names(df)]
  if (length(preds) == 0 || !(response %in% names(df))) {
    return(NULL)
  }

  fit_df <- df %>%
    select(all_of(c(response, preds))) %>%
    na.omit()

  if (nrow(fit_df) < 10) {
    return(NULL)
  }

  fm <- as.formula(paste(response, "~", paste(preds, collapse = " + ")))
  lm(fm, data = fit_df)
}

scatter_plots <- list()
diag_panels <- list()

for (i in seq_len(nrow(model_summary))) {
  response <- model_summary$Response[i]
  predictor_string <- model_summary$Predictors_Final[i]

  lm_fit <- fit_from_summary(analysis_data, response, predictor_string)
  if (is.null(lm_fit)) {
    next
  }

  preds <- attr(terms(lm_fit), "term.labels")

  top_preds <- preds[1:min(3, length(preds))]
  for (pred in top_preds) {
    p_scatter <- ggplot(analysis_data, aes(x = .data[[pred]], y = .data[[response]])) +
      geom_point(alpha = 0.5, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
      labs(title = paste(response, "vs", pred), x = pred, y = response) +
      theme_classic(base_size = 10)

    scatter_plots[[paste(response, pred, sep = "__")]] <- p_scatter
  }

  fit_df <- model.frame(lm_fit)
  fit_df$fitted <- fitted(lm_fit)
  fit_df$residuals <- residuals(lm_fit)

  p_resid <- ggplot(fit_df, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Residuals vs Fitted:", response), x = "Fitted", y = "Residuals") +
    theme_classic(base_size = 10)

  p_qq <- ggplot(fit_df, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = paste("Q-Q Plot:", response)) +
    theme_classic(base_size = 10)

  diag_panels[[response]] <- p_resid + p_qq
}

if (length(scatter_plots) > 0) {
  ggsave(
    file.path(plot_dir, "QA_Storage_Thermal_Scatterplots.png"),
    wrap_plots(scatter_plots, ncol = 3),
    width = 14,
    height = 12,
    dpi = 300
  )
}

if (length(diag_panels) > 0) {
  ggsave(
    file.path(plot_dir, "QA_Model_Diagnostics.png"),
    wrap_plots(diag_panels, ncol = 1),
    width = 12,
    height = 10,
    dpi = 300
  )
}
