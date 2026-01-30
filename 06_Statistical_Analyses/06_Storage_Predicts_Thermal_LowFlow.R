# =============================================================================
# Storage Metrics Predict Thermal & Low-Flow Responses
# =============================================================================
# Purpose: Test whether hydrometric storage metrics can predict ecological
#          thermal and low-flow responses
#
# Hypotheses:
#   H1: Greater storage capacity → lower maximum stream temperatures
#   H2: Greater storage capacity → higher minimum discharge
#   H3: Greater storage capacity → lower temperature at minimum discharge
#
# Workflow:
#   1. Load annual storage metrics
#   2. Load thermal/low-flow metrics (from Stream_Temperature_LowFlow_Metrics.R)
#   3. Merge datasets by site and water year
#   4. For each response variable:
#      a. Explore bivariate relationships with storage metrics
#      b. Build MLR models with stepwise selection
#      c. Test model assumptions
#      d. Visualize relationships
#
# Inputs:
#   - HJA_StorageMetrics_Annual_All.csv: Annual storage metrics
#   - stream_thermal_lowflow_metrics_annual.csv: Response variables
#
# Outputs:
#   - Storage_Thermal_LowFlow_Models.csv: Model results
#   - QA_Storage_Thermal_Scatterplots.png: Bivariate relationships
#   - QA_Model_Diagnostics.png: Residual plots
#
# Author: Sidney Bush
# Date: 2026-01-23
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(MASS)         # for stepAIC()
library(car)          # for vif()
library(patchwork)
library(GGally)       # for ggpairs()

theme_set(theme_classic(base_size = 12))

# Clear environment
rm(list = ls())

# Source configuration (paths, site definitions, water year range)
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir) || script_dir == "") script_dir <- getwd()
config_path <- file.path(dirname(script_dir), "config.R")
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found. Please ensure config.R exists in the repo root.")
}

# Use configuration values
site_order <- SITE_ORDER_HYDROMETRIC
base_dir   <- BASE_DATA_DIR
output_dir <- OUTPUT_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 2. LOAD STORAGE & THERMAL METRICS (ANNUAL)
# =============================================================================
# This file already contains storage metrics AND thermal/low-flow metrics
# (merged by 06_Aggregate_All_Metrics.R)

merged_data <- read_csv(
  file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv"),
  show_col_types = FALSE
) %>%
  filter(!site %in% c("GSLOOK_FULL", "GSWSMA", "GSWSMF", "GSMACK")) %>%  # Exclude non-analysis sites
  filter(!is.na(max_temp_7d_C))  # Keep only years with thermal data

# =============================================================================
# 5. DEFINE PREDICTOR AND RESPONSE VARIABLES
# =============================================================================

# Storage metrics (predictors)
# Note: mean_bf excluded - no data overlap with thermal metrics
storage_predictors <- c(
  "recession_curve_slope",
  "RBI",
  "Q5norm",
  "CV_Q5norm",
  "fdc_slope",
  "S_annual_mm",
  "DS_sum"
)

# Thermal/low-flow responses
response_vars <- c(
  "max_temp_7d_C",           # Maximum 7-day average stream temperature
  "min_Q_7d_mm_d",           # Minimum 7-day average discharge
  "temp_at_min_Q_7d_C"       # Temperature at time of minimum discharge
)

# =============================================================================
# 6. EXPLORATORY CORRELATIONS
# =============================================================================

# Calculate correlation matrix
cor_data <- merged_data %>%
  dplyr::select(all_of(c(response_vars, storage_predictors)))

cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

# Extract correlations between responses and storage metrics
cor_response_storage <- cor_matrix[response_vars, storage_predictors]

write.csv(cor_response_storage,
          file.path(output_dir, "Correlations_Storage_Thermal_LowFlow.csv"))

# =============================================================================
# 7. BUILD MLR MODELS FOR EACH RESPONSE
# =============================================================================

model_results <- list()

for (response in response_vars) {

  # Create formula: response ~ all storage predictors
  formula_full <- as.formula(paste(response, "~", paste(storage_predictors, collapse = " + ")))

  # Fit full model
  data_for_model <- merged_data %>%
    dplyr::select(all_of(c(response, storage_predictors, "site", "year"))) %>%
    na.omit()

  if (nrow(data_for_model) < 20) {
    next  # Skip if too few observations
  }

  lm_full <- lm(formula_full, data = data_for_model)

  # Stepwise AIC (backward selection)
  lm_aic <- stepAIC(lm_full, direction = "backward", trace = 0)

  # Extract retained variables
  retained_vars <- names(coef(lm_aic))[-1]

  if (length(retained_vars) == 0) {
    next
  }

  # Check VIF
  if (length(retained_vars) > 1) {
    vif_vals <- vif(lm_aic)
    vif_df <- data.frame(
      variable = names(vif_vals),
      VIF = as.numeric(vif_vals)
    )
  } else {
    vif_df <- data.frame(
      variable = retained_vars,
      VIF = NA
    )
  }

  # Extract model summary
  lm_summary <- summary(lm_aic)
  coef_df <- as.data.frame(lm_summary$coefficients)
  coef_df$variable <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df <- coef_df %>% filter(variable != "(Intercept)")

  # Compute standardized beta coefficients
  data_scaled <- data_for_model %>%
    mutate(across(all_of(c(response, retained_vars)), scale))

  formula_scaled <- as.formula(paste(response, "~", paste(retained_vars, collapse = " + ")))
  lm_scaled <- lm(formula_scaled, data = data_scaled)
  beta_df <- data.frame(
    variable = names(coef(lm_scaled))[-1],
    beta_std = as.numeric(coef(lm_scaled)[-1])
  )

  # Merge results
  result_df <- coef_df %>%
    left_join(vif_df, by = "variable") %>%
    left_join(beta_df, by = "variable") %>%
    mutate(
      response = response,
      R2 = lm_summary$r.squared,
      R2_adj = lm_summary$adj.r.squared,
      n_obs = nrow(data_for_model)
    ) %>%
    dplyr::select(response, variable, beta_std, `Pr(>|t|)`, VIF, R2, R2_adj, n_obs)

  colnames(result_df) <- c("Response", "Predictor", "Beta_Std", "p_value", "VIF", "R2", "R2_adj", "n")

  model_results[[response]] <- result_df

  # Store fitted model for diagnostics
  assign(paste0("model_", response), lm_aic)
}

# Combine all results
model_results_combined <- bind_rows(model_results)

write.csv(model_results_combined,
          file.path(output_dir, "Storage_Thermal_LowFlow_Models.csv"),
          row.names = FALSE)

# =============================================================================
# 8. VISUALIZE: BIVARIATE RELATIONSHIPS
# =============================================================================

# Create scatterplots for each response vs top storage predictors
plot_list <- list()

for (response in response_vars) {

  # Skip if no model was built for this response
  if (is.null(model_results[[response]])) next

  # Get top 3 predictors for this response (by absolute beta)
  top_predictors <- model_results[[response]] %>%
    arrange(desc(abs(Beta_Std))) %>%
    slice(1:min(3, n())) %>%
    pull(Predictor)

  if (length(top_predictors) == 0) next

  for (pred in top_predictors) {

    p <- ggplot(merged_data, aes(x = .data[[pred]], y = .data[[response]])) +
      geom_point(alpha = 0.5, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
      labs(
        title = paste(response, "vs", pred),
        x = pred,
        y = response
      ) +
      theme_classic(base_size = 10)

    plot_list[[paste(response, pred, sep = "_")]] <- p
  }
}

if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = 3)

  ggsave(
    file.path(output_dir, "QA_Storage_Thermal_Scatterplots.png"),
    combined_plot, width = 14, height = 12, dpi = 300
  )
}

# =============================================================================
# 9. MODEL DIAGNOSTICS: RESIDUAL PLOTS
# =============================================================================

diag_plots <- list()

for (response in response_vars) {

  model_name <- paste0("model_", response)

  if (!exists(model_name)) next

  model <- get(model_name)

  # Residuals vs Fitted
  data_for_model <- merged_data %>%
    dplyr::select(all_of(c(response, storage_predictors))) %>%
    na.omit()

  data_for_model$fitted <- fitted(model)
  data_for_model$residuals <- residuals(model)

  p1 <- ggplot(data_for_model, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Residuals vs Fitted:", response),
         x = "Fitted values", y = "Residuals") +
    theme_classic(base_size = 10)

  # Q-Q plot
  p2 <- ggplot(data_for_model, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = paste("Q-Q Plot:", response)) +
    theme_classic(base_size = 10)

  diag_plots[[response]] <- p1 + p2
}

if (length(diag_plots) > 0) {
  combined_diag <- wrap_plots(diag_plots, ncol = 1)

  ggsave(
    file.path(output_dir, "QA_Model_Diagnostics.png"),
    combined_diag, width = 12, height = 10, dpi = 300
  )
}
