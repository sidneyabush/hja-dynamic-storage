# =============================================================================
# Multiple Linear Regression: Catchment Characteristics → Storage Metrics
# =============================================================================
# Purpose: Use stepwise regression to identify which catchment attributes
#          best predict hydrometric storage metrics
#
# Workflow:
#   1. Load site-averaged storage metrics and catchment attributes
#   2. For each storage metric:
#      a. Build full model with all catchment predictors
#      b. Perform backward stepwise AIC selection
#      c. Check VIF for multicollinearity (retain if VIF ≤ 10)
#      d. Extract standardized beta coefficients and p-values
#      e. Report model R² and adjusted R²
#   3. Summarize results across all storage metrics
#
# Inputs:
#   - HJA_Ave_StorageMetrics_CatCharacter.csv (from Correlations_Metrics.R)
#
# Outputs:
#   - MLR_Storage_Catchment_Results.csv: Beta coefficients, p-values, VIF, R²
#   - QA_MLR_model_fits.png: Predicted vs observed for each metric
#
# Author: Based on Pamela Sullivan/Keira Johnson code, adapted by Sidney Bush
# Date: 2026-01-23
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(MASS)         # for stepAIC()
library(car)          # for vif()
library(tidyr)
library(patchwork)

theme_set(theme_classic(base_size = 12))

# Clear environment
rm(list = ls())

# =============================================================================
# 1. SETUP: Directories
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir  <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 2. LOAD SITE-AVERAGED DATA
# =============================================================================

HJA_Ave <- read_csv(
  file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
  show_col_types = FALSE
)

# =============================================================================
# 3. DEFINE OUTCOME VARIABLES (STORAGE METRICS)
# =============================================================================

outcome_vars <- c(
  "recession_curve_slope_mean",
  "RBI_mean",
  "Q5norm_mean",
  "CV_Q5norm_mean",
  "mean_bf_mean",
  "fdc_slope_mean",
  "S_annual_mm_mean",
  "DR_Overall"
)

# =============================================================================
# 4. DEFINE PREDICTOR VARIABLES (CATCHMENT ATTRIBUTES)
# =============================================================================

predictor_vars <- c(
  "Area_km2",
  "Elevation_mean_m",
  "Slope_mean",
  "Harvest",
  "Landslide_Young",
  "Landslide_Total",
  "Lava1_per",
  "Lava2_per",
  "Ash_Per",
  "Pyro_per"
)

# =============================================================================
# 5. RUN STEPWISE REGRESSION FOR EACH STORAGE METRIC
# =============================================================================

results_list <- list()

for (outcome in outcome_vars) {

  # Create formula: outcome ~ all predictors
  formula_full <- as.formula(paste(outcome, "~", paste(predictor_vars, collapse = " + ")))

  # Fit full model
  data_for_model <- HJA_Ave %>%
    select(all_of(c(outcome, predictor_vars))) %>%
    na.omit()

  if (nrow(data_for_model) < 5) {
    next  # Skip if too few observations
  }

  lm_full <- lm(formula_full, data = data_for_model)

  # Stepwise AIC (backward selection)
  lm_aic <- stepAIC(lm_full, direction = "backward", trace = 0)

  # Extract retained variables (excluding intercept)
  retained_vars <- names(coef(lm_aic))[-1]

  if (length(retained_vars) == 0) {
    next  # Skip if no variables retained
  }

  # Check VIF (only if more than 1 predictor)
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

  # Remove intercept from coefficient table
  coef_df <- coef_df %>% filter(variable != "(Intercept)")

  # Compute standardized beta coefficients
  # Scale predictors and outcome, then refit
  data_scaled <- data_for_model %>%
    mutate(across(all_of(c(outcome, retained_vars)), scale))

  formula_scaled <- as.formula(paste(outcome, "~", paste(retained_vars, collapse = " + ")))
  lm_scaled <- lm(formula_scaled, data = data_scaled)
  beta_df <- data.frame(
    variable = names(coef(lm_scaled))[-1],
    beta_std = as.numeric(coef(lm_scaled)[-1])
  )

  # Merge coefficient table with VIF and beta
  result_df <- coef_df %>%
    left_join(vif_df, by = "variable") %>%
    left_join(beta_df, by = "variable") %>%
    mutate(
      outcome = outcome,
      R2 = lm_summary$r.squared,
      R2_adj = lm_summary$adj.r.squared
    ) %>%
    select(outcome, variable, beta_std, `Pr(>|t|)`, VIF, R2, R2_adj)

  colnames(result_df) <- c("Outcome", "Predictor", "Beta_Std", "p_value", "VIF", "R2", "R2_adj")

  results_list[[outcome]] <- result_df
}

# =============================================================================
# 6. COMBINE RESULTS
# =============================================================================

results_combined <- bind_rows(results_list)

# Save results
write.csv(results_combined,
          file.path(output_dir, "MLR_Storage_Catchment_Results.csv"),
          row.names = FALSE)

# =============================================================================
# 7. QA PLOTS: PREDICTED VS OBSERVED
# =============================================================================

plot_list <- list()

for (outcome in outcome_vars) {

  formula_full <- as.formula(paste(outcome, "~", paste(predictor_vars, collapse = " + ")))

  data_for_model <- HJA_Ave %>%
    select(all_of(c(outcome, predictor_vars, "site"))) %>%
    na.omit()

  if (nrow(data_for_model) < 5) {
    next
  }

  lm_full <- lm(formula_full, data = data_for_model)
  lm_aic <- stepAIC(lm_full, direction = "backward", trace = 0)

  data_for_model$predicted <- predict(lm_aic, data_for_model)

  r2_adj <- round(summary(lm_aic)$adj.r.squared, 2)

  p <- ggplot(data_for_model, aes(x = .data[[outcome]], y = predicted)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
      title = outcome,
      subtitle = paste0("R² adj = ", r2_adj),
      x = "Observed",
      y = "Predicted"
    ) +
    theme_classic(base_size = 10)

  plot_list[[outcome]] <- p
}

# Combine plots
if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = 3)

  ggsave(
    file.path(output_dir, "QA_MLR_model_fits.png"),
    combined_plot, width = 12, height = 10, dpi = 300
  )
}
