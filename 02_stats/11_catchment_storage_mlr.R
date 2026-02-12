# -----------------------------------------------------------------------------
# Multiple Linear Regression: Catchment Characteristics -> Storage Metrics
# -----------------------------------------------------------------------------
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
#   - master_site.csv (from 03b_aggregate_metrics.R)
#
# Outputs:
#   - MLR_Storage_Catchment_Results.csv: Beta coefficients, p-values, VIF, R²
#   - MLR_Storage_Catchment_ModelSummary.csv: Final predictors and model fit stats
#   - MLR_Storage_Catchment_BetaMatrix.csv: Legacy-style beta matrix for plotting
#   - QA_MLR_model_fits.png: Predicted vs observed for each metric
#
# Author: Based on Pamela Sullivan/Keira Johnson code, adapted by Sidney Bush
# Date: 2026-01-23
# -----------------------------------------------------------------------------

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

# Source configuration (paths, site definitions, water year range)
# Get script directory (works with source() and Rscript)
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

# Use configuration values
site_order <- SITE_ORDER_HYDROMETRIC
base_dir   <- BASE_DATA_DIR
output_dir <- OUTPUT_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Modeling switches
# Keep defaults aligned with legacy behavior until you decide otherwise.
MODEL_USE_SCALED_PREDICTORS <- FALSE
MODEL_USE_ITERATIVE_VIF <- FALSE
VIF_THRESHOLD <- 10

# -----------------------------------------------------------------------------
# 2. LOAD SITE-AVERAGED DATA
# -----------------------------------------------------------------------------

site_file <- file.path(output_dir, MASTER_SITE_FILE)
if (!file.exists(site_file)) {
  site_file <- file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv")
}

HJA_Ave <- read_csv(
  site_file,
  show_col_types = FALSE
) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

# -----------------------------------------------------------------------------
# 3. DEFINE OUTCOME VARIABLES (STORAGE METRICS)
# -----------------------------------------------------------------------------
# NOTE: Q5norm and CV_Q5norm are NOT storage metrics - they are response variables
# Storage metrics by type (using method abbreviations):
#   Dynamic: RBI, RCS, FDC, SD
#   Mobile: CHS, MTT, Fyw, DR
#   Extended Dynamic: WB

outcome_vars <- c(
  "RCS_mean",   # RCS - Recession Curve Slope - Dynamic
  "RBI_mean",   # RBI - Richards-Baker Index - Dynamic
  "FDC_mean",   # FDC - Flow Duration Curve - Dynamic
  "SD_mean",    # SD  - Storage-Discharge - Dynamic
  "CHS_mean",   # CHS - Chemical Hydrograph Separation - Mobile
  "MTT",        # MTT - Mean Transit Time - Mobile
  "Fyw",        # Fyw - Young Water Fraction - Mobile
  "DR",         # DR  - Damping Ratio - Mobile
  "WB_mean"     # WB  - Water Balance - Extended Dynamic
)

# -----------------------------------------------------------------------------
# 4. DEFINE PREDICTOR VARIABLES (CATCHMENT ATTRIBUTES)
# -----------------------------------------------------------------------------
# Note: Reduced predictor set to avoid overfitting (N=10 sites)
# Selected key predictors based on theoretical importance

predictor_vars <- c(
  "Slope_mean",
  "Harvest",
  "Landslide_Total",
  "Landslide_Young",
  "Lava1_per",
  "Lava2_per",
  "Ash_Per",
  "Pyro_per"
)

# -----------------------------------------------------------------------------
# 5. RUN STEPWISE REGRESSION FOR EACH STORAGE METRIC
# -----------------------------------------------------------------------------

results_list <- list()
model_summary_list <- list()
model_objects <- list()

fit_mlr_with_vif <- function(data_in, outcome, predictors,
                             use_scaled_predictors = FALSE,
                             use_iterative_vif = FALSE,
                             vif_threshold = 10) {
  predictor_keep <- predictors[predictors %in% names(data_in)]
  if (!(outcome %in% names(data_in)) || length(predictor_keep) == 0) {
    return(NULL)
  }

  model_df <- data_in %>%
    dplyr::select(all_of(c(outcome, predictor_keep))) %>%
    na.omit()

  if (nrow(model_df) < 5) {
    return(NULL)
  }

  if (use_scaled_predictors) {
    model_df <- model_df %>%
      mutate(across(all_of(predictor_keep), ~ as.numeric(scale(.x))))
  }

  formula_full <- as.formula(paste(outcome, "~", paste(predictor_keep, collapse = " + ")))
  lm_full <- tryCatch(lm(formula_full, data = model_df), error = function(e) NULL)
  if (is.null(lm_full)) {
    return(NULL)
  }

  lm_aic <- tryCatch(
    stepAIC(lm_full, direction = "backward", trace = 0),
    error = function(e) NULL
  )
  if (is.null(lm_aic)) {
    return(NULL)
  }

  if (use_iterative_vif) {
    # iterative VIF pruning until all retained predictors are <= threshold
    repeat {
      retained_vars <- setdiff(names(coef(lm_aic)), "(Intercept)")
      if (length(retained_vars) <= 1) {
        break
      }

      vif_vals <- tryCatch(vif(lm_aic), error = function(e) NULL)
      if (is.null(vif_vals)) {
        break
      }
      if (max(vif_vals, na.rm = TRUE) <= vif_threshold) {
        break
      }

      drop_var <- names(which.max(vif_vals))
      retained_next <- setdiff(retained_vars, drop_var)
      if (length(retained_next) == 0) {
        return(NULL)
      }

      lm_next <- lm(
        as.formula(paste(outcome, "~", paste(retained_next, collapse = " + "))),
        data = model_df
      )
      lm_aic <- tryCatch(
        stepAIC(lm_next, direction = "backward", trace = 0),
        error = function(e) lm_next
      )
    }
  }

  retained_vars <- setdiff(names(coef(lm_aic)), "(Intercept)")
  if (length(retained_vars) == 0) {
    return(NULL)
  }

  if (length(retained_vars) > 1) {
    vif_vals <- tryCatch(vif(lm_aic), error = function(e) NULL)
    if (is.null(vif_vals)) {
      vif_df <- data.frame(variable = retained_vars, VIF = NA_real_)
    } else {
      vif_df <- data.frame(variable = names(vif_vals), VIF = as.numeric(vif_vals))
    }
  } else {
    vif_df <- data.frame(variable = retained_vars, VIF = NA_real_)
  }

  lm_summary <- summary(lm_aic)
  coef_df <- as.data.frame(lm_summary$coefficients)
  coef_df$variable <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df <- coef_df %>% filter(variable != "(Intercept)")

  beta_df_data <- model_df
  beta_df_data[[outcome]] <- as.numeric(scale(beta_df_data[[outcome]]))
  lm_beta <- lm(
    as.formula(paste(outcome, "~", paste(retained_vars, collapse = " + "))),
    data = beta_df_data
  )
  beta_df <- data.frame(
    variable = names(coef(lm_beta))[-1],
    beta_std = as.numeric(coef(lm_beta)[-1])
  )

  result_df <- coef_df %>%
    left_join(vif_df, by = "variable") %>%
    left_join(beta_df, by = "variable") %>%
    mutate(
      outcome = outcome,
      R2 = lm_summary$r.squared,
      R2_adj = lm_summary$adj.r.squared,
      RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE))
    ) %>%
    dplyr::select(outcome, variable, beta_std, `Pr(>|t|)`, VIF, R2, R2_adj, RMSE)

  names(result_df) <- c("Outcome", "Predictor", "Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE")

  summary_df <- tibble(
    Outcome = outcome,
    N = nrow(model_df),
    Predictors_Final = paste(retained_vars, collapse = "; "),
    R2 = lm_summary$r.squared,
    R2_adj = lm_summary$adj.r.squared,
    RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
    AIC = AIC(lm_aic)
  )

  list(
    coefficients = result_df,
    summary = summary_df,
    model = lm_aic,
    data = model_df
  )
}

for (outcome in outcome_vars) {
  fit_obj <- fit_mlr_with_vif(
    HJA_Ave,
    outcome,
    predictor_vars,
    use_scaled_predictors = MODEL_USE_SCALED_PREDICTORS,
    use_iterative_vif = MODEL_USE_ITERATIVE_VIF,
    vif_threshold = VIF_THRESHOLD
  )
  if (is.null(fit_obj)) {
    next
  }
  results_list[[outcome]] <- fit_obj$coefficients
  model_summary_list[[outcome]] <- fit_obj$summary
  model_objects[[outcome]] <- fit_obj
}

# -----------------------------------------------------------------------------
# 6. COMBINE RESULTS
# -----------------------------------------------------------------------------

results_combined <- bind_rows(results_list)
model_summary <- bind_rows(model_summary_list)
beta_matrix <- results_combined %>%
  dplyr::select(Predictor, Outcome, Beta_Std) %>%
  tidyr::pivot_wider(names_from = Outcome, values_from = Beta_Std)

# Save results
write.csv(results_combined,
          file.path(output_dir, "MLR_Storage_Catchment_Results.csv"),
          row.names = FALSE)
write.csv(model_summary,
          file.path(output_dir, "MLR_Storage_Catchment_ModelSummary.csv"),
          row.names = FALSE)
write.csv(beta_matrix,
          file.path(output_dir, "MLR_Storage_Catchment_BetaMatrix.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. QA PLOTS: PREDICTED VS OBSERVED
# -----------------------------------------------------------------------------

plot_list <- list()

for (outcome in outcome_vars) {
  if (!(outcome %in% names(model_objects))) {
    next
  }

  model_obj <- model_objects[[outcome]]$model
  data_for_model <- model_objects[[outcome]]$data
  data_for_model$predicted <- predict(model_obj, data_for_model)

  r2_adj <- round(summary(model_obj)$adj.r.squared, 2)

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
