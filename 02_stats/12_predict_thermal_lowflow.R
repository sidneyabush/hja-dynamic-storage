# -----------------------------------------------------------------------------
# Predict Thermal and Low-Flow Metrics from Storage and Catchment Predictors
# -----------------------------------------------------------------------------
# Purpose: Fit stepwise-AIC linear models for annual eco response metrics.
#
# Response metrics:
#   1) max_temp_7d_C
#   2) min_Q_7d_mm_d
#   3) temp_during_min_Q_7d_C
#
# Inputs:
#   - master_annual.csv
#
# Outputs:
#   - Correlations_Storage_Thermal_LowFlow.csv
#   - Storage_Thermal_LowFlow_Models.csv
#   - Storage_Thermal_LowFlow_ModelSummary.csv
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(MASS)
library(car)

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
base_dir <- BASE_DATA_DIR
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# model switches
MODEL_USE_SCALED_PREDICTORS <- FALSE
MODEL_MIN_N <- 20
ALLOW_LAVA_AND_ASH_TOGETHER <- FALSE
ALLOW_TOTAL_AND_YOUNG_LS_TOGETHER <- FALSE

annual_file <- file.path(output_dir, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
}

merged_data <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

if (!("temp_during_min_Q_7d_C" %in% names(merged_data)) && ("temp_at_min_Q_7d_C" %in% names(merged_data))) {
  merged_data <- merged_data %>%
    mutate(temp_during_min_Q_7d_C = temp_at_min_Q_7d_C)
}

response_vars <- c("max_temp_7d_C", "min_Q_7d_mm_d", "temp_during_min_Q_7d_C")
response_vars <- response_vars[response_vars %in% names(merged_data)]

storage_predictors_all <- c("RCS", "RBI", "FDC", "SD", "WB", "CHS", "MTT", "Fyw", "DR")
storage_predictors_all <- storage_predictors_all[storage_predictors_all %in% names(merged_data)]

geology_predictors <- c("Lava1_per", "Lava2_per", "Ash_Per")
geology_predictors <- geology_predictors[geology_predictors %in% names(merged_data)]

landslide_predictors <- c("Landslide_Total", "Landslide_Young")
landslide_predictors <- landslide_predictors[landslide_predictors %in% names(merged_data)]

if (length(response_vars) == 0) {
  stop("No response variables found in master annual dataset")
}
if (length(storage_predictors_all) == 0) {
  stop("No storage predictors found in master annual dataset")
}

# Candidate sets follow your current modeling rule:
# - geology terms are not used together unless switched on
# - landslide total and young are not used together unless switched on
if (ALLOW_LAVA_AND_ASH_TOGETHER) {
  geology_options <- list(character(0), geology_predictors)
} else {
  geology_options <- c(list(character(0)), lapply(geology_predictors, function(x) x))
}

if (ALLOW_TOTAL_AND_YOUNG_LS_TOGETHER) {
  ls_options <- list(character(0), landslide_predictors)
} else {
  ls_options <- c(list(character(0)), lapply(landslide_predictors, function(x) x))
}

candidate_predictor_sets <- list()
set_idx <- 1
for (geo_set in geology_options) {
  for (ls_set in ls_options) {
    candidate_predictor_sets[[set_idx]] <- unique(c(storage_predictors_all, geo_set, ls_set))
    set_idx <- set_idx + 1
  }
}

fit_one_model <- function(df_in, response, predictors, use_scaled_predictors = FALSE, min_n = 20) {
  keep_cols <- unique(c(response, predictors, "site", "year"))
  keep_cols <- keep_cols[keep_cols %in% names(df_in)]

  model_df <- df_in %>%
    dplyr::select(all_of(keep_cols)) %>%
    na.omit()

  if (nrow(model_df) < min_n) {
    return(NULL)
  }

  predictors_use <- predictors[predictors %in% names(model_df)]
  if (length(predictors_use) == 0) {
    return(NULL)
  }

  if (use_scaled_predictors) {
    model_df <- model_df %>%
      mutate(across(all_of(predictors_use), ~ as.numeric(scale(.x))))
  }

  formula_full <- as.formula(paste(response, "~", paste(predictors_use, collapse = " + ")))

  lm_full <- tryCatch(lm(formula_full, data = model_df), error = function(e) NULL)
  if (is.null(lm_full)) {
    return(NULL)
  }

  lm_aic <- tryCatch(stepAIC(lm_full, direction = "backward", trace = 0), error = function(e) NULL)
  if (is.null(lm_aic)) {
    return(NULL)
  }

  retained <- setdiff(names(coef(lm_aic)), "(Intercept)")
  if (length(retained) == 0) {
    return(NULL)
  }

  lm_summary <- summary(lm_aic)

  vif_df <- if (length(retained) > 1) {
    vif_vals <- tryCatch(vif(lm_aic), error = function(e) NULL)
    if (is.null(vif_vals)) {
      data.frame(variable = retained, VIF = NA_real_)
    } else {
      data.frame(variable = names(vif_vals), VIF = as.numeric(vif_vals))
    }
  } else {
    data.frame(variable = retained, VIF = NA_real_)
  }

  coef_df <- as.data.frame(lm_summary$coefficients)
  coef_df$variable <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df <- coef_df %>% filter(variable != "(Intercept)")

  beta_df_data <- model_df
  beta_df_data[[response]] <- as.numeric(scale(beta_df_data[[response]]))
  lm_beta <- lm(as.formula(paste(response, "~", paste(retained, collapse = " + "))), data = beta_df_data)
  beta_df <- data.frame(
    variable = names(coef(lm_beta))[-1],
    beta_std = as.numeric(coef(lm_beta)[-1])
  )

  coef_out <- coef_df %>%
    left_join(vif_df, by = "variable") %>%
    left_join(beta_df, by = "variable") %>%
    mutate(
      response = response,
      R2 = lm_summary$r.squared,
      R2_adj = lm_summary$adj.r.squared,
      RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
      AIC = AIC(lm_aic),
      n = nrow(model_df)
    ) %>%
    dplyr::select(response, variable, beta_std, `Pr(>|t|)`, VIF, R2, R2_adj, RMSE, AIC, n)

  names(coef_out) <- c("Response", "Predictor", "Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "n")

  summary_out <- tibble(
    Response = response,
    Predictors_Final = paste(retained, collapse = "; "),
    R2 = lm_summary$r.squared,
    R2_adj = lm_summary$adj.r.squared,
    RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
    AIC = AIC(lm_aic),
    n = nrow(model_df)
  )

  list(model = lm_aic, data = model_df, coef = coef_out, summary = summary_out)
}

model_results <- list()
model_summary <- list()

for (response in response_vars) {
  candidate_fits <- list()

  for (i in seq_along(candidate_predictor_sets)) {
    fit_obj <- fit_one_model(
      merged_data,
      response,
      candidate_predictor_sets[[i]],
      use_scaled_predictors = MODEL_USE_SCALED_PREDICTORS,
      min_n = MODEL_MIN_N
    )

    if (!is.null(fit_obj)) {
      fit_obj$coef <- fit_obj$coef %>% mutate(Candidate_Set = i)
      fit_obj$summary <- fit_obj$summary %>% mutate(Candidate_Set = i)
      candidate_fits[[length(candidate_fits) + 1]] <- fit_obj
    }
  }

  if (length(candidate_fits) == 0) {
    next
  }

  best_idx <- which.min(sapply(candidate_fits, function(x) x$summary$AIC))
  best_fit <- candidate_fits[[best_idx]]

  model_results[[response]] <- best_fit$coef
  model_summary[[response]] <- best_fit$summary
}

model_results_combined <- bind_rows(model_results)
model_summary_combined <- bind_rows(model_summary)

cor_data <- merged_data %>%
  dplyr::select(all_of(unique(c(response_vars, storage_predictors_all))))

if (ncol(cor_data) >= 2) {
  cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")
  cor_response_storage <- cor_matrix[response_vars, storage_predictors_all, drop = FALSE]
  write.csv(cor_response_storage,
            file.path(output_dir, "Correlations_Storage_Thermal_LowFlow.csv"),
            row.names = TRUE)
}

write.csv(model_results_combined,
          file.path(output_dir, "Storage_Thermal_LowFlow_Models.csv"),
          row.names = FALSE)

write.csv(model_summary_combined,
          file.path(output_dir, "Storage_Thermal_LowFlow_ModelSummary.csv"),
          row.names = FALSE)
