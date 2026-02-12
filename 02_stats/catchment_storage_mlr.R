# -----------------------------------------------------------------------------
# Multiple Linear Regression: Watershed Characteristics to Storage Metrics
# -----------------------------------------------------------------------------
# Purpose: Use stepwise regression to identify which watershed attributes
#          best predict hydrometric storage metrics
#
# Workflow:
#   1. Load site-averaged storage metrics and watershed attributes
#   2. For each storage metric:
#      a. Build full model with all watershed predictors
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
output_dir <- OUT_STATS_MLR_CATCH_CHARS_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
file_prefix <- "catch_chars_storage_mlr"

# Remove deprecated PCA-screening output for catchment-characteristics MLR.
unlink(file.path(output_dir, paste0(file_prefix, "_pca_screening.csv")))

VIF_THRESHOLD <- 10
LOW_N_THRESHOLD_CATCH <- 10

env_vif_threshold <- suppressWarnings(as.numeric(Sys.getenv("HJA_MLR_VIF_THRESHOLD", unset = "")))
if (!is.na(env_vif_threshold) && is.finite(env_vif_threshold) && env_vif_threshold > 0) {
  VIF_THRESHOLD <- env_vif_threshold
}

# -----------------------------------------------------------------------------
# 2. LOAD SITE-AVERAGED DATA
# -----------------------------------------------------------------------------

site_file <- file.path(OUT_MASTER_DIR, MASTER_SITE_FILE)
if (!file.exists(site_file)) {
  site_file <- file.path(OUTPUT_DIR, MASTER_SITE_FILE)
}

HJA_Ave <- read_csv(
  site_file,
  show_col_types = FALSE
) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

if (!("basin_slope" %in% names(HJA_Ave)) && ("Slope_mean" %in% names(HJA_Ave))) {
  HJA_Ave <- HJA_Ave %>% mutate(basin_slope = Slope_mean)
}

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
  "basin_slope",
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

fit_mlr_with_vif <- function(data_in, outcome, predictors,
                             use_scaled_predictors = FALSE,
                             use_iterative_vif = FALSE,
                             vif_threshold = 10,
                             enforce_correlated_exclusions = FALSE) {
  calc_loocv_stats <- function(model_formula, df) {
    n <- nrow(df)
    if (n < 6) {
      return(list(rmse_loocv = NA_real_, r2_loocv = NA_real_))
    }
    response <- all.vars(model_formula)[1]
    obs <- df[[response]]
    preds <- rep(NA_real_, n)

    for (i in seq_len(n)) {
      fit_i <- tryCatch(lm(model_formula, data = df[-i, , drop = FALSE]), error = function(e) NULL)
      if (!is.null(fit_i)) {
        preds[i] <- tryCatch(as.numeric(predict(fit_i, newdata = df[i, , drop = FALSE])), error = function(e) NA_real_)
      }
    }

    valid <- is.finite(obs) & is.finite(preds)
    if (sum(valid) < 3) {
      return(list(rmse_loocv = NA_real_, r2_loocv = NA_real_))
    }
    rmse_loocv <- sqrt(mean((obs[valid] - preds[valid])^2, na.rm = TRUE))
    sst <- sum((obs[valid] - mean(obs[valid], na.rm = TRUE))^2, na.rm = TRUE)
    sse <- sum((obs[valid] - preds[valid])^2, na.rm = TRUE)
    r2_loocv <- ifelse(sst > 0, 1 - (sse / sst), NA_real_)
    list(rmse_loocv = rmse_loocv, r2_loocv = r2_loocv)
  }

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

  # Optional enforcement of known high-correlation predictor exclusions (Perry et al., 2025).
  if (enforce_correlated_exclusions) {
    fit_candidate_from_retained <- function(retained_set) {
      if (length(retained_set) == 0) return(NULL)
      lm_next <- tryCatch(
        lm(as.formula(paste(outcome, "~", paste(retained_set, collapse = " + "))), data = model_df),
        error = function(e) NULL
      )
      if (is.null(lm_next)) return(NULL)
      lm_next <- tryCatch(stepAIC(lm_next, direction = "backward", trace = 0), error = function(e) lm_next)
      retained_next <- setdiff(names(coef(lm_next)), "(Intercept)")
      if (length(retained_next) == 0) return(NULL)
      n_obs_next <- nrow(model_df)
      k_params_next <- length(coef(lm_next)) + 1
      aic_next <- AIC(lm_next)
      aicc_next <- if ((n_obs_next - k_params_next - 1) > 0) {
        aic_next + (2 * k_params_next * (k_params_next + 1)) / (n_obs_next - k_params_next - 1)
      } else {
        Inf
      }
      list(model = lm_next, aicc = aicc_next)
    }

    repeat {
      retained_vars <- setdiff(names(coef(lm_aic)), "(Intercept)")
      has_ash <- "Ash_Per" %in% retained_vars
      has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained_vars)
      has_ls_both <- all(c("Landslide_Total", "Landslide_Young") %in% retained_vars)

      if (!(has_ash && has_lava) && !has_ls_both) {
        break
      }

      if (has_ash && has_lava) {
        # Resolve Ash/Lava collinearity by selecting the lowest-AICc candidate.
        drop_options <- c("Ash_Per")
        if ("Lava1_per" %in% retained_vars) drop_options <- c(drop_options, "Lava1_per")
        if ("Lava2_per" %in% retained_vars) drop_options <- c(drop_options, "Lava2_per")

        cand_models <- lapply(drop_options, function(v) {
          retained_next <- setdiff(retained_vars, v)
          fit_candidate_from_retained(retained_next)
        })
        valid_idx <- which(vapply(cand_models, function(x) !is.null(x), logical(1)))
        if (length(valid_idx) == 0) break
        best_idx <- valid_idx[which.min(vapply(cand_models[valid_idx], function(x) x$aicc, numeric(1)))]
        lm_aic <- cand_models[[best_idx]]$model
      } else if (has_ls_both) {
        # Resolve landslide-pair collinearity by selecting the lowest-AICc candidate.
        drop_options <- c("Landslide_Total", "Landslide_Young")
        cand_models <- lapply(drop_options, function(v) {
          retained_next <- setdiff(retained_vars, v)
          fit_candidate_from_retained(retained_next)
        })
        valid_idx <- which(vapply(cand_models, function(x) !is.null(x), logical(1)))
        if (length(valid_idx) == 0) break
        best_idx <- valid_idx[which.min(vapply(cand_models[valid_idx], function(x) x$aicc, numeric(1)))]
        lm_aic <- cand_models[[best_idx]]$model
      }
    }
  }

  retained_vars <- setdiff(names(coef(lm_aic)), "(Intercept)")
  if (length(retained_vars) == 0) {
    return(NULL)
  }

  has_ash <- "Ash_Per" %in% retained_vars
  has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained_vars)
  flag_lava_ash <- has_ash && has_lava
  flag_ls_total_young <- all(c("Landslide_Total", "Landslide_Young") %in% retained_vars)
  constraint_flag <- flag_lava_ash || flag_ls_total_young

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

  n_obs <- nrow(model_df)
  k_params <- length(coef(lm_aic)) + 1
  aic_val <- AIC(lm_aic)
  aicc_val <- if ((n_obs - k_params - 1) > 0) {
    aic_val + (2 * k_params * (k_params + 1)) / (n_obs - k_params - 1)
  } else {
    NA_real_
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
      outcome = ifelse(outcome == "DR", "DR_mean", outcome),
      R2 = lm_summary$r.squared,
      R2_adj = lm_summary$adj.r.squared,
      RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
      AIC = aic_val,
      AICc = aicc_val
    ) %>%
    dplyr::select(outcome, variable, beta_std, `Pr(>|t|)`, VIF, R2, R2_adj, RMSE, AIC, AICc)

  names(result_df) <- c("Outcome", "Predictor", "Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "AICc")

  summary_df <- tibble(
    Outcome = ifelse(outcome == "DR", "DR_mean", outcome),
    N = nrow(model_df),
    PCA_PCs_Retained = NA_integer_,
    Predictors_PCA = NA_character_,
    Predictors_Final = paste(retained_vars, collapse = "; "),
    R2 = lm_summary$r.squared,
    R2_adj = lm_summary$adj.r.squared,
    RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
    AIC = aic_val,
    AICc = aicc_val,
    flag_lava_and_ash = flag_lava_ash,
    flag_landslide_total_and_young = flag_ls_total_young,
    constraint_flag = constraint_flag
  )

  loocv <- calc_loocv_stats(formula(lm_aic), model_df)
  summary_df <- summary_df %>%
    mutate(
      RMSE_LOOCV = loocv$rmse_loocv,
      R2_LOOCV = loocv$r2_loocv,
      delta_RMSE_LOOCV_minus_model = RMSE_LOOCV - RMSE
    )

  list(
    coefficients = result_df,
    summary = summary_df,
    model = lm_aic,
    data = model_df,
    selection_path = tryCatch(as.data.frame(lm_aic$anova), error = function(e) NULL)
  )
}

run_method <- function(method_name, use_scaled_predictors, use_iterative_vif, enforce_correlated_exclusions = FALSE) {
  results_list <- list()
  model_summary_list <- list()
  selection_path_list <- list()
  model_selection_list <- list()
  model_avg_list <- list()

  # Candidate predictor sets for model comparison and averaging.
  candidate_predictor_sets <- unlist(
    lapply(seq_along(predictor_vars), function(k) {
      combn(predictor_vars, k, simplify = FALSE)
    }),
    recursive = FALSE
  )

  for (outcome in outcome_vars) {
    candidate_fits <- list()
    candidate_summaries <- list()
    candidate_coefs <- list()

    for (i in seq_along(candidate_predictor_sets)) {
      fit_obj <- fit_mlr_with_vif(
        HJA_Ave,
        outcome,
        candidate_predictor_sets[[i]],
        use_scaled_predictors = use_scaled_predictors,
        use_iterative_vif = use_iterative_vif,
        vif_threshold = VIF_THRESHOLD,
        enforce_correlated_exclusions = enforce_correlated_exclusions
      )
      if (is.null(fit_obj)) next

      fit_obj$summary <- fit_obj$summary %>% mutate(Candidate_Set = i, Method = method_name)
      fit_obj$coefficients <- fit_obj$coefficients %>% mutate(Candidate_Set = i, Method = method_name)

      candidate_fits[[length(candidate_fits) + 1]] <- fit_obj
      candidate_summaries[[length(candidate_summaries) + 1]] <- fit_obj$summary
      candidate_coefs[[length(candidate_coefs) + 1]] <- fit_obj$coefficients
    }

    if (length(candidate_fits) == 0) next

    cand_tbl <- bind_rows(candidate_summaries) %>%
      mutate(delta_AICc = AICc - min(AICc, na.rm = TRUE)) %>%
      arrange(AICc)
    model_selection_list[[outcome]] <- cand_tbl

    best_idx <- which.min(sapply(candidate_fits, function(x) x$summary$AICc))
    best_fit <- candidate_fits[[best_idx]]
    results_list[[outcome]] <- best_fit$coefficients

    if (!is.null(best_fit$selection_path) && nrow(best_fit$selection_path) > 0) {
      sp <- best_fit$selection_path
      sp$step <- seq_len(nrow(sp))
      sp$Outcome <- ifelse(outcome == "DR", "DR_mean", outcome)
      sp$Method <- method_name
      sp$Candidate_Set <- best_fit$summary$Candidate_Set[1]
      selection_path_list[[outcome]] <- sp
    }

    # Model-averaged standardized coefficients for delta AICc <= 2.
    cand_weights <- cand_tbl %>%
      filter(is.finite(AICc), delta_AICc <= 2) %>%
      mutate(w_raw = exp(-0.5 * delta_AICc),
             weight = w_raw / sum(w_raw, na.rm = TRUE)) %>%
      dplyr::select(Candidate_Set, Outcome, weight)

    best_summary <- best_fit$summary %>%
      mutate(
        R2_weighted_AICc_LTE2 = NA_real_,
        R2_adj_weighted_AICc_LTE2 = NA_real_
      )

    if (nrow(cand_weights) > 0) {
      weighted_fit <- cand_tbl %>%
        inner_join(cand_weights, by = c("Candidate_Set", "Outcome")) %>%
        summarise(
          R2_weighted_AICc_LTE2 = sum(R2 * weight, na.rm = TRUE),
          R2_adj_weighted_AICc_LTE2 = sum(R2_adj * weight, na.rm = TRUE),
          .groups = "drop"
        )
      if (nrow(weighted_fit) == 1) {
        best_summary <- best_summary %>%
          mutate(
            R2_weighted_AICc_LTE2 = weighted_fit$R2_weighted_AICc_LTE2[1],
            R2_adj_weighted_AICc_LTE2 = weighted_fit$R2_adj_weighted_AICc_LTE2[1]
          )
      }
    }
    model_summary_list[[outcome]] <- best_summary

    if (nrow(cand_weights) > 0) {
      coef_tbl <- bind_rows(candidate_coefs)
      predictor_universe <- sort(unique(coef_tbl$Predictor))
      avg_tbl <- tibble(Predictor = predictor_universe) %>%
        left_join(
          coef_tbl %>%
            inner_join(cand_weights, by = c("Candidate_Set", "Outcome")) %>%
            group_by(Predictor) %>%
            summarise(
              Beta_Std_Avg = sum(Beta_Std * weight, na.rm = TRUE),
              Inclusion_Weight = sum(weight, na.rm = TRUE),
              Models_Present = n_distinct(Candidate_Set),
              .groups = "drop"
            ),
          by = "Predictor"
        ) %>%
        mutate(
          Outcome = ifelse(outcome == "DR", "DR_mean", outcome),
          Beta_Std_Avg = ifelse(is.na(Beta_Std_Avg), 0, Beta_Std_Avg),
          Inclusion_Weight = ifelse(is.na(Inclusion_Weight), 0, Inclusion_Weight),
          Models_DeltaAICc_LTE2 = nrow(cand_weights)
        ) %>%
        dplyr::select(Outcome, Predictor, Beta_Std_Avg, Inclusion_Weight, Models_Present, Models_DeltaAICc_LTE2)

      model_avg_list[[outcome]] <- avg_tbl
    }
  }

  list(
    results = bind_rows(results_list),
    summary = bind_rows(model_summary_list),
    selection = bind_rows(selection_path_list),
    model_selection = bind_rows(model_selection_list),
    model_avg = bind_rows(model_avg_list)
  )
}

strict <- run_method(
  method_name = "strict",
  use_scaled_predictors = TRUE,
  use_iterative_vif = TRUE,
  enforce_correlated_exclusions = TRUE
)

strict$summary <- strict$summary %>%
  mutate(
    low_n_flag = N < LOW_N_THRESHOLD_CATCH,
    confidence_note = ifelse(low_n_flag, "lower confidence (small n)", "standard confidence")
  )

beta_matrix_strict <- strict$results %>%
  dplyr::select(Predictor, Outcome, Beta_Std) %>%
  tidyr::pivot_wider(names_from = Outcome, values_from = Beta_Std)

# Save strict outputs as primary for downstream plotting
write.csv(strict$results,
          file.path(output_dir, paste0(file_prefix, "_results.csv")),
          row.names = FALSE)
write.csv(strict$summary,
          file.path(output_dir, paste0(file_prefix, "_summary.csv")),
          row.names = FALSE)
write.csv(beta_matrix_strict,
          file.path(output_dir, paste0(file_prefix, "_beta_matrix.csv")),
          row.names = FALSE)

write.csv(strict$results,
          file.path(output_dir, paste0(file_prefix, "_results_strict.csv")),
          row.names = FALSE)
write.csv(strict$summary,
          file.path(output_dir, paste0(file_prefix, "_summary_strict.csv")),
          row.names = FALSE)
write.csv(strict$selection,
          file.path(output_dir, paste0(file_prefix, "_stepwise_path.csv")),
          row.names = FALSE)
write.csv(strict$model_selection,
          file.path(output_dir, paste0(file_prefix, "_model_selection.csv")),
          row.names = FALSE)
write.csv(strict$model_avg,
          file.path(output_dir, paste0(file_prefix, "_model_avg_coef.csv")),
          row.names = FALSE)
flags <- strict$summary %>%
  filter(constraint_flag) %>%
  arrange(Method, Outcome)
write.csv(flags,
          file.path(output_dir, paste0(file_prefix, "_corr_flags.csv")),
          row.names = FALSE)
