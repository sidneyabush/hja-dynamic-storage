# Use stepwise regression to identify which watershed attributes.
# Inputs: No direct CSV file reads in this script.
# Author: Based on Pamela Sullivan/Keira Johnson code, adapted by Sidney Bush
# Date: 2026-01-23

library(dplyr)
library(readr)
library(ggplot2)
library(MASS)         # for stepAIC()
library(car)          # for vif()
library(tidyr)

# Clear environment
rm(list = ls())

# Source configuration (paths, site definitions, water year range)
# Get script directory (works with source() and Rscript)
# Load project config
source("config.R")

theme_set(theme_pub(base_size = 12))

# Use configuration values
site_order <- SITE_ORDER_HYDROMETRIC
base_dir   <- BASE_DATA_DIR
output_dir <- OUT_MODELS_WATERSHED_CHAR_STORAGE_MLR_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
file_prefix <- "watershed_char_storage_mlr"

# Remove older PCA screening output for catchment-characteristics MLR.
unlink(file.path(output_dir, paste0(file_prefix, "_pca_screening.csv")))

VIF_THRESHOLD <- 10
LOW_N_THRESHOLD_CATCH <- 10
USE_MUMIN_LOOCV <- TRUE
STRICT_METHOD_LABEL <- "strict"
STRICT_USE_SCALED_PREDICTORS <- TRUE
STRICT_USE_ITERATIVE_VIF <- TRUE
STRICT_ENFORCE_CORRELATED_EXCLUSIONS <- TRUE
EXCLUDE_SITES_BY_OUTCOME <- list(
  "CHS_mean" = CHS_EXCLUDE_SITES
)

env_vif_threshold <- suppressWarnings(as.numeric(Sys.getenv("HJA_MLR_VIF_THRESHOLD", unset = "")))
if (!is.na(env_vif_threshold) && is.finite(env_vif_threshold) && env_vif_threshold > 0) {
  VIF_THRESHOLD <- env_vif_threshold
}

if (USE_MUMIN_LOOCV && !requireNamespace("MuMIn", quietly = TRUE)) {
  stop("MuMIn is required for catchment-model LOOCV. Install with: install.packages('MuMIn')")
}

# LOAD SITE-AVERAGED DATA

master_dir <- file.path(OUTPUT_DIR, "master")
site_file <- file.path(master_dir, MASTER_SITE_FILE)

HJA_Ave <- read_csv(
  site_file,
  show_col_types = FALSE
) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

if (!("basin_slope" %in% names(HJA_Ave)) && ("Slope_mean" %in% names(HJA_Ave))) {
  HJA_Ave <- HJA_Ave %>% mutate(basin_slope = Slope_mean)
}

# DEFINE OUTCOME VARIABLES (STORAGE METRICS)
# Note: Q5norm and CV_Q5norm are NOT storage metrics - they are response variables
# Storage metrics by type (using method abbreviations):
#   Dynamic: RBI, RCS, FDC, SD
#   Mobile: CHS, MTT1, MTT2, Fyw, DR
#   Extended Dynamic: WB

outcome_vars <- c(
  "RCS_mean",   # RCS - Recession Curve Slope - Dynamic
  "RBI_mean",   # RBI - Richards-Baker Index - Dynamic
  "FDC_mean",   # FDC - Flow Duration Curve - Dynamic
  "SD_mean",    # SD  - Storage-Discharge - Dynamic
  "CHS_mean",   # CHS - Chemical Hydrograph Separation - Mobile
  "MTT1",       # MTT1 - Mean Transit Time period 1 - Mobile
  "MTT2",       # MTT2 - Mean Transit Time period 2 - Mobile
  "Fyw",        # Fyw - Young Water Fraction - Mobile
  "DR",         # DR  - Damping Ratio - Mobile
  "WB_mean"     # WB  - Water Balance - Extended Dynamic
)

# DEFINE PREDICTOR VARIABLES (CATCHMENT ATTRIBUTES)
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

# RUN STEPWISE REGRESSION FOR EACH STORAGE METRIC

fit_mlr_with_vif <- function(data_in, outcome, predictors,
                             use_scaled_predictors = FALSE,
                             use_iterative_vif = FALSE,
                             vif_threshold = 10,
                             enforce_correlated_exclusions = FALSE) {
  calc_loocv_stats <- function(model_formula, df) {
    n <- nrow(df)
    if (n < 6) {
      return(list(
        rmse_loocv = NA_real_,
        r2_loocv = NA_real_,
        rmse_loocv_mean_runs = NA_real_
      ))
    }

    rmse_loocv <- NA_real_
    if (USE_MUMIN_LOOCV) {
      model_obj <- tryCatch(lm(model_formula, data = df), error = function(e) NULL)
      if (!is.null(model_obj)) {
        rmse_loocv <- tryCatch(
          as.numeric(MuMIn::loo(model_obj, type = "rmse")),
          error = function(e) NA_real_
        )
      }
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
      return(list(
        rmse_loocv = NA_real_,
        r2_loocv = NA_real_,
        rmse_loocv_mean_runs = NA_real_
      ))
    }
    errs <- obs[valid] - preds[valid]
    if (!is.finite(rmse_loocv)) {
      rmse_loocv <- sqrt(mean(errs^2, na.rm = TRUE))
    }
    rmse_loocv_mean_runs <- mean(abs(errs), na.rm = TRUE)
    sst <- sum((obs[valid] - mean(obs[valid], na.rm = TRUE))^2, na.rm = TRUE)
    sse <- sum(errs^2, na.rm = TRUE)
    r2_loocv <- ifelse(sst > 0, 1 - (sse / sst), NA_real_)
    list(
      rmse_loocv = rmse_loocv,
      r2_loocv = r2_loocv,
      rmse_loocv_mean_runs = rmse_loocv_mean_runs
    )
  }

  model_input <- data_in
  if ("site" %in% names(model_input)) {
    excluded_sites <- EXCLUDE_SITES_BY_OUTCOME[[outcome]]
    if (!is.null(excluded_sites) && length(excluded_sites) > 0) {
      model_input <- model_input %>%
        filter(!(site %in% excluded_sites))
    }
  }

  predictor_keep <- predictors[predictors %in% names(model_input)]
  if (!(outcome %in% names(model_input)) || length(predictor_keep) == 0) {
    return(NULL)
  }

  model_df <- model_input %>%
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
    # Refit while removing the most collinear predictor until VIF is acceptable.
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

  # Optional rule set: avoid known predictor pairs that are strongly correlated.
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
        # If ash and lava are both present, keep the lower-AICc alternative.
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
        # If both landslide metrics are present, keep the lower-AICc alternative.
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
      RMSE_LOOCV_MEAN_RUNS = loocv$rmse_loocv_mean_runs,
      R2_LOOCV = loocv$r2_loocv,
      delta_RMSE_LOOCV_minus_model = RMSE_LOOCV - RMSE,
      delta_RMSE_LOOCV_mean_runs_minus_model = RMSE_LOOCV_MEAN_RUNS - RMSE
    )

  list(
    coefficients = result_df,
    summary = summary_df,
    model = lm_aic,
    data = model_df,
    selection_path = tryCatch(as.data.frame(lm_aic$anova), error = function(e) NULL)
  )
}

compute_residual_diagnostics <- function(model_obj, model_df) {
  pull_scalar <- function(obj, key) {
    val <- tryCatch(obj[[key]], error = function(e) NULL)
    num <- suppressWarnings(as.numeric(val))
    if (length(num) >= 1 && is.finite(num[1])) {
      return(num[1])
    }
    NA_real_
  }

  resid_vals <- residuals(model_obj)
  resid_vals <- resid_vals[is.finite(resid_vals)]

  shapiro_w <- NA_real_
  shapiro_p <- NA_real_
  if (length(resid_vals) >= 3 && length(resid_vals) <= 5000) {
    sh <- tryCatch(shapiro.test(resid_vals), error = function(e) NULL)
    if (!is.null(sh)) {
      shapiro_w <- suppressWarnings(as.numeric(unname(sh$statistic)))
      shapiro_p <- suppressWarnings(as.numeric(sh$p.value))
    }
  }

  ncv_chisq <- NA_real_
  ncv_p <- NA_real_
  ncv <- tryCatch(car::ncvTest(model_obj), error = function(e) NULL)
  if (!is.null(ncv)) {
    ncv_chisq <- pull_scalar(ncv, "Chisquare")
    if (!is.finite(ncv_chisq)) {
      ncv_chisq <- pull_scalar(ncv, "ChiSquare")
    }
    ncv_p <- pull_scalar(ncv, "p")
  }

  tibble(
    n_residuals = length(resid_vals),
    shapiro_W = shapiro_w,
    shapiro_p = shapiro_p,
    ncv_chisq = ncv_chisq,
    ncv_p = ncv_p,
    normality_pass_p05 = ifelse(is.finite(shapiro_p), shapiro_p > 0.05, NA),
    homoscedasticity_pass_p05 = ifelse(is.finite(ncv_p), ncv_p > 0.05, NA)
  )
}

run_strict_method <- function() {
  results_list <- list()
  model_summary_list <- list()
  selection_path_list <- list()
  model_selection_list <- list()
  diagnostics_list <- list()

  # Candidate predictor sets for model comparison and best-model selection.
  candidate_predictor_sets <- unlist(
    lapply(seq_along(predictor_vars), function(k) {
      combn(predictor_vars, k, simplify = FALSE)
    }),
    recursive = FALSE
  )

  for (outcome in outcome_vars) {
    candidate_fits <- list()
    candidate_summaries <- list()

    for (i in seq_along(candidate_predictor_sets)) {
      fit_obj <- fit_mlr_with_vif(
        HJA_Ave,
        outcome,
        candidate_predictor_sets[[i]],
        use_scaled_predictors = STRICT_USE_SCALED_PREDICTORS,
        use_iterative_vif = STRICT_USE_ITERATIVE_VIF,
        vif_threshold = VIF_THRESHOLD,
        enforce_correlated_exclusions = STRICT_ENFORCE_CORRELATED_EXCLUSIONS
      )
      if (is.null(fit_obj)) next

      fit_obj$summary <- fit_obj$summary %>% mutate(Candidate_Set = i, Method = STRICT_METHOD_LABEL)
      fit_obj$coefficients <- fit_obj$coefficients %>% mutate(Candidate_Set = i, Method = STRICT_METHOD_LABEL)

      candidate_fits[[length(candidate_fits) + 1]] <- fit_obj
      candidate_summaries[[length(candidate_summaries) + 1]] <- fit_obj$summary
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
      sp$Method <- STRICT_METHOD_LABEL
      sp$Candidate_Set <- best_fit$summary$Candidate_Set[1]
      selection_path_list[[outcome]] <- sp
    }
    model_summary_list[[outcome]] <- best_fit$summary

    diag_tbl <- compute_residual_diagnostics(best_fit$model, best_fit$data) %>%
      mutate(
        Outcome = ifelse(outcome == "DR", "DR_mean", outcome),
        Predictors_Final = best_fit$summary$Predictors_Final[1],
        N = best_fit$summary$N[1],
        Method = STRICT_METHOD_LABEL
      ) %>%
      dplyr::select(
        Outcome, Predictors_Final, N, n_residuals,
        shapiro_W, shapiro_p, ncv_chisq, ncv_p,
        normality_pass_p05, homoscedasticity_pass_p05, Method
      )
    diagnostics_list[[outcome]] <- diag_tbl
  }

  list(
    results = bind_rows(results_list),
    summary = bind_rows(model_summary_list),
    selection = bind_rows(selection_path_list),
    model_selection = bind_rows(model_selection_list),
    diagnostics = bind_rows(diagnostics_list)
  )
}

model_run <- run_strict_method()

model_run$summary <- model_run$summary %>%
  mutate(
    low_n_flag = N < LOW_N_THRESHOLD_CATCH,
    confidence_note = ifelse(low_n_flag, "lower confidence (small n)", "standard confidence")
  )

beta_matrix <- model_run$results %>%
  dplyr::select(Predictor, Outcome, Beta_Std) %>%
  tidyr::pivot_wider(names_from = Outcome, values_from = Beta_Std)

# Save outputs as primary for downstream plotting
write.csv(model_run$results,
          file.path(output_dir, paste0(file_prefix, "_results.csv")),
          row.names = FALSE)
write.csv(model_run$summary,
          file.path(output_dir, paste0(file_prefix, "_summary.csv")),
          row.names = FALSE)
write.csv(beta_matrix,
          file.path(output_dir, paste0(file_prefix, "_beta_matrix.csv")),
          row.names = FALSE)
write.csv(model_run$selection,
          file.path(output_dir, paste0(file_prefix, "_stepwise_path.csv")),
          row.names = FALSE)
write.csv(model_run$model_selection,
          file.path(output_dir, paste0(file_prefix, "_model_selection.csv")),
          row.names = FALSE)
diagnostics_out <- model_run$diagnostics %>%
  left_join(
    model_run$summary %>% dplyr::select(Outcome, low_n_flag, confidence_note),
    by = "Outcome"
  )
write.csv(
  diagnostics_out,
  file.path(output_dir, paste0(file_prefix, "_diagnostics.csv")),
  row.names = FALSE
)
aicc_lt2 <- model_run$model_selection %>%
  filter(is.finite(delta_AICc), delta_AICc <= 2) %>%
  mutate(Predictors_Final_Label = label_catchment_predictor_list(Predictors_Final)) %>%
  arrange(Outcome, delta_AICc, AICc, Candidate_Set)
write.csv(
  aicc_lt2,
  file.path(output_dir, paste0(file_prefix, "_aicc_lt2.csv")),
  row.names = FALSE
)
old_model_avg_file <- file.path(output_dir, paste0(file_prefix, "_model_avg_coef.csv"))
if (file.exists(old_model_avg_file)) unlink(old_model_avg_file)
flags <- model_run$summary %>%
  filter(constraint_flag) %>%
  arrange(Outcome)
write.csv(flags,
          file.path(output_dir, paste0(file_prefix, "_corr_flags.csv")),
          row.names = FALSE)

# Explicit LOOCV validation output (models + tables validation folders)
loocv_validation <- model_run$summary %>%
  transmute(
    model_family = "watershed_char_storage_mlr",
    outcome = Outcome,
    n = N,
    rmse_model = RMSE,
    rmse_loocv = RMSE_LOOCV,
    rmse_loocv_mean_runs = RMSE_LOOCV_MEAN_RUNS,
    r2_model = R2,
    r2_adj_model = R2_adj,
    r2_loocv = R2_LOOCV,
    delta_rmse_loocv_minus_model = delta_RMSE_LOOCV_minus_model,
    delta_rmse_loocv_mean_runs_minus_model = delta_RMSE_LOOCV_mean_runs_minus_model
  ) %>%
  arrange(outcome)

write.csv(
  loocv_validation,
  file.path(OUT_STATS_VALIDATION_DIR, "watershed_char_storage_mlr_loocv_validation.csv"),
  row.names = FALSE
)

if (!dir.exists(file.path(OUT_TABLES_DIR, "validation"))) {
  dir.create(file.path(OUT_TABLES_DIR, "validation"), recursive = TRUE, showWarnings = FALSE)
}
write.csv(
  loocv_validation,
  file.path(OUT_TABLES_DIR, "validation", "watershed_char_storage_mlr_loocv_validation.csv"),
  row.names = FALSE
)
if (!dir.exists(OUT_TABLES_MLR_DIR)) {
  dir.create(OUT_TABLES_MLR_DIR, recursive = TRUE, showWarnings = FALSE)
}
write.csv(
  aicc_lt2,
  file.path(OUT_TABLES_MLR_DIR, "watershed_char_storage_mlr_aicc_lt2.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics_out,
  file.path(OUT_TABLES_MLR_DIR, "watershed_char_storage_mlr_diagnostics.csv"),
  row.names = FALSE
)
