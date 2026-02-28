# run catchment controls mlr and export only manuscript tables

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(MASS)
  library(car)
})

rm(list = ls())
source("config.R")

output_dir <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
file_prefix <- "catchment_char_storage_mlr"

VIF_THRESHOLD <- 10
LOW_N_THRESHOLD_CATCH <- 10

site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
if (!file.exists(site_file)) {
  stop("Missing required file: ", site_file)
}

site_df <- read_csv(site_file, show_col_types = FALSE) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

if (!("basin_slope" %in% names(site_df)) && ("Slope_mean" %in% names(site_df))) {
  site_df <- site_df %>% mutate(basin_slope = Slope_mean)
}

outcome_map <- c(
  "RBI" = "RBI_mean",
  "RCS" = "RCS_mean",
  "FDC" = "FDC_mean",
  "SD" = "SD_mean",
  "WB" = "WB_mean",
  "CHS" = "CHS_mean",
  "DR" = "DR",
  "Fyw" = "Fyw",
  "MTT1" = "MTT1",
  "MTT2" = "MTT2"
)
outcome_vars <- unname(outcome_map[STORAGE_METRIC_ORDER])

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

candidate_predictor_sets <- unlist(
  lapply(seq_along(predictor_vars), function(k) {
    combn(predictor_vars, k, simplify = FALSE)
  }),
  recursive = FALSE
)

calc_aicc <- function(model_obj, n_obs) {
  k_params <- length(coef(model_obj)) + 1
  aic_val <- AIC(model_obj)
  if ((n_obs - k_params - 1) <= 0) {
    return(NA_real_)
  }
  aic_val + (2 * k_params * (k_params + 1)) / (n_obs - k_params - 1)
}

calc_loocv_stats <- function(model_formula, df) {
  n <- nrow(df)
  if (n < 6) {
    return(list(
      rmse_loocv = NA_real_,
      rmse_loocv_mean_runs = NA_real_,
      r2_loocv = NA_real_
    ))
  }

  response <- all.vars(model_formula)[1]
  obs <- df[[response]]
  preds <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    fit_i <- tryCatch(lm(model_formula, data = df[-i, , drop = FALSE]), error = function(e) NULL)
    if (!is.null(fit_i)) {
      preds[i] <- tryCatch(
        as.numeric(predict(fit_i, newdata = df[i, , drop = FALSE])),
        error = function(e) NA_real_
      )
    }
  }

  valid <- is.finite(obs) & is.finite(preds)
  if (sum(valid) < 3) {
    return(list(
      rmse_loocv = NA_real_,
      rmse_loocv_mean_runs = NA_real_,
      r2_loocv = NA_real_
    ))
  }

  errs <- obs[valid] - preds[valid]
  rmse_loocv <- sqrt(mean(errs^2, na.rm = TRUE))
  rmse_loocv_mean_runs <- mean(abs(errs), na.rm = TRUE)
  sst <- sum((obs[valid] - mean(obs[valid], na.rm = TRUE))^2, na.rm = TRUE)
  sse <- sum(errs^2, na.rm = TRUE)
  r2_loocv <- ifelse(sst > 0, 1 - sse / sst, NA_real_)

  list(
    rmse_loocv = rmse_loocv,
    rmse_loocv_mean_runs = rmse_loocv_mean_runs,
    r2_loocv = r2_loocv
  )
}

iterative_vif <- function(model_obj, outcome, model_df, threshold) {
  lm_current <- model_obj
  repeat {
    retained <- setdiff(names(coef(lm_current)), "(Intercept)")
    if (length(retained) <= 1) {
      break
    }

    vif_vals <- tryCatch(vif(lm_current), error = function(e) NULL)
    if (is.null(vif_vals)) {
      break
    }

    if (max(vif_vals, na.rm = TRUE) <= threshold) {
      break
    }

    drop_var <- names(which.max(vif_vals))
    retained_next <- setdiff(retained, drop_var)
    if (length(retained_next) == 0) {
      return(NULL)
    }

    lm_next <- lm(
      as.formula(paste(outcome, "~", paste(retained_next, collapse = " + "))),
      data = model_df
    )
    lm_current <- tryCatch(stepAIC(lm_next, direction = "backward", trace = 0), error = function(e) lm_next)
  }
  lm_current
}

refit_drop_and_score <- function(model_obj, outcome, model_df, drop_var) {
  retained <- setdiff(setdiff(names(coef(model_obj)), "(Intercept)"), drop_var)
  if (length(retained) == 0) {
    return(NULL)
  }

  fit <- tryCatch(
    lm(as.formula(paste(outcome, "~", paste(retained, collapse = " + "))), data = model_df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NULL)
  }
  fit <- tryCatch(stepAIC(fit, direction = "backward", trace = 0), error = function(e) fit)

  retained_fit <- setdiff(names(coef(fit)), "(Intercept)")
  if (length(retained_fit) == 0) {
    return(NULL)
  }

  list(model = fit, aicc = calc_aicc(fit, nrow(model_df)))
}

enforce_correlated_exclusions <- function(model_obj, outcome, model_df) {
  lm_current <- model_obj

  repeat {
    retained <- setdiff(names(coef(lm_current)), "(Intercept)")
    has_ash <- "Ash_Per" %in% retained
    has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained)
    has_both_landslide <- all(c("Landslide_Total", "Landslide_Young") %in% retained)

    if (!(has_ash && has_lava) && !has_both_landslide) {
      break
    }

    candidate_models <- list()

    if (has_ash && has_lava) {
      drop_options <- c("Ash_Per")
      if ("Lava1_per" %in% retained) drop_options <- c(drop_options, "Lava1_per")
      if ("Lava2_per" %in% retained) drop_options <- c(drop_options, "Lava2_per")
      candidate_models <- lapply(drop_options, function(v) refit_drop_and_score(lm_current, outcome, model_df, v))
    } else if (has_both_landslide) {
      drop_options <- c("Landslide_Total", "Landslide_Young")
      candidate_models <- lapply(drop_options, function(v) refit_drop_and_score(lm_current, outcome, model_df, v))
    }

    valid_idx <- which(vapply(candidate_models, function(x) !is.null(x), logical(1)))
    if (length(valid_idx) == 0) {
      break
    }

    best_idx <- valid_idx[which.min(vapply(candidate_models[valid_idx], function(x) x$aicc, numeric(1)))]
    lm_current <- candidate_models[[best_idx]]$model
  }

  lm_current
}

compute_residual_diagnostics <- function(model_obj) {
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

fit_candidate <- function(data_in, outcome, predictors) {
  model_input <- data_in
  if ("site" %in% names(model_input) && identical(outcome, "CHS_mean")) {
    model_input <- model_input %>% filter(!(site %in% CHS_EXCLUDE_SITES))
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

  model_df <- model_df %>%
    mutate(across(all_of(predictor_keep), ~ as.numeric(scale(.x))))

  lm_full <- tryCatch(
    lm(as.formula(paste(outcome, "~", paste(predictor_keep, collapse = " + "))), data = model_df),
    error = function(e) NULL
  )
  if (is.null(lm_full)) {
    return(NULL)
  }

  lm_step <- tryCatch(stepAIC(lm_full, direction = "backward", trace = 0), error = function(e) NULL)
  if (is.null(lm_step)) {
    return(NULL)
  }

  lm_step <- iterative_vif(lm_step, outcome, model_df, VIF_THRESHOLD)
  if (is.null(lm_step)) {
    return(NULL)
  }

  lm_step <- enforce_correlated_exclusions(lm_step, outcome, model_df)
  retained <- setdiff(names(coef(lm_step)), "(Intercept)")
  if (length(retained) == 0) {
    return(NULL)
  }

  vif_df <- if (length(retained) > 1) {
    vif_vals <- tryCatch(vif(lm_step), error = function(e) NULL)
    if (is.null(vif_vals)) {
      tibble(variable = retained, VIF = NA_real_)
    } else {
      tibble(variable = names(vif_vals), VIF = as.numeric(vif_vals))
    }
  } else {
    tibble(variable = retained, VIF = NA_real_)
  }

  fit_sum <- summary(lm_step)
  aic_val <- AIC(lm_step)
  aicc_val <- calc_aicc(lm_step, nrow(model_df))

  coef_df <- as.data.frame(fit_sum$coefficients)
  coef_df$variable <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df <- coef_df %>% filter(variable != "(Intercept)")

  beta_df_data <- model_df
  beta_df_data[[outcome]] <- as.numeric(scale(beta_df_data[[outcome]]))
  lm_beta <- lm(as.formula(paste(outcome, "~", paste(retained, collapse = " + "))), data = beta_df_data)
  beta_df <- tibble(
    variable = names(coef(lm_beta))[-1],
    beta_std = as.numeric(coef(lm_beta)[-1])
  )

  has_ash <- "Ash_Per" %in% retained
  has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained)
  has_both_landslide <- all(c("Landslide_Total", "Landslide_Young") %in% retained)

  coef_out <- coef_df %>%
    left_join(vif_df, by = "variable") %>%
    left_join(beta_df, by = "variable") %>%
    mutate(
      outcome = ifelse(outcome == "DR", "DR_mean", outcome),
      R2 = fit_sum$r.squared,
      R2_adj = fit_sum$adj.r.squared,
      RMSE = sqrt(mean(residuals(lm_step)^2, na.rm = TRUE)),
      AIC = aic_val,
      AICc = aicc_val
    ) %>%
    dplyr::select(outcome, variable, beta_std, `Pr(>|t|)`, VIF, R2, R2_adj, RMSE, AIC, AICc)

  names(coef_out) <- c("Outcome", "Predictor", "Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "AICc")

  loocv <- calc_loocv_stats(formula(lm_step), model_df)

  summary_out <- tibble(
    Outcome = ifelse(outcome == "DR", "DR_mean", outcome),
    N = nrow(model_df),
    Predictors_Final = paste(retained, collapse = "; "),
    R2 = fit_sum$r.squared,
    R2_adj = fit_sum$adj.r.squared,
    RMSE = sqrt(mean(residuals(lm_step)^2, na.rm = TRUE)),
    AIC = aic_val,
    AICc = aicc_val,
    RMSE_LOOCV = loocv$rmse_loocv,
    RMSE_LOOCV_MEAN_RUNS = loocv$rmse_loocv_mean_runs,
    R2_LOOCV = loocv$r2_loocv,
    delta_RMSE_LOOCV_minus_model = loocv$rmse_loocv - sqrt(mean(residuals(lm_step)^2, na.rm = TRUE)),
    delta_RMSE_LOOCV_mean_runs_minus_model = loocv$rmse_loocv_mean_runs - sqrt(mean(residuals(lm_step)^2, na.rm = TRUE)),
    flag_lava_and_ash = has_ash && has_lava,
    flag_landslide_total_and_young = has_both_landslide,
    constraint_flag = (has_ash && has_lava) || has_both_landslide
  )

  diagnostics_out <- compute_residual_diagnostics(lm_step) %>%
    mutate(
      Outcome = ifelse(outcome == "DR", "DR_mean", outcome),
      Predictors_Final = paste(retained, collapse = "; "),
      N = nrow(model_df)
    ) %>%
    dplyr::select(
      Outcome, Predictors_Final, N, n_residuals,
      shapiro_W, shapiro_p, ncv_chisq, ncv_p,
      normality_pass_p05, homoscedasticity_pass_p05
    )

  list(
    coefficients = coef_out,
    summary = summary_out,
    diagnostics = diagnostics_out
  )
}

results_list <- list()
summary_list <- list()
diagnostics_list <- list()
model_selection_list <- list()

for (outcome in outcome_vars) {
  candidate_fits <- list()
  candidate_summaries <- list()

  for (i in seq_along(candidate_predictor_sets)) {
    fit_obj <- fit_candidate(site_df, outcome, candidate_predictor_sets[[i]])
    if (is.null(fit_obj)) {
      next
    }

    fit_obj$summary <- fit_obj$summary %>% mutate(Candidate_Set = i)
    candidate_fits[[length(candidate_fits) + 1]] <- fit_obj
    candidate_summaries[[length(candidate_summaries) + 1]] <- fit_obj$summary
  }

  if (length(candidate_fits) == 0) {
    next
  }

  cand_tbl <- bind_rows(candidate_summaries) %>%
    mutate(delta_AICc = AICc - min(AICc, na.rm = TRUE)) %>%
    arrange(AICc)

  model_selection_list[[outcome]] <- cand_tbl
  best_idx <- which.min(vapply(candidate_fits, function(x) x$summary$AICc[1], numeric(1)))
  best_fit <- candidate_fits[[best_idx]]

  results_list[[outcome]] <- best_fit$coefficients
  summary_list[[outcome]] <- best_fit$summary %>% dplyr::select(-Candidate_Set)
  diagnostics_list[[outcome]] <- best_fit$diagnostics
}

model_run <- list(
  results = bind_rows(results_list),
  summary = bind_rows(summary_list),
  diagnostics = bind_rows(diagnostics_list),
  model_selection = bind_rows(model_selection_list)
)

model_run$summary <- model_run$summary %>%
  mutate(
    low_n_flag = N < LOW_N_THRESHOLD_CATCH,
    confidence_note = ifelse(low_n_flag, "lower confidence (small n)", "standard confidence")
  )

format_export_outcome <- function(x) {
  out <- gsub("_mean$", "", as.character(x))
  gsub("_", "", out, fixed = TRUE)
}

round_export_cols <- function(df, cols, digits = 3) {
  keep <- intersect(cols, names(df))
  if (length(keep) == 0) {
    return(df)
  }
  df %>% mutate(across(all_of(keep), ~ signif(.x, digits)))
}

results_export <- model_run$results %>%
  mutate(
    Outcome = format_export_outcome(Outcome),
    Predictor = label_catchment_predictor(Predictor)
  ) %>%
  round_export_cols(c("Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "AICc"))

summary_export <- model_run$summary %>%
  mutate(
    Outcome = format_export_outcome(Outcome),
    Predictors_Final = label_catchment_predictor_list(Predictors_Final)
  ) %>%
  round_export_cols(
    c(
      "R2", "R2_adj", "RMSE", "AIC", "AICc",
      "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
      "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model"
    )
  )

diagnostics_export <- model_run$diagnostics %>%
  mutate(
    Outcome = format_export_outcome(Outcome),
    Predictors_Final = label_catchment_predictor_list(Predictors_Final)
  ) %>%
  round_export_cols(c("shapiro_W", "shapiro_p", "ncv_chisq", "ncv_p"))

aicc_lt2_export <- model_run$model_selection %>%
  filter(is.finite(delta_AICc), delta_AICc <= 2) %>%
  mutate(
    Outcome = format_export_outcome(Outcome),
    Predictors_Final = label_catchment_predictor_list(Predictors_Final)
  ) %>%
  round_export_cols(
    c(
      "R2", "R2_adj", "RMSE", "AIC", "AICc",
      "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
      "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model",
      "delta_AICc"
    )
  ) %>%
  arrange(Outcome, delta_AICc, AICc, Candidate_Set)

write.csv(results_export, file.path(output_dir, paste0(file_prefix, "_results.csv")), row.names = FALSE)
write.csv(summary_export, file.path(output_dir, paste0(file_prefix, "_summary.csv")), row.names = FALSE)
write.csv(diagnostics_export, file.path(output_dir, paste0(file_prefix, "_diagnostics.csv")), row.names = FALSE)
write.csv(aicc_lt2_export, file.path(output_dir, paste0(file_prefix, "_aicc_lt2.csv")), row.names = FALSE)

if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  dir.create(OUT_TABLES_MLR_DIR, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    aicc_lt2_export,
    file.path(OUT_TABLES_MLR_DIR, "catchment_char_storage_mlr_aicc_lt2.csv"),
    row.names = FALSE
  )
  write.csv(
    diagnostics_export,
    file.path(OUT_TABLES_MLR_DIR, "catchment_char_storage_mlr_diagnostics.csv"),
    row.names = FALSE
  )
}
