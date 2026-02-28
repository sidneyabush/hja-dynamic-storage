# pooled eco response mlr models

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(MASS)
  library(car)
})

rm(list = ls())
source("config.R")
source("helpers/mlr_utils.R")

output_dir <- OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
file_prefix <- "storage_eco_response_mlr"

MODEL_MIN_N <- 20
LOW_N_THRESHOLD_ECO <- 25
VIF_THRESHOLD <- 10

annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  stop("Missing required file: ", annual_file)
}

data_all <- read_csv(annual_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

if (!("T_7DMax" %in% names(data_all)) && ("max_temp_7d_C" %in% names(data_all))) {
  data_all <- data_all %>% mutate(T_7DMax = max_temp_7d_C)
}
if (!("Q_7Q5" %in% names(data_all)) && ("q5_7d_mm_d" %in% names(data_all))) {
  data_all <- data_all %>% mutate(Q_7Q5 = q5_7d_mm_d)
}
if (!("Q_7Q5" %in% names(data_all)) && ("min_Q_7d_mm_d" %in% names(data_all))) {
  data_all <- data_all %>% mutate(Q_7Q5 = min_Q_7d_mm_d)
}

if (!("Pws" %in% names(data_all))) {
  data_all$Pws <- NA_real_
}
if ("P_WetSeason" %in% names(data_all)) {
  data_all$Pws <- dplyr::coalesce(data_all$Pws, data_all$P_WetSeason)
}
if ("precip_nov_may_mm" %in% names(data_all)) {
  data_all$Pws <- dplyr::coalesce(data_all$Pws, data_all$precip_nov_may_mm)
}
if ("P_NovJan" %in% names(data_all)) {
  data_all$Pws <- dplyr::coalesce(data_all$Pws, data_all$P_NovJan)
}
if ("precip_nov_jan_mm" %in% names(data_all)) {
  data_all$Pws <- dplyr::coalesce(data_all$Pws, data_all$precip_nov_jan_mm)
}

response_vars <- c("Q_7Q5", "T_7DMax")
missing_responses <- setdiff(response_vars, names(data_all))
if (length(missing_responses) > 0) {
  stop("Missing eco response variable(s): ", paste(missing_responses, collapse = ", "))
}

mandatory_predictors <- c("Pws")
storage_predictors <- c("RBI", "RCS", "FDC", "SD", "WB", "CHS")
storage_predictors <- storage_predictors[storage_predictors %in% names(data_all)]
if (!("Pws" %in% names(data_all))) {
  stop("Missing required predictor: Pws")
}

candidate_sets <- list(c("Pws"))
if (length(storage_predictors) > 0) {
  storage_combos <- unlist(
    lapply(seq_along(storage_predictors), function(k) combn(storage_predictors, k, simplify = FALSE)),
    recursive = FALSE
  )
  candidate_sets <- c(candidate_sets, lapply(storage_combos, function(x) c("Pws", x)))
}

fit_candidate <- function(df_in, response, predictors) {
  predictors <- unique(predictors[predictors %in% names(df_in)])
  if (!("Pws" %in% predictors)) {
    return(NULL)
  }

  model_df <- df_in %>%
    dplyr::select(any_of(c("site", "year", response, predictors))) %>%
    na.omit()

  if (nrow(model_df) < MODEL_MIN_N) {
    return(NULL)
  }

  pred_sd <- sapply(predictors, function(p) suppressWarnings(sd(model_df[[p]], na.rm = TRUE)))
  predictors <- predictors[is.finite(pred_sd) & pred_sd > 0]
  if (!("Pws" %in% predictors)) {
    return(NULL)
  }

  model_df <- model_df %>%
    mutate(across(all_of(predictors), ~ as.numeric(scale(.x))))

  model_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  fit <- tryCatch(lm(model_formula, data = model_df), error = function(e) NULL)
  if (is.null(fit)) {
    return(NULL)
  }

  fit <- apply_vif_filter(
    model_obj = fit,
    outcome = response,
    model_df = model_df,
    threshold = VIF_THRESHOLD,
    protected = mandatory_predictors
  )
  if (is.null(fit)) {
    return(NULL)
  }

  retained <- setdiff(names(coef(fit)), "(Intercept)")
  if (!("Pws" %in% retained)) {
    retained <- unique(c("Pws", retained))
    fit <- tryCatch(
      lm(as.formula(paste(response, "~", paste(retained, collapse = " + "))), data = model_df),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
  }

  retained <- setdiff(names(coef(fit)), "(Intercept)")
  if (length(retained) == 0) {
    return(NULL)
  }

  fit_sum <- summary(fit)
  aic_val <- AIC(fit)
  aicc_val <- calc_aicc(fit, nrow(model_df))

  fstat <- suppressWarnings(as.numeric(fit_sum$fstatistic))
  model_p_global <- if (!is.null(fstat) && length(fstat) >= 3 && all(is.finite(fstat[1:3]))) {
    pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  } else {
    NA_real_
  }

  vif_df <- if (length(retained) > 1) {
    vif_vals <- tryCatch(vif(fit), error = function(e) NULL)
    if (is.null(vif_vals)) tibble(Predictor = retained, VIF = NA_real_)
    else tibble(Predictor = names(vif_vals), VIF = as.numeric(vif_vals))
  } else {
    tibble(Predictor = retained, VIF = NA_real_)
  }

  coef_df <- as.data.frame(fit_sum$coefficients)
  coef_df$Predictor <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df <- coef_df %>% filter(Predictor != "(Intercept)")

  beta_df_data <- model_df
  beta_df_data[[response]] <- as.numeric(scale(beta_df_data[[response]]))
  fit_beta <- lm(as.formula(paste(response, "~", paste(retained, collapse = " + "))), data = beta_df_data)
  beta_df <- tibble(
    Predictor = names(coef(fit_beta))[-1],
    Beta_Std = as.numeric(coef(fit_beta)[-1])
  )

  coef_out <- coef_df %>%
    transmute(
      Predictor,
      p_value = `Pr(>|t|)`
    ) %>%
    left_join(vif_df, by = "Predictor") %>%
    left_join(beta_df, by = "Predictor") %>%
    mutate(
      Site = "all_sites",
      Response = response,
      R2 = fit_sum$r.squared,
      R2_adj = fit_sum$adj.r.squared,
      RMSE = sqrt(mean(residuals(fit)^2, na.rm = TRUE)),
      AIC = aic_val,
      AICc = aicc_val,
      n = nrow(model_df)
    ) %>%
    dplyr::select(Site, Response, Predictor, Beta_Std, p_value, VIF, R2, R2_adj, RMSE, AIC, AICc, n)

  loocv <- calc_loocv_stats(formula(fit), model_df)

  summary_out <- tibble(
    Site = "all_sites",
    Response = response,
    Predictors_Final = paste(retained, collapse = "; "),
    R2 = fit_sum$r.squared,
    R2_adj = fit_sum$adj.r.squared,
    model_p_global = model_p_global,
    RMSE = sqrt(mean(residuals(fit)^2, na.rm = TRUE)),
    AIC = aic_val,
    AICc = aicc_val,
    n = nrow(model_df),
    RMSE_LOOCV = loocv$rmse,
    RMSE_LOOCV_MEAN_RUNS = loocv$mae,
    R2_LOOCV = loocv$r2,
    delta_RMSE_LOOCV_minus_model = loocv$rmse - sqrt(mean(residuals(fit)^2, na.rm = TRUE)),
    delta_RMSE_LOOCV_mean_runs_minus_model = loocv$mae - sqrt(mean(residuals(fit)^2, na.rm = TRUE))
  )

  pred_obs_out <- model_df %>%
    mutate(
      Site = if ("site" %in% names(model_df)) as.character(site) else "all_sites",
      Response = response,
      Predicted = as.numeric(predict(fit, newdata = model_df)),
      Residual = .data[[response]] - Predicted
    ) %>%
    transmute(
      Site,
      Response,
      year = if ("year" %in% names(model_df)) as.numeric(year) else NA_real_,
      Observed = .data[[response]],
      Predicted,
      Residual
    )

  diagnostics_out <- compute_residual_diagnostics(fit) %>%
    mutate(
      Site = "all_sites",
      Response = response,
      Predictors_Final = paste(retained, collapse = "; "),
      n = nrow(model_df)
    ) %>%
    dplyr::select(
      Site, Response, Predictors_Final, n, n_residuals,
      shapiro_W, shapiro_p, ncv_chisq, ncv_p,
      normality_pass_p05, homoscedasticity_pass_p05
    )

  list(
    coefficients = coef_out,
    summary = summary_out,
    pred_obs = pred_obs_out,
    diagnostics = diagnostics_out
  )
}

results_list <- list()
summary_list <- list()
pred_obs_list <- list()
diagnostics_list <- list()
selection_list <- list()

for (response in response_vars) {
  fits <- list()
  summaries <- list()

  for (i in seq_along(candidate_sets)) {
    fit_obj <- fit_candidate(data_all, response, candidate_sets[[i]])
    if (is.null(fit_obj)) next

    fit_obj$summary <- fit_obj$summary %>% mutate(Candidate_Set = i)
    fits[[length(fits) + 1]] <- fit_obj
    summaries[[length(summaries) + 1]] <- fit_obj$summary
  }

  if (length(fits) == 0) next

  selection_tbl <- bind_rows(summaries) %>%
    mutate(delta_AICc = AICc - min(AICc, na.rm = TRUE)) %>%
    arrange(AICc)
  selection_list[[response]] <- selection_tbl

  best_idx <- which.min(vapply(fits, function(x) x$summary$AICc[1], numeric(1)))
  best_fit <- fits[[best_idx]]

  results_list[[response]] <- best_fit$coefficients
  summary_list[[response]] <- best_fit$summary %>% dplyr::select(-Candidate_Set)
  pred_obs_list[[response]] <- best_fit$pred_obs
  diagnostics_list[[response]] <- best_fit$diagnostics
}

if (length(summary_list) == 0) {
  stop("No pooled eco models could be fit")
}

results <- bind_rows(results_list)
summary_tbl <- bind_rows(summary_list) %>%
  mutate(
    low_n_flag = n < LOW_N_THRESHOLD_ECO,
    confidence_note = ifelse(low_n_flag, "lower confidence (small n)", "standard confidence")
  )
pred_obs <- bind_rows(pred_obs_list)
diagnostics <- bind_rows(diagnostics_list) %>%
  left_join(
    summary_tbl %>% dplyr::select(Site, Response, low_n_flag, confidence_note),
    by = c("Site", "Response")
  )
selection <- bind_rows(selection_list)

results_export <- results %>%
  mutate(
    Response = format_export_outcome(Response, strip_mean = FALSE, drop_underscores = TRUE),
    Predictor = as.character(Predictor)
  ) %>%
  round_export_cols(c("Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "AICc"))

summary_export <- summary_tbl %>%
  mutate(
    Response = format_export_outcome(Response, strip_mean = FALSE, drop_underscores = TRUE),
    Predictors_Final = as.character(Predictors_Final)
  ) %>%
  round_export_cols(
    c(
      "R2", "R2_adj", "model_p_global", "RMSE", "AIC", "AICc",
      "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
      "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model"
    )
  )

diagnostics_export <- diagnostics %>%
  mutate(
    Response = format_export_outcome(Response, strip_mean = FALSE, drop_underscores = TRUE),
    Predictors_Final = as.character(Predictors_Final)
  ) %>%
  round_export_cols(c("shapiro_W", "shapiro_p", "ncv_chisq", "ncv_p"))

pred_obs_export <- pred_obs %>%
  mutate(Response = format_export_outcome(Response, strip_mean = FALSE, drop_underscores = TRUE)) %>%
  round_export_cols(c("Observed", "Predicted", "Residual"), digits = 4)

aicc_lt2_export <- selection %>%
  filter(is.finite(delta_AICc), delta_AICc <= 2) %>%
  mutate(
    Response = format_export_outcome(Response, strip_mean = FALSE, drop_underscores = TRUE),
    Predictors_Final = as.character(Predictors_Final)
  ) %>%
  round_export_cols(
    c(
      "R2", "R2_adj", "model_p_global", "RMSE", "AIC", "AICc",
      "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
      "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model",
      "delta_AICc"
    )
  )

write.csv(results_export, file.path(output_dir, paste0(file_prefix, "_results.csv")), row.names = FALSE)
write.csv(summary_export, file.path(output_dir, paste0(file_prefix, "_summary.csv")), row.names = FALSE)
write.csv(diagnostics_export, file.path(output_dir, paste0(file_prefix, "_diagnostics.csv")), row.names = FALSE)
write.csv(pred_obs_export, file.path(output_dir, paste0(file_prefix, "_predicted_observed.csv")), row.names = FALSE)
write.csv(aicc_lt2_export, file.path(output_dir, paste0(file_prefix, "_aicc_lt2.csv")), row.names = FALSE)

if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  dir.create(OUT_TABLES_MLR_DIR, recursive = TRUE, showWarnings = FALSE)
  write.csv(aicc_lt2_export, file.path(OUT_TABLES_MLR_DIR, paste0(file_prefix, "_aicc_lt2.csv")), row.names = FALSE)
  write.csv(diagnostics_export, file.path(OUT_TABLES_MLR_DIR, paste0(file_prefix, "_diagnostics.csv")), row.names = FALSE)
}
