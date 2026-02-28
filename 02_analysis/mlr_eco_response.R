# pooled eco response mlr models

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(MASS)
  library(car)
  library(tidyr)
})

rm(list = ls())
source("config.R")

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
    return(list(rmse_loocv = NA_real_, rmse_loocv_mean_runs = NA_real_, r2_loocv = NA_real_))
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
    return(list(rmse_loocv = NA_real_, rmse_loocv_mean_runs = NA_real_, r2_loocv = NA_real_))
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

iterative_vif <- function(model_obj, response, model_df, threshold) {
  lm_current <- model_obj

  repeat {
    retained <- setdiff(names(coef(lm_current)), "(Intercept)")
    non_mandatory <- setdiff(retained, mandatory_predictors)
    if (length(non_mandatory) == 0 || length(retained) <= 1) {
      break
    }

    vif_vals <- tryCatch(vif(lm_current), error = function(e) NULL)
    if (is.null(vif_vals)) {
      break
    }

    if (max(vif_vals, na.rm = TRUE) <= threshold) {
      break
    }

    ordered <- names(sort(vif_vals, decreasing = TRUE))
    drop_var <- ordered[ordered %in% non_mandatory][1]
    if (is.na(drop_var) || !nzchar(drop_var)) {
      break
    }

    retained_next <- setdiff(retained, drop_var)
    lm_current <- lm(
      as.formula(paste(response, "~", paste(retained_next, collapse = " + "))),
      data = model_df
    )
  }

  lm_current
}

compute_residual_diagnostics <- function(model_obj) {
  pull_scalar <- function(obj, key) {
    val <- tryCatch(obj[[key]], error = function(e) NULL)
    num <- suppressWarnings(as.numeric(val))
    if (length(num) >= 1 && is.finite(num[1])) return(num[1])
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

  ncv <- tryCatch(car::ncvTest(model_obj), error = function(e) NULL)
  ncv_chisq <- if (!is.null(ncv)) pull_scalar(ncv, "Chisquare") else NA_real_
  if (!is.finite(ncv_chisq) && !is.null(ncv)) ncv_chisq <- pull_scalar(ncv, "ChiSquare")
  ncv_p <- if (!is.null(ncv)) pull_scalar(ncv, "p") else NA_real_

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

  fit <- iterative_vif(fit, response, model_df, VIF_THRESHOLD)

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
    RMSE_LOOCV = loocv$rmse_loocv,
    RMSE_LOOCV_MEAN_RUNS = loocv$rmse_loocv_mean_runs,
    R2_LOOCV = loocv$r2_loocv,
    delta_RMSE_LOOCV_minus_model = loocv$rmse_loocv - sqrt(mean(residuals(fit)^2, na.rm = TRUE)),
    delta_RMSE_LOOCV_mean_runs_minus_model = loocv$rmse_loocv_mean_runs - sqrt(mean(residuals(fit)^2, na.rm = TRUE))
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

to_export_response <- function(x) {
  out <- as.character(x)
  out <- gsub("_", "", out, fixed = TRUE)
  out
}

to_export_predictor <- function(x) {
  as.character(x)
}

round_export_cols <- function(df, cols, digits = 3) {
  keep <- intersect(cols, names(df))
  if (length(keep) == 0) return(df)
  df %>% mutate(across(all_of(keep), ~ signif(.x, digits)))
}

results_export <- results %>%
  mutate(
    Response = to_export_response(Response),
    Predictor = to_export_predictor(Predictor)
  ) %>%
  round_export_cols(c("Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "AICc"))

summary_export <- summary_tbl %>%
  mutate(
    Response = to_export_response(Response),
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
    Response = to_export_response(Response),
    Predictors_Final = as.character(Predictors_Final)
  ) %>%
  round_export_cols(c("shapiro_W", "shapiro_p", "ncv_chisq", "ncv_p"))

pred_obs_export <- pred_obs %>%
  mutate(Response = to_export_response(Response)) %>%
  round_export_cols(c("Observed", "Predicted", "Residual"), digits = 4)

aicc_lt2_export <- selection %>%
  filter(is.finite(delta_AICc), delta_AICc <= 2) %>%
  mutate(
    Response = to_export_response(Response),
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
