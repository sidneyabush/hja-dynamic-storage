# catchment and ecological response MLR models

# inputs:
# outputs/master/master_site.csv
# outputs/master/master_annual.csv
# outputs/metrics/mobile/isotope_metrics_site.csv

# outputs:
# outputs/models/catchment_char_storage_mlr/*.csv
# outputs/models/storage_eco_response_mlr/*.csv

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(dplyr, readr, MASS, car, cran_repo = "https://cloud.r-project.org")

rm(list = ls())
source("config.R")

# catchment characteristic models
file_prefix <- "catchment_char_storage_mlr"
output_dir <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

MIN_N <- 5
LOW_N_THRESHOLD_CATCH <- 10
VIF_THRESHOLD <- 10

site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
site_df <- read_csv(site_file, show_col_types = FALSE) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

# use the combined site level MTT values from isotope_metrics_site
isotope_site_mean_file <- file.path(OUT_MET_MOBILE_DIR, "isotope_metrics_site.csv")
mtt_site_mean <- read_csv(isotope_site_mean_file, show_col_types = FALSE) %>%
  mutate(
    site = if ("site" %in% names(.)) site else SITECODE,
    site = standardize_site_code(site),
    MTT_site_combined = suppressWarnings(as.numeric(MTT))
  ) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC, !is.na(site), site != "") %>%
  group_by(site) %>%
  summarise(
    MTT_site_combined = ifelse(any(is.finite(MTT_site_combined)), mean(MTT_site_combined, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )

site_df <- site_df %>%
  mutate(site = standardize_site_code(site)) %>%
  left_join(mtt_site_mean, by = "site") %>%
  mutate(
    MTT = dplyr::coalesce(MTT_site_combined, suppressWarnings(as.numeric(MTT)))
  ) %>%
  dplyr::select(-MTT_site_combined)

if (!("basin_slope" %in% names(site_df)) && ("Slope_mean" %in% names(site_df))) {
  site_df <- site_df %>% mutate(basin_slope = Slope_mean)
}

outcome_map <- c(
  "RBI" = "RBI_mean",
  "RCS" = "RCS_mean",
  "FDC" = "FDC_mean",
  "SD" = "SD_mean",
  "WB" = "WB_mean",
  "BF" = "BF_mean",
  "DR" = "DR",
  "Fyw" = "Fyw",
  "MTT" = "MTT"
)

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

candidate_sets <- unlist(
  lapply(seq_along(predictor_vars), function(k) {
    combn(predictor_vars, k, simplify = FALSE)
  }),
  recursive = FALSE
)

# fit one catchment predictor set
fit_candidate <- function(df_in, outcome, predictors) {
  model_input <- df_in

  # remove any sites excluded from baseflow modeling
  if (identical(outcome, "BF_mean") && "site" %in% names(model_input)) {
    model_input <- model_input %>% filter(!(site %in% BF_EXCLUDE_SITES))
  }

  predictors <- predictors[predictors %in% names(model_input)]

  # skip candidates missing the response or all predictors
  if (!(outcome %in% names(model_input)) || length(predictors) == 0) {
    return(NULL)
  }

  model_df <- model_input %>%
    dplyr::select(all_of(c(outcome, predictors))) %>%
    na.omit()

  # require at least five complete sites for catchment models
  if (nrow(model_df) < MIN_N) {
    return(NULL)
  }

  # drop predictors that are constant across complete cases
  pred_sd <- sapply(predictors, function(x) suppressWarnings(sd(model_df[[x]], na.rm = TRUE)))
  predictors <- predictors[is.finite(pred_sd) & pred_sd > 0]
  if (length(predictors) == 0) {
    return(NULL)
  }

  model_df <- model_df %>%
    mutate(across(all_of(predictors), ~ as.numeric(scale(.x))))

  fit <- tryCatch(
    lm(as.formula(paste(outcome, "~", paste(predictors, collapse = " + "))), data = model_df),
    error = function(e) NULL
  )

  # skip candidates that cannot be fit by lm
  if (is.null(fit)) {
    return(NULL)
  }

  # let AIC remove unhelpful predictors before the shared filters are applied
  fit <- tryCatch(stepAIC(fit, direction = "backward", trace = 0), error = function(e) fit)

  fit <- apply_vif_filter(fit, outcome, model_df, VIF_THRESHOLD)

  # skip models that lose all predictors during VIF filtering
  if (is.null(fit)) {
    return(NULL)
  }

  fit <- apply_correlated_predictor_rules_catchment(fit, outcome, model_df)
  retained <- setdiff(names(coef(fit)), "(Intercept)")

  # skip intercept-only models
  if (length(retained) == 0) {
    return(NULL)
  }

  # collect fit statistics for the selected candidate
  fit_sum <- summary(fit)
  rmse <- sqrt(mean(residuals(fit)^2, na.rm = TRUE))
  aic_val <- AIC(fit)
  aicc_val <- calc_aicc(fit, nrow(model_df))
  loocv <- calc_loocv_stats(formula(fit), model_df)

  fstat <- suppressWarnings(as.numeric(fit_sum$fstatistic))
  model_p <- if (!is.null(fstat) && length(fstat) >= 3 && all(is.finite(fstat[1:3]))) {
    pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  } else {
    NA_real_
  }

  coef_df <- as.data.frame(fit_sum$coefficients)
  coef_df$Predictor <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df <- coef_df %>% filter(Predictor != "(Intercept)")

  vif_tbl <- if (length(retained) > 1) {
    vif_vals <- tryCatch(car::vif(fit), error = function(e) NULL)
    if (is.null(vif_vals)) tibble(Predictor = retained, VIF = NA_real_)
    else tibble(Predictor = names(vif_vals), VIF = as.numeric(vif_vals))
  } else {
    tibble(Predictor = retained, VIF = NA_real_)
  }

  beta_df <- model_df
  beta_df[[outcome]] <- as.numeric(scale(beta_df[[outcome]]))
  fit_beta <- lm(as.formula(paste(outcome, "~", paste(retained, collapse = " + "))), data = beta_df)
  beta_tbl <- tibble(
    Predictor = names(coef(fit_beta))[-1],
    Beta_Std = as.numeric(coef(fit_beta)[-1])
  )

  coef_out <- coef_df %>%
    transmute(
      Outcome = ifelse(outcome == "DR", "DR_mean", outcome),
      Predictor = Predictor,
      Beta_Std = NA_real_,
      p_value = `Pr(>|t|)`,
      VIF = NA_real_,
      R2 = fit_sum$r.squared,
      R2_adj = fit_sum$adj.r.squared,
      RMSE = rmse,
      AIC = aic_val,
      AICc = aicc_val
    ) %>%
    left_join(beta_tbl, by = "Predictor", suffix = c("", "_calc")) %>%
    mutate(Beta_Std = dplyr::coalesce(Beta_Std_calc, Beta_Std)) %>%
    dplyr::select(-Beta_Std_calc) %>%
    left_join(vif_tbl, by = "Predictor", suffix = c("", "_calc")) %>%
    mutate(VIF = dplyr::coalesce(VIF_calc, VIF)) %>%
    dplyr::select(-VIF_calc)

  summary_out <- tibble(
    Outcome = ifelse(outcome == "DR", "DR_mean", outcome),
    N = nrow(model_df),
    Predictors_Final = paste(retained, collapse = ", "),
    R2 = fit_sum$r.squared,
    R2_adj = fit_sum$adj.r.squared,
    model_p_global = model_p,
    RMSE = rmse,
    AIC = aic_val,
    AICc = aicc_val,
    RMSE_LOOCV = loocv$rmse,
    RMSE_LOOCV_MEAN_RUNS = loocv$mae,
    R2_LOOCV = loocv$r2,
    delta_RMSE_LOOCV_minus_model = loocv$rmse - rmse,
    delta_RMSE_LOOCV_mean_runs_minus_model = loocv$mae - rmse
  )

  diagnostics_out <- compute_residual_diagnostics(fit) %>%
    mutate(
      Outcome = ifelse(outcome == "DR", "DR_mean", outcome),
      Predictors_Final = paste(retained, collapse = ", "),
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

all_results <- list()
all_summary <- list()
all_diag <- list()
all_select <- list()

# each storage metric gets the same candidate-model search
# the selected model is the lowest AICc fit left after VIF and correlated-predictor filters
for (metric in STORAGE_METRIC_ORDER) {
  outcome <- unname(outcome_map[metric])
  if (!is.character(outcome) || length(outcome) == 0 || !(outcome %in% names(site_df))) {
    next
  }

  candidate_fits <- list()
  candidate_rows <- list()

  # keep candidate models that pass the data and VIF filters
  for (i in seq_along(candidate_sets)) {
    fit_obj <- fit_candidate(site_df, outcome, candidate_sets[[i]])
    if (is.null(fit_obj)) {
      next
    }

    row_i <- fit_obj$summary %>%
      mutate(Candidate_Set = i)
    candidate_fits[[length(candidate_fits) + 1]] <- fit_obj
    candidate_rows[[length(candidate_rows) + 1]] <- row_i
  }

  if (length(candidate_fits) == 0) {
    next
  }

  # rank retained candidate models by AICc
  selection_tbl <- bind_rows(candidate_rows) %>%
    mutate(delta_AICc = AICc - min(AICc, na.rm = TRUE)) %>%
    arrange(AICc)
  all_select[[metric]] <- selection_tbl

  best_idx <- which.min(vapply(candidate_fits, function(x) x$summary$AICc[1], numeric(1)))
  best_fit <- candidate_fits[[best_idx]]

  all_results[[metric]] <- best_fit$coefficients
  all_summary[[metric]] <- best_fit$summary
  all_diag[[metric]] <- best_fit$diagnostics
}

# stop if no catchment storage models can be fit
if (length(all_summary) == 0) {
  stop("No catchment characteristic models could be fit")
}

results_tbl <- bind_rows(all_results)
summary_tbl <- bind_rows(all_summary) %>%
  mutate(
    low_n_flag = N < LOW_N_THRESHOLD_CATCH,
    confidence_note = ifelse(low_n_flag, "lower confidence (small n)", "standard confidence")
  )
diag_tbl <- bind_rows(all_diag) %>%
  left_join(
    summary_tbl %>% dplyr::select(Outcome, low_n_flag, confidence_note),
    by = "Outcome"
  )
select_tbl <- bind_rows(all_select)

# format labels and round numeric columns before writing catchment model files
results_export <- results_tbl %>%
  mutate(
    Outcome = format_export_outcome(Outcome),
    Predictor = label_catchment_predictor(Predictor)
  ) %>%
  round_export_cols(c("Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "AICc"))

summary_export <- summary_tbl %>%
  mutate(
    Outcome = format_export_outcome(Outcome),
    Predictors_Final = label_catchment_predictor_list(Predictors_Final)
  ) %>%
  round_export_cols(
    c(
      "R2", "R2_adj", "model_p_global", "RMSE", "AIC", "AICc",
      "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
      "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model"
    )
  )

diag_export <- diag_tbl %>%
  mutate(
    Outcome = format_export_outcome(Outcome),
    Predictors_Final = label_catchment_predictor_list(Predictors_Final)
  ) %>%
  round_export_cols(c("shapiro_W", "shapiro_p", "ncv_chisq", "ncv_p"))

aicc_lt2_export <- select_tbl %>%
  filter(is.finite(delta_AICc), delta_AICc <= 2) %>%
  mutate(
    Outcome = format_export_outcome(Outcome),
    Predictors_Final = label_catchment_predictor_list(Predictors_Final)
  ) %>%
  round_export_cols(
    c(
      "R2", "R2_adj", "model_p_global", "RMSE", "AIC", "AICc",
      "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
      "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model",
      "delta_AICc"
    )
  ) %>%
  arrange(Outcome, delta_AICc, AICc, Candidate_Set)

# write catchment model outputs for Figure 4 and Tables S8 to S10
write.csv(
  results_export,
  file.path(output_dir, paste0(file_prefix, "_results.csv")),
  row.names = FALSE
)
write.csv(
  summary_export,
  file.path(output_dir, paste0(file_prefix, "_summary.csv")),
  row.names = FALSE
)
write.csv(
  diag_export,
  file.path(output_dir, paste0(file_prefix, "_diagnostics.csv")),
  row.names = FALSE
)
write.csv(
  aicc_lt2_export,
  file.path(output_dir, paste0(file_prefix, "_aicc_lt2.csv")),
  row.names = FALSE
)

# ecological response models
output_dir <- OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
file_prefix <- "storage_eco_response_mlr"

MODEL_MIN_N <- 20
LOW_N_THRESHOLD_ECO <- 25
VIF_THRESHOLD <- 10

annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
data_all <- read_csv(annual_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

# use site level isotope means (DR, Fyw, MTT) as eco model predictors
isotope_site_mean_file <- file.path(OUT_MET_MOBILE_DIR, "isotope_metrics_site.csv")
isotope_site_mean <- read_csv(isotope_site_mean_file, show_col_types = FALSE) %>%
  mutate(
    site = if ("site" %in% names(.)) site else SITECODE,
    site = standardize_site_code(site),
    DR_site_mean = suppressWarnings(as.numeric(DR)),
    Fyw_site_mean = suppressWarnings(as.numeric(Fyw)),
    MTT_site_mean = suppressWarnings(as.numeric(MTT))
  ) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC, !is.na(site), site != "") %>%
  group_by(site) %>%
  summarise(
    DR_site_mean = ifelse(any(is.finite(DR_site_mean)), mean(DR_site_mean, na.rm = TRUE), NA_real_),
    Fyw_site_mean = ifelse(any(is.finite(Fyw_site_mean)), mean(Fyw_site_mean, na.rm = TRUE), NA_real_),
    MTT_site_mean = ifelse(any(is.finite(MTT_site_mean)), mean(MTT_site_mean, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )

data_all <- data_all %>%
  mutate(site = standardize_site_code(site)) %>%
  left_join(isotope_site_mean, by = "site") %>%
  mutate(
    DR = DR_site_mean,
    Fyw = Fyw_site_mean,
    MTT = MTT_site_mean
  ) %>%
  dplyr::select(-DR_site_mean, -Fyw_site_mean, -MTT_site_mean)

response_vars <- c("Q_7Q5", "T_7DMax")
missing_responses <- setdiff(response_vars, names(data_all))

# stop if the ecological response table changed structure
if (length(missing_responses) > 0) {
  stop("Missing eco response variable(s): ", paste(missing_responses, collapse = ", "))
}

mandatory_predictors <- c("Pws")
storage_predictors <- c("RBI", "RCS", "FDC", "SD", "WB", "BF", "DR", "Fyw", "MTT")
storage_predictors <- storage_predictors[storage_predictors %in% names(data_all)]

# wet-season precipitation stays in every ecological response model
if (!("Pws" %in% names(data_all))) {
  stop("Missing predictor: Pws")
}

candidate_sets <- list(c("Pws"))

# add every storage predictor combination while keeping Pws in the model
if (length(storage_predictors) > 0) {
  storage_combos <- unlist(
    lapply(seq_along(storage_predictors), function(k) combn(storage_predictors, k, simplify = FALSE)),
    recursive = FALSE
  )
  candidate_sets <- c(candidate_sets, lapply(storage_combos, function(x) c("Pws", x)))
}

# fit one ecological response predictor set
fit_candidate <- function(df_in, response, predictors) {
  predictors <- unique(predictors[predictors %in% names(df_in)])

  # every ecological-response model includes wet-season precipitation
  if (!("Pws" %in% predictors)) {
    return(NULL)
  }

  model_df <- df_in %>%
    dplyr::select(any_of(c("site", "year", response, predictors))) %>%
    na.omit()

  # require at least 20 complete site years for ecological response models
  if (nrow(model_df) < MODEL_MIN_N) {
    return(NULL)
  }

  # drop predictors that are constant across complete cases
  pred_sd <- sapply(predictors, function(p) suppressWarnings(sd(model_df[[p]], na.rm = TRUE)))
  predictors <- predictors[is.finite(pred_sd) & pred_sd > 0]

  # skip candidates where Pws became unusable
  if (!("Pws" %in% predictors)) {
    return(NULL)
  }

  model_df <- model_df %>%
    mutate(across(all_of(predictors), ~ as.numeric(scale(.x))))

  model_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  fit <- tryCatch(lm(model_formula, data = model_df), error = function(e) NULL)

  # skip candidates that cannot be fit by lm
  if (is.null(fit)) {
    return(NULL)
  }

  # let AIC remove unhelpful predictors before the shared filters are applied
  fit <- tryCatch(stepAIC(fit, direction = "backward", trace = 0), error = function(e) fit)

  fit <- apply_vif_filter(
    model_obj = fit,
    outcome = response,
    model_df = model_df,
    threshold = VIF_THRESHOLD,
    protected = mandatory_predictors
  )

  # skip models that lose all predictors during VIF filtering
  if (is.null(fit)) {
    return(NULL)
  }

  retained <- setdiff(names(coef(fit)), "(Intercept)")

  # refit if stepwise selection removed protected Pws
  if (!("Pws" %in% retained)) {
    retained <- unique(c("Pws", retained))
    fit <- tryCatch(
      lm(as.formula(paste(response, "~", paste(retained, collapse = " + "))), data = model_df),
      error = function(e) NULL
    )

    if (is.null(fit)) return(NULL)
  }

  retained <- setdiff(names(coef(fit)), "(Intercept)")

  # skip intercept-only models
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

  # report VIF values only when more than one predictor is retained
  vif_df <- if (length(retained) > 1) {
    vif_vals <- tryCatch(vif(fit), error = function(e) NULL)
    if (is.null(vif_vals)) tibble(Predictor = retained, VIF = NA_real_)
    else tibble(Predictor = names(vif_vals), VIF = as.numeric(vif_vals))
  } else {
    tibble(Predictor = retained, VIF = NA_real_)
  }

  # standardize the response to calculate comparable beta values
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

  # summarize model fit and leave one out prediction error
  summary_out <- tibble(
    Site = "all_sites",
    Response = response,
    Predictors_Final = paste(retained, collapse = ", "),
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

  # store observed, predicted, and residual values for Figure 6
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

  # keep residual diagnostics with the same labels as the summary table
  diagnostics_out <- compute_residual_diagnostics(fit) %>%
    mutate(
      Site = "all_sites",
      Response = response,
      Predictors_Final = paste(retained, collapse = ", "),
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

# Q7Q5 and T7DMax use the same candidate-model search
# the selected model is the lowest AICc fit left after the shared filters
for (response in response_vars) {
  fits <- list()
  summaries <- list()

  # keep candidate models that pass the data and VIF filters
  for (i in seq_along(candidate_sets)) {
    fit_obj <- fit_candidate(data_all, response, candidate_sets[[i]])
    if (is.null(fit_obj)) next

    fit_obj$summary <- fit_obj$summary %>% mutate(Candidate_Set = i)
    fits[[length(fits) + 1]] <- fit_obj
    summaries[[length(summaries) + 1]] <- fit_obj$summary
  }

  if (length(fits) == 0) next

  # rank retained candidate models by AICc
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

# stop if no ecological response models can be fit
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

# format labels and round numeric columns before writing ecological response files
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

# write ecological response model outputs for Figures 5 and 6 and Tables S10 to S12
write.csv(results_export, file.path(output_dir, paste0(file_prefix, "_results.csv")), row.names = FALSE)
write.csv(summary_export, file.path(output_dir, paste0(file_prefix, "_summary.csv")), row.names = FALSE)
write.csv(diagnostics_export, file.path(output_dir, paste0(file_prefix, "_diagnostics.csv")), row.names = FALSE)
write.csv(pred_obs_export, file.path(output_dir, paste0(file_prefix, "_predicted_observed.csv")), row.names = FALSE)
write.csv(aicc_lt2_export, file.path(output_dir, paste0(file_prefix, "_aicc_lt2.csv")), row.names = FALSE)
