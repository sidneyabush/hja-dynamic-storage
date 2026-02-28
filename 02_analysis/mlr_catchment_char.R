# catchment characteristics mlr

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(MASS)
  library(car)
})

rm(list = ls())
source("config.R")
source("helpers/mlr_utils.R")

file_prefix <- "catchment_char_storage_mlr"
output_dir <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

MIN_N <- 5
LOW_N_THRESHOLD_CATCH <- 10
VIF_THRESHOLD <- 10

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

fit_candidate <- function(df_in, outcome, predictors) {
  model_input <- df_in
  if (identical(outcome, "CHS_mean") && "site" %in% names(model_input)) {
    model_input <- model_input %>% filter(!(site %in% CHS_EXCLUDE_SITES))
  }

  predictors <- predictors[predictors %in% names(model_input)]
  if (!(outcome %in% names(model_input)) || length(predictors) == 0) {
    return(NULL)
  }

  model_df <- model_input %>%
    dplyr::select(all_of(c(outcome, predictors))) %>%
    na.omit()

  if (nrow(model_df) < MIN_N) {
    return(NULL)
  }

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
  if (is.null(fit)) {
    return(NULL)
  }
  fit <- tryCatch(stepAIC(fit, direction = "backward", trace = 0), error = function(e) fit)

  fit <- apply_vif_filter(fit, outcome, model_df, VIF_THRESHOLD)
  if (is.null(fit)) {
    return(NULL)
  }

  fit <- apply_correlated_predictor_rules_catchment(fit, outcome, model_df)
  retained <- setdiff(names(coef(fit)), "(Intercept)")
  if (length(retained) == 0) {
    return(NULL)
  }

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
    Predictors_Final = paste(retained, collapse = "; "),
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

all_results <- list()
all_summary <- list()
all_diag <- list()
all_select <- list()

for (metric in STORAGE_METRIC_ORDER) {
  outcome <- unname(outcome_map[metric])
  if (!is.character(outcome) || length(outcome) == 0 || !(outcome %in% names(site_df))) {
    next
  }

  candidate_fits <- list()
  candidate_rows <- list()

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

if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  dir.create(OUT_TABLES_MLR_DIR, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    aicc_lt2_export,
    file.path(OUT_TABLES_MLR_DIR, "catchment_char_storage_mlr_aicc_lt2.csv"),
    row.names = FALSE
  )
  write.csv(
    diag_export,
    file.path(OUT_TABLES_MLR_DIR, "catchment_char_storage_mlr_diagnostics.csv"),
    row.names = FALSE
  )
}
