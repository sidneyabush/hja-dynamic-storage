# Fit pooled (all-sites) stepwise-AIC linear models for annual eco response metrics.
# Inputs: No direct CSV file reads in this script.
# Author: Sidney Bush
# Date: 2026-02-24

library(dplyr)
library(readr)
library(MASS)
library(car)
library(tidyr)

rm(list = ls())

# Load project config
source("config.R")

output_dir <- OUT_MODELS_STORAGE_ECOVAR_MLR_DIR
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
file_prefix <- "storage_ecovar_mlr_all_sites"
ALL_SITES_ID <- "all_sites"

MODEL_MIN_N <- 20
VIF_THRESHOLD <- 10
LOW_N_THRESHOLD_ECO <- 25
USE_MUMIN_LOOCV <- TRUE
STRICT_USE_SCALED_PREDICTORS <- TRUE
STRICT_USE_ITERATIVE_VIF <- TRUE

env_vif_threshold <- suppressWarnings(as.numeric(Sys.getenv("HJA_MLR_VIF_THRESHOLD", unset = "")))
if (!is.na(env_vif_threshold) && is.finite(env_vif_threshold) && env_vif_threshold > 0) {
  VIF_THRESHOLD <- env_vif_threshold
}

if (USE_MUMIN_LOOCV && !requireNamespace("MuMIn", quietly = TRUE)) {
  stop("MuMIn is required for eco-model LOOCV. Install with: install.packages('MuMIn')")
}

master_dir <- file.path(OUTPUT_DIR, "master")
annual_file <- file.path(master_dir, MASTER_ANNUAL_FILE)

merged_data <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

if (!("T_7DMax" %in% names(merged_data)) && ("max_temp_7d_C" %in% names(merged_data))) {
  merged_data <- merged_data %>%
    mutate(T_7DMax = max_temp_7d_C)
}

if (!("Q_7Q5" %in% names(merged_data)) && ("q5_7d_mm_d" %in% names(merged_data))) {
  merged_data <- merged_data %>%
    mutate(Q_7Q5 = q5_7d_mm_d)
}
if (!("Q_7Q5" %in% names(merged_data)) && ("min_Q_7d_mm_d" %in% names(merged_data))) {
  merged_data <- merged_data %>%
    mutate(Q_7Q5 = min_Q_7d_mm_d)
}

if (!("T_Q7Q5" %in% names(merged_data)) && ("temp_during_q5_7d_C" %in% names(merged_data))) {
  merged_data <- merged_data %>%
    mutate(T_Q7Q5 = temp_during_q5_7d_C)
}
if (!("T_Q7Q5" %in% names(merged_data)) && ("temp_during_min_Q_7d_C" %in% names(merged_data))) {
  merged_data <- merged_data %>%
    mutate(T_Q7Q5 = temp_during_min_Q_7d_C)
}

if (!("T_Q7Q5" %in% names(merged_data)) && ("temp_at_min_Q_7d_C" %in% names(merged_data))) {
  merged_data <- merged_data %>%
    mutate(T_Q7Q5 = temp_at_min_Q_7d_C)
}

if (!("P_WetSeason" %in% names(merged_data))) {
  merged_data$P_WetSeason <- NA_real_
}
if ("precip_nov_may_mm" %in% names(merged_data)) {
  merged_data$P_WetSeason <- dplyr::coalesce(merged_data$P_WetSeason, merged_data$precip_nov_may_mm)
}
if ("P_NovJan" %in% names(merged_data)) {
  merged_data$P_WetSeason <- dplyr::coalesce(merged_data$P_WetSeason, merged_data$P_NovJan)
}
if ("precip_nov_jan_mm" %in% names(merged_data)) {
  merged_data$P_WetSeason <- dplyr::coalesce(merged_data$P_WetSeason, merged_data$precip_nov_jan_mm)
}

response_vars_required <- c("T_7DMax", "Q_7Q5", "T_Q7Q5")
missing_response_vars <- setdiff(response_vars_required, names(merged_data))
if (length(missing_response_vars) > 0) {
  stop(
    paste(
      "Missing required eco response variable(s):",
      paste(missing_response_vars, collapse = ", ")
    )
  )
}
response_vars <- response_vars_required

required_predictors <- c("P_WetSeason")
missing_required_predictors <- setdiff(required_predictors, names(merged_data))
if (length(missing_required_predictors) > 0) {
  stop(
    paste(
      "Missing required eco predictor variable(s):",
      paste(missing_required_predictors, collapse = ", ")
    )
  )
}

eco_predictors_all <- c(
  "P_WetSeason",
  STORAGE_METRIC_ORDER[STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD", "WB", "CHS")]
)
eco_predictors_all <- eco_predictors_all[eco_predictors_all %in% names(merged_data)]

if (length(eco_predictors_all) == 0) {
  stop("No eco predictors found in master annual dataset")
}

# Guardrail: eco models must not include static watershed-characteristic predictors.
disallowed_watershed_predictors <- c(
  "basin_slope", "Harvest",
  "Landslide_Total", "Landslide_Young",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
)
if (any(disallowed_watershed_predictors %in% eco_predictors_all)) {
  stop(
    "Eco predictor set contains watershed-characteristic predictors; ",
    "remove them from eco MLR."
  )
}

# Eco models use only annual storage metrics plus wet-season precipitation.
candidate_predictor_sets <- list(eco_predictors_all)

calc_loocv_stats <- function(model_obj, model_formula, df) {
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
    rmse_loocv <- tryCatch(
      as.numeric(MuMIn::loo(model_obj, type = "rmse")),
      error = function(e) NA_real_
    )
  }

  response_var <- all.vars(model_formula)[1]
  obs <- df[[response_var]]
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
  # Each LOOCV test run has one held-out observation; RMSE per run is |error|.
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

fit_one_model <- function(df_in, response, predictors,
                          use_scaled_predictors = FALSE,
                          use_iterative_vif = FALSE,
                          vif_threshold = 10,
                          min_n = 20,
                          model_id = ALL_SITES_ID) {
  model_df_in <- df_in

  predictors_use <- predictors[predictors %in% names(model_df_in)]
  if (length(predictors_use) == 0) {
    return(NULL)
  }

  # Keep only predictors with enough overlap with the response and real variation.
  pred_stats <- lapply(predictors_use, function(p) {
    x <- model_df_in[[p]]
    y <- model_df_in[[response]]
    n_pair <- sum(is.finite(x) & is.finite(y))
    sdx <- suppressWarnings(sd(x, na.rm = TRUE))
    tibble(Predictor = p, n_pair = n_pair, sdx = sdx)
  }) %>% bind_rows()

  predictors_use <- pred_stats %>%
    filter(is.finite(sdx), sdx > 0, n_pair >= min_n) %>%
    arrange(desc(n_pair)) %>%
    pull(Predictor)
  if (length(predictors_use) == 0) {
    return(NULL)
  }

  # Build the site-level modeling table using only rows with complete data.
  # If we do not have enough years, drop the sparsest predictor and try again.
  build_model_df <- function(preds) {
    keep_cols <- unique(c("site", "year", response, preds))
    keep_cols <- keep_cols[keep_cols %in% names(model_df_in)]
    model_df_in %>%
      dplyr::select(all_of(keep_cols)) %>%
      na.omit()
  }

  model_df <- build_model_df(predictors_use)
  while (nrow(model_df) < min_n && length(predictors_use) > 1) {
    keep_counts <- sapply(predictors_use, function(p) {
      sum(is.finite(model_df_in[[p]]) & is.finite(model_df_in[[response]]))
    })
    drop_var <- names(which.min(keep_counts))
    predictors_use <- setdiff(predictors_use, drop_var)
    model_df <- build_model_df(predictors_use)
  }

  if (nrow(model_df) < min_n || length(predictors_use) == 0) {
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
  fallback_full_model <- FALSE
  if (length(setdiff(names(coef(lm_aic)), "(Intercept)")) == 0) {
    # If stepwise AIC drops everything, keep the screened full model instead.
    # This avoids returning a blank model when data are otherwise usable.
    lm_aic <- lm_full
    fallback_full_model <- TRUE
  }

  if (use_iterative_vif) {
    repeat {
      retained_iter <- setdiff(names(coef(lm_aic)), "(Intercept)")
      if (length(retained_iter) <= 1) {
        break
      }

      vif_vals_iter <- tryCatch(vif(lm_aic), error = function(e) NULL)
      if (is.null(vif_vals_iter)) {
        break
      }
      if (max(vif_vals_iter, na.rm = TRUE) <= vif_threshold) {
        break
      }

      drop_var <- names(which.max(vif_vals_iter))
      retained_next <- setdiff(retained_iter, drop_var)
      if (length(retained_next) == 0) {
        return(NULL)
      }

      lm_next <- lm(
        as.formula(paste(response, "~", paste(retained_next, collapse = " + "))),
        data = model_df
      )
      lm_aic <- tryCatch(
        stepAIC(lm_next, direction = "backward", trace = 0),
        error = function(e) lm_next
      )
    }
  }

  retained <- setdiff(names(coef(lm_aic)), "(Intercept)")
  if (length(retained) == 0) {
    return(NULL)
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
  fstat <- suppressWarnings(as.numeric(lm_summary$fstatistic))
  model_p_global <- if (!is.null(fstat) && length(fstat) >= 3 && all(is.finite(fstat[1:3]))) {
    pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  } else {
    NA_real_
  }

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
      Site = model_id,
      Response = response,
      R2 = lm_summary$r.squared,
      R2_adj = lm_summary$adj.r.squared,
      RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
      AIC = aic_val,
      AICc = aicc_val,
      n = n_obs
    ) %>%
    mutate(selection_fallback_full_model = fallback_full_model) %>%
    dplyr::select(Site, Response, Predictor = variable, Beta_Std = beta_std, p_value = `Pr(>|t|)`, VIF, R2, R2_adj, RMSE, AIC, AICc, n)

  summary_out <- tibble(
    Site = model_id,
    Response = response,
    Predictors_Final = paste(retained, collapse = "; "),
    R2 = lm_summary$r.squared,
    R2_adj = lm_summary$adj.r.squared,
    model_p_global = model_p_global,
    RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
    AIC = aic_val,
    AICc = aicc_val,
    n = n_obs,
    selection_fallback_full_model = fallback_full_model
  )

  loocv <- calc_loocv_stats(lm_aic, formula(lm_aic), model_df)
  summary_out <- summary_out %>%
    mutate(
      RMSE_LOOCV = loocv$rmse_loocv,
      RMSE_LOOCV_MEAN_RUNS = loocv$rmse_loocv_mean_runs,
      R2_LOOCV = loocv$r2_loocv,
      delta_RMSE_LOOCV_minus_model = RMSE_LOOCV - RMSE,
      delta_RMSE_LOOCV_mean_runs_minus_model = RMSE_LOOCV_MEAN_RUNS - RMSE
    )

  list(
    model = lm_aic,
    data = model_df,
    coef = coef_out,
    summary = summary_out
  )
}

diagnose_pooled_response <- function(df_in, response, predictors, min_n = 20) {
  if (!(response %in% names(df_in))) {
    return(list(reason = "response_missing", n_response = 0, usable_predictors = 0))
  }

  n_response <- sum(is.finite(df_in[[response]]))
  if (n_response < min_n) {
    return(list(reason = "insufficient_response_n", n_response = n_response, usable_predictors = 0))
  }

  predictors_use <- predictors[predictors %in% names(df_in)]
  if (length(predictors_use) == 0) {
    return(list(reason = "no_predictors_present", n_response = n_response, usable_predictors = 0))
  }

  pred_stats <- lapply(predictors_use, function(p) {
    x <- df_in[[p]]
    y <- df_in[[response]]
    tibble(
      predictor = p,
      n_pair = sum(is.finite(x) & is.finite(y)),
      sdx = suppressWarnings(sd(x, na.rm = TRUE))
    )
  }) %>% bind_rows()

  usable <- pred_stats %>%
    filter(is.finite(sdx), sdx > 0, n_pair >= min_n)

  if (nrow(usable) == 0) {
    return(list(reason = "no_usable_predictors_after_screening", n_response = n_response, usable_predictors = 0))
  }

  keep_cols <- unique(c("site", "year", response, usable$predictor))
  model_df <- df_in %>%
    dplyr::select(all_of(keep_cols[keep_cols %in% names(df_in)])) %>%
    na.omit()
  if (nrow(model_df) < min_n) {
    return(list(reason = "insufficient_complete_cases", n_response = n_response, usable_predictors = nrow(usable)))
  }

  list(reason = "fit_failed_after_selection", n_response = n_response, usable_predictors = nrow(usable))
}

run_strict_method <- function() {
  model_results <- list()
  model_summary <- list()
  model_selection <- list()
  model_coverage <- list()
  model_diagnostics <- list()

  for (response in response_vars) {
    candidate_fits <- list()
    candidate_summary <- list()

    for (i in seq_along(candidate_predictor_sets)) {
      fit_obj <- fit_one_model(
        merged_data,
        response,
        candidate_predictor_sets[[i]],
        use_scaled_predictors = STRICT_USE_SCALED_PREDICTORS,
        use_iterative_vif = STRICT_USE_ITERATIVE_VIF,
        vif_threshold = VIF_THRESHOLD,
        min_n = MODEL_MIN_N,
        model_id = ALL_SITES_ID
      )

      if (!is.null(fit_obj)) {
        fit_obj$coef <- fit_obj$coef %>% mutate(Candidate_Set = i)
        fit_obj$summary <- fit_obj$summary %>% mutate(Candidate_Set = i)
        candidate_fits[[length(candidate_fits) + 1]] <- fit_obj
        candidate_summary[[length(candidate_summary) + 1]] <- fit_obj$summary
      }
    }

    if (length(candidate_fits) == 0) {
      diag <- diagnose_pooled_response(
        merged_data,
        response,
        eco_predictors_all,
        min_n = MODEL_MIN_N
      )
      model_coverage[[response]] <- tibble(
        Site = ALL_SITES_ID,
        Response = response,
        model_status = "not_fit",
        reason_not_fit = diag$reason,
        n_response = diag$n_response,
        usable_predictors = diag$usable_predictors
      )
      next
    }

    cand_tbl <- bind_rows(candidate_summary) %>%
      mutate(delta_AICc = AICc - min(AICc, na.rm = TRUE)) %>%
      arrange(AICc)
    model_selection[[response]] <- cand_tbl

    best_idx <- which.min(sapply(candidate_fits, function(x) x$summary$AICc))
    best_fit <- candidate_fits[[best_idx]]

    model_results[[response]] <- best_fit$coef
    model_summary[[response]] <- best_fit$summary
    model_diagnostics[[response]] <- compute_residual_diagnostics(best_fit$model, best_fit$data) %>%
      mutate(
        Site = ALL_SITES_ID,
        Response = response,
        Predictors_Final = best_fit$summary$Predictors_Final[1],
        n = best_fit$summary$n[1]
      ) %>%
      dplyr::select(
        Site, Response, Predictors_Final, n, n_residuals,
        shapiro_W, shapiro_p, ncv_chisq, ncv_p,
        normality_pass_p05, homoscedasticity_pass_p05
      )
    n_response_all <- merged_data %>%
      summarise(n_response = sum(is.finite(.data[[response]]), na.rm = TRUE)) %>%
      pull(n_response)
    model_coverage[[response]] <- tibble(
      Site = ALL_SITES_ID,
      Response = response,
      model_status = "fit",
      reason_not_fit = NA_character_,
      n_response = n_response_all,
      usable_predictors = NA_integer_
    )
  }

  list(
    results = bind_rows(model_results),
    summary = bind_rows(model_summary),
    selection = bind_rows(model_selection),
    coverage = bind_rows(model_coverage),
    diagnostics = bind_rows(model_diagnostics)
  )
}

model_run <- run_strict_method()

if (nrow(model_run$summary) == 0) {
  stop("No pooled eco models could be fit. Check predictor availability and MODEL_MIN_N.")
}

model_run$summary <- model_run$summary %>%
  mutate(
    low_n_flag = .data$n < LOW_N_THRESHOLD_ECO,
    confidence_note = ifelse(low_n_flag, "lower confidence (small n)", "standard confidence")
  )

format_export_response <- function(x) {
  gsub("_", "", as.character(x), fixed = TRUE)
}

format_export_predictor_text <- function(x) {
  out <- as.character(x)
  out <- gsub("P_WetSeason", "Pws", out, fixed = TRUE)
  out <- gsub("_", "", out, fixed = TRUE)
  out
}

round_export_cols <- function(df, cols, digits = 3) {
  keep <- intersect(cols, names(df))
  if (length(keep) == 0) return(df)
  df %>%
    mutate(across(all_of(keep), ~ signif(.x, digits)))
}

model_results_combined <- model_run$results %>%
  arrange(match(Response, response_vars), Predictor)
model_summary_combined <- model_run$summary %>%
  arrange(match(Response, response_vars))
selection_combined <- model_run$selection %>%
  arrange(match(Response, response_vars), Candidate_Set, AICc)
coverage_combined <- model_run$coverage %>%
  arrange(match(Response, response_vars))
diagnostics_combined <- model_run$diagnostics %>%
  left_join(
    model_run$summary %>% dplyr::select(Site, Response, low_n_flag, confidence_note),
    by = c("Site", "Response")
  ) %>%
  arrange(match(Response, response_vars))
aicc_lt2 <- selection_combined %>%
  filter(is.finite(delta_AICc), delta_AICc <= 2) %>%
  arrange(match(Response, response_vars), delta_AICc, AICc, Candidate_Set)

model_results_export <- model_results_combined %>%
  mutate(
    Response = format_export_response(Response),
    Predictor = format_export_predictor_text(Predictor)
  ) %>%
  round_export_cols(
    c("Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "AICc")
  )

model_summary_export <- model_summary_combined %>%
  mutate(
    Response = format_export_response(Response),
    Predictors_Final = format_export_predictor_text(Predictors_Final)
  ) %>%
  round_export_cols(
    c(
      "R2", "R2_adj", "model_p_global", "RMSE", "AIC", "AICc",
      "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
      "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model"
    )
  )

selection_export <- selection_combined %>%
  mutate(
    Response = format_export_response(Response),
    Predictors_Final = format_export_predictor_text(Predictors_Final)
  ) %>%
  round_export_cols(
    c(
      "R2", "R2_adj", "model_p_global", "RMSE", "AIC", "AICc",
      "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
      "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model",
      "delta_AICc"
    )
  )

coverage_export <- coverage_combined %>%
  mutate(Response = format_export_response(Response))

diagnostics_export <- diagnostics_combined %>%
  mutate(
    Response = format_export_response(Response),
    Predictors_Final = format_export_predictor_text(Predictors_Final)
  ) %>%
  round_export_cols(c("shapiro_W", "shapiro_p", "ncv_chisq", "ncv_p"))

aicc_lt2_export <- aicc_lt2 %>%
  mutate(
    Response = format_export_response(Response),
    Predictors_Final = format_export_predictor_text(Predictors_Final)
  ) %>%
  round_export_cols(
    c(
      "R2", "R2_adj", "model_p_global", "RMSE", "AIC", "AICc",
      "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
      "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model",
      "delta_AICc"
    )
  )

cor_data <- merged_data %>%
  dplyr::select(all_of(unique(c(response_vars, eco_predictors_all))))

if (ncol(cor_data) >= 2) {
  cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")
  cor_response_predictors <- cor_matrix[response_vars, eco_predictors_all, drop = FALSE]
  write.csv(cor_response_predictors,
            file.path(output_dir, paste0(file_prefix, "_corr_matrix.csv")),
            row.names = TRUE)
}

write.csv(model_results_export,
          file.path(output_dir, paste0(file_prefix, "_results.csv")),
          row.names = FALSE)

write.csv(model_summary_export,
          file.path(output_dir, paste0(file_prefix, "_summary.csv")),
          row.names = FALSE)
write.csv(selection_export,
          file.path(output_dir, paste0(file_prefix, "_model_selection.csv")),
          row.names = FALSE)
write.csv(
  diagnostics_export,
  file.path(output_dir, paste0(file_prefix, "_diagnostics.csv")),
  row.names = FALSE
)
write.csv(
  aicc_lt2_export,
  file.path(output_dir, paste0(file_prefix, "_aicc_lt2.csv")),
  row.names = FALSE
)
write.csv(coverage_export,
          file.path(output_dir, paste0(file_prefix, "_coverage.csv")),
          row.names = FALSE)

# Cleanup legacy strict-suffixed outputs; strict is now the default workflow.
legacy_strict_files <- c(
  file.path(output_dir, paste0(file_prefix, "_summary_strict.csv")),
  file.path(output_dir, paste0(file_prefix, "_results_strict.csv")),
  file.path(output_dir, paste0(file_prefix, "_corr_flags.csv"))
)
for (legacy_file in legacy_strict_files) {
  if (file.exists(legacy_file)) unlink(legacy_file)
}

# Explicit LOOCV validation output
loocv_validation <- model_run$summary %>%
  transmute(
    model_family = file_prefix,
    site = Site,
    response = Response,
    n = n,
    rmse_model = RMSE,
    rmse_loocv = RMSE_LOOCV,
    rmse_loocv_mean_runs = RMSE_LOOCV_MEAN_RUNS,
    r2_model = R2,
    r2_adj_model = R2_adj,
    r2_loocv = R2_LOOCV,
    delta_rmse_loocv_minus_model = delta_RMSE_LOOCV_minus_model,
    delta_rmse_loocv_mean_runs_minus_model = delta_RMSE_LOOCV_mean_runs_minus_model
  ) %>%
  arrange(match(response, response_vars))

write.csv(
  loocv_validation,
  file.path(OUT_STATS_VALIDATION_DIR, paste0(file_prefix, "_loocv_validation.csv")),
  row.names = FALSE
)
if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  if (!dir.exists(file.path(OUT_TABLES_DIR, "validation"))) {
    dir.create(file.path(OUT_TABLES_DIR, "validation"), recursive = TRUE, showWarnings = FALSE)
  }
  write.csv(
    loocv_validation,
    file.path(OUT_TABLES_DIR, "validation", paste0(file_prefix, "_loocv_validation.csv")),
    row.names = FALSE
  )
  if (!dir.exists(OUT_TABLES_MLR_DIR)) {
    dir.create(OUT_TABLES_MLR_DIR, recursive = TRUE, showWarnings = FALSE)
  }
  write.csv(
    aicc_lt2,
    file.path(OUT_TABLES_MLR_DIR, paste0(file_prefix, "_aicc_lt2.csv")),
    row.names = FALSE
  )
  write.csv(
    diagnostics_combined,
    file.path(OUT_TABLES_MLR_DIR, paste0(file_prefix, "_diagnostics.csv")),
    row.names = FALSE
  )
}
