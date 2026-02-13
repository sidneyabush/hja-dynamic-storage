# -----------------------------------------------------------------------------
# Predict Thermal and Low-Flow Metrics from Storage and Catchment Predictors
# -----------------------------------------------------------------------------
# Purpose: Fit stepwise-AIC linear models for annual eco response metrics,
#          separately for each site and response.
#
# Response metrics:
#   1) T_7DMax
#   2) Q_7Q5
#   3) T_Q7Q5
#
# Inputs:
#   - master_annual.csv
#
# Outputs:
#   - storage_ecovar_mlr_results.csv
#   - storage_ecovar_mlr_summary.csv
#   - storage_ecovar_mlr_model_selection.csv
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(MASS)
library(car)
library(tidyr)

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

output_dir <- OUT_MODELS_STORAGE_ECOVAR_MLR_DIR
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
file_prefix <- "storage_ecovar_mlr"

MODEL_MIN_N <- 20
VIF_THRESHOLD <- 10
LOW_N_THRESHOLD_ECO <- 25
ALLOW_LAVA_AND_ASH_TOGETHER <- FALSE
ALLOW_TOTAL_AND_YOUNG_LS_TOGETHER <- FALSE
USE_MUMIN_LOOCV <- TRUE

env_vif_threshold <- suppressWarnings(as.numeric(Sys.getenv("HJA_MLR_VIF_THRESHOLD", unset = "")))
if (!is.na(env_vif_threshold) && is.finite(env_vif_threshold) && env_vif_threshold > 0) {
  VIF_THRESHOLD <- env_vif_threshold
}

if (USE_MUMIN_LOOCV && !requireNamespace("MuMIn", quietly = TRUE)) {
  stop("MuMIn is required for eco-model LOOCV. Install with: install.packages('MuMIn')")
}

annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUTPUT_DIR, MASTER_ANNUAL_FILE)
}

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

storage_predictors_all <- c("RCS", "RBI", "FDC", "SD", "WB", "CHS")
storage_predictors_all <- storage_predictors_all[storage_predictors_all %in% names(merged_data)]

geology_predictors <- c("Lava1_per", "Lava2_per", "Ash_Per")
geology_predictors <- geology_predictors[geology_predictors %in% names(merged_data)]

landslide_predictors <- c("Landslide_Total", "Landslide_Young")
landslide_predictors <- landslide_predictors[landslide_predictors %in% names(merged_data)]

if (length(storage_predictors_all) == 0) {
  stop("No storage predictors found in master annual dataset")
}

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

fit_one_model <- function(df_in, site_id, response, predictors,
                          use_scaled_predictors = FALSE,
                          use_iterative_vif = FALSE,
                          vif_threshold = 10,
                          min_n = 20,
                          enforce_correlated_exclusions = FALSE) {
  site_df <- df_in %>%
    filter(site == site_id)

  predictors_use <- predictors[predictors %in% names(site_df)]
  if (length(predictors_use) == 0) {
    return(NULL)
  }

  # Drop predictors that have too little support or are constant within site.
  pred_stats <- lapply(predictors_use, function(p) {
    x <- site_df[[p]]
    y <- site_df[[response]]
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

  # Build a complete-case model frame and iteratively prune sparse predictors
  # until the requested minimum annual sample size is met.
  build_model_df <- function(preds) {
    keep_cols <- unique(c("site", "year", response, preds))
    keep_cols <- keep_cols[keep_cols %in% names(site_df)]
    site_df %>%
      dplyr::select(all_of(keep_cols)) %>%
      na.omit()
  }

  model_df <- build_model_df(predictors_use)
  while (nrow(model_df) < min_n && length(predictors_use) > 1) {
    keep_counts <- sapply(predictors_use, function(p) {
      sum(is.finite(site_df[[p]]) & is.finite(site_df[[response]]))
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
    # If AIC selection collapses to intercept-only, retain the screened full model
    # so each site-response with usable data still yields coefficients.
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

  if (enforce_correlated_exclusions) {
    fit_candidate_from_retained <- function(retained_set) {
      if (length(retained_set) == 0) return(NULL)
      lm_next <- tryCatch(
        lm(as.formula(paste(response, "~", paste(retained_set, collapse = " + "))), data = model_df),
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
      retained_iter <- setdiff(names(coef(lm_aic)), "(Intercept)")
      has_ash <- "Ash_Per" %in% retained_iter
      has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained_iter)
      has_ls_both <- all(c("Landslide_Total", "Landslide_Young") %in% retained_iter)
      if (!(has_ash && has_lava) && !has_ls_both) break

      if (has_ash && has_lava) {
        drop_options <- c("Ash_Per")
        if ("Lava1_per" %in% retained_iter) drop_options <- c(drop_options, "Lava1_per")
        if ("Lava2_per" %in% retained_iter) drop_options <- c(drop_options, "Lava2_per")
        cand_models <- lapply(drop_options, function(v) {
          retained_next <- setdiff(retained_iter, v)
          fit_candidate_from_retained(retained_next)
        })
        valid_idx <- which(vapply(cand_models, function(x) !is.null(x), logical(1)))
        if (length(valid_idx) == 0) break
        best_idx <- valid_idx[which.min(vapply(cand_models[valid_idx], function(x) x$aicc, numeric(1)))]
        lm_aic <- cand_models[[best_idx]]$model
      } else if (has_ls_both) {
        drop_options <- c("Landslide_Total", "Landslide_Young")
        cand_models <- lapply(drop_options, function(v) {
          retained_next <- setdiff(retained_iter, v)
          fit_candidate_from_retained(retained_next)
        })
        valid_idx <- which(vapply(cand_models, function(x) !is.null(x), logical(1)))
        if (length(valid_idx) == 0) break
        best_idx <- valid_idx[which.min(vapply(cand_models[valid_idx], function(x) x$aicc, numeric(1)))]
        lm_aic <- cand_models[[best_idx]]$model
      }
    }
  }

  retained <- setdiff(names(coef(lm_aic)), "(Intercept)")
  if (length(retained) == 0) {
    return(NULL)
  }

  has_ash <- "Ash_Per" %in% retained
  has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained)
  flag_lava_ash <- has_ash && has_lava
  flag_ls_total_young <- all(c("Landslide_Total", "Landslide_Young") %in% retained)
  constraint_flag <- flag_lava_ash || flag_ls_total_young

  n_obs <- nrow(model_df)
  k_params <- length(coef(lm_aic)) + 1
  aic_val <- AIC(lm_aic)
  aicc_val <- if ((n_obs - k_params - 1) > 0) {
    aic_val + (2 * k_params * (k_params + 1)) / (n_obs - k_params - 1)
  } else {
    NA_real_
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
      Site = site_id,
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
    Site = site_id,
    Response = response,
    Predictors_Final = paste(retained, collapse = "; "),
    R2 = lm_summary$r.squared,
    R2_adj = lm_summary$adj.r.squared,
    RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
    AIC = aic_val,
    AICc = aicc_val,
    n = n_obs,
    selection_fallback_full_model = fallback_full_model,
    flag_lava_and_ash = flag_lava_ash,
    flag_landslide_total_and_young = flag_ls_total_young,
    constraint_flag = constraint_flag
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

diagnose_site_response <- function(df_in, site_id, response, predictors, min_n = 20) {
  site_df <- df_in %>%
    filter(site == site_id)

  if (!(response %in% names(site_df))) {
    return(list(reason = "response_missing", n_response = 0, usable_predictors = 0))
  }

  n_response <- sum(is.finite(site_df[[response]]))
  if (n_response < min_n) {
    return(list(reason = "insufficient_response_n", n_response = n_response, usable_predictors = 0))
  }

  predictors_use <- predictors[predictors %in% names(site_df)]
  if (length(predictors_use) == 0) {
    return(list(reason = "no_predictors_present", n_response = n_response, usable_predictors = 0))
  }

  pred_stats <- lapply(predictors_use, function(p) {
    x <- site_df[[p]]
    y <- site_df[[response]]
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
  model_df <- site_df %>%
    dplyr::select(all_of(keep_cols[keep_cols %in% names(site_df)])) %>%
    na.omit()
  if (nrow(model_df) < min_n) {
    return(list(reason = "insufficient_complete_cases", n_response = n_response, usable_predictors = nrow(usable)))
  }

  list(reason = "fit_failed_after_selection", n_response = n_response, usable_predictors = nrow(usable))
}

run_method <- function(method_name, use_scaled_predictors, use_iterative_vif, enforce_correlated_exclusions = FALSE) {
  model_results <- list()
  model_summary <- list()
  model_selection <- list()
  model_coverage <- list()

  site_order_for_model <- SITE_ORDER_HYDROMETRIC[SITE_ORDER_HYDROMETRIC %in% unique(merged_data$site)]

  for (site_id in site_order_for_model) {
    for (response in response_vars) {
      candidate_fits <- list()
      candidate_summary <- list()

      for (i in seq_along(candidate_predictor_sets)) {
        fit_obj <- fit_one_model(
          merged_data,
          site_id,
          response,
          candidate_predictor_sets[[i]],
          use_scaled_predictors = use_scaled_predictors,
          use_iterative_vif = use_iterative_vif,
          vif_threshold = VIF_THRESHOLD,
          min_n = MODEL_MIN_N,
          enforce_correlated_exclusions = enforce_correlated_exclusions
        )

        if (!is.null(fit_obj)) {
          fit_obj$coef <- fit_obj$coef %>% mutate(Candidate_Set = i, Method = method_name)
          fit_obj$summary <- fit_obj$summary %>% mutate(Candidate_Set = i, Method = method_name)
          candidate_fits[[length(candidate_fits) + 1]] <- fit_obj
          candidate_summary[[length(candidate_summary) + 1]] <- fit_obj$summary
        }
      }

      if (length(candidate_fits) == 0) {
        diag <- diagnose_site_response(
          merged_data,
          site_id,
          response,
          storage_predictors_all,
          min_n = MODEL_MIN_N
        )
        model_coverage[[paste(site_id, response, sep = "_")]] <- tibble(
          Site = site_id,
          Response = response,
          Method = method_name,
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
      model_selection[[paste(site_id, response, sep = "_")]] <- cand_tbl

      best_idx <- which.min(sapply(candidate_fits, function(x) x$summary$AICc))
      best_fit <- candidate_fits[[best_idx]]

      model_results[[paste(site_id, response, sep = "_")]] <- best_fit$coef
      model_summary[[paste(site_id, response, sep = "_")]] <- best_fit$summary
      n_response_site <- merged_data %>%
        filter(site == site_id) %>%
        summarise(n_response = sum(is.finite(.data[[response]]), na.rm = TRUE)) %>%
        pull(n_response)
      model_coverage[[paste(site_id, response, sep = "_")]] <- tibble(
        Site = site_id,
        Response = response,
        Method = method_name,
        model_status = "fit",
        reason_not_fit = NA_character_,
        n_response = n_response_site,
        usable_predictors = NA_integer_
      )
    }
  }

  list(
    results = bind_rows(model_results),
    summary = bind_rows(model_summary),
    selection = bind_rows(model_selection),
    coverage = bind_rows(model_coverage)
  )
}

strict <- run_method(
  method_name = "strict",
  use_scaled_predictors = TRUE,
  use_iterative_vif = TRUE,
  enforce_correlated_exclusions = TRUE
)

if (nrow(strict$summary) == 0) {
  stop("No site-response eco models could be fit. Check predictor availability and MODEL_MIN_N.")
}

strict$summary <- strict$summary %>%
  mutate(
    low_n_flag = .data$n < LOW_N_THRESHOLD_ECO,
    confidence_note = ifelse(low_n_flag, "lower confidence (small n)", "standard confidence")
  )

model_results_combined <- strict$results
model_summary_combined <- strict$summary
selection_combined <- strict$selection
coverage_combined <- strict$coverage %>%
  arrange(match(Site, SITE_ORDER_HYDROMETRIC), match(Response, response_vars))

cor_data <- merged_data %>%
  dplyr::select(all_of(unique(c(response_vars, storage_predictors_all))))

if (ncol(cor_data) >= 2) {
  cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")
  cor_response_storage <- cor_matrix[response_vars, storage_predictors_all, drop = FALSE]
  write.csv(cor_response_storage,
            file.path(output_dir, "storage_ecovar_mlr_corr_matrix.csv"),
            row.names = TRUE)
}

write.csv(model_results_combined,
          file.path(output_dir, paste0(file_prefix, "_results.csv")),
          row.names = FALSE)

write.csv(model_summary_combined,
          file.path(output_dir, paste0(file_prefix, "_summary.csv")),
          row.names = FALSE)

write.csv(strict$results,
          file.path(output_dir, paste0(file_prefix, "_results_strict.csv")),
          row.names = FALSE)
write.csv(strict$summary,
          file.path(output_dir, paste0(file_prefix, "_summary_strict.csv")),
          row.names = FALSE)
write.csv(selection_combined,
          file.path(output_dir, paste0(file_prefix, "_model_selection.csv")),
          row.names = FALSE)
write.csv(coverage_combined,
          file.path(output_dir, paste0(file_prefix, "_coverage.csv")),
          row.names = FALSE)

flags <- strict$summary %>%
  filter(constraint_flag) %>%
  arrange(Method, Site, Response)
write.csv(flags,
          file.path(output_dir, paste0(file_prefix, "_corr_flags.csv")),
          row.names = FALSE)

# Explicit LOOCV validation output (models + tables validation folders)
loocv_validation <- strict$summary %>%
  transmute(
    model_family = "storage_ecovar_mlr",
    method = Method,
    site = Site,
    response = Response,
    n = n,
    rmse_model = RMSE,
    rmse_loocv = RMSE_LOOCV,
    rmse_loocv_mean_runs = RMSE_LOOCV_MEAN_RUNS,
    r2_model = R2,
    r2_loocv = R2_LOOCV,
    delta_rmse_loocv_minus_model = delta_RMSE_LOOCV_minus_model,
    delta_rmse_loocv_mean_runs_minus_model = delta_RMSE_LOOCV_mean_runs_minus_model
  ) %>%
  arrange(match(site, SITE_ORDER_HYDROMETRIC), response)

write.csv(
  loocv_validation,
  file.path(OUT_STATS_VALIDATION_DIR, "storage_ecovar_mlr_loocv_validation.csv"),
  row.names = FALSE
)

if (!dir.exists(file.path(OUT_TABLES_DIR, "validation"))) {
  dir.create(file.path(OUT_TABLES_DIR, "validation"), recursive = TRUE, showWarnings = FALSE)
}
write.csv(
  loocv_validation,
  file.path(OUT_TABLES_DIR, "validation", "storage_ecovar_mlr_loocv_validation.csv"),
  row.names = FALSE
)
