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

output_dir <- OUT_STATS_MLR_ECO_DIR
base_dir <- BASE_DATA_DIR
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
file_prefix <- "storage_eco_mlr"

MODEL_MIN_N <- 20
VIF_THRESHOLD <- 10
ALLOW_LAVA_AND_ASH_TOGETHER <- FALSE
ALLOW_TOTAL_AND_YOUNG_LS_TOGETHER <- FALSE

env_vif_threshold <- suppressWarnings(as.numeric(Sys.getenv("HJA_MLR_VIF_THRESHOLD", unset = "")))
if (!is.na(env_vif_threshold) && is.finite(env_vif_threshold) && env_vif_threshold > 0) {
  VIF_THRESHOLD <- env_vif_threshold
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

if (!("temp_during_min_Q_7d_C" %in% names(merged_data)) && ("temp_at_min_Q_7d_C" %in% names(merged_data))) {
  merged_data <- merged_data %>%
    mutate(temp_during_min_Q_7d_C = temp_at_min_Q_7d_C)
}

response_vars_required <- c("max_temp_7d_C", "min_Q_7d_mm_d", "temp_during_min_Q_7d_C")
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

storage_predictors_all <- c("RCS", "RBI", "FDC", "SD", "WB", "CHS", "MTT", "Fyw", "DR")
storage_predictors_all <- storage_predictors_all[storage_predictors_all %in% names(merged_data)]

geology_predictors <- c("Lava1_per", "Lava2_per", "Ash_Per")
geology_predictors <- geology_predictors[geology_predictors %in% names(merged_data)]

landslide_predictors <- c("Landslide_Total", "Landslide_Young")
landslide_predictors <- landslide_predictors[landslide_predictors %in% names(merged_data)]

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

fit_one_model <- function(df_in, response, predictors,
                          use_scaled_predictors = FALSE,
                          use_iterative_vif = FALSE,
                          vif_threshold = 10,
                          min_n = 20,
                          enforce_correlated_exclusions = FALSE) {
  pca_screen_predictors <- function(df, predictor_names, conf_level = 0.99, collinear_cutoff = 0.7) {
    p_use <- predictor_names[predictor_names %in% names(df)]
    if (length(p_use) <= 1) {
      return(list(kept = p_use, detail = tibble(), retained_pcs = NA_integer_))
    }

    x <- df %>% dplyr::select(all_of(p_use)) %>% na.omit()
    if (nrow(x) < 6) {
      return(list(kept = p_use, detail = tibble(), retained_pcs = NA_integer_))
    }

    x_scaled <- scale(x)
    pca_obj <- prcomp(x_scaled, center = TRUE, scale. = TRUE)
    eig <- pca_obj$sdev^2
    var_exp <- eig / sum(eig)

    k_idx <- seq_along(eig)
    line_y <- eig[1] + (eig[length(eig)] - eig[1]) * (k_idx - 1) / (length(eig) - 1)
    elbow_k <- which.max(eig - line_y)
    elbow_k <- max(1, elbow_k)

    pc_scores <- pca_obj$x[, seq_len(elbow_k), drop = FALSE]
    p_cut <- 1 - conf_level

    sig_rows <- lapply(p_use, function(v) {
      x_vec <- as.numeric(x_scaled[, v])
      p_vals <- sapply(seq_len(elbow_k), function(i) {
        suppressWarnings(cor.test(x_vec, pc_scores[, i])$p.value)
      })
      tibble(
        Predictor = v,
        min_p = min(p_vals, na.rm = TRUE),
        significant_any_pc = any(p_vals < p_cut, na.rm = TRUE)
      )
    })
    sig_tbl <- bind_rows(sig_rows)
    keep_sig <- sig_tbl %>%
      filter(significant_any_pc) %>%
      pull(Predictor)
    if (length(keep_sig) == 0) keep_sig <- p_use

    loadings <- pca_obj$rotation[keep_sig, seq_len(elbow_k), drop = FALSE]
    w <- var_exp[seq_len(elbow_k)]
    score_tbl <- tibble(
      Predictor = keep_sig,
      pca_score = as.numeric(loadings^2 %*% w)
    )

    keep_final <- keep_sig
    drop_tbl <- tibble()
    repeat {
      if (length(keep_final) <= 1) break
      cm <- suppressWarnings(cor(x[, keep_final, drop = FALSE], use = "pairwise.complete.obs"))
      cm[lower.tri(cm, diag = TRUE)] <- NA_real_
      max_corr <- max(abs(cm), na.rm = TRUE)
      if (!is.finite(max_corr) || max_corr < collinear_cutoff) break

      idx <- which(abs(cm) == max_corr, arr.ind = TRUE)[1, ]
      a <- colnames(cm)[idx[2]]
      b <- rownames(cm)[idx[1]]
      sa <- score_tbl$pca_score[match(a, score_tbl$Predictor)]
      sb <- score_tbl$pca_score[match(b, score_tbl$Predictor)]
      drop_var <- ifelse(sa >= sb, b, a)
      keep_var <- ifelse(sa >= sb, a, b)
      drop_tbl <- bind_rows(drop_tbl, tibble(dropped = drop_var, kept = keep_var, corr_abs = max_corr, reason = "collinear_pair"))
      keep_final <- setdiff(keep_final, drop_var)
    }

    dropped_vars <- if ("dropped" %in% names(drop_tbl)) drop_tbl$dropped else character(0)
    detail <- sig_tbl %>%
      left_join(score_tbl, by = "Predictor") %>%
      mutate(
        retained_after_pca = Predictor %in% keep_final,
        drop_reason = ifelse(Predictor %in% dropped_vars, "collinear_pair", NA_character_)
      )
    list(kept = keep_final, detail = detail, retained_pcs = elbow_k)
  }

  calc_loocv_stats <- function(model_formula, df) {
    n <- nrow(df)
    if (n < 6) return(list(rmse_loocv = NA_real_, r2_loocv = NA_real_))
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
    if (sum(valid) < 3) return(list(rmse_loocv = NA_real_, r2_loocv = NA_real_))
    rmse_loocv <- sqrt(mean((obs[valid] - preds[valid])^2, na.rm = TRUE))
    sst <- sum((obs[valid] - mean(obs[valid], na.rm = TRUE))^2, na.rm = TRUE)
    sse <- sum((obs[valid] - preds[valid])^2, na.rm = TRUE)
    r2_loocv <- ifelse(sst > 0, 1 - (sse / sst), NA_real_)
    list(rmse_loocv = rmse_loocv, r2_loocv = r2_loocv)
  }

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

  pca_filter <- pca_screen_predictors(model_df, predictors_use)
  predictors_use <- pca_filter$kept
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
    repeat {
      retained_iter <- setdiff(names(coef(lm_aic)), "(Intercept)")
      has_ash <- "Ash_Per" %in% retained_iter
      has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained_iter)
      has_ls_both <- all(c("Landslide_Total", "Landslide_Young") %in% retained_iter)
      if (!(has_ash && has_lava) && !has_ls_both) break

      drop_var <- NULL
      if (has_ash && has_lava) {
        drop_var <- "Ash_Per"
      } else if (has_ls_both) {
        vif_vals_iter <- tryCatch(vif(lm_aic), error = function(e) NULL)
        if (!is.null(vif_vals_iter) && all(c("Landslide_Total", "Landslide_Young") %in% names(vif_vals_iter))) {
          drop_var <- names(which.max(vif_vals_iter[c("Landslide_Total", "Landslide_Young")]))
        } else {
          drop_var <- "Landslide_Total"
        }
      }

      if (is.null(drop_var)) break
      retained_next <- setdiff(retained_iter, drop_var)
      if (length(retained_next) == 0) return(NULL)
      lm_next <- lm(
        as.formula(paste(response, "~", paste(retained_next, collapse = " + "))),
        data = model_df
      )
      lm_aic <- tryCatch(stepAIC(lm_next, direction = "backward", trace = 0), error = function(e) lm_next)
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
      response = response,
      R2 = lm_summary$r.squared,
      R2_adj = lm_summary$adj.r.squared,
      RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
      AIC = aic_val,
      AICc = aicc_val,
      n = n_obs
    ) %>%
    dplyr::select(response, variable, beta_std, `Pr(>|t|)`, VIF, R2, R2_adj, RMSE, AIC, AICc, n)

  names(coef_out) <- c("Response", "Predictor", "Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "AICc", "n")

  summary_out <- tibble(
    Response = response,
    PCA_PCs_Retained = pca_filter$retained_pcs,
    Predictors_PCA = paste(predictors_use, collapse = "; "),
    Predictors_Final = paste(retained, collapse = "; "),
    R2 = lm_summary$r.squared,
    R2_adj = lm_summary$adj.r.squared,
    RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
    AIC = aic_val,
    AICc = aicc_val,
    n = n_obs,
    flag_lava_and_ash = flag_lava_ash,
    flag_landslide_total_and_young = flag_ls_total_young,
    constraint_flag = constraint_flag
  )

  loocv <- calc_loocv_stats(formula(lm_aic), model_df)
  summary_out <- summary_out %>%
    mutate(
      RMSE_LOOCV = loocv$rmse_loocv,
      R2_LOOCV = loocv$r2_loocv,
      delta_RMSE_LOOCV_minus_model = RMSE_LOOCV - RMSE
    )

  list(
    model = lm_aic,
    data = model_df,
    coef = coef_out,
    summary = summary_out,
    pca_screening = pca_filter$detail %>% mutate(Response = response)
  )
}

run_method <- function(method_name, use_scaled_predictors, use_iterative_vif, enforce_correlated_exclusions = FALSE) {
  model_results <- list()
  model_summary <- list()
  model_selection <- list()
  pca_screening <- list()

  for (response in response_vars) {
    candidate_fits <- list()
    candidate_summary <- list()

    for (i in seq_along(candidate_predictor_sets)) {
      fit_obj <- fit_one_model(
        merged_data,
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
        pca_screening[[paste(response, i, sep = "_")]] <- fit_obj$pca_screening %>%
          mutate(Candidate_Set = i, Method = method_name)
      }
    }

    if (length(candidate_fits) == 0) {
      next
    }

    cand_tbl <- bind_rows(candidate_summary) %>%
      mutate(delta_AICc = AICc - min(AICc, na.rm = TRUE)) %>%
      arrange(AICc)
    model_selection[[response]] <- cand_tbl

    best_idx <- which.min(sapply(candidate_fits, function(x) x$summary$AIC))
    best_fit <- candidate_fits[[best_idx]]

    model_results[[response]] <- best_fit$coef
    model_summary[[response]] <- best_fit$summary
  }

  list(
    results = bind_rows(model_results),
    summary = bind_rows(model_summary),
    selection = bind_rows(model_selection),
    pca_screening = bind_rows(pca_screening)
  )
}

non_strict <- run_method(
  method_name = "non_strict",
  use_scaled_predictors = FALSE,
  use_iterative_vif = FALSE,
  enforce_correlated_exclusions = FALSE
)
strict <- run_method(
  method_name = "strict",
  use_scaled_predictors = TRUE,
  use_iterative_vif = TRUE,
  enforce_correlated_exclusions = TRUE
)

model_results_combined <- strict$results
model_summary_combined <- strict$summary
selection_combined <- bind_rows(non_strict$selection, strict$selection)
pca_screening_combined <- bind_rows(non_strict$pca_screening, strict$pca_screening)

cor_data <- merged_data %>%
  dplyr::select(all_of(unique(c(response_vars, storage_predictors_all))))

if (ncol(cor_data) >= 2) {
  cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")
  cor_response_storage <- cor_matrix[response_vars, storage_predictors_all, drop = FALSE]
  write.csv(cor_response_storage,
            file.path(output_dir, "storage_eco_corr_matrix.csv"),
            row.names = TRUE)
}

write.csv(model_results_combined,
          file.path(output_dir, paste0(file_prefix, "_results.csv")),
          row.names = FALSE)

write.csv(model_summary_combined,
          file.path(output_dir, paste0(file_prefix, "_summary.csv")),
          row.names = FALSE)

write.csv(non_strict$results,
          file.path(output_dir, paste0(file_prefix, "_results_non_strict.csv")),
          row.names = FALSE)
write.csv(non_strict$summary,
          file.path(output_dir, paste0(file_prefix, "_summary_non_strict.csv")),
          row.names = FALSE)
write.csv(strict$results,
          file.path(output_dir, paste0(file_prefix, "_results_strict.csv")),
          row.names = FALSE)
write.csv(strict$summary,
          file.path(output_dir, paste0(file_prefix, "_summary_strict.csv")),
          row.names = FALSE)

write.csv(bind_rows(non_strict$results, strict$results),
          file.path(output_dir, paste0(file_prefix, "_results_all_methods.csv")),
          row.names = FALSE)
write.csv(bind_rows(non_strict$summary, strict$summary),
          file.path(output_dir, paste0(file_prefix, "_summary_all_methods.csv")),
          row.names = FALSE)
write.csv(selection_combined,
          file.path(output_dir, paste0(file_prefix, "_model_selection_all_methods.csv")),
          row.names = FALSE)
write.csv(pca_screening_combined,
          file.path(output_dir, paste0(file_prefix, "_pca_screening_all_methods.csv")),
          row.names = FALSE)

comparison <- non_strict$summary %>%
  dplyr::select(Response, Predictors_PCA, Predictors_Final, R2_adj, RMSE, RMSE_LOOCV, AIC, AICc, n) %>%
  rename_with(~ paste0(.x, "_non_strict"), -Response) %>%
  full_join(
    strict$summary %>%
      dplyr::select(Response, Predictors_PCA, Predictors_Final, R2_adj, RMSE, RMSE_LOOCV, AIC, AICc, n) %>%
      rename_with(~ paste0(.x, "_strict"), -Response),
    by = "Response"
  ) %>%
  mutate(
    delta_R2_adj = R2_adj_strict - R2_adj_non_strict,
    delta_RMSE = RMSE_strict - RMSE_non_strict,
    delta_RMSE_LOOCV = RMSE_LOOCV_strict - RMSE_LOOCV_non_strict,
    delta_AIC = AIC_strict - AIC_non_strict,
    delta_AICc = AICc_strict - AICc_non_strict
  ) %>%
  arrange(Response)

write.csv(comparison,
          file.path(output_dir, paste0(file_prefix, "_strict_vs_non_strict.csv")),
          row.names = FALSE)

flags <- bind_rows(non_strict$summary, strict$summary) %>%
  filter(constraint_flag) %>%
  arrange(Method, Response)
write.csv(flags,
          file.path(output_dir, paste0(file_prefix, "_corr_flags.csv")),
          row.names = FALSE)
