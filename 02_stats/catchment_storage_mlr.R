# -----------------------------------------------------------------------------
# Multiple Linear Regression: Watershed Characteristics -> Storage Metrics
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

VIF_THRESHOLD <- 10

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
  "Slope_mean",
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

    # Scree elbow via max distance to line between first and last eigenvalue.
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
    if (length(keep_sig) == 0) {
      keep_sig <- p_use
    }

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

      drop_tbl <- bind_rows(drop_tbl, tibble(
        dropped = drop_var,
        kept = keep_var,
        corr_abs = max_corr,
        reason = "collinear_pair"
      ))
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

  pca_filter <- pca_screen_predictors(model_df, predictor_keep)
  predictor_keep <- pca_filter$kept
  if (length(predictor_keep) == 0) {
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
    repeat {
      retained_vars <- setdiff(names(coef(lm_aic)), "(Intercept)")
      has_ash <- "Ash_Per" %in% retained_vars
      has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained_vars)
      has_ls_both <- all(c("Landslide_Total", "Landslide_Young") %in% retained_vars)

      if (!(has_ash && has_lava) && !has_ls_both) {
        break
      }

      drop_var <- NULL
      if (has_ash && has_lava) {
        # Prefer lava representation over ash when both are retained.
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
      retained_next <- setdiff(retained_vars, drop_var)
      if (length(retained_next) == 0) return(NULL)
      lm_next <- lm(
        as.formula(paste(outcome, "~", paste(retained_next, collapse = " + "))),
        data = model_df
      )
      lm_aic <- tryCatch(stepAIC(lm_next, direction = "backward", trace = 0), error = function(e) lm_next)
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
      outcome = outcome,
      R2 = lm_summary$r.squared,
      R2_adj = lm_summary$adj.r.squared,
      RMSE = sqrt(mean(residuals(lm_aic)^2, na.rm = TRUE)),
      AIC = aic_val,
      AICc = aicc_val
    ) %>%
    dplyr::select(outcome, variable, beta_std, `Pr(>|t|)`, VIF, R2, R2_adj, RMSE, AIC, AICc)

  names(result_df) <- c("Outcome", "Predictor", "Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AIC", "AICc")

  summary_df <- tibble(
    Outcome = outcome,
    N = nrow(model_df),
    PCA_PCs_Retained = pca_filter$retained_pcs,
    Predictors_PCA = paste(predictor_keep, collapse = "; "),
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
    selection_path = tryCatch(as.data.frame(lm_aic$anova), error = function(e) NULL),
    pca_screening = pca_filter$detail %>% mutate(Outcome = outcome)
  )
}

run_method <- function(method_name, use_scaled_predictors, use_iterative_vif, enforce_correlated_exclusions = FALSE) {
  results_list <- list()
  model_summary_list <- list()
  selection_path_list <- list()
  pca_screening_list <- list()

  for (outcome in outcome_vars) {
    fit_obj <- fit_mlr_with_vif(
      HJA_Ave,
      outcome,
      predictor_vars,
      use_scaled_predictors = use_scaled_predictors,
      use_iterative_vif = use_iterative_vif,
      vif_threshold = VIF_THRESHOLD,
      enforce_correlated_exclusions = enforce_correlated_exclusions
    )
    if (is.null(fit_obj)) {
      next
    }
    results_list[[outcome]] <- fit_obj$coefficients %>% mutate(Method = method_name)
    model_summary_list[[outcome]] <- fit_obj$summary %>% mutate(Method = method_name)

    if (!is.null(fit_obj$selection_path) && nrow(fit_obj$selection_path) > 0) {
      sp <- fit_obj$selection_path
      sp$step <- seq_len(nrow(sp))
      sp$Outcome <- outcome
      sp$Method <- method_name
      selection_path_list[[outcome]] <- sp
    }
    if (!is.null(fit_obj$pca_screening) && nrow(fit_obj$pca_screening) > 0) {
      pca_screening_list[[outcome]] <- fit_obj$pca_screening %>% mutate(Method = method_name)
    }
  }

  list(
    results = bind_rows(results_list),
    summary = bind_rows(model_summary_list),
    selection = bind_rows(selection_path_list),
    pca_screening = bind_rows(pca_screening_list)
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

results_all_methods <- bind_rows(non_strict$results, strict$results)
summary_all_methods <- bind_rows(non_strict$summary, strict$summary)
selection_all_methods <- bind_rows(non_strict$selection, strict$selection)
pca_screening_all_methods <- bind_rows(non_strict$pca_screening, strict$pca_screening)

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

# Save method-specific outputs
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

# Save all-method stacks and selection paths
write.csv(results_all_methods,
          file.path(output_dir, paste0(file_prefix, "_results_all_methods.csv")),
          row.names = FALSE)
write.csv(summary_all_methods,
          file.path(output_dir, paste0(file_prefix, "_summary_all_methods.csv")),
          row.names = FALSE)
write.csv(selection_all_methods,
          file.path(output_dir, paste0(file_prefix, "_stepwise_path_all_methods.csv")),
          row.names = FALSE)
write.csv(pca_screening_all_methods,
          file.path(output_dir, paste0(file_prefix, "_pca_screening_all_methods.csv")),
          row.names = FALSE)

# Method comparison table for side-by-side review
comparison <- non_strict$summary %>%
  dplyr::select(Outcome, Predictors_PCA, Predictors_Final, R2_adj, RMSE, RMSE_LOOCV, AIC, AICc, N) %>%
  rename_with(~ paste0(.x, "_non_strict"), -Outcome) %>%
  full_join(
    strict$summary %>%
      dplyr::select(Outcome, Predictors_PCA, Predictors_Final, R2_adj, RMSE, RMSE_LOOCV, AIC, AICc, N) %>%
      rename_with(~ paste0(.x, "_strict"), -Outcome),
    by = "Outcome"
  ) %>%
  mutate(
    delta_R2_adj = R2_adj_strict - R2_adj_non_strict,
    delta_RMSE = RMSE_strict - RMSE_non_strict,
    delta_RMSE_LOOCV = RMSE_LOOCV_strict - RMSE_LOOCV_non_strict,
    delta_AIC = AIC_strict - AIC_non_strict,
    delta_AICc = AICc_strict - AICc_non_strict
  ) %>%
  arrange(Outcome)

write.csv(comparison,
          file.path(output_dir, paste0(file_prefix, "_strict_vs_non_strict.csv")),
          row.names = FALSE)

flags <- summary_all_methods %>%
  filter(constraint_flag) %>%
  arrange(Method, Outcome)
write.csv(flags,
          file.path(output_dir, paste0(file_prefix, "_corr_flags.csv")),
          row.names = FALSE)
