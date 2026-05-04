suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(MASS)
  library(car)
})

rm(list = ls())

source("config.R")
source("helpers/mlr_utils.R")

output_dir <- file.path(REPO_DIR, "outputs", "MTT_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

VIF_THRESHOLD <- 10
MODEL_MIN_N_ECO <- 20
MIN_N_CATCH <- 5

annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
isotope_file <- file.path(ISOTOPE_DIR, "MTT_FYW.csv")
damping_file <- file.path(ISOTOPE_DIR, "DampingRatios_2025-07-07.csv")

required_files <- c(annual_file, site_file, isotope_file, damping_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required file(s): ", paste(missing_files, collapse = "; "))
}

annual_df <- read_csv(annual_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

site_df_base <- read_csv(site_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC, !site %in% SITE_EXCLUDE_STANDARD)

raw_mtt <- read_csv(isotope_file, show_col_types = FALSE) %>%
  mutate(
    site = standardize_site_code(trimws(site)),
    MTT_early = suppressWarnings(as.numeric(MTT1)),
    MTT_late_mean = dplyr::coalesce(
      suppressWarnings(as.numeric(MTTM)),
      rowMeans(cbind(
        suppressWarnings(as.numeric(MTT2L)),
        suppressWarnings(as.numeric(MTT2H))
      ), na.rm = TRUE)
    ),
    Fyw = suppressWarnings(as.numeric(FYWM))
  ) %>%
  mutate(MTT_late_mean = ifelse(is.nan(MTT_late_mean), NA_real_, MTT_late_mean)) %>%
  dplyr::select(site, MTT_early, MTT_late_mean, Fyw) %>%
  filter(!is.na(site), site != "", site %in% SITE_ORDER_HYDROMETRIC)

damping <- read_csv(damping_file, show_col_types = FALSE) %>%
  mutate(
    site = standardize_site_code(trimws(site)),
    DR = suppressWarnings(as.numeric(DR_Overall))
  ) %>%
  dplyr::select(site, DR) %>%
  filter(!is.na(site), site != "", site %in% SITE_ORDER_HYDROMETRIC)

mtt_definitions <- raw_mtt %>%
  dplyr::transmute(
    site,
    early_only = MTT_early,
    late_only = MTT_late_mean,
    combined_mean = rowMeans(cbind(MTT_early, MTT_late_mean), na.rm = TRUE)
  ) %>%
  mutate(combined_mean = ifelse(is.nan(combined_mean), NA_real_, combined_mean)) %>%
  pivot_longer(
    cols = c(early_only, late_only, combined_mean),
    names_to = "mtt_definition",
    values_to = "MTT"
  ) %>%
  arrange(mtt_definition, site)

write_csv(mtt_definitions, file.path(output_dir, "mtt_site_definitions.csv"))

mtt_availability <- mtt_definitions %>%
  group_by(mtt_definition) %>%
  summarise(
    n_sites = sum(is.finite(MTT)),
    min_mtt = ifelse(any(is.finite(MTT)), min(MTT, na.rm = TRUE), NA_real_),
    max_mtt = ifelse(any(is.finite(MTT)), max(MTT, na.rm = TRUE), NA_real_),
    mean_mtt = ifelse(any(is.finite(MTT)), mean(MTT, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )

write_csv(mtt_availability, file.path(output_dir, "mtt_definition_summary.csv"))

label_mtt_definition <- function(x) {
  dplyr::case_when(
    x == "early_only" ~ "MTT1 (McGuire et al., 2005; 2001-2003)",
    x == "late_only" ~ "MTT2 (Segura, 2015-2018)",
    x == "combined_mean" ~ "Combined MTT (mean of MTT1 and MTT2)",
    x == "none" ~ "No MTT",
    TRUE ~ as.character(x)
  )
}

clean_predictor_labels <- function(x) {
  out <- as.character(x)
  out <- gsub("\\bCHS\\b", "BF", out)
  out <- gsub("\\bbasin_slope\\b", "Basin slope", out)
  out <- gsub("\\bPyro_per\\b", "Pyroclastic", out)
  out
}

metric_map_site <- c(
  "RBI" = "RBI_mean",
  "RCS" = "RCS_mean",
  "FDC" = "FDC_mean",
  "SD" = "SD_mean",
  "WB" = "WB_mean",
  "CHS" = "CHS_mean"
)

pairwise_corr <- function(df, x, y) {
  keep <- is.finite(df[[x]]) & is.finite(df[[y]])
  n_keep <- sum(keep)
  if (n_keep < 3) {
    return(tibble(
      n = n_keep,
      r = NA_real_,
      p_value = NA_real_
    ))
  }
  ct <- suppressWarnings(cor.test(df[[x]][keep], df[[y]][keep]))
  tibble(
    n = n_keep,
    r = unname(ct$estimate),
    p_value = ct$p.value
  )
}

site_corrs <- bind_rows(lapply(unique(mtt_definitions$mtt_definition), function(defn) {
  iso_tbl <- mtt_definitions %>%
    filter(mtt_definition == defn) %>%
    dplyr::select(site, MTT) %>%
    left_join(raw_mtt %>% dplyr::select(site, Fyw), by = "site") %>%
    left_join(damping, by = "site")

  site_tbl <- site_df_base %>%
    dplyr::select(site, all_of(unname(metric_map_site))) %>%
    left_join(iso_tbl, by = "site")

  bind_rows(lapply(c("RBI", "RCS", "FDC", "SD", "WB", "CHS", "DR", "Fyw"), function(metric) {
    metric_col <- if (metric %in% names(metric_map_site)) metric_map_site[[metric]] else metric
    pairwise_corr(site_tbl, metric_col, "MTT") %>%
      mutate(
        mtt_definition = defn,
        comparison_metric = metric
      )
  }))
}))

write_csv(
  site_corrs %>%
    dplyr::select(mtt_definition, comparison_metric, n, r, p_value),
  file.path(output_dir, "mtt_site_level_correlations.csv")
)

candidate_sets_catch <- unlist(
  lapply(seq_len(8), function(k) {
    combn(
      c("basin_slope", "Harvest", "Landslide_Total", "Landslide_Young",
        "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"),
      k,
      simplify = FALSE
    )
  }),
  recursive = FALSE
)

fit_catchment_mtt <- function(site_df_in) {
  fit_candidate <- function(df_in, predictors) {
    model_df <- df_in %>%
      dplyr::select(MTT, all_of(predictors)) %>%
      na.omit()

    if (nrow(model_df) < MIN_N_CATCH) {
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
      lm(as.formula(paste("MTT ~", paste(predictors, collapse = " + "))), data = model_df),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    fit <- apply_vif_filter(fit, "MTT", model_df, VIF_THRESHOLD)
    if (is.null(fit)) {
      return(NULL)
    }
    fit <- apply_correlated_predictor_rules_catchment(fit, "MTT", model_df)

    retained <- setdiff(names(coef(fit)), "(Intercept)")
    if (length(retained) == 0) {
      return(NULL)
    }

    fit_sum <- summary(fit)
    rmse <- sqrt(mean(residuals(fit)^2, na.rm = TRUE))
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

    beta_df <- model_df
    beta_df$MTT <- as.numeric(scale(beta_df$MTT))
    fit_beta <- lm(as.formula(paste("MTT ~", paste(retained, collapse = " + "))), data = beta_df)
    beta_tbl <- tibble(
      Predictor = names(coef(fit_beta))[-1],
      Beta_Std = as.numeric(coef(fit_beta)[-1])
    )

    diagnostics <- compute_residual_diagnostics(fit)

    list(
      summary = tibble(
        Predictors_Final = paste(retained, collapse = "; "),
        N = nrow(model_df),
        R2 = fit_sum$r.squared,
        R2_adj = fit_sum$adj.r.squared,
        model_p_global = model_p,
        RMSE = rmse,
        RMSE_LOOCV = loocv$rmse,
        R2_LOOCV = loocv$r2,
        AIC = AIC(fit),
        AICc = calc_aicc(fit, nrow(model_df))
      ) %>% bind_cols(diagnostics),
      coefficients = coef_df %>%
        dplyr::transmute(Predictor, estimate = Estimate, p_value = `Pr(>|t|)`) %>%
        left_join(beta_tbl, by = "Predictor")
    )
  }

  fits <- lapply(candidate_sets_catch, function(preds) fit_candidate(site_df_in, preds))
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) {
    return(NULL)
  }

  aicc_vals <- vapply(fits, function(x) x$summary$AICc[1], numeric(1))
  best_idx <- which.min(aicc_vals)
  fits[[best_idx]]
}

catchment_results <- bind_rows(lapply(unique(mtt_definitions$mtt_definition), function(defn) {
  site_tbl <- site_df_base %>%
    mutate(
      basin_slope = if ("basin_slope" %in% names(.)) basin_slope else Slope_mean
    ) %>%
    dplyr::select(site, basin_slope, Harvest, Landslide_Total, Landslide_Young,
           Lava1_per, Lava2_per, Ash_Per, Pyro_per) %>%
    left_join(mtt_definitions %>% filter(mtt_definition == defn) %>% dplyr::select(site, MTT), by = "site")

  fit_obj <- fit_catchment_mtt(site_tbl)
  if (is.null(fit_obj)) {
    return(tibble(
      mtt_definition = defn,
      Predictors_Final = NA_character_,
      N = NA_integer_,
      R2 = NA_real_,
      R2_adj = NA_real_,
      model_p_global = NA_real_,
      RMSE = NA_real_,
      RMSE_LOOCV = NA_real_,
      R2_LOOCV = NA_real_,
      AIC = NA_real_,
      AICc = NA_real_,
      n_residuals = NA_real_,
      shapiro_W = NA_real_,
      shapiro_p = NA_real_,
      ncv_chisq = NA_real_,
      ncv_p = NA_real_,
      normality_pass_p05 = NA,
      homoscedasticity_pass_p05 = NA
    ))
  }
  fit_obj$summary %>% mutate(mtt_definition = defn, .before = 1)
}))

catchment_coefficients <- bind_rows(lapply(unique(mtt_definitions$mtt_definition), function(defn) {
  site_tbl <- site_df_base %>%
    mutate(
      basin_slope = if ("basin_slope" %in% names(.)) basin_slope else Slope_mean
    ) %>%
    dplyr::select(site, basin_slope, Harvest, Landslide_Total, Landslide_Young,
           Lava1_per, Lava2_per, Ash_Per, Pyro_per) %>%
    left_join(mtt_definitions %>% filter(mtt_definition == defn) %>% dplyr::select(site, MTT), by = "site")
  fit_obj <- fit_catchment_mtt(site_tbl)
  if (is.null(fit_obj)) return(NULL)
  fit_obj$coefficients %>% mutate(mtt_definition = defn, .before = 1)
}))

write_csv(catchment_results, file.path(output_dir, "mtt_catchment_model_summary.csv"))
write_csv(catchment_coefficients, file.path(output_dir, "mtt_catchment_model_coefficients.csv"))

fit_eco_models <- function(annual_df_in, include_mtt = TRUE) {
  mandatory_predictors <- "Pws"
  storage_predictors <- c("RBI", "RCS", "FDC", "SD", "WB", "CHS", "DR", "Fyw", "MTT")
  if (!include_mtt) {
    storage_predictors <- setdiff(storage_predictors, "MTT")
  }

  candidate_sets <- list(c("Pws"))
  storage_combos <- unlist(
    lapply(seq_along(storage_predictors), function(k) combn(storage_predictors, k, simplify = FALSE)),
    recursive = FALSE
  )
  candidate_sets <- c(candidate_sets, lapply(storage_combos, function(x) c("Pws", x)))

  fit_candidate <- function(df_in, response, predictors) {
    predictors <- unique(predictors[predictors %in% names(df_in)])
    if (!("Pws" %in% predictors)) {
      return(NULL)
    }

    model_df <- df_in %>%
      dplyr::select(any_of(c("site", "year", response, predictors))) %>%
      na.omit()

    if (nrow(model_df) < MODEL_MIN_N_ECO) {
      return(NULL)
    }

    pred_sd <- sapply(predictors, function(p) suppressWarnings(sd(model_df[[p]], na.rm = TRUE)))
    predictors <- predictors[is.finite(pred_sd) & pred_sd > 0]
    if (!("Pws" %in% predictors)) {
      return(NULL)
    }

    model_df <- model_df %>%
      mutate(across(all_of(predictors), ~ as.numeric(scale(.x))))

    fit <- tryCatch(
      lm(as.formula(paste(response, "~", paste(predictors, collapse = " + "))), data = model_df),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    fit <- apply_vif_filter(fit, response, model_df, VIF_THRESHOLD, protected = mandatory_predictors)
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
    fit_sum <- summary(fit)
    loocv <- calc_loocv_stats(formula(fit), model_df)
    diagnostics <- compute_residual_diagnostics(fit)

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

    beta_df_data <- model_df
    beta_df_data[[response]] <- as.numeric(scale(beta_df_data[[response]]))
    fit_beta <- lm(as.formula(paste(response, "~", paste(retained, collapse = " + "))), data = beta_df_data)
    beta_df <- tibble(
      Predictor = names(coef(fit_beta))[-1],
      Beta_Std = as.numeric(coef(fit_beta)[-1])
    )

    list(
      summary = tibble(
        Response = response,
        Predictors_Final = paste(retained, collapse = "; "),
        n = nrow(model_df),
        R2 = fit_sum$r.squared,
        R2_adj = fit_sum$adj.r.squared,
        model_p_global = model_p,
        RMSE = sqrt(mean(residuals(fit)^2, na.rm = TRUE)),
        RMSE_LOOCV = loocv$rmse,
        R2_LOOCV = loocv$r2,
        AIC = AIC(fit),
        AICc = calc_aicc(fit, nrow(model_df))
      ) %>% bind_cols(diagnostics),
      coefficients = coef_df %>%
        dplyr::transmute(Response = response, Predictor, estimate = Estimate, p_value = `Pr(>|t|)`) %>%
        left_join(beta_df, by = "Predictor")
    )
  }

  bind_rows(lapply(c("Q_7Q5", "T_7DMax"), function(response) {
    fits <- lapply(candidate_sets, function(preds) fit_candidate(annual_df_in, response, preds))
    fits <- Filter(Negate(is.null), fits)
    if (length(fits) == 0) return(NULL)
    best_idx <- which.min(vapply(fits, function(x) x$summary$AICc[1], numeric(1)))
    fits[[best_idx]]$summary
  }))
}

fit_eco_models_full <- function(annual_df_in, include_mtt = TRUE) {
  mandatory_predictors <- "Pws"
  storage_predictors <- c("RBI", "RCS", "FDC", "SD", "WB", "CHS", "DR", "Fyw", "MTT")
  if (!include_mtt) {
    storage_predictors <- setdiff(storage_predictors, "MTT")
  }

  candidate_sets <- list(c("Pws"))
  storage_combos <- unlist(
    lapply(seq_along(storage_predictors), function(k) combn(storage_predictors, k, simplify = FALSE)),
    recursive = FALSE
  )
  candidate_sets <- c(candidate_sets, lapply(storage_combos, function(x) c("Pws", x)))

  fit_candidate <- function(df_in, response, predictors) {
    predictors <- unique(predictors[predictors %in% names(df_in)])
    if (!("Pws" %in% predictors)) {
      return(NULL)
    }

    model_df <- df_in %>%
      dplyr::select(any_of(c("site", "year", response, predictors))) %>%
      na.omit()

    if (nrow(model_df) < MODEL_MIN_N_ECO) {
      return(NULL)
    }

    pred_sd <- sapply(predictors, function(p) suppressWarnings(sd(model_df[[p]], na.rm = TRUE)))
    predictors <- predictors[is.finite(pred_sd) & pred_sd > 0]
    if (!("Pws" %in% predictors)) {
      return(NULL)
    }

    model_df <- model_df %>%
      mutate(across(all_of(predictors), ~ as.numeric(scale(.x))))

    fit <- tryCatch(
      lm(as.formula(paste(response, "~", paste(predictors, collapse = " + "))), data = model_df),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    fit <- apply_vif_filter(fit, response, model_df, VIF_THRESHOLD, protected = mandatory_predictors)
    if (is.null(fit)) return(NULL)

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
    fit_sum <- summary(fit)
    loocv <- calc_loocv_stats(formula(fit), model_df)
    diagnostics <- compute_residual_diagnostics(fit)

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

    beta_df_data <- model_df
    beta_df_data[[response]] <- as.numeric(scale(beta_df_data[[response]]))
    fit_beta <- lm(as.formula(paste(response, "~", paste(retained, collapse = " + "))), data = beta_df_data)
    beta_df <- tibble(
      Predictor = names(coef(fit_beta))[-1],
      Beta_Std = as.numeric(coef(fit_beta)[-1])
    )

    list(
      summary = tibble(
        Response = response,
        Predictors_Final = paste(retained, collapse = "; "),
        n = nrow(model_df),
        R2 = fit_sum$r.squared,
        R2_adj = fit_sum$adj.r.squared,
        model_p_global = model_p,
        RMSE = sqrt(mean(residuals(fit)^2, na.rm = TRUE)),
        RMSE_LOOCV = loocv$rmse,
        R2_LOOCV = loocv$r2,
        AIC = AIC(fit),
        AICc = calc_aicc(fit, nrow(model_df))
      ) %>% bind_cols(diagnostics),
      coefficients = coef_df %>%
        dplyr::transmute(Response = response, Predictor, estimate = Estimate, p_value = `Pr(>|t|)`) %>%
        left_join(beta_df, by = "Predictor")
    )
  }

  summaries <- list()
  coefficients <- list()
  for (response in c("Q_7Q5", "T_7DMax")) {
    fits <- lapply(candidate_sets, function(preds) fit_candidate(annual_df_in, response, preds))
    fits <- Filter(Negate(is.null), fits)
    if (length(fits) == 0) next
    best_idx <- which.min(vapply(fits, function(x) x$summary$AICc[1], numeric(1)))
    summaries[[response]] <- fits[[best_idx]]$summary
    coefficients[[response]] <- fits[[best_idx]]$coefficients
  }

  list(
    summary = bind_rows(summaries),
    coefficients = bind_rows(coefficients)
  )
}

build_annual_model_df <- function(mtt_defn = NULL, include_mtt = TRUE) {
  iso_tbl <- raw_mtt %>%
    dplyr::select(site, Fyw) %>%
    left_join(damping, by = "site")

  if (include_mtt) {
    iso_tbl <- iso_tbl %>%
      left_join(
        mtt_definitions %>%
          filter(mtt_definition == mtt_defn) %>%
          dplyr::select(site, MTT),
        by = "site"
      )
  }

  annual_df %>%
    left_join(iso_tbl, by = "site")
}

eco_with_mtt <- lapply(unique(mtt_definitions$mtt_definition), function(defn) {
  fit_eco_models_full(build_annual_model_df(defn, include_mtt = TRUE), include_mtt = TRUE)
})
names(eco_with_mtt) <- unique(mtt_definitions$mtt_definition)

eco_without_mtt <- fit_eco_models_full(build_annual_model_df(include_mtt = FALSE), include_mtt = FALSE)

eco_summaries <- bind_rows(lapply(names(eco_with_mtt), function(defn) {
  eco_with_mtt[[defn]]$summary %>%
    mutate(mtt_definition = defn, model_variant = "with_mtt", .before = 1)
})) %>%
  bind_rows(
    eco_without_mtt$summary %>%
      mutate(mtt_definition = "none", model_variant = "without_mtt", .before = 1)
  )

eco_coefficients <- bind_rows(lapply(names(eco_with_mtt), function(defn) {
  eco_with_mtt[[defn]]$coefficients %>%
    mutate(mtt_definition = defn, model_variant = "with_mtt", .before = 1)
})) %>%
  bind_rows(
    eco_without_mtt$coefficients %>%
      mutate(mtt_definition = "none", model_variant = "without_mtt", .before = 1)
  )

write_csv(eco_summaries, file.path(output_dir, "mtt_eco_model_summary.csv"))
write_csv(eco_coefficients, file.path(output_dir, "mtt_eco_model_coefficients.csv"))

key_results <- eco_summaries %>%
  dplyr::select(mtt_definition, model_variant, Response, Predictors_Final, R2_adj, RMSE_LOOCV, R2_LOOCV, n) %>%
  arrange(Response, model_variant, mtt_definition)

write_csv(key_results, file.path(output_dir, "mtt_sensitivity_key_results.csv"))

safe_z <- function(x) {
  x_num <- as.numeric(x)
  s <- suppressWarnings(sd(x_num, na.rm = TRUE))
  if (!is.finite(s) || s <= 0) {
    return(rep(NA_real_, length(x_num)))
  }
  (x_num - mean(x_num, na.rm = TRUE)) / s
}

row_mean_min <- function(df_like, min_non_na = 1L) {
  mat <- as.matrix(df_like)
  n_non_na <- rowSums(is.finite(mat))
  out <- rowMeans(mat, na.rm = TRUE)
  out[n_non_na < min_non_na] <- NA_real_
  tibble(mean = out, n = n_non_na)
}

conceptual_axis_tbl <- bind_rows(lapply(unique(mtt_definitions$mtt_definition), function(defn) {
  site_tbl <- site_df_base %>%
    dplyr::select(site, RBI_mean, RCS_mean, FDC_mean, SD_mean, WB_mean, CHS_mean) %>%
    left_join(damping, by = "site") %>%
    left_join(raw_mtt %>% dplyr::select(site, Fyw), by = "site") %>%
    left_join(
      mtt_definitions %>% filter(mtt_definition == defn) %>% dplyr::select(site, MTT),
      by = "site"
    )

  metric_oriented <- site_tbl %>%
    dplyr::transmute(
      site,
      RBI_inv = -RBI_mean,
      RCS = RCS_mean,
      FDC = FDC_mean,
      SD = SD_mean,
      WB_depletion_mag = WB_mean,
      DR_inv = -DR,
      Fyw_inv = -Fyw,
      MTT = MTT,
      CHS = CHS_mean
    )

  metric_z <- metric_oriented %>%
    mutate(across(-site, safe_z, .names = "{.col}_z"))

  dynamic_vals <- row_mean_min(
    metric_z[, c("RBI_inv_z", "RCS_z", "FDC_z", "SD_z", "WB_depletion_mag_z")],
    min_non_na = 4L
  )

  mobile_vals <- row_mean_min(
    tibble(
      MTT = metric_z$MTT_z,
      DR = metric_z$DR_inv_z,
      Fyw = metric_z$Fyw_inv_z,
      CHS = metric_z$CHS_z
    ),
    min_non_na = 2L
  )

  tibble(
    mtt_definition = defn,
    site = metric_oriented$site,
    dynamic_storage_strength = dynamic_vals$mean,
    mobile_mixing_with_chs = mobile_vals$mean,
    n_mobile_components = mobile_vals$n,
    eligible_panel_b = is.finite(metric_oriented$RBI_inv) &
      is.finite(metric_oriented$RCS) &
      is.finite(metric_oriented$FDC) &
      is.finite(metric_oriented$SD) &
      is.finite(metric_oriented$WB_depletion_mag) &
      is.finite(metric_oriented$CHS) &
      is.finite(metric_oriented$DR_inv) &
      is.finite(metric_oriented$Fyw_inv) &
      is.finite(metric_oriented$MTT)
  )
}))

write_csv(conceptual_axis_tbl, file.path(output_dir, "mtt_conceptual_axis_sensitivity.csv"))

panel_b_summary <- conceptual_axis_tbl %>%
  filter(eligible_panel_b) %>%
  group_by(mtt_definition) %>%
  summarise(
    n_panel_b_sites = n(),
    sites_panel_b = paste(site, collapse = "; "),
    .groups = "drop"
  )

write_csv(panel_b_summary, file.path(output_dir, "mtt_conceptual_panel_b_summary.csv"))

scope_tbl <- tibble(
  analysis_component = c(
    "Dynamic/extended-dynamic PCA (Figure 3)",
    "ANOVA + Tukey HSD for annual site differences (Figure 2, BF in Figure 4a)",
    "Mobile storage descriptive MTT values (Figure 4d, Table 1)",
    "Dynamic-mobile correlation matrix (Figure 5)",
    "Catchment-control model for MTT (Figure 6, Table 4)",
    "Low-flow eco model Q7Q5 (Figure 7, Table 5)",
    "Temperature eco model T7DMax (Figure 7, Table 5)",
    "Conceptual framework mobile axis / panel b eligibility (Figure 9)"
  ),
  affected_by_mtt_definition = c(
    "No",
    "No",
    "Yes",
    "Yes",
    "Yes",
    "No",
    "Yes",
    "Yes"
  ),
  note = c(
    "PCA uses RBI, RCS, FDC, SD, and WB only.",
    "ANOVA/Tukey only use annual RBI, RCS, FDC, SD, WB, and CHS.",
    "Reported site-level MTT values depend directly on the chosen site-level definition.",
    "Correlations involving MTT change across early, late, and combined definitions.",
    "Selected predictor and model support for MTT vary across definitions.",
    "Selected Q7Q5 model is identical with early-only, late-only, combined, and no-MTT variants.",
    "Selected T7DMax model and fit vary across MTT definitions and when MTT is omitted.",
    "Mobile-axis values and the complete-case site set can change when MTT availability changes."
  )
)

write_csv(scope_tbl, file.path(output_dir, "mtt_sensitivity_scope.csv"))

supp_table <- eco_summaries %>%
  mutate(
    response = dplyr::case_when(
      Response == "Q_7Q5" ~ "Q7Q5",
      Response == "T_7DMax" ~ "T7DMax",
      TRUE ~ Response
    ),
    `MTT Definition` = label_mtt_definition(mtt_definition)
  ) %>%
  dplyr::select(
    `Section` = response,
    `MTT Definition`,
    `Selected Predictor(s)` = Predictors_Final,
    `n` = n,
    `Adj R2` = R2_adj,
    `Model p-value` = model_p_global,
    `LOOCV R2` = R2_LOOCV,
    `LOOCV RMSE` = RMSE_LOOCV
  ) %>%
  mutate(`Selected Predictor(s)` = clean_predictor_labels(`Selected Predictor(s)`)) %>%
  filter(`Section` == "T7DMax") %>%
  arrange(`Section`, `MTT Definition`)

catchment_supp_table <- catchment_results %>%
  mutate(
    `MTT Definition` = label_mtt_definition(mtt_definition)
  ) %>%
  dplyr::select(
    `Section` = mtt_definition,
    `MTT Definition`,
    `Selected Predictor(s)` = Predictors_Final,
    `n` = N,
    `Adj R2` = R2_adj,
    `Model p-value` = model_p_global,
    `LOOCV R2` = R2_LOOCV,
    `LOOCV RMSE` = RMSE_LOOCV
  ) %>%
  mutate(`Selected Predictor(s)` = clean_predictor_labels(`Selected Predictor(s)`)) %>%
  mutate(`Section` = "Catchment MTT model") %>%
  arrange(`MTT Definition`)

supp_table_combined <- bind_rows(
  supp_table,
  catchment_supp_table
) %>%
  transmute(
    `Response/model` = `Section`,
    `MTT definition` = `MTT Definition`,
    `Selected predictor(s)` = `Selected Predictor(s)`,
    `Sample size (n)` = n,
    `Adj R2` = `Adj R2`,
    `LOOCV R2` = `LOOCV R2`,
    `LOOCV RMSE` = `LOOCV RMSE`,
    `Model p-value` = `Model p-value`
  ) %>%
  mutate(
    `Adj R2` = signif(`Adj R2`, 9),
    `Model p-value` = signif(`Model p-value`, 9),
    `LOOCV R2` = signif(`LOOCV R2`, 9),
    `LOOCV RMSE` = signif(`LOOCV RMSE`, 9)
  ) %>%
  arrange(
    factor(`Response/model`, levels = c("T7DMax", "Catchment MTT model")),
    factor(
      `MTT definition`,
      levels = c(
        "MTT1 (McGuire et al., 2005; 2001-2003)",
        "MTT2 (Segura, 2015-2018)",
        "Combined MTT (mean of MTT1 and MTT2)",
        "No MTT"
      )
    )
  )

write_csv(supp_table_combined, file.path(output_dir, "TableS5_MTT_sensitivity.csv"))

write_csv(
  supp_table_combined,
  file.path(MS_TABLES_SUPP_DIR, "TableS5_MTT_sensitivity.csv")
)

summary_lines <- c(
  "# MTT Sensitivity Summary",
  "",
  "This analysis compares manuscript conclusions under alternative site-level MTT definitions.",
  "",
  "## MTT Definitions",
  paste0(
    "- ", mtt_availability$mtt_definition,
    ": n_sites=", mtt_availability$n_sites,
    ", mean=", signif(mtt_availability$mean_mtt, 3)
  ),
  "",
  "## Key Eco-Model Results",
  paste0(
    "- ", key_results$Response,
    " [", key_results$mtt_definition, "/", key_results$model_variant, "]: ",
    key_results$Predictors_Final,
    " | adj_R2=", signif(key_results$R2_adj, 3),
    " | LOOCV_R2=", signif(key_results$R2_LOOCV, 3)
  )
)

writeLines(summary_lines, file.path(output_dir, "README_MTT_sensitivity.md"))
