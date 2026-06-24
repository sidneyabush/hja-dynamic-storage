# test whether alternate MTT definitions change the main model conclusions

# inputs:
# outputs/master/master_annual.csv
# outputs/master/master_site.csv
# isotope_dir/MTT_FYW.csv
# isotope_dir/DampingRatios_2025-07-07.csv
# outputs/metrics/mobile/isotope_metrics_site.csv
# outputs/models/storage_eco_response_mlr/storage_eco_response_mlr_summary.csv

# outputs:
# outputs/MTT_sensitivity/TableS7_MTT_sensitivity.csv

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(dplyr, readr, tidyr, tibble, MASS, car, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

# write the sensitivity table separately from the main model outputs
output_dir <- file.path(OUTPUT_DIR, "MTT_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

VIF_THRESHOLD <- 10
MODEL_MIN_N_ECO <- 20
MIN_N_CATCH <- 5

annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
isotope_file <- file.path(ISOTOPE_DIR, "MTT_FYW.csv")
damping_file <- file.path(ISOTOPE_DIR, "DampingRatios_2025-07-07.csv")
main_eco_summary_file <- file.path(OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR, "storage_eco_response_mlr_summary.csv")

# load annual site year values and site level catchment predictors
annual_df <- read_csv(annual_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

site_df_base <- read_csv(site_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC, !site %in% SITE_EXCLUDE_STANDARD)

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

# damping ratio and Fyw stay fixed while MTT definitions are varied
damping_raw <- read_csv(damping_file, show_col_types = FALSE)
dr_col <- dplyr::case_when(
  "DR_Overall" %in% names(damping_raw) ~ "DR_Overall",
  "DR" %in% names(damping_raw) ~ "DR",
  TRUE ~ NA_character_
)

# stop if the damping ratio source table changed column names
if (is.na(dr_col)) {
  stop("No damping-ratio column found in: ", damping_file)
}
damping <- damping_raw %>%
  mutate(
    site = if ("site" %in% names(.)) site else SITECODE,
    site = standardize_site_code(site),
    DR = suppressWarnings(as.numeric(.data[[dr_col]]))
  ) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC, !is.na(site), site != "") %>%
  group_by(site) %>%
  summarise(
    DR = ifelse(any(is.finite(DR)), mean(DR, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )

# build MTT1, MTT2, and combined MTT alternatives from the source isotope table
raw_mtt <- read_csv(isotope_file, show_col_types = FALSE) %>%
  mutate(
    site = standardize_site_code(trimws(site)),
    MTT_early = suppressWarnings(as.numeric(MTT1)),
    MTT_second_period = rowMeans(cbind(
      suppressWarnings(as.numeric(MTT2L)),
      suppressWarnings(as.numeric(MTT2H))
    ), na.rm = TRUE
    )
  ) %>%
  mutate(MTT_second_period = ifelse(is.nan(MTT_second_period), NA_real_, MTT_second_period)) %>%
  dplyr::select(site, MTT_early, MTT_second_period) %>%
  filter(!is.na(site), site != "", site %in% SITE_ORDER_HYDROMETRIC)

mtt_definitions <- raw_mtt %>%
  dplyr::transmute(
    site,
    early_only = MTT_early,
    mtt2_only = MTT_second_period,
    combined_mean = NA_real_
  ) %>%
  left_join(
    isotope_site_mean %>% dplyr::select(site, MTT_site_mean),
    by = "site"
  ) %>%
  mutate(combined_mean = MTT_site_mean) %>%
  pivot_longer(
    cols = c(early_only, mtt2_only, combined_mean),
    names_to = "mtt_definition",
    values_to = "MTT"
  ) %>%
  arrange(mtt_definition, site)

label_mtt_definition <- function(x) {
  dplyr::case_when(
    x == "early_only" ~ "MTT1 (McGuire et al., 2005; 2001-2003)",
    x == "mtt2_only" ~ "MTT2 (Segura, 2015-2018)",
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

# catchment MTT models use the same predictor search as the main catchment models
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
  # fit one catchment predictor set
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
        Predictors_Final = paste(retained, collapse = ", "),
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

# run the catchment MTT model for each alternate MTT definition
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

fit_eco_models_full <- function(annual_df_in, include_mtt = TRUE) {
  # ecological response models always keep wet season precipitation
  mandatory_predictors <- "Pws"
  storage_predictors <- c("RBI", "RCS", "FDC", "SD", "WB", "BF", "DR", "Fyw", "MTT")
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
    # fit one ecological response predictor set
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

    # calculate the overall model p value from the fitted F statistic
    fstat <- suppressWarnings(as.numeric(fit_sum$fstatistic))
    model_p <- if (!is.null(fstat) && length(fstat) >= 3 && all(is.finite(fstat[1:3]))) {
      pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
    } else {
      NA_real_
    }

    # use standardized coefficients so Table S7 matches the main model tables
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

    # store the AICc summary and selected coefficients for each MTT definition
    list(
      summary = tibble(
        Response = response,
        Predictors_Final = paste(retained, collapse = ", "),
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
  # add site isotope metrics to each annual site year record
  iso_tbl <- isotope_site_mean %>%
    dplyr::select(site, DR_site_mean, Fyw_site_mean)

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
    left_join(iso_tbl, by = "site") %>%
    mutate(
      DR = DR_site_mean,
      Fyw = Fyw_site_mean
    ) %>%
    { if (include_mtt) mutate(., MTT = MTT) else . } %>%
  dplyr::select(-DR_site_mean, -Fyw_site_mean)
}

# fit ecological response models with alternate MTT definitions and without MTT
eco_mtt_variants <- setdiff(unique(mtt_definitions$mtt_definition), "combined_mean")

eco_with_mtt <- lapply(eco_mtt_variants, function(defn) {
  fit_eco_models_full(build_annual_model_df(defn, include_mtt = TRUE), include_mtt = TRUE)
})
names(eco_with_mtt) <- eco_mtt_variants

eco_without_mtt <- fit_eco_models_full(build_annual_model_df(include_mtt = FALSE), include_mtt = FALSE)

eco_summaries <- bind_rows(lapply(names(eco_with_mtt), function(defn) {
  eco_with_mtt[[defn]]$summary %>%
    mutate(mtt_definition = defn, model_variant = "with_mtt", .before = 1)
})) %>%
  bind_rows(
    eco_without_mtt$summary %>%
      mutate(mtt_definition = "none", model_variant = "without_mtt", .before = 1)
  )

main_eco_summary <- read_csv(main_eco_summary_file, show_col_types = FALSE) %>%
  mutate(Response = as.character(Response)) %>%
  filter(Response == "T7DMax") %>%
  mutate(
    mtt_definition = "combined_mean",
    model_variant = "with_mtt"
  ) %>%
  dplyr::select(
    mtt_definition,
    model_variant,
    Response,
    Predictors_Final,
    n,
    R2,
    R2_adj,
    model_p_global,
    RMSE,
    RMSE_LOOCV,
    R2_LOOCV,
    AIC,
    AICc
  )

# use the already selected main T7DMax model for the combined MTT case
eco_summaries <- eco_summaries %>%
  filter(!(mtt_definition == "combined_mean" & Response == "T_7DMax")) %>%
  bind_rows(main_eco_summary)

# assemble the final Table S7 rows and keep only T7DMax sensitivity results
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
    `LOOCV R2` = R2_LOOCV,
    `LOOCV RMSE` = RMSE_LOOCV,
    `Model p-value` = model_p_global
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
    `LOOCV R2` = R2_LOOCV,
    `LOOCV RMSE` = RMSE_LOOCV,
    `Model p-value` = model_p_global
  ) %>%
  mutate(`Selected Predictor(s)` = clean_predictor_labels(`Selected Predictor(s)`)) %>%
  mutate(`Section` = "Catchment MTT model") %>%
  arrange(`MTT Definition`)

# combine ecological response and catchment MTT sensitivity results for Table S7
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

write_csv(supp_table_combined, file.path(output_dir, "TableS7_MTT_sensitivity.csv"))

try(grDevices::graphics.off(), silent = TRUE)
quit(save = "no", status = 0, runLast = FALSE)
