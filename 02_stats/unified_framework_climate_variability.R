# Quantify interannual climate variability around unified storage state.
# Inputs: OUTPUT_DIR/master/master_annual.csv; unified framework site axes.
# Author: Sidney Bush
# Date: 2026-02-21

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

rm(list = ls())

# Load project config
source("config.R")

find_read_root <- function() {
  primary <- OUTPUT_DIR
  fallback <- file.path(REPO_DIR, "outputs")
  if (file.exists(file.path(primary, "master", MASTER_ANNUAL_FILE))) {
    return(primary)
  }
  fallback
}

is_writable_root <- function(root_dir) {
  if (!dir.exists(root_dir)) {
    ok_dir <- tryCatch(
      {
        dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
        TRUE
      },
      error = function(e) FALSE
    )
    if (!ok_dir) return(FALSE)
  }
  probe <- file.path(root_dir, ".write_probe")
  ok <- tryCatch(
    {
      file.create(probe)
      file.remove(probe)
      TRUE
    },
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
  isTRUE(ok)
}

find_write_root <- function() {
  primary <- OUTPUT_DIR
  fallback <- file.path(REPO_DIR, "outputs")
  if (is_writable_root(primary)) return(primary)
  if (is_writable_root(fallback)) return(fallback)
  stop("No writable output root found for unified-framework climate analysis.")
}

find_axes_file <- function(read_root, write_root) {
  candidates <- c(
    file.path(read_root, "models", "unified_framework", "unified_framework_site_axes.csv"),
    file.path(write_root, "models", "unified_framework", "unified_framework_site_axes.csv"),
    file.path(REPO_DIR, "outputs", "models", "unified_framework", "unified_framework_site_axes.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Missing unified framework site axes file. Run 02_stats/unified_framework.R first.")
  }
  hit[1]
}

safe_z <- function(x) {
  x_num <- as.numeric(x)
  s <- suppressWarnings(sd(x_num, na.rm = TRUE))
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x_num)))
  (x_num - mean(x_num, na.rm = TRUE)) / s
}

row_mean_min <- function(df_like, min_non_na = 1L) {
  mat <- as.matrix(df_like)
  n_non_na <- rowSums(is.finite(mat))
  out <- rowMeans(mat, na.rm = TRUE)
  out[n_non_na < min_non_na] <- NA_real_
  list(mean = out, n = n_non_na)
}

calc_sen_slope <- function(x, y, min_pairs = 10L) {
  valid <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  n <- length(x)
  if (n < 3L) return(NA_real_)

  slopes <- numeric(0)
  for (i in seq_len(n - 1L)) {
    dx <- x[(i + 1L):n] - x[i]
    dy <- y[(i + 1L):n] - y[i]
    keep <- is.finite(dx) & is.finite(dy) & (dx != 0)
    if (any(keep)) {
      slopes <- c(slopes, dy[keep] / dx[keep])
    }
  }
  if (length(slopes) < min_pairs) return(NA_real_)
  stats::median(slopes, na.rm = TRUE)
}

extract_lm_summary <- function(df, formula_txt, min_n = 30L) {
  vars <- all.vars(stats::as.formula(formula_txt))
  dat <- df %>%
    select(all_of(vars)) %>%
    filter(if_all(everything(), is.finite))

  if (nrow(dat) < min_n) {
    return(list(
      model = NULL,
      n = nrow(dat),
      r2 = NA_real_,
      r2_adj = NA_real_,
      rmse = NA_real_,
      coefs = tibble(term = character(), estimate = numeric(), p_value = numeric())
    ))
  }

  fit <- tryCatch(stats::lm(stats::as.formula(formula_txt), data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(
      model = NULL,
      n = nrow(dat),
      r2 = NA_real_,
      r2_adj = NA_real_,
      rmse = NA_real_,
      coefs = tibble(term = character(), estimate = numeric(), p_value = numeric())
    ))
  }

  s <- summary(fit)
  coef_df <- as.data.frame(s$coefficients)
  coef_df$term <- rownames(coef_df)
  coefs <- coef_df %>%
    transmute(
      term = term,
      estimate = as.numeric(Estimate),
      p_value = as.numeric(`Pr(>|t|)`)
    )

  res <- stats::residuals(fit)
  rmse <- sqrt(mean(res^2, na.rm = TRUE))

  list(
    model = fit,
    n = nrow(dat),
    r2 = as.numeric(s$r.squared),
    r2_adj = as.numeric(s$adj.r.squared),
    rmse = rmse,
    coefs = coefs
  )
}

read_root <- find_read_root()
write_root <- find_write_root()
axes_file <- find_axes_file(read_root, write_root)

model_out_dir <- file.path(write_root, "models", "unified_framework")
table_out_dir <- file.path(write_root, "tables", "unified_framework")
for (d in c(model_out_dir, table_out_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

annual_file <- file.path(read_root, "master", MASTER_ANNUAL_FILE)
annual <- read_csv(annual_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site, year)

axes_site <- read_csv(axes_file, show_col_types = FALSE) %>%
  select(site, unified_state_index, dynamic_storage_strength_z, mobile_mixing_z, flow_path_partitioning_z)

annual2 <- annual %>%
  left_join(axes_site, by = "site") %>%
  mutate(
    WB_drawdown = -WB,
    RBI_inv = -RBI
  )

# Annual dynamic strength from annual hydrometric metrics.
annual_dyn <- annual2 %>%
  select(site, year, RBI_inv, RCS, FDC, SD, WB_drawdown)

annual_dyn_z <- annual_dyn %>%
  mutate(
    RBI_inv_z = safe_z(RBI_inv),
    RCS_z = safe_z(RCS),
    FDC_z = safe_z(FDC),
    SD_z = safe_z(SD),
    WB_drawdown_z = safe_z(WB_drawdown)
  )

dyn_strength_vals <- row_mean_min(
  annual_dyn_z %>% select(RBI_inv_z, RCS_z, FDC_z, SD_z, WB_drawdown_z),
  min_non_na = 4L
)

annual2 <- annual2 %>%
  mutate(dynamic_storage_strength_annual = dyn_strength_vals$mean)

# Annual dynamic-storage PCA (site-year space) for trajectory visualization.
pca_df <- annual_dyn_z %>%
  select(site, year, RBI_inv_z, RCS_z, FDC_z, SD_z, WB_drawdown_z) %>%
  mutate(
    RBI_inv_z = ifelse(is.na(RBI_inv_z), mean(RBI_inv_z, na.rm = TRUE), RBI_inv_z),
    RCS_z = ifelse(is.na(RCS_z), mean(RCS_z, na.rm = TRUE), RCS_z),
    FDC_z = ifelse(is.na(FDC_z), mean(FDC_z, na.rm = TRUE), FDC_z),
    SD_z = ifelse(is.na(SD_z), mean(SD_z, na.rm = TRUE), SD_z),
    WB_drawdown_z = ifelse(is.na(WB_drawdown_z), mean(WB_drawdown_z, na.rm = TRUE), WB_drawdown_z)
  )

pca_fit <- stats::prcomp(
  pca_df %>% select(RBI_inv_z, RCS_z, FDC_z, SD_z, WB_drawdown_z),
  center = TRUE,
  scale. = TRUE
)

pca_scores <- as.data.frame(pca_fit$x[, 1:2])
names(pca_scores) <- c("dynamic_pc1", "dynamic_pc2")
pca_scores <- bind_cols(
  pca_df %>% select(site, year),
  pca_scores
)

pca_loadings <- as.data.frame(pca_fit$rotation[, 1:2])
pca_loadings$metric <- rownames(pca_loadings)
pca_var <- tibble(
  pc = paste0("PC", seq_along(pca_fit$sdev)),
  variance_explained = (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
)

annual2 <- annual2 %>%
  left_join(pca_scores, by = c("site", "year"))

# Per-site annual anomalies around local climatology/state.
annual_anom <- annual2 %>%
  group_by(site) %>%
  mutate(
    WB_drawdown_anom = WB_drawdown - mean(WB_drawdown, na.rm = TRUE),
    P_NovJan_anom = P_NovJan - mean(P_NovJan, na.rm = TRUE),
    Q_7Q5_anom = Q_7Q5 - mean(Q_7Q5, na.rm = TRUE),
    T_Q7Q5_anom = T_Q7Q5 - mean(T_Q7Q5, na.rm = TRUE),
    T_7DMax_anom = T_7DMax - mean(T_7DMax, na.rm = TRUE),
    dynamic_storage_strength_annual_anom =
      dynamic_storage_strength_annual - mean(dynamic_storage_strength_annual, na.rm = TRUE)
  ) %>%
  ungroup()

# Main vs supplementary interaction models.
responses <- c("Q_7Q5_anom", "T_Q7Q5_anom", "T_7DMax_anom")

# Main model set: annual dynamic storage anomaly + winter input anomaly.
main_formula_by_response <- setNames(
  paste0(
    responses,
    " ~ dynamic_storage_strength_annual_anom * unified_state_index + P_NovJan_anom"
  ),
  responses
)

# Supplementary model set: WB-based annual depletion anomaly substitution.
supp_formula_by_response <- setNames(
  paste0(
    responses,
    " ~ WB_drawdown_anom * unified_state_index + P_NovJan_anom"
  ),
  responses
)

fit_formula_record <- function(df, response, formula_txt, model_set, model_id, min_n = 40L) {
  fit_out <- extract_lm_summary(df, formula_txt, min_n = min_n)
  model_row <- tibble(
    model_set = model_set,
    model_id = model_id,
    response = response,
    formula = formula_txt,
    n = fit_out$n,
    r2 = fit_out$r2,
    r2_adj = fit_out$r2_adj,
    rmse = fit_out$rmse
  )

  coef_rows <- fit_out$coefs %>%
    mutate(
      model_set = model_set,
      model_id = model_id,
      response = response,
      formula = formula_txt,
      n = fit_out$n,
      .before = 1
    )

  list(model_row = model_row, coef_rows = coef_rows)
}

model_rows <- list()
coef_rows <- list()
m_idx <- 1L

for (resp in responses) {
  main_fit <- fit_formula_record(
    annual_anom,
    response = resp,
    formula_txt = main_formula_by_response[[resp]],
    model_set = "main",
    model_id = "dynamic_plus_pwinter"
  )
  model_rows[[m_idx]] <- main_fit$model_row
  coef_rows[[m_idx]] <- main_fit$coef_rows
  m_idx <- m_idx + 1L

  supp_fit <- fit_formula_record(
    annual_anom,
    response = resp,
    formula_txt = supp_formula_by_response[[resp]],
    model_set = "supplementary",
    model_id = "wb_plus_pwinter"
  )
  model_rows[[m_idx]] <- supp_fit$model_row
  coef_rows[[m_idx]] <- supp_fit$coef_rows
  m_idx <- m_idx + 1L
}

interaction_models <- bind_rows(model_rows)
interaction_coefs <- bind_rows(coef_rows)

main_models <- interaction_models %>% filter(model_set == "main")
main_coefs <- interaction_coefs %>% filter(model_set == "main")
supp_models <- interaction_models %>% filter(model_set == "supplementary")
supp_coefs <- interaction_coefs %>% filter(model_set == "supplementary")

# Site-level slope extraction and relation to unified state.
# Keep both primary and supplementary predictor choices.
slope_rows <- list()
s_idx <- 1L
predictors <- c("dynamic_storage_strength_annual_anom", "P_NovJan_anom", "WB_drawdown_anom")
for (resp in responses) {
  for (pred in predictors) {
    for (st in SITE_ORDER_HYDROMETRIC) {
      d <- annual_anom %>%
        filter(site == st) %>%
        select(all_of(c(resp, pred))) %>%
        filter(if_all(everything(), is.finite))
      if (nrow(d) < 15L) next

      fit <- tryCatch(stats::lm(stats::as.formula(paste0(resp, " ~ ", pred)), data = d), error = function(e) NULL)
      if (is.null(fit)) next
      s <- summary(fit)
      if (!(pred %in% rownames(s$coefficients))) next

      slope_rows[[s_idx]] <- tibble(
        site = st,
        response = resp,
        predictor = pred,
        n = nrow(d),
        ols_slope = as.numeric(s$coefficients[pred, "Estimate"]),
        ols_slope_p = as.numeric(s$coefficients[pred, "Pr(>|t|)"]),
        sen_slope = calc_sen_slope(d[[pred]], d[[resp]], min_pairs = 10L),
        sen_tau = suppressWarnings(cor(d[[pred]], d[[resp]], method = "kendall", use = "pairwise.complete.obs")),
        sen_tau_p = tryCatch(
          stats::cor.test(d[[pred]], d[[resp]], method = "kendall")$p.value,
          error = function(e) NA_real_
        ),
        r2 = as.numeric(s$r.squared),
        r2_adj = as.numeric(s$adj.r.squared)
      )
      s_idx <- s_idx + 1L
    }
  }
}

site_slopes <- bind_rows(slope_rows) %>%
  left_join(axes_site %>% select(site, unified_state_index), by = "site") %>%
  mutate(
    slope = ols_slope,
    slope_p = ols_slope_p
  ) %>%
  mutate(
    predictor_group = case_when(
      predictor == "dynamic_storage_strength_annual_anom" ~ "main_storage_anomaly",
      predictor == "P_NovJan_anom" ~ "main_climate_input",
      predictor == "WB_drawdown_anom" ~ "supp_wb_depletion",
      TRUE ~ "other"
    )
  )

slope_state_corr_ols <- site_slopes %>%
  group_by(response, predictor, predictor_group) %>%
  summarise(
    n_sites = sum(is.finite(ols_slope) & is.finite(unified_state_index)),
    pearson_r = suppressWarnings(cor(ols_slope, unified_state_index, use = "pairwise.complete.obs", method = "pearson")),
    spearman_rho = suppressWarnings(cor(ols_slope, unified_state_index, use = "pairwise.complete.obs", method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(method = "ols", .before = 1)

slope_state_corr_sen <- site_slopes %>%
  group_by(response, predictor, predictor_group) %>%
  summarise(
    n_sites = sum(is.finite(sen_slope) & is.finite(unified_state_index)),
    pearson_r = suppressWarnings(cor(sen_slope, unified_state_index, use = "pairwise.complete.obs", method = "pearson")),
    spearman_rho = suppressWarnings(cor(sen_slope, unified_state_index, use = "pairwise.complete.obs", method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(method = "sen", .before = 1)

slope_state_corr <- bind_rows(slope_state_corr_ols, slope_state_corr_sen)

site_slopes_sen <- site_slopes %>%
  transmute(
    site,
    response,
    predictor,
    predictor_group,
    n,
    slope = sen_slope,
    slope_p = sen_tau_p,
    sen_tau,
    r2 = NA_real_,
    r2_adj = NA_real_,
    unified_state_index
  )

site_slopes_ols <- site_slopes %>%
  transmute(
    site,
    response,
    predictor,
    predictor_group,
    n,
    slope = ols_slope,
    slope_p = ols_slope_p,
    sen_tau = NA_real_,
    r2,
    r2_adj,
    unified_state_index
  )

slope_method_comparison <- site_slopes %>%
  group_by(response, predictor, predictor_group) %>%
  summarise(
    n_sites = sum(is.finite(ols_slope) & is.finite(sen_slope)),
    spearman_ols_vs_sen = suppressWarnings(cor(ols_slope, sen_slope, use = "pairwise.complete.obs", method = "spearman")),
    pearson_ols_vs_sen = suppressWarnings(cor(ols_slope, sen_slope, use = "pairwise.complete.obs", method = "pearson")),
    mean_abs_diff = mean(abs(ols_slope - sen_slope), na.rm = TRUE),
    .groups = "drop"
  )

# Save outputs
write_csv(
  annual_anom,
  file.path(model_out_dir, "unified_framework_annual_anomalies.csv")
)
write_csv(
  pca_loadings,
  file.path(model_out_dir, "unified_framework_dynamic_pca_loadings.csv")
)
write_csv(
  pca_var,
  file.path(model_out_dir, "unified_framework_dynamic_pca_variance.csv")
)
write_csv(
  main_models,
  file.path(model_out_dir, "unified_framework_climate_main_models.csv")
)
write_csv(
  main_coefs,
  file.path(model_out_dir, "unified_framework_climate_main_coefficients.csv")
)
write_csv(
  supp_models,
  file.path(model_out_dir, "unified_framework_climate_supp_wb_models.csv")
)
write_csv(
  supp_coefs,
  file.path(model_out_dir, "unified_framework_climate_supp_wb_coefficients.csv")
)
write_csv(
  interaction_models,
  file.path(model_out_dir, "unified_framework_climate_interaction_models.csv")
)
write_csv(
  interaction_coefs,
  file.path(model_out_dir, "unified_framework_climate_interaction_coefficients.csv")
)
write_csv(
  site_slopes,
  file.path(model_out_dir, "unified_framework_site_sensitivity_slopes.csv")
)
write_csv(
  site_slopes_ols,
  file.path(model_out_dir, "unified_framework_site_sensitivity_slopes_ols.csv")
)
write_csv(
  site_slopes_sen,
  file.path(model_out_dir, "unified_framework_site_sensitivity_slopes_sen.csv")
)
write_csv(
  slope_method_comparison,
  file.path(model_out_dir, "unified_framework_slope_method_comparison.csv")
)
write_csv(
  slope_state_corr_ols,
  file.path(model_out_dir, "unified_framework_slope_state_correlations_ols.csv")
)
write_csv(
  slope_state_corr_sen,
  file.path(model_out_dir, "unified_framework_slope_state_correlations_sen.csv")
)
write_csv(
  slope_state_corr,
  file.path(model_out_dir, "unified_framework_slope_state_correlations.csv")
)

write_csv(
  main_models,
  file.path(table_out_dir, "unified_framework_climate_main_models.csv")
)
write_csv(
  main_coefs,
  file.path(table_out_dir, "unified_framework_climate_main_coefficients.csv")
)
write_csv(
  supp_models,
  file.path(table_out_dir, "unified_framework_climate_supp_wb_models.csv")
)
write_csv(
  supp_coefs,
  file.path(table_out_dir, "unified_framework_climate_supp_wb_coefficients.csv")
)
write_csv(
  interaction_models,
  file.path(table_out_dir, "unified_framework_climate_interaction_models.csv")
)
write_csv(
  interaction_coefs,
  file.path(table_out_dir, "unified_framework_climate_interaction_coefficients.csv")
)
write_csv(
  site_slopes,
  file.path(table_out_dir, "unified_framework_site_sensitivity_slopes.csv")
)
write_csv(
  site_slopes_ols,
  file.path(table_out_dir, "unified_framework_site_sensitivity_slopes_ols.csv")
)
write_csv(
  site_slopes_sen,
  file.path(table_out_dir, "unified_framework_site_sensitivity_slopes_sen.csv")
)
write_csv(
  slope_method_comparison,
  file.path(table_out_dir, "unified_framework_slope_method_comparison.csv")
)
write_csv(
  slope_state_corr_ols,
  file.path(table_out_dir, "unified_framework_slope_state_correlations_ols.csv")
)
write_csv(
  slope_state_corr_sen,
  file.path(table_out_dir, "unified_framework_slope_state_correlations_sen.csv")
)
write_csv(
  slope_state_corr,
  file.path(table_out_dir, "unified_framework_slope_state_correlations.csv")
)
