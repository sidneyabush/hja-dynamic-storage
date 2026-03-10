# calculate mobile-storage metrics (chs + isotope metrics) by site and water year.
# inputs: discharge_dir/hf00402_v14.csv; ec_dir/cf01201_v4.txt; ec_dir/cf00201_v7.csv; isotope_dir/mtt_fyw.csv; isotope_dir/dampingratios_2025-07-07.csv.
# author: keira johnson (original chs), sidney bush
# date: 2026-03-09

library(dplyr)
library(readr)
library(lubridate)

rm(list = ls())

# load project config
source("config.R")

output_dir <- OUT_MET_MOBILE_DIR
discharge_dir <- DISCHARGE_DIR
ec_dir <- EC_DIR
isotope_dir <- ISOTOPE_DIR

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- part 1: chs from chemistry tracers and discharge ----

discharge <- read.csv(file.path(discharge_dir, "HF00402_v14.csv")) %>%
  mutate(
    date = as.Date(DATE, "%m/%d/%Y"),
    SITECODE = standardize_site_code(SITECODE)
  ) %>%
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    WATERYEAR >= WY_START,
    WATERYEAR <= WY_END
  ) %>%
  select(SITECODE, date, MEAN_Q, WATERYEAR)

calc_chs_from_tracer <- function(tracer_daily, discharge_tbl, min_obs_per_wy) {
  tracer_q <- tracer_daily %>%
    left_join(discharge_tbl %>% select(SITECODE, date, MEAN_Q), by = c("SITECODE", "date")) %>%
    filter(is.finite(tracer_value), is.finite(MEAN_Q), MEAN_Q > 0) %>%
    mutate(waterYear = get_water_year(date)) %>%
    filter(waterYear >= WY_START, waterYear <= WY_END)

  goodyears <- tracer_q %>%
    group_by(SITECODE, waterYear) %>%
    summarise(num_obs = n_distinct(date), .groups = "drop") %>%
    filter(num_obs >= min_obs_per_wy)

  tracer_q_kept <- tracer_q %>%
    semi_join(goodyears, by = c("SITECODE", "waterYear"))

  if (nrow(tracer_q_kept) == 0) {
    annual <- tibble(
      SITECODE = character(),
      waterYear = integer(),
      CHS = numeric(),
      median_bf = numeric(),
      sd_bf = numeric(),
      n_obs = integer()
    )
    return(list(annual = annual, goodyears = goodyears))
  }

  tracer_q_kept <- tracer_q_kept %>%
    group_by(SITECODE) %>%
    mutate(
      tracer_runoff = quantile(tracer_value, 0.01, na.rm = TRUE),
      tracer_groundwater = quantile(tracer_value, 0.99, na.rm = TRUE),
      tracer_span = tracer_groundwater - tracer_runoff
    ) %>%
    ungroup() %>%
    mutate(
      Q_baseflow = ifelse(
        is.finite(tracer_span) & tracer_span > 0,
        MEAN_Q * (tracer_value - tracer_runoff) / tracer_span,
        NA_real_
      ),
      GW_prop = Q_baseflow / MEAN_Q,
      GW_prop = pmax(0, pmin(1, GW_prop))
    )

  annual <- tracer_q_kept %>%
    group_by(SITECODE, waterYear) %>%
    summarise(
      CHS = mean(GW_prop, na.rm = TRUE),
      median_bf = median(GW_prop, na.rm = TRUE),
      sd_bf = sd(GW_prop, na.rm = TRUE),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(
      CHS = ifelse(is.nan(CHS), NA_real_, CHS),
      median_bf = ifelse(is.nan(median_bf), NA_real_, median_bf),
      sd_bf = ifelse(is.nan(sd_bf), NA_real_, sd_bf)
    )

  list(annual = annual, goodyears = goodyears)
}

ec_file <- file.path(ec_dir, "CF01201_v4.txt")
if (!file.exists(ec_file)) {
  stop("Missing required EC file: ", ec_file)
}
ec_inst <- read_csv(ec_file, show_col_types = FALSE)

ec_daily <- ec_inst %>%
  mutate(
    date = as.Date(DATE_TIME),
    SITECODE = standardize_site_code(SITECODE),
    tracer_value = suppressWarnings(as.numeric(EC_INST))
  ) %>%
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    is.finite(tracer_value)
  ) %>%
  group_by(SITECODE, date) %>%
  summarise(tracer_value = mean(tracer_value, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(tracer_value))

ec_chs <- calc_chs_from_tracer(ec_daily, discharge, CHS_MIN_DAYS_PER_WY)
annual_bf_prop <- ec_chs$annual %>%
  rename(n_days = n_obs)

write.csv(
  annual_bf_prop,
  file.path(output_dir, "annual_gw_prop.csv"),
  row.names = FALSE
)

write.csv(
  ec_chs$goodyears %>% arrange(SITECODE, waterYear),
  file.path(OUT_MET_SUPPORT_DIR, "chs_wy_day_counts_kept.csv"),
  row.names = FALSE
)

chem_file <- file.path(ec_dir, "CF00201_v7.csv")
if (!file.exists(chem_file)) {
  stop("Missing required chemistry file: ", chem_file)
}
chem <- read_csv(chem_file, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(DATE_TIME),
    SITECODE = standardize_site_code(SITECODE),
    COND = suppressWarnings(as.numeric(COND)),
    CA = suppressWarnings(as.numeric(CA)),
    waterYear = get_water_year(date)
  ) %>%
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    waterYear >= WY_START,
    waterYear <= WY_END
  )

chem_ec_daily <- chem %>%
  group_by(SITECODE, date) %>%
  summarise(tracer_value = mean(COND, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(tracer_value))

chem_ca_daily <- chem %>%
  group_by(SITECODE, date) %>%
  summarise(tracer_value = mean(CA, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(tracer_value))

ec_chem_chs <- calc_chs_from_tracer(chem_ec_daily, discharge, CHS_MIN_OBS_PER_WY_CHEM)
ca_chs <- calc_chs_from_tracer(chem_ca_daily, discharge, CHS_MIN_OBS_PER_WY_CHEM)

annual_bf_prop_ec_chem <- ec_chem_chs$annual
annual_bf_prop_ca <- ca_chs$annual

write.csv(
  annual_bf_prop_ec_chem,
  file.path(output_dir, "annual_gw_prop_ec_chem.csv"),
  row.names = FALSE
)

write.csv(
  annual_bf_prop_ca,
  file.path(output_dir, "annual_gw_prop_ca.csv"),
  row.names = FALSE
)

write.csv(
  ec_chem_chs$goodyears %>% arrange(SITECODE, waterYear),
  file.path(OUT_MET_SUPPORT_DIR, "chs_wy_obs_counts_ec_chem_kept.csv"),
  row.names = FALSE
)

write.csv(
  ca_chs$goodyears %>% arrange(SITECODE, waterYear),
  file.path(OUT_MET_SUPPORT_DIR, "chs_wy_obs_counts_ca_kept.csv"),
  row.names = FALSE
)

comparison <- annual_bf_prop %>%
  select(
    SITECODE,
    waterYear,
    CHS_EC = CHS,
    n_days_ec = n_days
  ) %>%
  full_join(
    annual_bf_prop_ec_chem %>%
      select(
        SITECODE,
        waterYear,
        CHS_EC_CHEM = CHS,
        n_obs_ec_chem = n_obs
      ),
    by = c("SITECODE", "waterYear")
  ) %>%
  full_join(
    annual_bf_prop_ca %>%
      select(
        SITECODE,
        waterYear,
        CHS_CA = CHS,
        n_obs_ca = n_obs
      ),
    by = c("SITECODE", "waterYear")
  ) %>%
  mutate(
    diff_ca_minus_ec = CHS_CA - CHS_EC,
    diff_ca_minus_ec_chem = CHS_CA - CHS_EC_CHEM
  ) %>%
  arrange(SITECODE, waterYear)

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  cor(x[keep], y[keep])
}

safe_cor_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  cor(x[keep], y[keep], method = "spearman")
}

safe_mean <- function(x) {
  if (!any(is.finite(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

safe_rmse <- function(x, y) {
  d <- x - y
  keep <- is.finite(d)
  if (sum(keep) == 0) {
    return(NA_real_)
  }
  sqrt(mean(d[keep]^2))
}

safe_lm_intercept <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  coef(lm(y[keep] ~ x[keep]))[1]
}

safe_lm_slope <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  coef(lm(y[keep] ~ x[keep]))[2]
}

site_summary <- comparison %>%
  group_by(SITECODE) %>%
  summarise(
    n_years_ec_ca = sum(is.finite(CHS_EC) & is.finite(CHS_CA)),
    corr_ec_vs_ca = safe_cor(CHS_EC, CHS_CA),
    mean_diff_ca_minus_ec = safe_mean(diff_ca_minus_ec),
    rmse_ca_vs_ec = safe_rmse(CHS_CA, CHS_EC),
    n_years_ec_chem_ca = sum(is.finite(CHS_EC_CHEM) & is.finite(CHS_CA)),
    corr_ec_chem_vs_ca = safe_cor(CHS_EC_CHEM, CHS_CA),
    mean_diff_ca_minus_ec_chem = safe_mean(diff_ca_minus_ec_chem),
    rmse_ca_vs_ec_chem = safe_rmse(CHS_CA, CHS_EC_CHEM),
    .groups = "drop"
  )

overall_n_years_ec_ca <- sum(is.finite(comparison$CHS_EC) & is.finite(comparison$CHS_CA))
overall_n_years_ec_chem_ca <- sum(is.finite(comparison$CHS_EC_CHEM) & is.finite(comparison$CHS_CA))

overall_summary <- tibble(
  n_years_ec_ca = overall_n_years_ec_ca,
  corr_ec_vs_ca = safe_cor(comparison$CHS_EC, comparison$CHS_CA),
  mean_diff_ca_minus_ec = safe_mean(comparison$diff_ca_minus_ec),
  rmse_ca_vs_ec = safe_rmse(comparison$CHS_CA, comparison$CHS_EC),
  n_years_ec_chem_ca = overall_n_years_ec_chem_ca,
  corr_ec_chem_vs_ca = safe_cor(comparison$CHS_EC_CHEM, comparison$CHS_CA),
  mean_diff_ca_minus_ec_chem = safe_mean(comparison$diff_ca_minus_ec_chem),
  rmse_ca_vs_ec_chem = safe_rmse(comparison$CHS_CA, comparison$CHS_EC_CHEM)
)

site_year_pairs_ec_ca <- comparison %>%
  filter(is.finite(CHS_EC_CHEM), is.finite(CHS_CA)) %>%
  transmute(
    SITECODE,
    waterYear,
    CHS_EC = CHS_EC_CHEM,
    CHS_CA,
    diff_ca_minus_ec = CHS_CA - CHS_EC_CHEM,
    abs_diff_ca_minus_ec = abs(CHS_CA - CHS_EC_CHEM),
    pct_diff_ca_vs_ec = ifelse(
      abs(CHS_EC_CHEM) > 0,
      100 * (CHS_CA - CHS_EC_CHEM) / CHS_EC_CHEM,
      NA_real_
    )
  ) %>%
  arrange(SITECODE, waterYear)

site_by_site_stats <- site_year_pairs_ec_ca %>%
  group_by(SITECODE) %>%
  summarise(
    n_years_paired = n(),
    corr_pearson = safe_cor(CHS_EC, CHS_CA),
    corr_spearman = safe_cor_spearman(CHS_EC, CHS_CA),
    mean_diff_ca_minus_ec = safe_mean(diff_ca_minus_ec),
    median_diff_ca_minus_ec = median(diff_ca_minus_ec, na.rm = TRUE),
    rmse_ca_vs_ec = safe_rmse(CHS_CA, CHS_EC),
    lm_intercept_ca_on_ec = safe_lm_intercept(CHS_EC, CHS_CA),
    lm_slope_ca_on_ec = safe_lm_slope(CHS_EC, CHS_CA),
    .groups = "drop"
  ) %>%
  arrange(SITECODE)

write.csv(
  comparison,
  file.path(output_dir, "annual_gw_prop_ec_ca_comparison.csv"),
  row.names = FALSE
)

write.csv(
  site_summary,
  file.path(output_dir, "annual_gw_prop_ec_ca_site_summary.csv"),
  row.names = FALSE
)

write.csv(
  overall_summary,
  file.path(output_dir, "annual_gw_prop_ec_ca_overall_summary.csv"),
  row.names = FALSE
)

write.csv(
  site_year_pairs_ec_ca,
  file.path(output_dir, "annual_gw_prop_ec_ca_site_year_pairs.csv"),
  row.names = FALSE
)

write.csv(
  site_by_site_stats,
  file.path(output_dir, "annual_gw_prop_ec_ca_site_by_site_stats.csv"),
  row.names = FALSE
)

# ---- part 2: isotope-derived site metrics (mtt/fyw/dr) ----

mtt_fyw <- read_csv(
  file.path(isotope_dir, "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  mutate(
    MTT2L_val = if ("MTT2L" %in% names(.)) suppressWarnings(as.numeric(MTT2L)) else NA_real_,
    MTT2H_val = if ("MTT2H" %in% names(.)) suppressWarnings(as.numeric(MTT2H)) else NA_real_,
    MTT2M_val = if ("MTT2M" %in% names(.)) suppressWarnings(as.numeric(MTT2M)) else NA_real_
  ) %>%
  mutate(
    MTT1 = suppressWarnings(as.numeric(MTT1)),
    MTT2 = suppressWarnings(as.numeric(dplyr::coalesce(
      MTT2M_val,
      rowMeans(cbind(
        MTT2L_val,
        MTT2H_val
      ), na.rm = TRUE)
    ))),
    MTT2 = ifelse(is.nan(MTT2), NA_real_, MTT2),
    Fyw = suppressWarnings(as.numeric(FYWM))
  ) %>%
  select(site, MTT1, MTT2, Fyw) %>%
  filter(!is.na(site), site != "")

damping <- read_csv(
  file.path(isotope_dir, "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  select(site, DR = DR_Overall, DR_err = DR__err)

isotope_metrics <- mtt_fyw %>%
  full_join(damping, by = "site") %>%
  mutate(site = factor(site, levels = SITE_ORDER_ALL)) %>%
  arrange(site)

write.csv(
  isotope_metrics,
  file.path(output_dir, "isotope_metrics_site.csv"),
  row.names = FALSE
)
