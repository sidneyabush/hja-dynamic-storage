# Assemble all annual/site-level metrics into master analysis tables.
# Inputs: dynamic_dir/RBI_RecessionCurve_Annual.csv; dynamic_dir/StorageDischarge_FDC_Annual.csv; mobile_dir/Annual_GW_Prop.csv; extended_dir/DS_drawdown_annual.csv; eco_dir/stream_thermal_lowflow_metrics_annual.csv; isotope_dir/MTT_FYW.csv; +1 more CSV files.
# Author: Sidney Bush
# Date: 2026-02-13

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(GGally)
library(ggcorrplot)

# Start clean
rm(list = ls())

# Load config (works from Rscript or source)
# Load project config
source("config.R")


theme_set(theme_pub(base_size = 12))

# Key input/output directories

base_dir    <- BASE_DATA_DIR
dynamic_dir <- OUT_MET_DYNAMIC_DIR
mobile_dir  <- OUT_MET_MOBILE_DIR
extended_dir <- OUT_MET_EXTENDED_DIR
eco_dir <- OUT_MET_ECO_DIR
master_dir <- OUT_MASTER_DIR
isotope_dir <- ISOTOPE_DIR
catchment_dir <- CATCHMENT_CHARACTERISTICS_DIR

# Make sure output directory exists
if (!dir.exists(master_dir)) dir.create(master_dir, recursive = TRUE)

# Helper to fail fast if a join key is not unique

assert_unique_keys <- function(df, keys, df_name) {
  dupes <- df %>%
    count(across(all_of(keys)), name = "n") %>%
    filter(n > 1)
  if (nrow(dupes) > 0) {
    stop(
      paste0(
        "Non-unique join keys detected in ", df_name, " for keys (",
        paste(keys, collapse = ", "), ")."
      )
    )
  }
}

# Dynamic storage metrics (annual)

# RBI and recession slope
rbi_path <- file.path(dynamic_dir, "rbi_rcs_annual.csv")
rbi_recession <- read_csv(
  rbi_path,
  show_col_types = FALSE
) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, RCS, RBI)
assert_unique_keys(rbi_recession, c("site", "year"), "rbi_recession")

# Storage-discharge, FDC, and supporting flow quantiles
fdc_path <- file.path(dynamic_dir, "storage_discharge_fdc_annual.csv")
storage_fdc <- read_csv(
  fdc_path,
  show_col_types = FALSE
) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, SD, FDC, Q99, Q50, Q01, Q5norm, CV_Q5norm)
assert_unique_keys(storage_fdc, c("site", "year"), "storage_fdc")

# Mobile storage metric (annual CHS)

# CHS = annual mean baseflow fraction
chs_path <- file.path(mobile_dir, "annual_gw_prop.csv")
baseflow <- read_csv(
  chs_path,
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, CHS)
assert_unique_keys(baseflow, c("site", "year"), "baseflow")

# Extended dynamic storage metric (annual WB)

# Water-balance drawdown
wb_path <- file.path(extended_dir, "ds_drawdown_annual.csv")
wb_storage <- read_csv(
  wb_path,
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, WB)
assert_unique_keys(wb_storage, c("site", "year"), "wb_storage")

# Ecological response variables (annual)

# Thermal and low-flow response metrics
thermal_lowflow <- read_csv(
  file.path(eco_dir, "stream_thermal_lowflow_metrics_annual.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = if ("site" %in% names(.)) site else SITECODE,
    year = if ("year" %in% names(.)) year else wateryear
  ) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(
    site, year, T_7DMax, Q_7Q5, T_at_Q7Q5, T_Q7Q5,
    max_temp_7d_C,
    q5_7d_mm_d, temp_at_q5_7d_C, temp_during_q5_7d_C,
    min_Q_7d_mm_d, temp_at_min_Q_7d_C, temp_during_min_Q_7d_C,
    Q5_CV
  ) %>%
  mutate(
    T_7DMax = ifelse(is.na(T_7DMax), max_temp_7d_C, T_7DMax),
    Q_7Q5 = ifelse(is.na(Q_7Q5), q5_7d_mm_d, Q_7Q5),
    T_at_Q7Q5 = ifelse(is.na(T_at_Q7Q5), temp_at_q5_7d_C, T_at_Q7Q5),
    T_Q7Q5 = ifelse(is.na(T_Q7Q5), temp_during_q5_7d_C, T_Q7Q5)
  )
assert_unique_keys(thermal_lowflow, c("site", "year"), "thermal_lowflow")

# Merge annual tables into one site-year master table

HJA_annual <- rbi_recession %>%
  full_join(storage_fdc, by = c("site", "year")) %>%
  full_join(baseflow, by = c("site", "year")) %>%
  full_join(wb_storage, by = c("site", "year")) %>%
  full_join(thermal_lowflow, by = c("site", "year")) %>%
  # Keep only hydrometric sites used in analysis
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  # Keep site order consistent across outputs/plots
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site, year)
assert_unique_keys(HJA_annual, c("site", "year"), "HJA_annual")

# Save annual site-year master table
write.csv(HJA_annual,
          file.path(master_dir, MASTER_ANNUAL_FILE),
          row.names = FALSE)

# Calculate site-level means across available years

HJA_avg <- HJA_annual %>%
  group_by(site) %>%
  summarise(
    n_years = n(),
    across(
      where(is.numeric),
      list(mean = ~mean(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# Add isotope-based site-level metrics (not annual)

# MTT and Fyw
mtt_fyw <- read_csv(
  file.path(isotope_dir, "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  select(site, MTT = MTTM, Fyw = FYWM) %>%
  filter(!is.na(site), site != "", site %in% SITE_ORDER_HYDROMETRIC)

# Damping ratio
damping_ratios <- read_csv(
  file.path(isotope_dir, "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  select(site, DR = DR_Overall, DR_err = DR__err) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

# Combine isotope metrics into one site-level table
isotope_metrics <- mtt_fyw %>%
  full_join(damping_ratios, by = "site")
assert_unique_keys(isotope_metrics, c("site"), "isotope_metrics")

# Join isotope metrics into site-average table
HJA_avg <- HJA_avg %>%
  left_join(isotope_metrics, by = "site")

# Add static catchment characteristics

catchment_chars <- read_csv(
  resolve_catchment_characteristics_file(),
  show_col_types = FALSE
) %>%
  rename(site = Site) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)
assert_unique_keys(catchment_chars, c("site"), "catchment_chars")

HJA_avg <- HJA_avg %>%
  left_join(catchment_chars, by = "site")

# Save final site-level master table
write.csv(HJA_avg,
          file.path(master_dir, MASTER_SITE_FILE),
          row.names = FALSE)

# Track sample sizes by metric/site
# Storage and response metric names used below
dynamic_metrics <- DYNAMIC_METRICS
mobile_metrics_annual <- MOBILE_METRICS_ANNUAL
mobile_metrics_site <- MOBILE_METRICS_SITE
extended_metrics <- EXTENDED_DYNAMIC_METRICS
response_metrics <- c("T_7DMax", "Q_7Q5", "T_Q7Q5")

# Count available annual observations per site
annual_sample_sizes <- HJA_annual %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    n_years_total = n(),
    n_RBI = sum(!is.na(RBI)),
    n_RCS = sum(!is.na(RCS)),
    n_FDC = sum(!is.na(FDC)),
    n_SD = sum(!is.na(SD)),
    n_CHS = sum(!is.na(CHS)),
    n_WB = sum(!is.na(WB)),
    n_t_7dmax = sum(!is.na(T_7DMax)),
    n_q_7q5 = sum(!is.na(Q_7Q5)),
    .groups = "drop"
  )

# Add isotope availability flags
site_isotope <- HJA_avg %>%
  select(site, MTT, Fyw, DR) %>%
  mutate(
    has_MTT = !is.na(MTT),
    has_Fyw = !is.na(Fyw),
    has_DR = !is.na(DR)
  ) %>%
  select(site, has_MTT, has_Fyw, has_DR)

sample_sizes <- annual_sample_sizes %>%
  left_join(site_isotope, by = "site")

# Save sample-size summary table
write.csv(
  sample_sizes,
  file.path(master_dir, "metric_sample_sizes_by_site.csv"),
  row.names = FALSE
)

# Build per-site summary statistics for all metrics
annual_metric_cols <- intersect(
  c(
    dynamic_metrics,
    mobile_metrics_annual,
    extended_metrics,
    response_metrics
  ),
  names(HJA_annual)
)

annual_long <- HJA_annual %>%
  mutate(site = as.character(site)) %>%
  select(site, all_of(annual_metric_cols)) %>%
  pivot_longer(
    cols = all_of(annual_metric_cols),
    names_to = "metric",
    values_to = "value"
  )

site_level_long <- isotope_metrics %>%
  select(site, any_of(mobile_metrics_site)) %>%
  pivot_longer(
    cols = any_of(mobile_metrics_site),
    names_to = "metric",
    values_to = "value"
  )

site_metric_summary_stats <- bind_rows(annual_long, site_level_long) %>%
  group_by(site, metric) %>%
  summarise(
    n = sum(is.finite(value)),
    mean = ifelse(n > 0, mean(value, na.rm = TRUE), NA_real_),
    sd = ifelse(n > 1, stats::sd(value, na.rm = TRUE), NA_real_),
    min = ifelse(n > 0, min(value, na.rm = TRUE), NA_real_),
    median = ifelse(n > 0, stats::median(value, na.rm = TRUE), NA_real_),
    max = ifelse(n > 0, max(value, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    metric_group = case_when(
      metric %in% dynamic_metrics ~ "dynamic",
      metric %in% mobile_metrics_annual ~ "mobile_annual",
      metric %in% mobile_metrics_site ~ "mobile_site",
      metric %in% extended_metrics ~ "extended_dynamic",
      metric %in% response_metrics ~ "ecological_response",
      TRUE ~ "other"
    ),
    site = factor(site, levels = SITE_ORDER_HYDROMETRIC)
  ) %>%
  arrange(site, metric) %>%
  mutate(site = as.character(site))

write.csv(
  site_metric_summary_stats,
  file.path(master_dir, "master_site_metric_summary_stats.csv"),
  row.names = FALSE
)
