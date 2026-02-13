# -----------------------------------------------------------------------------
# Stream Temperature & Low-Flow Metrics for Storage Manuscript
# -----------------------------------------------------------------------------
# Purpose: Calculate ecologically-relevant thermal and low-flow metrics:
#   1. Maximum 7-day moving average stream temperature (per water year)
#   2. Q5 of 7-day moving average discharge (per water year)
#   3. Stream temperature during Q5 low-flow period
#   4. Q5_CV: Coefficient of variation of stream temp during low-flow period
#
# Timeline: Water Years 1997-2020
#
# Inputs:
#   - HT00451_v10.txt stream temperature data
#   - HF00402_v14.csv discharge data (daily)
#   - drainage_area.csv for normalization
#
# Outputs:
#   - stream_thermal_lowflow_metrics_annual.csv
#
# Author: Sidney Bush
# Date: 2026-01-23
# -----------------------------------------------------------------------------

# Load libraries
library(dplyr)
library(lubridate)
library(readr)
library(tidyr)
library(zoo)        # for rollmean()
library(ggplot2)
library(patchwork)

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

theme_set(theme_pub(base_size = 12))

# -----------------------------------------------------------------------------
# 1. SETUP: Directories and site list (from config.R)
# -----------------------------------------------------------------------------

base_dir      <- BASE_DATA_DIR
output_dir    <- OUT_MET_ECO_DIR
temp_dir      <- STREAM_TEMP_DIR
discharge_dir <- DISCHARGE_DIR
catchment_dir <- CATCHMENT_CHARACTERISTICS_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Target sites and years (from config.R)
sites_keep <- SITE_ORDER_HYDROMETRIC
target_years <- WY_START:WY_END

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

# -----------------------------------------------------------------------------
# 2. LOAD & PROCESS STREAM TEMPERATURE DATA
# -----------------------------------------------------------------------------

# Load daily stream temperature file
temp_file <- file.path(temp_dir, "HT00451_v10.txt")
if (!file.exists(temp_file)) {
  stop("Missing stream temperature file: ", temp_file)
}

temp_raw <- read_csv(temp_file, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(DATE_TIME),
    site = standardize_site_code(SITECODE)
  ) %>%
  select(site, date, temp_mean_C = WATERTEMP_MEAN) %>%
  filter(site %in% sites_keep)

# Aggregate to daily mean stream temperature per site-date
temp_daily <- temp_raw %>%
  group_by(site, date) %>%
  summarise(
    temp_mean_C = mean(temp_mean_C, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(temp_mean_C))

if (nrow(temp_daily) == 0 || n_distinct(temp_daily$site) < 3) {
  available_streams <- temp_raw %>%
    distinct(site) %>%
    pull(site) %>%
    unique()

  stop(
    paste0(
      "Stream temperature inputs do not cover hydrometric analysis sites. ",
      "Found stream IDs: ", paste(available_streams, collapse = ", "), ". ",
      "Expected site coverage includes: ", paste(sites_keep, collapse = ", "), "."
    )
  )
}

# WS09 has no stream-temperature records in HT00451 source.
# Keep WS09 in downstream annual masters; temp responses remain NA for WS09.

# -----------------------------------------------------------------------------
# 3. LOAD & PROCESS DISCHARGE DATA
# -----------------------------------------------------------------------------

# Load drainage areas
da_df <- read_csv(resolve_drainage_area_file(),
                  show_col_types = FALSE)

# Load discharge
discharge <- read_csv(file.path(discharge_dir, "HF00402_v14.csv"),
                      show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(SITECODE)) %>%
  left_join(da_df, by = "SITECODE") %>%
  filter(!is.na(DA_M2), site %in% sites_keep) %>%
  mutate(
    date = as.Date(DATE, format = "%m/%d/%Y"),
    Q_cms = MEAN_Q * 0.02831683199881,  # Convert to m³/s
    Q_mm_d = (Q_cms * 86400) / DA_M2 * 1000  # Convert to mm/day
  ) %>%
  select(site, date, Q_cms, Q_mm_d, WATERYEAR) %>%
  arrange(site, date)

# -----------------------------------------------------------------------------
# 4. CALCULATE 7-DAY MOVING AVERAGES
# -----------------------------------------------------------------------------

# Function to calculate 7-day rolling mean
calc_7day_rolling <- function(df, value_col) {
  df %>%
    arrange(date) %>%
    mutate(
      rolling_7d = zoo::rollmean(.data[[value_col]], k = 7,
                                  fill = NA, align = "center")
    )
}

# 7-day rolling average stream temperature (by stream)
temp_rolling <- temp_daily %>%
  group_by(site) %>%
  calc_7day_rolling("temp_mean_C") %>%
  ungroup() %>%
  mutate(year = get_water_year(date)) %>%
  rename(temp_7d_avg_C = rolling_7d)

# 7-day rolling average discharge (by site)
discharge_rolling <- discharge %>%
  group_by(site) %>%
  calc_7day_rolling("Q_mm_d") %>%
  ungroup() %>%
  mutate(year = get_water_year(date)) %>%
  rename(Q_7d_avg_mm_d = rolling_7d)

# -----------------------------------------------------------------------------
# 5. EXTRACT ANNUAL METRICS (PER WATER YEAR)
# -----------------------------------------------------------------------------

# 5.1 Maximum 7-day average temperature per water year
t_7dmax <- temp_rolling %>%
  filter(year %in% target_years) %>%
  group_by(site, year) %>%
  filter(!is.na(temp_7d_avg_C)) %>%
  slice_max(temp_7d_avg_C, n = 1, with_ties = FALSE) %>%
  select(
    site,
    year,
    date_t_7dmax = date,
    T_7DMax = temp_7d_avg_C
  ) %>%
  ungroup()
assert_unique_keys(t_7dmax, c("site", "year"), "t_7dmax")

# 5.2 Q5 of 7-day average discharge per water year
q_7q5 <- discharge_rolling %>%
  filter(year %in% target_years) %>%
  group_by(site, year) %>%
  summarise(
    Q_7Q5 = quantile(Q_7d_avg_mm_d, probs = 0.05, na.rm = TRUE),
    .groups = "drop"
  )
assert_unique_keys(q_7q5, c("site", "year"), "q_7q5")

# Representative date at Q5 threshold (closest 7-day Q to annual Q5)
q_7q5_date <- discharge_rolling %>%
  filter(year %in% target_years, !is.na(Q_7d_avg_mm_d)) %>%
  inner_join(q_7q5, by = c("site", "year")) %>%
  group_by(site, year) %>%
  slice_min(abs(Q_7d_avg_mm_d - Q_7Q5), n = 1, with_ties = FALSE) %>%
  select(
    site,
    year,
    date_q_7q5 = date
  ) %>%
  ungroup()
assert_unique_keys(q_7q5_date, c("site", "year"), "q_7q5_date")

# 5.3 Temperature at and during Q5 low-flow period
# - T_at_Q7Q5: temperature on representative Q5 date
# - T_Q7Q5: mean temperature across all days with 7-day Q <= annual Q5
temp_at_q5 <- q_7q5_date %>%
  left_join(
    temp_daily %>%
      select(site, date, temp_mean_C) %>%
      rename(T_at_Q7Q5 = temp_mean_C),
    by = c("site", "date_q_7q5" = "date")
  )
assert_unique_keys(temp_at_q5, c("site", "year"), "temp_at_q5")

q5_period_temp <- discharge_rolling %>%
  filter(year %in% target_years, !is.na(Q_7d_avg_mm_d)) %>%
  inner_join(q_7q5, by = c("site", "year")) %>%
  filter(Q_7d_avg_mm_d <= Q_7Q5) %>%
  select(site, year, date) %>%
  inner_join(temp_daily %>% select(site, date, temp_mean_C), by = c("site", "date")) %>%
  group_by(site, year) %>%
  summarise(
    T_Q7Q5 = mean(temp_mean_C, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    T_Q7Q5 = if_else(
      is.nan(T_Q7Q5),
      as.numeric(NA),
      T_Q7Q5
    )
  )
assert_unique_keys(q5_period_temp, c("site", "year"), "q5_period_temp")

# 5.4 Q5_CV: Coefficient of variation of stream temp during low-flow period
# CV of daily mean stream temperature during Aug-Oct (late-summer recession)
# Lower values indicate greater thermal buffering by subsurface storage
temp_cv_lowflow <- temp_daily %>%
  mutate(
    year = get_water_year(date),
    month_num = month(date)
  ) %>%
  filter(year %in% target_years, month_num %in% 8:10) %>%  # Aug-Oct
  group_by(site, year) %>%
  summarise(
    temp_mean = mean(temp_mean_C, na.rm = TRUE),
    temp_sd = sd(temp_mean_C, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  ) %>%
  mutate(Q5_CV = temp_sd / temp_mean) %>%
  select(site, year, Q5_CV)
assert_unique_keys(temp_cv_lowflow, c("site", "year"), "temp_cv_lowflow")

# -----------------------------------------------------------------------------
# 6. COMBINE METRICS INTO MASTER TABLE
# -----------------------------------------------------------------------------

# Merge all metrics by standardized site-year keys
master_metrics <- t_7dmax %>%
  full_join(q_7q5, by = c("site", "year")) %>%
  left_join(q_7q5_date, by = c("site", "year")) %>%
  left_join(temp_at_q5, by = c("site", "year", "date_q_7q5")) %>%
  left_join(q5_period_temp, by = c("site", "year")) %>%
  left_join(temp_cv_lowflow, by = c("site", "year")) %>%
  # Keep compatibility aliases while names transition.
  mutate(
    max_temp_7d_C = T_7DMax,
    q5_7d_mm_d = Q_7Q5,
    temp_during_q5_7d_C = T_Q7Q5,
    min_Q_7d_mm_d = Q_7Q5,
    temp_at_q5_7d_C = T_at_Q7Q5,
    temp_at_min_Q_7d_C = T_at_Q7Q5,
    temp_during_min_Q_7d_C = T_Q7Q5
  ) %>%
  mutate(
    SITECODE = site,
    wateryear = year
  ) %>%
  arrange(site, year)
assert_unique_keys(master_metrics, c("site", "year"), "master_metrics")

# Save output
output_file <- file.path(output_dir, "stream_thermal_lowflow_metrics_annual.csv")
write_csv(master_metrics, output_file)

# -----------------------------------------------------------------------------
# 7. SUMMARY STATISTICS
# -----------------------------------------------------------------------------

summary_stats <- master_metrics %>%
  group_by(site) %>%
  summarise(
    n_years = n(),
    t_7dmax_mean = mean(T_7DMax, na.rm = TRUE),
    t_7dmax_sd   = sd(T_7DMax, na.rm = TRUE),
    q_7q5_mean   = mean(Q_7Q5, na.rm = TRUE),
    q_7q5_sd     = sd(Q_7Q5, na.rm = TRUE),
    t_q7q5_mean  = mean(T_Q7Q5, na.rm = TRUE),
    t_q7q5_sd    = sd(T_Q7Q5, na.rm = TRUE),
    Q5_CV_mean = mean(Q5_CV, na.rm = TRUE),
    Q5_CV_sd   = sd(Q5_CV, na.rm = TRUE),
    .groups = "drop"
  )
