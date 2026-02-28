# calculate ecologically-relevant thermal, low-flow, and seasonal-precip metrics.
# inputs: discharge_dir/hf00402_v14.csv; out_met_support_dir/catchments_met_q.csv.
# author: sidney bush
# date: 2026-01-23

library(dplyr)
library(lubridate)
library(readr)
library(tidyr)
library(zoo)        # for rollmean()
library(ggplot2)
library(patchwork)

# clear environment
rm(list = ls())

# source configuration (paths, site definitions, water year range)
# get script directory (works with source() and rscript)
# load project config
source("config.R")


theme_set(theme_pub(base_size = 12))

# setup: directories and site list (from config.r)

base_dir      <- BASE_DATA_DIR
output_dir    <- OUT_MET_ECO_DIR
temp_dir      <- STREAM_TEMP_DIR

# create output directory if needed
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# target sites and years (from config.r)
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

# load & process stream temperature data

# load daily stream temperature file
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

# aggregate to daily mean stream temperature per site-date
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

# ws09 has no stream-temperature records in ht00451 source.
# keep ws09 in downstream annual masters; temp responses remain na for ws09.

# load & process preprocessed daily met data

met_support_file <- file.path(OUT_MET_SUPPORT_DIR, "catchments_met_q.csv")
if (!file.exists(met_support_file)) {
  stop(
    "Missing paired met support file: ", met_support_file,
    ". Run 00_data_preprocessing/create_hydromet_master.R first."
  )
}

met_support <- read_csv(met_support_file, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%m/%d/%y")),
    site = standardize_site_code(SITECODE),
    P_mm_d = suppressWarnings(as.numeric(P_mm_d)),
    Q_mm_d = suppressWarnings(as.numeric(Q_mm_d)),
    WATERYEAR = get_water_year(date)
  ) %>%
  filter(site %in% sites_keep, !is.na(date))
assert_unique_keys(met_support, c("site", "date"), "met_support")

# use standardized q_mm_d from creation-step support table.
# no downstream unit conversion should occur in this script.
discharge <- met_support %>%
  select(site, date, Q_mm_d, WATERYEAR) %>%
  arrange(site, date)

met_daily <- met_support %>%
  select(site, date, P_mm_d)

# wet-season precipitation (nov-may) by season year.
# example: nov-dec 2019 + jan-may 2020 are assigned to year 2020.
precip_nov_may <- met_daily %>%
  mutate(
    month_num = month(date),
    season_year = if_else(month_num %in% c(11L, 12L), year(date) + 1L, year(date))
  ) %>%
  filter(month_num %in% c(11L, 12L, 1L, 2L, 3L, 4L, 5L), season_year %in% target_years) %>%
  group_by(site, year = season_year) %>%
  summarise(
    n_days = sum(is.finite(P_mm_d)),
    Pws = if_else(n_days > 0, sum(P_mm_d, na.rm = TRUE), as.numeric(NA)),
    .groups = "drop"
  ) %>%
  mutate(
    precip_nov_may_mm = Pws
  ) %>%
  select(site, year, Pws, precip_nov_may_mm)
assert_unique_keys(precip_nov_may, c("site", "year"), "precip_nov_may")

# calculate 7-day moving averages

# function to calculate 7-day rolling mean
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

# extract annual metrics (per water year)

# maximum 7-day average temperature per water year
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

# q5 of 7-day average discharge per water year
q_7q5 <- discharge_rolling %>%
  filter(year %in% target_years) %>%
  group_by(site, year) %>%
  summarise(
    Q_7Q5 = quantile(Q_7d_avg_mm_d, probs = 0.05, na.rm = TRUE),
    .groups = "drop"
  )
assert_unique_keys(q_7q5, c("site", "year"), "q_7q5")

# representative date at q5 threshold (closest 7-day q to annual q5)
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

# temperature at q5 low-flow date
# - t_at_q7q5: temperature on representative q5 date
temp_at_q5 <- q_7q5_date %>%
  left_join(
    temp_daily %>%
      select(site, date, temp_mean_C) %>%
      rename(T_at_Q7Q5 = temp_mean_C),
    by = c("site", "date_q_7q5" = "date")
  )
assert_unique_keys(temp_at_q5, c("site", "year"), "temp_at_q5")

# q5_cv: coefficient of variation of stream temp during low-flow period
# cv of daily mean stream temperature during aug-oct (late-summer recession)
# lower values indicate greater thermal buffering by subsurface storage
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

# combine metrics into master table

# merge all metrics by standardized site-year keys
master_metrics <- t_7dmax %>%
  full_join(q_7q5, by = c("site", "year")) %>%
  left_join(q_7q5_date, by = c("site", "year")) %>%
  left_join(temp_at_q5, by = c("site", "year", "date_q_7q5")) %>%
  left_join(precip_nov_may, by = c("site", "year")) %>%
  left_join(temp_cv_lowflow, by = c("site", "year")) %>%
  # keep compatibility aliases while names transition.
  mutate(
    max_temp_7d_C = T_7DMax,
    q5_7d_mm_d = Q_7Q5,
    min_Q_7d_mm_d = Q_7Q5,
    temp_at_q5_7d_C = T_at_Q7Q5,
    temp_at_min_Q_7d_C = T_at_Q7Q5,
    precip_nov_may_mm = ifelse(is.na(precip_nov_may_mm), Pws, precip_nov_may_mm)
  ) %>%
  mutate(
    SITECODE = site,
    wateryear = year
  ) %>%
  arrange(site, year)
assert_unique_keys(master_metrics, c("site", "year"), "master_metrics")

# save output
output_file <- file.path(output_dir, "stream_thermal_lowflow_metrics_annual.csv")
write_csv(master_metrics, output_file)

# summary statistics

summary_stats <- master_metrics %>%
  group_by(site) %>%
  summarise(
    n_years = n(),
    t_7dmax_mean = mean(T_7DMax, na.rm = TRUE),
    t_7dmax_sd   = sd(T_7DMax, na.rm = TRUE),
    q_7q5_mean   = mean(Q_7Q5, na.rm = TRUE),
    q_7q5_sd     = sd(Q_7Q5, na.rm = TRUE),
    p_wetseason_mean = mean(Pws, na.rm = TRUE),
    p_wetseason_sd   = sd(Pws, na.rm = TRUE),
    Q5_CV_mean = mean(Q5_CV, na.rm = TRUE),
    Q5_CV_sd   = sd(Q5_CV, na.rm = TRUE),
    .groups = "drop"
  )
