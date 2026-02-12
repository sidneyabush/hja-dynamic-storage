# -----------------------------------------------------------------------------
# Stream Temperature & Low-Flow Metrics for Storage Manuscript
# -----------------------------------------------------------------------------
# Purpose: Calculate ecologically-relevant thermal and low-flow metrics:
#   1. Maximum 7-day moving average stream temperature (per calendar year)
#   2. Minimum 7-day moving average discharge (per calendar year)
#   3. Stream temperature at time of minimum 7-day discharge
#   4. Q5_CV: Coefficient of variation of stream temp during low-flow period
#
# Timeline: Calendar Years 1997-2020
#
# Inputs:
#   - HT00401_v8.csv daily stream temperature data
#   - HF00402_v14.csv discharge data (daily)
#   - drainage_area.csv for normalization
#
# Outputs:
#   - stream_thermal_lowflow_metrics_annual.csv
#   - QA plots showing time series and distributions
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

# -----------------------------------------------------------------------------
# 1. SETUP: Directories and site list (from config.R)
# -----------------------------------------------------------------------------

base_dir      <- BASE_DATA_DIR
output_dir    <- OUTPUT_DIR
temp_dir      <- STREAM_TEMP_DIR
discharge_dir <- DISCHARGE_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Target sites and years (from config.R)
sites_keep <- SITE_ORDER_HYDROMETRIC
target_years <- WY_START:WY_END

# -----------------------------------------------------------------------------
# 2. LOAD & PROCESS STREAM TEMPERATURE DATA
# -----------------------------------------------------------------------------

# Load daily stream temperature file
temp_file <- file.path(temp_dir, "HT00401_v8.csv")
if (!file.exists(temp_file)) {
  stop("Missing stream temperature file: ", temp_file)
}

temp_raw <- read_csv(temp_file, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(DATE),
    site = standardize_site_code(SITECODE)
  ) %>%
  select(site, date, temp_mean_C = AIRTEMP_MEAN_DAY) %>%
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

# -----------------------------------------------------------------------------
# 3. LOAD & PROCESS DISCHARGE DATA
# -----------------------------------------------------------------------------

# Load drainage areas
da_df <- read_csv(file.path(discharge_dir, "drainage_area.csv"),
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
  mutate(year = year(date)) %>%
  rename(temp_7d_avg_C = rolling_7d)

# 7-day rolling average discharge (by site)
discharge_rolling <- discharge %>%
  group_by(site) %>%
  calc_7day_rolling("Q_mm_d") %>%
  ungroup() %>%
  mutate(year = year(date)) %>%
  rename(Q_7d_avg_mm_d = rolling_7d)

# -----------------------------------------------------------------------------
# 5. EXTRACT ANNUAL METRICS (PER CALENDAR YEAR)
# -----------------------------------------------------------------------------

# 5.1 Maximum 7-day average temperature per calendar year
max_temp_7d <- temp_rolling %>%
  filter(year %in% target_years) %>%
  group_by(site, year) %>%
  filter(!is.na(temp_7d_avg_C)) %>%
  slice_max(temp_7d_avg_C, n = 1, with_ties = FALSE) %>%
  select(
    site,
    year,
    date_max_temp_7d = date,
    max_temp_7d_C = temp_7d_avg_C
  ) %>%
  ungroup()

# 5.2 Minimum 7-day average discharge per calendar year
min_Q_7d <- discharge_rolling %>%
  filter(year %in% target_years) %>%
  group_by(site, year) %>%
  filter(!is.na(Q_7d_avg_mm_d)) %>%
  slice_min(Q_7d_avg_mm_d, n = 1, with_ties = FALSE) %>%
  select(
    site,
    year,
    date_min_Q_7d = date,
    min_Q_7d_mm_d = Q_7d_avg_mm_d
  ) %>%
  ungroup()

# 5.3 Temperature during minimum 7-day discharge window
# - temp_at_min_Q_7d_C: center-day temperature at minimum 7-day Q
# - temp_during_min_Q_7d_C: mean temperature across the same 7-day window
temp_at_min_Q <- min_Q_7d %>%
  left_join(
    temp_daily %>%
      select(site, date, temp_mean_C) %>%
      rename(temp_at_min_Q_7d_C = temp_mean_C),
    by = c("site", "date_min_Q_7d" = "date")
  ) %>%
  rowwise() %>%
  mutate(
    temp_during_min_Q_7d_C = mean(
      temp_daily$temp_mean_C[
        temp_daily$site == site &
          temp_daily$date >= (date_min_Q_7d - 3) &
          temp_daily$date <= (date_min_Q_7d + 3)
      ],
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  mutate(
    temp_during_min_Q_7d_C = if_else(
      is.nan(temp_during_min_Q_7d_C),
      as.numeric(NA),
      temp_during_min_Q_7d_C
    )
  ) %>%
  select(site, year, date_min_Q_7d, temp_at_min_Q_7d_C, temp_during_min_Q_7d_C)

# 5.4 Q5_CV: Coefficient of variation of stream temp during low-flow period
# CV of daily mean stream temperature during Aug-Oct (late-summer recession)
# Lower values indicate greater thermal buffering by subsurface storage
temp_cv_lowflow <- temp_daily %>%
  mutate(
    year = year(date),
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

# -----------------------------------------------------------------------------
# 6. COMBINE METRICS INTO MASTER TABLE
# -----------------------------------------------------------------------------

# Merge all metrics by standardized site-year keys
master_metrics <- max_temp_7d %>%
  full_join(min_Q_7d, by = c("site", "year")) %>%
  left_join(temp_at_min_Q, by = c("site", "year", "date_min_Q_7d")) %>%
  left_join(temp_cv_lowflow, by = c("site", "year")) %>%
  mutate(
    SITECODE = site,
    wateryear = year
  ) %>%
  arrange(site, year)

# Save output
output_file <- file.path(output_dir, "stream_thermal_lowflow_metrics_annual.csv")
write_csv(master_metrics, output_file)

# -----------------------------------------------------------------------------
# 7. QA PLOTS
# -----------------------------------------------------------------------------

# Plot 1: Time series of 7-day rolling metrics
p1_temp_data <- temp_rolling %>%
  filter(year %in% target_years)

if (nrow(p1_temp_data) > 0) {
  p1_temp <- p1_temp_data %>%
    ggplot(aes(x = date, y = temp_7d_avg_C, color = site)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~site, scales = "free_y", ncol = 2) +
    labs(
      title = "7-Day Rolling Average Stream Temperature",
      subtitle = paste("Calendar Years", min(target_years), "-", max(target_years)),
      x = "Date",
      y = "Temperature (°C)",
      color = "Site"
    ) +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "QA_7day_temp_timeseries.png"),
         p1_temp, width = 12, height = 10, dpi = 300)
}

p1_Q_data <- discharge_rolling %>%
  filter(year %in% target_years, site %in% sites_keep)

if (nrow(p1_Q_data) > 0) {
  p1_Q <- p1_Q_data %>%
    ggplot(aes(x = date, y = Q_7d_avg_mm_d, color = site)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~site, scales = "free_y", ncol = 2) +
    labs(
      title = "7-Day Rolling Average Discharge",
      subtitle = paste("Calendar Years", min(target_years), "-", max(target_years)),
      x = "Date",
      y = "Discharge (mm/d)",
      color = "Site"
    ) +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "QA_7day_Q_timeseries.png"),
         p1_Q, width = 12, height = 10, dpi = 300)
}

# Plot 2: Annual metrics distributions
p2_data <- master_metrics %>%
  pivot_longer(
    cols = c(max_temp_7d_C, min_Q_7d_mm_d, temp_during_min_Q_7d_C, Q5_CV),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

if (nrow(p2_data) > 0) {
  p2 <- p2_data %>%
    ggplot(aes(x = year, y = value, color = site)) +
    geom_line() +
    geom_point() +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    labs(
      title = "Annual Thermal & Low-Flow Metrics",
      x = "Calendar Year",
      y = "Value",
      color = "Site"
    ) +
    theme_minimal()

  ggsave(file.path(output_dir, "QA_annual_metrics_timeseries.png"),
         p2, width = 12, height = 10, dpi = 300)
}

# -----------------------------------------------------------------------------
# 8. SUMMARY STATISTICS
# -----------------------------------------------------------------------------

summary_stats <- master_metrics %>%
  group_by(site) %>%
  summarise(
    n_years = n(),
    max_temp_7d_mean = mean(max_temp_7d_C, na.rm = TRUE),
    max_temp_7d_sd   = sd(max_temp_7d_C, na.rm = TRUE),
    min_Q_7d_mean    = mean(min_Q_7d_mm_d, na.rm = TRUE),
    min_Q_7d_sd      = sd(min_Q_7d_mm_d, na.rm = TRUE),
    temp_at_min_Q_mean = mean(temp_during_min_Q_7d_C, na.rm = TRUE),
    temp_at_min_Q_sd   = sd(temp_during_min_Q_7d_C, na.rm = TRUE),
    Q5_CV_mean = mean(Q5_CV, na.rm = TRUE),
    Q5_CV_sd   = sd(Q5_CV, na.rm = TRUE),
    .groups = "drop"
  )
