# =============================================================================
# Stream Temperature & Low-Flow Metrics for Storage Manuscript
# =============================================================================
# Purpose: Calculate ecologically-relevant thermal and low-flow metrics:
#   1. Maximum 7-day moving average stream temperature (per water year)
#   2. Minimum 7-day moving average discharge (per water year)
#   3. Stream temperature at time of minimum 7-day discharge
#   4. Q5_CV: Coefficient of variation of stream temp during low-flow period
#
# Timeline: Water Years 1997-2020
#
# Inputs:
#   - HT002* stream temperature files (30-min data)
#   - HF00402_v14.csv discharge data (daily)
#   - drainage_area.csv for normalization
#
# Outputs:
#   - stream_thermal_lowflow_metrics_annual.csv
#   - QA plots showing time series and distributions
#
# Author: Sidney Bush
# Date: 2026-01-23
# =============================================================================

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

config_path <- file.path(dirname(script_dir), "config.R")
if (!file.exists(config_path)) {
  config_path <- file.path(getwd(), "config.R")
}
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found. Please ensure config.R exists in the repo root.")
}

# =============================================================================
# 1. SETUP: Directories and site list (from config.R)
# =============================================================================

base_dir      <- BASE_DATA_DIR
output_dir    <- OUTPUT_DIR
temp_dir      <- STREAM_TEMP_DIR
discharge_dir <- DISCHARGE_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Target sites and years (from config.R)
sites_keep   <- SITE_ORDER_HYDROMETRIC
target_years <- WY_START:WY_END

# =============================================================================
# 2. LOAD & PROCESS STREAM TEMPERATURE DATA
# =============================================================================

# Find all HT002 stream temperature files
temp_files <- list.files(temp_dir, pattern = "^HT002.*\\.(csv|txt)$",
                         full.names = TRUE, ignore.case = TRUE)

if (length(temp_files) == 0) {
  stop("No HT002 stream temperature files found in: ", temp_dir)
}

# Read and combine all temperature files
temp_raw <- lapply(temp_files, function(f) {
  tryCatch({
    read_csv(f, show_col_types = FALSE) %>%
      select(DATE_TIME, STREAM, TEMP_C, MEDIA) %>%
      filter(MEDIA == "Water") %>%  # Stream water only
      mutate(
        datetime = as.POSIXct(DATE_TIME, tz = "America/Los_Angeles"),
        date     = as.Date(datetime)
      )
  }, error = function(e) {
    NULL
  })
}) %>%
  bind_rows()

# Aggregate to daily mean stream temperature
temp_daily <- temp_raw %>%
  group_by(STREAM, date) %>%
  summarise(
    temp_mean_C = mean(TEMP_C, na.rm = TRUE),
    n_obs       = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(temp_mean_C))

# =============================================================================
# 3. LOAD & PROCESS DISCHARGE DATA
# =============================================================================

# Load drainage areas
da_df <- read_csv(file.path(discharge_dir, "drainage_area.csv"),
                  show_col_types = FALSE)

# Load discharge
discharge <- read_csv(file.path(discharge_dir, "HF00402_v14.csv"),
                      show_col_types = FALSE) %>%
  left_join(da_df, by = "SITECODE") %>%
  filter(!is.na(DA_M2), SITECODE %in% sites_keep) %>%
  mutate(
    date = as.Date(DATE, format = "%m/%d/%Y"),
    Q_cms = MEAN_Q * 0.02831683199881,  # Convert to m³/s
    Q_mm_d = (Q_cms * 86400) / DA_M2 * 1000  # Convert to mm/day
  ) %>%
  select(SITECODE, date, Q_cms, Q_mm_d, WATERYEAR) %>%
  arrange(SITECODE, date)

# =============================================================================
# 4. CALCULATE 7-DAY MOVING AVERAGES
# =============================================================================

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
  group_by(STREAM) %>%
  calc_7day_rolling("temp_mean_C") %>%
  ungroup() %>%
  mutate(
    wateryear = if_else(month(date) >= 10,
                       year(date) + 1,
                       year(date))
  ) %>%
  rename(temp_7d_avg_C = rolling_7d)

# 7-day rolling average discharge (by site)
discharge_rolling <- discharge %>%
  group_by(SITECODE) %>%
  calc_7day_rolling("Q_mm_d") %>%
  ungroup() %>%
  rename(Q_7d_avg_mm_d = rolling_7d)

# =============================================================================
# 5. EXTRACT ANNUAL METRICS (PER WATER YEAR)
# =============================================================================

# 5.1 Maximum 7-day average temperature per water year
max_temp_7d <- temp_rolling %>%
  filter(wateryear %in% target_years) %>%
  group_by(STREAM, wateryear) %>%
  filter(!is.na(temp_7d_avg_C)) %>%
  slice_max(temp_7d_avg_C, n = 1, with_ties = FALSE) %>%
  select(
    STREAM,
    wateryear,
    date_max_temp_7d = date,
    max_temp_7d_C = temp_7d_avg_C
  ) %>%
  ungroup()

# 5.2 Minimum 7-day average discharge per water year
min_Q_7d <- discharge_rolling %>%
  filter(WATERYEAR %in% target_years) %>%
  group_by(SITECODE, WATERYEAR) %>%
  filter(!is.na(Q_7d_avg_mm_d)) %>%
  slice_min(Q_7d_avg_mm_d, n = 1, with_ties = FALSE) %>%
  select(
    SITECODE,
    wateryear = WATERYEAR,
    date_min_Q_7d = date,
    min_Q_7d_mm_d = Q_7d_avg_mm_d
  ) %>%
  ungroup()

# 5.3 Temperature at time of minimum 7-day discharge
# Join temperature data to min discharge timing
temp_at_min_Q <- min_Q_7d %>%
  left_join(
    temp_daily %>% select(date, STREAM, temp_mean_C),
    by = c("date_min_Q_7d" = "date")
  ) %>%
  rename(temp_at_min_Q_7d_C = temp_mean_C) %>%
  select(SITECODE, STREAM, wateryear, date_min_Q_7d, temp_at_min_Q_7d_C)

# 5.4 Q5_CV: Coefficient of variation of stream temp during low-flow period
# CV of daily mean stream temperature during Aug-Oct (late-summer recession)
# Lower values indicate greater thermal buffering by subsurface storage
temp_cv_lowflow <- temp_daily %>%
  mutate(
    wateryear = if_else(month(date) >= 10,
                       year(date) + 1,
                       year(date)),
    month_num = month(date)
  ) %>%
  filter(wateryear %in% target_years, month_num %in% 8:10) %>%  # Aug-Oct
  group_by(STREAM, wateryear) %>%
  summarise(
    temp_mean = mean(temp_mean_C, na.rm = TRUE),
    temp_sd = sd(temp_mean_C, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  ) %>%
  mutate(Q5_CV = temp_sd / temp_mean) %>%
  select(STREAM, wateryear, Q5_CV)

# =============================================================================
# 6. COMBINE METRICS INTO MASTER TABLE
# =============================================================================

# Merge all metrics
# Note: Stream names and SITECODE may not match perfectly, handle carefully
master_metrics <- max_temp_7d %>%
  full_join(
    min_Q_7d,
    by = c("wateryear")
  ) %>%
  left_join(
    temp_at_min_Q %>% select(SITECODE, wateryear, temp_at_min_Q_7d_C),
    by = c("SITECODE", "wateryear")
  ) %>%
  left_join(
    temp_cv_lowflow,
    by = c("STREAM", "wateryear")
  ) %>%
  arrange(SITECODE, wateryear)

# Save output
output_file <- file.path(output_dir, "stream_thermal_lowflow_metrics_annual.csv")
write_csv(master_metrics, output_file)

# =============================================================================
# 7. QA PLOTS
# =============================================================================

# Plot 1: Time series of 7-day rolling metrics
p1_temp <- temp_rolling %>%
  filter(wateryear %in% target_years) %>%
  ggplot(aes(x = date, y = temp_7d_avg_C, color = STREAM)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~STREAM, scales = "free_y", ncol = 2) +
  labs(
    title = "7-Day Rolling Average Stream Temperature",
    subtitle = paste("Water Years", min(target_years), "-", max(target_years)),
    x = "Date",
    y = "Temperature (°C)",
    color = "Stream"
  ) +
  theme(legend.position = "none")

p1_Q <- discharge_rolling %>%
  filter(WATERYEAR %in% target_years, SITECODE %in% sites_keep) %>%
  ggplot(aes(x = date, y = Q_7d_avg_mm_d, color = SITECODE)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~SITECODE, scales = "free_y", ncol = 2) +
  labs(
    title = "7-Day Rolling Average Discharge",
    subtitle = paste("Water Years", min(target_years), "-", max(target_years)),
    x = "Date",
    y = "Discharge (mm/d)",
    color = "Site"
  ) +
  theme(legend.position = "none")

# Save plots
ggsave(file.path(output_dir, "QA_7day_temp_timeseries.png"),
       p1_temp, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "QA_7day_Q_timeseries.png"),
       p1_Q, width = 12, height = 10, dpi = 300)

# Plot 2: Annual metrics distributions
p2 <- master_metrics %>%
  pivot_longer(
    cols = c(max_temp_7d_C, min_Q_7d_mm_d, temp_at_min_Q_7d_C, Q5_CV),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = wateryear, y = value, color = SITECODE)) +
  geom_line() +
  geom_point() +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  labs(
    title = "Annual Thermal & Low-Flow Metrics",
    x = "Water Year",
    y = "Value",
    color = "Site"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "QA_annual_metrics_timeseries.png"),
       p2, width = 12, height = 10, dpi = 300)

# =============================================================================
# 8. SUMMARY STATISTICS
# =============================================================================

summary_stats <- master_metrics %>%
  group_by(SITECODE) %>%
  summarise(
    n_years = n(),
    max_temp_7d_mean = mean(max_temp_7d_C, na.rm = TRUE),
    max_temp_7d_sd   = sd(max_temp_7d_C, na.rm = TRUE),
    min_Q_7d_mean    = mean(min_Q_7d_mm_d, na.rm = TRUE),
    min_Q_7d_sd      = sd(min_Q_7d_mm_d, na.rm = TRUE),
    temp_at_min_Q_mean = mean(temp_at_min_Q_7d_C, na.rm = TRUE),
    temp_at_min_Q_sd   = sd(temp_at_min_Q_7d_C, na.rm = TRUE),
    Q5_CV_mean = mean(Q5_CV, na.rm = TRUE),
    Q5_CV_sd   = sd(Q5_CV, na.rm = TRUE),
    .groups = "drop"
  )
