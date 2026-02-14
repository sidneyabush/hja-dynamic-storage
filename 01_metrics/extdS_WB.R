# Calculate dynamic storage drawdown during summer recession using.
# Inputs: discharge_dir/HF00402_v14.csv.
# Author: Keira Johnson (original), Sidney Bush (adapted)
# Date: 2026-02-13

library(pracma)
library(dplyr)
library(ggplot2)
library(zoo)
library(tibble)
library(lubridate)

# Clear environment
rm(list = ls())

# Source configuration (paths, site definitions, water year range)
# Get script directory (works with source() and Rscript)
# Load project config
source("config.R")


theme_set(theme_pub(base_size = 12))

# SETUP: Directories (from config.R)

base_dir    <- BASE_DATA_DIR
output_dir  <- OUT_MET_EXTENDED_DIR

discharge_dir <- DISCHARGE_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# LOAD & PROCESS DISCHARGE DATA

discharge <- read.csv(file.path(discharge_dir, "HF00402_v14.csv")) %>%
  mutate(
    SITECODE = standardize_site_code(SITECODE),
    date = as.Date(DATE, "%m/%d/%Y"),
    waterYear = get_water_year(date)
  ) %>%
  # Filter to hydrometric sites and water year range
  filter(SITECODE %in% SITE_ORDER_HYDROMETRIC,
         waterYear >= WY_START, waterYear <= WY_END)

# Keep only water years with near-complete daily records.
goodyears <- discharge %>%
  group_by(SITECODE, waterYear) %>%
  summarise(num_days = n_distinct(date), .groups = "drop") %>%
  filter(num_days >= 365)

discharge <- discharge %>%
  semi_join(goodyears, by = c("SITECODE", "waterYear")) %>%
  # Set factor levels for consistent ordering
  mutate(SITECODE = factor(SITECODE, levels = SITE_ORDER_HYDROMETRIC))

# Smooth discharge a bit so peak detection is less noisy.
discharge <- discharge %>%
  group_by(SITECODE) %>%
  arrange(date) %>%
  mutate(Q_smoothed = rollmean(MEAN_Q, k = 7, fill = NA, align = "right")) %>%
  ungroup()

# Drop rows where smoothing could not be computed.
discharge <- discharge[complete.cases(discharge$Q_smoothed), ]

# DEFINE FUNCTION TO FIND LAST SIGNIFICANT PEAK

find_last_peak <- function(data, threshold_pct = 0.08) {
  time_series <- data$MEAN_Q
  max_peak_discharge <- max(time_series, na.rm = TRUE)
  threshold_value <- max_peak_discharge * threshold_pct

  peaks <- tryCatch({
    findpeaks(time_series)
  }, error = function(e) {
    return(NULL)
  })

  if (!is.null(peaks)) {
    peaks_df <- as_tibble(peaks) %>%
      rename(peak_height = V1, peak_index = V2) %>%
      filter(peak_height >= threshold_value) %>%
      mutate(
        date = data$date[peak_index],
        wyd = get_water_year_day(date)
      ) %>%
      filter(wyd < 300) %>%  # Before day 300 of water year
      arrange(peak_index)

    last_valid_peak <- peaks_df %>%
      slice_tail(n = 1)

    if (nrow(last_valid_peak) > 0) {
      return(tibble(
        last_peak_date = last_valid_peak$date,
        last_peak_value = time_series[last_valid_peak$peak_index]
      ))
    }
  }

  return(tibble(last_peak_date = NA, last_peak_value = NA))
}

# IDENTIFY LAST PEAK FOR EACH SITE-YEAR

last_peak <- discharge %>%
  group_by(SITECODE, waterYear) %>%
  do(find_last_peak(.)) %>%
  ungroup() %>%
  mutate(wyd = get_water_year_day(last_peak_date))

# Save last peak dates
write.csv(last_peak,
          file.path(output_dir, "ds_drawdown_date.csv"),
          row.names = FALSE)

# LOAD WATER BALANCE DATA

DS_dat <- read.csv(resolve_water_balance_daily_file()) %>%
  mutate(
    SITECODE = standardize_site_code(SITECODE),
    DATE = as.Date(DATE, "%m/%d/%y"),
    waterYear = get_water_year(DATE)
  ) %>%
  filter(waterYear >= WY_START, waterYear <= WY_END)

# MERGE WITH LAST PEAK DATES

DS_dat <- left_join(DS_dat, last_peak, by = c("SITECODE", "waterYear"))
DS_dat <- DS_dat[complete.cases(DS_dat$last_peak_date), ]

# CALCULATE DYNAMIC STORAGE DRAWDOWN

# Keep only dates after the seasonal peak for each site-year.
DS_dat_cropped <- DS_dat %>%
  group_by(SITECODE, waterYear) %>%
  filter(DATE >= last_peak_date) %>%
  mutate(max_q = Q_mm_d[last_peak_date == DATE])

# Calculate daily storage change and then sum through the recession.
DS_compute <- DS_dat_cropped %>%
  group_by(SITECODE, waterYear) %>%
  mutate(
    DS_daily = P_mm_d - Q_mm_d - ET_mm_d,
    WB = cumsum(DS_daily)
  ) %>%
  ungroup()

# EXTRACT MAXIMUM DRAWDOWN (MINIMUM DS_SUM)

DS_max <- DS_compute %>%
  group_by(SITECODE, waterYear) %>%
  slice_min(WB, n = 1) %>%
  select(SITECODE, waterYear, WB) %>%
  ungroup()

# Save annual drawdown
write.csv(DS_max,
          file.path(output_dir, "ds_drawdown_annual.csv"),
          row.names = FALSE)
