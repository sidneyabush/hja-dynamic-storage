# =============================================================================
# Dynamic Storage Drawdown (Water Balance Method)
# =============================================================================
# Purpose: Calculate dynamic storage drawdown during summer recession using
#          the water balance approach (P - Q - ET)
#
# Method:
#   1. Identify last hydrograph peak (>8% of annual max, before WY day 300)
#   2. From that date forward, compute daily dS = P - Q - ET
#   3. Calculate cumulative drawdown (ΔS_sum)
#   4. Extract minimum (maximum drawdown) for each site-year
#
# Timeline: Water Years 1998-2019
#
# Inputs:
#   - HF00402_v14.csv: Daily discharge
#   - daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv: P, Q, ET
#
# Outputs:
#   - DS_Drawdown_Date.csv: Timing of last peak for each site-year
#   - DS_drawdown_annual.csv: Annual maximum drawdown (mm) per site
#   - QA plots: Last peak visualization, drawdown timeseries
#
# Author: Keira Johnson (original), Sidney Bush (adapted)
# Date: 2026-01-23
# =============================================================================

# Load libraries
library(pracma)
library(dplyr)
library(ggplot2)
library(zoo)
library(tibble)
library(lubridate)

theme_set(theme_classic(base_size = 12))

# Clear environment
rm(list = ls())

# Helper functions to replace EflowStats
get_waterYear <- function(date) {
  ifelse(month(date) >= 10, year(date) + 1, year(date))
}

get_waterYearDay <- function(date) {
  wy <- get_waterYear(date)
  wy_start <- as.Date(paste0(wy - 1, "-10-01"))
  as.numeric(date - wy_start) + 1
}

# =============================================================================
# 1. SETUP: Directories
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir  <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"

discharge_dir <- file.path(base_dir, "Q")
storage_dir   <- file.path(base_dir, "DynamicStorage")

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 2. LOAD & PROCESS DISCHARGE DATA
# =============================================================================

discharge <- read.csv(file.path(discharge_dir, "HF00402_v14.csv")) %>%
  mutate(date = as.Date(DATE, "%m/%d/%Y"))

# Filter to complete water years only
goodyears <- discharge %>%
  filter(!(SITECODE %in% c("GSWSMC", "GSWSMF"))) %>%
  mutate(waterYear = get_waterYear(date)) %>%
  group_by(SITECODE, waterYear) %>%
  summarise(num_days = n_distinct(date), .groups = "drop") %>%
  filter(num_days >= 365)

discharge <- discharge %>%
  mutate(waterYear = get_waterYear(date)) %>%
  semi_join(goodyears, by = c("SITECODE", "waterYear"))

# Add 7-day smoothed discharge
discharge <- discharge %>%
  group_by(SITECODE) %>%
  arrange(date) %>%
  mutate(Q_smoothed = rollmean(MEAN_Q, k = 7, fill = NA, align = "right")) %>%
  ungroup()

# Remove rows with missing smoothed Q
discharge <- discharge[complete.cases(discharge$Q_smoothed), ]

# =============================================================================
# 3. DEFINE FUNCTION TO FIND LAST SIGNIFICANT PEAK
# =============================================================================

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
        wyd = get_waterYearDay(date)
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

# =============================================================================
# 4. IDENTIFY LAST PEAK FOR EACH SITE-YEAR
# =============================================================================

last_peak <- discharge %>%
  group_by(SITECODE, waterYear) %>%
  do(find_last_peak(.)) %>%
  ungroup() %>%
  mutate(wyd = get_waterYearDay(last_peak_date))

# Save last peak dates
write.csv(last_peak,
          file.path(output_dir, "DS_Drawdown_Date.csv"),
          row.names = FALSE)

# =============================================================================
# 5. LOAD WATER BALANCE DATA
# =============================================================================

DS_dat <- read.csv(file.path(storage_dir, "daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv")) %>%
  mutate(
    SITECODE = case_when(
      SITECODE == "GSLOOK_FULL" ~ "GSLOOK",
      SITECODE == "GSMACK" ~ "GSWSMA",
      .default = SITECODE
    ),
    DATE = as.Date(DATE, "%m/%d/%y"),
    waterYear = get_waterYear(DATE)
  ) %>%
  filter(waterYear > 1997 & waterYear < 2020)

# =============================================================================
# 6. MERGE WITH LAST PEAK DATES
# =============================================================================

DS_dat <- left_join(DS_dat, last_peak, by = c("SITECODE", "waterYear"))
DS_dat <- DS_dat[complete.cases(DS_dat$last_peak_date), ]

# =============================================================================
# 7. CALCULATE DYNAMIC STORAGE DRAWDOWN
# =============================================================================

# Crop to recession period (from last peak onward)
DS_dat_cropped <- DS_dat %>%
  group_by(SITECODE, waterYear) %>%
  filter(DATE >= last_peak_date) %>%
  mutate(max_q = Q_mm_d[last_peak_date == DATE])

# Compute daily storage change and cumulative drawdown
DS_compute <- DS_dat_cropped %>%
  group_by(SITECODE, waterYear) %>%
  mutate(
    DS_daily = P_mm_d - Q_mm_d - ET_mm_d,
    DS_sum = cumsum(DS_daily)
  ) %>%
  ungroup()

# =============================================================================
# 8. EXTRACT MAXIMUM DRAWDOWN (MINIMUM DS_SUM)
# =============================================================================

DS_max <- DS_compute %>%
  group_by(SITECODE, waterYear) %>%
  slice_min(DS_sum, n = 1) %>%
  select(SITECODE, waterYear, DS_sum) %>%
  ungroup()

# Save annual drawdown
write.csv(DS_max,
          file.path(output_dir, "DS_drawdown_annual.csv"),
          row.names = FALSE)

# =============================================================================
# 9. QA PLOTS
# =============================================================================

# Plot 1: Last peak timing distribution by site
p1 <- ggplot(last_peak, aes(x = wyd, y = SITECODE)) +
  geom_boxplot() +
  labs(
    title = "Timing of Last Significant Peak (>8% of annual max)",
    x = "Water Year Day",
    y = ""
  ) +
  theme_bw(base_size = 14)

ggsave(
  file.path(output_dir, "QA_LastPeak_DSDrawdown_Timing.png"),
  p1, width = 10, height = 6, dpi = 300
)

# Plot 2: Drawdown timeseries by site and year
p2 <- ggplot(DS_compute, aes(x = get_waterYearDay(DATE), y = DS_sum)) +
  geom_line(aes(color = waterYear, group = waterYear)) +
  facet_wrap(~SITECODE, ncol = 3, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient(low = "black", high = "grey88") +
  labs(
    title = "Dynamic Storage Drawdown During Summer Recession",
    x = "Water Year Day",
    y = "Cumulative Storage Change (mm)",
    color = "Water Year"
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(output_dir, "QA_DS_Drawdown_Timeseries.png"),
  p2, width = 14, height = 10, dpi = 300
)

# Plot 3: Annual maximum drawdown distribution
p3 <- ggplot(DS_max, aes(x = SITECODE, y = DS_sum)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
  labs(
    title = "Annual Maximum Dynamic Storage Drawdown by Site",
    x = "Site",
    y = "Maximum Drawdown (mm)"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(output_dir, "QA_DS_Drawdown_Distribution.png"),
  p3, width = 10, height = 6, dpi = 300
)
