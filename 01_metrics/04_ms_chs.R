# -----------------------------------------------------------------------------
# Chemical Hydrograph Separation - Baseflow Proportion Calculation
# -----------------------------------------------------------------------------
# Purpose: Use specific conductance (EC) to estimate annual mean baseflow
#          proportion via two-component hydrograph separation
#
# Method: Chemical hydrograph separation using specific conductance as tracer
#   - Groundwater endmember = 99th percentile SC (high SC = deep groundwater)
#   - Runoff endmember = 1st percentile SC (low SC = quick runoff)
#   - Daily baseflow proportion = (SC_obs - SC_runoff) / (SC_gw - SC_runoff)
#
# Timeline: Water Years with ≥365 days of concurrent SC and discharge data
#
# Inputs:
#   - CF01201_v3.txt: Continuous specific conductance (EC) data
#   - HF00402_v14.csv: Daily discharge data
#
# Outputs:
#   - Annual_GW_Prop.csv: Annual mean baseflow proportion by site
#
# Author: Keira Johnson (original), Sidney Bush (adapted)
# -----------------------------------------------------------------------------

# Load libraries
library(dplyr)
library(ggplot2)
library(lubridate)

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
# SETUP: Directories (from config.R)
# -----------------------------------------------------------------------------

base_dir <- BASE_DATA_DIR
output_dir <- OUT_MET_MOBILE_DIR

discharge_dir <- DISCHARGE_DIR
ec_dir <- EC_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# LOAD DISCHARGE DATA
# -----------------------------------------------------------------------------

discharge <- read.csv(file.path(discharge_dir, "HF00402_v14.csv")) %>%
  mutate(
    date = as.Date(DATE, "%m/%d/%Y"),
    SITECODE = standardize_site_code(SITECODE)
  ) %>%
  # Filter to chemistry sites and water year range
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    WATERYEAR >= WY_START,
    WATERYEAR <= WY_END
  ) %>%
  select(SITECODE, date, MEAN_Q, WATERYEAR)

# -----------------------------------------------------------------------------
# LOAD & PROCESS SPECIFIC CONDUCTANCE DATA
# -----------------------------------------------------------------------------

EC <- read.delim(file.path(ec_dir, "CF01201_v3.txt"), sep = ",")

# Aggregate to daily mean SC
EC_daily <- EC %>%
  mutate(date = as.Date(DATE_TIME, "%Y-%m-%d")) %>%
  mutate(SITECODE = standardize_site_code(SITECODE)) %>%
  # Filter to chemistry sites
  filter(SITECODE %in% SITE_ORDER_CHEMISTRY) %>%
  group_by(SITECODE, date) %>%
  summarise(daily_SC = mean(EC_INST, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(daily_SC))

# -----------------------------------------------------------------------------
# MERGE SC AND DISCHARGE
# -----------------------------------------------------------------------------

EC_Q <- left_join(
  EC_daily,
  discharge %>% select(SITECODE, date, MEAN_Q, WATERYEAR),
  by = c("SITECODE", "date")
) %>%
  filter(!is.na(MEAN_Q), !is.na(daily_SC))

# -----------------------------------------------------------------------------
# FILTER TO COMPLETE WATER YEARS
# -----------------------------------------------------------------------------

# Get water year for each record
EC_Q <- EC_Q %>%
  mutate(waterYear = get_water_year(date))

# Identify complete water years (≥365 days of data)
goodyears <- EC_Q %>%
  group_by(SITECODE, waterYear) %>%
  summarise(num_days = n_distinct(date), .groups = "drop") %>%
  filter(num_days >= 365)

# Filter to keep only complete water years
EC_Q <- EC_Q %>%
  semi_join(goodyears, by = c("SITECODE", "waterYear"))

# -----------------------------------------------------------------------------
# CALCULATE BASEFLOW PROPORTION
# -----------------------------------------------------------------------------

# Define endmembers (per site, across all years)
EC_Q <- EC_Q %>%
  group_by(SITECODE) %>%
  mutate(
    SC_runoff = quantile(daily_SC, 0.01, na.rm = TRUE), # 1st percentile = quick runoff
    SC_groundwater = quantile(daily_SC, 0.99, na.rm = TRUE) # 99th percentile = groundwater
  ) %>%
  ungroup()

# Two-component mixing model
# Q_baseflow = Q_total × (SC_obs - SC_runoff) / (SC_gw - SC_runoff)
EC_Q <- EC_Q %>%
  mutate(
    Q_baseflow = MEAN_Q * (daily_SC - SC_runoff) / (SC_groundwater - SC_runoff),
    GW_prop = Q_baseflow / MEAN_Q
  )

# Clip to realistic range [0, 1]
EC_Q <- EC_Q %>%
  mutate(GW_prop = pmax(0, pmin(1, GW_prop)))

# -----------------------------------------------------------------------------
# AGGREGATE TO ANNUAL MEAN BASEFLOW
# -----------------------------------------------------------------------------

# CHS = Chemical Hydrograph Separation (method name)
# The metric measured is mean baseflow fraction
annual_bf_prop <- EC_Q %>%
  group_by(SITECODE, waterYear) %>%
  summarise(
    CHS = mean(GW_prop, na.rm = TRUE),  # CHS = mean baseflow fraction
    median_bf = median(GW_prop, na.rm = TRUE),
    sd_bf = sd(GW_prop, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  ) %>%
  # Match legacy output site codes and site set
  mutate(SITECODE = to_legacy_hydro_site_code(SITECODE)) %>%
  filter(SITECODE != "GSWSMC")

# -----------------------------------------------------------------------------
# SAVE OUTPUTS
# -----------------------------------------------------------------------------

# Save annual baseflow proportions
output_file <- file.path(output_dir, "annual_gw_prop.csv")
write.csv(annual_bf_prop, output_file, row.names = FALSE)

# -----------------------------------------------------------------------------
# SUMMARY STATISTICS
# -----------------------------------------------------------------------------

summary_stats <- annual_bf_prop %>%
  group_by(SITECODE) %>%
  summarise(
    n_years = n(),
    CHS_overall = mean(CHS, na.rm = TRUE),
    sd_CHS = sd(CHS, na.rm = TRUE),
    min_CHS = min(CHS, na.rm = TRUE),
    max_CHS = max(CHS, na.rm = TRUE),
    .groups = "drop"
  )
