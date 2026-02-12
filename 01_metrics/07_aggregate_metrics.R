# -----------------------------------------------------------------------------
# Aggregate All Storage Metrics - Master Table Creation
# -----------------------------------------------------------------------------
# Purpose: Combine all hydrometric storage metrics into comprehensive annual
#          and site-averaged master tables for statistical analysis
#
# Storage Metrics by Type:
#   DYNAMIC (from hydrometric data):
#     - RBI: Richards-Baker Index (flashiness)
#     - RCS: Recession Curve Slope
#     - FDC: Flow Duration Curve slope
#     - SD: Storage-Discharge (S_annual_mm)
#
#   MOBILE (from chemistry/isotopes):
#     - CHS: Chemical Hydrograph Separation (mean_bf)
#     - MTT: Mean Transit Time
#     - Fyw: Young Water Fraction
#     - DR: Isotopic Damping Ratio
#
#   EXTENDED DYNAMIC (from water balance):
#     - WB: Water Balance drawdown (DS_sum)
#
# Inputs (all from previous scripts):
#   - RBI_RecessionCurve_Annual.csv: RBI & recession slope
#   - StorageDischarge_FDC_Annual.csv: S_annual, FDC slopes, Q percentiles
#   - Annual_GW_Prop.csv: Mean baseflow from chemical separation
#   - DS_drawdown_annual.csv: Dynamic storage drawdown
#   - stream_thermal_lowflow_metrics_annual.csv: Temperature & low-flow metrics
#   - Catchment_Charc.csv: Topography, geology, landslides
#   - DampingRatios_2025-07-07.csv: Isotope damping ratios
#   - MTT_FYW.csv: Mean transit time and young water fraction
#
# Outputs:
#   - HJA_StorageMetrics_Annual_All.csv: All metrics, annual, all sites
#   - HJA_Ave_StorageMetrics_CatCharacter.csv: Site-averaged metrics + catchment
#   - master_annual.csv: Canonical annual master for all analyses
#   - master_site.csv: Canonical site-level master for all analyses
#   - HJA_Stor_Temp_Yr.csv: Legacy-compatible annual master alias
#   - Sample_Size_by_Metric.csv: Sample size tracking for each metric
#   - QA correlation plots and scatterplot matrices
# -----------------------------------------------------------------------------

# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(GGally)
library(ggcorrplot)

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
# 1. SETUP: Directories (from config.R)
# -----------------------------------------------------------------------------

base_dir    <- BASE_DATA_DIR
dynamic_dir <- OUT_MET_DYNAMIC_DIR
mobile_dir  <- OUT_MET_MOBILE_DIR
extended_dir <- OUT_MET_EXTENDED_DIR
eco_dir <- OUT_MET_ECO_DIR
master_dir <- OUT_MASTER_DIR

# Create output directory if needed
if (!dir.exists(master_dir)) dir.create(master_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. LOAD ALL METRIC FILES
# -----------------------------------------------------------------------------

# --- DYNAMIC STORAGE METRICS ---

# RBI & Recession Slope (RBI, RCS)
rbi_recession <- read_csv(
  file.path(dynamic_dir, "RBI_RecessionCurve_Annual.csv"),
  show_col_types = FALSE
) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, RCS, RBI)

# Storage-Discharge, FDC, Q percentiles (SD, FDC)
storage_fdc <- read_csv(
  file.path(dynamic_dir, "StorageDischarge_FDC_Annual.csv"),
  show_col_types = FALSE
) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, SD, FDC, Q99, Q50, Q01, Q5norm, CV_Q5norm)

# --- MOBILE STORAGE METRICS ---

# Chemical hydrograph separation (CHS = mean baseflow fraction)
baseflow <- read_csv(
  file.path(mobile_dir, "Annual_GW_Prop.csv"),
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, CHS)

# --- EXTENDED DYNAMIC STORAGE METRICS ---

# Water Balance (WB = extended dynamic storage)
wb_storage <- read_csv(
  file.path(extended_dir, "DS_drawdown_annual.csv"),
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, WB)

# --- RESPONSE VARIABLES ---

# Stream temperature & low-flow metrics
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
  select(site, year, max_temp_7d_C, min_Q_7d_mm_d, temp_at_min_Q_7d_C, temp_during_min_Q_7d_C, Q5_CV)

# -----------------------------------------------------------------------------
# 3. MERGE ALL ANNUAL METRICS
# -----------------------------------------------------------------------------

HJA_annual <- rbi_recession %>%
  full_join(storage_fdc, by = c("site", "year")) %>%
  full_join(baseflow, by = c("site", "year")) %>%
  full_join(wb_storage, by = c("site", "year")) %>%
  full_join(thermal_lowflow, by = c("site", "year")) %>%
  # Filter to keep only sites in our analysis
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  # Set factor levels for consistent ordering
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site, year)

# Save annual metrics
write.csv(HJA_annual,
          file.path(master_dir, MASTER_ANNUAL_FILE),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. CALCULATE SITE AVERAGES
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 5. ADD ISOTOPE METRICS (MTT, Fyw, DR) - SITE-LEVEL
# -----------------------------------------------------------------------------

# Mean Transit Time and Young Water Fraction
mtt_fyw <- read_csv(
  file.path(base_dir, "Isotopes", "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  select(site, MTT = MTTM, Fyw = FYWM) %>%
  filter(!is.na(site), site != "")

# Damping ratios
damping_ratios <- read_csv(
  file.path(base_dir, "Isotopes", "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  select(site, DR = DR_Overall, DR_err = DR__err)

# Merge isotope metrics
isotope_metrics <- mtt_fyw %>%
  full_join(damping_ratios, by = "site")

# Add isotope metrics to site averages
HJA_avg <- HJA_avg %>%
  left_join(isotope_metrics, by = "site")

# -----------------------------------------------------------------------------
# 6. ADD CATCHMENT CHARACTERISTICS
# -----------------------------------------------------------------------------

catchment_chars <- read_csv(
  file.path(base_dir, "DynamicStorage", "Catchment_Charc.csv"),
  show_col_types = FALSE
) %>%
  rename(site = Site) %>%
  mutate(site = standardize_site_code(site))

HJA_avg <- HJA_avg %>%
  left_join(catchment_chars, by = "site")

# Save site-averaged metrics with all data
write.csv(HJA_avg,
          file.path(master_dir, MASTER_SITE_FILE),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. SAMPLE SIZE TRACKING
# -----------------------------------------------------------------------------

# Define all storage metrics by type (using method abbreviations)
dynamic_metrics <- c("RBI", "RCS", "FDC", "SD")
mobile_metrics <- c("CHS", "MTT", "Fyw", "DR")
extended_metrics <- c("WB")
response_metrics <- c("max_temp_7d_C", "min_Q_7d_mm_d", "temp_during_min_Q_7d_C")

# Count non-NA values for annual metrics
annual_sample_sizes <- HJA_annual %>%
  group_by(site) %>%
  summarise(
    n_years_total = n(),
    n_RBI = sum(!is.na(RBI)),
    n_RCS = sum(!is.na(RCS)),
    n_FDC = sum(!is.na(FDC)),
    n_SD = sum(!is.na(SD)),
    n_CHS = sum(!is.na(CHS)),
    n_WB = sum(!is.na(WB)),
    n_max_temp = sum(!is.na(max_temp_7d_C)),
    n_min_Q = sum(!is.na(min_Q_7d_mm_d)),
    .groups = "drop"
  )

# Add site-level isotope availability
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

# Save sample size table

# -----------------------------------------------------------------------------
# 8. QA: CORRELATION MATRIX (STORAGE METRICS ONLY)
# -----------------------------------------------------------------------------

# Storage metrics by type (Q5norm, CV_Q5norm are NOT storage metrics)
#   Dynamic: RBI, RCS, FDC, SD
#   Mobile: CHS, MTT, Fyw, DR
#   Extended Dynamic: WB
# -----------------------------------------------------------------------------
# 8. FINAL SUMMARY
# -----------------------------------------------------------------------------

cat("\n=== AGGREGATION COMPLETE ===\n")
cat(sprintf("Water years: %d - %d\n", WY_START, WY_END))
cat(sprintf("Hydrometric sites: %d\n", length(unique(HJA_annual$site))))
cat(sprintf("Annual observations: %d\n", nrow(HJA_annual)))
cat(sprintf("Site averages: %d\n", nrow(HJA_avg)))
