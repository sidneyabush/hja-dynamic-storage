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

base_dir    <- BASE_DATA_DIR
dynamic_dir <- OUT_MET_DYNAMIC_DIR
mobile_dir  <- OUT_MET_MOBILE_DIR
extended_dir <- OUT_MET_EXTENDED_DIR
eco_dir <- OUT_MET_ECO_DIR
master_dir <- OUT_MASTER_DIR
isotope_dir <- ISOTOPE_DIR
catchment_dir <- CATCHMENT_CHARACTERISTICS_DIR

# Create output directory if needed
if (!dir.exists(master_dir)) dir.create(master_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# LOAD ALL METRIC FILES
# -----------------------------------------------------------------------------

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

# --- DYNAMIC STORAGE METRICS ---

# RBI & Recession Slope (RBI, RCS)
rbi_path <- file.path(dynamic_dir, "rbi_rcs_annual.csv")
if (!file.exists(rbi_path)) {
  rbi_path <- file.path(dynamic_dir, "RBI_RecessionCurve_Annual.csv")
}
rbi_recession <- read_csv(
  rbi_path,
  show_col_types = FALSE
) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, RCS, RBI)
assert_unique_keys(rbi_recession, c("site", "year"), "rbi_recession")

# Storage-Discharge, FDC, Q percentiles (SD, FDC)
fdc_path <- file.path(dynamic_dir, "storage_discharge_fdc_annual.csv")
if (!file.exists(fdc_path)) {
  fdc_path <- file.path(dynamic_dir, "StorageDischarge_FDC_Annual.csv")
}
storage_fdc <- read_csv(
  fdc_path,
  show_col_types = FALSE
) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, SD, FDC, Q99, Q50, Q01, Q5norm, CV_Q5norm)
assert_unique_keys(storage_fdc, c("site", "year"), "storage_fdc")

# --- MOBILE STORAGE METRICS ---

# Chemical hydrograph separation (CHS = mean baseflow fraction)
chs_path <- file.path(mobile_dir, "annual_gw_prop.csv")
if (!file.exists(chs_path)) {
  chs_path <- file.path(mobile_dir, "Annual_GW_Prop.csv")
}
baseflow <- read_csv(
  chs_path,
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, CHS)
assert_unique_keys(baseflow, c("site", "year"), "baseflow")

# --- EXTENDED DYNAMIC STORAGE METRICS ---

# Water Balance (WB = extended dynamic storage)
wb_path <- file.path(extended_dir, "ds_drawdown_annual.csv")
if (!file.exists(wb_path)) {
  wb_path <- file.path(extended_dir, "DS_drawdown_annual.csv")
}
wb_storage <- read_csv(
  wb_path,
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, WB)
assert_unique_keys(wb_storage, c("site", "year"), "wb_storage")

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

# -----------------------------------------------------------------------------
# MERGE ALL ANNUAL METRICS
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
assert_unique_keys(HJA_annual, c("site", "year"), "HJA_annual")

# Save annual metrics
write.csv(HJA_annual,
          file.path(master_dir, MASTER_ANNUAL_FILE),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# CALCULATE SITE AVERAGES
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
# ADD ISOTOPE METRICS (MTT, Fyw, DR) - SITE-LEVEL
# -----------------------------------------------------------------------------

# Mean Transit Time and Young Water Fraction
mtt_fyw <- read_csv(
  file.path(isotope_dir, "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  select(site, MTT = MTTM, Fyw = FYWM) %>%
  filter(!is.na(site), site != "", site %in% SITE_ORDER_HYDROMETRIC)

# Damping ratios
damping_ratios <- read_csv(
  file.path(isotope_dir, "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  select(site, DR = DR_Overall, DR_err = DR__err) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

# Merge isotope metrics
isotope_metrics <- mtt_fyw %>%
  full_join(damping_ratios, by = "site")
assert_unique_keys(isotope_metrics, c("site"), "isotope_metrics")

# Add isotope metrics to site averages
HJA_avg <- HJA_avg %>%
  left_join(isotope_metrics, by = "site")

# -----------------------------------------------------------------------------
# ADD CATCHMENT CHARACTERISTICS
# -----------------------------------------------------------------------------

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

# Save site-averaged metrics with all data
write.csv(HJA_avg,
          file.path(master_dir, MASTER_SITE_FILE),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# SAMPLE SIZE TRACKING
# -----------------------------------------------------------------------------

# Define all storage metrics by type (using method abbreviations)
dynamic_metrics <- c("RBI", "RCS", "FDC", "SD")
mobile_metrics <- c("CHS", "MTT", "Fyw", "DR")
extended_metrics <- c("WB")
response_metrics <- c("T_7DMax", "Q_7Q5", "T_Q7Q5")

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
    n_t_7dmax = sum(!is.na(T_7DMax)),
    n_q_7q5 = sum(!is.na(Q_7Q5)),
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
# QA: CORRELATION MATRIX (STORAGE METRICS ONLY)
# -----------------------------------------------------------------------------

# Storage metrics by type (Q5norm, CV_Q5norm are NOT storage metrics)
#   Dynamic: RBI, RCS, FDC, SD
#   Mobile: CHS, MTT, Fyw, DR
#   Extended Dynamic: WB
# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------

