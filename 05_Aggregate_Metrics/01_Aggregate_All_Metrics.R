# =============================================================================
# Aggregate All Storage Metrics - Master Table Creation
# =============================================================================
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
#   - Sample_Size_by_Metric.csv: Sample size tracking for each metric
#   - QA correlation plots and scatterplot matrices
# =============================================================================

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
# 1. SETUP: Directories (from config.R)
# =============================================================================

base_dir    <- BASE_DATA_DIR
output_dir  <- OUTPUT_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 2. LOAD ALL METRIC FILES
# =============================================================================

# --- DYNAMIC STORAGE METRICS ---

# RBI & Recession Slope (RBI, RCS)
rbi_recession <- read_csv(
  file.path(output_dir, "RBI_RecessionCurve_Annual.csv"),
  show_col_types = FALSE
) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, recession_curve_slope, RBI)

# Storage-Discharge, FDC, Q percentiles (SD, FDC)
storage_fdc <- read_csv(
  file.path(output_dir, "StorageDischarge_FDC_Annual.csv"),
  show_col_types = FALSE
) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, S_annual_mm, fdc_slope, Q99, Q50, Q01, Q5norm, CV_Q5norm)

# --- MOBILE STORAGE METRICS ---

# Chemical hydrograph separation (CHS -> mean_bf)
baseflow <- read_csv(
  file.path(output_dir, "Annual_GW_Prop.csv"),
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, mean_bf)

# --- EXTENDED DYNAMIC STORAGE METRICS ---

# Dynamic storage drawdown (WB -> DS_sum)
ds_drawdown <- read_csv(
  file.path(output_dir, "DS_drawdown_annual.csv"),
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, DS_sum)

# --- RESPONSE VARIABLES ---

# Stream temperature & low-flow metrics
thermal_lowflow <- read_csv(
  file.path(output_dir, "stream_thermal_lowflow_metrics_annual.csv"),
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = wateryear) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, max_temp_7d_C, min_Q_7d_mm_d, temp_at_min_Q_7d_C, Q5_CV)

# =============================================================================
# 3. MERGE ALL ANNUAL METRICS
# =============================================================================

HJA_annual <- rbi_recession %>%
  full_join(storage_fdc, by = c("site", "year")) %>%
  full_join(baseflow, by = c("site", "year")) %>%
  full_join(ds_drawdown, by = c("site", "year")) %>%
  full_join(thermal_lowflow, by = c("site", "year")) %>%
  # Filter to keep only sites in our analysis
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  # Set factor levels for consistent ordering
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site, year)

# Save annual metrics
write.csv(HJA_annual,
          file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv"),
          row.names = FALSE)

cat("Saved: HJA_StorageMetrics_Annual_All.csv\n")

# =============================================================================
# 4. CALCULATE SITE AVERAGES
# =============================================================================

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

# =============================================================================
# 5. ADD ISOTOPE METRICS (MTT, Fyw, DR) - SITE-LEVEL
# =============================================================================

# Mean Transit Time and Young Water Fraction
mtt_fyw <- read_csv(
  file.path(base_dir, "Isotopes", "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = case_when(
      site == "MCRAEC" ~ "MR",        # McRae Creek
      site == "GSLOOK " ~ "GSLOOK",   # Remove trailing space
      TRUE ~ site
    )
  ) %>%
  select(site, MTT = MTTM, Fyw = FYWM) %>%
  filter(!is.na(site), site != "")

# Damping ratios
damping_ratios <- read_csv(
  file.path(base_dir, "Isotopes", "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = case_when(
      site == "GSMACK" ~ "GSWSMC",
      TRUE ~ site
    )
  ) %>%
  select(site, DR = DR_Overall, DR_err = DR__err)

# Merge isotope metrics
isotope_metrics <- mtt_fyw %>%
  full_join(damping_ratios, by = "site")

# Add isotope metrics to site averages
HJA_avg <- HJA_avg %>%
  left_join(isotope_metrics, by = "site")

# =============================================================================
# 6. ADD CATCHMENT CHARACTERISTICS
# =============================================================================

catchment_chars <- read_csv(
  file.path(base_dir, "DynamicStorage", "Catchment_Charc.csv"),
  show_col_types = FALSE
) %>%
  rename(site = Site)

HJA_avg <- HJA_avg %>%
  left_join(catchment_chars, by = "site")

# Save site-averaged metrics with all data
write.csv(HJA_avg,
          file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
          row.names = FALSE)

cat("Saved: HJA_Ave_StorageMetrics_CatCharacter.csv\n")

# =============================================================================
# 7. SAMPLE SIZE TRACKING
# =============================================================================

# Define all storage metrics by type
dynamic_metrics <- c("RBI", "recession_curve_slope", "fdc_slope", "S_annual_mm")
mobile_metrics <- c("mean_bf", "MTT", "Fyw", "DR")
extended_metrics <- c("DS_sum")
response_metrics <- c("max_temp_7d_C", "min_Q_7d_mm_d", "temp_at_min_Q_7d_C")

# Count non-NA values for annual metrics
annual_sample_sizes <- HJA_annual %>%
  group_by(site) %>%
  summarise(
    n_years_total = n(),
    n_RBI = sum(!is.na(RBI)),
    n_RCS = sum(!is.na(recession_curve_slope)),
    n_FDC = sum(!is.na(fdc_slope)),
    n_SD = sum(!is.na(S_annual_mm)),
    n_CHS = sum(!is.na(mean_bf)),
    n_WB = sum(!is.na(DS_sum)),
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
write.csv(sample_sizes,
          file.path(output_dir, "Sample_Size_by_Metric.csv"),
          row.names = FALSE)

cat("Saved: Sample_Size_by_Metric.csv\n")

# Print summary
cat("\n=== SAMPLE SIZE SUMMARY ===\n\n")
cat("Storage Type | Metric | Sites | Total Obs\n")
cat("-------------|--------|-------|----------\n")
cat(sprintf("Dynamic      | RBI    | %d     | %d\n",
            sum(sample_sizes$n_RBI > 0), sum(sample_sizes$n_RBI)))
cat(sprintf("Dynamic      | RCS    | %d     | %d\n",
            sum(sample_sizes$n_RCS > 0), sum(sample_sizes$n_RCS)))
cat(sprintf("Dynamic      | FDC    | %d     | %d\n",
            sum(sample_sizes$n_FDC > 0), sum(sample_sizes$n_FDC)))
cat(sprintf("Dynamic      | SD     | %d     | %d\n",
            sum(sample_sizes$n_SD > 0), sum(sample_sizes$n_SD)))
cat(sprintf("Mobile       | CHS    | %d     | %d\n",
            sum(sample_sizes$n_CHS > 0), sum(sample_sizes$n_CHS)))
cat(sprintf("Mobile       | MTT    | %d     | (site-level)\n",
            sum(sample_sizes$has_MTT)))
cat(sprintf("Mobile       | Fyw    | %d     | (site-level)\n",
            sum(sample_sizes$has_Fyw)))
cat(sprintf("Mobile       | DR     | %d     | (site-level)\n",
            sum(sample_sizes$has_DR)))
cat(sprintf("Extended     | WB     | %d     | %d\n",
            sum(sample_sizes$n_WB > 0), sum(sample_sizes$n_WB)))

# =============================================================================
# 8. QA: CORRELATION MATRIX (STORAGE METRICS ONLY)
# =============================================================================

# Storage metrics by type (Q5norm, CV_Q5norm are NOT storage metrics)
#   Dynamic: RBI, RCS, FDC, SD
#   Mobile: CHS (mean_bf), MTT, Fyw, DR
#   Extended Dynamic: WB (DS_sum)
storage_vars <- c(
  "recession_curve_slope_mean",  # RCS - Dynamic
  "RBI_mean",                    # RBI - Dynamic
  "fdc_slope_mean",              # FDC - Dynamic
  "S_annual_mm_mean",            # SD  - Dynamic
  "mean_bf_mean",                # CHS - Mobile
  "MTT",                         # MTT - Mobile (site-level)
  "Fyw",                         # Fyw - Mobile (site-level)
  "DR",                          # DR  - Mobile (site-level)
  "DS_sum_mean"                  # WB  - Extended Dynamic
)

# Filter to available columns
available_vars <- storage_vars[storage_vars %in% colnames(HJA_avg)]

cor_data <- HJA_avg %>%
  select(all_of(available_vars)) %>%
  filter(complete.cases(.))

if (nrow(cor_data) >= 3) {
  cor_storage <- cor(cor_data, use = "pairwise.complete.obs")

  p_cor <- ggcorrplot(cor_storage,
             hc.order = FALSE,
             type = "lower",
             outline.col = "white",
             lab = TRUE,
             lab_size = 2.5) +
    labs(title = "Storage Metrics Correlation Matrix (Site Averages)")

  ggsave(
    file.path(output_dir, "QA_Storage_Metrics_CorrPlot.png"),
    p_cor, width = 12, height = 12, dpi = 300
  )
  cat("\nSaved: QA_Storage_Metrics_CorrPlot.png\n")
}

# =============================================================================
# 9. QA: SCATTERPLOT MATRIX (KEY STORAGE METRICS)
# =============================================================================

no_se_smoother <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_smooth(se = FALSE, method = "lm", color = "black", ...) +
    geom_point()
}

# Subset to key storage metrics (exclude _mean suffix for cleaner labels)
key_vars <- c("RBI_mean", "fdc_slope_mean", "S_annual_mm_mean",
              "DS_sum_mean", "mean_bf_mean", "DR")
key_vars <- key_vars[key_vars %in% colnames(HJA_avg)]

key_metrics <- HJA_avg %>%
  select(all_of(key_vars)) %>%
  rename_with(~gsub("_mean$", "", .x))

if (ncol(key_metrics) >= 3) {
  p_pairs <- ggpairs(
    key_metrics,
    lower = list(continuous = no_se_smoother)
  ) +
    theme_bw()

  ggsave(
    file.path(output_dir, "QA_Storage_Metrics_Scatterplot_Matrix.png"),
    p_pairs, width = 14, height = 14, dpi = 300
  )
  cat("Saved: QA_Storage_Metrics_Scatterplot_Matrix.png\n")
}

# =============================================================================
# 10. FINAL SUMMARY
# =============================================================================

cat("\n=== AGGREGATION COMPLETE ===\n")
cat(sprintf("Water years: %d - %d\n", WY_START, WY_END))
cat(sprintf("Hydrometric sites: %d\n", length(unique(HJA_annual$site))))
cat(sprintf("Annual observations: %d\n", nrow(HJA_annual)))
cat(sprintf("Site averages: %d\n", nrow(HJA_avg)))
