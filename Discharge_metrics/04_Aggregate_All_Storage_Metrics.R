# =============================================================================
# Aggregate All Storage Metrics - Master Table Creation
# =============================================================================
# Purpose: Combine all hydrometric storage metrics into comprehensive annual
#          and site-averaged master tables for statistical analysis
#
# Inputs (all from previous scripts):
#   - RBI_RecessionCurve_Annual.csv: RBI & recession slope
#   - StorageDischarge_FDC_Annual.csv: S_annual, FDC slopes, Q percentiles
#   - Annual_GW_Prop.csv: Mean baseflow from chemical separation
#   - DS_drawdown_annual.csv: Dynamic storage drawdown
#   - stream_thermal_lowflow_metrics_annual.csv: Temperature & low-flow metrics
#   - Catchment_Charc.csv: Topography, geology, landslides
#   - DampingRatios_2025-07-07.csv: Isotope damping ratios
#
# Outputs:
#   - HJA_StorageMetrics_Annual_All.csv: All metrics, annual, all sites
#   - HJA_StorageMetrics_Average_All.csv: Site-averaged metrics
#   - QA correlation plots and scatterplot matrices
#
# Author: Keira Johnson (original), Sidney Bush (adapted)
# Date: 2026-01-23
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(GGally)
library(ggcorrplot)

theme_set(theme_classic(base_size = 12))

# Clear environment
rm(list = ls())

# =============================================================================
# 1. SETUP: Directories
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir  <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 2. LOAD ALL METRIC FILES
# =============================================================================

# RBI & Recession Slope
rbi_recession <- read_csv(
  file.path(output_dir, "RBI_RecessionCurve_Annual.csv"),
  show_col_types = FALSE
) %>%
  select(site, year, recession_curve_slope, RBI)

# Storage-Discharge, FDC, Q percentiles
storage_fdc <- read_csv(
  file.path(output_dir, "StorageDischarge_FDC_Annual.csv"),
  show_col_types = FALSE
) %>%
  select(site, year, S_annual_mm, fdc_slope, Q99, Q50, Q01, Q5norm, CV_Q5norm)

# Chemical hydrograph separation (baseflow)
baseflow <- read_csv(
  file.path(output_dir, "Annual_GW_Prop.csv"),
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  select(site, year, mean_bf)

# Dynamic storage drawdown
ds_drawdown <- read_csv(
  file.path(output_dir, "DS_drawdown_annual.csv"),
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  select(site, year, DS_sum)

# Stream temperature & low-flow metrics
thermal_lowflow <- read_csv(
  file.path(output_dir, "stream_thermal_lowflow_metrics_annual.csv"),
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = wateryear) %>%
  select(site, year, max_temp_7d_C, min_Q_7d_mm_d, temp_at_min_Q_7d_C, Q5_CV)

# =============================================================================
# 3. MERGE ALL ANNUAL METRICS
# =============================================================================

HJA_annual <- rbi_recession %>%
  full_join(storage_fdc, by = c("site", "year")) %>%
  full_join(baseflow, by = c("site", "year")) %>%
  full_join(ds_drawdown, by = c("site", "year")) %>%
  full_join(thermal_lowflow, by = c("site", "year")) %>%
  arrange(site, year)

# Save annual metrics
write.csv(HJA_annual,
          file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv"),
          row.names = FALSE)

# =============================================================================
# 4. CALCULATE SITE AVERAGES
# =============================================================================

HJA_avg <- HJA_annual %>%
  group_by(site) %>%
  summarise(
    across(
      where(is.numeric),
      list(mean = ~mean(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# =============================================================================
# 5. ADD CATCHMENT CHARACTERISTICS
# =============================================================================

# Topography
catchment_chars <- read_csv(
  file.path(base_dir, "DynamicStorage", "Catchment_Charc.csv"),
  show_col_types = FALSE
) %>%
  rename(site = Site)

HJA_avg <- HJA_avg %>%
  left_join(catchment_chars, by = "site")

# Damping ratios (isotope data)
damping_ratios <- read_csv(
  file.path(base_dir, "Isotopes", "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = if_else(site == "GSMACK", "GSWSMC", site))

HJA_avg <- HJA_avg %>%
  left_join(damping_ratios, by = "site")

# Save site-averaged metrics
write.csv(HJA_avg,
          file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
          row.names = FALSE)

# =============================================================================
# 6. QA: CORRELATION MATRIX (STORAGE METRICS ONLY)
# =============================================================================

storage_vars <- c(
  "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean",
  "mean_bf_mean", "fdc_slope_mean", "S_annual_mm_mean", "DS_sum_mean"
)

cor_storage <- cor(HJA_avg[storage_vars], use = "pairwise.complete.obs")

p_cor <- ggcorrplot(cor_storage,
           hc.order = FALSE,
           type = "lower",
           outline.col = "white",
           lab = TRUE) +
  labs(title = "Storage Metrics Correlation Matrix (Site Averages)")

ggsave(
  file.path(output_dir, "QA_Storage_Metrics_CorrPlot.png"),
  p_cor, width = 10, height = 10, dpi = 300
)

# =============================================================================
# 7. QA: SCATTERPLOT MATRIX (STORAGE METRICS)
# =============================================================================

no_se_smoother <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_smooth(se = FALSE, method = "lm", color = "black", ...) +
    geom_point()
}

# Subset to key storage metrics
key_metrics <- HJA_avg %>%
  select(all_of(storage_vars)) %>%
  rename_with(~gsub("_mean$", "", .x))  # Remove "_mean" suffix for cleaner labels

p_pairs <- ggpairs(
  key_metrics,
  lower = list(continuous = no_se_smoother)
) +
  theme_bw()

ggsave(
  file.path(output_dir, "QA_Storage_Metrics_Scatterplot_Matrix.png"),
  p_pairs, width = 14, height = 14, dpi = 300
)

# =============================================================================
# 8. SUMMARY STATISTICS
# =============================================================================

summary_stats <- HJA_annual %>%
  group_by(site) %>%
  summarise(
    n_years = n(),
    RBI_mean = mean(RBI, na.rm = TRUE),
    RBI_sd = sd(RBI, na.rm = TRUE),
    Q5norm_mean = mean(Q5norm, na.rm = TRUE),
    Q5norm_sd = sd(Q5norm, na.rm = TRUE),
    mean_bf_mean = mean(mean_bf, na.rm = TRUE),
    mean_bf_sd = sd(mean_bf, na.rm = TRUE),
    .groups = "drop"
  )
