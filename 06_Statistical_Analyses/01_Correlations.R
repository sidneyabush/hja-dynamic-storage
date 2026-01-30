# =============================================================================
# Correlation Analysis: Storage Metrics & Catchment Attributes
# =============================================================================
# Purpose: Generate correlation matrices for storage metrics and catchment
#          attributes to explore relationships
#
# Analyses:
#   1. Storage metrics only - identify relationships among hydrometric indicators
#   2. Catchment attributes only - examine landscape controls
#   3. Combined analysis - explore storage-landscape relationships
#
# Inputs:
#   - HJA_Ave_StorageMetrics_CatCharacter.csv: Site-averaged metrics
#     (created by 06_Aggregate_All_Metrics.R)
#
# Outputs:
#   - QA_Storage_Metrics_CorrPlot.png: Storage metrics correlation matrix
#   - QA_Catchment_Attributes_CorrPlot.png: Catchment correlation matrix
#   - QA_Storage_Catchment_Combined_CorrPlot.png: Combined correlation matrix
#
# Author: Sidney Bush
# Date: 2026-01-24
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(ggcorrplot)

theme_set(theme_classic(base_size = 12))

# Clear environment
rm(list = ls())

# Source configuration (paths, site definitions, water year range)
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir) || script_dir == "") script_dir <- getwd()
config_path <- file.path(dirname(script_dir), "config.R")
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found. Please ensure config.R exists in the repo root.")
}

# Use configuration values
site_order <- SITE_ORDER_HYDROMETRIC
output_dir <- OUTPUT_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 2. LOAD SITE-AVERAGED DATA
# =============================================================================
# This file was created by 06_Aggregate_All_Metrics.R and contains:
# - Storage metrics (averaged across years)
# - Catchment characteristics
# - Temperature metrics
# - Damping ratios

HJA_Ave <- read_csv(
  file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
  show_col_types = FALSE
) %>%
  filter(!site %in% c("GSLOOK_FULL", "GSWSMA", "GSWSMF", "GSMACK"))  # Exclude non-analysis sites

# =============================================================================
# 3. CORRELATION ANALYSIS: STORAGE METRICS ONLY
# =============================================================================

storage_vars <- c(
  "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean",
  "mean_bf_mean", "fdc_slope_mean", "S_annual_mm_mean", "DS_sum_mean"
)

cor_storage <- cor(HJA_Ave[storage_vars], use = "pairwise.complete.obs")

p_storage <- ggcorrplot(cor_storage,
           hc.order = FALSE,
           type = "lower",
           outline.col = "white",
           lab = TRUE) +
  labs(title = "Storage Metrics Correlation Matrix (Site Averages)")

ggsave(
  file.path(output_dir, "QA_Storage_Metrics_CorrPlot.png"),
  p_storage, width = 10, height = 10, dpi = 300
)

# =============================================================================
# 4. CORRELATION ANALYSIS: CATCHMENT ATTRIBUTES ONLY
# =============================================================================

vars_catchment <- c(
  "Area_km2", "Elevation_mean_m", "Slope_mean", "Aspec_Mean_deg",
  "Harvest", "Landslide_Young", "Landslide_Mod", "Landslide_Old",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
)

cor_catchment <- cor(HJA_Ave[vars_catchment], use = "pairwise.complete.obs")

p_catchment <- ggcorrplot(cor_catchment,
           hc.order = FALSE,
           type = "lower",
           outline.col = "white",
           lab = TRUE) +
  labs(title = "Catchment Attributes Correlation Matrix")

ggsave(
  file.path(output_dir, "QA_Catchment_Attributes_CorrPlot.png"),
  p_catchment, width = 10, height = 10, dpi = 300
)

# =============================================================================
# 5. CORRELATION ANALYSIS: COMBINED (STORAGE + CATCHMENT)
# =============================================================================

vars_combined <- c(storage_vars, vars_catchment)

cor_combined <- cor(HJA_Ave[vars_combined], use = "pairwise.complete.obs")

p_combined <- ggcorrplot(cor_combined,
           hc.order = FALSE,
           type = "lower",
           outline.col = "white",
           lab = FALSE) +  # Too many variables for labels
  labs(title = "Storage Metrics & Catchment Attributes Correlation Matrix")

ggsave(
  file.path(output_dir, "QA_Storage_Catchment_Combined_CorrPlot.png"),
  p_combined, width = 14, height = 14, dpi = 300
)
