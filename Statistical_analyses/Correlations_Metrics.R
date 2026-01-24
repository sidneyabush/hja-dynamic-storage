# =============================================================================
# Correlation Analysis: Storage Metrics & Catchment Attributes
# =============================================================================
# Purpose: Develop correlations and relationships among storage metrics
#          and with catchment attributes
#
# Analyses:
#   1. Compute site-averaged storage metrics from annual data
#   2. Merge with catchment characteristics and temperature data
#   3. Generate correlation matrices for:
#      - Storage metrics only
#      - Catchment attributes only
#      - Combined storage + catchment analysis
#   4. Create correlation plots (ggcorrplot)
#
# Inputs:
#   - HJA_StorageMetrics_Annual_All.csv: Annual storage metrics
#   - Catchment_Charc.csv: Catchment characteristics
#   - mean_july.temp.csv: Mean July stream temperature
#   - mean_july_airtemp.csv: Mean July air temperature
#   - DampingRatios_2025-07-07.csv: Isotope damping ratios
#
# Outputs:
#   - HJA_Stor_Temp_Yr.csv: Annual storage + temperature data
#   - HJA_Ave_StorageMetrics_CatCharacter.csv: Site-averaged metrics
#   - HJA_Storage_CorrPlot_CatchM.pdf: Catchment correlation plot
#   - HJA_Storage_CorrPlot.pdf: Full correlation plot
#
# Author: Pamela Sullivan (original), Sidney Bush (adapted)
# Date: 2026-01-23
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(colorspace)
library(tidyr)
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
# 2. LOAD CATCHMENT ATTRIBUTES
# =============================================================================

Catt <- read_csv(
  file.path(base_dir, "DynamicStorage", "Catchment_Charc.csv"),
  show_col_types = FALSE
) %>%
  rename(site = Site)

# =============================================================================
# 3. LOAD TEMPERATURE DATA
# =============================================================================

# Mean stream water temperature in July
WtempM <- read_csv(
  file.path(base_dir, "Stream_T/Output", "mean_july.temp.csv"),
  show_col_types = FALSE
) %>%
  rename(year = Year) %>%
  pivot_longer(
    cols = -year,
    names_to = "site",
    values_to = "JulyM_ST"
  ) %>%
  mutate(site = if_else(site == "GSMACK", "GSWSMC", site))

# Mean air temperature in July
AtempM <- read_csv(
  file.path(base_dir, "Stream_T/Output", "mean_july_airtemp.csv"),
  show_col_types = FALSE
) %>%
  rename(year = Year) %>%
  pivot_longer(
    cols = -year,
    names_to = "site",
    values_to = "JulyM_AT"
  ) %>%
  mutate(site = if_else(site == "GSMACK", "GSWSMC", site))

# Combine and calculate stream/air temperature ratio
HJA_Temp <- left_join(WtempM, AtempM, by = c("site", "year")) %>%
  mutate(JST_AT = JulyM_ST / JulyM_AT)

# =============================================================================
# 4. LOAD STORAGE METRICS
# =============================================================================

HJA_storage <- read_csv(
  file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv"),
  show_col_types = FALSE,
  col_select = -1
)

# =============================================================================
# 5. MERGE STORAGE + TEMPERATURE (ANNUAL)
# =============================================================================

HJA_Stor_Temp_Yr <- left_join(HJA_storage, HJA_Temp, by = c("site", "year"))

# Save yearly data
write.csv(HJA_Stor_Temp_Yr,
          file.path(output_dir, "HJA_Stor_Temp_Yr.csv"),
          row.names = FALSE)

# =============================================================================
# 6. LOAD DAMPING RATIOS (ISOTOPES)
# =============================================================================

Damp_R <- read_csv(
  file.path(base_dir, "Isotopes", "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = if_else(site == "GSMACK", "GSWSMC", site))

# =============================================================================
# 7. CALCULATE SITE-AVERAGED STORAGE METRICS
# =============================================================================

# Columns to average
cols <- c(
  "recession_curve_slope", "RBI", "Q5norm", "CV_Q5norm", "mean_bf",
  "fdc_slope", "S_annual_mm", "DS_sum", "JulyM_ST", "JulyM_AT", "JST_AT"
)

HJA_storage_ave <- HJA_Stor_Temp_Yr %>%
  group_by(site) %>%
  summarise(
    across(
      all_of(cols),
      list(mean = ~mean(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# =============================================================================
# 8. MERGE WITH DAMPING RATIOS & CATCHMENT ATTRIBUTES
# =============================================================================

HJA_storage_ave <- HJA_storage_ave %>%
  left_join(Damp_R, by = "site") %>%
  left_join(Catt, by = "site")

# Save overall site-based values
write.csv(HJA_storage_ave,
          file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
          row.names = FALSE)

# =============================================================================
# 9. CORRELATION ANALYSIS: CATCHMENT ATTRIBUTES ONLY
# =============================================================================

vars_catchment <- c(
  "Area_km2", "Elevation_mean_m", "Slope_mean", "Aspec_Mean_deg",
  "Harvest", "Landslide_Young", "Landslide_Mod", "Landslide_Old",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
)

cor_catchment <- cor(HJA_storage_ave[vars_catchment],
                     use = "pairwise.complete.obs")

# Plot catchment correlation matrix
p1 <- ggcorrplot(cor_catchment,
           hc.order = FALSE,
           type = "lower",
           outline.col = "white",
           lab = TRUE) +
  labs(title = "Catchment Attributes Correlation Matrix")

ggsave(
  file.path(output_dir, "QA_HJA_Storage_CorrPlot_CatchM.png"),
  p1, width = 10, height = 10, dpi = 300
)

# =============================================================================
# 10. CORRELATION ANALYSIS: STORAGE + CATCHMENT (FULL)
# =============================================================================

vars_all <- c(
  "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean",
  "mean_bf_mean", "fdc_slope_mean", "S_annual_mm_mean", "DS_sum_mean",
  "DR_Overall", "JST_AT_mean", "Area_km2", "Elevation_mean_m",
  "Slope_mean", "Aspec_Mean_deg", "Harvest", "Landslide_Young",
  "Landslide_Mod", "Landslide_Old", "Lava1_per", "Lava2_per",
  "Ash_Per", "Pyro_per"
)

cor_all <- cor(HJA_storage_ave[vars_all], use = "pairwise.complete.obs")

# Plot full correlation matrix
p2 <- ggcorrplot(cor_all,
           hc.order = FALSE,
           type = "lower",
           outline.col = "white",
           lab = TRUE) +
  labs(title = "Storage Metrics & Catchment Attributes Correlation Matrix")

ggsave(
  file.path(output_dir, "QA_HJA_Storage_CorrPlot_Full.png"),
  p2, width = 14, height = 14, dpi = 300
)

# =============================================================================
# 11. CORRELATION ANALYSIS: STORAGE METRICS ONLY
# =============================================================================

vars_storage <- c(
  "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean",
  "mean_bf_mean", "fdc_slope_mean", "S_annual_mm_mean", "DS_sum_mean",
  "DR_Overall", "JST_AT_mean"
)

cor_storage <- cor(HJA_storage_ave[vars_storage], use = "pairwise.complete.obs")

# Plot storage metrics correlation matrix
p3 <- ggcorrplot(cor_storage,
           hc.order = FALSE,
           type = "lower",
           outline.col = "white",
           lab = TRUE) +
  labs(title = "Storage Metrics Correlation Matrix")

ggsave(
  file.path(output_dir, "QA_HJA_Storage_CorrPlot_Metrics.png"),
  p3, width = 12, height = 12, dpi = 300
)
