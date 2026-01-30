# =============================================================================
# Organized Storage Metrics Figures
# =============================================================================
# Purpose: Create properly organized faceted plots for storage metrics
#
# Groups:
#   1. Dynamic Storage (RBI, RCS, FDC, SD) - faceted plot
#   2. Extended Dynamic (WB) - single panel
#
# Mobile storage plots are separate:
#   - Mobile_Isotope_Metrics.png (DR, Fyw, MTT1, MTT2)
#   - Mobile_CHS_Boxplot.png (CHS/baseflow)
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)

# Clear environment
rm(list = ls())

# =============================================================================
# SOURCE CONFIGURATION
# =============================================================================

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
  script_dir <- file.path(getwd(), "07_Plots")
}

config_path <- file.path(dirname(script_dir), "config.R")
if (!file.exists(config_path)) {
  config_path <- file.path(getwd(), "config.R")
}
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found.")
}

# =============================================================================
# SETUP
# =============================================================================

base_dir <- BASE_DATA_DIR
output_dir <- file.path(FIGURES_DIR, "Hydrometric")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Site colors and order from config
site_colors <- SITE_COLORS
site_order <- SITE_ORDER_HYDROMETRIC

# Metric labels (using method abbreviations)
metric_labels <- c(
  "RBI" = "Richards-Baker Index (RBI)",
  "RCS" = "Recession Curve Slope (RCS)",
  "FDC" = "Flow Duration Curve (FDC)",
  "SD" = "Storage-Discharge (SD, mm)",
  "WB" = "Water Balance (WB, mm)",
  "CHS" = "Baseflow Fraction (CHS)"
)

# Publication theme
theme_pub <- theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    legend.position = "none"
  )

theme_set(theme_pub)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading storage metrics data...\n")

# Load master annual file
annual_file <- file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUTPUT_DIR, "HJA_StorageMetrics_Annual_All.csv")
}

annual_data <- read_csv(annual_file, show_col_types = FALSE) %>%
  # Rename old column names to new abbreviations
  rename_with(~ case_when(
    .x == "recession_curve_slope" ~ "RCS",
    .x == "fdc_slope" ~ "FDC",
    .x == "S_annual_mm" ~ "SD",
    .x == "mean_bf" ~ "CHS",
    .x == "DS_sum" ~ "WB",
    TRUE ~ .x
  )) %>%
  # Standardize site names
  mutate(
    site = case_when(
      site %in% c("GSWSMC", "GSMACK") ~ "Mack",
      TRUE ~ site
    )
  ) %>%
  filter(site %in% site_order) %>%
  mutate(site = factor(site, levels = site_order))

cat("Loaded", nrow(annual_data), "rows\n")

# =============================================================================
# 1. DYNAMIC STORAGE METRICS (RBI, RCS, FDC, SD)
# =============================================================================

cat("\nCreating Dynamic Storage faceted plot...\n")

dynamic_metrics <- c("RBI", "RCS", "FDC", "SD")

dynamic_long <- annual_data %>%
  select(site, year, all_of(dynamic_metrics)) %>%
  pivot_longer(
    cols = all_of(dynamic_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric, levels = dynamic_metrics))

# Calculate summary stats
dynamic_summary <- dynamic_long %>%
  group_by(site, metric) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Faceted mean ± SD plot
p_dynamic <- ggplot(dynamic_summary, aes(x = site, y = mean_val, color = site)) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.3, linewidth = 0.5
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2,
             labeller = labeller(metric = metric_labels)) +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(
    title = "Dynamic Storage Metrics",
    x = NULL,
    y = "Mean ± 1 SD"
  )

ggsave(
  file.path(output_dir, "Dynamic_Storage_Faceted.png"),
  p_dynamic, width = 10, height = 8, dpi = 300
)
ggsave(
  file.path(output_dir, "Dynamic_Storage_Faceted.pdf"),
  p_dynamic, width = 10, height = 8
)

cat("Saved: Dynamic_Storage_Faceted.png/pdf\n")

# =============================================================================
# 2. EXTENDED DYNAMIC STORAGE (WB only)
# =============================================================================

cat("\nCreating Extended Dynamic (WB) plot...\n")

wb_data <- annual_data %>%
  select(site, year, WB) %>%
  filter(!is.na(WB))

wb_summary <- wb_data %>%
  group_by(site) %>%
  summarise(
    mean_val = mean(WB, na.rm = TRUE),
    sd_val = sd(WB, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_wb <- ggplot(wb_summary, aes(x = site, y = mean_val, color = site)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.3, linewidth = 0.6
  ) +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(
    title = "Extended Dynamic Storage - Water Balance (WB)",
    subtitle = "Cumulative summer recession storage change (mm)",
    x = NULL,
    y = "Mean ± 1 SD (mm)"
  )

ggsave(
  file.path(output_dir, "Extended_Dynamic_WB.png"),
  p_wb, width = 8, height = 5, dpi = 300
)
ggsave(
  file.path(output_dir, "Extended_Dynamic_WB.pdf"),
  p_wb, width = 8, height = 5
)

cat("Saved: Extended_Dynamic_WB.png/pdf\n")

# =============================================================================
# 3. TIME SERIES - DYNAMIC STORAGE METRICS
# =============================================================================

cat("\nCreating Dynamic Storage time series...\n")

p_dynamic_ts <- ggplot(dynamic_long, aes(x = year, y = value, color = site)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  facet_grid(metric ~ site, scales = "free_y",
             labeller = labeller(metric = metric_labels)) +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(
    title = "Dynamic Storage Metrics - Time Series",
    x = "Water Year",
    y = "Value"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text.y = element_text(angle = 0, size = 8),
    strip.text.x = element_text(size = 8)
  )

ggsave(
  file.path(output_dir, "Dynamic_Storage_TimeSeries.png"),
  p_dynamic_ts, width = 14, height = 10, dpi = 300
)
ggsave(
  file.path(output_dir, "Dynamic_Storage_TimeSeries.pdf"),
  p_dynamic_ts, width = 14, height = 10
)

cat("Saved: Dynamic_Storage_TimeSeries.png/pdf\n")

# =============================================================================
# 4. TIME SERIES - EXTENDED DYNAMIC (WB)
# =============================================================================

cat("\nCreating WB time series...\n")

p_wb_ts <- ggplot(wb_data, aes(x = year, y = WB, color = site)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~site, ncol = 2) +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(
    title = "Extended Dynamic Storage (WB) - Time Series",
    x = "Water Year",
    y = "Water Balance (mm)"
  )

ggsave(
  file.path(output_dir, "Extended_Dynamic_WB_TimeSeries.png"),
  p_wb_ts, width = 10, height = 8, dpi = 300
)
ggsave(
  file.path(output_dir, "Extended_Dynamic_WB_TimeSeries.pdf"),
  p_wb_ts, width = 10, height = 8
)

cat("Saved: Extended_Dynamic_WB_TimeSeries.png/pdf\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== ORGANIZED PLOTS COMPLETE ===\n")
cat("Output directory:", output_dir, "\n")
cat("\nFigures created:\n")
cat("  - Dynamic_Storage_Faceted.png/pdf (RBI, RCS, FDC, SD)\n")
cat("  - Dynamic_Storage_TimeSeries.png/pdf\n")
cat("  - Extended_Dynamic_WB.png/pdf (WB)\n")
cat("  - Extended_Dynamic_WB_TimeSeries.png/pdf\n")
cat("\nMobile storage figures in Publication folder:\n")
cat("  - Mobile_Isotope_Metrics.png/pdf (DR, Fyw, MTT1, MTT2)\n")
cat("  - Mobile_CHS_Boxplot.png/pdf (CHS)\n")
