# =============================================================================
# Meteorological Context Figure
# =============================================================================
# Purpose: Create figure showing mean monthly climate and annual anomalies
#          for precipitation, temperature, and SWE at CENMET station
#
# Output:
#   - Fig_Met_Context.png/pdf: 2x3 panel figure
#     Top row: Mean monthly P, T, SWE
#     Bottom row: Annual anomalies relative to study period mean
#
# Study Period: Water Years 1997-2020
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(patchwork)

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

met_dir <- file.path(BASE_DATA_DIR, "all_hydromet")
output_dir <- file.path(FIGURES_DIR, "Publication")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Study period
wy_start <- 1997
wy_end <- 2020

# Theme for publication
theme_pub <- theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

theme_set(theme_pub)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading meteorological data...\n")

# Temperature
temp_data <- read_csv(file.path(met_dir, "Temperature_filtered_post1997.csv"),
                      show_col_types = FALSE) %>%
  mutate(
    DATE = as.Date(DATE),
    year = year(DATE),
    month = month(DATE),
    water_year = ifelse(month >= 10, year + 1, year)
  ) %>%
  filter(water_year >= wy_start & water_year <= wy_end) %>%
  select(DATE, year, month, water_year, T_C = CENMET)

# Precipitation
precip_data <- read_csv(file.path(met_dir, "Precipitation_filtered_post1997.csv"),
                        show_col_types = FALSE) %>%
  mutate(
    DATE = as.Date(DATE),
    year = year(DATE),
    month = month(DATE),
    water_year = ifelse(month >= 10, year + 1, year)
  ) %>%
  filter(water_year >= wy_start & water_year <= wy_end) %>%
  select(DATE, year, month, water_year, P_mm = CENMET)

# SWE
swe_data <- read_csv(file.path(met_dir, "SWE_original_&_filled_1997_2023_v5.csv"),
                     show_col_types = FALSE) %>%
  mutate(
    DATE = as.Date(DATE, format = "%m/%d/%Y"),
    year = year(DATE),
    month = month(DATE),
    water_year = ifelse(month >= 10, year + 1, year)
  ) %>%
  filter(water_year >= wy_start & water_year <= wy_end) %>%
  select(DATE, year, month, water_year, SWE_mm = CENMET)

cat("  Temperature records:", nrow(temp_data), "\n")
cat("  Precipitation records:", nrow(precip_data), "\n")
cat("  SWE records:", nrow(swe_data), "\n")

# =============================================================================
# CALCULATE MONTHLY MEANS (CLIMATOLOGY)
# =============================================================================

cat("Calculating monthly climatology...\n")

# Monthly means across all years
temp_monthly <- temp_data %>%
  group_by(month) %>%
  summarise(
    mean_T = mean(T_C, na.rm = TRUE),
    sd_T = sd(T_C, na.rm = TRUE),
    .groups = "drop"
  )

precip_monthly <- precip_data %>%
  group_by(month) %>%
  summarise(
    mean_P = mean(P_mm, na.rm = TRUE),
    sd_P = sd(P_mm, na.rm = TRUE),
    .groups = "drop"
  )

# For SWE, use monthly mean (typically reported as point-in-time)
swe_monthly <- swe_data %>%
  group_by(month) %>%
  summarise(
    mean_SWE = mean(SWE_mm, na.rm = TRUE),
    sd_SWE = sd(SWE_mm, na.rm = TRUE),
    .groups = "drop"
  )

# Reorder months for water year (Oct-Sep)
month_order <- c(10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9)
month_labels <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar",
                  "Apr", "May", "Jun", "Jul", "Aug", "Sep")

temp_monthly <- temp_monthly %>%
  mutate(month_wy = factor(month, levels = month_order, labels = month_labels))

precip_monthly <- precip_monthly %>%
  mutate(month_wy = factor(month, levels = month_order, labels = month_labels))

swe_monthly <- swe_monthly %>%
  mutate(month_wy = factor(month, levels = month_order, labels = month_labels))

# =============================================================================
# CALCULATE ANNUAL ANOMALIES
# =============================================================================

cat("Calculating annual anomalies...\n")

# Annual means
temp_annual <- temp_data %>%
  group_by(water_year) %>%
  summarise(annual_T = mean(T_C, na.rm = TRUE), .groups = "drop")

precip_annual <- precip_data %>%
  group_by(water_year) %>%
  summarise(annual_P = sum(P_mm, na.rm = TRUE), .groups = "drop")

# For SWE, use peak SWE (April 1 or max)
swe_annual <- swe_data %>%
  group_by(water_year) %>%
  summarise(annual_SWE = max(SWE_mm, na.rm = TRUE), .groups = "drop")

# Calculate long-term means
mean_T <- mean(temp_annual$annual_T, na.rm = TRUE)
mean_P <- mean(precip_annual$annual_P, na.rm = TRUE)
mean_SWE <- mean(swe_annual$annual_SWE, na.rm = TRUE)

# Calculate anomalies
temp_annual <- temp_annual %>%
  mutate(anomaly_T = annual_T - mean_T)

precip_annual <- precip_annual %>%
  mutate(anomaly_P = annual_P - mean_P)

swe_annual <- swe_annual %>%
  mutate(anomaly_SWE = annual_SWE - mean_SWE)

cat("  Mean annual T:", round(mean_T, 1), "C\n")
cat("  Mean annual P:", round(mean_P, 0), "mm\n")
cat("  Mean peak SWE:", round(mean_SWE, 0), "mm\n")

# =============================================================================
# CREATE PLOTS
# =============================================================================

cat("Creating plots...\n")

# --- Top Row: Monthly Climatology ---

# Precipitation
p1 <- ggplot(precip_monthly, aes(x = month_wy, y = mean_P)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_P - sd_P, ymax = mean_P + sd_P),
                width = 0.3, linewidth = 0.3) +
  labs(title = "Mean Monthly Precipitation",
       x = NULL, y = "Precipitation (mm)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Temperature
p2 <- ggplot(temp_monthly, aes(x = month_wy, y = mean_T)) +
  geom_col(fill = "gray70", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_T - sd_T, ymax = mean_T + sd_T),
                width = 0.3, linewidth = 0.3) +
  labs(title = "Mean Monthly Temperature",
       x = NULL, y = "Temperature (\u00B0C)")

# SWE
p3 <- ggplot(swe_monthly, aes(x = month_wy, y = mean_SWE)) +
  geom_col(fill = "gray70", alpha = 0.8) +
  geom_errorbar(aes(ymin = pmax(0, mean_SWE - sd_SWE), ymax = mean_SWE + sd_SWE),
                width = 0.3, linewidth = 0.3) +
  labs(title = "Mean Monthly SWE",
       x = NULL, y = "SWE (mm)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# --- Bottom Row: Annual Anomalies ---

# Precipitation anomaly
p4 <- ggplot(precip_annual, aes(x = water_year, y = anomaly_P)) +
  geom_col(aes(fill = anomaly_P > 0), show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "coral")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Annual Precipitation Anomaly",
       x = "Water Year", y = "Anomaly (mm)") +
  scale_x_continuous(breaks = seq(1998, 2020, by = 4))

# Temperature anomaly
p5 <- ggplot(temp_annual, aes(x = water_year, y = anomaly_T)) +
  geom_col(aes(fill = anomaly_T > 0), show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Annual Temperature Anomaly",
       x = "Water Year", y = "Anomaly (\u00B0C)") +
  scale_x_continuous(breaks = seq(1998, 2020, by = 4))

# SWE anomaly
p6 <- ggplot(swe_annual, aes(x = water_year, y = anomaly_SWE)) +
  geom_col(aes(fill = anomaly_SWE > 0), show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "dodgerblue3", "FALSE" = "coral")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Annual Peak SWE Anomaly",
       x = "Water Year", y = "Anomaly (mm)") +
  scale_x_continuous(breaks = seq(1998, 2020, by = 4))

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

# Combine into 2x3 grid
fig_combined <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(
    title = paste0("Meteorological Context: CENMET Station (WY ", wy_start, "-", wy_end, ")"),
    caption = paste0("Top: Mean monthly values (\u00B1 1 SD). Bottom: Annual anomalies relative to ",
                     wy_start, "-", wy_end, " mean (dashed line)."),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  )

# Save
ggsave(file.path(output_dir, "Fig_Met_Context.png"),
       fig_combined, width = 12, height = 7, dpi = 300)

ggsave(file.path(output_dir, "Fig_Met_Context.pdf"),
       fig_combined, width = 12, height = 7)

cat("\n=== FIGURE COMPLETE ===\n")
cat("Saved to:", file.path(output_dir, "Fig_Met_Context.png"), "\n")
cat("Study period: WY", wy_start, "-", wy_end, "\n")
