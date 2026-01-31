# =============================================================================
# Meteorological Context Figure
# =============================================================================
# Purpose: Create figure showing monthly climatology (median with IQR; Jan–Dec)
#          and annual anomalies (departures from WY 1997–2020 mean; dashed line = 0)
#          for precipitation, temperature, and SWE at the CENMET station.
#
# Output:
#   - Fig_Met_Context.png/pdf: 2x3 panel figure
#     Top row: Monthly climatology (median with IQR)
#     Bottom row: Annual anomalies relative to WY mean
#   - Panel labels: a)–f)
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

script_dir <- tryCatch(
  {
    dirname(sys.frame(1)$ofile)
  },
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      getwd()
    }
  }
)

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
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

wy_start <- 1997
wy_end <- 2020

# =============================================================================
# MUTED, COHESIVE PALETTE (consistent across rows)
# =============================================================================
# Precipitation — darker, calm blue
col_precip <- "#4A6FA5"

# Temperature — subdued orange (not yellow, not red)
col_temp <- "#D08C2F"

# SWE — darker glacier teal
col_swe <- "#4F9A8A"

# Negative anomaly tints (same hues, much lighter)
col_precip_neg <- "#C9D7EB"
col_temp_neg <- "#F0D9B5"
col_swe_neg <- "#CFE5E1"

# Top row (variable identity by lightness)
col_precip <- "#4D4D4D" # dark gray
col_temp <- "#7A7A7A" # medium gray
col_swe <- "#B0B0B0" # light gray

# Bottom row anomalies (same variable identity)
# Positive = same gray as top row
# Negative = lighter tint of same gray
col_precip_neg <- "#BFBFBF"
col_temp_neg <- "#D9D9D9"
col_swe_neg <- "#E6E6E6"


# =============================================================================
# THEME (no bold; centered titles; small panel tags)
# =============================================================================

theme_pub <- theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    plot.title = element_text(face = "plain", size = 10, hjust = 0.5),
    axis.title = element_text(face = "plain", size = 9),
    axis.text = element_text(size = 8)
  )

theme_set(theme_pub)

# =============================================================================
# LOAD DATA
# =============================================================================

temp_data <- read_csv(
  file.path(met_dir, "Temperature_filtered_post1997.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    DATE = as.Date(DATE),
    month = month(DATE),
    year = year(DATE),
    water_year = ifelse(month >= 10, year + 1, year)
  ) %>%
  filter(water_year >= wy_start & water_year <= wy_end) %>%
  select(month, water_year, T_C = CENMET)

precip_data <- read_csv(
  file.path(met_dir, "Precipitation_filtered_post1997.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    DATE = as.Date(DATE),
    month = month(DATE),
    year = year(DATE),
    water_year = ifelse(month >= 10, year + 1, year)
  ) %>%
  filter(water_year >= wy_start & water_year <= wy_end) %>%
  select(month, water_year, P_mm = CENMET)

swe_data <- read_csv(
  file.path(met_dir, "SWE_original_&_filled_1997_2023_v5.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    DATE = as.Date(DATE, format = "%m/%d/%Y"),
    month = month(DATE),
    year = year(DATE),
    water_year = ifelse(month >= 10, year + 1, year)
  ) %>%
  filter(water_year >= wy_start & water_year <= wy_end) %>%
  select(month, water_year, SWE_mm = CENMET)

# =============================================================================
# MONTHLY CLIMATOLOGY (median + IQR), ORDERED BY CALENDAR MONTH
# =============================================================================

temp_monthly <- temp_data %>%
  group_by(month) %>%
  summarise(
    med = median(T_C, na.rm = TRUE),
    q25 = quantile(T_C, 0.25, na.rm = TRUE),
    q75 = quantile(T_C, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

precip_monthly <- precip_data %>%
  group_by(month) %>%
  summarise(
    med = median(P_mm, na.rm = TRUE),
    q25 = quantile(P_mm, 0.25, na.rm = TRUE),
    q75 = quantile(P_mm, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

swe_monthly <- swe_data %>%
  group_by(month) %>%
  summarise(
    med = median(SWE_mm, na.rm = TRUE),
    q25 = quantile(SWE_mm, 0.25, na.rm = TRUE),
    q75 = quantile(SWE_mm, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

month_labels <- c(
  "Jan",
  "Feb",
  "Mar",
  "Apr",
  "May",
  "Jun",
  "Jul",
  "Aug",
  "Sep",
  "Oct",
  "Nov",
  "Dec"
)

temp_monthly$month <- factor(
  temp_monthly$month,
  levels = 1:12,
  labels = month_labels
)
precip_monthly$month <- factor(
  precip_monthly$month,
  levels = 1:12,
  labels = month_labels
)
swe_monthly$month <- factor(
  swe_monthly$month,
  levels = 1:12,
  labels = month_labels
)

# =============================================================================
# ANNUAL ANOMALIES (WY aggregates)
# =============================================================================

temp_annual <- temp_data %>%
  group_by(water_year) %>%
  summarise(val = mean(T_C, na.rm = TRUE), .groups = "drop") %>%
  mutate(anom = val - mean(val, na.rm = TRUE))

precip_annual <- precip_data %>%
  group_by(water_year) %>%
  summarise(val = sum(P_mm, na.rm = TRUE), .groups = "drop") %>%
  mutate(anom = val - mean(val, na.rm = TRUE))

swe_annual <- swe_data %>%
  group_by(water_year) %>%
  summarise(val = max(SWE_mm, na.rm = TRUE), .groups = "drop") %>%
  mutate(anom = val - mean(val, na.rm = TRUE))

# =============================================================================
# PLOTS
# =============================================================================

# Top row: Monthly climatology (median + IQR); ONLY these have titles
p1 <- ggplot(precip_monthly, aes(month, med)) +
  geom_col(fill = col_precip, color = "black", linewidth = 0.15, alpha = 0.8) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, linewidth = 0.3) +
  labs(title = "Precipitation", y = "Precipitation (mm)", x = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

p2 <- ggplot(temp_monthly, aes(month, med)) +
  geom_col(fill = col_temp, color = "black", linewidth = 0.15, alpha = 0.8) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, linewidth = 0.3) +
  labs(title = "Temperature", y = "Temperature (°C)", x = NULL)

p3 <- ggplot(swe_monthly, aes(month, med)) +
  geom_col(fill = col_swe, color = "black", linewidth = 0.15, alpha = 0.8) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, linewidth = 0.3) +
  labs(title = "SWE", y = "SWE (mm)", x = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Bottom row: Anomalies (no titles; no legends)
p4 <- ggplot(precip_annual, aes(water_year, anom)) +
  geom_col(
    aes(fill = anom > 0),
    color = "black",
    linewidth = 0.15,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("TRUE" = col_precip, "FALSE" = col_precip_neg)) +
  guides(fill = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = NULL, x = "Water Year", y = "Anomaly (mm)") +
  scale_x_continuous(breaks = seq(1998, 2020, by = 4))

p5 <- ggplot(temp_annual, aes(water_year, anom)) +
  geom_col(
    aes(fill = anom > 0),
    color = "black",
    linewidth = 0.15,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("TRUE" = col_temp, "FALSE" = col_temp_neg)) +
  guides(fill = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = NULL, x = "Water Year", y = "Anomaly (°C)") +
  scale_x_continuous(breaks = seq(1998, 2020, by = 4))

p6 <- ggplot(swe_annual, aes(water_year, anom)) +
  geom_col(
    aes(fill = anom > 0),
    color = "black",
    linewidth = 0.15,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("TRUE" = col_swe, "FALSE" = col_swe_neg)) +
  guides(fill = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = NULL, x = "Water Year", y = "Anomaly (mm)") +
  scale_x_continuous(breaks = seq(1998, 2020, by = 4))

# =============================================================================
# COMBINE AND SAVE (no overall title; no caption; panel labels a)–f))
# =============================================================================

fig_combined <- (p1 | p2 | p3) /
  (p4 | p5 | p6) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(size = 7, face = "plain"), # smaller, journal-like
      plot.tag.position = c(0.02, 0.985)
    )
  )

ggsave(
  file.path(output_dir, "Fig_Met_Context.png"),
  fig_combined,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(output_dir, "Fig_Met_Context.pdf"),
  fig_combined,
  width = 12,
  height = 7
)
