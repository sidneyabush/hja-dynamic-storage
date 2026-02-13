# -----------------------------------------------------------------------------
# Meteorological Context Figure
# -----------------------------------------------------------------------------
# Purpose: Create figure showing monthly climatology (median with IQR; Jan–Dec)
#          and annual anomalies (departures from WY 1997–2020 mean; dashed line = 0)
#          for precipitation, temperature, and SWE at the CENMET station.
#
# Output:
#   - met_context.png/pdf: 2x3 panel figure
#     Top row: Monthly climatology (median with IQR)
#     Bottom row: Annual anomalies relative to WY mean
#   - Panel labels: a)–f)
#
# Study Period: Water Years 1997-2020
#
# Author: Sidney Bush
# Date: 2026-01-30
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(patchwork)

rm(list = ls())

# -----------------------------------------------------------------------------
# SOURCE CONFIGURATION
# -----------------------------------------------------------------------------

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
  stop("config.R not found.")
}

# -----------------------------------------------------------------------------
# SETUP
# -----------------------------------------------------------------------------

met_dir <- file.path(BASE_DATA_DIR, "all_hydromet")
output_dir <- file.path(FIGURES_DIR, "main")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

wy_start <- 1997
wy_end <- 2020
wy_breaks <- seq(2000, wy_end, by = 5)

# -----------------------------------------------------------------------------
# GREYSCALE PALETTE (print-safe, manuscript-ready)
# -----------------------------------------------------------------------------
col_precip <- "#4D4D4D" # dark gray
col_temp <- "#7A7A7A" # medium gray
col_swe <- "#B0B0B0" # light gray


# -----------------------------------------------------------------------------
# THEME
# -----------------------------------------------------------------------------

theme_pub_base <- theme_pub() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    plot.caption = element_blank(),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE)
  )

theme_set(theme_pub_base)

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

temp_data <- read_csv(
  file.path(met_dir, "Temperature_original_&_filled_1979_2023_v2.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    DATE = as.Date(DATE, tryFormats = c("%m/%d/%Y", "%m/%d/%y", "%Y-%m-%d")),
    month = month(DATE),
    year = year(DATE),
    water_year = ifelse(month >= 10, year + 1, year),
    T_C = coalesce(CENMET_inter, CENMET)
  ) %>%
  filter(water_year >= wy_start & water_year <= wy_end) %>%
  select(month, water_year, T_C)

precip_data <- read_csv(
  file.path(met_dir, "Precipitation_original_&_filled_1979_2023.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    DATE = as.Date(DATE, tryFormats = c("%m/%d/%Y", "%m/%d/%y", "%Y-%m-%d")),
    month = month(DATE),
    year = year(DATE),
    water_year = ifelse(month >= 10, year + 1, year),
    P_mm = coalesce(CENMET_inter, CENMET)
  ) %>%
  filter(water_year >= wy_start & water_year <= wy_end) %>%
  select(month, water_year, P_mm)

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

# -----------------------------------------------------------------------------
# MONTHLY CLIMATOLOGY (median + IQR), ORDERED BY CALENDAR MONTH
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# ANNUAL ANOMALIES (WY aggregates)
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# PLOTS
# -----------------------------------------------------------------------------

# Top row: Monthly climatology (median + IQR)
p1 <- ggplot(precip_monthly, aes(month, med)) +
  geom_col(fill = col_precip, color = "black", linewidth = 0.15, alpha = 0.8) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, linewidth = 0.3) +
  labs(y = "Precipitation (mm)", x = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(temp_monthly, aes(month, med)) +
  geom_col(fill = col_temp, color = "black", linewidth = 0.15, alpha = 0.8) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, linewidth = 0.3) +
  labs(y = "Temperature (°C)", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(swe_monthly, aes(month, med)) +
  geom_col(fill = col_swe, color = "black", linewidth = 0.15, alpha = 0.8) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, linewidth = 0.3) +
  labs(y = "SWE (mm)", x = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Bottom row: Anomalies (no titles; no legends)
p4 <- ggplot(precip_annual, aes(water_year, anom)) +
  geom_col(fill = col_precip, color = "black", linewidth = 0.15, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Water Year", y = "Anomaly (mm)") +
  scale_x_continuous(breaks = wy_breaks) +
  coord_cartesian(xlim = c(wy_start - 0.5, wy_end + 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p5 <- ggplot(temp_annual, aes(water_year, anom)) +
  geom_col(fill = col_temp, color = "black", linewidth = 0.15, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Water Year", y = "Anomaly (°C)") +
  scale_x_continuous(breaks = wy_breaks) +
  coord_cartesian(xlim = c(wy_start - 0.5, wy_end + 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p6 <- ggplot(swe_annual, aes(water_year, anom)) +
  geom_col(fill = col_swe, color = "black", linewidth = 0.15, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Water Year", y = "Anomaly (mm)") +
  scale_x_continuous(breaks = wy_breaks) +
  coord_cartesian(xlim = c(wy_start - 0.5, wy_end + 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -----------------------------------------------------------------------------
# COMBINE AND SAVE (no overall title; no caption; panel labels a)–f))
# -----------------------------------------------------------------------------

fig_combined <- (p1 | p2 | p3) /
  (p4 | p5 | p6) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(size = FIG_STRIP_TEXT_SIZE, face = "plain"),
      plot.tag.position = c(0.02, 0.985)
    )
  )

ggsave(
  file.path(output_dir, "met_context.png"),
  fig_combined,
  width = 13 * FIG_WIDTH_SCALE,
  height = 8.5 * FIG_HEIGHT_SCALE,
  dpi = 300
)

ggsave(
  file.path(output_dir, "met_context.pdf"),
  fig_combined,
  width = 13 * FIG_WIDTH_SCALE,
  height = 8.5 * FIG_HEIGHT_SCALE
)
