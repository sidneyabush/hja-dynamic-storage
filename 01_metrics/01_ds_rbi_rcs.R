# -----------------------------------------------------------------------------
# Dynamic Storage: RBI and Recession Curve Slope
# -----------------------------------------------------------------------------
# Purpose: Calculate Richards-Baker Flashiness Index (RBI) and recession curve
#          slope for each site and water year
#
# Methods:
#   - RBI: sum(|dQ|) / sum(Q) for each water year
#   - Recession slope: log-log linear regression of -dQ/dt vs Q
#
# Inputs: Daily discharge data (HF00402_v14.csv)
# Outputs: Annual RBI and recession slopes
# -----------------------------------------------------------------------------

library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(scales)
library(colorspace)
library(tidyr)

# Clear environment
rm(list = ls())

# Source configuration (paths, site definitions, water year range)
# Get script directory for relative sourcing
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

theme_set(theme_pub(base_size = 14))

# Use configuration values
output_dir <- OUT_MET_DYNAMIC_DIR
discharge_dir <- DISCHARGE_DIR
catchment_dir <- CATCHMENT_CHARACTERISTICS_DIR
sites_keep <- SITE_ORDER_HYDROMETRIC

# Read & prep data
da_df <- read_csv(resolve_drainage_area_file(), show_col_types = FALSE) %>%
  mutate(SITECODE = standardize_site_code(SITECODE))
discharge <- read_csv(file.path(discharge_dir, "HF00402_v14.csv"), show_col_types = FALSE) %>%
  # Standardize site codes (e.g., GSWSMC -> Mack, GSLOOK -> Look)
  mutate(SITECODE = standardize_site_code(SITECODE)) %>%
  # Filter to water years 1997-2020 and hydrometric sites
  filter(WATERYEAR >= WY_START, WATERYEAR <= WY_END, SITECODE %in% sites_keep) %>%
  left_join(da_df, by = "SITECODE") %>%
  filter(!is.na(DA_M2)) %>%
  mutate(
    Date = as.Date(DATE, "%m/%d/%Y"),
    Q = MEAN_Q * 0.02831683199881 # m³/s
  ) %>%
  arrange(SITECODE, Date)

# Need to update this function to take the log-log:
# lm_model <- lm(log(recession_slope) ~ log(Q), data = recession_data)

# 1) Recession slope
calc_recession <- function(df) {
  tmp <- df %>%
    mutate(
      dQ = Q - lag(Q),
      dQ_dt = dQ / as.numeric(Date - lag(Date)),
      change_ratio = Q / lag(Q)
    ) %>%
    filter(!is.na(dQ_dt), change_ratio >= 0.7, dQ < 0) %>%
    mutate(recession_slope = -dQ_dt)

  tibble(slope = coef(lm(log(recession_slope) ~ log(Q), data = tmp))[2])
}

# 2) RBI
calc_RBI <- function(df) {
  tmp <- df %>% mutate(dQ = Q - lag(Q)) %>% filter(!is.na(dQ))
  total_Q <- sum(df$Q, na.rm = TRUE)
  props <- abs(tmp$dQ) / total_Q

  tibble(RBI = sum(props, na.rm = TRUE))
}

# Compute annual metrics
annual_metrics <- discharge %>%
  group_by(SITECODE, WATERYEAR) %>%
  group_map(
    ~ bind_cols(
      tibble(site = .y$SITECODE, year = .y$WATERYEAR),
      calc_recession(.x),
      calc_RBI(.x)
    )
  ) %>%
  bind_rows() %>%
  mutate(site = factor(site, levels = sites_keep))

# Prepare full‐record recession data in mm/day
recession_clean <- discharge %>%
  mutate(
    dQ = Q - lag(Q),
    dQ_dt = dQ / as.numeric(Date - lag(Date)), # m³/s per day
    recession_slope = -dQ_dt, # m³/s per day
    # convert to mm/day: (m³/s) ÷ basin area (m²) × 86400 s/day × 1000 mm/m
    Q_mm_day = Q / DA_M2 * 86400 * 1000,
    slope_mm_day = recession_slope / DA_M2 * 86400 * 1000
  ) %>%
  filter(!is.na(Q_mm_day), !is.na(slope_mm_day)) %>%
  transmute(
    site = factor(SITECODE, levels = sites_keep),
    Q_mm_day,
    slope_mm_day
  )

# Annotation positions (force positive slopes)
slopes_df <- recession_clean %>%
  group_by(site) %>%
  summarise(
    slope_mm = abs(round(coef(lm(slope_mm_day ~ Q_mm_day))[2], 2)),
    Q_pos_mm = quantile(Q_mm_day, 0.9, na.rm = TRUE),
    slope_pos_mm = quantile(slope_mm_day, 0.85, na.rm = TRUE),
    .groups = "drop"
  )

# Color palette → lighten by 10%
palette10 <- c(
  "#AA4499",
  "#882255",
  "#CC6677",
  "#DDCC77",
  "#999933",
  "#117733",
  "#44AA99",
  "#88CCEE",
  "#6699CC",
  "#332288"
)
site_cols <- setNames(lighten(palette10, amount = 0.1), sites_keep)

# Save annual metrics for aggregation script
# RCS = Recession Curve Slope (method name)
annual_metrics %>%
  rename(RCS = slope) %>%
  select(site, year, RCS, RBI) %>%
  write_csv(file.path(output_dir, "rbi_rcs_annual.csv"))

# End of script
