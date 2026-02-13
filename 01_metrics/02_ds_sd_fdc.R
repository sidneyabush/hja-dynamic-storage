# -----------------------------------------------------------------------------
# Dynamic Storage: Storage-Discharge Method + FDC Slope (SD metric)
# -----------------------------------------------------------------------------
# Purpose: Calculate dynamic storage using storage-discharge relationship and
#          flow duration curve (FDC) slopes for each site and water year
#
# Methods:
#   - Fit recession law: dQ/dt = k*Q^p using rain-free falling limbs
#   - Storage depth: S = integral of Q/dQ/dt from Qmin to Qmax
#   - FDC slope: slope of log10(Q) vs exceedance probability (5th-95th)
#
# Inputs: Daily water balance data with P, Q, ET
# Outputs: Annual storage depths and FDC slopes
# -----------------------------------------------------------------------------

# Load libraries & clear environment
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
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

# Directories & input
output_dir <- OUT_MET_DYNAMIC_DIR
input_file <- resolve_water_balance_daily_file()
stopifnot(file.exists(input_file))

# Read data, filter sites, parse date & compute water-year
df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(
    # robustly parse DATE into a Date object
    date = parse_date_time(
      DATE,
      orders = c("Ymd", "Y-m-d", "mdy", "dmy")
    ) %>% as.Date(),

    # water-year: Oct–Dec → next calendar year, Jan–Sep → same year
    wateryear = if_else(
      month(date) >= 10,
      year(date) + 1,
      year(date)
    )
  ) %>%
  rename(
    site = SITECODE,
    P    = P_mm_d,
    Q    = Q_mm_d,
    ET   = ET_mm_d
  ) %>%
  # Standardize site codes (e.g., GSLOOK_FULL -> Look, GSWSMC -> Mack)
  mutate(site = standardize_site_code(site)) %>%
  # Filter to hydrometric sites and water year range
  filter(site %in% SITE_ORDER_HYDROMETRIC,
         wateryear >= WY_START, wateryear <= WY_END) %>%
  # Set factor levels for consistent ordering
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site, date)


# Helper functions
compute_fdc <- function(Q) {
  Qpos <- Q[Q > 0]
  Qs   <- sort(Qpos, decreasing = TRUE)
  n    <- length(Qs)
  r    <- seq_len(n)
  exc  <- (r - 0.44) / (n + 0.12) * 100
  tibble(exceedance = exc, Q = Qs)
}

get_Q <- function(fdc, P_exc) {
  if (nrow(fdc) < 2) return(NA_real_)
  approx(fdc$exceedance, fdc$Q, xout = P_exc)$y
}

deltaS <- function(Qu, Ql, k, p) {
  (Qu^(2 - p) - Ql^(2 - p)) / (k * (2 - p))
}

# Site-water-year analysis function
analyze <- function(df_sub) {
  # Flow-duration thresholds
  fdc  <- compute_fdc(df_sub$Q)
  Q99  <- get_Q(fdc, 99)
  Q50  <- get_Q(fdc, 50)
  Q01  <- get_Q(fdc, 1)
  if (any(is.na(c(Q99, Q50, Q01)))) {
    return(tibble(
      k = NA_real_, p = NA_real_,
      Q_max = NA_real_, Q_min = NA_real_,
      Q99 = Q99, Q50 = Q50, Q01 = Q01,
      S_annual_mm   = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm  = NA_real_
    ))
  }
  
  # Recession data (rain-free falling limbs)
  tmp <- df_sub %>%
    arrange(date) %>%
    mutate(
      dt      = as.numeric(difftime(date, lag(date), "days")),
      dQ      = (lag(Q) - Q) / dt,
      is_rain = P > 0
    )
  rec <- tmp %>% filter(!is_rain, !is.na(dQ), dQ > 0, Q > 0)
  if (nrow(rec) < 10) {
    return(tibble(
      k = NA_real_, p = NA_real_,
      Q_max = NA_real_, Q_min = NA_real_,
      Q99 = Q99, Q50 = Q50, Q01 = Q01,
      S_annual_mm   = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm  = NA_real_
    ))
  }
  
  # Fit recession law
  fit   <- lm(log(dQ) ~ log(Q), data = rec)
  p_est <- coef(fit)["log(Q)"]
  k_est <- exp(coef(fit)["(Intercept)"])
  
  # Compute annual-mean extremes by water-year
  yearly <- df_sub %>%
    group_by(wateryear) %>%
    summarize(
      yr_max = if (all(is.na(Q))) NA_real_ else max(Q, na.rm = TRUE),
      yr_min = if (all(Q <= 0 | is.na(Q))) NA_real_ else min(Q[Q > 0], na.rm = TRUE),
      .groups = "drop"
    )
  Qmax_avg <- mean(yearly$yr_max, na.rm = TRUE)
  Qmin_avg <- mean(yearly$yr_min, na.rm = TRUE)
  
  # Compute dynamic storage depths
  tibble(
    k             = k_est,
    p             = p_est,
    Q_max         = Qmax_avg,
    Q_min         = Qmin_avg,
    Q99           = Q99,
    Q50           = Q50,
    Q01           = Q01,
    S_annual_mm   = deltaS(Qmax_avg, Qmin_avg, k_est, p_est),
    S_high_med_mm = deltaS(Q01, Q50, k_est, p_est),
    S_med_low_mm  = deltaS(Q50, Q99, k_est, p_est)
  )
}

# Compute overall results (across all water-years)
overall <- df %>%
  group_by(site) %>%
  group_map(~ analyze(.x) %>% mutate(site = .y$site), .keep = TRUE) %>%
  bind_rows()

# Compute per-water-year results
annual <- df %>%
  group_by(site, wateryear) %>%
  group_map(~ analyze(.x) %>% 
              mutate(site = .y$site, wateryear = .y$wateryear),
            .keep = TRUE) %>%
  bind_rows()

# Add catchment areas & convert mm → m³
areas_ha <- tibble(
  SITECODE = c("GSWS01","GSWS02","GSWS03","GSWS06","GSWS07",
               "GSWS08","GSWS09","GSWS10", "Mack", "Look"),
  area_ha  = c(96, 60, 101, 13.0, 15.4, 21.4, 8.5, 10.2, 580, 6242)
)
add_vol <- function(df_in) {
  df_in %>%
    left_join(areas_ha, by = c("site" = "SITECODE")) %>%
    mutate(
      area_m2        = area_ha * 10000,
      S_annual_m3    = (S_annual_mm   / 1000) * area_m2,
      S_high_med_m3  = (S_high_med_mm / 1000) * area_m2,
      S_med_low_m3   = (S_med_low_mm  / 1000) * area_m2
    )
}

overall_vol <- add_vol(overall)
annual_vol  <- add_vol(annual)

# Calculate FDC slope (Q5-Q95), Q5norm, and CV_Q5norm
# FDC slope: slope of log10(Q) ~ exceedance probability for 5th-95th percentile
# Q5norm: 5th percentile discharge during Aug-Oct low-flow period (mm/day)
# CV_Q5norm: coefficient of variation of annual Q5norm values

# Calculate FDC slopes (overall per site, 5th-95th percentile)
fdc_slopes <- df %>%
  filter(Q > 0) %>%  # Remove zero/negative values before log transform
  group_by(site) %>%
  arrange(desc(Q)) %>%
  mutate(
    rank = row_number(),
    n = n(),
    exceedance_prob = 100 * rank / (n + 1)
  ) %>%
  filter(exceedance_prob >= 5, exceedance_prob <= 95) %>%
  summarise(
    fdc_slope = coef(lm(log10(Q) ~ exceedance_prob))[2],
    .groups = "drop"
  )

# Legacy-compatible FDC exports:
# - FDC_slopes_overall.csv (one slope per site)
# - FDC_slopes_WY.csv (one slope per site-water year)
fdc_curves_overall <- df %>%
  filter(Q > 0) %>%
  group_by(site) %>%
  arrange(desc(Q), .by_group = TRUE) %>%
  mutate(
    Rank = row_number(),
    N = n(),
    ExceedanceProbability = 100 * Rank / (N + 1)
  ) %>%
  filter(ExceedanceProbability >= 5, ExceedanceProbability <= 95) %>%
  summarise(
    Slope = coef(lm(log10(Q) ~ ExceedanceProbability))[2],
    .groups = "drop"
  )

fdc_curves_WY <- df %>%
  filter(Q > 0) %>%
  group_by(site, wateryear) %>%
  arrange(desc(Q), .by_group = TRUE) %>%
  mutate(
    Rank = row_number(),
    N = n(),
    ExceedanceProbability = 100 * Rank / (N + 1)
  ) %>%
  filter(ExceedanceProbability >= 5, ExceedanceProbability <= 95) %>%
  summarise(
    Slope = coef(lm(log10(Q) ~ ExceedanceProbability))[2],
    .groups = "drop"
  ) %>%
  rename(WaterYear = wateryear)

# Calculate Q5norm: 5th percentile during Aug-Oct (months 8-10)
Q5_annual <- df %>%
  filter(month(date) >= 8, month(date) <= 10) %>%
  group_by(site, wateryear) %>%
  summarise(
    Q5norm = quantile(Q, 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate CV_Q5norm: coefficient of variation across years
CV_Q5norm <- Q5_annual %>%
  group_by(site) %>%
  summarise(
    CV_Q5norm = sd(Q5norm, na.rm = TRUE) / mean(Q5norm, na.rm = TRUE),
    .groups = "drop"
  )

# Merge FDC metrics into annual data
annual <- annual %>%
  left_join(fdc_slopes, by = "site") %>%
  left_join(Q5_annual, by = c("site", "wateryear")) %>%
  left_join(CV_Q5norm, by = "site")

# Save outputs (with wateryear in filenames for clarity)
write.csv(overall_vol,
          file = file.path(output_dir, "storage_overall_per_site.csv"),
          row.names = FALSE)
write.csv(annual,
          file = file.path(output_dir,
                           "storage_discharge_method_annual_mm_metrics_per_site_wateryear.csv"),
          row.names = FALSE)
write.csv(annual_vol,
          file = file.path(output_dir,
                           "storage_discharge_method_annual_vol_per_site_wateryear.csv"),
          row.names = FALSE)
write.csv(fdc_curves_overall,
          file = file.path(output_dir, "fdc_slopes_overall.csv"),
          row.names = FALSE)
write.csv(fdc_curves_WY,
          file = file.path(output_dir, "fdc_slopes_wy.csv"),
          row.names = FALSE)

# Save annual FDC metrics for aggregation script
# SD = Storage-Discharge method, FDC = Flow Duration Curve method
write.csv(
  annual %>%
    rename(SD = S_annual_mm, FDC = fdc_slope) %>%
    select(site, year = wateryear, SD, FDC, Q99, Q50, Q01, Q5norm, CV_Q5norm),
  file.path(output_dir, "storage_discharge_fdc_annual.csv"),
  row.names = FALSE
)
