# calibrate Hamon ET to PT ET methods and produce daily water balance ET inputs

# inputs:
# out_met_support_dir/PT_ET_methods_timeseries.csv

# outputs:
# out_met_support_dir/daily_water_balance_et_hamon_zhang_coeff_interp.csv

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(readr, dplyr, lubridate, purrr, zoo, cran_repo = "https://cloud.r-project.org")

source("config.R")

# use the shared support folder for PT inputs and Hamon output
input_dir <- OUT_MET_SUPPORT_DIR
output_dir <- OUT_MET_SUPPORT_DIR
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# load PT ET timeseries and compute uncalibrated Hamon
pt_file <- file.path(input_dir, "PT_ET_methods_timeseries.csv")
pt <- read_csv(pt_file, show_col_types = FALSE) %>%
  mutate(DATE = as.Date(DATE))

calculate_et_hamon <- function(temp_c, date_vec, lat, coeff) {
  doy <- yday(date_vec)
  delta <- 0.4093 * sin((2 * pi * doy / 365) - 1.405)
  phi_rad <- lat * pi / 180
  cos_sun <- pmax(-1, pmin(1, -tan(phi_rad) * tan(delta)))
  omega <- acos(cos_sun)
  dayl <- (24 / pi) * omega
  es <- 0.6108 * exp(17.27 * temp_c / (temp_c + 237.3))
  rho_v <- 216.7 * es / (temp_c + 273.16)
  et_raw <- coeff * (dayl / 12) * rho_v
  ifelse(is.na(temp_c), NA_real_, pmax(0, et_raw))
}

study_lat <- (44.20127400 + 44.28226000) / 2

data <- pt %>%
  mutate(
    ET_Hamon_uncalibrated = calculate_et_hamon(T_C, DATE, study_lat, 0.1651),
    MONTH = month(DATE)
  )

calibration_start <- as.Date("2013-01-01")
calibration_end <- as.Date(sprintf("%d-09-30", WY_END))

# derive monthly median coefficients for the calibration period
cal_window_raw <- data %>%
  filter(DATE >= calibration_start, DATE <= calibration_end)

cal_window_complete <- cal_window_raw %>%
  filter(!is.na(ET_PT_zhang), !is.na(ET_Hamon_uncalibrated))

# stop if there is no overlap for the ET calibration
if (nrow(cal_window_complete) == 0) {
  stop(
    "No overlapping PT and uncalibrated Hamon ET values in calibration window: ",
    as.character(calibration_start),
    " to ",
    as.character(calibration_end)
  )
}

monthly_coef_zhang <- cal_window_raw %>%
  filter(!is.na(ET_PT_zhang), !is.na(ET_Hamon_uncalibrated)) %>%
  group_by(SITECODE, MONTH) %>%
  summarise(
    coef_zhang_month = median(ET_PT_zhang / ET_Hamon_uncalibrated),
    .groups = "drop"
  )

# fallback to site mean and global median coefficients
site_mean_zhang <- monthly_coef_zhang %>%
  group_by(SITECODE) %>%
  summarise(coef_zhang_site = mean(coef_zhang_month), .groups = "drop")

global_coef_zhang <- mean(monthly_coef_zhang$coef_zhang_month, na.rm = TRUE)

# merge coefficients and compute monthly calibrated Hamon
data <- data %>%
  left_join(monthly_coef_zhang, by = c("SITECODE", "MONTH")) %>%
  left_join(site_mean_zhang, by = "SITECODE") %>%
  mutate(
    coef_zhang_final = coalesce(
      coef_zhang_month,
      coef_zhang_site,
      global_coef_zhang
    ),
    ET_Hamon_pt_zhang_monthly = ifelse(
      is.na(ET_Hamon_uncalibrated),
      NA_real_,
      pmax(0, ET_Hamon_uncalibrated * coef_zhang_final)
    )
  )

# define calibration window using the new columns
cal_window <- data %>%
  filter(DATE >= calibration_start, DATE <= calibration_end)

# fit per site regressions that require at least 10 points
reg_sites_zhang <- cal_window %>%
  filter(!is.na(ET_PT_zhang), !is.na(ET_Hamon_pt_zhang_monthly)) %>%
  pull(SITECODE) %>%
  unique()

regressions_zhang <- map(reg_sites_zhang, function(site) {
  df <- cal_window %>% filter(SITECODE == site)
  if (nrow(df) >= 10) {
    lm(ET_Hamon_pt_zhang_monthly ~ ET_PT_zhang, data = df)
  } else {
    NULL
  }
}) %>%
  set_names(reg_sites_zhang)

# predict from regression where monthly calibrated ET is missing but PT exists
pred_zhang <- map_dfr(reg_sites_zhang, function(site) {
  m <- regressions_zhang[[site]]
  df <- data %>%
    filter(
      SITECODE == site,
      is.na(ET_Hamon_pt_zhang_monthly),
      !is.na(ET_PT_zhang)
    )
  if (!is.null(m) && nrow(df) > 0) {
    df %>%
      transmute(
        SITECODE,
        DATE,
        ET_Hamon_pt_zhang_pred = pmax(0, predict(m, newdata = df))
      )
  } else {
    tibble(
      SITECODE = character(),
      DATE = as.Date(character()),
      ET_Hamon_pt_zhang_pred = double()
    )
  }
})

data <- data %>%
  left_join(pred_zhang, by = c("SITECODE", "DATE"))

# create initial series from coefficient and regression estimates
data <- data %>%
  mutate(
    ET_Hamon_pt_zhang_interpolated = coalesce(
      ET_Hamon_pt_zhang_monthly,
      ET_Hamon_pt_zhang_pred
    )
  )

# fill remaining missing values by interpolation within each site
data <- data %>%
  arrange(SITECODE, DATE) %>%
  group_by(SITECODE) %>%
  mutate(
    # linear interpolate interior gaps
    ET_Hamon_pt_zhang_interp_full = na.approx(
      ET_Hamon_pt_zhang_interpolated,
      x = DATE,
      na.rm = FALSE
    ),

    # carry forward to fill remaining leading missing values
    ET_Hamon_pt_zhang_interp_full = na.locf(
      ET_Hamon_pt_zhang_interp_full,
      na.rm = FALSE
    ),

    # carry backward to fill remaining trailing missing values
    ET_Hamon_pt_zhang_interp_full = na.locf(
      ET_Hamon_pt_zhang_interp_full,
      fromLast = TRUE
    )
  ) %>%
  ungroup()

# drop temporary coefficient and prediction columns
data <- data %>%
  select(
    -contains("coef"),
    -contains("alpha"),
    -contains("interpolated"),
    -contains("pred"),
    -MONTH
  )

# export ET series used by WB dynamic storage
wb_data_final_et <- data %>%
  select(
    -any_of(c("RH_d_pct", "NR_Wm2_d", "VPD_kPa")),
    -contains("monthly"),
    -contains("uncalibrated"),
    -any_of("ET_PT_zhang")
  ) %>%
  rename(ET_mm_d = ET_Hamon_pt_zhang_interp_full)

write_csv(
  wb_data_final_et,
  file.path(output_dir, "daily_water_balance_et_hamon_zhang_coeff_interp.csv")
)
