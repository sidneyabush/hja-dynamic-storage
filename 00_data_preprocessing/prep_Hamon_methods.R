# calibrate hamon et to pt et methods and produce daily water-balance et inputs.
# inputs: out_met_support_dir/pt_et_methods_timeseries.csv
# outputs: out_met_support_dir/daily_water_balance_et_hamon_zhang_coeff_interp.csv

library(readr)
library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)
library(zoo)

source("config.R")

# --- directories ---
input_dir <- OUT_MET_SUPPORT_DIR
output_dir <- OUT_MET_SUPPORT_DIR
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. load pt et timeseries and compute uncalibrated hamon ---
pt_file <- file.path(input_dir, "PT_ET_methods_timeseries.csv")
if (!file.exists(pt_file)) {
  stop("Missing required PT ET timeseries file: ", pt_file)
}

pt <- read_csv(pt_file, show_col_types = FALSE) %>%
  mutate(DATE = as.Date(DATE))

calculate_et_hamon <- function(temp_c, date_vec, lat, coeff) {
  doy      <- yday(date_vec)
  delta    <- 0.4093 * sin((2*pi*doy/365) - 1.405)
  phi_rad  <- lat * pi/180
  cos_sun  <- pmax(-1, pmin(1, -tan(phi_rad)*tan(delta)))
  omega    <- acos(cos_sun)
  dayl     <- (24/pi)*omega
  es       <- 0.6108 * exp(17.27*temp_c/(temp_c+237.3))
  rho_v    <- 216.7 * es/(temp_c+273.16)
  et_raw   <- coeff * (dayl/12) * rho_v
  ifelse(is.na(temp_c), NA_real_, pmax(0, et_raw))
}

study_lat <- (44.20127400 + 44.28226000)/2

data <- pt %>%
  mutate(
    ET_Hamon_uncalibrated = calculate_et_hamon(T_C, DATE, study_lat, 0.1651),
    MONTH                 = month(DATE)
  )

calibration_start <- as.Date("2013-01-01")
calibration_end <- as.Date(sprintf("%d-09-30", WY_END))

# --- 2. derive monthly median coefficients (target: 2013 through wy_end) ---
cal_window_raw <- data %>%
  filter(DATE >= calibration_start, DATE <= calibration_end)

cal_window_complete <- cal_window_raw %>%
  filter(!is.na(ET_PT_zhang), !is.na(ET_Hamon_uncalibrated))

if (nrow(cal_window_complete) == 0) {
  stop(
    "No overlapping PT and uncalibrated Hamon ET values in calibration window: ",
    as.character(calibration_start),
    " to ",
    as.character(calibration_end)
  )
}

calibration_window_summary <- cal_window_complete %>%
  mutate(
    waterYear = get_water_year(DATE)
  ) %>%
  summarise(
    calibration_window_start = min(DATE, na.rm = TRUE),
    calibration_window_end = max(DATE, na.rm = TRUE),
    calibration_wy_start = min(waterYear, na.rm = TRUE),
    calibration_wy_end = max(waterYear, na.rm = TRUE),
    n_days_complete = n()
  )

write_csv(
  calibration_window_summary,
  file.path(output_dir, "et_calibration_window_summary.csv")
)

monthly_coef_zhang <- cal_window_raw %>%
  filter(!is.na(ET_PT_zhang), !is.na(ET_Hamon_uncalibrated)) %>%
  group_by(SITECODE, MONTH) %>%
  summarise(coef_zhang_month = median(ET_PT_zhang / ET_Hamon_uncalibrated),
            .groups="drop")

monthly_coef_szilagyi <- cal_window_raw %>%
  filter(!is.na(ET_PT_szilagyi), !is.na(ET_Hamon_uncalibrated)) %>%
  group_by(SITECODE, MONTH) %>%
  summarise(coef_szilagyi_month = median(ET_PT_szilagyi / ET_Hamon_uncalibrated),
            .groups="drop")

# --- 3. fallback: site‐mean & global median coefficients ---
site_mean_zhang    <- monthly_coef_zhang    %>% group_by(SITECODE) %>% summarise(coef_zhang_site    = mean(coef_zhang_month),    .groups="drop")
site_mean_szilagyi <- monthly_coef_szilagyi %>% group_by(SITECODE) %>% summarise(coef_szilagyi_site = mean(coef_szilagyi_month), .groups="drop")

global_coef_zhang    <- mean(monthly_coef_zhang$coef_zhang_month,    na.rm=TRUE)
global_coef_szilagyi <- mean(monthly_coef_szilagyi$coef_szilagyi_month, na.rm=TRUE)

# --- 4. merge in coefficients & compute monthly‐calibrated hamon ---
data <- data %>%
  left_join(monthly_coef_zhang,    by = c("SITECODE","MONTH")) %>%
  left_join(site_mean_zhang,       by = "SITECODE") %>%
  left_join(monthly_coef_szilagyi, by = c("SITECODE","MONTH")) %>%
  left_join(site_mean_szilagyi,    by = "SITECODE") %>%
  mutate(
    coef_zhang_final    = coalesce(coef_zhang_month,    coef_zhang_site,    global_coef_zhang),
    coef_szilagyi_final = coalesce(coef_szilagyi_month, coef_szilagyi_site, global_coef_szilagyi),
    ET_Hamon_pt_zhang_monthly   = ifelse(is.na(ET_Hamon_uncalibrated),
                                         NA_real_,
                                         pmax(0, ET_Hamon_uncalibrated * coef_zhang_final)),
    ET_Hamon_pt_szilagyi_monthly = ifelse(is.na(ET_Hamon_uncalibrated),
                                          NA_real_,
                                          pmax(0, ET_Hamon_uncalibrated * coef_szilagyi_final))
  )

# --- 5. re‐define calibration window on the _new_ columns ---
cal_window <- data %>%
  filter(DATE >= calibration_start, DATE <= calibration_end)

# --- 6. fit per‐site regressions (require ≥10 pts) ---
reg_sites_zhang    <- cal_window %>%
  filter(!is.na(ET_PT_zhang), !is.na(ET_Hamon_pt_zhang_monthly)) %>%
  pull(SITECODE) %>% unique()

reg_sites_szilagyi <- cal_window %>%
  filter(!is.na(ET_PT_szilagyi), !is.na(ET_Hamon_pt_szilagyi_monthly)) %>%
  pull(SITECODE) %>% unique()

regressions_zhang <- map(reg_sites_zhang, function(site) {
  df <- cal_window %>% filter(SITECODE == site)
  if(nrow(df) >= 10) lm(ET_Hamon_pt_zhang_monthly ~ ET_PT_zhang, data = df) else NULL
}) %>% set_names(reg_sites_zhang)

regressions_szilagyi <- map(reg_sites_szilagyi, function(site) {
  df <- cal_window %>% filter(SITECODE == site)
  if(nrow(df) >= 10) lm(ET_Hamon_pt_szilagyi_monthly ~ ET_PT_szilagyi, data = df) else NULL
}) %>% set_names(reg_sites_szilagyi)

# --- 7. predict via regression where monthly‐calibrated is na but pt exists ---
pred_zhang <- map_dfr(reg_sites_zhang, function(site) {
  m <- regressions_zhang[[site]]
  df <- data %>%
    filter(SITECODE == site,
           is.na(ET_Hamon_pt_zhang_monthly),
           !is.na(ET_PT_zhang))
  if(!is.null(m) && nrow(df)>0) {
    df %>% transmute(SITECODE, DATE,
                     ET_Hamon_pt_zhang_pred = pmax(0, predict(m, newdata = df)))
  } else {
    tibble(SITECODE=character(), DATE=as.Date(character()), ET_Hamon_pt_zhang_pred=double())
  }
})

pred_szilagyi <- map_dfr(reg_sites_szilagyi, function(site) {
  m <- regressions_szilagyi[[site]]
  df <- data %>%
    filter(SITECODE == site,
           is.na(ET_Hamon_pt_szilagyi_monthly),
           !is.na(ET_PT_szilagyi))
  if(!is.null(m) && nrow(df)>0) {
    df %>% transmute(SITECODE, DATE,
                     ET_Hamon_pt_szilagyi_pred = pmax(0, predict(m, newdata = df)))
  } else {
    tibble(SITECODE=character(), DATE=as.Date(character()), ET_Hamon_pt_szilagyi_pred=double())
  }
})

data <- data %>%
  left_join(pred_zhang,    by = c("SITECODE","DATE")) %>%
  left_join(pred_szilagyi, by = c("SITECODE","DATE"))

# --- 8. create initial interpolated‐by‐coefficient‐and‐regression series ---
data <- data %>%
  mutate(
    ET_Hamon_pt_zhang_interpolated   = coalesce(ET_Hamon_pt_zhang_monthly,
                                                ET_Hamon_pt_zhang_pred),
    ET_Hamon_pt_szilagyi_interpolated = coalesce(ET_Hamon_pt_szilagyi_monthly,
                                                 ET_Hamon_pt_szilagyi_pred)
  )

# --- 9. fill any remaining na by time‐series interpolation within each site ---
data <- data %>%
  arrange(SITECODE, DATE) %>%
  group_by(SITECODE) %>%
  mutate(
    # 9a) linear interpolate interior gaps
    ET_Hamon_pt_zhang_interp_full   = na.approx(ET_Hamon_pt_zhang_interpolated,
                                                x = DATE, na.rm = FALSE),
    ET_Hamon_pt_szilagyi_interp_full = na.approx(ET_Hamon_pt_szilagyi_interpolated,
                                                 x = DATE, na.rm = FALSE),
    
    # 9b) carry‐forward to fill any remaining leading nas
    ET_Hamon_pt_zhang_interp_full   = na.locf(ET_Hamon_pt_zhang_interp_full,
                                              na.rm = FALSE),
    ET_Hamon_pt_szilagyi_interp_full = na.locf(ET_Hamon_pt_szilagyi_interp_full,
                                               na.rm = FALSE),
    
    # 9c) carry‐backward to fill any remaining trailing nas
    ET_Hamon_pt_zhang_interp_full   = na.locf(ET_Hamon_pt_zhang_interp_full,
                                              fromLast = TRUE),
    ET_Hamon_pt_szilagyi_interp_full = na.locf(ET_Hamon_pt_szilagyi_interp_full,
                                               fromLast = TRUE)
  ) %>%
  ungroup()

# --- 10. drop all the temporary coefficient & prediction columns ---
data <- data %>%
  select(-contains("coef"),
         -contains("alpha"),
         -contains("interpolated"),
         -contains("pred"),
         -MONTH
         )

# --- 11. export full methods table (for diagnostics) ---
write_csv(
  data,
  file.path(output_dir, "daily_MET_water_balance_all_ET_methods_1997_present.csv")
)

# --- 12. export water-balance subset with all et methods ---
wb_data <- data %>%
  select(-any_of(c("RH_d_pct", "NR_Wm2_d", "VPD_kPa")))

write_csv(
  wb_data,
  file.path(output_dir, "daily_water_balance_all_ET_methods_1997_present.csv")
)

# --- 13. export default et series used by wb dynamic storage ---
wb_data_final_et <- data %>%
  select(
    -any_of(c("RH_d_pct", "NR_Wm2_d", "VPD_kPa")),
    -contains("szilagyi"),
    -contains("monthly"),
    -contains("uncalibrated"),
    -any_of("ET_PT_zhang")
  ) %>%
  rename(ET_mm_d = ET_Hamon_pt_zhang_interp_full)

write_csv(
  wb_data_final_et,
  file.path(output_dir, "daily_water_balance_et_hamon_zhang_coeff_interp.csv")
)
