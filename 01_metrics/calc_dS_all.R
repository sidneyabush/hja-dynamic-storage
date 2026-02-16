# Calculate hydrometric dynamic-storage metrics (RBI/RCS + SD/FDC) by site and water year.
# Inputs: DISCHARGE_DIR/HF00402_v14.csv; CATCHMENT_CHARACTERISTICS_DIR/drainage_area.csv; OUT_MET_SUPPORT_DIR/daily_water_balance_et_hamon_zhang_coeff_interp.csv.
# Author: Sidney Bush
# Date: 2026-02-14

library(dplyr)
library(lubridate)
library(readr)
library(tidyr)

rm(list = ls())

# Load project config
source("config.R")

output_dir <- OUT_MET_DYNAMIC_DIR
discharge_dir <- DISCHARGE_DIR
sites_keep <- SITE_ORDER_HYDROMETRIC

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- PART 1: RBI and recession slope (RCS) from discharge ----

da_df <- read_csv(resolve_drainage_area_file(), show_col_types = FALSE) %>%
  mutate(SITECODE = standardize_site_code(SITECODE))

discharge <- read_csv(
  file.path(discharge_dir, "HF00402_v14.csv"),
  show_col_types = FALSE
) %>%
  mutate(SITECODE = standardize_site_code(SITECODE)) %>%
  filter(WATERYEAR >= WY_START, WATERYEAR <= WY_END, SITECODE %in% sites_keep) %>%
  left_join(da_df, by = "SITECODE") %>%
  filter(!is.na(DA_M2)) %>%
  mutate(
    Date = as.Date(DATE, "%m/%d/%Y"),
    Q = MEAN_Q * 0.02831683199881
  ) %>%
  arrange(SITECODE, Date)

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

calc_RBI <- function(df) {
  tmp <- df %>% mutate(dQ = Q - lag(Q)) %>% filter(!is.na(dQ))
  total_Q <- sum(df$Q, na.rm = TRUE)
  props <- abs(tmp$dQ) / total_Q
  tibble(RBI = sum(props, na.rm = TRUE))
}

rbi_rcs_annual <- discharge %>%
  group_by(SITECODE, WATERYEAR) %>%
  group_map(
    ~ bind_cols(
      tibble(site = .y$SITECODE, year = .y$WATERYEAR),
      calc_recession(.x),
      calc_RBI(.x)
    )
  ) %>%
  bind_rows() %>%
  rename(RCS = slope) %>%
  mutate(site = factor(site, levels = sites_keep)) %>%
  select(site, year, RCS, RBI)

write_csv(rbi_rcs_annual, file.path(output_dir, "rbi_rcs_annual.csv"))

# ---- PART 2: Storage-discharge and FDC metrics from water balance ----

wb_daily_file <- resolve_water_balance_daily_file()
stopifnot(file.exists(wb_daily_file))

wb_df <- read.csv(wb_daily_file, stringsAsFactors = FALSE) %>%
  mutate(
    date = parse_date_time(DATE, orders = c("Ymd", "Y-m-d", "mdy", "dmy")) %>%
      as.Date(),
    wateryear = if_else(month(date) >= 10, year(date) + 1, year(date))
  ) %>%
  rename(
    site = SITECODE,
    P = P_mm_d,
    Q = Q_mm_d,
    ET = ET_mm_d
  ) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(
    site %in% SITE_ORDER_HYDROMETRIC,
    wateryear >= WY_START,
    wateryear <= WY_END
  ) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site, date)

compute_fdc <- function(Q) {
  Qpos <- Q[Q > 0]
  Qs <- sort(Qpos, decreasing = TRUE)
  n <- length(Qs)
  r <- seq_len(n)
  exc <- (r - 0.44) / (n + 0.12) * 100
  tibble(exceedance = exc, Q = Qs)
}

get_Q <- function(fdc, P_exc) {
  if (nrow(fdc) < 2) {
    return(NA_real_)
  }
  approx(fdc$exceedance, fdc$Q, xout = P_exc)$y
}

deltaS <- function(Qu, Ql, k, p) {
  (Qu^(2 - p) - Ql^(2 - p)) / (k * (2 - p))
}

analyze <- function(df_sub) {
  fdc <- compute_fdc(df_sub$Q)
  Q99 <- get_Q(fdc, 99)
  Q50 <- get_Q(fdc, 50)
  Q01 <- get_Q(fdc, 1)

  if (any(is.na(c(Q99, Q50, Q01)))) {
    return(tibble(
      k = NA_real_,
      p = NA_real_,
      Q_max = NA_real_,
      Q_min = NA_real_,
      Q99 = Q99,
      Q50 = Q50,
      Q01 = Q01,
      S_annual_mm = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm = NA_real_
    ))
  }

  tmp <- df_sub %>%
    arrange(date) %>%
    mutate(
      dt = as.numeric(difftime(date, lag(date), "days")),
      dQ = (lag(Q) - Q) / dt,
      is_rain = P > 0
    )

  rec <- tmp %>% filter(!is_rain, !is.na(dQ), dQ > 0, Q > 0)
  if (nrow(rec) < 10) {
    return(tibble(
      k = NA_real_,
      p = NA_real_,
      Q_max = NA_real_,
      Q_min = NA_real_,
      Q99 = Q99,
      Q50 = Q50,
      Q01 = Q01,
      S_annual_mm = NA_real_,
      S_high_med_mm = NA_real_,
      S_med_low_mm = NA_real_
    ))
  }

  fit <- lm(log(dQ) ~ log(Q), data = rec)
  p_est <- coef(fit)["log(Q)"]
  k_est <- exp(coef(fit)["(Intercept)"])

  yearly <- df_sub %>%
    group_by(wateryear) %>%
    summarize(
      yr_max = if (all(is.na(Q))) NA_real_ else max(Q, na.rm = TRUE),
      yr_min = if (all(Q <= 0 | is.na(Q))) {
        NA_real_
      } else {
        min(Q[Q > 0], na.rm = TRUE)
      },
      .groups = "drop"
    )

  Qmax_avg <- mean(yearly$yr_max, na.rm = TRUE)
  Qmin_avg <- mean(yearly$yr_min, na.rm = TRUE)

  tibble(
    k = k_est,
    p = p_est,
    Q_max = Qmax_avg,
    Q_min = Qmin_avg,
    Q99 = Q99,
    Q50 = Q50,
    Q01 = Q01,
    S_annual_mm = deltaS(Qmax_avg, Qmin_avg, k_est, p_est),
    S_high_med_mm = deltaS(Q01, Q50, k_est, p_est),
    S_med_low_mm = deltaS(Q50, Q99, k_est, p_est)
  )
}

overall <- wb_df %>%
  group_by(site) %>%
  group_map(~ analyze(.x) %>% mutate(site = .y$site), .keep = TRUE) %>%
  bind_rows()

annual <- wb_df %>%
  group_by(site, wateryear) %>%
  group_map(
    ~ analyze(.x) %>% mutate(site = .y$site, wateryear = .y$wateryear),
    .keep = TRUE
  ) %>%
  bind_rows()

areas_ha <- tibble(
  SITECODE = c(
    "GSWS01", "GSWS02", "GSWS03", "GSWS06", "GSWS07",
    "GSWS08", "GSWS09", "GSWS10", "Mack", "Look"
  ),
  area_ha = c(96, 60, 101, 13.0, 15.4, 21.4, 8.5, 10.2, 580, 6242)
)

add_vol <- function(df_in) {
  df_in %>%
    left_join(areas_ha, by = c("site" = "SITECODE")) %>%
    mutate(
      area_m2 = area_ha * 10000,
      S_annual_m3 = (S_annual_mm / 1000) * area_m2,
      S_high_med_m3 = (S_high_med_mm / 1000) * area_m2,
      S_med_low_m3 = (S_med_low_mm / 1000) * area_m2
    )
}

overall_vol <- add_vol(overall)
annual_vol <- add_vol(annual)

fdc_slopes <- wb_df %>%
  filter(Q > 0) %>%
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

fdc_curves_overall <- wb_df %>%
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

fdc_curves_wy <- wb_df %>%
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

q5_annual <- wb_df %>%
  filter(month(date) >= 8, month(date) <= 10) %>%
  group_by(site, wateryear) %>%
  summarise(
    Q5norm = quantile(Q, 0.05, na.rm = TRUE),
    .groups = "drop"
  )

cv_q5norm <- q5_annual %>%
  group_by(site) %>%
  summarise(
    CV_Q5norm = sd(Q5norm, na.rm = TRUE) / mean(Q5norm, na.rm = TRUE),
    .groups = "drop"
  )

annual <- annual %>%
  left_join(fdc_slopes, by = "site") %>%
  left_join(q5_annual, by = c("site", "wateryear")) %>%
  left_join(cv_q5norm, by = "site")

write.csv(
  overall_vol,
  file = file.path(output_dir, "storage_overall_per_site.csv"),
  row.names = FALSE
)
write.csv(
  annual,
  file = file.path(
    output_dir,
    "storage_discharge_method_annual_mm_metrics_per_site_wateryear.csv"
  ),
  row.names = FALSE
)
write.csv(
  annual_vol,
  file = file.path(
    output_dir,
    "storage_discharge_method_annual_vol_per_site_wateryear.csv"
  ),
  row.names = FALSE
)
write.csv(
  fdc_curves_overall,
  file = file.path(output_dir, "fdc_slopes_overall.csv"),
  row.names = FALSE
)
write.csv(
  fdc_curves_wy,
  file = file.path(output_dir, "fdc_slopes_wy.csv"),
  row.names = FALSE
)

# Canonical annual dynamic-storage output consumed by aggregation script
write.csv(
  annual %>%
    rename(SD = S_annual_mm, FDC = fdc_slope) %>%
    select(site, year = wateryear, SD, FDC, Q99, Q50, Q01, Q5norm, CV_Q5norm),
  file.path(output_dir, "storage_discharge_fdc_annual.csv"),
  row.names = FALSE
)
