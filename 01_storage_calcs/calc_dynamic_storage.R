# calculate hydrometric dynamic storage metrics (RBI/RCS + SD/FDC) by site and water year

# inputs:
# discharge_dir/HF00402_v14.csv
# catchment_characteristics_dir/drainage_area.csv
# out_met_support_dir/daily_water_balance_et_hamon_zhang_coeff_interp.csv

# outputs:
# outputs/metrics/dynamic/rbi_rcs_annual.csv
# outputs/metrics/dynamic/fdc_slopes_overall.csv
# outputs/metrics/dynamic/fdc_slopes_wy.csv
# outputs/metrics/dynamic/storage_discharge_fdc_annual.csv
# outputs/metrics/extended_dynamic/ds_depletion_annual.csv

# author: Sidney Bush
# date: 2026-02-14

librarian::shelf(dplyr, lubridate, readr, tidyr, tibble, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

output_dir <- OUT_MET_DYNAMIC_DIR
output_dir_ed <- OUT_MET_EXTENDED_DIR
discharge_dir <- DISCHARGE_DIR
sites_keep <- SITE_ORDER_HYDROMETRIC

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_ed, recursive = TRUE, showWarnings = FALSE)

# calculate RBI and recession slope (RCS) from Q

# pair discharge with drainage area for each gaged catchment
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
  # keep gradual flow declines and fit the recession exponent
  tmp <- df %>%
    mutate(
      dQ = Q - lag(Q),
      dQ_dt = dQ / as.numeric(Date - lag(Date)),
      change_ratio = Q / lag(Q)
    ) %>%
    filter(!is.na(dQ_dt), change_ratio >= 0.7, dQ < 0) %>%
    mutate(recession_slope = -dQ_dt)

  if (nrow(tmp) < 10) {
    return(tibble(slope = NA_real_))
  }

  tibble(slope = coef(lm(log(recession_slope) ~ log(Q), data = tmp))[2])
}

calc_RBI <- function(df) {
  # sum day to day flow changes relative to total annual flow
  tmp <- df %>% mutate(dQ = Q - lag(Q)) %>% filter(!is.na(dQ))
  total_Q <- sum(df$Q, na.rm = TRUE)
  props <- abs(tmp$dQ) / total_Q
  tibble(RBI = sum(props, na.rm = TRUE))
}

# calculate annual RBI and RCS for each catchment
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

# calculate storage discharge and FDC metrics from water balance data

wb_daily_file <- resolve_water_balance_daily_file()

# use the daily water balance table for FDC, storage discharge, and WB metrics
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
  # empirical flow duration curve for positive daily discharge
  Qpos <- Q[Q > 0]
  Qs <- sort(Qpos, decreasing = TRUE)
  n <- length(Qs)
  r <- seq_len(n)
  exc <- (r - 0.44) / (n + 0.12) * 100
  tibble(exceedance = exc, Q = Qs)
}

get_Q <- function(fdc, P_exc) {
  # interpolate Q99, Q50, or Q01 from the empirical FDC
  if (nrow(fdc) < 2) {
    return(NA_real_)
  }
  approx(fdc$exceedance, fdc$Q, xout = P_exc)$y
}

deltaS <- function(Qu, Ql, k, p) {
  # storage change implied by the fitted recession relationship
  (Qu^(2 - p) - Ql^(2 - p)) / (k * (2 - p))
}

analyze <- function(df_sub) {
  # calculate FDC quantiles and recession storage for one site or site year
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

  rec <- tmp %>%
    filter(
      !is_rain,
      !is.na(dQ),
      dQ > 0,
      Q > 0
    )
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

# annual storage values are used in Figure 2, ANOVA/Tukey tests, and ecological response models
annual <- wb_df %>%
  group_by(site, wateryear) %>%
  group_map(
    ~ analyze(.x) %>% mutate(site = .y$site, wateryear = .y$wateryear),
    .keep = TRUE
  ) %>%
  bind_rows()

# full period FDC slopes are used in master_site, catchment models, and Figure 7
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

# annual FDC slopes are used for Figure 2 and ANOVA/Tukey output
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

# full period FDC slope used in master_site
fdc_slopes_site <- fdc_slopes %>%
  transmute(
    site = site,
    fdc_slope = fdc_slope
  )

# Q5 used for eco response
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
  left_join(fdc_slopes_site, by = "site") %>%
  left_join(q5_annual, by = c("site", "wateryear")) %>%
  left_join(cv_q5norm, by = "site")

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

# annual dynamic storage output used by the aggregation script
write.csv(
  annual %>%
    rename(SD = S_annual_mm, FDC = fdc_slope) %>%
    select(site, year = wateryear, SD, FDC, Q99, Q50, Q01, Q5norm, CV_Q5norm),
  file.path(output_dir, "storage_discharge_fdc_annual.csv"),
  row.names = FALSE
)

# calculate annual water balance depletion

# keep only complete water years before calculating depletion
wb_daily <- read.csv(wb_daily_file, stringsAsFactors = FALSE) %>%
  mutate(
    SITECODE = standardize_site_code(SITECODE),
    DATE = parse_date_time(
      DATE,
      orders = c("Ymd", "Y-m-d", "mdy", "m/d/y", "dmy")
    ) %>% as.Date(),
    waterYear = get_water_year(DATE)
  ) %>%
  filter(
    SITECODE %in% SITE_ORDER_HYDROMETRIC,
    waterYear >= WY_START,
    waterYear <= WY_END,
    is.finite(P_mm_d),
    is.finite(Q_mm_d),
    is.finite(ET_mm_d)
  ) %>%
  arrange(SITECODE, waterYear, DATE)

goodyears_wb <- wb_daily %>%
  group_by(SITECODE, waterYear) %>%
  summarise(num_days = n_distinct(DATE), .groups = "drop") %>%
  filter(num_days >= 365)

wb_daily <- wb_daily %>%
  semi_join(goodyears_wb, by = c("SITECODE", "waterYear"))

calc_wb_depletion <- function(df_sub) {
  # annual WB is the largest within year drawdown from cumulative P minus Q minus ET
  df_sub <- df_sub %>% arrange(DATE)
  ds_daily <- df_sub$P_mm_d - df_sub$Q_mm_d - df_sub$ET_mm_d
  wb_cum <- cumsum(ds_daily)
  drawdown <- cummax(wb_cum) - wb_cum

  if (!any(is.finite(drawdown))) {
    return(tibble(
      peak_storage_date = as.Date(NA),
      depletion_date = as.Date(NA),
      WB = NA_real_
    ))
  }

  depletion_idx <- which.max(drawdown)[1]
  peak_idx <- which.max(wb_cum[seq_len(depletion_idx)])[1]

  tibble(
    peak_storage_date = df_sub$DATE[peak_idx],
    depletion_date = df_sub$DATE[depletion_idx],
    WB = as.numeric(drawdown[depletion_idx])
  )
}

wb_summary <- wb_daily %>%
  group_by(SITECODE, waterYear) %>%
  group_modify(~ calc_wb_depletion(.x)) %>%
  ungroup()

write.csv(
  wb_summary %>%
    select(SITECODE, waterYear, WB),
  file.path(output_dir_ed, "ds_depletion_annual.csv"),
  row.names = FALSE
)
