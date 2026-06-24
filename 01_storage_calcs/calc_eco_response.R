# calculate thermal, low flow, and wet season precipitation metrics

# inputs:
# stream_temp_dir/HT00451_v10.txt
# out_met_support_dir/catchments_met_q.csv

# outputs:
# outputs/metrics/eco/stream_thermal_lowflow_metrics_annual.csv

# author: Sidney Bush
# date: 2026-01-23

librarian::shelf(dplyr, lubridate, readr, zoo, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

base_dir      <- BASE_DATA_DIR
output_dir    <- OUT_MET_ECO_DIR
temp_dir      <- STREAM_TEMP_DIR

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sites_keep <- SITE_ORDER_HYDROMETRIC
target_years <- WY_START:WY_END

stop_if_repeated_rows <- function(df, columns, df_name) {
  dupes <- df %>%
    count(across(all_of(columns)), name = "n") %>%
    filter(n > 1)

  # stop if repeated rows would duplicate outputs
  if (nrow(dupes) > 0) {
    stop(
      paste0(
        "Repeated rows in ", df_name, " for columns (",
        paste(columns, collapse = ", "), ")."
      )
    )
  }
}

temp_file <- file.path(temp_dir, "HT00451_v10.txt")
temp_raw <- read_csv(temp_file, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(DATE_TIME),
    site = standardize_site_code(SITECODE)
  ) %>%
  select(site, date, temp_mean_C = WATERTEMP_MEAN) %>%
  filter(site %in% sites_keep)

temp_daily <- temp_raw %>%
  group_by(site, date) %>%
  summarise(
    temp_mean_C = mean(temp_mean_C, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(temp_mean_C))

# stop early if the stream temperature file does not include enough study sites
if (nrow(temp_daily) == 0 || n_distinct(temp_daily$site) < 3) {
  available_streams <- temp_raw %>%
    distinct(site) %>%
    pull(site) %>%
    unique()

  stop(
    paste0(
      "Not enough stream temperature data for the analysis sites. ",
      "Found stream IDs: ", paste(available_streams, collapse = ", "), ". ",
      "Sites needed: ", paste(sites_keep, collapse = ", "), "."
    )
  )
}

# WS09 has no stream temperature record
# keep WS09 in the annual master tables but leave the temperature outputs as NA

met_support_file <- file.path(OUT_MET_SUPPORT_DIR, "catchments_met_q.csv")
met_support <- read_csv(met_support_file, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%m/%d/%y")),
    site = standardize_site_code(SITECODE),
    P_mm_d = suppressWarnings(as.numeric(P_mm_d)),
    Q_mm_d = suppressWarnings(as.numeric(Q_mm_d)),
    WATERYEAR = get_water_year(date)
  ) %>%
  filter(site %in% sites_keep, !is.na(date))
stop_if_repeated_rows(met_support, c("site", "date"), "met_support")

# use Q_mm_d from catchments_met_q so low flow uses the same discharge depth as storage metrics
discharge <- met_support %>%
  select(site, date, Q_mm_d, WATERYEAR) %>%
  arrange(site, date)

met_daily <- met_support %>%
  select(site, date, P_mm_d)

# wet season precipitation from November through May
precip_nov_may <- met_daily %>%
  mutate(
    month_num = month(date),
    season_year = if_else(month_num %in% c(11L, 12L), year(date) + 1L, year(date))
  ) %>%
  filter(month_num %in% c(11L, 12L, 1L, 2L, 3L, 4L, 5L), season_year %in% target_years) %>%
  group_by(site, year = season_year) %>%
  summarise(
    n_days = sum(is.finite(P_mm_d)),
    Pws = if_else(n_days > 0, sum(P_mm_d, na.rm = TRUE), as.numeric(NA)),
    .groups = "drop"
  ) %>%
  mutate(
    precip_nov_may_mm = Pws
  ) %>%
  select(site, year, Pws, precip_nov_may_mm)
stop_if_repeated_rows(precip_nov_may, c("site", "year"), "precip_nov_may")

calc_7day_rolling <- function(df, value_col) {
  df %>%
    arrange(date) %>%
    mutate(
      rolling_7d = zoo::rollmean(.data[[value_col]], k = 7,
                                  fill = NA, align = "center")
    )
}

temp_rolling <- temp_daily %>%
  group_by(site) %>%
  calc_7day_rolling("temp_mean_C") %>%
  ungroup() %>%
  mutate(year = get_water_year(date)) %>%
  rename(temp_7d_avg_C = rolling_7d)

discharge_rolling <- discharge %>%
  group_by(site) %>%
  calc_7day_rolling("Q_mm_d") %>%
  ungroup() %>%
  mutate(year = get_water_year(date)) %>%
  rename(Q_7d_avg_mm_d = rolling_7d)

t_7dmax <- temp_rolling %>%
  filter(year %in% target_years) %>%
  group_by(site, year) %>%
  filter(!is.na(temp_7d_avg_C)) %>%
  slice_max(temp_7d_avg_C, n = 1, with_ties = FALSE) %>%
  select(
    site,
    year,
    date_t_7dmax = date,
    T_7DMax = temp_7d_avg_C
  ) %>%
  ungroup()
stop_if_repeated_rows(t_7dmax, c("site", "year"), "t_7dmax")

# Q5 of the 7-day average discharge in each water year
q_7q5 <- discharge_rolling %>%
  filter(year %in% target_years) %>%
  group_by(site, year) %>%
  summarise(
    Q_7Q5 = quantile(Q_7d_avg_mm_d, probs = 0.05, na.rm = TRUE),
    .groups = "drop"
  )
stop_if_repeated_rows(q_7q5, c("site", "year"), "q_7q5")

master_metrics <- t_7dmax %>%
  full_join(q_7q5, by = c("site", "year")) %>%
  left_join(precip_nov_may, by = c("site", "year")) %>%
  arrange(site, year)
stop_if_repeated_rows(master_metrics, c("site", "year"), "master_metrics")

output_file <- file.path(output_dir, "stream_thermal_lowflow_metrics_annual.csv")
write_csv(master_metrics, output_file)
