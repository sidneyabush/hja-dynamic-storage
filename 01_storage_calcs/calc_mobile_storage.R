# calculate mobile-storage metrics (chs + isotope metrics) by site and water year.
# inputs: discharge_dir/hf00402_v14.csv; ec_dir/cf01201_v4.txt; isotope_dir/mtt_fyw.csv; isotope_dir/dampingratios_2025-07-07.csv.
# author: keira johnson (original chs), sidney bush
# date: 2026-02-14

library(dplyr)
library(readr)
library(lubridate)

rm(list = ls())

# load project config
source("config.R")

output_dir <- OUT_MET_MOBILE_DIR
discharge_dir <- DISCHARGE_DIR
ec_dir <- EC_DIR
isotope_dir <- ISOTOPE_DIR

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- part 1: chs from specific conductance and discharge ----

discharge <- read.csv(file.path(discharge_dir, "HF00402_v14.csv")) %>%
  mutate(
    date = as.Date(DATE, "%m/%d/%Y"),
    SITECODE = standardize_site_code(SITECODE)
  ) %>%
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    WATERYEAR >= WY_START,
    WATERYEAR <= WY_END
  ) %>%
  select(SITECODE, date, MEAN_Q, WATERYEAR)

ec_file <- file.path(ec_dir, "CF01201_v4.txt")
if (!file.exists(ec_file)) {
  stop("Missing required EC file: ", ec_file)
}
EC <- read_csv(ec_file, show_col_types = FALSE)

EC_daily <- EC %>%
  mutate(date = as.Date(DATE_TIME)) %>%
  mutate(SITECODE = standardize_site_code(SITECODE)) %>%
  filter(SITECODE %in% SITE_ORDER_CHEMISTRY) %>%
  group_by(SITECODE, date) %>%
  summarise(daily_SC = mean(EC_INST, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(daily_SC))

EC_Q <- left_join(
  EC_daily,
  discharge %>% select(SITECODE, date, MEAN_Q, WATERYEAR),
  by = c("SITECODE", "date")
) %>%
  filter(!is.na(MEAN_Q), !is.na(daily_SC)) %>%
  mutate(waterYear = get_water_year(date))

goodyears <- EC_Q %>%
  group_by(SITECODE, waterYear) %>%
  summarise(num_days = n_distinct(date), .groups = "drop") %>%
  filter(num_days >= CHS_MIN_DAYS_PER_WY)

EC_Q <- EC_Q %>%
  semi_join(goodyears, by = c("SITECODE", "waterYear")) %>%
  group_by(SITECODE) %>%
  mutate(
    SC_runoff = quantile(daily_SC, 0.01, na.rm = TRUE),
    SC_groundwater = quantile(daily_SC, 0.99, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    Q_baseflow = MEAN_Q * (daily_SC - SC_runoff) / (SC_groundwater - SC_runoff),
    GW_prop = Q_baseflow / MEAN_Q,
    GW_prop = pmax(0, pmin(1, GW_prop))
  )

annual_bf_prop <- EC_Q %>%
  group_by(SITECODE, waterYear) %>%
  summarise(
    CHS = mean(GW_prop, na.rm = TRUE),
    median_bf = median(GW_prop, na.rm = TRUE),
    sd_bf = sd(GW_prop, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  )

write.csv(
  annual_bf_prop,
  file.path(output_dir, "annual_gw_prop.csv"),
  row.names = FALSE
)

write.csv(
  goodyears %>% arrange(SITECODE, waterYear),
  file.path(OUT_MET_SUPPORT_DIR, "chs_wy_day_counts_kept.csv"),
  row.names = FALSE
)

# ---- part 2: isotope-derived site metrics (mtt/fyw/dr) ----

mtt_fyw <- read_csv(
  file.path(isotope_dir, "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  mutate(
    MTT2L_val = if ("MTT2L" %in% names(.)) suppressWarnings(as.numeric(MTT2L)) else NA_real_,
    MTT2H_val = if ("MTT2H" %in% names(.)) suppressWarnings(as.numeric(MTT2H)) else NA_real_,
    MTT2M_val = if ("MTT2M" %in% names(.)) suppressWarnings(as.numeric(MTT2M)) else NA_real_
  ) %>%
  mutate(
    MTT1 = suppressWarnings(as.numeric(MTT1)),
    MTT2 = suppressWarnings(as.numeric(dplyr::coalesce(
      MTT2M_val,
      rowMeans(cbind(
        MTT2L_val,
        MTT2H_val
      ), na.rm = TRUE)
    ))),
    MTT2 = ifelse(is.nan(MTT2), NA_real_, MTT2),
    Fyw = suppressWarnings(as.numeric(FYWM))
  ) %>%
  select(site, MTT1, MTT2, Fyw) %>%
  filter(!is.na(site), site != "")

damping <- read_csv(
  file.path(isotope_dir, "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  select(site, DR = DR_Overall, DR_err = DR__err)

isotope_metrics <- mtt_fyw %>%
  full_join(damping, by = "site") %>%
  mutate(site = factor(site, levels = SITE_ORDER_ALL)) %>%
  arrange(site)

write.csv(
  isotope_metrics,
  file.path(output_dir, "isotope_metrics_site.csv"),
  row.names = FALSE
)
