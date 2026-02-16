# Calculate mobile-storage metrics (CHS + isotope metrics) by site and water year.
# Inputs: DISCHARGE_DIR/HF00402_v14.csv; EC_DIR/CF01201_v3.txt; ISOTOPE_DIR/MTT_FYW.csv; ISOTOPE_DIR/DampingRatios_2025-07-07.csv.
# Author: Keira Johnson (original CHS), Sidney Bush
# Date: 2026-02-14

library(dplyr)
library(readr)
library(lubridate)

rm(list = ls())

# Load project config
source("config.R")

output_dir <- OUT_MET_MOBILE_DIR
discharge_dir <- DISCHARGE_DIR
ec_dir <- EC_DIR
isotope_dir <- ISOTOPE_DIR

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- PART 1: CHS from specific conductance and discharge ----

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

# EC source is comma-delimited in current workflow input
EC <- read.delim(file.path(ec_dir, "CF01201_v4.txt"), sep = ",")

EC_daily <- EC %>%
  mutate(date = as.Date(DATE_TIME, "%Y-%m-%d")) %>%
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
  filter(num_days >= 365)

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

# ---- PART 2: Isotope-derived site metrics (MTT/Fyw/DR) ----

mtt_fyw <- read_csv(
  file.path(isotope_dir, "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  select(site, MTT = MTTM, Fyw = FYWM) %>%
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
