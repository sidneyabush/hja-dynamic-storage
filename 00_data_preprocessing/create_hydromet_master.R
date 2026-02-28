# build catchment-scale daily met + discharge inputs for wy 1997-2020.
# inputs:
# met_dir/ms00102_v9.csv
# met_dir/ms05025_v3.csv
# discharge_dir/hf00402_v14.csv

# author: sidney bush
# date: 2026-02-13

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

rm(list = ls())

# setup: load helpers and config
source("helpers/hydromet_utils.R")
source("config.R")

# paths and water-year date bounds

met_dir <- MET_DIR
discharge_dir <- DISCHARGE_DIR
catchment_dir <- CATCHMENT_CHARACTERISTICS_DIR
output_dir <- OUT_MET_SUPPORT_DIR
exploratory_plot_dir <- EXPLORATORY_ET_METHODS_DIR
wy_start_date <- as.Date(sprintf("%d-10-01", WY_START - 1))
wy_end_date <- as.Date(sprintf("%d-09-30", WY_END))

# make sure the output folder exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# site mapping: which met stations feed each catchment variable.
# one station = direct use, two stations = pair interpolation, three = multivariate interpolation.

site_mapping <- list(
  # lower elevation catchments - use primet
  "GSWS09" = list(
    temp = c("PRIMET"),
    precip = c("PRIMET"),
    rh = c("PRIMET", "CS2MET"),
    netrad = c("PRIMET")
  ),
  "GSWS10" = list(
    temp = c("PRIMET"),
    precip = c("PRIMET"),
    rh = c("PRIMET", "CS2MET"),
    netrad = c("PRIMET")
  ),
  "GSWS01" = list(
    temp = c("PRIMET"),
    precip = c("PRIMET"),
    rh = c("PRIMET", "CS2MET"),
    netrad = c("PRIMET")
  ),
  "GSWS02" = list(
    temp = c("PRIMET"),
    precip = c("PRIMET"),
    rh = c("PRIMET", "CS2MET"),
    netrad = c("PRIMET")
  ),
  "GSWS03" = list(
    temp = c("PRIMET"),
    precip = c("PRIMET"),
    rh = c("PRIMET", "CS2MET"),
    netrad = c("PRIMET")
  ),

  # mack creek - use cenmet/uplmet blend
  "GSMACK" = list(
    temp = c("CENMET", "UPLMET"),
    precip = c("GSMACK", "UPLMET"),
    rh = c("CENMET", "UPLMET"),
    netrad = c("VANMET")
  ),

  # upper elevation catchments - use h15met/vanmet blend
  "GSWS06" = list(
    temp = c("H15MET", "VANMET"),
    precip = c("H15MET"),
    rh = c("H15MET", "VANMET", "WS7MET"),
    netrad = c("VANMET")
  ),
  "GSWS07" = list(
    temp = c("H15MET", "VANMET"),
    precip = c("H15MET"),
    rh = c("H15MET", "VANMET", "WS7MET"),
    netrad = c("VANMET")
  ),
  "GSWS08" = list(
    temp = c("H15MET", "VANMET"),
    precip = c("H15MET"),
    rh = c("H15MET", "VANMET", "WS7MET"),
    netrad = c("VANMET")
  ),

  # lookout creek tributaries
  "LONGER" = list(
    temp = c("CENMET"),
    precip = c("CENMET"),
    rh = c("CENMET"),
    netrad = c("VANMET")
  ),
  "COLD" = list(
    temp = c("CENMET", "UPLMET"),
    precip = c("CENMET", "UPLMET"),
    rh = c("CENMET", "UPLMET"),
    netrad = c("VANMET")
  )
)

# load daily met inputs

# temperature
Temp <- make_inter_long(
  "Temperature_original_&_filled_1979_2023_v2.csv",
  "Temp",
  met_dir,
  date_start = wy_start_date,
  date_end = wy_end_date
) %>%
  select(DATE, SITECODE, Temp) %>%
  rename(T_C = Temp)

# precipitation
Precip <- make_inter_long(
  "Precipitation_original_&_filled_1979_2023.csv",
  "Precip",
  met_dir,
  date_start = wy_start_date,
  date_end = wy_end_date
) %>%
  select(DATE, SITECODE, Precip) %>%
  rename(P_mm_d = Precip)

# add mack creek precipitation from its dedicated file
MACK_Precip <- read_mack_precip(
  "MS00403_v2.csv",
  met_dir,
  recode_map = SITECODE_RECODE_TO_GSMACK
)
Precip <- bind_rows(Precip, MACK_Precip)

# relative humidity
RH <- read_csv(file.path(met_dir, "MS00102_v9.csv"), show_col_types = FALSE) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  select(SITECODE, DATE, RELHUM_MEAN_DAY) %>%
  rename(RH_d_pct = RELHUM_MEAN_DAY)

# net radiation
NetRad <- read_csv(
  file.path(met_dir, "MS05025_v3.csv"),
  show_col_types = FALSE
) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  select(SITECODE, DATE, NR_TOT_MEAN_DAY) %>%
  rename(NR_Wm2_d = NR_TOT_MEAN_DAY)

# join met variables into one daily table and keep target dates

combined_met <- Temp %>%
  full_join(Precip, by = c("DATE", "SITECODE")) %>%
  full_join(RH, by = c("DATE", "SITECODE")) %>%
  full_join(NetRad, by = c("DATE", "SITECODE")) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  arrange(SITECODE, DATE)

# fill gaps using the station-group interpolation rules above

# pull pair/triplet groups used by the interpolation helper
station_groups <- extract_station_groups(site_mapping)

for (pair_name in names(station_groups$pairs)) {
  pair <- station_groups$pairs[[pair_name]]
  cat(sprintf("  - %s: %s and %s\n", pair_name, pair$site1, pair$site2))
}

for (triplet_name in names(station_groups$triplets)) {
  triplet <- station_groups$triplets[[triplet_name]]
  cat(sprintf(
    "  - %s: %s, %s, and %s\n",
    triplet_name,
    triplet$site1,
    triplet$site2,
    triplet$site3
  ))
}

# remove duplicate site-date rows before interpolation
combined_met_clean <- combined_met %>%
  group_by(DATE, SITECODE) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# variables passed to interpolation (vpd is computed inside helper)
variables <- c("T_C", "P_mm_d", "RH_d_pct", "NR_Wm2_d")

# run interpolation and keep filled output
results <- process_station_groups(
  combined_met_clean,
  station_groups,
  variables,
  plot_dir = exploratory_plot_dir
)

interpolated_data <- results$data

# build per-catchment daily met datasets

# include vpd in catchment outputs
catchment_variables <- c(variables, "VPD_kPa")

# aggregate station-level data to catchment-level series
catchment_datasets <- create_catchment_datasets(
  interpolated_data,
  site_mapping,
  catchment_variables
)

# combine all catchment tables into one master table
all_catchments_data <- bind_rows(catchment_datasets)

# add discharge (and area-normalized q) to each catchment

# drainage area lookup
da_df <- read_csv(resolve_drainage_area_file(), show_col_types = FALSE) %>%
  mutate(SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK)))

# load discharge, standardize site codes, and keep target dates/sites
discharge <- read_csv(
  file.path(discharge_dir, "HF00402_v14.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    DATE = parse_my_date(DATE),
    SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK))
  ) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  filter(!SITECODE %in% SITE_EXCLUDE_RAW) %>%
  group_by(DATE, SITECODE) %>%
  summarise(MEAN_Q = sum(MEAN_Q, na.rm = TRUE), .groups = "drop")

# join discharge into each catchment dataset
catchment_datasets <- add_discharge_to_catchments(
  catchment_datasets,
  discharge,
  da_df
)

# rebuild the combined table after adding discharge
all_catchments_data <- bind_rows(catchment_datasets)

# clean up join artifact column names
if ("SITECODE.x" %in% names(all_catchments_data)) {
  all_catchments_data <- all_catchments_data %>%
    select(-any_of(c("SITECODE", "SITECODE.y"))) %>%
    rename(SITECODE = SITECODE.x) %>%
    select(DATE, SITECODE, everything())
}

# build gslook composite from configured component catchments

# component catchments used for gslook
gslook_components <- GSLOOK_COMPOSITE_COMPONENT_SITES

# average component met variables by date to create gslook met series
gslook_full_df <- all_catchments_data %>%
  filter(SITECODE %in% gslook_components) %>%
  group_by(DATE) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(SITECODE = "GSLOOK") %>%
  select(DATE, SITECODE, everything())

# pull gslook discharge and convert to mm/day
gslook_q_full <- discharge %>%
  filter(SITECODE == "GSLOOK") %>%
  left_join(da_df, by = "SITECODE") %>%
  mutate(
    Q_mm_d = MEAN_Q * 0.0283168 * 86400 / DA_M2 * 1000,
    SITECODE = "GSLOOK"
  ) %>%
  select(DATE, SITECODE, Q_mm_d)

# attach discharge to the gslook composite table
gslook_full_df <- gslook_full_df %>%
  select(-any_of("Q_mm_d")) %>%
  left_join(gslook_q_full, by = c("DATE", "SITECODE"))

# replace any existing gslook rows with the composite rows
all_catchments_data <- bind_rows(
  all_catchments_data %>% filter(SITECODE != "GSLOOK"),
  gslook_full_df
)

# resolve any duplicate q_mm_d columns after join/bind steps
if ("Q_mm_d.x" %in% names(all_catchments_data)) {
  all_catchments_data <- all_catchments_data %>%
    mutate(Q_mm_d = coalesce(Q_mm_d.x, Q_mm_d)) %>%
    select(-any_of(c("Q_mm_d.x", "Q_mm_d.y")))
}

# keep gslook in the per-catchment list for downstream use
catchment_datasets[["GSLOOK"]] <- gslook_full_df

# save the final daily catchment met+q table

output_file <- file.path(output_dir, "catchments_met_q.csv")
write_csv(all_catchments_data, output_file)
