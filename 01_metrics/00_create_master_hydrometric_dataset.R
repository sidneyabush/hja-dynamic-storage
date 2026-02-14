# Build watershed-scale daily met + discharge inputs for WY 1997-2020.
# Inputs: met_dir/MS00102_v9.csv; met_dir/MS05025_v3.csv; discharge_dir/HF00402_v14.csv.
# Author: Sidney Bush
# Date: 2026-02-13

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

rm(list = ls())

# Setup: load helpers and config
source("helpers/hydromet_utils.R")
source("config.R")

# Paths and water-year date bounds

met_dir <- MET_DIR
discharge_dir <- DISCHARGE_DIR
catchment_dir <- CATCHMENT_CHARACTERISTICS_DIR
output_dir <- OUT_MET_SUPPORT_DIR
exploratory_plot_dir <- EXPLORATORY_ET_METHODS_DIR
wy_start_date <- as.Date(sprintf("%d-10-01", WY_START - 1))
wy_end_date <- as.Date(sprintf("%d-09-30", WY_END))

# Make sure the output folder exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Site mapping: which met stations feed each watershed variable.
# One station = direct use, two stations = pair interpolation, three = multivariate interpolation.

site_mapping <- list(
  # Lower elevation watersheds - use PRIMET
  "GSWS09" = list(temp = c("PRIMET"), precip = c("PRIMET"),
                  rh = c("PRIMET", "CS2MET"), netrad = c("PRIMET")),
  "GSWS10" = list(temp = c("PRIMET"), precip = c("PRIMET"),
                  rh = c("PRIMET", "CS2MET"), netrad = c("PRIMET")),
  "GSWS01" = list(temp = c("PRIMET"), precip = c("PRIMET"),
                  rh = c("PRIMET", "CS2MET"), netrad = c("PRIMET")),
  "GSWS02" = list(temp = c("PRIMET"), precip = c("PRIMET"),
                  rh = c("PRIMET", "CS2MET"), netrad = c("PRIMET")),
  "GSWS03" = list(temp = c("PRIMET"), precip = c("PRIMET"),
                  rh = c("PRIMET", "CS2MET"), netrad = c("PRIMET")),

  # Mack Creek - use CENMET/UPLMET blend
  "GSMACK" = list(temp = c("CENMET", "UPLMET"), precip = c("GSMACK", "UPLMET"),
                  rh = c("CENMET", "UPLMET"), netrad = c("VANMET")),

  # Upper elevation watersheds - use H15MET/VANMET blend
  "GSWS06" = list(temp = c("H15MET", "VANMET"), precip = c("H15MET"),
                  rh = c("H15MET", "VANMET", "WS7MET"), netrad = c("VANMET")),
  "GSWS07" = list(temp = c("H15MET", "VANMET"), precip = c("H15MET"),
                  rh = c("H15MET", "VANMET", "WS7MET"), netrad = c("VANMET")),
  "GSWS08" = list(temp = c("H15MET", "VANMET"), precip = c("H15MET"),
                  rh = c("H15MET", "VANMET", "WS7MET"), netrad = c("VANMET")),

  # Lookout Creek tributaries
  "LONGER" = list(temp = c("CENMET"), precip = c("CENMET"),
                  rh = c("CENMET"), netrad = c("VANMET")),
  "COLD"   = list(temp = c("CENMET", "UPLMET"), precip = c("CENMET", "UPLMET"),
                  rh = c("CENMET", "UPLMET"), netrad = c("VANMET"))
)

# Load daily met inputs

# Temperature
Temp <- make_inter_long(
  "Temperature_original_&_filled_1979_2023_v2.csv",
  "Temp",
  met_dir,
  date_start = wy_start_date,
  date_end = wy_end_date
) %>%
  select(DATE, SITECODE, Temp) %>%
  rename(T_C = Temp)

# Precipitation
Precip <- make_inter_long(
  "Precipitation_original_&_filled_1979_2023.csv",
  "Precip",
  met_dir,
  date_start = wy_start_date,
  date_end = wy_end_date
) %>%
  select(DATE, SITECODE, Precip) %>%
  rename(P_mm_d = Precip)

# Add Mack Creek precipitation from its dedicated file
MACK_Precip <- read_mack_precip(
  "MS00403_v2.csv",
  met_dir,
  recode_map = SITECODE_RECODE_TO_GSMACK
)
Precip <- bind_rows(Precip, MACK_Precip)

# Relative humidity
RH <- read_csv(file.path(met_dir, "MS00102_v9.csv"), show_col_types = FALSE) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  select(SITECODE, DATE, RELHUM_MEAN_DAY) %>%
  rename(RH_d_pct = RELHUM_MEAN_DAY)

# Net radiation
NetRad <- read_csv(file.path(met_dir, "MS05025_v3.csv"), show_col_types = FALSE) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  select(SITECODE, DATE, NR_TOT_MEAN_DAY) %>%
  rename(NR_Wm2_d = NR_TOT_MEAN_DAY)

# Join met variables into one daily table and keep target dates

combined_met <- Temp %>%
  full_join(Precip, by = c("DATE", "SITECODE")) %>%
  full_join(RH, by = c("DATE", "SITECODE")) %>%
  full_join(NetRad, by = c("DATE", "SITECODE")) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  arrange(SITECODE, DATE)

# Fill gaps using the station-group interpolation rules above

# Pull pair/triplet groups used by the interpolation helper
station_groups <- extract_station_groups(site_mapping)

for (pair_name in names(station_groups$pairs)) {
  pair <- station_groups$pairs[[pair_name]]
  cat(sprintf("  - %s: %s and %s\n", pair_name, pair$site1, pair$site2))
}

for (triplet_name in names(station_groups$triplets)) {
  triplet <- station_groups$triplets[[triplet_name]]
  cat(sprintf("  - %s: %s, %s, and %s\n", triplet_name,
              triplet$site1, triplet$site2, triplet$site3))
}

# Remove duplicate site-date rows before interpolation
combined_met_clean <- combined_met %>%
  group_by(DATE, SITECODE) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# Variables passed to interpolation (VPD is computed inside helper)
variables <- c("T_C", "P_mm_d", "RH_d_pct", "NR_Wm2_d")

# Run interpolation and keep filled output
results <- process_station_groups(
  combined_met_clean,
  station_groups,
  variables,
  plot_dir = exploratory_plot_dir
)

interpolated_data <- results$data

# Build per-watershed daily met datasets

# Include VPD in watershed outputs
watershed_variables <- c(variables, "VPD_kPa")

# Aggregate station-level data to watershed-level series
watershed_datasets <- create_watershed_datasets(interpolated_data, site_mapping, watershed_variables)

# Combine all watershed tables into one master table
all_watersheds_data <- bind_rows(watershed_datasets)

# Add discharge (and area-normalized Q) to each watershed

# Drainage area lookup
da_df <- read_csv(resolve_drainage_area_file(), show_col_types = FALSE) %>%
  mutate(SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK)))

# Load discharge, standardize site codes, and keep target dates/sites
discharge <- read_csv(file.path(discharge_dir, "HF00402_v14.csv"), show_col_types = FALSE) %>%
  mutate(
    DATE = parse_my_date(DATE),
    SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK))
  ) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  filter(!SITECODE %in% SITE_EXCLUDE_RAW) %>%
  group_by(DATE, SITECODE) %>%
  summarise(MEAN_Q = sum(MEAN_Q, na.rm = TRUE), .groups = "drop")

# Join discharge into each watershed dataset
watershed_datasets <- add_discharge_to_watersheds(watershed_datasets, discharge, da_df)

# Rebuild the combined table after adding discharge
all_watersheds_data <- bind_rows(watershed_datasets)

# Clean up join artifact column names
if ("SITECODE.x" %in% names(all_watersheds_data)) {
  all_watersheds_data <- all_watersheds_data %>%
    select(-any_of(c("SITECODE", "SITECODE.y"))) %>%
    rename(SITECODE = SITECODE.x) %>%
    select(DATE, SITECODE, everything())
}

# Build GSLOOK composite from configured component watersheds

# Component watersheds used for GSLOOK
gslook_components <- GSLOOK_COMPOSITE_COMPONENT_SITES

# Average component met variables by date to create GSLOOK met series
gslook_full_df <- all_watersheds_data %>%
  filter(SITECODE %in% gslook_components) %>%
  group_by(DATE) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(SITECODE = "GSLOOK") %>%
  select(DATE, SITECODE, everything())

# Pull GSLOOK discharge and convert to mm/day
gslook_q_full <- discharge %>%
  filter(SITECODE == "GSLOOK") %>%
  left_join(da_df, by = "SITECODE") %>%
  mutate(
    Q_mm_d = MEAN_Q * 0.0283168 * 86400 / DA_M2,
    SITECODE = "GSLOOK"
  ) %>%
  select(DATE, SITECODE, Q_mm_d)

# Attach discharge to the GSLOOK composite table
gslook_full_df <- gslook_full_df %>%
  select(-any_of("Q_mm_d")) %>%
  left_join(gslook_q_full, by = c("DATE", "SITECODE"))

# Replace any existing GSLOOK rows with the composite rows
all_watersheds_data <- bind_rows(
  all_watersheds_data %>% filter(SITECODE != "GSLOOK"),
  gslook_full_df
)

# Resolve any duplicate Q_mm_d columns after join/bind steps
if ("Q_mm_d.x" %in% names(all_watersheds_data)) {
  all_watersheds_data <- all_watersheds_data %>%
    mutate(Q_mm_d = coalesce(Q_mm_d.x, Q_mm_d)) %>%
    select(-any_of(c("Q_mm_d.x", "Q_mm_d.y")))
}

# Keep GSLOOK in the per-watershed list for downstream use
watershed_datasets[["GSLOOK"]] <- gslook_full_df

# Save the final daily watershed met+Q table

output_file <- file.path(output_dir, "watersheds_met_q.csv")
write_csv(all_watersheds_data, output_file)
