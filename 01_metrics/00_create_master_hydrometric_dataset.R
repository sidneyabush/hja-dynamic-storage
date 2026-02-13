# -----------------------------------------------------------------------------
# Create Master Hydrometric Dataset
# -----------------------------------------------------------------------------
# Purpose: Combine and interpolate meteorological station data to create
#          watershed-level daily meteorological datasets with P, T, RH, NR, VPD, Q
#
# Workflow:
#   Load meteorological data from multiple stations
#   Interpolate missing values using OLS (pairs) or MLR (triplets)
#   Aggregate station data to watershed level
#   Add discharge data
#   Create GSLOOK composite from component watersheds
#
# Inputs (from all_hydromet/):
#   - Temperature_original_&_filled_1979_2023_v2.csv
#   - Precipitation_original_&_filled_1979_2023.csv
#   - MS00102_v9.csv (relative humidity)
#   - MS05025_v3.csv (net radiation)
#   - MS00403_v2.csv (Mack Creek precipitation)
#   - HF00402_v14.csv (discharge)
#   - drainage_area.csv
#
# Outputs:
#   - watersheds_met_q.csv: Daily P, T, RH, NR, VPD, Q for all watersheds
#
# Methods preserved from original 1440-line script. Helper functions are in
# helpers/hydromet_utils.R for maintainability.
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

rm(list = ls())

# -----------------------------------------------------------------------------
# SOURCE UTILITIES AND CONFIGURATION
# -----------------------------------------------------------------------------

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

# Source helper functions
helpers_path <- file.path(script_dir, "helpers", "hydromet_utils.R")
if (!file.exists(helpers_path)) {
  helpers_path <- file.path(dirname(script_dir), "helpers", "hydromet_utils.R")
}
source(helpers_path)

# Source configuration (paths, sites, water year range)
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

# -----------------------------------------------------------------------------
# DIRECTORIES
# -----------------------------------------------------------------------------

met_dir <- MET_DIR
discharge_dir <- DISCHARGE_DIR
catchment_dir <- CATCHMENT_CHARACTERISTICS_DIR
output_dir <- OUT_MET_SUPPORT_DIR
exploratory_plot_dir <- EXPLORATORY_ET_METHODS_DIR
wy_start_date <- as.Date(sprintf("%d-10-01", WY_START - 1))
wy_end_date <- as.Date(sprintf("%d-09-30", WY_END))

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# SITE MAPPING: Watershed -> Met Stations
# -----------------------------------------------------------------------------
# Each watershed is mapped to one or more met stations for each variable.
# - Single station: Use that station's data
# - Two stations: OLS interpolation between the pair
# - Three stations: Multiple regression interpolation

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

# -----------------------------------------------------------------------------
# LOAD METEOROLOGICAL DATA
# -----------------------------------------------------------------------------


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

# Add Mack Creek precipitation
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

# -----------------------------------------------------------------------------
# COMBINE AND FILTER DATA
# -----------------------------------------------------------------------------


combined_met <- Temp %>%
  full_join(Precip, by = c("DATE", "SITECODE")) %>%
  full_join(RH, by = c("DATE", "SITECODE")) %>%
  full_join(NetRad, by = c("DATE", "SITECODE")) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  arrange(SITECODE, DATE)

# -----------------------------------------------------------------------------
# INTERPOLATE MISSING VALUES
# -----------------------------------------------------------------------------


# Extract station groups that need interpolation
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

# Clean duplicates
combined_met_clean <- combined_met %>%
  group_by(DATE, SITECODE) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# Variables to process (VPD computed internally)
variables <- c("T_C", "P_mm_d", "RH_d_pct", "NR_Wm2_d")

# Run interpolation
results <- process_station_groups(
  combined_met_clean,
  station_groups,
  variables,
  plot_dir = exploratory_plot_dir
)

interpolated_data <- results$data


# -----------------------------------------------------------------------------
# CREATE WATERSHED DATASETS
# -----------------------------------------------------------------------------


# Variables including VPD
watershed_variables <- c(variables, "VPD_kPa")

# Aggregate to watersheds
watershed_datasets <- create_watershed_datasets(interpolated_data, site_mapping, watershed_variables)

# Combine all watersheds
all_watersheds_data <- bind_rows(watershed_datasets)

# -----------------------------------------------------------------------------
# ADD DISCHARGE DATA
# -----------------------------------------------------------------------------


# Load drainage areas
da_df <- read_csv(resolve_drainage_area_file(), show_col_types = FALSE) %>%
  mutate(SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK)))

# Load and process discharge
discharge <- read_csv(file.path(discharge_dir, "HF00402_v14.csv"), show_col_types = FALSE) %>%
  mutate(
    DATE = parse_my_date(DATE),
    SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK))
  ) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  filter(!SITECODE %in% SITE_EXCLUDE_RAW) %>%
  group_by(DATE, SITECODE) %>%
  summarise(MEAN_Q = sum(MEAN_Q, na.rm = TRUE), .groups = "drop")

# Add discharge to watersheds
watershed_datasets <- add_discharge_to_watersheds(watershed_datasets, discharge, da_df)

# Rebuild master table
all_watersheds_data <- bind_rows(watershed_datasets)

# Clean up column names
if ("SITECODE.x" %in% names(all_watersheds_data)) {
  all_watersheds_data <- all_watersheds_data %>%
    select(-any_of(c("SITECODE", "SITECODE.y"))) %>%
    rename(SITECODE = SITECODE.x) %>%
    select(DATE, SITECODE, everything())
}

# -----------------------------------------------------------------------------
# CREATE GSLOOK COMPOSITE
# -----------------------------------------------------------------------------


# Component watersheds for Lookout
gslook_components <- GSLOOK_COMPOSITE_COMPONENT_SITES

# Build composite met variables (average of component watersheds)
gslook_full_df <- all_watersheds_data %>%
  filter(SITECODE %in% gslook_components) %>%
  group_by(DATE) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(SITECODE = "GSLOOK") %>%
  select(DATE, SITECODE, everything())

# Get GSLOOK discharge for the composite
gslook_q_full <- discharge %>%
  filter(SITECODE == "GSLOOK") %>%
  left_join(da_df, by = "SITECODE") %>%
  mutate(
    Q_mm_d = MEAN_Q * 0.0283168 * 86400 / DA_M2,
    SITECODE = "GSLOOK"
  ) %>%
  select(DATE, SITECODE, Q_mm_d)

# Attach discharge to composite
gslook_full_df <- gslook_full_df %>%
  select(-any_of("Q_mm_d")) %>%
  left_join(gslook_q_full, by = c("DATE", "SITECODE"))

# Add to master table
all_watersheds_data <- bind_rows(
  all_watersheds_data %>% filter(SITECODE != "GSLOOK"),
  gslook_full_df
)

# Clean up any duplicate Q_mm_d columns
if ("Q_mm_d.x" %in% names(all_watersheds_data)) {
  all_watersheds_data <- all_watersheds_data %>%
    mutate(Q_mm_d = coalesce(Q_mm_d.x, Q_mm_d)) %>%
    select(-any_of(c("Q_mm_d.x", "Q_mm_d.y")))
}

# Update watershed datasets list
watershed_datasets[["GSLOOK"]] <- gslook_full_df

# -----------------------------------------------------------------------------
# SAVE OUTPUT
# -----------------------------------------------------------------------------

output_file <- file.path(output_dir, "watersheds_met_q.csv")
write_csv(all_watersheds_data, output_file)
