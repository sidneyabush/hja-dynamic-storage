# build daily catchment met and discharge inputs for WY 1997-2020

# inputs:
# met_dir/Temperature_original_&_filled_1979_2023_v2.csv
# met_dir/Precipitation_original_&_filled_1979_2023.csv
# met_dir/SWE_original_&_filled_1997_2023_v5.csv
# met_dir/MS00102_v9.csv
# met_dir/MS05025_v3.csv
# discharge_dir/HF00402_v14.csv
# catchment_characteristics_dir/drainage_area.csv

# outputs:
# out_met_support_dir/catchments_met_q.csv

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(readr, dplyr, lubridate, tidyr, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

met_dir <- MET_DIR
discharge_dir <- DISCHARGE_DIR
catchment_dir <- CATCHMENT_CHARACTERISTICS_DIR
output_dir <- OUT_MET_SUPPORT_DIR
wy_start_date <- as.Date(sprintf("%d-10-01", WY_START - 1))
wy_end_date <- as.Date(sprintf("%d-09-30", WY_END))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# assign met stations by catchment and variable
# single stations are used directly, paired stations fill each other's gaps
site_mapping <- get_met_station_assignments()

# met inputs start from the Attias (2025) gap-filled station files
Temp <- make_inter_long(
  "Temperature_original_&_filled_1979_2023_v2.csv",
  "Temp",
  met_dir,
  date_start = wy_start_date,
  date_end = wy_end_date
) %>%
  select(DATE, SITECODE, Temp) %>%
  rename(T_C = Temp)

Precip <- make_inter_long(
  "Precipitation_original_&_filled_1979_2023.csv",
  "Precip",
  met_dir,
  date_start = wy_start_date,
  date_end = wy_end_date
) %>%
  select(DATE, SITECODE, Precip) %>%
  rename(P_mm_d = Precip)

MACK_Precip <- read_mack_precip(
  "MS00403_v2.csv",
  met_dir,
  recode_map = SITECODE_RECODE_TO_GSMACK
)
Precip <- bind_rows(Precip, MACK_Precip)

RH <- read_csv(file.path(met_dir, "MS00102_v9.csv"), show_col_types = FALSE) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  select(SITECODE, DATE, RELHUM_MEAN_DAY) %>%
  rename(RH_d_pct = RELHUM_MEAN_DAY)

NetRad <- read_csv(
  file.path(met_dir, "MS05025_v3.csv"),
  show_col_types = FALSE
) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  select(SITECODE, DATE, NR_TOT_MEAN_DAY) %>%
  rename(NR_Wm2_d = NR_TOT_MEAN_DAY)

combined_met <- Temp %>%
  full_join(Precip, by = c("DATE", "SITECODE")) %>%
  full_join(RH, by = c("DATE", "SITECODE")) %>%
  full_join(NetRad, by = c("DATE", "SITECODE")) %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date) %>%
  arrange(SITECODE, DATE)

station_groups <- extract_station_groups(site_mapping)

combined_met_clean <- combined_met %>%
  group_by(DATE, SITECODE) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# VPD is calculated after temperature and relative humidity have been filled
variables <- c("T_C", "P_mm_d", "RH_d_pct", "NR_Wm2_d")

results <- process_station_groups(
  combined_met_clean,
  station_groups,
  variables
)

interpolated_data <- results$data

catchment_variables <- c(variables, "VPD_kPa")

catchment_datasets <- create_catchment_datasets(
  interpolated_data,
  site_mapping,
  catchment_variables
)

all_catchments_data <- bind_rows(catchment_datasets)

da_df <- read_csv(resolve_drainage_area_file(), show_col_types = FALSE) %>%
  mutate(SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK)))

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

catchment_datasets <- add_discharge_to_catchments(
  catchment_datasets,
  discharge,
  da_df
)

all_catchments_data <- bind_rows(catchment_datasets)

if ("SITECODE.x" %in% names(all_catchments_data)) {
  all_catchments_data <- all_catchments_data %>%
    select(-any_of(c("SITECODE", "SITECODE.y"))) %>%
    rename(SITECODE = SITECODE.x) %>%
    select(DATE, SITECODE, everything())
}

# GSLOOK is represented by a met composite of the configured tributary catchments, paired with the observed GSLOOK discharge record
gslook_components <- GSLOOK_COMPOSITE_COMPONENT_SITES

gslook_full_df <- all_catchments_data %>%
  filter(SITECODE %in% gslook_components) %>%
  group_by(DATE) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(SITECODE = "GSLOOK") %>%
  select(DATE, SITECODE, everything())

gslook_q_full <- discharge %>%
  filter(SITECODE == "GSLOOK") %>%
  left_join(da_df, by = "SITECODE") %>%
  mutate(
    Q_mm_d = MEAN_Q * 0.0283168 * 86400 / DA_M2 * 1000,
    SITECODE = "GSLOOK"
  ) %>%
  select(DATE, SITECODE, Q_mm_d)

gslook_full_df <- gslook_full_df %>%
  select(-any_of("Q_mm_d")) %>%
  left_join(gslook_q_full, by = c("DATE", "SITECODE"))

all_catchments_data <- bind_rows(
  all_catchments_data %>% filter(SITECODE != "GSLOOK"),
  gslook_full_df
)

if ("Q_mm_d.x" %in% names(all_catchments_data)) {
  all_catchments_data <- all_catchments_data %>%
    mutate(Q_mm_d = coalesce(Q_mm_d.x, Q_mm_d)) %>%
    select(-any_of(c("Q_mm_d.x", "Q_mm_d.y")))
}

catchment_datasets[["GSLOOK"]] <- gslook_full_df

output_file <- file.path(output_dir, "catchments_met_q.csv")
write_csv(all_catchments_data, output_file)
