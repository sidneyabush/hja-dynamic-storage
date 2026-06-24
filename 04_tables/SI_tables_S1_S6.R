# write supporting information Tables S1 to S6

# inputs:
# discharge_dir/HF00402_v14.csv
# ec_dir/CF00201_v7.csv
# stream_temp_dir/HT00451_v10.txt
# out_met_support_dir/catchments_met_q.csv
# outputs/metrics/mobile/isotope_metrics_site.csv
# catchment_characteristics_dir/catchment_char.csv
# met_dir/Temperature_original_&_filled_1979_2023_v2.csv
# met_dir/Precipitation_original_&_filled_1979_2023.csv
# met_dir/SWE_original_&_filled_1997_2023_v5.csv
# met_dir/MS00102_v9.csv
# met_dir/MS05025_v3.csv
# met_dir/MS00403_v2.csv
# isotope_dir/MTT_FYW.csv
# isotope_dir/DampingRatios_2025-07-07.csv

# outputs:
# figs_tables_pub/supp/TableS1_data_periods_of_record.csv
# figs_tables_pub/supp/TableS2_catchment_physiography_land_use.csv
# figs_tables_pub/supp/TableS3_catchment_geology_landslides.csv
# figs_tables_pub/supp/TableS4_met_station_record_summary.csv
# figs_tables_pub/supp/TableS5_met_station_assignments.csv
# figs_tables_pub/supp/TableS6_isotope_metrics.csv

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(dplyr, readr, tidyr, tibble, lubridate, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

dir.create(MS_TABLES_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

# study period used to summarize data availability
study_start <- as.Date(sprintf("%d-10-01", WY_START - 1))
study_end <- as.Date(sprintf("%d-09-30", WY_END))

# formatting functions used by several SI tables
clean_num <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9.-]", "", as.character(x))))
}

site_factor <- function(x) {
  factor(as.character(x), levels = SITE_ORDER_HYDROMETRIC)
}

summarise_period <- function(df, site_col, date_col, value_col) {
  df %>%
    transmute(
      site = standardize_site_code(as.character({{ site_col }})),
      date = as.Date({{ date_col }}),
      value = suppressWarnings(as.numeric({{ value_col }}))
    ) %>%
    filter(
      site %in% SITE_ORDER_HYDROMETRIC,
      date >= study_start,
      date <= study_end,
      is.finite(value)
    ) %>%
    mutate(water_year = get_water_year(date)) %>%
    group_by(site) %>%
    summarise(
      first_date = min(date, na.rm = TRUE),
      last_date = max(date, na.rm = TRUE),
      n_records = n(),
      n_water_years = n_distinct(water_year),
      .groups = "drop"
    )
}

add_prefix <- function(df, prefix) {
  df %>%
    rename_with(~ paste0(prefix, "_", .x), -site)
}

read_csv_dates <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(DATE = parse_my_date(DATE))
}

summarise_wide_station_file <- function(path, variable, use_interpolated = FALSE) {
  raw <- read_csv_dates(path)
  value_cols <- setdiff(names(raw), c("DATE", "WY"))
  if (isTRUE(use_interpolated)) {
    value_cols <- value_cols[grepl("_inter$", value_cols)]
  }

  raw %>%
    select(DATE, all_of(value_cols)) %>%
    pivot_longer(
      cols = all_of(value_cols),
      names_to = "station",
      values_to = "value"
    ) %>%
    mutate(
      station = gsub("_inter$", "", station),
      variable = variable,
      source_file = basename(path),
      date = as.Date(DATE),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    filter(is.finite(value)) %>%
    group_by(variable, source_file, station) %>%
    summarise(
      first_date = min(date, na.rm = TRUE),
      last_date = max(date, na.rm = TRUE),
      n_records = n(),
      .groups = "drop"
    )
}

summarise_long_station_file <- function(path, variable, value_col) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      date = parse_my_date(DATE),
      station = as.character(SITECODE),
      value = suppressWarnings(as.numeric(.data[[value_col]])),
      variable = variable,
      source_file = basename(path)
    ) %>%
    filter(is.finite(value)) %>%
    group_by(variable, source_file, station) %>%
    summarise(
      first_date = min(date, na.rm = TRUE),
      last_date = max(date, na.rm = TRUE),
      n_records = n(),
      .groups = "drop"
    )
}

# Table S1: data periods of record used in the analysis
discharge_period <- read_csv(file.path(DISCHARGE_DIR, "HF00402_v14.csv"), show_col_types = FALSE) %>%
  mutate(
    SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK)),
    DATE = parse_my_date(DATE)
  ) %>%
  summarise_period(SITECODE, DATE, MEAN_Q) %>%
  add_prefix("discharge")

chemistry_period <- read_csv(file.path(EC_DIR, "CF00201_v7.csv"), show_col_types = FALSE) %>%
  mutate(DATE = as.Date(DATE_TIME)) %>%
  summarise_period(SITECODE, DATE, CA) %>%
  add_prefix("chemistry_ca")

temperature_period <- read_csv(file.path(STREAM_TEMP_DIR, "HT00451_v10.txt"), show_col_types = FALSE) %>%
  mutate(DATE = as.Date(DATE_TIME)) %>%
  summarise_period(SITECODE, DATE, WATERTEMP_MEAN) %>%
  add_prefix("stream_temperature")

met_period <- read_csv(file.path(OUT_MET_SUPPORT_DIR, "catchments_met_q.csv"), show_col_types = FALSE) %>%
  mutate(DATE = parse_my_date(DATE)) %>%
  summarise_period(SITECODE, DATE, P_mm_d) %>%
  add_prefix("met_forcing")

isotope_metrics <- read_csv(
  file.path(OUT_MET_MOBILE_DIR, "isotope_metrics_site.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  transmute(
    site,
    has_damping_ratio = is.finite(DR),
    has_young_water_fraction = is.finite(Fyw),
    has_mean_transit_time = is.finite(MTT)
  )

table_s1 <- tibble(site = SITE_ORDER_HYDROMETRIC) %>%
  left_join(discharge_period, by = "site") %>%
  left_join(chemistry_period, by = "site") %>%
  left_join(temperature_period, by = "site") %>%
  left_join(met_period, by = "site") %>%
  left_join(isotope_metrics, by = "site") %>%
  arrange(site_factor(site))

write_csv(table_s1, file.path(MS_TABLES_SUPP_DIR, "TableS1_data_periods_of_record.csv"))

# Tables S2 and S3: catchment physiography, geology, and landslides
catchment_chars <- read_csv(
  resolve_catchment_characteristics_file(),
  show_col_types = FALSE
) %>%
  mutate(
    site = standardize_site_code(Site),
    Area_km2 = clean_num(Area_km2)
  ) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

table_s2 <- catchment_chars %>%
  transmute(
    site,
    area_ha = Area_km2 * 100,
    elevation_min_m = clean_num(Elevation_min_m),
    elevation_max_m = clean_num(Elevation_max_m),
    elevation_mean_m = clean_num(Elevation_mean_m),
    basin_slope_mean_deg = clean_num(Slope_mean),
    basin_slope_sd_deg = clean_num(Slope_Std),
    aspect_mean_deg = clean_num(Aspect_Mean_deg),
    harvest_pct = clean_num(Harvest),
    stand_age_yr = clean_num(Age)
  ) %>%
  arrange(site_factor(site))

write_csv(table_s2, file.path(MS_TABLES_SUPP_DIR, "TableS2_catchment_physiography_land_use.csv"))

table_s3 <- catchment_chars %>%
  transmute(
    site,
    lava_1_pct = clean_num(Lava1_per),
    lava_2_pct = clean_num(Lava2_per),
    ash_pct = clean_num(Ash_Per),
    pyroclastic_pct = clean_num(Pyro_per),
    landslide_young_pct = clean_num(Landslide_Young),
    landslide_mod_pct = clean_num(Landslide_Mod),
    landslide_old_pct = clean_num(Landslide_Old),
    landslide_total_pct = clean_num(Landslide_Total)
  ) %>%
  arrange(site_factor(site))

write_csv(table_s3, file.path(MS_TABLES_SUPP_DIR, "TableS3_catchment_geology_landslides.csv"))

# Table S4: meteorological and climate source station periods
table_s4 <- bind_rows(
  summarise_wide_station_file(
    file.path(MET_DIR, "Temperature_original_&_filled_1979_2023_v2.csv"),
    "air_temperature",
    use_interpolated = TRUE
  ),
  summarise_wide_station_file(
    file.path(MET_DIR, "Precipitation_original_&_filled_1979_2023.csv"),
    "precipitation",
    use_interpolated = TRUE
  ),
  summarise_wide_station_file(
    file.path(MET_DIR, "SWE_original_&_filled_1997_2023_v5.csv"),
    "snow_water_equivalent",
    use_interpolated = TRUE
  ),
  summarise_long_station_file(
    file.path(MET_DIR, "MS00102_v9.csv"),
    "relative_humidity",
    "RELHUM_MEAN_DAY"
  ),
  summarise_long_station_file(
    file.path(MET_DIR, "MS05025_v3.csv"),
    "net_radiation",
    "NR_TOT_MEAN_DAY"
  ),
  summarise_long_station_file(
    file.path(MET_DIR, "MS00403_v2.csv"),
    "mack_creek_precipitation",
    "PRECIP_TOT_DAY"
  )
) %>%
  arrange(variable, station)

write_csv(table_s4, file.path(MS_TABLES_SUPP_DIR, "TableS4_met_station_record_summary.csv"))

# Table S5: catchment meteorological station assignments
format_station_list <- function(x) {
  paste(x, collapse = ", ")
}

site_mapping <- get_met_station_assignments()
table_s5 <- bind_rows(lapply(names(site_mapping), function(source_site) {
  site_info <- site_mapping[[source_site]]
  tibble(
    source_site = source_site,
    analysis_site = standardize_site_code(source_site),
    temperature_station = format_station_list(site_info$temp),
    precipitation_station = format_station_list(site_info$precip),
    relative_humidity_station = format_station_list(site_info$rh),
    net_radiation_station = format_station_list(site_info$netrad),
    note = ifelse(
      source_site %in% c("LONGER", "COLD"),
      "Used only in the Lookout Creek composite forcing record",
      NA_character_
    )
  )
})) %>%
  bind_rows(
    tibble(
      source_site = "GSLOOK",
      analysis_site = "Look",
      temperature_station = paste(GSLOOK_COMPOSITE_COMPONENT_SITES, collapse = ", "),
      precipitation_station = paste(GSLOOK_COMPOSITE_COMPONENT_SITES, collapse = ", "),
      relative_humidity_station = paste(GSLOOK_COMPOSITE_COMPONENT_SITES, collapse = ", "),
      net_radiation_station = paste(GSLOOK_COMPOSITE_COMPONENT_SITES, collapse = ", "),
      note = "Composite forcing record from configured Lookout Creek component catchments"
    )
  ) %>%
  mutate(
    analysis_site = factor(analysis_site, levels = c(SITE_ORDER_HYDROMETRIC, "LONGER", "COLD"))
  ) %>%
  arrange(analysis_site, source_site) %>%
  mutate(analysis_site = as.character(analysis_site))

write_csv(table_s5, file.path(MS_TABLES_SUPP_DIR, "TableS5_met_station_assignments.csv"))

# Table S6: isotope source values and site metrics used in models
# MTT and Fyw source values are shown beside the model values
mtt_fyw_raw <- read_csv(file.path(ISOTOPE_DIR, "MTT_FYW.csv"), show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  transmute(
    site,
    MTT1 = clean_num(MTT1),
    MTT1_SD = clean_num(MTT1_SD),
    MTT2L = clean_num(MTT2L),
    MTT2L_SD = clean_num(MTT2L_SD),
    MTT2H = clean_num(MTT2H),
    MTT2H_SD = clean_num(MTT2H_SD),
    MTTM = clean_num(MTTM),
    FYWL = clean_num(FYWL),
    FYWL_SD = clean_num(FYWL_SD),
    FYWH = clean_num(FYWH),
    FYWH_SD = clean_num(FYWH_SD),
    FYWM = clean_num(FYWM)
  )

# damping ratio source columns show which studies contribute to DR_Overall
damping_raw <- read_csv(
  file.path(ISOTOPE_DIR, "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  transmute(
    site,
    DR_Segura = clean_num(Segu),
    DR_Segura_err = clean_num(Seg_err),
    DR_McGuire = clean_num(McGuire),
    DR_McGuire_err = clean_num(McGui_err),
    DR_Ortega = clean_num(Ortega),
    DR_Ortega_err = clean_num(Orte_err),
    DR_Overall = clean_num(DR_Overall),
    DR_Overall_err = clean_num(DR__err)
  )

# isotope_metrics_site contains the model values used by the analysis workflow
isotope_final <- read_csv(
  file.path(OUT_MET_MOBILE_DIR, "isotope_metrics_site.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  transmute(
    site,
    MTT_final = clean_num(MTT),
    Fyw_final = clean_num(Fyw),
    DR_final = clean_num(DR),
    DR_final_err = clean_num(DR_err)
  )

table_s6 <- tibble(site = SITE_ORDER_HYDROMETRIC) %>%
  left_join(isotope_final, by = "site") %>%
  left_join(mtt_fyw_raw, by = "site") %>%
  left_join(damping_raw, by = "site") %>%
  arrange(site_factor(site))

write_csv(table_s6, file.path(MS_TABLES_SUPP_DIR, "TableS6_isotope_metrics.csv"))
