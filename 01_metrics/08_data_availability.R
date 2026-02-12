# -----------------------------------------------------------------------------
# Data Availability Summary - Sites, Metrics, and Date Ranges
# -----------------------------------------------------------------------------
# Purpose: Generate comprehensive tables showing data availability for:
#   1. Which sites can be used for each storage metric
#   2. Date ranges for each metric/site combination
#   3. Meteorological variables and their date ranges
#
# Output Tables:
#   - Site_Metric_Availability.csv: Sites × Metrics matrix
#   - Metric_Date_Ranges.csv: Date ranges by site and metric
#   - Met_Variables_Summary.csv: Meteorological data coverage
#
# Author: Sidney Bush
# Date: 2026-01-30
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(lubridate)

rm(list = ls())

# Source configuration (paths, site definitions, water year range)
# Get script directory (works with source() and Rscript)
script_dir <- tryCatch(
  {
    dirname(sys.frame(1)$ofile)
  },
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      getwd()
    }
  }
)
if (is.null(script_dir) || script_dir == "" || script_dir == ".") {
  script_dir <- getwd()
}

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
# SITE DEFINITIONS (extends config.R definitions with metadata)
# -----------------------------------------------------------------------------
# Complete site list with descriptions
# Uses SITE_NAMES from config.R, with additional metadata
site_info <- tribble(
  ~site_code , ~site_name        , ~hydrometric , ~chemistry , ~isotopes , ~notes                      ,
  "GSWS09"   , "Watershed 09"    , TRUE         , TRUE       , TRUE      , ""                          ,
  "GSWS10"   , "Watershed 10"    , TRUE         , TRUE       , TRUE      , ""                          ,
  "GSWS01"   , "Watershed 01"    , TRUE         , TRUE       , TRUE      , ""                          ,
  "Look"     , "Lookout"         , TRUE         , TRUE       , TRUE      , ""                          ,
  "GSWS02"   , "Watershed 02"    , TRUE         , TRUE       , TRUE      , ""                          ,
  "GSWS03"   , "Watershed 03"    , TRUE         , TRUE       , TRUE      , ""                          ,
  "MR"       , "McRae"           , FALSE        , FALSE      , TRUE      , "Isotope only (2022-2023)"  ,
  "GSWS06"   , "Watershed 06"    , TRUE         , FALSE      , FALSE     , "No chemistry/isotope data" ,
  "GSWS07"   , "Watershed 07"    , TRUE         , TRUE       , TRUE      , ""                          ,
  "GSWS08"   , "Watershed 08"    , TRUE         , TRUE       , TRUE      , ""                          ,
  "NC"       , "Nostoc"          , FALSE        , FALSE      , TRUE      , "Isotope only (2022-2023)"  ,
  "Mack"     , "Mack"            , TRUE         , TRUE       , TRUE      , ""                          ,
  "LC"       , "Longer"          , FALSE        , FALSE      , TRUE      , "Isotope only (2022-2023)"  ,
  "LO2"      , "Upper Lookout 2" , FALSE        , FALSE      , TRUE      , "Isotope only (2022-2023)"  ,
  "CC"       , "Cold Creek"      , FALSE        , FALSE      , TRUE      , "Isotope only (2022-2023)"  ,
  "LO1"      , "Upper Lookout 1" , FALSE        , FALSE      , TRUE      , "Isotope only (2022-2023)"
)

site_info <- site_info %>%
  filter(!site_code %in% c("MR", "NC", "LC", "LO2", "CC", "LO1"))

# -----------------------------------------------------------------------------
# DIRECTORIES (from config.R)
# -----------------------------------------------------------------------------
base_dir <- BASE_DATA_DIR
output_dir <- OUTPUT_DIR
support_dir <- OUT_MET_SUPPORT_DIR
if (!dir.exists(support_dir)) dir.create(support_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# STORAGE METRICS AVAILABILITY TABLE
# -----------------------------------------------------------------------------
# Storage metrics definitions
# NOTE: Q5norm, CV_Q5norm, and temperature are NOT storage metrics
metrics_info <- tribble(
  ~storage_type      , ~method                    , ~abbreviation , ~variable_name          , ~requires     ,
  "Dynamic"          , "Richards-Baker Index"     , "RBI"         , "RBI"                   , "hydrometric" ,
  "Dynamic"          , "Recession Curve Slope"    , "RCS"         , "recession_curve_slope" , "hydrometric" ,
  "Dynamic"          , "Flow Duration Curve"      , "FDC"         , "fdc_slope"             , "hydrometric" ,
  "Dynamic"          , "Storage-Discharge"        , "SD"          , "S_annual_mm"           , "hydrometric" ,
  "Mobile"           , "Mean Transit Time"        , "MTT"         , "MTT"                   , "isotopes"    ,
  "Mobile"           , "Young Water Fraction"     , "Fyw"         , "Fyw"                   , "isotopes"    ,
  "Mobile"           , "Chemical Hydrograph Sep." , "CHS"         , "mean_bf"               , "chemistry"   ,
  "Mobile"           , "Isotopic Damping Ratio"   , "DR"          , "DR"                    , "isotopes"    ,
  "Extended Dynamic" , "Water Balance"            , "WB"          , "DS_sum"                , "hydrometric"
)

# Response variables (NOT storage metrics, but used in analyses)
response_vars <- tribble(
  ~type      , ~method                 , ~abbreviation , ~variable_name       , ~requires     ,
  "Response" , "Max 7-day Stream Temp" , "T_7DMax"     , "T_7DMax"            , "stream_temp" ,
  "Response" , "Q5 of 7-day Discharge" , "Q_7Q5"       , "Q_7Q5"              , "hydrometric" ,
  "Response" , "Temp during Q5 Period" , "T_Q7Q5"      , "T_Q7Q5"             , "stream_temp"
)

# Create site × metric availability matrix
site_metric_matrix <- site_info %>%
  select(site_code, site_name, hydrometric, chemistry, isotopes) %>%
  crossing(metrics_info %>% select(abbreviation, requires)) %>%
  mutate(
    available = case_when(
      requires == "hydrometric" & hydrometric ~ TRUE,
      requires == "chemistry" & chemistry ~ TRUE,
      requires == "isotopes" & isotopes ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  select(site_code, site_name, abbreviation, available) %>%
  pivot_wider(names_from = abbreviation, values_from = available) %>%
  # Reorder columns
  select(site_code, site_name, RBI, RCS, FDC, SD, MTT, Fyw, CHS, DR, WB)


# -----------------------------------------------------------------------------
# DATE RANGES BY METRIC
# -----------------------------------------------------------------------------
# Load discharge data for hydrometric date ranges
discharge <- read_csv(
  file.path(base_dir, "Q", "HF00402_v14.csv"),
  show_col_types = FALSE
) %>%
  mutate(Date = as.Date(DATE, "%m/%d/%Y"))

hydro_ranges <- discharge %>%
  group_by(SITECODE) %>%
  summarise(
    hydro_start = min(Date, na.rm = TRUE),
    hydro_end = max(Date, na.rm = TRUE),
    hydro_n_days = sum(!is.na(MEAN_Q)),
    .groups = "drop"
  ) %>%
  rename(site_code = SITECODE)

# Load EC data for chemistry date ranges
ec_file <- file.path(base_dir, "EC", "CF01201_v3.txt")
if (file.exists(ec_file)) {
  ec_data <- read_delim(ec_file, delim = "\t", show_col_types = FALSE)

  # Check column names
  if ("DATE" %in% names(ec_data) & "SITECODE" %in% names(ec_data)) {
    chem_ranges <- ec_data %>%
      mutate(Date = as.Date(DATE)) %>%
      group_by(SITECODE) %>%
      summarise(
        chem_start = min(Date, na.rm = TRUE),
        chem_end = max(Date, na.rm = TRUE),
        chem_n_days = n(),
        .groups = "drop"
      ) %>%
      rename(site_code = SITECODE)
  } else {
    chem_ranges <- tibble(site_code = character())
  }
} else {
  chem_ranges <- tibble(site_code = character())
}

# Load water balance data
wb_file <- file.path(
  base_dir,
  "DynamicStorage",
  "daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv"
)
if (file.exists(wb_file)) {
  wb_data <- read_csv(wb_file, show_col_types = FALSE) %>%
    mutate(Date = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y")))

  wb_ranges <- wb_data %>%
    group_by(SITECODE) %>%
    summarise(
      wb_start = min(Date, na.rm = TRUE),
      wb_end = max(Date, na.rm = TRUE),
      wb_n_days = sum(!is.na(P_mm_d) & !is.na(Q_mm_d) & !is.na(ET_mm_d)),
      .groups = "drop"
    ) %>%
    rename(site_code = SITECODE)
} else {
  wb_ranges <- tibble(site_code = character())
}

# Combine all date ranges
date_ranges <- site_info %>%
  select(site_code, site_name) %>%
  left_join(hydro_ranges, by = "site_code") %>%
  left_join(chem_ranges, by = "site_code") %>%
  left_join(wb_ranges, by = "site_code")


# -----------------------------------------------------------------------------
# METEOROLOGICAL VARIABLES SUMMARY
# -----------------------------------------------------------------------------
# Load master met dataset
met_file <- file.path(OUT_MET_SUPPORT_DIR, "watersheds_met_q.csv")
if (!file.exists(met_file)) {
  met_file <- file.path(base_dir, "MET", "watersheds_met_data_q.csv")
}

if (file.exists(met_file)) {
  met_data <- read_csv(met_file, show_col_types = FALSE)

  # Get list of met variables (exclude date/site columns)
  exclude_cols <- c(
    "DATE",
    "Date",
    "date",
    "SITECODE",
    "site",
    "WATERYEAR",
    "wateryear"
  )
  met_vars <- setdiff(names(met_data), exclude_cols)

  # Parse date
  if ("DATE" %in% names(met_data)) {
    met_data$Date <- as.Date(
      met_data$DATE,
      tryFormats = c("%Y-%m-%d", "%m/%d/%Y")
    )
  }

  # Calculate coverage for each variable separately
  met_summary <- data.frame(
    variable = met_vars,
    n_obs = sapply(met_vars, function(v) sum(!is.na(met_data[[v]]))),
    pct_complete = sapply(met_vars, function(v) {
      round(100 * sum(!is.na(met_data[[v]])) / nrow(met_data), 1)
    }),
    min_date = sapply(met_vars, function(v) {
      valid <- !is.na(met_data[[v]])
      if (any(valid)) {
        as.character(min(met_data$Date[valid], na.rm = TRUE))
      } else {
        NA_character_
      }
    }),
    max_date = sapply(met_vars, function(v) {
      valid <- !is.na(met_data[[v]])
      if (any(valid)) {
        as.character(max(met_data$Date[valid], na.rm = TRUE))
      } else {
        NA_character_
      }
    }),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(n_obs))

  # Save met summary
  write.csv(
    met_summary,
    file.path(support_dir, "met_variables_summary.csv"),
    row.names = FALSE
  )
} else {
  cat("Met data file not found.\n")

  # Check ET calculation outputs
  et_file <- file.path(
    output_dir,
    "ET",
    "daily_water_balance_all_ET_methods_1997_present.csv"
  )
  if (file.exists(et_file)) {
    et_data <- read_csv(et_file, show_col_types = FALSE)
    et_vars <- names(et_data)
    cat("Water balance data columns:\n")
    print(et_vars)

    # Get date range
    if ("DATE" %in% names(et_data)) {
      et_data$Date <- as.Date(
        et_data$DATE,
        tryFormats = c("%Y-%m-%d", "%m/%d/%Y")
      )
      cat(
        "\nDate range:",
        min(et_data$Date, na.rm = TRUE),
        "to",
        max(et_data$Date, na.rm = TRUE),
        "\n"
      )
    }
  }
}

# -----------------------------------------------------------------------------
# CREATE COMPREHENSIVE SUMMARY TABLE
# -----------------------------------------------------------------------------
# Merge all information
comprehensive_summary <- site_info %>%
  left_join(
    site_metric_matrix %>% select(-site_name),
    by = "site_code"
  ) %>%
  left_join(
    date_ranges %>% select(site_code, hydro_start, hydro_end, hydro_n_days),
    by = "site_code"
  )


# -----------------------------------------------------------------------------
# SAVE OUTPUT TABLES
# -----------------------------------------------------------------------------
write.csv(
  site_metric_matrix,
  file.path(support_dir, "site_metric_availability.csv"),
  row.names = FALSE
)

write.csv(
  date_ranges,
  file.path(support_dir, "metric_date_ranges.csv"),
  row.names = FALSE
)

write.csv(
  comprehensive_summary,
  file.path(support_dir, "comprehensive_data_summary.csv"),
  row.names = FALSE
)

# Save metrics info
write.csv(
  metrics_info,
  file.path(support_dir, "storage_metrics_definitions.csv"),
  row.names = FALSE
)
