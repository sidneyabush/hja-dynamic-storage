# Assemble all annual/site-level metrics into master analysis tables.
# Inputs: dynamic_dir/RBI_RecessionCurve_Annual.csv; dynamic_dir/StorageDischarge_FDC_Annual.csv; mobile_dir/Annual_GW_Prop.csv; extended_dir/DS_drawdown_annual.csv; eco_dir/stream_thermal_lowflow_metrics_annual.csv; isotope_dir/MTT_FYW.csv; +1 more CSV files.
# Author: Sidney Bush
# Date: 2026-02-13

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(GGally)
library(ggcorrplot)

# Start clean
rm(list = ls())

# Load config (works from Rscript or source)
# Load project config
source("config.R")


theme_set(theme_pub(base_size = 12))

# Key input/output directories

base_dir    <- BASE_DATA_DIR
dynamic_dir <- OUT_MET_DYNAMIC_DIR
mobile_dir  <- OUT_MET_MOBILE_DIR
extended_dir <- OUT_MET_EXTENDED_DIR
eco_dir <- OUT_MET_ECO_DIR
master_dir <- file.path(OUTPUT_DIR, "master")
isotope_dir <- ISOTOPE_DIR
catchment_dir <- CATCHMENT_CHARACTERISTICS_DIR

# Make sure output directory exists
if (!dir.exists(master_dir)) dir.create(master_dir, recursive = TRUE)

# Helper to fail fast if a join key is not unique

assert_unique_keys <- function(df, keys, df_name) {
  dupes <- df %>%
    count(across(all_of(keys)), name = "n") %>%
    filter(n > 1)
  if (nrow(dupes) > 0) {
    stop(
      paste0(
        "Non-unique join keys detected in ", df_name, " for keys (",
        paste(keys, collapse = ", "), ")."
      )
    )
  }
}

# Dynamic storage metrics (annual)

# RBI and recession slope
rbi_path <- file.path(dynamic_dir, "rbi_rcs_annual.csv")
rbi_recession <- read_csv(
  rbi_path,
  show_col_types = FALSE
) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, RCS, RBI)
assert_unique_keys(rbi_recession, c("site", "year"), "rbi_recession")

# Storage-discharge, FDC, and supporting flow quantiles
fdc_path <- file.path(dynamic_dir, "storage_discharge_fdc_annual.csv")
storage_fdc <- read_csv(
  fdc_path,
  show_col_types = FALSE
) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, SD, FDC, Q99, Q50, Q01, Q5norm, CV_Q5norm)
assert_unique_keys(storage_fdc, c("site", "year"), "storage_fdc")

# Mobile storage metric (annual CHS)

# CHS = annual mean baseflow fraction
chs_path <- file.path(mobile_dir, "annual_gw_prop.csv")
baseflow <- read_csv(
  chs_path,
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(year = as.integer(year)) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  mutate(
    CHS = suppressWarnings(as.numeric(CHS)),
    CHS = ifelse(site %in% CHS_EXCLUDE_SITES, NA_real_, CHS)
  ) %>%
  select(site, year, CHS)
assert_unique_keys(baseflow, c("site", "year"), "baseflow")

# Extended dynamic storage metric (annual WB)

# Water-balance drawdown
wb_path <- file.path(extended_dir, "ds_drawdown_annual.csv")
wb_storage <- read_csv(
  wb_path,
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(year = as.integer(year)) %>%
  mutate(site = standardize_site_code(site)) %>%
  mutate(WB = suppressWarnings(as.numeric(WB))) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, WB)
assert_unique_keys(wb_storage, c("site", "year"), "wb_storage")

# Ecological response variables (annual)

# Thermal and low-flow response metrics
thermal_lowflow <- read_csv(
  file.path(eco_dir, "stream_thermal_lowflow_metrics_annual.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = if ("site" %in% names(.)) site else SITECODE,
    year = if ("year" %in% names(.)) year else wateryear
  ) %>%
  mutate(
    site = standardize_site_code(site),
    year = as.integer(year)
  ) %>%
  filter(year >= WY_START, year <= WY_END)

thermal_cols_required <- c(
  "T_7DMax", "Q_7Q5", "T_at_Q7Q5", "T_Q7Q5",
  "max_temp_7d_C",
  "q5_7d_mm_d", "temp_at_q5_7d_C", "temp_during_q5_7d_C",
  "min_Q_7d_mm_d", "temp_at_min_Q_7d_C", "temp_during_min_Q_7d_C",
  "P_WetSeason", "precip_nov_may_mm", "P_NovJan", "precip_nov_jan_mm", "Q5_CV"
)
thermal_cols_output <- c(
  "T_7DMax", "Q_7Q5", "T_at_Q7Q5", "T_Q7Q5",
  "max_temp_7d_C",
  "q5_7d_mm_d", "temp_at_q5_7d_C", "temp_during_q5_7d_C",
  "min_Q_7d_mm_d", "temp_at_min_Q_7d_C", "temp_during_min_Q_7d_C",
  "P_WetSeason", "precip_nov_may_mm", "Q5_CV"
)
for (nm in thermal_cols_required) {
  if (!(nm %in% names(thermal_lowflow))) {
    thermal_lowflow[[nm]] <- NA_real_
  }
}

thermal_lowflow <- thermal_lowflow %>%
  select(
    site, year, all_of(thermal_cols_required)
  ) %>%
  mutate(
    T_7DMax = ifelse(is.na(T_7DMax), max_temp_7d_C, T_7DMax),
    Q_7Q5 = ifelse(is.na(Q_7Q5), q5_7d_mm_d, Q_7Q5),
    T_at_Q7Q5 = ifelse(is.na(T_at_Q7Q5), temp_at_q5_7d_C, T_at_Q7Q5),
    T_Q7Q5 = ifelse(is.na(T_Q7Q5), temp_during_q5_7d_C, T_Q7Q5),
    P_WetSeason = dplyr::coalesce(P_WetSeason, precip_nov_may_mm, P_NovJan, precip_nov_jan_mm),
    precip_nov_may_mm = dplyr::coalesce(precip_nov_may_mm, P_WetSeason)
  ) %>%
  select(
    site, year, all_of(thermal_cols_output)
  )
assert_unique_keys(thermal_lowflow, c("site", "year"), "thermal_lowflow")

# Merge annual tables into one site-year master table

HJA_annual <- rbi_recession %>%
  full_join(storage_fdc, by = c("site", "year")) %>%
  full_join(baseflow, by = c("site", "year")) %>%
  full_join(wb_storage, by = c("site", "year")) %>%
  full_join(thermal_lowflow, by = c("site", "year")) %>%
  # Keep only hydrometric sites used in analysis
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  # Keep site order consistent across outputs/plots
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site, year)
assert_unique_keys(HJA_annual, c("site", "year"), "HJA_annual")

# Save annual site-year master table
write.csv(HJA_annual,
          file.path(master_dir, MASTER_ANNUAL_FILE),
          row.names = FALSE)

# Calculate site-level means across available years

HJA_avg <- HJA_annual %>%
  group_by(site) %>%
  summarise(
    n_years = n(),
    across(
      where(is.numeric),
      list(mean = ~mean(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# Add isotope-based site-level metrics (not annual)

# MTT1, MTT2, and Fyw
mtt_fyw <- read_csv(
  file.path(isotope_dir, "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  mutate(
    MTT1 = suppressWarnings(as.numeric(MTT1)),
    MTT2L_val = if ("MTT2L" %in% names(.)) suppressWarnings(as.numeric(MTT2L)) else NA_real_,
    MTT2H_val = if ("MTT2H" %in% names(.)) suppressWarnings(as.numeric(MTT2H)) else NA_real_,
    MTT2M_val = if ("MTT2M" %in% names(.)) suppressWarnings(as.numeric(MTT2M)) else NA_real_,
    MTT2 = suppressWarnings(as.numeric(dplyr::coalesce(
      MTT2M_val,
      rowMeans(cbind(MTT2L_val, MTT2H_val), na.rm = TRUE)
    ))),
    MTT2 = ifelse(is.nan(MTT2), NA_real_, MTT2),
    Fyw = suppressWarnings(as.numeric(FYWM))
  ) %>%
  select(site, MTT1, MTT2, Fyw) %>%
  filter(!is.na(site), site != "", site %in% SITE_ORDER_HYDROMETRIC)

# Damping ratio
damping_ratios <- read_csv(
  file.path(isotope_dir, "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  select(site, DR = DR_Overall, DR_err = DR__err) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

# Combine isotope metrics into one site-level table
isotope_metrics <- mtt_fyw %>%
  full_join(damping_ratios, by = "site")
assert_unique_keys(isotope_metrics, c("site"), "isotope_metrics")

# Join isotope metrics into site-average table
HJA_avg <- HJA_avg %>%
  left_join(isotope_metrics, by = "site")

# Add static catchment characteristics

catchment_chars <- read_csv(
  resolve_catchment_characteristics_file(),
  show_col_types = FALSE
) %>%
  rename(site = Site) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)
assert_unique_keys(catchment_chars, c("site"), "catchment_chars")

HJA_avg <- HJA_avg %>%
  left_join(catchment_chars, by = "site")

# Save final site-level master table
write.csv(HJA_avg,
          file.path(master_dir, MASTER_SITE_FILE),
          row.names = FALSE)

# Track sample sizes by metric/site
# Storage and response metric names used below
dynamic_metrics <- DYNAMIC_METRICS
mobile_metrics_annual <- MOBILE_METRICS_ANNUAL
mobile_metrics_site <- MOBILE_METRICS_SITE
extended_metrics <- EXTENDED_DYNAMIC_METRICS
response_metrics <- c("T_7DMax", "Q_7Q5", "T_Q7Q5", "P_WetSeason")

# Count available annual observations per site
annual_sample_sizes <- HJA_annual %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    n_years_total = n(),
    n_RBI = sum(!is.na(RBI)),
    n_RCS = sum(!is.na(RCS)),
    n_FDC = sum(!is.na(FDC)),
    n_SD = sum(!is.na(SD)),
    n_CHS = sum(!is.na(CHS)),
    n_WB = sum(!is.na(WB)),
    n_t_7dmax = sum(!is.na(T_7DMax)),
    n_q_7q5 = sum(!is.na(Q_7Q5)),
    n_p_wetseason = sum(!is.na(P_WetSeason)),
    .groups = "drop"
  )

# Add isotope availability flags
site_isotope <- HJA_avg %>%
  select(site, MTT1, MTT2, Fyw, DR) %>%
  mutate(
    has_MTT1 = !is.na(MTT1),
    has_MTT2 = !is.na(MTT2),
    has_Fyw = !is.na(Fyw),
    has_DR = !is.na(DR)
  ) %>%
  select(site, has_MTT1, has_MTT2, has_Fyw, has_DR)

sample_sizes <- annual_sample_sizes %>%
  left_join(site_isotope, by = "site")

# Save sample-size summary table (optional auxiliary output)
if (isTRUE(WRITE_AUX_OUTPUTS)) {
  write.csv(
    sample_sizes,
    file.path(master_dir, "metric_sample_sizes_by_site.csv"),
    row.names = FALSE
  )
}

# Build per-site summary statistics for all metrics
storage_metric_cols <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c(dynamic_metrics, mobile_metrics_annual, extended_metrics)
]
annual_metric_cols <- intersect(c(storage_metric_cols, response_metrics), names(HJA_annual))

annual_long <- HJA_annual %>%
  mutate(site = as.character(site)) %>%
  mutate(across(all_of(annual_metric_cols), ~suppressWarnings(as.numeric(.x)))) %>%
  select(site, all_of(annual_metric_cols)) %>%
  pivot_longer(
    cols = all_of(annual_metric_cols),
    names_to = "metric",
    values_to = "value"
  )

site_level_long <- isotope_metrics %>%
  select(site, any_of(mobile_metrics_site)) %>%
  pivot_longer(
    cols = any_of(mobile_metrics_site),
    names_to = "metric",
    values_to = "value"
  )

site_metric_summary_stats <- bind_rows(annual_long, site_level_long) %>%
  group_by(site, metric) %>%
  summarise(
    n = sum(is.finite(value)),
    mean = ifelse(n > 0, mean(value, na.rm = TRUE), NA_real_),
    sd = ifelse(n > 1, stats::sd(value, na.rm = TRUE), NA_real_),
    min = ifelse(n > 0, min(value, na.rm = TRUE), NA_real_),
    median = ifelse(n > 0, stats::median(value, na.rm = TRUE), NA_real_),
    max = ifelse(n > 0, max(value, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    metric_group = case_when(
      metric %in% dynamic_metrics ~ "dynamic",
      metric %in% mobile_metrics_annual ~ "mobile_annual",
      metric %in% mobile_metrics_site ~ "mobile_site",
      metric %in% extended_metrics ~ "extended_dynamic",
      metric %in% response_metrics ~ "ecological_response",
      TRUE ~ "other"
    ),
    site = factor(site, levels = SITE_ORDER_HYDROMETRIC)
  ) %>%
  arrange(site, metric) %>%
  mutate(site = as.character(site))

if (isTRUE(WRITE_AUX_OUTPUTS)) {
  write.csv(
    site_metric_summary_stats,
    file.path(master_dir, "master_site_metric_summary_stats.csv"),
    row.names = FALSE
  )
}

# ---- Data availability summary tables (moved from stats into metrics) ----

support_dir <- OUT_MET_SUPPORT_DIR
if (!dir.exists(support_dir)) {
  dir.create(support_dir, recursive = TRUE)
}

isotope_sites <- c(
  "WS09", "WS10", "WS01", "Look", "WS02", "WS03", "WS07", "WS08", "Mack"
)
site_info <- tibble(
  site = SITE_ORDER_HYDROMETRIC,
  site_name = unname(SITE_NAMES[SITE_ORDER_HYDROMETRIC]),
  hydrometric = TRUE,
  chemistry = SITE_ORDER_HYDROMETRIC %in% SITE_ORDER_CHEMISTRY,
  isotopes = SITE_ORDER_HYDROMETRIC %in% isotope_sites
)

metrics_info <- tribble(
  ~storage_type      , ~method                    , ~abbreviation , ~variable_name          , ~requires     ,
  "Dynamic"          , "Richards-Baker Index"     , "RBI"         , "RBI"                   , "hydrometric" ,
  "Dynamic"          , "Recession Curve Slope"    , "RCS"         , "RCS"                   , "hydrometric" ,
  "Dynamic"          , "Flow Duration Curve"      , "FDC"         , "FDC"                   , "hydrometric" ,
  "Dynamic"          , "Storage-Discharge"        , "SD"          , "SD"                    , "hydrometric" ,
  "Extended Dynamic" , "Water Balance"            , "WB"          , "WB"                    , "hydrometric" ,
  "Mobile"           , "Chemical Hydrograph Sep." , "CHS"         , "CHS"                   , "chemistry"   ,
  "Mobile"           , "Isotopic Damping Ratio"   , "DR"          , "DR"                    , "isotopes"    ,
  "Mobile"           , "Young Water Fraction"     , "Fyw"         , "Fyw"                   , "isotopes"    ,
  "Mobile"           , "Mean Transit Time (MTT1)" , "MTT1"        , "MTT1"                  , "isotopes"    ,
  "Mobile"           , "Mean Transit Time (MTT2)" , "MTT2"        , "MTT2"                  , "isotopes"
)

annual_metric_presence <- HJA_annual %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    RBI = any(is.finite(RBI)),
    RCS = any(is.finite(RCS)),
    FDC = any(is.finite(FDC)),
    SD = any(is.finite(SD)),
    CHS = any(is.finite(CHS)),
    WB = any(is.finite(WB)),
    .groups = "drop"
  )

site_metric_presence <- isotope_metrics %>%
  mutate(site = as.character(site)) %>%
  transmute(
    site,
    MTT1 = is.finite(MTT1),
    MTT2 = is.finite(MTT2),
    Fyw = is.finite(Fyw),
    DR = is.finite(DR)
  )

metric_cols <- STORAGE_METRIC_ORDER

site_metric_matrix <- site_info %>%
  select(site, site_name) %>%
  left_join(annual_metric_presence, by = "site") %>%
  left_join(site_metric_presence, by = "site") %>%
  mutate(across(all_of(metric_cols), ~ ifelse(is.na(.x), FALSE, .x))) %>%
  select(site, site_name, all_of(metric_cols))

metric_site_counts <- tibble(metric = metric_cols) %>%
  rowwise() %>%
  mutate(
    n_sites = sum(site_metric_matrix[[metric]], na.rm = TRUE),
    sites = paste(site_metric_matrix$site[site_metric_matrix[[metric]]], collapse = ", ")
  ) %>%
  ungroup()

discharge_avail <- read_csv(
  file.path(DISCHARGE_DIR, "HF00402_v14.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Date = as.Date(DATE, "%m/%d/%Y"),
    site = standardize_site_code(SITECODE)
  )

hydro_ranges <- discharge_avail %>%
  group_by(site) %>%
  summarise(
    hydro_start = min(Date, na.rm = TRUE),
    hydro_end = max(Date, na.rm = TRUE),
    hydro_n_days = sum(!is.na(MEAN_Q)),
    .groups = "drop"
  )

ec_candidates <- c(
  file.path(EC_DIR, "CF01201_v4.txt"),
  file.path(EC_DIR, "CF01201_v3.txt")
)
ec_file <- ec_candidates[file.exists(ec_candidates)][1]
if (!is.na(ec_file) && file.exists(ec_file)) {
  ec_data <- read_delim(ec_file, delim = "\t", show_col_types = FALSE)
  if ("DATE" %in% names(ec_data) && "SITECODE" %in% names(ec_data)) {
    chem_ranges <- ec_data %>%
      mutate(Date = as.Date(DATE), site = standardize_site_code(SITECODE)) %>%
      group_by(site) %>%
      summarise(
        chem_start = min(Date, na.rm = TRUE),
        chem_end = max(Date, na.rm = TRUE),
        chem_n_days = n(),
        .groups = "drop"
      )
  } else {
    chem_ranges <- tibble(site = character())
  }
} else {
  chem_ranges <- tibble(site = character())
}

wb_daily_file <- resolve_water_balance_daily_file()
wb_data <- read_csv(wb_daily_file, show_col_types = FALSE) %>%
  mutate(
    Date = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y")),
    site = standardize_site_code(SITECODE)
  )
wb_ranges <- wb_data %>%
  group_by(site) %>%
  summarise(
    wb_start = min(Date, na.rm = TRUE),
    wb_end = max(Date, na.rm = TRUE),
    wb_n_days = sum(!is.na(P_mm_d) & !is.na(Q_mm_d) & !is.na(ET_mm_d)),
    .groups = "drop"
  )

date_ranges <- site_info %>%
  select(site, site_name) %>%
  left_join(hydro_ranges, by = "site") %>%
  left_join(chem_ranges, by = "site") %>%
  left_join(wb_ranges, by = "site")

met_file <- file.path(OUT_MET_SUPPORT_DIR, "watersheds_met_q.csv")
if (file.exists(met_file)) {
  met_data <- read_csv(met_file, show_col_types = FALSE)
  exclude_cols <- c("DATE", "Date", "date", "SITECODE", "site", "WATERYEAR", "wateryear")
  met_vars <- setdiff(names(met_data), exclude_cols)
  if ("DATE" %in% names(met_data)) {
    met_data$Date <- as.Date(met_data$DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y"))
  }
  met_summary <- data.frame(
    variable = met_vars,
    n_obs = sapply(met_vars, function(v) sum(!is.na(met_data[[v]]))),
    pct_complete = sapply(met_vars, function(v) round(100 * sum(!is.na(met_data[[v]])) / nrow(met_data), 1)),
    min_date = sapply(met_vars, function(v) {
      valid <- !is.na(met_data[[v]])
      if (any(valid)) as.character(min(met_data$Date[valid], na.rm = TRUE)) else NA_character_
    }),
    max_date = sapply(met_vars, function(v) {
      valid <- !is.na(met_data[[v]])
      if (any(valid)) as.character(max(met_data$Date[valid], na.rm = TRUE)) else NA_character_
    }),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(n_obs))
  if (isTRUE(WRITE_AUX_OUTPUTS)) {
    write.csv(
      met_summary,
      file.path(support_dir, "met_variables_summary.csv"),
      row.names = FALSE
    )
  }
}

comprehensive_summary <- site_info %>%
  left_join(site_metric_matrix %>% select(-site_name), by = "site") %>%
  left_join(date_ranges %>% select(site, hydro_start, hydro_end, hydro_n_days), by = "site")

eco_response_availability <- HJA_annual %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    n_wy_total = n_distinct(year),
    n_wy_T_7DMax = sum(is.finite(T_7DMax)),
    n_wy_Q_7Q5 = sum(is.finite(Q_7Q5)),
    n_wy_T_Q7Q5 = sum(is.finite(T_Q7Q5)),
    n_wy_P_WetSeason = sum(is.finite(P_WetSeason)),
    .groups = "drop"
  ) %>%
  right_join(tibble(site = SITE_ORDER_HYDROMETRIC), by = "site") %>%
  mutate(
    n_wy_total = ifelse(is.na(n_wy_total), 0L, n_wy_total),
    n_wy_T_7DMax = ifelse(is.na(n_wy_T_7DMax), 0L, n_wy_T_7DMax),
    n_wy_Q_7Q5 = ifelse(is.na(n_wy_Q_7Q5), 0L, n_wy_Q_7Q5),
    n_wy_T_Q7Q5 = ifelse(is.na(n_wy_T_Q7Q5), 0L, n_wy_T_Q7Q5),
    n_wy_P_WetSeason = ifelse(is.na(n_wy_P_WetSeason), 0L, n_wy_P_WetSeason),
    missing_reason = case_when(
      site == "WS09" & n_wy_T_7DMax == 0 ~ "No WS09 stream-temperature records in HT00451_v10.txt",
      n_wy_T_7DMax == 0 ~ "No stream-temperature WY records",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(factor(site, levels = SITE_ORDER_HYDROMETRIC))

eco_response_wy_coverage <- HJA_annual %>%
  mutate(site = as.character(site)) %>%
  select(site, year, T_7DMax, Q_7Q5, T_Q7Q5, P_WetSeason) %>%
  pivot_longer(
    cols = c(T_7DMax, Q_7Q5, T_Q7Q5, P_WetSeason),
    names_to = "response",
    values_to = "value"
  ) %>%
  group_by(site, response) %>%
  summarise(
    n_wy_with_data = sum(is.finite(value)),
    first_wy_with_data = ifelse(any(is.finite(value)), min(year[is.finite(value)]), NA_integer_),
    last_wy_with_data = ifelse(any(is.finite(value)), max(year[is.finite(value)]), NA_integer_),
    .groups = "drop"
  ) %>%
  arrange(factor(site, levels = SITE_ORDER_HYDROMETRIC), response)

ht004_temp_file <- file.path(STREAM_TEMP_DIR, "HT00451_v10.txt")
if (file.exists(ht004_temp_file)) {
  ht004_temp_raw <- read_csv(ht004_temp_file, show_col_types = FALSE) %>%
    mutate(
      site_raw = as.character(SITECODE),
      site = standardize_site_code(site_raw),
      datetime = as.POSIXct(DATE_TIME, tz = "UTC"),
      temp_val = WATERTEMP_MEAN
    )
  stream_temp_source_coverage <- ht004_temp_raw %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
    group_by(site_raw, site) %>%
    summarise(
      n_records = n(),
      n_non_na_temp = sum(is.finite(temp_val), na.rm = TRUE),
      first_datetime = min(datetime, na.rm = TRUE),
      last_datetime = max(datetime, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor(site, levels = SITE_ORDER_HYDROMETRIC))
  if (!("WS09" %in% stream_temp_source_coverage$site)) {
    stream_temp_source_coverage <- bind_rows(
      stream_temp_source_coverage,
      tibble(
        site_raw = "GSWS09",
        site = "WS09",
        n_records = 0L,
        n_non_na_temp = 0L,
        first_datetime = as.POSIXct(NA),
        last_datetime = as.POSIXct(NA)
      )
    ) %>%
      arrange(factor(site, levels = SITE_ORDER_HYDROMETRIC))
  }
} else {
  stream_temp_source_coverage <- tibble(
    site_raw = character(),
    site = character(),
    n_records = integer(),
    n_non_na_temp = integer(),
    first_datetime = as.POSIXct(character()),
    last_datetime = as.POSIXct(character())
  )
}

if (isTRUE(WRITE_AUX_OUTPUTS)) {
  write.csv(
    site_metric_matrix,
    file.path(support_dir, "site_metric_availability.csv"),
    row.names = FALSE
  )
  write.csv(
    metric_site_counts,
    file.path(support_dir, "metric_site_counts.csv"),
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
  write.csv(
    metrics_info,
    file.path(support_dir, "storage_metrics_definitions.csv"),
    row.names = FALSE
  )
  write.csv(
    eco_response_availability,
    file.path(support_dir, "eco_response_availability_by_site.csv"),
    row.names = FALSE
  )
  write.csv(
    eco_response_wy_coverage,
    file.path(support_dir, "eco_response_wy_coverage.csv"),
    row.names = FALSE
  )
  write.csv(
    stream_temp_source_coverage,
    file.path(support_dir, "stream_temp_source_coverage_ht00451.csv"),
    row.names = FALSE
  )
}
