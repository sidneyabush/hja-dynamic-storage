# -----------------------------------------------------------------------------
# HJA Dynamic Storage - Configuration File
# -----------------------------------------------------------------------------
# This file contains all paths, constants, and shared settings used across
# the analysis workflow. Source this file at the beginning of each script.
#
# Usage: source("config.R") at the top of each analysis script
#
# For users running this code:
#   1. Set USE_LOCAL_DATA = TRUE if data files are in the repo's data/ folder
#   2. Set USE_LOCAL_DATA = FALSE and update BOX_BASE_DIR if using Box storage
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# PATH CONFIGURATION
# -----------------------------------------------------------------------------

# Set to TRUE if data files are in the local repo data/ folder
# Set to FALSE if data files are in Box cloud storage
USE_LOCAL_DATA <- FALSE

# Get the repo root directory (where this config.R file lives)
repo_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(repo_ofile) || !is.character(repo_ofile) || length(repo_ofile) == 0 || is.na(repo_ofile[1]) || repo_ofile[1] == "") {
  REPO_DIR <- getwd()
} else {
  REPO_DIR <- normalizePath(dirname(repo_ofile[1]), mustWork = FALSE)
}

if (USE_LOCAL_DATA) {
  # Local data paths (relative to repo root)
  BASE_DATA_DIR <- file.path(REPO_DIR, "data")
  OUTPUT_DIR <- file.path(REPO_DIR, "outputs")
} else {
  # Box cloud storage paths (update these for your system)
  BOX_BASE_DIR <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript"
  BASE_DATA_DIR <- file.path(BOX_BASE_DIR, "03_Data")
  OUTPUT_DIR <- file.path(BOX_BASE_DIR, "05_Outputs")
}

# Keep current workflow outputs grouped under one subfolder in 05_Outputs.
USE_CONSOLIDATED_OUTPUT_SUBDIR <- TRUE
CONSOLIDATED_OUTPUT_SUBDIR <- "final_workflow"

if (USE_CONSOLIDATED_OUTPUT_SUBDIR) {
  OUTPUT_DIR <- file.path(OUTPUT_DIR, CONSOLIDATED_OUTPUT_SUBDIR)
}

# Keep all figures under outputs for one unified tree.
FIGURES_DIR <- file.path(OUTPUT_DIR, "figs")

# Optional path overrides for sandboxed/local runs
BASE_DATA_DIR_OVERRIDE <- Sys.getenv("HJA_BASE_DATA_DIR", unset = "")
OUTPUT_DIR_OVERRIDE <- Sys.getenv("HJA_OUTPUT_DIR", unset = "")
FIGURES_DIR_OVERRIDE <- Sys.getenv("HJA_FIGURES_DIR", unset = "")

if (nzchar(BASE_DATA_DIR_OVERRIDE)) BASE_DATA_DIR <- BASE_DATA_DIR_OVERRIDE
if (nzchar(OUTPUT_DIR_OVERRIDE)) OUTPUT_DIR <- OUTPUT_DIR_OVERRIDE
if (nzchar(FIGURES_DIR_OVERRIDE)) FIGURES_DIR <- FIGURES_DIR_OVERRIDE
if (!nzchar(FIGURES_DIR_OVERRIDE)) FIGURES_DIR <- file.path(OUTPUT_DIR, "figs")

# Subdirectories
DISCHARGE_DIR <- file.path(BASE_DATA_DIR, "Q")
ET_DIR <- file.path(BASE_DATA_DIR, "DynamicStorage")
EC_DIR <- file.path(BASE_DATA_DIR, "EC")
ISOTOPE_DIR <- file.path(BASE_DATA_DIR, "Isotopes")
STREAM_TEMP_DIR <- file.path(BASE_DATA_DIR, "Stream_T")
MET_DIR <- file.path(BASE_DATA_DIR, "MET")

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Organized output directories
OUT_METRICS_DIR <- file.path(OUTPUT_DIR, "metrics")
OUT_MET_DYNAMIC_DIR <- file.path(OUT_METRICS_DIR, "dynamic")
OUT_MET_MOBILE_DIR <- file.path(OUT_METRICS_DIR, "mobile")
OUT_MET_EXTENDED_DIR <- file.path(OUT_METRICS_DIR, "extended_dynamic")
OUT_MET_ECO_DIR <- file.path(OUT_METRICS_DIR, "eco")
OUT_MET_SUPPORT_DIR <- file.path(OUT_METRICS_DIR, "support")

OUT_MASTER_DIR <- file.path(OUTPUT_DIR, "master")

OUT_STATS_DIR <- file.path(OUTPUT_DIR, "stats")
OUT_STATS_ANOVA_DIR <- file.path(OUT_STATS_DIR, "anova_tukey")
OUT_STATS_PCA_DIR <- file.path(OUT_STATS_DIR, "pca")
OUT_STATS_MLR_CATCH_CHARS_DIR <- file.path(OUT_STATS_DIR, "mlr_catch_chars_storage")
OUT_STATS_MLR_ECO_DIR <- file.path(OUT_STATS_DIR, "mlr_storage_eco")

OUT_TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
OUT_TABLES_MLR_DIR <- file.path(OUT_TABLES_DIR, "mlr")

for (d in c(
  OUT_METRICS_DIR,
  OUT_MET_DYNAMIC_DIR,
  OUT_MET_MOBILE_DIR,
  OUT_MET_EXTENDED_DIR,
  OUT_MET_ECO_DIR,
  OUT_MET_SUPPORT_DIR,
  OUT_MASTER_DIR,
  OUT_STATS_DIR,
  OUT_STATS_ANOVA_DIR,
  OUT_STATS_PCA_DIR,
  OUT_STATS_MLR_CATCH_CHARS_DIR,
  OUT_STATS_MLR_ECO_DIR,
  OUT_TABLES_DIR,
  OUT_TABLES_MLR_DIR
)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

# -----------------------------------------------------------------------------
# SITE DEFINITIONS
# -----------------------------------------------------------------------------

# Hydrometric sites with continuous streamflow data (for Dynamic metrics)
# Order: WS06 grouped with 07 and 08 for plotting
SITE_ORDER_HYDROMETRIC <- c(
  "WS09", # Watershed 09
  "WS10", # Watershed 10
  "WS01", # Watershed 01
  "Look", # Lookout Creek
  "WS02", # Watershed 02
  "WS03", # Watershed 03
  "WS06", # Watershed 06 (no chemistry/isotope data)
  "WS07", # Watershed 07
  "WS08", # Watershed 08
  "Mack" # Mack Creek
)

# Sites with chemistry data (for CHS) - excludes WS06, WS09, Look
SITE_ORDER_CHEMISTRY <- c(
  "WS10", # Watershed 10
  "WS01", # Watershed 01
  "WS02", # Watershed 02
  "WS03", # Watershed 03
  "WS06", # Watershed 06 (limited: 2017-2019)
  "WS07", # Watershed 07
  "WS08", # Watershed 08
  "Mack" # Mack Creek
)

# All sites including isotope-only sites (for MTT, Fyw, DR)
SITE_ORDER_ALL <- c(
  "WS09", # Watershed 09
  "WS10", # Watershed 10
  "WS01", # Watershed 01
  "Look", # Lookout Creek
  "WS02", # Watershed 02
  "WS03", # Watershed 03
  "MR", # McRae Creek (isotope-only)
  "WS06", # Watershed 06 (hydrometric only)
  "WS07", # Watershed 07
  "WS08", # Watershed 08
  "NC", # Nostoc Creek (isotope-only)
  "Mack", # Mack Creek
  "LC", # Longer Creek (isotope-only)
  "LO2", # Upper Lookout 2 (isotope-only)
  "CC", # Cold Creek (isotope-only)
  "LO1" # Upper Lookout 1 (isotope-only)
)

# Site name lookup table
SITE_NAMES <- c(
  "WS09" = "Watershed 09",
  "WS10" = "Watershed 10",
  "WS01" = "Watershed 01",
  "Look" = "Lookout",
  "WS02" = "Watershed 02",
  "WS03" = "Watershed 03",
  "WS06" = "Watershed 06",
  "WS07" = "Watershed 07",
  "WS08" = "Watershed 08",
  "Mack" = "Mack",
  "MR" = "McRae",
  "NC" = "Nostoc",
  "LC" = "Longer",
  "LO2" = "Upper Lookout 2",
  "CC" = "Cold Creek",
  "LO1" = "Upper Lookout 1"
)

# -----------------------------------------------------------------------------
# WATER YEAR RANGE
# -----------------------------------------------------------------------------

WY_START <- 1997
WY_END <- 2020

# -----------------------------------------------------------------------------
# STORAGE METRICS DEFINITIONS
# -----------------------------------------------------------------------------

# Dynamic storage metrics (from hydrometric data - annual time series)
# RBI = Richards-Baker Index, RCS = Recession Curve Slope
# FDC = Flow Duration Curve, SD = Storage-Discharge
DYNAMIC_METRICS <- c("RBI", "RCS", "FDC", "SD")

# Mobile storage metrics
# CHS = Chemical Hydrograph Separation (mean baseflow fraction)
# MTT = Mean Transit Time, Fyw = Young Water Fraction, DR = Damping Ratio
MOBILE_METRICS_ANNUAL <- c("CHS")
MOBILE_METRICS_SITE <- c("MTT", "Fyw", "DR") # Site-level from isotopes

# Extended dynamic storage metrics (from water balance - annual)
# WB = Water Balance (extended dynamic storage)
EXTENDED_DYNAMIC_METRICS <- c("WB")

# All storage metrics
ALL_STORAGE_METRICS <- c(
  DYNAMIC_METRICS,
  MOBILE_METRICS_ANNUAL,
  MOBILE_METRICS_SITE,
  EXTENDED_DYNAMIC_METRICS
)

# -----------------------------------------------------------------------------
# SHARED FILE NAMES AND CODE MAPPINGS
# -----------------------------------------------------------------------------

MASTER_ANNUAL_FILE <- "master_annual.csv"
MASTER_SITE_FILE <- "master_site.csv"
LEGACY_ANNUAL_FILE <- "HJA_Stor_Temp_Yr.csv"
LEGACY_SITE_FILE <- "HJA_Ave_StorageMetrics_CatCharacter.csv"

# Raw site codes that should be excluded from analysis tables.
SITE_EXCLUDE_RAW <- c("GSWSMA", "GSWSMF")

# Site recode map used when bringing met/discharge records into watershed codes.
SITECODE_RECODE_TO_GSMACK <- c("GSWSMC" = "GSMACK")

# Component watersheds used to make GSLOOK met composites.
GSLOOK_COMPOSITE_COMPONENT_SITES <- c("GSWS01", "GSWS06", "LONGER", "COLD")

# Legacy names that still appear in old exports.
LEGACY_STORAGE_RENAME_MAP <- c(
  recession_curve_slope = "RCS",
  fdc_slope = "FDC",
  S_annual_mm = "SD",
  mean_bf = "CHS",
  DS_sum = "WB"
)

# -----------------------------------------------------------------------------
# COLOR PALETTE
# -----------------------------------------------------------------------------

# Global plot text size (used across plotting scripts).
FIG_BASE_SIZE <- 18
FIG_AXIS_TEXT_SIZE <- 16
FIG_AXIS_TITLE_SIZE <- 18
FIG_STRIP_TEXT_SIZE <- 16
FIG_ANNOT_TEXT_SIZE <- 5
FIG_TILE_TEXT_SIZE <- 6
FIG_POINT_SIZE_SMALL <- 1.5
FIG_POINT_SIZE_MED <- 2.5
FIG_POINT_SIZE_LARGE <- 3.0

# 10-color palette for hydrometric sites (colorblind-friendly)
SITE_COLORS <- c(
  "WS09" = "#882255",
  "WS10" = "#AA4499",
  "WS01" = "#CC6677",
  "Look" = "#DDCC77",
  "WS02" = "#999933",
  "WS03" = "#117733",
  "WS06" = "#44AA99",
  "WS07" = "#88CCEE",
  "WS08" = "#6699CC",
  "Mack" = "#332288"
)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

#' Get water year from date
#' @param date A Date object
#' @return Integer water year (Oct-Dec = next year)
get_water_year <- function(date) {
  ifelse(
    lubridate::month(date) >= 10,
    lubridate::year(date) + 1,
    lubridate::year(date)
  )
}

#' Get water year day (1 = Oct 1)
#' @param date A Date object
#' @return Integer day of water year (1-366)
get_water_year_day <- function(date) {
  wy <- get_water_year(date)
  wy_start <- as.Date(paste0(wy - 1, "-10-01"))
  as.numeric(date - wy_start) + 1
}

#' Standardize site codes from various source files
#' @param site_code Character vector of site codes
#' @return Standardized site codes
standardize_site_code <- function(site_code) {
  site_code <- trimws(site_code)
  dplyr::case_when(
    # Named sites
    site_code == "GSWSMC" ~ "Mack",
    site_code == "GSMACK" ~ "Mack",
    site_code == "GSLOOK_FULL" ~ "Look",
    site_code == "GSLOOK" ~ "Look",
    site_code == "MCRAEC" ~ "MR",
    site_code == "NCCREC" ~ "NC",
    site_code == "LCCREC" ~ "LC",
    site_code == "LO2CRE" ~ "LO2",
    site_code == "CCCREE" ~ "CC",
    site_code == "LO1CRE" ~ "LO1",
    # Numbered watersheds: GSWS## -> WS##
    grepl("^GSWS[0-9]+$", site_code) ~ gsub("^GSWS", "WS", site_code),
    TRUE ~ site_code
  )
}

# Reverse mapping for legacy hydrometric exports (GSWS##/GSLOOK/GSWSMC style).
to_legacy_hydro_site_code <- function(site_code) {
  site_code <- trimws(site_code)
  dplyr::case_when(
    site_code == "Mack" ~ "GSWSMC",
    site_code == "Look" ~ "GSLOOK",
    grepl("^WS[0-9]+$", site_code) ~ gsub("^WS", "GSWS", site_code),
    TRUE ~ site_code
  )
}

rename_legacy_storage_metrics <- function(df) {
  hits <- names(LEGACY_STORAGE_RENAME_MAP)[names(LEGACY_STORAGE_RENAME_MAP) %in% names(df)]
  if (length(hits) == 0) {
    return(df)
  }
  dplyr::rename(df, !!!setNames(as.list(hits), unname(LEGACY_STORAGE_RENAME_MAP[hits])))
}

SITE_EXCLUDE_STANDARD <- unique(standardize_site_code(SITE_EXCLUDE_RAW))
