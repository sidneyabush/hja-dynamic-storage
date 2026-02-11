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
REPO_DIR <- normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE)
if (is.na(REPO_DIR) || REPO_DIR == "") {
  # Fallback if sourced interactively
  REPO_DIR <- getwd()
}

if (USE_LOCAL_DATA) {
  # Local data paths (relative to repo root)
  BASE_DATA_DIR <- file.path(REPO_DIR, "data")
  OUTPUT_DIR <- file.path(REPO_DIR, "outputs")
  FIGURES_DIR <- file.path(REPO_DIR, "figures")
} else {
  # Box cloud storage paths (update these for your system)
  BOX_BASE_DIR <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript"
  BASE_DATA_DIR <- file.path(BOX_BASE_DIR, "03_Data")
  OUTPUT_DIR <- file.path(BOX_BASE_DIR, "05_Outputs")
  FIGURES_DIR <- file.path(BOX_BASE_DIR, "Final_Workflow", "03_Figures")
}

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
# COLOR PALETTE
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# PRINT CONFIGURATION SUMMARY
# -----------------------------------------------------------------------------

cat(
  "
================================================================================
HJA Dynamic Storage - Configuration Loaded
================================================================================
"
)
cat(
  "Data source:",
  ifelse(USE_LOCAL_DATA, "Local (data/)", "Box cloud storage"),
  "\n"
)
cat("Data directory:", BASE_DATA_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Water year range:", WY_START, "-", WY_END, "\n")
cat("Hydrometric sites:", length(SITE_ORDER_HYDROMETRIC), "\n")
cat(
  "================================================================================\n\n"
)
