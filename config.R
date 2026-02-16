# HJA Dynamic Storage - Configuration File.
# Author: Sidney Bush
# Date: 2026-02-13

# Default project paths (single supported layout)
REPO_DIR <- normalizePath(getwd(), mustWork = FALSE)
BOX_BASE_DIR <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript"
FINAL_WORKFLOW_ROOT <- file.path(BOX_BASE_DIR, "final_workflow")
BASE_DATA_DIR <- file.path(FINAL_WORKFLOW_ROOT, "inputs")
OUTPUT_DIR <- file.path(FINAL_WORKFLOW_ROOT, "outputs")
EXPLORATORY_PLOTS_DIR <- file.path(OUTPUT_DIR, "exploratory_plots")

# Keep all figures under outputs.
FIGURES_DIR <- file.path(OUTPUT_DIR, "figs")

# Input subdirectories
DISCHARGE_DIR <- file.path(BASE_DATA_DIR, "q")
EC_DIR <- file.path(BASE_DATA_DIR, "ec")
ISOTOPE_DIR <- file.path(BASE_DATA_DIR, "isotopes")
STREAM_TEMP_DIR <- file.path(BASE_DATA_DIR, "stream_t")
MET_DIR <- file.path(BASE_DATA_DIR, "all_hydromet")
CATCHMENT_CHARACTERISTICS_DIR <- file.path(
  BASE_DATA_DIR,
  "catchment_characteristics"
)

EXPLORATORY_ET_METHODS_DIR <- file.path(EXPLORATORY_PLOTS_DIR, "et_methods")

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

OUT_STATS_DIR <- file.path(OUTPUT_DIR, "models")
OUT_STATS_ANOVA_DIR <- file.path(OUT_STATS_DIR, "anova_tukey")
OUT_STATS_PCA_DIR <- file.path(OUT_STATS_DIR, "pca")
OUT_MODELS_WATERSHED_CHAR_STORAGE_MLR_DIR <- file.path(
  OUT_STATS_DIR,
  "watershed_char_storage_mlr"
)
OUT_MODELS_STORAGE_ECOVAR_MLR_DIR <- file.path(
  OUT_STATS_DIR,
  "storage_ecovar_mlr"
)
OUT_STATS_VALIDATION_DIR <- file.path(OUT_STATS_DIR, "validation")

OUT_TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
OUT_TABLES_MLR_DIR <- file.path(OUT_TABLES_DIR, "mlr")

for (d in c(
  OUT_METRICS_DIR,
  OUT_MET_DYNAMIC_DIR,
  OUT_MET_MOBILE_DIR,
  OUT_MET_EXTENDED_DIR,
  OUT_MET_ECO_DIR,
  OUT_MET_SUPPORT_DIR,
  file.path(OUTPUT_DIR, "master"),
  OUT_STATS_DIR,
  OUT_STATS_ANOVA_DIR,
  OUT_STATS_PCA_DIR,
  OUT_MODELS_WATERSHED_CHAR_STORAGE_MLR_DIR,
  OUT_MODELS_STORAGE_ECOVAR_MLR_DIR,
  OUT_STATS_VALIDATION_DIR,
  OUT_TABLES_DIR,
  OUT_TABLES_MLR_DIR,
  EXPLORATORY_PLOTS_DIR,
  EXPLORATORY_ET_METHODS_DIR
)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

# SITE DEFINITIONS

# Hydrometric sites with continuous streamflow data (for Dynamic metrics)
# Order: WS06 grouped with 07 and 08 for plotting
SITE_ORDER_HYDROMETRIC <- c(
  "WS09", # Watershed 09
  "WS10", # Watershed 10
  "WS01", # Watershed 01
  "Look", # Lookout Creek
  "WS02", # Watershed 02
  "WS03", # Watershed 03
  "Mack", # Mack Creek
  "WS06", # Watershed 06 (no chemistry/isotope data)
  "WS07", # Watershed 07
  "WS08" # Watershed 08
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

# Active analysis/plotting site set
SITE_ORDER_ALL <- SITE_ORDER_HYDROMETRIC

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

# WATER YEAR RANGE

WY_START <- 1997
WY_END <- 2020

# STORAGE METRICS DEFINITIONS

# Dynamic storage metrics (from streamflow data - year-by-year records)
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

# SHARED FILE NAMES AND CODE MAPPINGS

MASTER_ANNUAL_FILE <- "master_annual.csv"
MASTER_SITE_FILE <- "master_site.csv"

# Raw site codes that should be excluded from analysis tables.
SITE_EXCLUDE_RAW <- c("GSWSMA", "GSWSMF")

# Site recode map used when bringing met/discharge records into watershed codes.
SITECODE_RECODE_TO_GSMACK <- c("GSWSMC" = "GSMACK")

# Component watersheds used to make GSLOOK met composites.
GSLOOK_COMPOSITE_COMPONENT_SITES <- c("GSWS01", "GSWS06", "LONGER", "COLD")


# COLOR PALETTE

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
FIG_WIDTH_SCALE <- 1.35
FIG_HEIGHT_SCALE <- 1.35

# Global label/annotation behavior.
# Use these in plotting scripts so labels are consistently readable.
FIG_LABEL_CHECK_OVERLAP <- TRUE
FIG_LABEL_CLIP <- "off" # "off" prevents annotation clipping at panel bounds
FIG_LABEL_PLOT_MARGIN_PT <- 18 # extra margin so outer labels are not cut
FIG_MEAN_LINE_LINETYPES <- c(
  "solid",
  "dashed",
  "dotdash",
  "longdash",
  "twodash",
  "dotted",
  "22",
  "42",
  "F2"
)
FIG_MEAN_LABEL_DIGITS <- 2

# 10-color palette for streamflow sites (colorblind-friendly)
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

# HELPER FUNCTIONS

# Publication plot theme: no title/subtitle and no grid lines.
theme_pub <- function(base_size = FIG_BASE_SIZE) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(hjust = 0),
      strip.text.y = ggplot2::element_text(hjust = 0)
    )
}

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

make_panel_label_map <- function(values) {
  values <- as.character(values)
  if (length(values) == 0) {
    return(character())
  }
  idx <- seq_along(values)
  labels <- paste0(letters[((idx - 1) %% 26) + 1], ") ", values)
  stats::setNames(labels, values)
}

SITE_EXCLUDE_STANDARD <- unique(standardize_site_code(SITE_EXCLUDE_RAW))

# Water-balance daily input (current workflow path only).
resolve_water_balance_daily_file <- function() {
  path <- file.path(
    OUT_MET_SUPPORT_DIR,
    "daily_water_balance_et_hamon_zhang_coeff_interp.csv"
  )
  if (!file.exists(path)) {
    stop(
      paste0(
        "Missing required water-balance daily file: ",
        path
      )
    )
  }
  path
}

resolve_drainage_area_file <- function() {
  path <- file.path(CATCHMENT_CHARACTERISTICS_DIR, "drainage_area.csv")
  if (!file.exists(path)) {
    stop(
      paste0(
        "Missing required drainage area file: ",
        path
      )
    )
  }
  path
}

resolve_catchment_characteristics_file <- function() {
  canonical <- file.path(CATCHMENT_CHARACTERISTICS_DIR, "catchment_char.csv")
  if (!file.exists(canonical)) {
    stop(
      paste0(
        "Missing canonical catchment characteristics file: ",
        canonical
      )
    )
  }
  canonical
}
