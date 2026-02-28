# hja dynamic storage - configuration file.
# author: sidney bush
# date: 2026-02-13

# default project paths (single supported layout)
REPO_DIR <- normalizePath(getwd(), mustWork = FALSE)
BOX_BASE_DIR <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript"
FINAL_WORKFLOW_ROOT <- file.path(BOX_BASE_DIR, "final_workflow")
BASE_DATA_DIR <- file.path(FINAL_WORKFLOW_ROOT, "inputs")
OUTPUT_DIR <- file.path(FINAL_WORKFLOW_ROOT, "outputs")
MS_MATERIALS_DIR <- file.path(FINAL_WORKFLOW_ROOT, "ms_materials")
EXPLORATORY_PLOTS_DIR <- file.path(OUTPUT_DIR, "exploratory_plots")
CONCEPTUAL_DIAGRAM_DIR <- file.path(EXPLORATORY_PLOTS_DIR, "conceptual_diagram")
SUPP_EXPLORATORY_DIR <- file.path(EXPLORATORY_PLOTS_DIR, "supp")
SUPP_EXPLORATORY_PDF_DIR <- SUPP_EXPLORATORY_DIR
# alias names still used in some scripts
SUPP_LEGACY_DIR <- SUPP_EXPLORATORY_DIR
SUPP_LEGACY_PDF_DIR <- SUPP_EXPLORATORY_PDF_DIR

# manuscript-facing outputs.
# figures and tables now write directly to main/supp.
FIGURES_DIR <- MS_MATERIALS_DIR
MS_MAIN_DIR <- file.path(MS_MATERIALS_DIR, "main")
MS_SUPP_DIR <- file.path(MS_MATERIALS_DIR, "supp")
MS_FIG_MAIN_DIR <- MS_MAIN_DIR
MS_FIG_SUPP_DIR <- MS_SUPP_DIR
MS_FIG_MAIN_PDF_DIR <- file.path(MS_MAIN_DIR, "pdf")
MS_FIG_SUPP_PDF_DIR <- file.path(MS_SUPP_DIR, "pdf")
MS_TABLES_MAIN_DIR <- MS_MAIN_DIR
MS_TABLES_SUPP_DIR <- MS_SUPP_DIR

# input subdirectories
# default input layout keeps small input groups at the inputs root.
DISCHARGE_DIR <- BASE_DATA_DIR
EC_DIR <- BASE_DATA_DIR
ISOTOPE_DIR <- BASE_DATA_DIR
STREAM_TEMP_DIR <- BASE_DATA_DIR
MET_DIR <- file.path(BASE_DATA_DIR, "all_hydromet")
CATCHMENT_CHARACTERISTICS_DIR <- BASE_DATA_DIR

EXPLORATORY_ET_METHODS_DIR <- file.path(EXPLORATORY_PLOTS_DIR, "et_methods")

# output controls
# keep only workflow-essential csvs by default; set true when you need
# manuscript/helper table exports.
WRITE_TABLE_OUTPUTS <- FALSE
WRITE_AUX_OUTPUTS <- FALSE

# create output directory if it doesn't exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# organized output directories
OUT_METRICS_DIR <- file.path(OUTPUT_DIR, "metrics")
OUT_MET_DYNAMIC_DIR <- file.path(OUT_METRICS_DIR, "dynamic")
OUT_MET_MOBILE_DIR <- file.path(OUT_METRICS_DIR, "mobile")
OUT_MET_EXTENDED_DIR <- file.path(OUT_METRICS_DIR, "extended_dynamic")
OUT_MET_ECO_DIR <- file.path(OUT_METRICS_DIR, "eco")
OUT_MET_SUPPORT_DIR <- file.path(OUT_METRICS_DIR, "support")

OUT_STATS_DIR <- file.path(OUTPUT_DIR, "models")
OUT_STATS_ANOVA_DIR <- file.path(OUT_STATS_DIR, "anova_tukey")
OUT_STATS_PCA_DIR <- file.path(OUT_STATS_DIR, "pca")
OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR <- file.path(
  OUT_STATS_DIR,
  "catchment_char_storage_mlr"
)
# alias name still used in some scripts
OUT_MODELS_WATERSHED_CHAR_STORAGE_MLR_DIR <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR <- file.path(
  OUT_STATS_DIR,
  "storage_eco_response_mlr"
)
# alias name still used in some scripts
OUT_MODELS_STORAGE_ECOVAR_MLR_DIR <- OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR
OUT_STATS_VALIDATION_DIR <- file.path(OUT_STATS_DIR, "validation")

OUT_TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
OUT_TABLES_MLR_DIR <- file.path(OUT_TABLES_DIR, "mlr")

output_dirs <- c(
  OUT_METRICS_DIR,
  OUT_MET_DYNAMIC_DIR,
  OUT_MET_MOBILE_DIR,
  OUT_MET_EXTENDED_DIR,
  OUT_MET_ECO_DIR,
  OUT_MET_SUPPORT_DIR,
  file.path(OUTPUT_DIR, "master"),
  MS_MATERIALS_DIR,
  MS_MAIN_DIR,
  MS_SUPP_DIR,
  MS_FIG_MAIN_DIR,
  MS_FIG_SUPP_DIR,
  MS_FIG_MAIN_PDF_DIR,
  MS_FIG_SUPP_PDF_DIR,
  MS_TABLES_MAIN_DIR,
  MS_TABLES_SUPP_DIR,
  OUT_STATS_DIR,
  OUT_STATS_ANOVA_DIR,
  OUT_STATS_PCA_DIR,
  OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR,
  OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
  OUT_STATS_VALIDATION_DIR,
  EXPLORATORY_PLOTS_DIR,
  CONCEPTUAL_DIAGRAM_DIR,
  SUPP_LEGACY_DIR,
  EXPLORATORY_ET_METHODS_DIR
)

if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  output_dirs <- c(output_dirs, OUT_TABLES_DIR, OUT_TABLES_MLR_DIR)
}

for (d in output_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# site definitions

# hydrometric sites with continuous streamflow data (for dynamic metrics)
# order: ws06 grouped with 07 and 08 for plotting
SITE_ORDER_HYDROMETRIC <- c(
  "WS09", # Catchment 09
  "WS10", # Catchment 10
  "WS01", # Catchment 01
  "Look", # Lookout Creek
  "WS02", # Catchment 02
  "WS03", # Catchment 03
  "Mack", # Mack Creek
  "WS06", # Catchment 06 (no chemistry/isotope data)
  "WS07", # Catchment 07
  "WS08" # Catchment 08
)

# sites with chemistry data (for chs)
SITE_ORDER_CHEMISTRY <- c(
  "WS10", # Catchment 10
  "WS01", # Catchment 01
  "Look", # Lookout Creek
  "WS02", # Catchment 02
  "WS03", # Catchment 03
  "WS06", # Catchment 06
  "WS07", # Catchment 07
  "WS08", # Catchment 08
  "Mack" # Mack Creek
)

# active analysis/plotting site set
SITE_ORDER_ALL <- SITE_ORDER_HYDROMETRIC

# site name lookup table
SITE_NAMES <- c(
  "WS09" = "Catchment 09",
  "WS10" = "Catchment 10",
  "WS01" = "Catchment 01",
  "Look" = "Lookout",
  "WS02" = "Catchment 02",
  "WS03" = "Catchment 03",
  "WS06" = "Catchment 06",
  "WS07" = "Catchment 07",
  "WS08" = "Catchment 08",
  "Mack" = "Mack",
  "MR" = "McRae",
  "NC" = "Nostoc",
  "LC" = "Longer",
  "LO2" = "Upper Lookout 2",
  "CC" = "Cold Creek",
  "LO1" = "Upper Lookout 1"
)

# water year range

WY_START <- 1997
WY_END <- 2020

# minimum number of daily ec+q observations needed to keep a chs water year.
CHS_MIN_DAYS_PER_WY <- 300

# sites intentionally excluded from chs-based modeling summaries.
# note: gslook is standardized to look upstream.
CHS_EXCLUDE_SITES <- c("Look", "WS09")

# storage metrics definitions

# default storage-metric order across figures/tables:
# rbi, rcs, fdc, sd, wb, chs, dr, fyw, mtt1, mtt2
STORAGE_METRIC_ORDER <- c("RBI", "RCS", "FDC", "SD", "WB", "CHS", "DR", "Fyw", "MTT1", "MTT2")

# dynamic storage metrics (from streamflow data - year-by-year records)
# rbi = richards-baker index, rcs = recession curve slope
# fdc = flow duration curve, sd = storage-discharge
DYNAMIC_METRICS <- STORAGE_METRIC_ORDER[STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD")]

# mobile storage metrics
# chs = chemical hydrograph separation (mean baseflow fraction)
# mtt1/mtt2 = mean transit time period-specific values, fyw = young water fraction, dr = damping ratio
MOBILE_METRICS_ANNUAL <- STORAGE_METRIC_ORDER[STORAGE_METRIC_ORDER %in% c("CHS")]
MOBILE_METRICS_SITE <- STORAGE_METRIC_ORDER[STORAGE_METRIC_ORDER %in% c("DR", "Fyw", "MTT1", "MTT2")] # Site-level from isotopes

# extended dynamic storage metrics (from water balance - annual)
# wb = water balance (extended dynamic storage)
EXTENDED_DYNAMIC_METRICS <- STORAGE_METRIC_ORDER[STORAGE_METRIC_ORDER %in% c("WB")]

# shared metric display order for plotting.
# keep isotope mobile metrics split as mtt1 and mtt2.
PLOT_ORDER_DYNAMIC_STORAGE <- STORAGE_METRIC_ORDER[STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD", "WB")]
PLOT_ORDER_MOBILE_STORAGE <- STORAGE_METRIC_ORDER[STORAGE_METRIC_ORDER %in% c("CHS", "DR", "Fyw", "MTT1", "MTT2")]
PLOT_MOBILE_STORAGE_SITE_COLS <- c(
  "CHS" = "CHS_mean",
  "DR" = "DR",
  "Fyw" = "Fyw",
  "MTT1" = "MTT1",
  "MTT2" = "MTT2"
)

# all storage metrics
ALL_STORAGE_METRICS <- STORAGE_METRIC_ORDER

# shared file names and code mappings

MASTER_ANNUAL_FILE <- "master_annual.csv"
MASTER_SITE_FILE <- "master_site.csv"

# raw site codes that should be excluded from analysis tables.
SITE_EXCLUDE_RAW <- c("GSWSMA", "GSWSMF")

# site recode map used when bringing met/discharge records into catchment codes.
SITECODE_RECODE_TO_GSMACK <- c("GSWSMC" = "GSMACK")

# component catchments used to make gslook met composites.
GSLOOK_COMPOSITE_COMPONENT_SITES <- c("GSWS01", "GSWS06", "LONGER", "COLD")


# color palette

# global plot text size (used across plotting scripts).
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

# global label/annotation behavior.
# use these in plotting scripts so labels are consistently readable.
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

# helper functions

# publication plot theme: no title/subtitle and no grid lines.
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
    # named sites
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
    # numbered catchments: gsws## -> ws##
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

CATCHMENT_PREDICTOR_LABELS <- c(
  "basin_slope" = "Basin Slope",
  "Harvest" = "Harvest",
  "Landslide_Total" = "Total Landslide",
  "Landslide_Young" = "Young Landslide",
  "Lava 1 (%)" = "Lava-1 (%)",
  "Lava1_per" = "Lava-1 (%)",
  "Lava 2 (%)" = "Lava-2 (%)",
  "Lava2_per" = "Lava-2 (%)",
  "Ash_Per" = "Ash (%)",
  "Pyro_per" = "Pyroclastic (%)"
)

label_catchment_predictor <- function(x) {
  x_chr <- as.character(x)
  out <- unname(CATCHMENT_PREDICTOR_LABELS[x_chr])
  miss <- is.na(out)
  out[miss] <- x_chr[miss]
  out
}

label_catchment_predictor_list <- function(x, sep = ";") {
  x_chr <- as.character(x)
  vapply(
    x_chr,
    function(val) {
      if (is.na(val) || !nzchar(val)) return(val)
      parts <- trimws(strsplit(val, sep, fixed = TRUE)[[1]])
      parts <- parts[nzchar(parts)]
      if (length(parts) == 0) return(val)
      paste(label_catchment_predictor(parts), collapse = "; ")
    },
    FUN.VALUE = character(1)
  )
}

SITE_EXCLUDE_STANDARD <- unique(standardize_site_code(SITE_EXCLUDE_RAW))

# water-balance daily input file
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
  default_file <- file.path(CATCHMENT_CHARACTERISTICS_DIR, "catchment_char.csv")
  if (!file.exists(default_file)) {
    stop(
      paste0(
        "Missing default catchment characteristics file: ",
        default_file
      )
    )
  }
  default_file
}
