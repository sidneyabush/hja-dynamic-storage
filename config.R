# settings used across the HJA dynamic storage analysis
# author: Sidney Bush
# date: 2026-02-13

# default paths
# by default the code looks for inputs/ and writes outputs inside the repo
# you can still override these with environment variables
REPO_DIR <- normalizePath(
  Sys.getenv("HJA_REPO_DIR", unset = getwd()),
  mustWork = FALSE
)
default_final_workflow_root <- REPO_DIR
FINAL_WORKFLOW_ROOT <- Sys.getenv(
  "HJA_FINAL_WORKFLOW_ROOT",
  unset = default_final_workflow_root
)
BASE_DATA_DIR <- Sys.getenv(
  "HJA_BASE_DATA_DIR",
  unset = file.path(FINAL_WORKFLOW_ROOT, "inputs")
)
OUTPUT_DIR <- Sys.getenv(
  "HJA_OUTPUT_DIR",
  unset = file.path(FINAL_WORKFLOW_ROOT, "outputs")
)
MS_MATERIALS_DIR <- Sys.getenv(
  "HJA_MS_MATERIALS_DIR",
  unset = file.path(FINAL_WORKFLOW_ROOT, "ms_materials")
)
EXPLORATORY_PLOTS_DIR <- file.path(OUTPUT_DIR, "exploratory_plots")
CONCEPTUAL_DIAGRAM_DIR <- file.path(EXPLORATORY_PLOTS_DIR, "conceptual_diagram")
SUPP_EXPLORATORY_DIR <- file.path(EXPLORATORY_PLOTS_DIR, "supp")
SUPP_EXPLORATORY_PDF_DIR <- SUPP_EXPLORATORY_DIR
# keep older names that still appear in some scripts
SUPP_LEGACY_DIR <- SUPP_EXPLORATORY_DIR
SUPP_LEGACY_PDF_DIR <- SUPP_EXPLORATORY_PDF_DIR

# paper figures and tables
# write figures and tables straight to main/ and supp/
FIGURES_DIR <- MS_MATERIALS_DIR
MS_MAIN_DIR <- file.path(MS_MATERIALS_DIR, "main")
MS_SUPP_DIR <- file.path(MS_MATERIALS_DIR, "supp")
MS_FIG_MAIN_DIR <- MS_MAIN_DIR
MS_FIG_SUPP_DIR <- MS_SUPP_DIR
MS_FIG_MAIN_PDF_DIR <- file.path(MS_MAIN_DIR, "pdf")
MS_FIG_SUPP_PDF_DIR <- file.path(MS_SUPP_DIR, "pdf")
MS_TABLES_MAIN_DIR <- MS_MAIN_DIR
MS_TABLES_SUPP_DIR <- MS_SUPP_DIR

# input folders
# small input groups live at the top of inputs/ by default
DISCHARGE_DIR <- BASE_DATA_DIR
EC_DIR <- BASE_DATA_DIR
ISOTOPE_DIR <- BASE_DATA_DIR
STREAM_TEMP_DIR <- BASE_DATA_DIR
MET_DIR <- file.path(BASE_DATA_DIR, "all_hydromet")
CATCHMENT_CHARACTERISTICS_DIR <- BASE_DATA_DIR

EXPLORATORY_ET_METHODS_DIR <- file.path(EXPLORATORY_PLOTS_DIR, "et_methods")

# output options
# leave extra helper exports off unless you need them
WRITE_TABLE_OUTPUTS <- FALSE
WRITE_AUX_OUTPUTS <- FALSE

# create the main output folder if needed
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# output folders
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
# keep older name used in some scripts
OUT_MODELS_WATERSHED_CHAR_STORAGE_MLR_DIR <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR <- file.path(
  OUT_STATS_DIR,
  "storage_eco_response_mlr"
)
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

# site lists

# sites with continuous streamflow used for dynamic storage
# ws06 stays grouped with ws07 and ws08 for plotting
SITE_ORDER_HYDROMETRIC <- c(
  "WS09", # WS09
  "WS10", # WS10
  "WS01", # WS01
  "Look", # Look
  "WS02", # WS02
  "WS03", # WS03
  "Mack", # Mack
  "WS06", # WS06 (no chemistry/isotope data)
  "WS07", # WS07
  "WS08" # WS08
)

# sites with chemistry data
SITE_ORDER_CHEMISTRY <- c(
  "WS09", # WS09
  "WS10", # WS10
  "WS01", # WS01
  "Look", # Look
  "WS02", # WS02
  "WS03", # WS03
  "WS06", # WS06
  "WS07", # WS07
  "WS08", # WS08
  "Mack" # Mack
)

# site order used in most outputs
SITE_ORDER_ALL <- SITE_ORDER_HYDROMETRIC

# site name lookup
SITE_NAMES <- c(
  "WS09" = "WS09",
  "WS10" = "WS10",
  "WS01" = "WS01",
  "Look" = "Look",
  "WS02" = "WS02",
  "WS03" = "WS03",
  "WS06" = "WS06",
  "WS07" = "WS07",
  "WS08" = "WS08",
  "Mack" = "Mack",
  "MR" = "MR",
  "NC" = "NC",
  "LC" = "LC",
  "LO2" = "LO2",
  "CC" = "CC",
  "LO1" = "LO1"
)

# water year range

WY_START <- 1997
WY_END <- 2020

# minimum number of daily EC + Q values needed to keep a BF water year
BF_MIN_DAYS_PER_WY <- 300
CHS_MIN_DAYS_PER_WY <- BF_MIN_DAYS_PER_WY

# minimum number of chemistry samples per water year for sample-based BF
# based on the CF002 chemistry records such as cond and Ca
BF_MIN_OBS_PER_WY_CHEM <- 10
CHS_MIN_OBS_PER_WY_CHEM <- BF_MIN_OBS_PER_WY_CHEM

# isotope metrics are used as site averages in this run
# there is no annual isotope branch in the main code path

# sites intentionally excluded from BF-based summaries
BF_EXCLUDE_SITES <- character(0)

# storage metric groups

# default storage metric order across figures and tables
# rbi, rcs, fdc, sd, wb, bf, dr, fyw, mtt
STORAGE_METRIC_ORDER <- c(
  "RBI",
  "RCS",
  "FDC",
  "SD",
  "WB",
  "BF",
  "DR",
  "Fyw",
  "MTT"
)

# dynamic storage metrics
# rbi = Richards-Baker index, rcs = recession-curve slope
# fdc = full-period site-level flow-duration-curve slope
# annual site-year fdc is used only for Figure 2 and ANOVA/Tukey output
# sd = storage-discharge
DYNAMIC_METRICS <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD")
]

# mobile storage metrics
# bf = baseflow fraction from chemical hydrograph separation
# mtt = mean transit time, fyw = young water fraction, dr = damping ratio
MOBILE_METRICS_ANNUAL <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("BF")
]
MOBILE_METRICS_SITE <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("DR", "Fyw", "MTT")
] # site-level isotope values

# extended-dynamic storage metric from the annual water balance
# wb = water-balance deficit
EXTENDED_DYNAMIC_METRICS <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("WB")
]

# plotting order
PLOT_ORDER_DYNAMIC_STORAGE <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD", "WB")
]
PLOT_ORDER_MOBILE_STORAGE <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("BF", "DR", "Fyw", "MTT")
]
PLOT_MOBILE_STORAGE_SITE_COLS <- c(
  "BF" = "BF_mean",
  "DR" = "DR",
  "Fyw" = "Fyw",
  "MTT" = "MTT"
)

# all storage metrics
ALL_STORAGE_METRICS <- STORAGE_METRIC_ORDER

# shared file names and code mappings

MASTER_ANNUAL_FILE <- "master_annual.csv"
MASTER_SITE_FILE <- "master_site.csv"

# raw site codes to leave out of analysis tables
SITE_EXCLUDE_RAW <- c("GSWSMA", "GSWSMF")

# site recode map used when bringing met and discharge records into site codes
SITECODE_RECODE_TO_GSMACK <- c("GSWSMC" = "GSMACK")

# component sites used to make GSLOOK met composites
GSLOOK_COMPOSITE_COMPONENT_SITES <- c("GSWS01", "GSWS06", "LONGER", "COLD")


# plotting defaults ----

# shared plot text size
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

# shared label settings
# use these so labels stay readable across figures
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

# 10-color palette for streamflow sites
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

# plot theme used across the paper
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
    # numbered sites: gsws## -> ws##
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

STORAGE_METRIC_FULL_NAMES <- c(
  "RBI" = "Richards-Baker Index",
  "RCS" = "Recession Curve Slope",
  "FDC" = "Flow Duration Curve Slope",
  "SD" = "Storage-Discharge",
  "WB" = "Water-balance deficit",
  "CHS" = "Baseflow Fraction",
  "DR" = "Damping Ratio",
  "Fyw" = "Young Water Fraction",
  "MTT" = "Mean Transit Time"
)

METRIC_ABBREV_DISPLAY <- c(
  "RBI" = "RBI",
  "RCS" = "RCS",
  "FDC" = "FDC",
  "SD" = "SD",
  "WB" = "WB",
  "CHS" = "BF",
  "DR" = "DR",
  "Fyw" = "Fyw",
  "MTT" = "MTT"
)

ECO_RESPONSE_FULL_NAMES <- c(
  "Q7Q5" = "Seven-Day Low-Flow Discharge (Q7Q5)",
  "T7DMax" = "Seven-Day Maximum Stream Temperature (T7DMax)"
)

ECO_PREDICTOR_FULL_NAMES <- c(
  "Pws" = "Wet-Season Precipitation (Pws)",
  "RBI" = "Richards-Baker Index (RBI)",
  "RCS" = "Recession Curve Slope (RCS)",
  "FDC" = "Flow Duration Curve Slope (FDC)",
  "SD" = "Storage-Discharge (SD)",
  "WB" = "Water-balance deficit (WB)",
  "CHS" = "Baseflow Fraction (BF)",
  "DR" = "Damping Ratio (DR)",
  "Fyw" = "Young Water Fraction (Fyw)",
  "MTT" = "Mean Transit Time (MTT)"
)

label_metric_abbrev <- function(x) {
  x_chr <- as.character(x)
  out <- unname(METRIC_ABBREV_DISPLAY[x_chr])
  miss <- is.na(out)
  out[miss] <- x_chr[miss]
  out
}

label_storage_metric <- function(x, include_abbrev = TRUE) {
  x_chr <- as.character(x)
  long <- unname(STORAGE_METRIC_FULL_NAMES[x_chr])
  abbr <- label_metric_abbrev(x_chr)
  miss <- is.na(long)
  long[miss] <- x_chr[miss]
  if (!isTRUE(include_abbrev)) {
    return(long)
  }
  out <- paste0(long, " (", abbr, ")")
  out[miss] <- x_chr[miss]
  out
}

label_eco_response <- function(x) {
  x_chr <- as.character(x)
  out <- unname(ECO_RESPONSE_FULL_NAMES[x_chr])
  miss <- is.na(out)
  out[miss] <- x_chr[miss]
  out
}

label_eco_predictor <- function(x) {
  x_chr <- as.character(x)
  out <- unname(ECO_PREDICTOR_FULL_NAMES[x_chr])
  miss <- is.na(out)
  out[miss] <- x_chr[miss]
  out
}

wrap_plot_label <- function(x, width = 22) {
  x_chr <- as.character(x)
  vapply(
    x_chr,
    function(lbl) {
      parts <- strwrap(lbl, width = width)
      if (length(parts) == 0) {
        return("")
      }
      paste(parts, collapse = "\n")
    },
    FUN.VALUE = character(1)
  )
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
      if (is.na(val) || !nzchar(val)) {
        return(val)
      }
      parts <- trimws(strsplit(val, sep, fixed = TRUE)[[1]])
      parts <- parts[nzchar(parts)]
      if (length(parts) == 0) {
        return(val)
      }
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
