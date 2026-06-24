# shared settings used by the analysis scripts

# inputs:
# HJA_REPO_DIR and HJA_FINAL_WORKFLOW_ROOT environment variables, if set

# outputs:
# shared workflow paths, constants, and shared functions

# author: Sidney Bush
# date: 2026-02-13

# paths
# inputs and workflow outputs live outside the Git repository
REPO_DIR <- normalizePath(
  Sys.getenv("HJA_REPO_DIR", unset = getwd()),
  winslash = "/",
  mustWork = FALSE
)
DEFAULT_WORKFLOW_ROOT <- file.path(dirname(REPO_DIR), "hja-dynamic-storage-workflow")
FINAL_WORKFLOW_ROOT <- normalizePath(
  Sys.getenv("HJA_FINAL_WORKFLOW_ROOT", unset = DEFAULT_WORKFLOW_ROOT),
  winslash = "/",
  mustWork = FALSE
)
BASE_DATA_DIR <- file.path(FINAL_WORKFLOW_ROOT, "inputs")
OUTPUT_DIR <- file.path(FINAL_WORKFLOW_ROOT, "outputs")
FIGS_TABLES_PUB_DIR <- file.path(FINAL_WORKFLOW_ROOT, "figs_tables_pub")

# published figure and table files
# run_all.R creates the main and supporting information folders before scripts write files
MS_MAIN_DIR <- file.path(FIGS_TABLES_PUB_DIR, "main")
MS_SUPP_DIR <- file.path(FIGS_TABLES_PUB_DIR, "supp")
MS_FIG_MAIN_DIR <- MS_MAIN_DIR
MS_FIG_SUPP_DIR <- MS_SUPP_DIR
MS_FIG_MAIN_PDF_DIR <- file.path(MS_MAIN_DIR, "pdf")
MS_FIG_SUPP_PDF_DIR <- file.path(MS_SUPP_DIR, "pdf")
MS_FIG_MAIN_TIFF_DIR <- file.path(MS_MAIN_DIR, "tiff")
MS_FIG_SUPP_TIFF_DIR <- file.path(MS_SUPP_DIR, "tiff")
MS_TABLES_MAIN_DIR <- MS_MAIN_DIR
MS_TABLES_SUPP_DIR <- MS_SUPP_DIR

# input folders
# most input files are in inputs, hydromet files are in inputs/all_hydromet
DISCHARGE_DIR <- BASE_DATA_DIR
EC_DIR <- BASE_DATA_DIR
ISOTOPE_DIR <- BASE_DATA_DIR
STREAM_TEMP_DIR <- BASE_DATA_DIR
MET_DIR <- file.path(BASE_DATA_DIR, "all_hydromet")
CATCHMENT_CHARACTERISTICS_DIR <- BASE_DATA_DIR

# output folders
# run_all.R creates the metrics and model folders before scripts write files
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
OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR <- file.path(
  OUT_STATS_DIR,
  "storage_eco_response_mlr"
)

# sites
# order used for streamflow metrics, figures, and most tables
SITE_ORDER_HYDROMETRIC <- c(
  "WS09",
  "WS10",
  "WS01",
  "Look",
  "WS02",
  "WS03",
  "Mack",
  "WS06",
  "WS07",
  "WS08"
)

# order used by chemistry and baseflow scripts
SITE_ORDER_CHEMISTRY <- c(
  "WS09",
  "WS10",
  "WS01",
  "Look",
  "WS02",
  "WS03",
  "WS06",
  "WS07",
  "WS08",
  "Mack"
)
# water years
WY_START <- 1997
WY_END <- 2020

# baseflow fraction filters
BF_MIN_OBS_PER_WY_CHEM <- 10

# isotope metrics are site averages for the manuscript models
# annual isotope files are not written

BF_EXCLUDE_SITES <- character(0)

# storage metrics
# metric order for plots and tables
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

# dynamic metrics are based on discharge
DYNAMIC_METRICS <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD")
]

# BF is annual, isotope metrics are site averages
MOBILE_METRICS_ANNUAL <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("BF")
]
MOBILE_METRICS_SITE <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("DR", "Fyw", "MTT")
]

# WB is kept separate because it is based on the annual water balance
EXTENDED_DYNAMIC_METRICS <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("WB")
]

# shared file names and code mappings
MASTER_ANNUAL_FILE <- "master_annual.csv"
MASTER_SITE_FILE <- "master_site.csv"

# source site codes to leave out of analysis tables
SITE_EXCLUDE_RAW <- c("GSWSMA", "GSWSMF")

# recode source site names before standardizing names
SITECODE_RECODE_TO_GSMACK <- c("GSWSMC" = "GSMACK")

# met sites used to build the Lookout Creek composite record
GSLOOK_COMPOSITE_COMPONENT_SITES <- c("GSWS01", "GSWS06", "LONGER", "COLD")


# shared functions
# config stays limited to paths and constants
source(file.path(REPO_DIR, "figure_functions.R"))
source(file.path(REPO_DIR, "workflow_functions.R"))
source(file.path(REPO_DIR, "hydromet_functions.R"))
source(file.path(REPO_DIR, "model_functions.R"))
