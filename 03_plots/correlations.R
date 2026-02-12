# -----------------------------------------------------------------------------
# Correlation and Scatter-Matrix Plots
# -----------------------------------------------------------------------------
# This script makes the core correlation heatmaps used in the workflow.
#
# Inputs:
#   - master_site.csv
#   - master_annual.csv
#
# Outputs:
#   - catch_chars_storage_mlr_corr.png
#   - storage_eco_model_corrplot.png
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(ggcorrplot)

theme_set(theme_classic(base_size = 12))

rm(list = ls())

# get script directory (works with source() and Rscript)
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    getwd()
  }
})
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

output_dir <- OUT_MASTER_DIR
plot_dir <- file.path(FIGURES_DIR, "supp", "stats", "correlations")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Clear prior correlation files so the folder reflects current workflow outputs.
old_corr_files <- list.files(plot_dir, pattern = "corrplot\\.png$", full.names = TRUE)
if (length(old_corr_files) > 0) {
  unlink(old_corr_files)
}

# -----------------------------------------------------------------------------
# 1. load data
# -----------------------------------------------------------------------------

site_file <- file.path(output_dir, MASTER_SITE_FILE)
if (!file.exists(site_file)) {
  site_file <- file.path(OUTPUT_DIR, MASTER_SITE_FILE)
}

HJA_Ave <- read_csv(
  site_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

# -----------------------------------------------------------------------------
# 2. catch_chars MLR family correlation matrix
# -----------------------------------------------------------------------------

watershed_predictors <- c(
  "Slope_mean", "Harvest", "Landslide_Total", "Landslide_Young",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
)
watershed_predictors <- watershed_predictors[watershed_predictors %in% names(HJA_Ave)]

if (length(watershed_predictors) >= 2) {
  cor_catchment <- cor(HJA_Ave[watershed_predictors], use = "pairwise.complete.obs")

  p_catchment <- ggcorrplot(
    cor_catchment,
    hc.order = FALSE,
    type = "lower",
    outline.col = "white",
    lab = TRUE
  ) + labs(title = "Watershed Characteristics MLR Predictor Correlation Matrix")

  ggsave(
    file.path(plot_dir, "catch_chars_storage_mlr_corr.png"),
    p_catchment,
    width = 9,
    height = 9,
    dpi = 300
  )
}

# -----------------------------------------------------------------------------
# 3. storage->eco MLR family correlation matrix
# -----------------------------------------------------------------------------

annual_file <- file.path(output_dir, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUTPUT_DIR, MASTER_ANNUAL_FILE)
}
if (!file.exists(annual_file)) {
  annual_file <- file.path(BASE_DATA_DIR, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
}

HJA_Yr <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

if (!("q5_7d_mm_d" %in% names(HJA_Yr)) && ("min_Q_7d_mm_d" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(q5_7d_mm_d = min_Q_7d_mm_d)
}
if (!("T_7DMax" %in% names(HJA_Yr)) && ("max_temp_7d_C" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(T_7DMax = max_temp_7d_C)
}
if (!("Q_7Q5" %in% names(HJA_Yr)) && ("q5_7d_mm_d" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(Q_7Q5 = q5_7d_mm_d)
}
if (!("temp_during_q5_7d_C" %in% names(HJA_Yr)) && ("temp_during_min_Q_7d_C" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(temp_during_q5_7d_C = temp_during_min_Q_7d_C)
}
if (!("temp_during_q5_7d_C" %in% names(HJA_Yr)) && ("temp_at_min_Q_7d_C" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(temp_during_q5_7d_C = temp_at_min_Q_7d_C)
}
if (!("T_Q7Q5" %in% names(HJA_Yr)) && ("temp_during_q5_7d_C" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(T_Q7Q5 = temp_during_q5_7d_C)
}

eco_response_vars <- c("T_7DMax", "Q_7Q5", "T_Q7Q5")
eco_response_vars <- eco_response_vars[eco_response_vars %in% names(HJA_Yr)]

storage_predictor_vars <- c("RCS", "RBI", "FDC", "SD", "WB", "CHS", "MTT", "Fyw", "DR")
storage_predictor_vars <- storage_predictor_vars[storage_predictor_vars %in% names(HJA_Yr)]

eco_corr_vars <- unique(c(eco_response_vars, storage_predictor_vars))
eco_corr_vars <- eco_corr_vars[eco_corr_vars %in% names(HJA_Yr)]

if (length(eco_corr_vars) >= 2) {
  cor_eco <- cor(HJA_Yr[eco_corr_vars], use = "pairwise.complete.obs")

  p_eco <- ggcorrplot(
    cor_eco,
    hc.order = FALSE,
    type = "lower",
    outline.col = "white",
    lab = FALSE
  ) + labs(title = "Storage->Eco MLR Predictor/Response Correlation Matrix")

  ggsave(
    file.path(plot_dir, "storage_eco_mlr_corr.png"),
    p_eco,
    width = 11,
    height = 11,
    dpi = 300
  )
}
