#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
})

find_repo_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(10)) {
    has_config <- file.exists(file.path(cur, "config.R"))
    has_metrics <- dir.exists(file.path(cur, "01_metrics"))
    if (has_config && has_metrics) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(start_dir, winslash = "/", mustWork = FALSE)
}

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
if (is.na(script_path) || script_path == "") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
}
start_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_root <- find_repo_root(start_dir)

source(file.path(repo_root, "config.R"))

required_files <- c(
  file.path(BASE_DATA_DIR, "all_hydromet", "Temperature_original_&_filled_1979_2023_v2.csv"),
  file.path(BASE_DATA_DIR, "all_hydromet", "Precipitation_original_&_filled_1979_2023.csv"),
  file.path(BASE_DATA_DIR, "all_hydromet", "MS00102_v9.csv"),
  file.path(BASE_DATA_DIR, "all_hydromet", "MS05025_v3.csv"),
  file.path(BASE_DATA_DIR, "all_hydromet", "MS00403_v2.csv"),
  file.path(BASE_DATA_DIR, "all_hydromet", "HF00402_v14.csv"),
  file.path(BASE_DATA_DIR, "all_hydromet", "drainage_area.csv"),
  file.path(BASE_DATA_DIR, "Q", "HF00402_v14.csv"),
  file.path(BASE_DATA_DIR, "Q", "drainage_area.csv"),
  file.path(BASE_DATA_DIR, "DynamicStorage", "daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv"),
  file.path(BASE_DATA_DIR, "DynamicStorage", "Catchment_Charc.csv"),
  file.path(BASE_DATA_DIR, "Isotopes", "MTT_FYW.csv"),
  file.path(BASE_DATA_DIR, "Isotopes", "DampingRatios_2025-07-07.csv"),
  file.path(BASE_DATA_DIR, "Stream_T", "HT00401_v8.csv")
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required input file(s):\n- ",
      paste(missing_files, collapse = "\n- ")
    )
  )
}

# Accept any known EC filename used across legacy/current variants.
ec_candidates <- c(
  file.path(BASE_DATA_DIR, "EC", "CF01201_v3.txt"),
  file.path(BASE_DATA_DIR, "EC", "CF01203_v8.csv"),
  file.path(BASE_DATA_DIR, "EC", "CF01203_v9.csv")
)
if (!any(file.exists(ec_candidates))) {
  stop(
    paste0(
      "Missing EC input file. Expected one of:\\n- ",
      paste(ec_candidates, collapse = "\\n- ")
    )
  )
}

# Minimal schema checks on key files.
check_columns <- function(path, required_cols) {
  cols <- names(suppressMessages(read_csv(path, n_max = 0, show_col_types = FALSE)))
  missing <- setdiff(required_cols, cols)
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing required columns in ", path, ": ",
        paste(missing, collapse = ", ")
      )
    )
  }
}

check_columns(file.path(BASE_DATA_DIR, "Q", "HF00402_v14.csv"), c("DATE", "SITECODE", "MEAN_Q"))
temp_cols <- names(suppressMessages(read_csv(file.path(BASE_DATA_DIR, "Stream_T", "HT00401_v8.csv"), n_max = 0, show_col_types = FALSE)))
if (!all(c("DATE", "SITECODE") %in% temp_cols)) {
  stop("Missing required columns in Stream_T/HT00401_v8.csv: DATE and/or SITECODE")
}
if (!any(c("AIRTEMP_MEAN_DAY", "WATERTEMP_MEAN") %in% temp_cols)) {
  stop("Missing stream temperature column in Stream_T/HT00401_v8.csv (expected AIRTEMP_MEAN_DAY or WATERTEMP_MEAN)")
}
check_columns(file.path(BASE_DATA_DIR, "DynamicStorage", "Catchment_Charc.csv"), c("Site"))

message("Input check passed.")
