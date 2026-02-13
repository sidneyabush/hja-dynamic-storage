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
  file.path(MET_DIR, "Temperature_original_&_filled_1979_2023_v2.csv"),
  file.path(MET_DIR, "Precipitation_original_&_filled_1979_2023.csv"),
  file.path(MET_DIR, "MS00102_v9.csv"),
  file.path(MET_DIR, "MS05025_v3.csv"),
  file.path(MET_DIR, "MS00403_v2.csv"),
  file.path(DISCHARGE_DIR, "HF00402_v14.csv"),
  resolve_drainage_area_file(),
  resolve_catchment_characteristics_file(),
  file.path(ISOTOPE_DIR, "MTT_FYW.csv"),
  file.path(ISOTOPE_DIR, "DampingRatios_2025-07-07.csv"),
  file.path(STREAM_TEMP_DIR, "HT00451_v10.txt")
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
  file.path(EC_DIR, "CF01201_v3.txt"),
  file.path(EC_DIR, "CF01203_v8.csv"),
  file.path(EC_DIR, "CF01203_v9.csv")
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

check_columns(file.path(DISCHARGE_DIR, "HF00402_v14.csv"), c("DATE", "SITECODE", "MEAN_Q"))
temp_cols <- names(suppressMessages(read_csv(file.path(STREAM_TEMP_DIR, "HT00451_v10.txt"), n_max = 0, show_col_types = FALSE)))
if (!all(c("DATE_TIME", "SITECODE") %in% temp_cols)) {
  stop("Missing required columns in Stream_T/HT00451_v10.txt: DATE_TIME and/or SITECODE")
}
if (!("WATERTEMP_MEAN" %in% temp_cols)) {
  stop("Missing stream temperature column in Stream_T/HT00451_v10.txt (expected WATERTEMP_MEAN)")
}
check_columns(resolve_catchment_characteristics_file(), c("Site"))

# Derived file check (must come from workflow outputs)
invisible(resolve_water_balance_daily_file())

