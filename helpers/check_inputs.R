# inputs: no direct csv file reads in this script.
# author: sidney bush
# date: 2026-02-13

suppressPackageStartupMessages({
  library(readr)
})

find_repo_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(10)) {
    has_config <- file.exists(file.path(cur, "config.R"))
    has_metrics <- dir.exists(file.path(cur, "01_storage_calcs"))
    if (has_config && has_metrics) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      break
    }
    cur <- parent
  }
  normalizePath(start_dir, winslash = "/", mustWork = FALSE)
}

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
  error = function(e) NA_character_
)
if (is.na(script_path) || script_path == "") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(
      sub("^--file=", "", file_arg[1]),
      winslash = "/",
      mustWork = FALSE
    )
  }
}
start_dir <- if (!is.na(script_path) && nzchar(script_path)) {
  dirname(script_path)
} else {
  getwd()
}
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

# ec file used by the current workflow.
ec_file <- file.path(EC_DIR, "CF01201_v4.txt")
if (!file.exists(ec_file)) {
  stop("Missing required EC input file: ", ec_file)
}

# minimal schema checks on key files.
check_columns <- function(path, required_cols) {
  cols <- names(suppressMessages(read_csv(
    path,
    n_max = 0,
    show_col_types = FALSE
  )))
  missing <- setdiff(required_cols, cols)
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing required columns in ",
        path,
        ": ",
        paste(missing, collapse = ", ")
      )
    )
  }
}

check_columns(
  file.path(DISCHARGE_DIR, "HF00402_v14.csv"),
  c("DATE", "SITECODE", "MEAN_Q")
)
temp_file <- file.path(STREAM_TEMP_DIR, "HT00451_v10.txt")
temp_cols <- names(suppressMessages(read_csv(
  temp_file,
  n_max = 0,
  show_col_types = FALSE
)))
if (!all(c("DATE_TIME", "SITECODE") %in% temp_cols)) {
  stop("Missing required columns in ", temp_file, ": DATE_TIME and/or SITECODE")
}
if (!("WATERTEMP_MEAN" %in% temp_cols)) {
  stop(
    "Missing stream temperature column in ",
    temp_file,
    " (expected WATERTEMP_MEAN)"
  )
}
check_columns(resolve_catchment_characteristics_file(), c("Site"))
