# =============================================================================
# Clean Final Workflow Data - Remove Q5 columns and update metric names
# =============================================================================
# Purpose: Remove Q5norm and CV_Q5norm (response variables, not storage metrics)
#          and update column names to new method abbreviations
#
# Column name changes:
#   recession_curve_slope -> RCS
#   fdc_slope -> FDC
#   S_annual_mm -> SD
#   mean_bf -> CHS
#   DS_sum -> WB
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

library(dplyr)
library(readr)

# Clear environment
rm(list = ls())

# =============================================================================
# SOURCE CONFIGURATION
# =============================================================================

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
  script_dir <- file.path(getwd(), "07_Plots")
}

config_path <- file.path(dirname(script_dir), "config.R")
if (!file.exists(config_path)) {
  config_path <- file.path(getwd(), "config.R")
}
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found.")
}

# =============================================================================
# SETUP PATHS
# =============================================================================

final_workflow_dir <- file.path(BOX_BASE_DIR, "Final_Workflow", "02_Processed_Data")

# Columns to remove (response variables, not storage metrics)
q5_cols <- c("Q5norm", "CV_Q5norm", "Q5_CV", "Q5norm_mean", "CV_Q5norm_mean", "Q5_CV_mean")

# Column renaming map
rename_map <- c(
  "recession_curve_slope" = "RCS",
  "fdc_slope" = "FDC",
  "S_annual_mm" = "SD",
  "mean_bf" = "CHS",
  "DS_sum" = "WB",
  "recession_curve_slope_mean" = "RCS_mean",
  "recession_curve_slope_sd" = "RCS_sd",
  "recession_curve_slope_min" = "RCS_min",
  "recession_curve_slope_max" = "RCS_max",
  "recession_curve_slope_n" = "RCS_n",
  "fdc_slope_mean" = "FDC_mean",
  "fdc_slope_sd" = "FDC_sd",
  "fdc_slope_min" = "FDC_min",
  "fdc_slope_max" = "FDC_max",
  "fdc_slope_n" = "FDC_n",
  "S_annual_mm_mean" = "SD_mean",
  "S_annual_mm_sd" = "SD_sd",
  "S_annual_mm_min" = "SD_min",
  "S_annual_mm_max" = "SD_max",
  "S_annual_mm_n" = "SD_n",
  "mean_bf_mean" = "CHS_mean",
  "mean_bf_sd" = "CHS_sd",
  "mean_bf_min" = "CHS_min",
  "mean_bf_max" = "CHS_max",
  "mean_bf_n" = "CHS_n",
  "DS_sum_mean" = "WB_mean",
  "DS_sum_sd" = "WB_sd",
  "DS_sum_min" = "WB_min",
  "DS_sum_max" = "WB_max",
  "DS_sum_n" = "WB_n"
)

# Helper function to clean a CSV file
clean_csv <- function(file_path, remove_cols, rename_cols, filter_metric_col = NULL) {
  if (!file.exists(file_path)) {
    cat("  File not found:", file_path, "\n")
    return(FALSE)
  }

  cat("  Processing:", basename(file_path), "\n")

  df <- read_csv(file_path, show_col_types = FALSE)
  original_cols <- ncol(df)
  original_rows <- nrow(df)

  # Remove Q5 columns
  cols_to_remove <- intersect(names(df), remove_cols)
  if (length(cols_to_remove) > 0) {
    df <- df %>% select(-all_of(cols_to_remove))
    cat("    Removed columns:", paste(cols_to_remove, collapse = ", "), "\n")
  }

  # Rename columns
  cols_to_rename <- intersect(names(df), names(rename_cols))
  if (length(cols_to_rename) > 0) {
    df <- df %>% rename_with(~ rename_cols[.x], all_of(cols_to_rename))
    cat("    Renamed columns:", paste(cols_to_rename, "->", rename_cols[cols_to_rename], collapse = ", "), "\n")
  }

  # Filter rows if metric column specified (for long format files)
  if (!is.null(filter_metric_col) && filter_metric_col %in% names(df)) {
    q5_metrics <- c("Q5norm", "CV_Q5norm", "Q5_CV")
    rows_before <- nrow(df)
    df <- df %>% filter(!(.data[[filter_metric_col]] %in% q5_metrics))
    rows_removed <- rows_before - nrow(df)
    if (rows_removed > 0) {
      cat("    Removed", rows_removed, "rows with Q5 metrics\n")
    }

    # Also update metric names in the metric column
    df <- df %>%
      mutate(!!filter_metric_col := case_when(
        .data[[filter_metric_col]] == "recession_curve_slope" ~ "RCS",
        .data[[filter_metric_col]] == "fdc_slope" ~ "FDC",
        .data[[filter_metric_col]] == "S_annual_mm" ~ "SD",
        .data[[filter_metric_col]] == "mean_bf" ~ "CHS",
        .data[[filter_metric_col]] == "DS_sum" ~ "WB",
        TRUE ~ .data[[filter_metric_col]]
      ))
  }

  # Write back
  write_csv(df, file_path)
  cat("    Saved: ", ncol(df), "cols (was", original_cols, "),", nrow(df), "rows (was", original_rows, ")\n")

  return(TRUE)
}

# =============================================================================
# CLEAN STORAGE_METRICS FILES
# =============================================================================

cat("\n=== Cleaning Storage_Metrics files ===\n")

clean_csv(
  file.path(final_workflow_dir, "Storage_Metrics", "HJA_StorageMetrics_Annual_All.csv"),
  q5_cols, rename_map
)

clean_csv(
  file.path(final_workflow_dir, "Storage_Metrics", "StorageDischarge_FDC_Annual.csv"),
  q5_cols, rename_map
)

clean_csv(
  file.path(final_workflow_dir, "Storage_Metrics", "HJA_Ave_StorageMetrics_CatCharacter.csv"),
  q5_cols, rename_map
)

# Also clean the raw storage_discharge file if it exists
clean_csv(
  file.path(final_workflow_dir, "Storage_Metrics", "storage_discharge_method_annual_mm_metrics_per_site_wateryear.csv"),
  q5_cols, rename_map
)

# =============================================================================
# CLEAN SUMMARY_STATS FILES
# =============================================================================

cat("\n=== Cleaning Summary_Stats files ===\n")

clean_csv(
  file.path(final_workflow_dir, "Summary_Stats", "Site_Summary_Statistics.csv"),
  q5_cols, rename_map
)

clean_csv(
  file.path(final_workflow_dir, "Summary_Stats", "Site_Metric_Summary_Long.csv"),
  q5_cols, rename_map, filter_metric_col = "metric"
)

# =============================================================================
# CLEAN STATISTICAL_RESULTS FILES
# =============================================================================

cat("\n=== Cleaning Statistical_Results files ===\n")

# Correlations file - need special handling for row names
corr_file <- file.path(final_workflow_dir, "Statistical_Results", "Correlations_Storage_Thermal_LowFlow.csv")
if (file.exists(corr_file)) {
  cat("  Processing:", basename(corr_file), "\n")
  df <- read_csv(corr_file, show_col_types = FALSE)

  # Remove Q5 columns
  cols_to_remove <- intersect(names(df), q5_cols)
  if (length(cols_to_remove) > 0) {
    df <- df %>% select(-all_of(cols_to_remove))
    cat("    Removed columns:", paste(cols_to_remove, collapse = ", "), "\n")
  }

  # Rename columns
  cols_to_rename <- intersect(names(df), names(rename_map))
  if (length(cols_to_rename) > 0) {
    df <- df %>% rename_with(~ rename_map[.x], all_of(cols_to_rename))
    cat("    Renamed columns:", paste(cols_to_rename, collapse = ", "), "\n")
  }

  write_csv(df, corr_file)
  cat("    Saved\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== DATA CLEANUP COMPLETE ===\n")
cat("Removed Q5norm and CV_Q5norm data from all files\n")
cat("Updated column names to method abbreviations (RCS, FDC, SD, CHS, WB)\n")
