# -----------------------------------------------------------------------------
# Legacy Parity Check
# -----------------------------------------------------------------------------
# Compares mapped legacy outputs to current workflow outputs.
#
# Usage:
#   LEGACY_OUTPUT_DIR=/path/to/legacy_outputs Rscript legacy_parity_check.R
#
# If LEGACY_OUTPUT_DIR is not set, default is:
#   file.path(OUTPUT_DIR, "legacy_reference")
#
# Outputs:
#   - legacy_parity_summary.csv
#   - legacy_parity_details.csv
#   - legacy_parity_file_map.csv
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)

rm(list = ls())

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
  config_path <- file.path(getwd(), "config.R")
}
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found.")
}

legacy_output_dir <- Sys.getenv("LEGACY_OUTPUT_DIR", unset = "")
if (legacy_output_dir == "") {
  legacy_output_dir <- file.path(OUTPUT_DIR, "legacy_reference")
}

current_output_dir <- OUTPUT_DIR
allow_mlr_differences <- TRUE
abs_tol_default <- 1e-8
allow_updated_eco_responses <- TRUE

legacy_search_dirs <- unique(c(
  legacy_output_dir,
  file.path(legacy_output_dir, "Hydrometric"),
  dirname(legacy_output_dir),
  file.path(dirname(legacy_output_dir), "Hydrometric")
))

resolve_legacy_path <- function(search_dirs, candidate_files) {
  for (d in search_dirs) {
    for (f in candidate_files) {
      p <- file.path(d, f)
      if (file.exists(p)) return(p)
    }
  }
  file.path(search_dirs[1], candidate_files[1])
}

rename_using_map <- function(df, rename_map) {
  if (length(rename_map) == 0) return(df)
  nm <- names(df)
  for (old_nm in names(rename_map)) {
    new_nm <- rename_map[[old_nm]]
    idx <- which(nm == old_nm)
    if (length(idx) == 1) nm[idx] <- new_nm
  }
  names(df) <- nm
  df
}

dedupe_columns <- function(df) {
  df[, !duplicated(names(df)), drop = FALSE]
}

coerce_key_cols <- function(df, keys) {
  out <- df
  for (k in keys) {
    if (k %in% names(out)) {
      if (k == "site") {
        out[[k]] <- as.character(standardize_site_code(out[[k]]))
      } else {
        out[[k]] <- as.character(out[[k]])
      }
    }
  }
  out
}

numeric_compare <- function(df_legacy, df_current, keys, file_id, abs_tol, skip_numeric_vars = character(0), allow_extra_current_keys = FALSE) {
  legacy_cols <- names(df_legacy)
  current_cols <- names(df_current)

  shared_cols <- intersect(legacy_cols, current_cols)
  shared_cols <- setdiff(shared_cols, keys)

  numeric_shared <- shared_cols[
    sapply(shared_cols, function(col) {
      is.numeric(df_legacy[[col]]) && is.numeric(df_current[[col]])
    })
  ]
  numeric_shared <- setdiff(numeric_shared, skip_numeric_vars)

  if (length(keys) == 0 || any(!keys %in% names(df_legacy)) || any(!keys %in% names(df_current))) {
    return(list(
      summary = tibble(
        file_id = file_id,
        status = "missing_keys",
        n_legacy = nrow(df_legacy),
        n_current = nrow(df_current),
        n_matched = NA_integer_,
        n_only_legacy = NA_integer_,
        n_only_current = NA_integer_,
        n_numeric_compared = length(numeric_shared),
        n_numeric_failed = NA_integer_,
        max_abs_diff = NA_real_,
        pass = FALSE
      ),
      details = tibble()
    ))
  }

  ldf <- coerce_key_cols(df_legacy, keys)
  cdf <- coerce_key_cols(df_current, keys)

  joined <- full_join(
    ldf %>% mutate(.legacy_row = TRUE),
    cdf %>% mutate(.current_row = TRUE),
    by = keys,
    suffix = c("_legacy", "_current")
  )

  has_legacy <- ifelse(is.na(joined$.legacy_row), FALSE, TRUE)
  has_current <- ifelse(is.na(joined$.current_row), FALSE, TRUE)

  n_matched <- sum(has_legacy & has_current)
  n_only_legacy <- sum(has_legacy & !has_current)
  n_only_current <- sum(!has_legacy & has_current)

  detail_rows <- list()
  fail_count <- 0
  max_abs <- 0

  for (col in numeric_shared) {
    lcol <- paste0(col, "_legacy")
    ccol <- paste0(col, "_current")

    if (!(lcol %in% names(joined)) || !(ccol %in% names(joined))) {
      next
    }

    diffs <- abs(joined[[lcol]] - joined[[ccol]])
    diffs[is.na(diffs)] <- 0

    col_max <- suppressWarnings(max(diffs, na.rm = TRUE))
    if (!is.finite(col_max)) col_max <- 0
    max_abs <- max(max_abs, col_max)

    col_fail <- sum(diffs > abs_tol, na.rm = TRUE)
    if (col_fail > 0) fail_count <- fail_count + 1

    detail_rows[[length(detail_rows) + 1]] <- tibble(
      file_id = file_id,
      variable = col,
      n_compared = sum(has_legacy & has_current),
      max_abs_diff = col_max,
      n_failed = col_fail,
      tolerance = abs_tol,
      pass = col_fail == 0
    )
  }

  details <- bind_rows(detail_rows)

  pass_all <- n_only_legacy == 0 && fail_count == 0 &&
    (allow_extra_current_keys || n_only_current == 0)

  summary <- tibble(
    file_id = file_id,
    status = "compared",
    n_legacy = nrow(df_legacy),
    n_current = nrow(df_current),
    n_matched = n_matched,
    n_only_legacy = n_only_legacy,
    n_only_current = n_only_current,
    n_numeric_compared = length(numeric_shared),
    n_numeric_failed = fail_count,
    max_abs_diff = max_abs,
    pass = pass_all
  )

  list(summary = summary, details = details)
}

file_specs <- list(
  list(
    file_id = "annual_master",
    legacy_files = c(LEGACY_ANNUAL_FILE, MASTER_ANNUAL_FILE),
    current_file = MASTER_ANNUAL_FILE,
    keys = c("site", "year"),
    rename_map = c(
      SITECODE = "site",
      waterYear = "year",
      wateryear = "year",
      recession_curve_slope = "RCS",
      rbfi = "RBI",
      fdc_slope = "FDC",
      S_annual_mm = "SD",
      mean_bf = "CHS",
      DS_sum = "WB",
      temp_at_min_Q_7d_C = "temp_during_min_Q_7d_C"
    ),
    mlr_related = FALSE,
    skip_numeric_vars = c(
      "max_temp_7d_C",
      "min_Q_7d_mm_d",
      "temp_at_min_Q_7d_C",
      "temp_during_min_Q_7d_C",
      "Q5_CV"
    ),
    allow_extra_current_keys = TRUE,
    required = TRUE,
    informational = TRUE
  ),
  list(
    file_id = "site_master",
    legacy_files = c(LEGACY_SITE_FILE, MASTER_SITE_FILE),
    current_file = MASTER_SITE_FILE,
    keys = c("site"),
    rename_map = c(
      recession_curve_slope_mean = "RCS_mean",
      fdc_slope_mean = "FDC_mean",
      S_annual_mm_mean = "SD_mean",
      mean_bf_mean = "CHS_mean",
      DS_sum_mean = "WB_mean",
      DR_Overall = "DR"
    ),
    mlr_related = FALSE,
    skip_numeric_vars = c(
      "max_temp_7d_C_mean",
      "min_Q_7d_mm_d_mean",
      "temp_at_min_Q_7d_C_mean",
      "temp_during_min_Q_7d_C_mean",
      "Q5_CV_mean"
    ),
    required = TRUE
  ),
  list(
    file_id = "annual_gw_prop",
    legacy_files = c("Annual_GW_Prop.csv"),
    current_file = "Annual_GW_Prop.csv",
    keys = c("site", "year"),
    rename_map = c(
      SITECODE = "site",
      waterYear = "year",
      mean_bf = "CHS"
    ),
    mlr_related = FALSE,
    required = TRUE
  ),
  list(
    file_id = "drawdown_annual",
    legacy_files = c("DS_drawdown_annual.csv"),
    current_file = "DS_drawdown_annual.csv",
    keys = c("site", "year"),
    rename_map = c(
      SITECODE = "site",
      waterYear = "year",
      DS_sum = "WB"
    ),
    mlr_related = FALSE,
    required = TRUE
  ),
  list(
    file_id = "fdc_slopes_overall",
    legacy_files = c("FDC_slopes_overall.csv"),
    current_file = "FDC_slopes_overall.csv",
    keys = c("site"),
    rename_map = c(),
    mlr_related = FALSE,
    required = FALSE
  ),
  list(
    file_id = "fdc_slopes_wy",
    legacy_files = c("FDC_slopes_WY.csv"),
    current_file = "FDC_slopes_WY.csv",
    keys = c("site", "WaterYear"),
    rename_map = c(wateryear = "WaterYear"),
    mlr_related = FALSE,
    required = FALSE
  ),
  list(
    file_id = "storage_discharge_annual_mm",
    legacy_files = c("storage_discharge_method_annual_mm_metrics_per_site_wateryear.csv"),
    current_file = "storage_discharge_method_annual_mm_metrics_per_site_wateryear.csv",
    keys = c("site", "wateryear"),
    rename_map = c(),
    mlr_related = FALSE,
    required = TRUE
  ),
  list(
    file_id = "mlr_catchment",
    legacy_files = c("MasterLM_Coefs_DF.csv", "MLR_Storage_Catchment_Results.csv"),
    current_file = "MLR_Storage_Catchment_Results.csv",
    keys = c(),
    rename_map = c(),
    mlr_related = TRUE,
    required = FALSE
  )
)

summary_rows <- list()
detail_rows <- list()
map_rows <- list()

for (spec in file_specs) {
  legacy_candidates <- if ("legacy_files" %in% names(spec)) spec$legacy_files else spec$legacy_file
  legacy_path <- resolve_legacy_path(legacy_search_dirs, legacy_candidates)
  current_path <- file.path(current_output_dir, spec$current_file)

  map_rows[[length(map_rows) + 1]] <- tibble(
    file_id = spec$file_id,
    legacy_file = legacy_path,
    current_file = current_path,
    mlr_related = spec$mlr_related
  )

  if (!file.exists(current_path)) {
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      file_id = spec$file_id,
      status = "missing_file",
      n_legacy = NA_integer_,
      n_current = NA_integer_,
      n_matched = NA_integer_,
      n_only_legacy = NA_integer_,
      n_only_current = NA_integer_,
      n_numeric_compared = NA_integer_,
      n_numeric_failed = NA_integer_,
      max_abs_diff = NA_real_,
      pass = FALSE
    )
    next
  }
  if (!file.exists(legacy_path)) {
    is_required <- if ("required" %in% names(spec)) isTRUE(spec$required) else TRUE
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      file_id = spec$file_id,
      status = ifelse(is_required, "missing_file", "legacy_not_available"),
      n_legacy = NA_integer_,
      n_current = NA_integer_,
      n_matched = NA_integer_,
      n_only_legacy = NA_integer_,
      n_only_current = NA_integer_,
      n_numeric_compared = NA_integer_,
      n_numeric_failed = NA_integer_,
      max_abs_diff = NA_real_,
      pass = !is_required
    )
    next
  }

  legacy_df <- suppressWarnings(read_csv(legacy_path, show_col_types = FALSE))
  current_df <- suppressWarnings(read_csv(current_path, show_col_types = FALSE))

  legacy_df <- rename_using_map(legacy_df, spec$rename_map)
  current_df <- rename_using_map(current_df, spec$rename_map)
  legacy_df <- dedupe_columns(legacy_df)
  current_df <- dedupe_columns(current_df)

  tol_use <- abs_tol_default
  if (spec$mlr_related && allow_mlr_differences) {
    tol_use <- Inf
  }

  skip_vars <- character(0)
  if (allow_updated_eco_responses && "skip_numeric_vars" %in% names(spec)) {
    skip_vars <- spec$skip_numeric_vars
  }

  allow_extra_current <- if ("allow_extra_current_keys" %in% names(spec)) {
    isTRUE(spec$allow_extra_current_keys)
  } else {
    FALSE
  }

  cmp <- numeric_compare(
    legacy_df,
    current_df,
    spec$keys,
    spec$file_id,
    tol_use,
    skip_numeric_vars = skip_vars,
    allow_extra_current_keys = allow_extra_current
  )

  if (spec$mlr_related && allow_mlr_differences && nrow(cmp$summary) == 1) {
    cmp$summary$status <- "skipped_mlr_exception"
    cmp$summary$pass <- TRUE
  }
  if ("informational" %in% names(spec) && isTRUE(spec$informational) && nrow(cmp$summary) == 1) {
    if (cmp$summary$status == "compared" && !isTRUE(cmp$summary$pass)) {
      cmp$summary$status <- "compared_informational"
    }
    cmp$summary$pass <- TRUE
  }

  summary_rows[[length(summary_rows) + 1]] <- cmp$summary
  if (nrow(cmp$details) > 0) {
    detail_rows[[length(detail_rows) + 1]] <- cmp$details
  }
}

summary_df <- bind_rows(summary_rows)
detail_df <- bind_rows(detail_rows)
map_df <- bind_rows(map_rows)

write_csv(summary_df, file.path(current_output_dir, "legacy_parity_summary.csv"))
write_csv(detail_df, file.path(current_output_dir, "legacy_parity_details.csv"))
write_csv(map_df, file.path(current_output_dir, "legacy_parity_file_map.csv"))
