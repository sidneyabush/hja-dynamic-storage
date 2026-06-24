# compare workflow tables so unit mistakes are caught before modeling

# inputs:
# out_met_support_dir/catchments_met_q.csv
# out_met_support_dir/daily_water_balance_et_hamon_zhang_coeff_interp.csv
# discharge_dir/HF00402_v14.csv
# catchment_characteristics_dir/drainage_area.csv
# out_met_dynamic_dir/storage_discharge_fdc_annual.csv
# output_dir/master/master_annual.csv
# output_dir/master/master_site.csv
# out_met_mobile_dir/isotope_metrics_site.csv

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(dplyr, readr, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

check_log <- tibble(
  check = character(),
  status = character(),
  detail = character()
)

failed_checks <- character()

add_check <- function(name, ok, detail) {
  check_log <<- bind_rows(
    check_log,
    tibble(
      check = name,
      status = ifelse(ok, "PASS", "FAIL"),
      detail = detail
    )
  )

  if (!ok) {
    failed_checks <<- c(failed_checks, paste0(name, ": ", detail))
  }
}

read_workflow_csv <- function(path, name) {
  ok <- file.exists(path)
  add_check(
    paste0("file_", name),
    ok,
    ifelse(ok, paste0("Found ", path), paste0("Missing ", path))
  )

  if (!ok) {
    return(tibble())
  }

  read_csv(path, show_col_types = FALSE)
}

require_columns <- function(df, cols, name) {
  missing <- setdiff(cols, names(df))
  add_check(
    paste0("columns_", name),
    length(missing) == 0,
    ifelse(
      length(missing) == 0,
      paste0(name, " has required columns"),
      paste0(name, " missing columns: ", paste(missing, collapse = ", "))
    )
  )
}

require_unique_rows <- function(df, cols, name) {
  if (!all(cols %in% names(df))) {
    return(invisible(NULL))
  }

  n_duplicate_groups <- df %>%
    filter(if_all(all_of(cols), ~ !is.na(.x))) %>%
    count(across(all_of(cols)), name = "n") %>%
    filter(n > 1) %>%
    nrow()

  add_check(
    paste0("unique_rows_", name),
    n_duplicate_groups == 0,
    ifelse(
      n_duplicate_groups == 0,
      paste0(name, " has one row per ", paste(cols, collapse = " and ")),
      paste0(name, " has ", n_duplicate_groups, " repeated row groups")
    )
  )
}

check_range <- function(df, col, low, high, name) {
  if (!(col %in% names(df))) {
    return(invisible(NULL))
  }

  x <- suppressWarnings(as.numeric(df[[col]]))
  bad <- sum(is.finite(x) & (x < low | x > high), na.rm = TRUE)
  add_check(
    paste0("range_", name, "_", col),
    bad == 0,
    ifelse(
      bad == 0,
      paste0(col, " is within ", low, " to ", high),
      paste0(col, " has ", bad, " values outside ", low, " to ", high)
    )
  )
}

max_abs_diff <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (!any(keep)) {
    return(NA_real_)
  }
  max(abs(x[keep] - y[keep]), na.rm = TRUE)
}

max_rel_diff <- function(x, y) {
  keep <- is.finite(x) & is.finite(y) & abs(y) > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  max(abs(x[keep] - y[keep]) / abs(y[keep]), na.rm = TRUE)
}

check_table_match <- function(df_left, df_right, by, cols, left_name, right_name, tolerance = 1e-8) {
  joined <- df_left %>%
    select(all_of(by), all_of(cols)) %>%
    inner_join(
      df_right %>% select(all_of(by), all_of(cols)),
      by = by,
      suffix = c("_left", "_right")
    )

  add_check(
    paste0("rows_", left_name, "_vs_", right_name),
    nrow(joined) > 0,
    paste0("Joined ", nrow(joined), " rows")
  )

  for (col in cols) {
    diff <- max_abs_diff(joined[[paste0(col, "_left")]], joined[[paste0(col, "_right")]])
    add_check(
      paste0(left_name, "_vs_", right_name, "_", col),
      is.finite(diff) && diff <= tolerance,
      paste0(col, " max absolute difference = ", signif(diff, 6))
    )
  }
}

met_q <- read_workflow_csv(
  file.path(OUT_MET_SUPPORT_DIR, "catchments_met_q.csv"),
  "catchments_met_q"
) %>%
  mutate(
    DATE = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y")),
    site = standardize_site_code(SITECODE)
  )

wb_daily <- read_workflow_csv(
  resolve_water_balance_daily_file(),
  "daily_water_balance"
) %>%
  mutate(
    DATE = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y")),
    site = standardize_site_code(SITECODE)
  )

discharge <- read_workflow_csv(
  file.path(DISCHARGE_DIR, "HF00402_v14.csv"),
  "discharge"
) %>%
  mutate(
    DATE = as.Date(DATE, "%m/%d/%Y"),
    SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK)),
    site = standardize_site_code(SITECODE)
  ) %>%
  filter(
    DATE >= as.Date(sprintf("%d-10-01", WY_START - 1)),
    DATE <= as.Date(sprintf("%d-09-30", WY_END)),
    !SITECODE %in% SITE_EXCLUDE_RAW
  ) %>%
  group_by(DATE, site) %>%
  summarise(MEAN_Q = sum(MEAN_Q, na.rm = TRUE), .groups = "drop")

drainage_area <- read_workflow_csv(
  resolve_drainage_area_file(),
  "drainage_area"
) %>%
  mutate(
    SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK)),
    site = standardize_site_code(SITECODE)
  )

dynamic_annual <- read_workflow_csv(
  file.path(OUT_MET_DYNAMIC_DIR, "storage_discharge_fdc_annual.csv"),
  "storage_discharge_fdc_annual"
) %>%
  mutate(site = as.character(site), year = as.integer(year))

master_annual <- read_workflow_csv(
  file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE),
  "master_annual"
) %>%
  mutate(site = as.character(site), year = as.integer(year))

master_site <- read_workflow_csv(
  file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE),
  "master_site"
) %>%
  mutate(site = as.character(site))

isotope_site <- read_workflow_csv(
  file.path(OUT_MET_MOBILE_DIR, "isotope_metrics_site.csv"),
  "isotope_metrics_site"
) %>%
  mutate(site = standardize_site_code(if ("site" %in% names(.)) site else SITECODE))

require_columns(met_q, c("DATE", "site", "T_C", "P_mm_d", "Q_mm_d", "RH_d_pct", "VPD_kPa"), "catchments_met_q")
require_columns(wb_daily, c("DATE", "site", "T_C", "P_mm_d", "Q_mm_d", "ET_mm_d"), "daily_water_balance")
require_columns(discharge, c("DATE", "site", "MEAN_Q"), "discharge")
require_columns(drainage_area, c("site", "DA_M2"), "drainage_area")
require_columns(dynamic_annual, c("site", "year", "SD", "FDC", "Q99", "Q50", "Q01"), "storage_discharge_fdc_annual")
require_columns(master_annual, c("site", "year", "RBI", "RCS", "FDC", "SD", "BF", "WB"), "master_annual")
require_columns(master_site, c("site", "RBI_mean", "RCS_mean", "FDC_mean", "SD_mean", "BF_mean", "WB_mean", "MTT", "Fyw", "DR"), "master_site")
require_columns(isotope_site, c("site", "MTT", "Fyw", "DR"), "isotope_metrics_site")

require_unique_rows(met_q, c("DATE", "site"), "catchments_met_q")
require_unique_rows(wb_daily, c("DATE", "site"), "daily_water_balance")
require_unique_rows(dynamic_annual, c("site", "year"), "storage_discharge_fdc_annual")
require_unique_rows(master_annual, c("site", "year"), "master_annual")
require_unique_rows(master_site, c("site"), "master_site")
require_unique_rows(isotope_site, c("site"), "isotope_metrics_site")

check_range(met_q, "P_mm_d", 0, Inf, "catchments_met_q")
check_range(met_q, "Q_mm_d", 0, Inf, "catchments_met_q")
check_range(met_q, "RH_d_pct", 0, 100, "catchments_met_q")
check_range(met_q, "T_C", -40, 60, "catchments_met_q")
check_range(met_q, "VPD_kPa", 0, Inf, "catchments_met_q")
check_range(wb_daily, "ET_mm_d", 0, Inf, "daily_water_balance")
check_range(master_site, "BF_mean", 0, 1, "master_site")
check_range(master_site, "Fyw", 0, 1, "master_site")

# confirm that the water balance table uses the same daily forcing values
check_table_match(
  met_q,
  wb_daily,
  by = c("DATE", "site"),
  cols = c("T_C", "P_mm_d", "Q_mm_d"),
  left_name = "met",
  right_name = "water_balance"
)

# rebuild Q_mm_d from source discharge and drainage area
q_formula <- discharge %>%
  left_join(drainage_area %>% select(site, DA_M2), by = "site") %>%
  filter(is.finite(DA_M2), DA_M2 > 0) %>%
  mutate(Q_mm_formula = MEAN_Q * 0.0283168 * 86400 / DA_M2 * 1000) %>%
  select(DATE, site, Q_mm_formula)

q_check <- met_q %>%
  select(DATE, site, Q_mm_d) %>%
  inner_join(q_formula, by = c("DATE", "site")) %>%
  filter(is.finite(Q_mm_d), is.finite(Q_mm_formula))

q_max_rel <- max_rel_diff(q_check$Q_mm_d, q_check$Q_mm_formula)
add_check(
  "q_mm_d_formula",
  nrow(q_check) > 0 && is.finite(q_max_rel) && q_max_rel <= 1e-6,
  paste0("Joined ", nrow(q_check), " rows, max relative difference = ", signif(q_max_rel, 6))
)

# confirm that derived tables feed master tables without changing values
check_table_match(
  dynamic_annual,
  master_annual,
  by = c("site", "year"),
  cols = c("SD", "FDC", "Q99", "Q50", "Q01"),
  left_name = "dynamic",
  right_name = "master_annual"
)

check_table_match(
  isotope_site,
  master_site,
  by = "site",
  cols = c("MTT", "Fyw", "DR"),
  left_name = "isotope",
  right_name = "master_site"
)

if (length(failed_checks) > 0) {
  stop(
    paste(
      c("Unit consistency checks failed:", paste0("- ", failed_checks)),
      collapse = "\n"
    ),
    call. = FALSE
  )
}
