# validate unit consistency across core workflow outputs.
# inputs:
# - out_met_support_dir/catchments_met_q.csv
# - resolve_water_balance_daily_file()
# - discharge_dir/hf00402_v14.csv
# - resolve_drainage_area_file()
# - out_met_dynamic_dir/storage_discharge_fdc_annual.csv
# - output_dir/master/master_annual.csv
# - output_dir/master/master_site.csv
# - out_met_mobile_dir/isotope_metrics_site.csv
# output (best effort):
# - out_met_support_dir/unit_consistency_report.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

find_repo_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(10)) {
    has_config <- file.exists(file.path(cur, "config.R"))
    has_metrics <- dir.exists(file.path(cur, "01_storage_calcs"))
    if (has_config && has_metrics) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
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
start_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_root <- find_repo_root(start_dir)

source(file.path(repo_root, "config.R"))

check_log <- tibble(
  check = character(),
  status = character(),
  detail = character()
)
fatal_msgs <- character()
warn_msgs <- character()

record_check <- function(name, ok, detail, fatal = TRUE) {
  status <- if (ok) "PASS" else if (fatal) "FAIL" else "WARN"
  check_log <<- bind_rows(
    check_log,
    tibble(check = as.character(name), status = status, detail = as.character(detail))
  )
  if (!ok) {
    if (fatal) {
      fatal_msgs <<- c(fatal_msgs, paste0(name, ": ", detail))
    } else {
      warn_msgs <<- c(warn_msgs, paste0(name, ": ", detail))
    }
  }
}

safe_read_csv <- function(path) {
  read_csv(path, show_col_types = FALSE)
}

check_required_columns <- function(df, required_cols, table_name) {
  missing <- setdiff(required_cols, names(df))
  ok <- length(missing) == 0
  detail <- if (ok) {
    paste0(table_name, " has required columns.")
  } else {
    paste0("Missing columns in ", table_name, ": ", paste(missing, collapse = ", "))
  }
  record_check(paste0("columns_", table_name), ok, detail, fatal = TRUE)
}

check_unique_key <- function(df, key_cols, table_name) {
  missing <- setdiff(key_cols, names(df))
  if (length(missing) > 0) {
    record_check(
      paste0("unique_key_", table_name),
      FALSE,
      paste0("Cannot test key; missing: ", paste(missing, collapse = ", ")),
      fatal = TRUE
    )
    return(invisible(NULL))
  }

  key_complete <- complete.cases(df[, key_cols, drop = FALSE])
  n_missing_keys <- sum(!key_complete)
  if (n_missing_keys > 0) {
    record_check(
      paste0("missing_key_rows_", table_name),
      FALSE,
      paste0(n_missing_keys, " rows have missing key values and were excluded from uniqueness check."),
      fatal = FALSE
    )
  }

  df_key <- df[key_complete, , drop = FALSE]
  dup_n <- df_key %>%
    count(across(all_of(key_cols)), name = "n") %>%
    filter(n > 1) %>%
    nrow()
  ok <- dup_n == 0
  detail <- if (ok) {
    paste0("Unique key ", paste(key_cols, collapse = "+"), " in ", table_name, ".")
  } else {
    paste0(dup_n, " duplicated key rows in ", table_name, ".")
  }
  record_check(paste0("unique_key_", table_name), ok, detail, fatal = TRUE)
}

max_diff_stats <- function(obs, ref, eps = 1e-12) {
  valid <- is.finite(obs) & is.finite(ref)
  if (!any(valid)) {
    return(list(n = 0L, max_abs = NA_real_, max_rel = NA_real_))
  }
  d <- abs(obs[valid] - ref[valid])
  rel <- d / pmax(abs(ref[valid]), eps)
  list(
    n = as.integer(sum(valid)),
    max_abs = max(d, na.rm = TRUE),
    max_rel = max(rel, na.rm = TRUE)
  )
}

check_nonnegative <- function(df, col, table_name) {
  if (!(col %in% names(df))) {
    record_check(
      paste0("nonnegative_", table_name, "_", col),
      FALSE,
      paste0("Column not found: ", col),
      fatal = TRUE
    )
    return(invisible(NULL))
  }
  x <- suppressWarnings(as.numeric(df[[col]]))
  bad <- sum(is.finite(x) & x < 0, na.rm = TRUE)
  ok <- bad == 0
  detail <- if (ok) {
    paste0(col, " is nonnegative in ", table_name, ".")
  } else {
    paste0(col, " has ", bad, " negative values in ", table_name, ".")
  }
  record_check(paste0("nonnegative_", table_name, "_", col), ok, detail, fatal = TRUE)
}

check_range <- function(df, col, low, high, table_name, fatal = TRUE) {
  if (!(col %in% names(df))) {
    record_check(
      paste0("range_", table_name, "_", col),
      FALSE,
      paste0("Column not found: ", col),
      fatal = fatal
    )
    return(invisible(NULL))
  }
  x <- suppressWarnings(as.numeric(df[[col]]))
  bad <- sum(is.finite(x) & (x < low | x > high), na.rm = TRUE)
  ok <- bad == 0
  x_min <- if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_
  x_max <- if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_
  detail <- if (ok) {
    paste0(col, " in ", table_name, " within [", low, ", ", high, "].")
  } else {
    paste0(
      col, " has ", bad, " values outside [", low, ", ", high, "] in ",
      table_name, " (min=", signif(x_min, 6), ", max=", signif(x_max, 6), ")."
    )
  }
  record_check(paste0("range_", table_name, "_", col), ok, detail, fatal = fatal)
}

check_site_median_ratio <- function(df, value_col, table_name, ratio_threshold = 500, min_n = 30) {
  if (!(value_col %in% names(df)) || !("site" %in% names(df))) {
    record_check(
      paste0("site_ratio_", table_name, "_", value_col),
      FALSE,
      "Missing site/value column for ratio check.",
      fatal = TRUE
    )
    return(invisible(NULL))
  }
  by_site <- df %>%
    transmute(site = as.character(site), value = abs(.data[[value_col]])) %>%
    filter(is.finite(value), value > 0) %>%
    group_by(site) %>%
    summarise(
      n = n(),
      median_abs = median(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n >= min_n, is.finite(median_abs), median_abs > 0)

  if (nrow(by_site) < 2) {
    record_check(
      paste0("site_ratio_", table_name, "_", value_col),
      TRUE,
      paste0("Insufficient site coverage for ", value_col, " ratio check."),
      fatal = TRUE
    )
    return(invisible(NULL))
  }

  ratio <- max(by_site$median_abs, na.rm = TRUE) / min(by_site$median_abs, na.rm = TRUE)
  ok <- is.finite(ratio) && ratio <= ratio_threshold
  detail <- if (ok) {
    paste0(value_col, " site-median max/min ratio=", signif(ratio, 6), " in ", table_name, ".")
  } else {
    worst <- by_site %>% arrange(median_abs)
    paste0(
      value_col, " site-median ratio=", signif(ratio, 6), " exceeds threshold ",
      ratio_threshold, " in ", table_name,
      " (min site=", worst$site[1], ", max site=", worst$site[nrow(worst)], ")."
    )
  }
  record_check(paste0("site_ratio_", table_name, "_", value_col), ok, detail, fatal = TRUE)
}

master_dir <- file.path(OUTPUT_DIR, "master")
wb_daily_file <- tryCatch(resolve_water_balance_daily_file(), error = function(e) NA_character_)

required_paths <- c(
  met_q = file.path(OUT_MET_SUPPORT_DIR, "catchments_met_q.csv"),
  wb_daily = wb_daily_file,
  discharge = file.path(DISCHARGE_DIR, "HF00402_v14.csv"),
  drainage_area = resolve_drainage_area_file(),
  dynamic_annual = file.path(OUT_MET_DYNAMIC_DIR, "storage_discharge_fdc_annual.csv"),
  master_annual = file.path(master_dir, MASTER_ANNUAL_FILE),
  master_site = file.path(master_dir, MASTER_SITE_FILE),
  isotope_site = file.path(OUT_MET_MOBILE_DIR, "isotope_metrics_site.csv")
)

for (nm in names(required_paths)) {
  p <- required_paths[[nm]]
  ok <- is.character(p) && nzchar(p) && file.exists(p)
  detail <- if (ok) paste0("Found: ", p) else paste0("Missing required file: ", p)
  record_check(paste0("file_", nm), ok, detail, fatal = TRUE)
}

if (length(fatal_msgs) > 0) {
  stop(
    paste(
      c("Unit consistency checks FAILED before data load:", paste0("- ", fatal_msgs)),
      collapse = "\n"
    )
  )
}

met_q <- safe_read_csv(required_paths[["met_q"]]) %>%
  mutate(
    DATE = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y")),
    site = standardize_site_code(SITECODE)
  )

wb_daily <- safe_read_csv(required_paths[["wb_daily"]]) %>%
  mutate(
    DATE = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y")),
    site = standardize_site_code(SITECODE)
  )

discharge <- safe_read_csv(required_paths[["discharge"]]) %>%
  mutate(
    DATE = as.Date(DATE, "%m/%d/%Y"),
    SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK))
  ) %>%
  filter(
    DATE >= as.Date(sprintf("%d-10-01", WY_START - 1)),
    DATE <= as.Date(sprintf("%d-09-30", WY_END)),
    !SITECODE %in% SITE_EXCLUDE_RAW
  ) %>%
  group_by(DATE, SITECODE) %>%
  summarise(MEAN_Q = sum(MEAN_Q, na.rm = TRUE), .groups = "drop") %>%
  mutate(site = standardize_site_code(SITECODE))

da_df <- safe_read_csv(required_paths[["drainage_area"]]) %>%
  mutate(
    SITECODE = recode(SITECODE, !!!as.list(SITECODE_RECODE_TO_GSMACK)),
    site = standardize_site_code(SITECODE)
  )

dynamic_annual <- safe_read_csv(required_paths[["dynamic_annual"]]) %>%
  mutate(
    DATE = NA,
    site = as.character(site),
    year = as.integer(year)
  )

master_annual <- safe_read_csv(required_paths[["master_annual"]]) %>%
  mutate(
    site = as.character(site),
    year = as.integer(year)
  )

master_site <- safe_read_csv(required_paths[["master_site"]]) %>%
  mutate(site = as.character(site))

isotope_site <- safe_read_csv(required_paths[["isotope_site"]])
isotope_site_col <- if ("site" %in% names(isotope_site)) {
  "site"
} else if ("SITECODE" %in% names(isotope_site)) {
  "SITECODE"
} else {
  NA_character_
}
if (is.na(isotope_site_col)) {
  record_check("columns_isotope_site_site_id", FALSE, "isotope_metrics_site missing site identifier column.", fatal = TRUE)
  isotope_site <- isotope_site %>% mutate(site = NA_character_)
} else {
  isotope_site <- isotope_site %>% mutate(site = standardize_site_code(.data[[isotope_site_col]]))
}

# required schema checks.
check_required_columns(
  met_q,
  c("DATE", "SITECODE", "site", "T_C", "P_mm_d", "RH_d_pct", "NR_Wm2_d", "VPD_kPa", "Q_mm_d"),
  "catchments_met_q"
)
check_required_columns(
  wb_daily,
  c("DATE", "SITECODE", "site", "T_C", "P_mm_d", "Q_mm_d", "ET_mm_d"),
  "daily_water_balance"
)
check_required_columns(
  dynamic_annual,
  c("site", "year", "SD", "FDC", "Q99", "Q50", "Q01"),
  "storage_discharge_fdc_annual"
)
check_required_columns(
  master_annual,
  c("site", "year", "RBI", "RCS", "FDC", "SD", "CHS", "WB"),
  "master_annual"
)
check_required_columns(
  master_site,
  c("site", "RBI_mean", "RCS_mean", "FDC_mean", "SD_mean", "CHS_mean", "WB_mean"),
  "master_site"
)
check_required_columns(
  isotope_site,
  c("site", "MTT1", "MTT2", "Fyw", "DR"),
  "isotope_metrics_site"
)

# uniqueness checks.
check_unique_key(met_q, c("DATE", "site"), "catchments_met_q")
check_unique_key(wb_daily, c("DATE", "site"), "daily_water_balance")
check_unique_key(dynamic_annual, c("site", "year"), "storage_discharge_fdc_annual")
check_unique_key(master_annual, c("site", "year"), "master_annual")
check_unique_key(master_site, c("site"), "master_site")
check_unique_key(isotope_site, c("site"), "isotope_metrics_site")

# core range checks.
check_nonnegative(met_q, "P_mm_d", "catchments_met_q")
check_nonnegative(met_q, "Q_mm_d", "catchments_met_q")
check_range(met_q, "RH_d_pct", low = 0, high = 100, table_name = "catchments_met_q", fatal = TRUE)
check_range(met_q, "T_C", low = -40, high = 60, table_name = "catchments_met_q", fatal = TRUE)
check_range(met_q, "NR_Wm2_d", low = -300, high = 1000, table_name = "catchments_met_q", fatal = FALSE)
check_nonnegative(met_q, "VPD_kPa", "catchments_met_q")

check_nonnegative(wb_daily, "P_mm_d", "daily_water_balance")
check_nonnegative(wb_daily, "Q_mm_d", "daily_water_balance")
check_nonnegative(wb_daily, "ET_mm_d", "daily_water_balance")
check_range(wb_daily, "T_C", low = -40, high = 60, table_name = "daily_water_balance", fatal = TRUE)

check_range(master_site, "CHS_mean", low = 0, high = 1, table_name = "master_site", fatal = FALSE)
check_range(master_site, "Fyw", low = 0, high = 1, table_name = "master_site", fatal = FALSE)
check_nonnegative(master_site, "RBI_mean", "master_site")

# compare shared columns between met support and water-balance daily.
common_cols <- intersect(c("T_C", "P_mm_d", "Q_mm_d"), intersect(names(met_q), names(wb_daily)))
met_wb_cmp <- met_q %>%
  select(DATE, site, all_of(common_cols)) %>%
  inner_join(
    wb_daily %>% select(DATE, site, all_of(common_cols)),
    by = c("DATE", "site"),
    suffix = c("_met", "_wb")
  )

record_check(
  "rows_met_vs_wb_join",
  nrow(met_wb_cmp) > 0,
  paste0("Joined ", nrow(met_wb_cmp), " overlapping DATE/site rows."),
  fatal = TRUE
)

for (col in common_cols) {
  s <- max_diff_stats(met_wb_cmp[[paste0(col, "_met")]], met_wb_cmp[[paste0(col, "_wb")]])
  ok <- s$n > 0 && is.finite(s$max_abs) && s$max_abs <= 1e-8
  detail <- if (ok) {
    paste0(col, " consistent between met and wb daily (n=", s$n, ", max_abs=", signif(s$max_abs, 6), ").")
  } else {
    paste0(col, " mismatch between met and wb daily (n=", s$n, ", max_abs=", signif(s$max_abs, 6), ").")
  }
  record_check(paste0("met_vs_wb_", col), ok, detail, fatal = TRUE)
}

# independent q_mm_d formula check against discharge + drainage area.
q_expected <- discharge %>%
  left_join(da_df %>% select(site, DA_M2), by = "site") %>%
  filter(is.finite(DA_M2), DA_M2 > 0) %>%
  mutate(Q_mm_expected = MEAN_Q * 0.0283168 * 86400 / DA_M2 * 1000) %>%
  select(DATE, site, Q_mm_expected)

q_cmp <- met_q %>%
  select(DATE, site, Q_mm_d) %>%
  inner_join(q_expected, by = c("DATE", "site")) %>%
  filter(is.finite(Q_mm_d), is.finite(Q_mm_expected), Q_mm_expected > 0)

record_check(
  "rows_q_formula_join",
  nrow(q_cmp) > 0,
  paste0("Joined ", nrow(q_cmp), " Q rows for independent formula check."),
  fatal = TRUE
)

q_stats <- max_diff_stats(q_cmp$Q_mm_d, q_cmp$Q_mm_expected)
q_row_ok <- q_stats$n > 0 && is.finite(q_stats$max_rel) && q_stats$max_rel <= 1e-6
record_check(
  "q_formula_row_match",
  q_row_ok,
  paste0(
    "Q_mm_d vs expected: n=", q_stats$n,
    ", max_abs=", signif(q_stats$max_abs, 6),
    ", max_rel=", signif(q_stats$max_rel, 6), "."
  ),
  fatal = TRUE
)

q_site_ratio <- q_cmp %>%
  mutate(ratio = Q_mm_d / Q_mm_expected) %>%
  group_by(site) %>%
  summarise(
    n = n(),
    median_ratio = median(ratio, na.rm = TRUE),
    p05 = quantile(ratio, 0.05, na.rm = TRUE),
    p95 = quantile(ratio, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

q_bad_sites <- q_site_ratio %>%
  filter(!is.finite(median_ratio) | abs(log10(median_ratio)) > 0.03)

record_check(
  "q_formula_site_ratio",
  nrow(q_bad_sites) == 0,
  if (nrow(q_bad_sites) == 0) {
    "All site median Q ratio checks are within tolerance."
  } else {
    paste0(
      "Sites outside tolerance: ",
      paste(q_bad_sites$site, collapse = ", ")
    )
  },
  fatal = TRUE
)

# potential cross-site unit mismatch checks for mm/day variables.
check_site_median_ratio(met_q, "P_mm_d", "catchments_met_q", ratio_threshold = 500)
check_site_median_ratio(met_q, "Q_mm_d", "catchments_met_q", ratio_threshold = 500)
check_site_median_ratio(wb_daily, "ET_mm_d", "daily_water_balance", ratio_threshold = 500)

# dynamic annual table should match master_annual for shared fields.
dyn_cmp <- dynamic_annual %>%
  select(site, year, SD, FDC, Q99, Q50, Q01) %>%
  inner_join(
    master_annual %>% select(site, year, SD, FDC, Q99, Q50, Q01),
    by = c("site", "year"),
    suffix = c("_dyn", "_master")
  )
record_check(
  "rows_dynamic_vs_master_annual_join",
  nrow(dyn_cmp) > 0,
  paste0("Joined ", nrow(dyn_cmp), " dynamic annual rows to master_annual."),
  fatal = TRUE
)

for (col in c("SD", "FDC", "Q99", "Q50", "Q01")) {
  s <- max_diff_stats(dyn_cmp[[paste0(col, "_dyn")]], dyn_cmp[[paste0(col, "_master")]])
  ok <- s$n > 0 && is.finite(s$max_abs) && s$max_abs <= 1e-8
  detail <- if (ok) {
    paste0(col, " matches between dynamic annual and master_annual (max_abs=", signif(s$max_abs, 6), ").")
  } else {
    paste0(col, " mismatch between dynamic annual and master_annual (max_abs=", signif(s$max_abs, 6), ").")
  }
  record_check(paste0("dynamic_vs_master_annual_", col), ok, detail, fatal = TRUE)
}

# master_site *_mean columns should equal means from master_annual where base columns exist.
mean_cols <- grep("_mean$", names(master_site), value = TRUE)
mean_pairs <- tibble(mean_col = mean_cols) %>%
  mutate(base_col = sub("_mean$", "", mean_col)) %>%
  filter(base_col %in% names(master_annual))

if (nrow(mean_pairs) == 0) {
  record_check(
    "master_site_mean_columns",
    FALSE,
    "No *_mean columns in master_site matched master_annual columns.",
    fatal = TRUE
  )
} else {
  annual_bases <- unique(mean_pairs$base_col)
  annual_means <- master_annual %>%
    group_by(site) %>%
    summarise(
      across(
        all_of(annual_bases),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  for (i in seq_len(nrow(mean_pairs))) {
    mcol <- mean_pairs$mean_col[i]
    bcol <- mean_pairs$base_col[i]
    cmp <- master_site %>%
      select(site, site_val = all_of(mcol)) %>%
      inner_join(annual_means %>% select(site, annual_val = all_of(bcol)), by = "site")
    s <- max_diff_stats(cmp$site_val, cmp$annual_val)
    ok <- s$n > 0 && is.finite(s$max_abs) && s$max_abs <= 1e-8
    detail <- if (ok) {
      paste0(mcol, " matches mean(", bcol, ") from master_annual (max_abs=", signif(s$max_abs, 6), ").")
    } else {
      paste0(mcol, " mismatch vs mean(", bcol, ") from master_annual (max_abs=", signif(s$max_abs, 6), ").")
    }
    record_check(paste0("master_site_mean_", mcol), ok, detail, fatal = TRUE)
  }
}

# isotope site metrics should pass through to master_site unchanged.
iso_cmp <- isotope_site %>%
  select(site, MTT1, MTT2, Fyw, DR) %>%
  inner_join(master_site %>% select(site, MTT1, MTT2, Fyw, DR), by = "site", suffix = c("_iso", "_site"))

record_check(
  "rows_isotope_vs_master_site_join",
  nrow(iso_cmp) > 0,
  paste0("Joined ", nrow(iso_cmp), " isotope site rows to master_site."),
  fatal = TRUE
)

for (col in c("MTT1", "MTT2", "Fyw", "DR")) {
  s <- max_diff_stats(iso_cmp[[paste0(col, "_iso")]], iso_cmp[[paste0(col, "_site")]])
  ok <- s$n > 0 && is.finite(s$max_abs) && s$max_abs <= 1e-8
  detail <- if (ok) {
    paste0(col, " matches between isotope_metrics_site and master_site.")
  } else {
    paste0(col, " mismatch between isotope_metrics_site and master_site (max_abs=", signif(s$max_abs, 6), ").")
  }
  record_check(paste0("isotope_vs_master_site_", col), ok, detail, fatal = TRUE)
}

# write detailed report when possible.
report_path <- file.path(OUT_MET_SUPPORT_DIR, "unit_consistency_report.csv")
tryCatch(
  write_csv(check_log, report_path),
  error = function(e) {
    record_check(
      "write_unit_consistency_report",
      FALSE,
      paste0("Could not write report to ", report_path, ": ", conditionMessage(e)),
      fatal = FALSE
    )
  }
)

if (length(warn_msgs) > 0) {
  warning(
    paste(
      c("Unit consistency checks produced warnings:", paste0("- ", warn_msgs)),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

if (length(fatal_msgs) > 0) {
  stop(
    paste(
      c("Unit consistency checks FAILED:", paste0("- ", fatal_msgs)),
      collapse = "\n"
    )
  )
}

message("Unit consistency checks passed: ", sum(check_log$status == "PASS"), " checks.")
