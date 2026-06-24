# shared workflow checks, folders, and site utilities

# inputs:
# none

# outputs:
# workflow functions loaded by config.R

# author: Sidney Bush
# date: 2026-02-13

# water year from date
# water years start on October 1
get_water_year <- function(date) {
  ifelse(
    lubridate::month(date) >= 10,
    lubridate::year(date) + 1,
    lubridate::year(date)
  )
}

# day of water year from date
# October 1 is day 1
get_water_year_day <- function(date) {
  wy <- get_water_year(date)
  wy_start <- as.Date(paste0(wy - 1, "-10-01"))
  as.numeric(date - wy_start) + 1
}

# standard site codes used in outputs
# source files use several naming systems for the same sites
standardize_site_code <- function(site_code) {
  site_code <- trimws(site_code)
  dplyr::case_when(
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
    grepl("^GSWS[0-9]+$", site_code) ~ gsub("^GSWS", "WS", site_code),
    TRUE ~ site_code
  )
}

# full storage metric names used in tables and figure labels
STORAGE_METRIC_FULL_NAMES <- c(
  "RBI" = "Richards-Baker Index",
  "RCS" = "Recession Curve Slope",
  "FDC" = "Flow Duration Curve Slope",
  "SD" = "Storage-Discharge",
  "WB" = "Water-balance deficit",
  "BF" = "Baseflow Fraction",
  "DR" = "Damping Ratio",
  "Fyw" = "Young Water Fraction",
  "MTT" = "Mean Transit Time"
)

# short metric names used when space is limited
METRIC_ABBREV_DISPLAY <- c(
  "RBI" = "RBI",
  "RCS" = "RCS",
  "FDC" = "FDC",
  "SD" = "SD",
  "WB" = "WB",
  "BF" = "BF",
  "DR" = "DR",
  "Fyw" = "Fyw",
  "MTT" = "MTT"
)

# return a short metric name
label_metric_abbrev <- function(x) {
  x_chr <- as.character(x)
  out <- unname(METRIC_ABBREV_DISPLAY[x_chr])
  miss <- is.na(out)
  out[miss] <- x_chr[miss]
  out
}

# return a full storage metric label
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

# catchment predictor names used in model figures and tables
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

# return a clean catchment predictor label
label_catchment_predictor <- function(x) {
  x_chr <- as.character(x)
  out <- unname(CATCHMENT_PREDICTOR_LABELS[x_chr])
  miss <- is.na(out)
  out[miss] <- x_chr[miss]
  out
}

# return clean labels for comma separated predictor lists
label_catchment_predictor_list <- function(x) {
  x_chr <- as.character(x)
  vapply(
    x_chr,
    function(val) {
      if (is.na(val) || !nzchar(val)) {
        return(val)
      }
      parts <- trimws(strsplit(val, "[[:space:]]*[,;][[:space:]]*")[[1]])
      parts <- parts[nzchar(parts)]
      if (length(parts) == 0) {
        return(val)
      }
      paste(label_catchment_predictor(parts), collapse = ", ")
    },
    FUN.VALUE = character(1)
  )
}

# standardized site codes to leave out of analysis tables
SITE_EXCLUDE_STANDARD <- unique(standardize_site_code(SITE_EXCLUDE_RAW))

# files used by several scripts
resolve_water_balance_daily_file <- function() {
  file.path(
    OUT_MET_SUPPORT_DIR,
    "daily_water_balance_et_hamon_zhang_coeff_interp.csv"
  )
}

resolve_drainage_area_file <- function() {
  file.path(CATCHMENT_CHARACTERISTICS_DIR, "drainage_area.csv")
}

resolve_catchment_characteristics_file <- function() {
  file.path(CATCHMENT_CHARACTERISTICS_DIR, "catchment_char.csv")
}

# check input columns before running the workflow
check_inputs <- function() {
  # require readr for lightweight column checks
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "Missing R package: readr. Run Rscript install_packages.R once, ",
      "then rerun Rscript run_all.R.",
      call. = FALSE
    )
  }

  check_columns <- function(path, cols) {
    file_cols <- names(suppressMessages(readr::read_csv(
      path,
      n_max = 0,
      show_col_types = FALSE
    )))
    missing_cols <- setdiff(cols, file_cols)

    # missing columns usually mean a source file changed
    if (length(missing_cols) > 0) {
      stop(
        "Missing columns in ",
        path,
        ": ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  check_columns(file.path(DISCHARGE_DIR, "HF00402_v14.csv"), c("DATE", "SITECODE", "MEAN_Q"))
  check_columns(file.path(EC_DIR, "CF00201_v7.csv"), c("DATE_TIME", "SITECODE", "CA"))
  check_columns(file.path(STREAM_TEMP_DIR, "HT00451_v10.txt"), c("DATE_TIME", "SITECODE", "WATERTEMP_MEAN"))
  check_columns(file.path(CATCHMENT_CHARACTERISTICS_DIR, "catchment_char.csv"), c("Site"))

  invisible(TRUE)
}

# make folders used by the full workflow
make_output_dirs <- function() {
  dirs <- c(
    OUT_METRICS_DIR,
    OUT_MET_DYNAMIC_DIR,
    OUT_MET_MOBILE_DIR,
    OUT_MET_EXTENDED_DIR,
    OUT_MET_ECO_DIR,
    OUT_MET_SUPPORT_DIR,
    file.path(OUTPUT_DIR, "master"),
    FIGS_TABLES_PUB_DIR,
    MS_MAIN_DIR,
    MS_SUPP_DIR,
    MS_FIG_MAIN_DIR,
    MS_FIG_SUPP_DIR,
    MS_FIG_MAIN_PDF_DIR,
    MS_FIG_SUPP_PDF_DIR,
    MS_FIG_MAIN_TIFF_DIR,
    MS_FIG_SUPP_TIFF_DIR,
    MS_TABLES_MAIN_DIR,
    MS_TABLES_SUPP_DIR,
    OUT_STATS_DIR,
    OUT_STATS_ANOVA_DIR,
    OUT_STATS_PCA_DIR,
    OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR,
    OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR
  )

  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  invisible(dirs)
}

# make sure the workflow wrote the files used by the manuscript
verify_outputs <- function() {
  # require readr for manuscript output checks
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "Missing R package: readr. Run Rscript install_packages.R once, ",
      "then rerun Rscript run_all.R.",
      call. = FALSE
    )
  }

  master_dir <- file.path(OUTPUT_DIR, "master")
  master_annual_file <- file.path(master_dir, MASTER_ANNUAL_FILE)
  letters_path <- file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv")

  figure_files <- function(bases, pdf_dir, tiff_dir) {
    unlist(
      lapply(bases, function(base) {
        c(
          file.path(pdf_dir, paste0(base, ".pdf")),
          file.path(tiff_dir, paste0(base, ".tiff"))
        )
      }),
      use.names = FALSE
    )
  }

  main_figure_bases <- c(
    "Fig2_dynamic_storage_pca",
    "Fig3_mobile_storage",
    "Fig4_catchment_controls",
    "Fig5_ecological_response_models",
    "Fig6_observed_predicted_ecological_responses",
    "Fig7_dynamic_mobile_framework"
  )

  supp_figure_bases <- c(
    "FigS1_met_context",
    "FigS2_dynamic_storage_corr",
    "FigS3_mobile_storage_corr",
    "FigS4_dynamic_mobile_corr"
  )

  si_table_files <- file.path(
    MS_TABLES_SUPP_DIR,
    paste0(
      c(
        "TableS1_data_periods_of_record",
        "TableS2_catchment_physiography_land_use",
        "TableS3_catchment_geology_landslides",
        "TableS4_met_station_record_summary",
        "TableS5_met_station_assignments",
        "TableS6_isotope_metrics",
        "TableS7_MTT_sensitivity",
        "TableS8_catchment_char_storage_mlr_model_stats",
        "TableS9_catchment_alt_models_unique_deltaAICc_le2_BF",
        "TableS10_mlr_model_diagnostics",
        "TableS11_storage_eco_response_mlr_model_stats",
        "TableS12_eco_alt_models_unique_deltaAICc_le2_BF"
      ),
      ".csv"
    )
  )

  required_outputs <- c(
    master_annual_file,
    letters_path,
    figure_files(main_figure_bases, MS_FIG_MAIN_PDF_DIR, MS_FIG_MAIN_TIFF_DIR),
    figure_files(supp_figure_bases, MS_FIG_SUPP_PDF_DIR, MS_FIG_SUPP_TIFF_DIR),
    si_table_files
  )

  missing_outputs <- required_outputs[!file.exists(required_outputs)]

  # missing manuscript outputs should stop the run
  if (length(missing_outputs) > 0) {
    stop(
      paste0(
        "Missing output file(s):\n- ",
        paste(missing_outputs, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }

  for (table_file in si_table_files) {
    table_df <- readr::read_csv(
      table_file,
      show_col_types = FALSE,
      name_repair = "minimal"
    )
    bad_names <- names(table_df) == "" | is.na(names(table_df))

    # blank headings make the supporting information tables hard to read
    if (any(bad_names)) {
      stop("Table has missing column heading(s): ", table_file, call. = FALSE)
    }
    duplicated_names <- names(table_df)[duplicated(names(table_df))]

    # repeated headings make the supporting information tables ambiguous
    if (length(duplicated_names) > 0) {
      stop(
        "Table has duplicated column heading(s): ",
        table_file,
        " (",
        paste(unique(duplicated_names), collapse = ", "),
        ")",
        call. = FALSE
      )
    }
  }

  annual <- readr::read_csv(master_annual_file, show_col_types = FALSE)

  # check the master annual table before downstream manuscript checks
  if (!all(c("site", "year") %in% names(annual))) {
    stop("master_annual is missing columns: site/year", call. = FALSE)
  }
  if (nrow(annual) == 0) {
    stop("master_annual has zero rows", call. = FALSE)
  }
  if (any(!is.finite(annual$year))) {
    stop("master_annual has non-finite years", call. = FALSE)
  }

  unknown_sites <- setdiff(unique(as.character(annual$site)), SITE_ORDER_HYDROMETRIC)

  # unknown site codes usually mean a source name was not standardized
  if (length(unknown_sites) > 0) {
    stop(
      "Unknown site code(s) in master_annual: ",
      paste(unknown_sites, collapse = ", "),
      call. = FALSE
    )
  }

  letters_df <- readr::read_csv(letters_path, show_col_types = FALSE)

  # check Tukey grouping columns before checking site order
  if (!all(c("metric", "site") %in% names(letters_df))) {
    stop("tukey_group_letters is missing columns: metric/site", call. = FALSE)
  }

  metric_groups <- split(as.character(letters_df$site), letters_df$metric)
  for (metric_name in names(metric_groups)) {
    site_values <- metric_groups[[metric_name]]
    allowed <- SITE_ORDER_HYDROMETRIC[SITE_ORDER_HYDROMETRIC %in% site_values]

    # keep Tukey letters in the manuscript site order
    if (!identical(site_values, allowed)) {
      stop(
        "Site order mismatch in tukey_group_letters for metric: ",
        metric_name,
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}
