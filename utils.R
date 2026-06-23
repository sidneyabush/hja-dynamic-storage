# shared functions, plot settings, and labels

# shared font sizes and export settings
FIG_BASE_SIZE <- 18
FIG_AXIS_TEXT_SIZE <- 16
FIG_AXIS_TITLE_SIZE <- 18
FIG_STRIP_TEXT_SIZE <- 16
FIG_ANNOT_TEXT_SIZE <- 5
FIG_TILE_TEXT_SIZE <- 6
FIG_POINT_SIZE_SMALL <- 1.5
FIG_POINT_SIZE_LARGE <- 3.0
FIG_WIDTH_SCALE <- 1.35
FIG_HEIGHT_SCALE <- 1.35
FIG_PREVIEW_DPI <- 300
FIG_PRODUCTION_DPI <- 300

# shared label spacing settings
FIG_LABEL_CLIP <- "off"
FIG_LABEL_PLOT_MARGIN_PT <- 18

# site colors used in figures
SITE_COLORS <- c(
  "WS09" = "#882255",
  "WS10" = "#AA4499",
  "WS01" = "#CC6677",
  "Look" = "#DDCC77",
  "WS02" = "#999933",
  "WS03" = "#117733",
  "WS06" = "#44AA99",
  "WS07" = "#88CCEE",
  "WS08" = "#6699CC",
  "Mack" = "#332288"
)

# plot theme
# individual figure scripts add any figure specific styling
theme_pub <- function(base_size = FIG_BASE_SIZE) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(hjust = 0),
      strip.text.y = ggplot2::element_text(hjust = 0)
    )
}

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

# return a file path and stop if the file is missing
resolve_file <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing ", label, ": ", path, call. = FALSE)
  }
  path
}

# files used by several scripts
resolve_water_balance_daily_file <- function() {
  resolve_file(
    file.path(
      OUT_MET_SUPPORT_DIR,
      "daily_water_balance_et_hamon_zhang_coeff_interp.csv"
    ),
    "water balance daily file"
  )
}

resolve_drainage_area_file <- function() {
  resolve_file(
    file.path(CATCHMENT_CHARACTERISTICS_DIR, "drainage_area.csv"),
    "drainage area file"
  )
}

resolve_catchment_characteristics_file <- function() {
  resolve_file(
    file.path(CATCHMENT_CHARACTERISTICS_DIR, "catchment_char.csv"),
    "catchment characteristics file"
  )
}

# check input files before running the workflow
check_inputs <- function() {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "Missing R package: readr. Run Rscript install_packages.R once, ",
      "then rerun Rscript run_all.R.",
      call. = FALSE
    )
  }

  input_files <- c(
    file.path(MET_DIR, "Temperature_original_&_filled_1979_2023_v2.csv"),
    file.path(MET_DIR, "Precipitation_original_&_filled_1979_2023.csv"),
    file.path(MET_DIR, "SWE_original_&_filled_1997_2023_v5.csv"),
    file.path(MET_DIR, "MS00102_v9.csv"),
    file.path(MET_DIR, "MS05025_v3.csv"),
    file.path(MET_DIR, "MS00403_v2.csv"),
    file.path(DISCHARGE_DIR, "HF00402_v14.csv"),
    file.path(CATCHMENT_CHARACTERISTICS_DIR, "drainage_area.csv"),
    file.path(CATCHMENT_CHARACTERISTICS_DIR, "catchment_char.csv"),
    file.path(ISOTOPE_DIR, "MTT_FYW.csv"),
    file.path(ISOTOPE_DIR, "DampingRatios_2025-07-07.csv"),
    file.path(STREAM_TEMP_DIR, "HT00451_v10.txt"),
    file.path(EC_DIR, "CF01201_v4.txt"),
    file.path(EC_DIR, "CF00201_v7.csv")
  )

  missing_files <- input_files[!file.exists(input_files)]
  if (length(missing_files) > 0) {
    stop(
      paste0(
        "Missing input file(s):\n- ",
        paste(missing_files, collapse = "\n- ")
      ),
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
  check_columns(file.path(EC_DIR, "CF01201_v4.txt"), c("DATE_TIME", "SITECODE", "EC_INST"))
  check_columns(file.path(EC_DIR, "CF00201_v7.csv"), c("DATE_TIME", "SITECODE", "COND", "CA"))
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
    MS_MATERIALS_DIR,
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
    OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
    EXPLORATORY_PLOTS_DIR,
    UNIFIED_FRAMEWORK_DIR,
    EXPLORATORY_ET_METHODS_DIR
  )

  if (isTRUE(WRITE_TABLE_OUTPUTS)) {
    dirs <- c(dirs, OUT_TABLES_DIR, OUT_TABLES_MLR_DIR)
  }

  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  invisible(dirs)
}

# check that the workflow wrote the files used by the manuscript
verify_outputs <- function() {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "Missing R package: readr. Run Rscript install_packages.R once, ",
      "then rerun Rscript run_all.R.",
      call. = FALSE
    )
  }

  master_dir <- file.path(OUTPUT_DIR, "master")

  required_outputs <- c(
    file.path(master_dir, MASTER_ANNUAL_FILE),
    file.path(master_dir, MASTER_SITE_FILE),
    file.path(OUT_STATS_ANOVA_DIR, "anova_results.csv"),
    file.path(OUT_STATS_ANOVA_DIR, "tukey_hsd_results.csv"),
    file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv"),
    file.path(OUT_STATS_ANOVA_DIR, "storage_metrics_summary_stats_by_site.csv"),
    file.path(OUT_STATS_PCA_DIR, "pca_loadings.csv"),
    file.path(OUT_STATS_PCA_DIR, "pca_variance_explained.csv"),
    file.path(OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR, "catchment_char_storage_mlr_summary.csv"),
    file.path(OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR, "storage_eco_response_mlr_summary.csv"),
    file.path(MS_FIG_MAIN_DIR, "Fig2_dynamic_storage_pca.png"),
    file.path(MS_FIG_MAIN_PDF_DIR, "Fig2_dynamic_storage_pca.pdf"),
    file.path(MS_FIG_MAIN_TIFF_DIR, "Fig2_dynamic_storage_pca.tiff"),
    file.path(MS_FIG_MAIN_DIR, "Fig3_mobile_storage.png"),
    file.path(MS_FIG_MAIN_PDF_DIR, "Fig3_mobile_storage.pdf"),
    file.path(MS_FIG_MAIN_TIFF_DIR, "Fig3_mobile_storage.tiff"),
    file.path(MS_FIG_MAIN_DIR, "Fig4_catchment_controls.png"),
    file.path(MS_FIG_MAIN_PDF_DIR, "Fig4_catchment_controls.pdf"),
    file.path(MS_FIG_MAIN_TIFF_DIR, "Fig4_catchment_controls.tiff"),
    file.path(MS_FIG_MAIN_DIR, "Fig5_ecological_response_models.png"),
    file.path(MS_FIG_MAIN_PDF_DIR, "Fig5_ecological_response_models.pdf"),
    file.path(MS_FIG_MAIN_TIFF_DIR, "Fig5_ecological_response_models.tiff"),
    file.path(MS_FIG_MAIN_DIR, "Fig6_observed_predicted_ecological_responses.png"),
    file.path(MS_FIG_MAIN_PDF_DIR, "Fig6_observed_predicted_ecological_responses.pdf"),
    file.path(MS_FIG_MAIN_TIFF_DIR, "Fig6_observed_predicted_ecological_responses.tiff"),
    file.path(MS_FIG_MAIN_DIR, "Fig7_dynamic_mobile_framework.png"),
    file.path(MS_FIG_MAIN_PDF_DIR, "Fig7_dynamic_mobile_framework.pdf"),
    file.path(MS_FIG_MAIN_TIFF_DIR, "Fig7_dynamic_mobile_framework.tiff"),
    file.path(MS_FIG_SUPP_DIR, "FigS1_met_context.png"),
    file.path(MS_FIG_SUPP_PDF_DIR, "FigS1_met_context.pdf"),
    file.path(MS_FIG_SUPP_TIFF_DIR, "FigS1_met_context.tiff"),
    file.path(MS_FIG_SUPP_DIR, "FigS2_dynamic_storage_corr.png"),
    file.path(MS_FIG_SUPP_PDF_DIR, "FigS2_dynamic_storage_corr.pdf"),
    file.path(MS_FIG_SUPP_TIFF_DIR, "FigS2_dynamic_storage_corr.tiff"),
    file.path(MS_FIG_SUPP_DIR, "FigS3_mobile_storage_corr.png"),
    file.path(MS_FIG_SUPP_PDF_DIR, "FigS3_mobile_storage_corr.pdf"),
    file.path(MS_FIG_SUPP_TIFF_DIR, "FigS3_mobile_storage_corr.tiff"),
    file.path(MS_FIG_SUPP_DIR, "FigS4_dynamic_mobile_corr.png"),
    file.path(MS_FIG_SUPP_PDF_DIR, "FigS4_dynamic_mobile_corr.pdf"),
    file.path(MS_FIG_SUPP_TIFF_DIR, "FigS4_dynamic_mobile_corr.tiff"),
    file.path(MS_TABLES_SUPP_DIR, "TableS7_MTT_sensitivity.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS8_catchment_char_storage_mlr_model_stats.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS9_catchment_alt_models_unique_deltaAICc_le2_BF.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS10_storage_eco_response_mlr_model_stats.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS11_eco_alt_models_unique_deltaAICc_le2_BF.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS12_mlr_model_diagnostics.csv")
  )

  if (isTRUE(WRITE_AUX_OUTPUTS)) {
    required_outputs <- c(
      required_outputs,
      file.path(master_dir, "master_site_metric_summary_stats.csv"),
      file.path(OUT_MET_SUPPORT_DIR, "site_metric_availability.csv")
    )
  }

  if (isTRUE(WRITE_TABLE_OUTPUTS)) {
    required_outputs <- c(
      required_outputs,
      file.path(OUT_TABLES_MLR_DIR, "mlr_main_results_table.csv")
    )
  }

  missing_outputs <- required_outputs[!file.exists(required_outputs)]
  if (length(missing_outputs) > 0) {
    stop(
      paste0(
        "Missing output file(s):\n- ",
        paste(missing_outputs, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }

  table_outputs <- c(
    file.path(MS_TABLES_SUPP_DIR, "TableS7_MTT_sensitivity.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS8_catchment_char_storage_mlr_model_stats.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS9_catchment_alt_models_unique_deltaAICc_le2_BF.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS10_storage_eco_response_mlr_model_stats.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS11_eco_alt_models_unique_deltaAICc_le2_BF.csv"),
    file.path(MS_TABLES_SUPP_DIR, "TableS12_mlr_model_diagnostics.csv")
  )

  for (table_file in table_outputs) {
    table_df <- readr::read_csv(
      table_file,
      show_col_types = FALSE,
      name_repair = "minimal"
    )
    bad_names <- names(table_df) == "" | is.na(names(table_df))
    if (any(bad_names)) {
      stop("Table has missing column heading(s): ", table_file, call. = FALSE)
    }
    duplicated_names <- names(table_df)[duplicated(names(table_df))]
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

  annual <- readr::read_csv(
    file.path(master_dir, MASTER_ANNUAL_FILE),
    show_col_types = FALSE
  )
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
  if (length(unknown_sites) > 0) {
    stop(
      "Unknown site code(s) in master_annual: ",
      paste(unknown_sites, collapse = ", "),
      call. = FALSE
    )
  }

  letters_path <- file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv")
  letters_df <- readr::read_csv(letters_path, show_col_types = FALSE)
  if (!all(c("metric", "site") %in% names(letters_df))) {
    stop("tukey_group_letters is missing columns: metric/site", call. = FALSE)
  }

  metric_groups <- split(as.character(letters_df$site), letters_df$metric)
  for (metric_name in names(metric_groups)) {
    site_values <- metric_groups[[metric_name]]
    allowed <- SITE_ORDER_HYDROMETRIC[SITE_ORDER_HYDROMETRIC %in% site_values]
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

# functions used to build daily meteorological and streamflow inputs
# this section does not read study data on its own


# date parsing

#' Parse dates robustly (handles Excel serials and various string formats)
#' @param d Vector of dates (character, numeric, or Date)
#' @return Vector of Date objects
parse_my_date <- function(d) {
  dd <- as.character(d)
  is_serial <- grepl("^[0-9]+$", dd)
  out <- rep(as.Date(NA), length(dd))

  # Excel serials → date

  if (any(is_serial)) {
    out[is_serial] <- as_date(
      as.numeric(dd[is_serial]),
      origin = "1899-12-30"
    )
  }

  # string dates → lubridate
  if (any(!is_serial)) {
    out[!is_serial] <- parse_date_time(
      dd[!is_serial],
      orders = c("Ymd", "Y-m-d", "mdy", "m/d/Y", "dmy", "d/m/Y", "dbY"),
      quiet = TRUE
    ) %>% as_date()
  }

  return(out)
}

# data loading helpers

#' Read interpolated meteorological data and convert to long format
#' @param fname Filename of CSV with *_inter columns
#' @param var Variable name for the value column
#' @param met_dir Directory containing the file
#' @return Long format tibble with DATE, SITECODE, and variable column
make_inter_long <- function(fname, var, met_dir, date_start = NULL, date_end = NULL) {
  raw <- read_csv(file.path(met_dir, fname), show_col_types = FALSE) %>%
    mutate(DATE = parse_my_date(DATE))

  if (!is.null(date_start)) {
    raw <- raw %>% filter(DATE >= as.Date(date_start))
  }
  if (!is.null(date_end)) {
    raw <- raw %>% filter(DATE <= as.Date(date_end))
  }

  long <- raw %>%
    pivot_longer(
      cols           = ends_with("_inter"),
      names_to       = "SITECODE",
      names_pattern  = "(.*)_inter$",
      values_to      = var,
      values_drop_na = TRUE
    )

  all_dates     <- seq.Date(min(raw$DATE, na.rm = TRUE),
                            max(raw$DATE, na.rm = TRUE),
                            by = "day")
  station_names <- unique(long$SITECODE)

  expand_grid(DATE = all_dates, SITECODE = station_names) %>%
    left_join(long, by = c("DATE", "SITECODE")) %>%
    arrange(SITECODE, DATE)
}

#' Read Mack Creek precipitation data
#' @param fname Filename
#' @param met_dir Directory containing the file
#' @return Tibble with DATE, SITECODE, P_mm_d
read_mack_precip <- function(fname, met_dir, recode_map = c("GSWSMC" = "GSMACK")) {
  read_csv(file.path(met_dir, fname), show_col_types = FALSE) %>%
    mutate(
      DATE     = parse_my_date(DATE),
      SITECODE = recode(SITECODE, !!!as.list(recode_map))
    ) %>%
    filter(SITECODE == "GSMACK") %>%
    select(DATE, SITECODE, PRECIP_TOT_DAY) %>%
    rename(P_mm_d = PRECIP_TOT_DAY)
}

# meteorological calculations

#' Calculate vapor pressure deficit (VPD) from temperature and relative humidity
#' @param temp_celsius Temperature in Celsius
#' @param rh_percent Relative humidity as percent (0-100)
#' @return VPD in kPa
calculate_vpd <- function(temp_celsius, rh_percent) {
  es <- 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))
  ea <- es * (rh_percent / 100)
  vpd <- es - ea
  return(vpd)
}

# pairwise relationship plot for interpolation checks
create_relationship_plot <- function(data, site1, site2, variable, r_squared, complete_count) {
  site1_col <- site1
  site2_col <- site2

  ggplot(data, aes(x = .data[[site1_col]], y = .data[[site2_col]])) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    labs(
      title = paste("Relationship for", variable, "between stations"),
      subtitle = paste(site1, "and", site2, "(n =", complete_count, ")"),
      x = site1,
      y = site2
    ) +
    annotate(
      "text",
      x = min(data[[site1_col]], na.rm = TRUE) + 0.8 * (max(data[[site1_col]], na.rm = TRUE) - min(data[[site1_col]], na.rm = TRUE)),
      y = min(data[[site2_col]], na.rm = TRUE) + 0.1 * (max(data[[site2_col]], na.rm = TRUE) - min(data[[site2_col]], na.rm = TRUE)),
      label = sprintf("R2 = %.3f", r_squared),
      hjust = 1,
      fontface = "bold"
    )
}

# model summary panel for triplet interpolation checks
create_multiple_regression_plot <- function(target_site, predictor_sites, variable, model_summary, complete_count) {
  title <- paste("Multiple Regression Model for", variable)
  subtitle <- paste(target_site, "predicted from", paste(predictor_sites, collapse = ", "))

  model_info <- data.frame(
    Metric = c("R2", "Adjusted R2", "F-statistic", "p-value", "Sample Size"),
    Value = c(
      round(model_summary$r.squared, 4),
      round(model_summary$adj.r.squared, 4),
      round(model_summary$fstatistic[1], 2),
      format.pval(
        pf(
          model_summary$fstatistic[1],
          model_summary$fstatistic[2],
          model_summary$fstatistic[3],
          lower.tail = FALSE
        ),
        digits = 3
      ),
      complete_count
    )
  )

  coef_info <- data.frame(
    Term = rownames(model_summary$coefficients),
    Estimate = model_summary$coefficients[, "Estimate"],
    `Pr(>|t|)` = model_summary$coefficients[, "Pr(>|t|)"]
  )

  ggplot() +
    theme_minimal(base_size = 12) +
    annotate("text", x = 0.5, y = 0.9, label = title, fontface = "bold", size = 5, hjust = 0.5) +
    annotate("text", x = 0.5, y = 0.85, label = subtitle, fontface = "italic", size = 4, hjust = 0.5) +
    annotate("text", x = 0.5, y = 0.75, label = "Model Metrics:", fontface = "bold", size = 4, hjust = 0.5) +
    annotate("text", x = 0.3, y = 0.7, label = paste(model_info$Metric, collapse = "\n"), hjust = 1, size = 3) +
    annotate("text", x = 0.7, y = 0.7, label = paste(model_info$Value, collapse = "\n"), hjust = 0, size = 3) +
    annotate("text", x = 0.5, y = 0.5, label = "Coefficients:", fontface = "bold", size = 4, hjust = 0.5) +
    annotate("text", x = 0.2, y = 0.45, label = paste(c("Term", coef_info$Term), collapse = "\n"), hjust = 0, fontface = "bold", size = 3) +
    annotate("text", x = 0.45, y = 0.45, label = paste(c("Estimate", round(coef_info$Estimate, 4)), collapse = "\n"), hjust = 0, size = 3) +
    annotate("text", x = 0.7, y = 0.45, label = paste(c("p-value", format.pval(coef_info$`Pr(>|t|)`, digits = 3)), collapse = "\n"), hjust = 0, size = 3) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(title = title, subtitle = subtitle)
}

# station interpolation

#' Extract station pairs and triplets that need interpolation from site mapping
#' @param site_mapping Named list of catchment -> station mappings
#' @return List with $pairs and $triplets
extract_station_groups <- function(site_mapping) {
  pairs <- list()
  triplets <- list()

  for (site in names(site_mapping)) {
    site_info <- site_mapping[[site]]

    for (var in c("temp", "precip", "rh", "netrad")) {
      stations <- site_info[[var]]

      if (length(stations) == 2) {
        pair_name <- paste(stations[1], stations[2], sep = "_")
        if (!pair_name %in% names(pairs)) {
          pairs[[pair_name]] <- list(site1 = stations[1], site2 = stations[2])
        }
      } else if (length(stations) == 3) {
        triplet_name <- paste(stations, collapse = "_")
        if (!triplet_name %in% names(triplets)) {
          triplets[[triplet_name]] <- list(
            site1 = stations[1],
            site2 = stations[2],
            site3 = stations[3]
          )
        }
      }
    }
  }

  return(list(pairs = pairs, triplets = triplets))
}

#' Interpolate missing values between two stations using linear regression
#' @param data Data frame with DATE, SITECODE, and variable columns
#' @param site1 First station code
#' @param site2 Second station code
#' @param variable Variable name to interpolate
#' @return Data frame with interpolated values added
interpolate_pair <- function(data, site1, site2, variable, plot_dir = NULL) {
  # check for duplicates
  check_dupes <- data %>%
    filter(SITECODE %in% c(site1, site2)) %>%
    group_by(DATE, SITECODE) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)

  if(nrow(check_dupes) > 0) {
    data <- data %>%
      group_by(DATE, SITECODE) %>%
      summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  }

  pair_data <- data %>%
    filter(SITECODE %in% c(site1, site2)) %>%
    select(DATE, SITECODE, !!sym(variable)) %>%
    pivot_wider(names_from = SITECODE, values_from = !!sym(variable), values_fn = mean)

  site1_col <- site1
  site2_col <- site2

  complete_rows <- complete.cases(pair_data[, c(site1_col, site2_col)])
  complete_count <- sum(complete_rows)

  if(complete_count < 5) {
    return(data)
  }

  # fit linear model
  model <- lm(formula = paste0("`", site2_col, "` ~ `", site1_col, "`"),
              data = pair_data[complete_rows, ])

  r_squared <- summary(model)$r.squared
  intercept <- coef(model)[1]
  slope <- coef(model)[2]

  if (!is.null(plot_dir) && nzchar(plot_dir)) {
    plot_data <- pair_data[complete_rows, ]
    if (nrow(plot_data) > 1) {
      p <- create_relationship_plot(plot_data, site1, site2, variable, r_squared, complete_count)
      ggsave(
        filename = file.path(plot_dir, paste0("relationship_", variable, "_", site1, "_", site2, ".png")),
        plot = p,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  }

  # interpolate site2 from site1
  s1_has_data_s2_missing <- !is.na(pair_data[[site1_col]]) & is.na(pair_data[[site2_col]])
  if(any(s1_has_data_s2_missing)) {
    predictions <- intercept + slope * pair_data[[site1_col]][s1_has_data_s2_missing]

    temp_df <- data.frame(
      DATE = pair_data$DATE[s1_has_data_s2_missing],
      SITECODE = site2,
      value = predictions
    )
    names(temp_df)[3] <- variable

    data <- bind_rows(data, temp_df) %>%
      distinct(DATE, SITECODE, .keep_all = TRUE)
  }

  # interpolate site1 from site2 (inverse relationship)
  s2_has_data_s1_missing <- is.na(pair_data[[site1_col]]) & !is.na(pair_data[[site2_col]])
  if(any(s2_has_data_s1_missing)) {
    inv_intercept <- -intercept/slope
    inv_slope <- 1/slope
    predictions <- inv_intercept + inv_slope * pair_data[[site2_col]][s2_has_data_s1_missing]

    temp_df <- data.frame(
      DATE = pair_data$DATE[s2_has_data_s1_missing],
      SITECODE = site1,
      value = predictions
    )
    names(temp_df)[3] <- variable

    data <- bind_rows(data, temp_df) %>%
      distinct(DATE, SITECODE, .keep_all = TRUE)
  }

  return(data)
}

#' Interpolate missing values using multiple regression from triplet of stations
#' @param data Data frame with DATE, SITECODE, and variable columns
#' @param site1 First station code
#' @param site2 Second station code
#' @param site3 Third station code
#' @param variable Variable name to interpolate
#' @return Data frame with interpolated values added
interpolate_triplet <- function(data, site1, site2, site3, variable, plot_dir = NULL) {
  # check for duplicates
  check_dupes <- data %>%
    filter(SITECODE %in% c(site1, site2, site3)) %>%
    group_by(DATE, SITECODE) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)

  if(nrow(check_dupes) > 0) {
    data <- data %>%
      group_by(DATE, SITECODE) %>%
      summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  }

  triplet_data <- data %>%
    filter(SITECODE %in% c(site1, site2, site3)) %>%
    select(DATE, SITECODE, !!sym(variable)) %>%
    pivot_wider(names_from = SITECODE, values_from = !!sym(variable), values_fn = mean)

  site1_col <- site1
  site2_col <- site2
  site3_col <- site3

  model1_rows <- complete.cases(triplet_data[, c(site1_col, site2_col, site3_col)])
  model1_count <- sum(model1_rows)

  if(model1_count < 10) {
    return(data)
  }

  # create three models (one for each site)
  model1 <- lm(formula = paste0("`", site1_col, "` ~ `", site2_col, "` + `", site3_col, "`"),
               data = triplet_data[model1_rows, ])
  model2 <- lm(formula = paste0("`", site2_col, "` ~ `", site1_col, "` + `", site3_col, "`"),
               data = triplet_data[model1_rows, ])
  model3 <- lm(formula = paste0("`", site3_col, "` ~ `", site1_col, "` + `", site2_col, "`"),
               data = triplet_data[model1_rows, ])

  if (!is.null(plot_dir) && nzchar(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    coef_value <- function(coef_tbl, pred_name, col_name) {
      idx <- which(rownames(coef_tbl) %in% c(pred_name, paste0("`", pred_name, "`")))
      if (length(idx) == 0) return(NA_real_)
      suppressWarnings(as.numeric(coef_tbl[idx[1], col_name]))
    }

    summarize_triplet_model <- function(target_site, pred_site1, pred_site2, model_obj, n_obs, variable_name) {
      sm <- summary(model_obj)
      coef_tbl <- sm$coefficients
      tibble(
        variable = variable_name,
        target_site = target_site,
        predictor_1 = pred_site1,
        predictor_2 = pred_site2,
        n_complete = n_obs,
        r2 = as.numeric(sm$r.squared),
        adj_r2 = as.numeric(sm$adj.r.squared),
        intercept = suppressWarnings(as.numeric(coef_tbl["(Intercept)", "Estimate"])),
        beta_1 = coef_value(coef_tbl, pred_site1, "Estimate"),
        p_beta_1 = coef_value(coef_tbl, pred_site1, "Pr(>|t|)"),
        beta_2 = coef_value(coef_tbl, pred_site2, "Estimate"),
        p_beta_2 = coef_value(coef_tbl, pred_site2, "Pr(>|t|)"),
        rmse = sqrt(mean(residuals(model_obj)^2, na.rm = TRUE))
      )
    }

    diag_rows <- bind_rows(
      summarize_triplet_model(site1, site2, site3, model1, model1_count, variable),
      summarize_triplet_model(site2, site1, site3, model2, model1_count, variable),
      summarize_triplet_model(site3, site1, site2, model3, model1_count, variable)
    )

    diag_file <- file.path(plot_dir, "multireg_triplet_diagnostics.csv")
    if (file.exists(diag_file)) {
      existing_diag <- read_csv(diag_file, show_col_types = FALSE)
      diag_rows <- bind_rows(existing_diag, diag_rows)
    }
    diag_rows <- diag_rows %>%
      distinct(variable, target_site, predictor_1, predictor_2, .keep_all = TRUE) %>%
      arrange(variable, target_site, predictor_1, predictor_2)

    write_csv(diag_rows, diag_file)
  }

  # case 1: only site1 is missing
  case1 <- is.na(triplet_data[[site1_col]]) & !is.na(triplet_data[[site2_col]]) & !is.na(triplet_data[[site3_col]])
  if(any(case1)) {
    temp_df <- triplet_data[case1, c("DATE", site2_col, site3_col)]
    names(temp_df) <- c("DATE", site2_col, site3_col)
    predictions <- predict(model1, newdata = temp_df)

    interp_df <- data.frame(
      DATE = triplet_data$DATE[case1],
      SITECODE = site1,
      value = predictions
    )
    names(interp_df)[3] <- variable

    data <- bind_rows(data, interp_df) %>%
      distinct(DATE, SITECODE, .keep_all = TRUE)
  }

  # case 2: only site2 is missing
  case2 <- !is.na(triplet_data[[site1_col]]) & is.na(triplet_data[[site2_col]]) & !is.na(triplet_data[[site3_col]])
  if(any(case2)) {
    temp_df <- triplet_data[case2, c("DATE", site1_col, site3_col)]
    names(temp_df) <- c("DATE", site1_col, site3_col)
    predictions <- predict(model2, newdata = temp_df)

    interp_df <- data.frame(
      DATE = triplet_data$DATE[case2],
      SITECODE = site2,
      value = predictions
    )
    names(interp_df)[3] <- variable

    data <- bind_rows(data, interp_df) %>%
      distinct(DATE, SITECODE, .keep_all = TRUE)
  }

  # case 3: only site3 is missing
  case3 <- !is.na(triplet_data[[site1_col]]) & !is.na(triplet_data[[site2_col]]) & is.na(triplet_data[[site3_col]])
  if(any(case3)) {
    temp_df <- triplet_data[case3, c("DATE", site1_col, site2_col)]
    names(temp_df) <- c("DATE", site1_col, site2_col)
    predictions <- predict(model3, newdata = temp_df)

    interp_df <- data.frame(
      DATE = triplet_data$DATE[case3],
      SITECODE = site3,
      value = predictions
    )
    names(interp_df)[3] <- variable

    data <- bind_rows(data, interp_df) %>%
      distinct(DATE, SITECODE, .keep_all = TRUE)
  }

  return(data)
}

#' Constrain interpolated values to physical limits
#' @param data Data frame with RH_d_pct and/or P_mm_d columns
#' @return Data frame with constrained values
constrain_interpolated_values <- function(data) {
  # cap rh at 100%
  if ("RH_d_pct" %in% names(data)) {
    over_100_count <- sum(data$RH_d_pct > 100, na.rm = TRUE)
    if (over_100_count > 0) {
      data <- data %>% mutate(RH_d_pct = ifelse(RH_d_pct > 100, 100, RH_d_pct))
    }
  }

  # cap precipitation at 0 (no negatives)
  if ("P_mm_d" %in% names(data)) {
    neg_precip_count <- sum(data$P_mm_d < 0, na.rm = TRUE)
    if (neg_precip_count > 0) {
      data <- data %>% mutate(P_mm_d = ifelse(P_mm_d < 0, 0, P_mm_d))
    }
  }

  return(data)
}

# main processing functions

#' Process all station pairs and triplets, apply interpolation
#' @param data Combined meteorological data
#' @param station_groups List with $pairs and $triplets from extract_station_groups()
#' @param variables Character vector of variable names to process
#' @return List with interpolated data and tracking info
process_station_groups <- function(data, station_groups, variables, plot_dir = NULL) {
  interpolated_data <- data
  interpolated_pairs <- list()
  interpolated_triplets <- list()

  # process pairs
  for (pair_name in names(station_groups$pairs)) {
    pair <- station_groups$pairs[[pair_name]]
    site1 <- pair$site1
    site2 <- pair$site2

    pair_can_be_interpolated <- FALSE

    for (var in variables) {
      pair_data <- data %>%
        filter(SITECODE %in% c(site1, site2)) %>%
        select(DATE, SITECODE, !!sym(var)) %>%
        pivot_wider(names_from = SITECODE, values_from = !!sym(var), values_fn = mean)

      complete_count <- sum(complete.cases(pair_data[, c(site1, site2)]))

      if (complete_count >= 5) {
        pair_can_be_interpolated <- TRUE
        break
      }
    }

    if (pair_can_be_interpolated) {
      for (var in variables) {
        interpolated_data <- interpolate_pair(interpolated_data, site1, site2, var, plot_dir = plot_dir)
      }
      interpolated_pairs[[pair_name]] <- pair
    }
  }

  # process triplets
  for (triplet_name in names(station_groups$triplets)) {
    triplet <- station_groups$triplets[[triplet_name]]
    site1 <- triplet$site1
    site2 <- triplet$site2
    site3 <- triplet$site3

    triplet_can_be_interpolated <- FALSE

    for (var in variables) {
      triplet_data <- data %>%
        filter(SITECODE %in% c(site1, site2, site3)) %>%
        select(DATE, SITECODE, !!sym(var)) %>%
        pivot_wider(names_from = SITECODE, values_from = !!sym(var), values_fn = mean)

      complete_count <- sum(complete.cases(triplet_data[, c(site1, site2, site3)]))

      if (complete_count >= 10) {
        triplet_can_be_interpolated <- TRUE
        break
      }
    }

    if (triplet_can_be_interpolated) {
      for (var in variables) {
        interpolated_data <- interpolate_triplet(interpolated_data, site1, site2, site3, var, plot_dir = plot_dir)
      }
      interpolated_triplets[[triplet_name]] <- triplet
    } else {
      # fall back to pairwise interpolation
      fallback_pairs <- list(
        list(site1 = site1, site2 = site2),
        list(site1 = site1, site2 = site3),
        list(site1 = site2, site2 = site3)
      )

      for (pair in fallback_pairs) {
        for (var in variables) {
          interpolated_data <- interpolate_pair(interpolated_data, pair$site1, pair$site2, var, plot_dir = plot_dir)
        }
      }
    }
  }

  # apply constraints and calculate vpd
  interpolated_data <- constrain_interpolated_values(interpolated_data)

  interpolated_data <- interpolated_data %>%
    rowwise() %>%
    mutate(VPD_kPa = if_else(!is.na(T_C) & !is.na(RH_d_pct),
                             calculate_vpd(T_C, RH_d_pct),
                             as.numeric(NA))) %>%
    ungroup()

  return(list(
    data = interpolated_data,
    interpolated_pairs = interpolated_pairs,
    interpolated_triplets = interpolated_triplets
  ))
}

#' Create catchment level datasets from interpolated station data
#' @param interpolated_data Interpolated station data
#' @param site_mapping Named list of catchment -> station mappings
#' @param variables Character vector of variable names
#' @return Named list of catchment datasets
create_catchment_datasets <- function(interpolated_data, site_mapping, variables) {
  site_datasets <- list()

  for (site_name in names(site_mapping)) {
    site_info <- site_mapping[[site_name]]
    site_data <- data.frame()

    for (var_idx in seq_along(variables)) {
      var <- variables[var_idx]

      # handle vpd separately (uses temperature stations)
      if (var == "VPD_kPa") {
        var_name <- "temp"
      } else if (var_idx > length(names(site_info))) {
        next
      } else {
        var_name <- names(site_info)[var_idx]
      }

      stations <- site_info[[var_name]]

      if (length(stations) == 1) {
        # single station
        single_data <- interpolated_data %>%
          filter(SITECODE == stations[1]) %>%
          select(DATE, !!sym(var))

        if (nrow(site_data) == 0) {
          site_data <- single_data %>% mutate(SITECODE = site_name)
        } else {
          site_data <- site_data %>% left_join(single_data, by = "DATE")
        }
      } else if (length(stations) >= 2) {
        # average multiple stations
        multi_data <- interpolated_data %>%
          filter(SITECODE %in% stations) %>%
          select(DATE, SITECODE, !!sym(var)) %>%
          pivot_wider(names_from = SITECODE, values_from = !!sym(var)) %>%
          mutate(avg_value = rowMeans(select(., all_of(stations)), na.rm = TRUE)) %>%
          select(DATE, avg_value)

        names(multi_data)[2] <- var

        if (nrow(site_data) == 0) {
          site_data <- multi_data %>% mutate(SITECODE = site_name)
        } else {
          site_data <- site_data %>% left_join(multi_data, by = "DATE")
        }
      }
    }

    site_datasets[[site_name]] <- site_data
  }

  return(site_datasets)
}

#' Add discharge data to catchment datasets
#' @param catchment_datasets Named list of catchment datasets
#' @param discharge Discharge data frame with DATE, SITECODE, MEAN_Q
#' @param da_df Drainage area data frame with SITECODE, DA_M2
#' @return Updated catchment datasets with Q_mm_d column
add_discharge_to_catchments <- function(catchment_datasets, discharge, da_df) {
  discharge_processed <- discharge %>%
    left_join(da_df, by = "SITECODE") %>%
    filter(!is.na(DA_M2)) %>%
    mutate(Q_mm_d = MEAN_Q * 0.0283168 * 86400 / DA_M2 * 1000) %>%
    select(DATE, SITECODE, Q_mm_d)

  for (site_name in names(catchment_datasets)) {
    site_dis <- discharge_processed %>% filter(SITECODE == site_name)
    if (nrow(site_dis) > 0) {
      catchment_datasets[[site_name]] <- catchment_datasets[[site_name]] %>%
        left_join(site_dis, by = "DATE")
    }
  }

  return(catchment_datasets)
}

# triplet RH diagnostics for WS7MET, VANMET, and H15MET
plot_triplet_station_comparisons <- function(interpolated_data, plot_dir) {
  triplet_stations <- c("WS7MET", "VANMET", "H15MET")

  triplet_rh_data <- interpolated_data %>%
    filter(SITECODE %in% triplet_stations) %>%
    select(DATE, SITECODE, RH_d_pct) %>%
    pivot_wider(names_from = SITECODE, values_from = RH_d_pct)

  complete_triplet_data <- triplet_rh_data %>%
    filter(!is.na(WS7MET) & !is.na(VANMET) & !is.na(H15MET))

  if (nrow(complete_triplet_data) == 0) {
    return(NULL)
  }

  p1 <- ggplot(complete_triplet_data, aes(x = WS7MET, y = VANMET)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_classic(base_size = 12) +
    labs(title = "Comparison of Relative Humidity", subtitle = "WS7MET vs VANMET", x = "WS7MET Relative Humidity (%)", y = "VANMET Relative Humidity (%)")

  model_ws7_van <- lm(VANMET ~ WS7MET, data = complete_triplet_data)
  p1 <- p1 + annotate("text", x = min(complete_triplet_data$WS7MET, na.rm = TRUE) + 0.8 * (max(complete_triplet_data$WS7MET, na.rm = TRUE) - min(complete_triplet_data$WS7MET, na.rm = TRUE)), y = min(complete_triplet_data$VANMET, na.rm = TRUE) + 0.1 * (max(complete_triplet_data$VANMET, na.rm = TRUE) - min(complete_triplet_data$VANMET, na.rm = TRUE)), label = sprintf("R2 = %.3f", summary(model_ws7_van)$r.squared), hjust = 1, fontface = "bold")

  p2 <- ggplot(complete_triplet_data, aes(x = WS7MET, y = H15MET)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_classic(base_size = 12) +
    labs(title = "Comparison of Relative Humidity", subtitle = "WS7MET vs H15MET", x = "WS7MET Relative Humidity (%)", y = "H15MET Relative Humidity (%)")

  model_ws7_h15 <- lm(H15MET ~ WS7MET, data = complete_triplet_data)
  p2 <- p2 + annotate("text", x = min(complete_triplet_data$WS7MET, na.rm = TRUE) + 0.8 * (max(complete_triplet_data$WS7MET, na.rm = TRUE) - min(complete_triplet_data$WS7MET, na.rm = TRUE)), y = min(complete_triplet_data$H15MET, na.rm = TRUE) + 0.1 * (max(complete_triplet_data$H15MET, na.rm = TRUE) - min(complete_triplet_data$H15MET, na.rm = TRUE)), label = sprintf("R2 = %.3f", summary(model_ws7_h15)$r.squared), hjust = 1, fontface = "bold")

  p3 <- ggplot(complete_triplet_data, aes(x = VANMET, y = H15MET)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_classic(base_size = 12) +
    labs(title = "Comparison of Relative Humidity", subtitle = "VANMET vs H15MET", x = "VANMET Relative Humidity (%)", y = "H15MET Relative Humidity (%)")

  model_van_h15 <- lm(H15MET ~ VANMET, data = complete_triplet_data)
  p3 <- p3 + annotate("text", x = min(complete_triplet_data$VANMET, na.rm = TRUE) + 0.8 * (max(complete_triplet_data$VANMET, na.rm = TRUE) - min(complete_triplet_data$VANMET, na.rm = TRUE)), y = min(complete_triplet_data$H15MET, na.rm = TRUE) + 0.1 * (max(complete_triplet_data$H15MET, na.rm = TRUE) - min(complete_triplet_data$H15MET, na.rm = TRUE)), label = sprintf("R2 = %.3f", summary(model_van_h15)$r.squared), hjust = 1, fontface = "bold")

  triplet_long <- triplet_rh_data %>%
    pivot_longer(cols = c("WS7MET", "VANMET", "H15MET"), names_to = "Station", values_to = "RH_pct")

  p4 <- ggplot(triplet_long, aes(x = DATE, y = RH_pct, color = Station)) +
    geom_line(linewidth = 0.5) +
    theme_classic(base_size = 12) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Time Series Comparison of Relative Humidity", subtitle = "Triplet Stations: WS7MET, VANMET, and H15MET", x = "Date", y = "Relative Humidity (%)")

  ggsave(file.path(plot_dir, "RH_comparison_WS7MET_vs_VANMET.png"), plot = p1, width = 8, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "RH_comparison_WS7MET_vs_H15MET.png"), plot = p2, width = 8, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "RH_comparison_VANMET_vs_H15MET.png"), plot = p3, width = 8, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "RH_timeseries_triplet_stations.png"), plot = p4, width = 12, height = 6, dpi = 300)

  list(
    model_ws7_van = model_ws7_van,
    model_ws7_h15 = model_ws7_h15,
    model_van_h15 = model_van_h15
  )
}

# model selection, AICc, VIF filtering, and diagnostics
# used by the catchment-control, ecological-response, and MTT sensitivity scripts

calc_aicc <- function(model_obj, n_obs) {
  k_params <- length(coef(model_obj)) + 1
  aic_val <- AIC(model_obj)
  if ((n_obs - k_params - 1) <= 0) {
    return(NA_real_)
  }
  aic_val + (2 * k_params * (k_params + 1)) / (n_obs - k_params - 1)
}

calc_loocv_stats <- function(model_formula, model_df, min_n = 6) {
  n <- nrow(model_df)
  if (n < min_n) {
    return(list(rmse = NA_real_, mae = NA_real_, r2 = NA_real_))
  }

  response <- all.vars(model_formula)[1]
  obs <- model_df[[response]]
  pred <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    fit_i <- tryCatch(lm(model_formula, data = model_df[-i, , drop = FALSE]), error = function(e) NULL)
    if (!is.null(fit_i)) {
      pred[i] <- tryCatch(
        as.numeric(predict(fit_i, newdata = model_df[i, , drop = FALSE])),
        error = function(e) NA_real_
      )
    }
  }

  valid <- is.finite(obs) & is.finite(pred)
  if (sum(valid) < 3) {
    return(list(rmse = NA_real_, mae = NA_real_, r2 = NA_real_))
  }

  errs <- obs[valid] - pred[valid]
  sst <- sum((obs[valid] - mean(obs[valid], na.rm = TRUE))^2, na.rm = TRUE)
  sse <- sum(errs^2, na.rm = TRUE)

  list(
    rmse = sqrt(mean(errs^2, na.rm = TRUE)),
    mae = mean(abs(errs), na.rm = TRUE),
    r2 = ifelse(sst > 0, 1 - sse / sst, NA_real_)
  )
}

compute_residual_diagnostics <- function(model_obj) {
  pull_scalar <- function(obj, key) {
    val <- tryCatch(obj[[key]], error = function(e) NULL)
    num <- suppressWarnings(as.numeric(val))
    if (length(num) >= 1 && is.finite(num[1])) {
      return(num[1])
    }
    NA_real_
  }

  resid_vals <- residuals(model_obj)
  resid_vals <- resid_vals[is.finite(resid_vals)]

  shapiro_w <- NA_real_
  shapiro_p <- NA_real_
  if (length(resid_vals) >= 3 && length(resid_vals) <= 5000) {
    sh <- tryCatch(shapiro.test(resid_vals), error = function(e) NULL)
    if (!is.null(sh)) {
      shapiro_w <- suppressWarnings(as.numeric(unname(sh$statistic)))
      shapiro_p <- suppressWarnings(as.numeric(sh$p.value))
    }
  }

  ncv <- tryCatch(car::ncvTest(model_obj), error = function(e) NULL)
  ncv_chisq <- if (!is.null(ncv)) pull_scalar(ncv, "Chisquare") else NA_real_
  if (!is.finite(ncv_chisq) && !is.null(ncv)) {
    ncv_chisq <- pull_scalar(ncv, "ChiSquare")
  }
  ncv_p <- if (!is.null(ncv)) pull_scalar(ncv, "p") else NA_real_

  tibble::tibble(
    n_residuals = length(resid_vals),
    shapiro_W = shapiro_w,
    shapiro_p = shapiro_p,
    ncv_chisq = ncv_chisq,
    ncv_p = ncv_p,
    normality_pass_p05 = ifelse(is.finite(shapiro_p), shapiro_p > 0.05, NA),
    homoscedasticity_pass_p05 = ifelse(is.finite(ncv_p), ncv_p > 0.05, NA)
  )
}

round_export_cols <- function(df, cols, digits = 3) {
  keep <- intersect(cols, names(df))
  if (length(keep) == 0) {
    return(df)
  }
  dplyr::mutate(df, dplyr::across(dplyr::all_of(keep), ~ signif(.x, digits)))
}

format_export_outcome <- function(x, strip_mean = TRUE, drop_underscores = TRUE) {
  out <- as.character(x)
  if (isTRUE(strip_mean)) {
    out <- gsub("_mean$", "", out)
  }
  if (isTRUE(drop_underscores)) {
    out <- gsub("_", "", out, fixed = TRUE)
  }
  out
}

apply_vif_filter <- function(model_obj, outcome, model_df, threshold, protected = character()) {
  fit <- model_obj
  protected <- as.character(protected)

  repeat {
    retained <- setdiff(names(coef(fit)), "(Intercept)")
    drop_pool <- setdiff(retained, protected)
    if (length(drop_pool) == 0 || length(retained) <= 1) {
      break
    }

    vif_vals <- tryCatch(car::vif(fit), error = function(e) NULL)
    if (is.null(vif_vals) || max(vif_vals, na.rm = TRUE) <= threshold) {
      break
    }

    ordered <- names(sort(vif_vals, decreasing = TRUE))
    drop_var <- ordered[ordered %in% drop_pool][1]
    if (is.na(drop_var) || !nzchar(drop_var)) {
      break
    }

    keep <- setdiff(retained, drop_var)
    if (length(keep) == 0) {
      return(NULL)
    }

    fit <- tryCatch(
      lm(as.formula(paste(outcome, "~", paste(keep, collapse = " + "))), data = model_df),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    fit <- tryCatch(MASS::stepAIC(fit, direction = "backward", trace = 0), error = function(e) fit)
  }

  fit
}

apply_correlated_predictor_rules_catchment <- function(model_obj, outcome, model_df) {
  refit_drop_var <- function(drop_var) {
    keep <- setdiff(setdiff(names(coef(model_obj_curr)), "(Intercept)"), drop_var)
    if (length(keep) == 0) {
      return(NULL)
    }
    fit <- tryCatch(
      lm(as.formula(paste(outcome, "~", paste(keep, collapse = " + "))), data = model_df),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    fit <- tryCatch(MASS::stepAIC(fit, direction = "backward", trace = 0), error = function(e) fit)
    keep_fit <- setdiff(names(coef(fit)), "(Intercept)")
    if (length(keep_fit) == 0) {
      return(NULL)
    }
    list(
      model = fit,
      aicc = calc_aicc(fit, nrow(model_df)),
      aic = tryCatch(AIC(fit), error = function(e) NA_real_)
    )
  }

  model_obj_curr <- model_obj
  repeat {
    retained <- setdiff(names(coef(model_obj_curr)), "(Intercept)")
    has_ash <- "Ash_Per" %in% retained
    has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained)
    has_both_landslide <- all(c("Landslide_Total", "Landslide_Young") %in% retained)

    if (!(has_ash && has_lava) && !has_both_landslide) {
      break
    }

    options <- list()
    if (has_ash && has_lava) {
      drop_candidates <- c("Ash_Per", intersect(c("Lava1_per", "Lava2_per"), retained))
      options <- lapply(drop_candidates, refit_drop_var)
    } else if (has_both_landslide) {
      options <- lapply(c("Landslide_Total", "Landslide_Young"), refit_drop_var)
    }

    valid <- which(vapply(options, function(x) !is.null(x), logical(1)))
    if (length(valid) == 0) {
      break
    }

    aicc_scores <- vapply(options[valid], function(x) x$aicc, numeric(1))
    if (any(is.finite(aicc_scores))) {
      best <- valid[which.min(aicc_scores)]
    } else {
      aic_scores <- vapply(options[valid], function(x) x$aic, numeric(1))
      if (any(is.finite(aic_scores))) {
        best <- valid[which.min(aic_scores)]
      } else {
        best <- valid[1]
      }
    }
    model_obj_curr <- options[[best]]$model
  }

  model_obj_curr
}
