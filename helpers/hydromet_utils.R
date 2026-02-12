# -----------------------------------------------------------------------------
# Hydrometric Data Processing Utilities
# -----------------------------------------------------------------------------
# Purpose: Helper functions for processing meteorological and hydrometric data
#
# Functions:
#   - parse_my_date(): Robust date parsing (handles Excel serials and strings)
#   - make_inter_long(): Convert wide interpolated data to long format
#   - read_mack_precip(): Read Mack Creek precipitation data
#   - calculate_vpd(): Vapor pressure deficit calculation
#   - interpolate_pair(): Interpolate missing values between two stations
#   - interpolate_triplet(): Interpolate using multiple regression (3 stations)
#   - constrain_interpolated_values(): Cap RH at 100%, precip at 0
#   - process_station_groups(): Main processing workflow for station pairs/triplets
#   - create_watershed_datasets(): Aggregate station data to watershed level
#
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# DATE PARSING
# -----------------------------------------------------------------------------

#' Parse dates robustly (handles Excel serials and various string formats)
#' @param d Vector of dates (character, numeric, or Date)
#' @return Vector of Date objects
parse_my_date <- function(d) {
  dd <- as.character(d)
  is_serial <- grepl("^[0-9]+$", dd)
  out <- rep(as.Date(NA), length(dd))

  # Excel serials → Date

  if (any(is_serial)) {
    out[is_serial] <- as_date(
      as.numeric(dd[is_serial]),
      origin = "1899-12-30"
    )
  }

  # String dates → lubridate
  if (any(!is_serial)) {
    out[!is_serial] <- parse_date_time(
      dd[!is_serial],
      orders = c("Ymd", "Y-m-d", "mdy", "m/d/Y", "dmy", "d/m/Y", "dbY"),
      quiet = TRUE
    ) %>% as_date()
  }

  return(out)
}

# -----------------------------------------------------------------------------
# DATA LOADING HELPERS
# -----------------------------------------------------------------------------

#' Read interpolated meteorological data and convert to long format
#' @param fname Filename of CSV with *_inter columns
#' @param var Variable name for the value column
#' @param met_dir Directory containing the file
#' @return Long-format tibble with DATE, SITECODE, and variable column
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

# -----------------------------------------------------------------------------
# METEOROLOGICAL CALCULATIONS
# -----------------------------------------------------------------------------

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

# Legacy-style pairwise relationship plot for interpolation diagnostics
create_relationship_plot <- function(data, site1, site2, variable, r_squared, complete_count) {
  site1_col <- site1
  site2_col <- site2

  ggplot(data, aes_string(x = paste0("`", site1_col, "`"), y = paste0("`", site2_col, "`"))) +
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

# Legacy-style model summary panel for triplet interpolation diagnostics
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

# -----------------------------------------------------------------------------
# STATION INTERPOLATION
# -----------------------------------------------------------------------------

#' Extract station pairs and triplets that need interpolation from site mapping
#' @param site_mapping Named list of watershed -> station mappings
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
  # Check for duplicates
  check_dupes <- data %>%
    filter(SITECODE %in% c(site1, site2)) %>%
    group_by(DATE, SITECODE) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)

  if(nrow(check_dupes) > 0) {
    data <- data %>%
      group_by(DATE, SITECODE) %>%
      summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
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

  # Fit linear model
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

  # Interpolate site2 from site1
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

  # Interpolate site1 from site2 (inverse relationship)
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
  # Check for duplicates
  check_dupes <- data %>%
    filter(SITECODE %in% c(site1, site2, site3)) %>%
    group_by(DATE, SITECODE) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)

  if(nrow(check_dupes) > 0) {
    data <- data %>%
      group_by(DATE, SITECODE) %>%
      summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
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

  # Create three models (one for each site)
  model1 <- lm(formula = paste0("`", site1_col, "` ~ `", site2_col, "` + `", site3_col, "`"),
               data = triplet_data[model1_rows, ])
  model2 <- lm(formula = paste0("`", site2_col, "` ~ `", site1_col, "` + `", site3_col, "`"),
               data = triplet_data[model1_rows, ])
  model3 <- lm(formula = paste0("`", site3_col, "` ~ `", site1_col, "` + `", site2_col, "`"),
               data = triplet_data[model1_rows, ])

  if (!is.null(plot_dir) && nzchar(plot_dir)) {
    model1_summary <- summary(model1)
    model2_summary <- summary(model2)
    model3_summary <- summary(model3)

    p1 <- create_multiple_regression_plot(site1, c(site2, site3), variable, model1_summary, model1_count)
    p2 <- create_multiple_regression_plot(site2, c(site1, site3), variable, model2_summary, model1_count)
    p3 <- create_multiple_regression_plot(site3, c(site1, site2), variable, model3_summary, model1_count)

    ggsave(file.path(plot_dir, paste0("multireg_", variable, "_", site1, "_from_", site2, "_", site3, ".png")), p1, width = 8, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, paste0("multireg_", variable, "_", site2, "_from_", site1, "_", site3, ".png")), p2, width = 8, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, paste0("multireg_", variable, "_", site3, "_from_", site1, "_", site2, ".png")), p3, width = 8, height = 6, dpi = 300)
  }

  # Case 1: Only site1 is missing
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

  # Case 2: Only site2 is missing
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

  # Case 3: Only site3 is missing
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
  # Cap RH at 100%
  if ("RH_d_pct" %in% names(data)) {
    over_100_count <- sum(data$RH_d_pct > 100, na.rm = TRUE)
    if (over_100_count > 0) {
      data <- data %>% mutate(RH_d_pct = ifelse(RH_d_pct > 100, 100, RH_d_pct))
    }
  }

  # Cap precipitation at 0 (no negatives)
  if ("P_mm_d" %in% names(data)) {
    neg_precip_count <- sum(data$P_mm_d < 0, na.rm = TRUE)
    if (neg_precip_count > 0) {
      data <- data %>% mutate(P_mm_d = ifelse(P_mm_d < 0, 0, P_mm_d))
    }
  }

  return(data)
}

# -----------------------------------------------------------------------------
# MAIN PROCESSING FUNCTIONS
# -----------------------------------------------------------------------------

#' Process all station pairs and triplets, apply interpolation
#' @param data Combined meteorological data
#' @param station_groups List with $pairs and $triplets from extract_station_groups()
#' @param variables Character vector of variable names to process
#' @return List with interpolated data and tracking info
process_station_groups <- function(data, station_groups, variables, plot_dir = NULL) {
  interpolated_data <- data
  interpolated_pairs <- list()
  interpolated_triplets <- list()

  # Process pairs
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

  # Process triplets
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
      # Fall back to pairwise interpolation
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

  # Apply constraints and calculate VPD
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

#' Create watershed-level datasets from interpolated station data
#' @param interpolated_data Interpolated station data
#' @param site_mapping Named list of watershed -> station mappings
#' @param variables Character vector of variable names
#' @return Named list of watershed datasets
create_watershed_datasets <- function(interpolated_data, site_mapping, variables) {
  site_datasets <- list()

  for (site_name in names(site_mapping)) {
    site_info <- site_mapping[[site_name]]
    site_data <- data.frame()

    for (var_idx in seq_along(variables)) {
      var <- variables[var_idx]

      # Handle VPD separately (uses temperature stations)
      if (var == "VPD_kPa") {
        var_name <- "temp"
      } else if (var_idx > length(names(site_info))) {
        next
      } else {
        var_name <- names(site_info)[var_idx]
      }

      stations <- site_info[[var_name]]

      if (length(stations) == 1) {
        # Single station
        single_data <- interpolated_data %>%
          filter(SITECODE == stations[1]) %>%
          select(DATE, !!sym(var))

        if (nrow(site_data) == 0) {
          site_data <- single_data %>% mutate(SITECODE = site_name)
        } else {
          site_data <- site_data %>% left_join(single_data, by = "DATE")
        }
      } else if (length(stations) >= 2) {
        # Average multiple stations
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

#' Add discharge data to watershed datasets
#' @param watershed_datasets Named list of watershed datasets
#' @param discharge Discharge data frame with DATE, SITECODE, MEAN_Q
#' @param da_df Drainage area data frame with SITECODE, DA_M2
#' @return Updated watershed datasets with Q_mm_d column
add_discharge_to_watersheds <- function(watershed_datasets, discharge, da_df) {
  discharge_processed <- discharge %>%
    left_join(da_df, by = "SITECODE") %>%
    filter(!is.na(DA_M2)) %>%
    mutate(Q_mm_d = MEAN_Q * 0.0283168 * 86400 / DA_M2 * 1000) %>%
    select(DATE, SITECODE, Q_mm_d)

  for (site_name in names(watershed_datasets)) {
    site_dis <- discharge_processed %>% filter(SITECODE == site_name)
    if (nrow(site_dis) > 0) {
      watershed_datasets[[site_name]] <- watershed_datasets[[site_name]] %>%
        left_join(site_dis, by = "DATE")
    }
  }

  return(watershed_datasets)
}

# Legacy triplet RH diagnostics (WS7MET, VANMET, H15MET)
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
