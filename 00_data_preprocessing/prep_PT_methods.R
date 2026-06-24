# calculate Priestley Taylor ET methods from daily catchment met+Q

# inputs:
# out_met_support_dir/catchments_met_q.csv

# outputs:
# out_met_support_dir/PT_ET_methods_timeseries.csv

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(readr, dplyr, lubridate, cran_repo = "https://cloud.r-project.org")

source("config.R")

input_file <- file.path(OUT_MET_SUPPORT_DIR, "catchments_met_q.csv")
output_dir <- OUT_MET_SUPPORT_DIR

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

all_catchments_data <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(DATE = as.Date(DATE))

required_cols <- c("DATE", "SITECODE", "T_C", "RH_d_pct", "NR_Wm2_d")
missing_cols <- setdiff(required_cols, names(all_catchments_data))

# the PT calculation needs these met columns
if (length(missing_cols) > 0) {
  stop(
    "Missing columns in catchments met file: ",
    paste(missing_cols, collapse = ", ")
  )
}

wy_start_date <- as.Date(sprintf("%d-10-01", WY_START - 1))
wy_end_date <- as.Date(sprintf("%d-09-30", WY_END))
input_min_date <- suppressWarnings(min(all_catchments_data$DATE, na.rm = TRUE))

# check date coverage before trimming to the study period
if (is.finite(input_min_date) && input_min_date > wy_start_date) {
  stop(
    "Input support data begins at ",
    as.character(input_min_date),
    ", but WY_START=",
    WY_START,
    " requires data from ",
    as.character(wy_start_date),
    ". Rebuild catchments_met_q.csv first."
  )
}

all_catchments_data <- all_catchments_data %>%
  filter(DATE >= wy_start_date, DATE <= wy_end_date)

calculate_alpha_zhang <- function(temp_celsius, rh_percent, pressure_kpa = 101.325) {
  q_specific <- 0.622 * 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3)) * (rh_percent/100) / (pressure_kpa - 0.378 * 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3)) * (rh_percent/100))
  delta <- 4098 * (0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3))) / ((temp_celsius + 237.3)^2)
  gamma <- 0.067
  Bo <- 1 # Bowen ratio term fixed
  alpha <- (1 + Bo) * (delta / (delta + gamma))
  return(alpha)
}
calculate_et_pt <- function(alpha, net_radiation_wm2, temp_celsius, rh_percent) {
  net_radiation_mjm2d <- net_radiation_wm2 * 0.0864
  delta <- 4098 * 0.6108 * exp(17.27 * temp_celsius / (temp_celsius + 237.3)) / ((temp_celsius + 237.3)^2)
  gamma <- 0.067
  lambda_v <- 2.501 - (0.002361 * temp_celsius)
  et_pt <- alpha * (delta / (delta + gamma)) * (net_radiation_mjm2d / lambda_v)
  pmax(0, et_pt)
}

all_catchments_data$alpha_zhang <- mapply(
  calculate_alpha_zhang,
  temp_celsius = all_catchments_data$T_C,
  rh_percent = all_catchments_data$RH_d_pct
)

all_catchments_data$ET_PT_zhang <- calculate_et_pt(
  alpha = all_catchments_data$alpha_zhang,
  net_radiation_wm2 = all_catchments_data$NR_Wm2_d,
  temp_celsius = all_catchments_data$T_C,
  rh_percent = all_catchments_data$RH_d_pct
)

write_csv(all_catchments_data, file.path(output_dir, "PT_ET_methods_timeseries.csv"))
