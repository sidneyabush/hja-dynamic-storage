# calculate mobile storage metrics (BF + isotope metrics)

# inputs:
# discharge_dir/HF00402_v14.csv
# ec_dir/CF00201_v7.csv
# isotope_dir/MTT_FYW.csv
# isotope_dir/DampingRatios_2025-07-07.csv

# outputs:
# outputs/metrics/mobile/annual_gw_prop_ca.csv
# outputs/metrics/mobile/isotope_metrics_site.csv

# author: Keira Johnson (original BF), Sidney Bush
# date: 2026-03-09 (updated by SB)

librarian::shelf(dplyr, readr, lubridate, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

output_dir <- OUT_MET_MOBILE_DIR
discharge_dir <- DISCHARGE_DIR
ec_dir <- EC_DIR
isotope_dir <- ISOTOPE_DIR

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# BF from chemistry tracers and discharge
# load daily discharge for the sites with chemistry records
discharge <- read.csv(file.path(discharge_dir, "HF00402_v14.csv")) %>%
  mutate(
    date = as.Date(DATE, "%m/%d/%Y"),
    SITECODE = standardize_site_code(SITECODE)
  ) %>%
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    WATERYEAR >= WY_START,
    WATERYEAR <= WY_END
  ) %>%
  select(SITECODE, date, MEAN_Q, WATERYEAR)

calc_bf_from_tracer <- function(tracer_daily, discharge_tbl, min_obs_per_wy) {
  # pair tracer samples with same day discharge and keep sampled water years
  tracer_q <- tracer_daily %>%
    left_join(discharge_tbl %>% select(SITECODE, date, MEAN_Q), by = c("SITECODE", "date")) %>%
    filter(is.finite(tracer_value), is.finite(MEAN_Q), MEAN_Q > 0) %>%
    mutate(waterYear = get_water_year(date)) %>%
    filter(waterYear >= WY_START, waterYear <= WY_END)

  goodyears <- tracer_q %>%
    group_by(SITECODE, waterYear) %>%
    summarise(num_obs = n_distinct(date), .groups = "drop") %>%
    filter(num_obs >= min_obs_per_wy)

  tracer_q_kept <- tracer_q %>%
    semi_join(goodyears, by = c("SITECODE", "waterYear"))

  # return an empty table when a tracer has no usable site years
  if (nrow(tracer_q_kept) == 0) {
    return(tibble(
      SITECODE = character(),
      waterYear = integer(),
      BF = numeric(),
      median_bf = numeric(),
      sd_bf = numeric(),
      n_obs = integer()
    ))
  }

  # define runoff and baseflow endmembers from each site's Ca distribution
  tracer_q_kept <- tracer_q_kept %>%
    group_by(SITECODE) %>%
    mutate(
      tracer_runoff = quantile(tracer_value, 0.01, na.rm = TRUE),
      tracer_groundwater = quantile(tracer_value, 0.99, na.rm = TRUE),
      tracer_span = tracer_groundwater - tracer_runoff
    ) %>%
    ungroup() %>%
    mutate(
      Q_baseflow = ifelse(
        is.finite(tracer_span) & tracer_span > 0,
        MEAN_Q * (tracer_value - tracer_runoff) / tracer_span,
        NA_real_
      ),
      GW_prop = Q_baseflow / MEAN_Q,
      GW_prop = pmax(0, pmin(1, GW_prop))
    )

  annual <- tracer_q_kept %>%
    group_by(SITECODE, waterYear) %>%
    summarise(
      BF = mean(GW_prop, na.rm = TRUE),
      median_bf = median(GW_prop, na.rm = TRUE),
      sd_bf = sd(GW_prop, na.rm = TRUE),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(
      BF = ifelse(is.nan(BF), NA_real_, BF),
      median_bf = ifelse(is.nan(median_bf), NA_real_, median_bf),
      sd_bf = ifelse(is.nan(sd_bf), NA_real_, sd_bf)
    )

  annual
}

# prepare daily calcium values for chemical hydrograph separation
chem_file <- file.path(ec_dir, "CF00201_v7.csv")
chem <- read_csv(chem_file, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(DATE_TIME),
    SITECODE = standardize_site_code(SITECODE),
    CA = suppressWarnings(as.numeric(CA)),
    waterYear = get_water_year(date)
  ) %>%
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    waterYear >= WY_START,
    waterYear <= WY_END
  )

chem_ca_daily <- chem %>%
  group_by(SITECODE, date) %>%
  summarise(tracer_value = mean(CA, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(tracer_value))

annual_bf_prop_ca <- calc_bf_from_tracer(
  chem_ca_daily,
  discharge,
  BF_MIN_OBS_PER_WY_CHEM
)

write.csv(
  annual_bf_prop_ca,
  file.path(output_dir, "annual_gw_prop_ca.csv"),
  row.names = FALSE
)

# site isotope metrics (MTT, Fyw, DR)

# read MTT and Fyw values from the isotope summary table
mtt_fyw <- read_csv(
  file.path(isotope_dir, "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  mutate(
    MTT2L_val = if ("MTT2L" %in% names(.)) suppressWarnings(as.numeric(MTT2L)) else NA_real_,
    MTT2H_val = if ("MTT2H" %in% names(.)) suppressWarnings(as.numeric(MTT2H)) else NA_real_,
    MTT2M_val = if ("MTT2M" %in% names(.)) suppressWarnings(as.numeric(MTT2M)) else NA_real_
  ) %>%
  mutate(
    MTT_early = suppressWarnings(as.numeric(MTT1)),
    MTT_late = suppressWarnings(as.numeric(dplyr::coalesce(
      MTT2M_val,
      rowMeans(cbind(
        MTT2L_val,
        MTT2H_val
      ), na.rm = TRUE)
    ))),
    MTT_late = ifelse(is.nan(MTT_late), NA_real_, MTT_late),
    # collapse period specific labels into one MTT value
    MTT = rowMeans(cbind(MTT_early, MTT_late), na.rm = TRUE),
    MTT = ifelse(is.nan(MTT), NA_real_, MTT),
    Fyw = suppressWarnings(as.numeric(FYWM))
  ) %>%
  select(site, MTT, Fyw) %>%
  filter(!is.na(site), site != "")

# read damping ratio values from the compiled isotope table
damping <- read_csv(
  file.path(isotope_dir, "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  select(site, DR = DR_Overall, DR_err = DR__err) %>%
  filter(!is.na(site), site != "", site %in% SITE_ORDER_HYDROMETRIC)

site_template <- tibble(site = SITE_ORDER_HYDROMETRIC)

# keep every hydrometric site in the isotope output, even when a metric is missing
isotope_metrics_mean <- site_template %>%
  left_join(mtt_fyw, by = "site") %>%
  left_join(damping, by = "site") %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site)

write.csv(
  isotope_metrics_mean,
  file.path(output_dir, "isotope_metrics_site.csv"),
  row.names = FALSE
)
