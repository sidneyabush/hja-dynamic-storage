# assemble the annual and site level metrics into master tables

# inputs:
# dynamic_dir/rbi_rcs_annual.csv
# dynamic_dir/storage_discharge_fdc_annual.csv
# dynamic_dir/fdc_slopes_overall.csv
# mobile_dir/annual_gw_prop_ca.csv
# extended_dir/ds_depletion_annual.csv
# eco_dir/stream_thermal_lowflow_metrics_annual.csv
# mobile_dir/isotope_metrics_site.csv
# catchment_characteristics_dir/catchment_char.csv

# outputs:
# outputs/master/master_annual.csv
# outputs/master/master_site.csv

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(dplyr, readr, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

dynamic_dir <- OUT_MET_DYNAMIC_DIR
mobile_dir  <- OUT_MET_MOBILE_DIR
extended_dir <- OUT_MET_EXTENDED_DIR
eco_dir <- OUT_MET_ECO_DIR
master_dir <- file.path(OUTPUT_DIR, "master")

dir.create(master_dir, recursive = TRUE, showWarnings = FALSE)

assert_unique_keys <- function(df, keys, df_name) {
  dupes <- df %>%
    count(across(all_of(keys)), name = "n") %>%
    filter(n > 1)

  # stop if repeated rows would duplicate outputs
  if (nrow(dupes) > 0) {
    stop(
      paste0(
        "Repeated rows in ", df_name, " for columns (",
        paste(keys, collapse = ", "), ")."
      )
    )
  }
}

rbi_path <- file.path(dynamic_dir, "rbi_rcs_annual.csv")
rbi_recession <- read_csv(
  rbi_path,
  show_col_types = FALSE
) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, RCS, RBI)
assert_unique_keys(rbi_recession, c("site", "year"), "rbi_recession")

fdc_path <- file.path(dynamic_dir, "storage_discharge_fdc_annual.csv")
storage_fdc <- read_csv(
  fdc_path,
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, SD, FDC, Q99, Q50, Q01, Q5norm, CV_Q5norm)
assert_unique_keys(storage_fdc, c("site", "year"), "storage_fdc")

# use the full period site level FDC in the later analysis steps
# annual site year FDC is used only for Figure 2 and the ANOVA/Tukey output
fdc_site_path <- file.path(dynamic_dir, "fdc_slopes_overall.csv")
fdc_site <- read_csv(
  fdc_site_path,
  show_col_types = FALSE
) %>%
  mutate(
    site = standardize_site_code(if ("site" %in% names(.)) site else SITECODE),
    FDC_site = suppressWarnings(as.numeric(
      if ("Slope" %in% names(.)) Slope else if ("fdc_slope" %in% names(.)) fdc_slope else if ("FDC" %in% names(.)) FDC else NA_real_
    ))
  ) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  select(site, FDC_site)
assert_unique_keys(fdc_site, c("site"), "fdc_site")

storage_fdc <- storage_fdc %>%
  left_join(fdc_site, by = "site") %>%
  mutate(FDC = dplyr::coalesce(FDC_site, FDC)) %>%
  select(-FDC_site)

# BF = annual mean baseflow fraction from calcium based hydrograph separation
chs_path <- file.path(mobile_dir, "annual_gw_prop_ca.csv")
baseflow <- read_csv(
  chs_path,
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(year = as.integer(year)) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  mutate(
    BF = suppressWarnings(as.numeric(BF)),
    BF = ifelse(site %in% BF_EXCLUDE_SITES, NA_real_, BF)
  ) %>%
  select(site, year, BF)
assert_unique_keys(baseflow, c("site", "year"), "baseflow")

wb_path <- file.path(extended_dir, "ds_depletion_annual.csv")
wb_storage <- read_csv(
  wb_path,
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  mutate(year = as.integer(year)) %>%
  mutate(site = standardize_site_code(site)) %>%
  mutate(WB = suppressWarnings(as.numeric(WB))) %>%
  filter(year >= WY_START, year <= WY_END) %>%
  select(site, year, WB)
assert_unique_keys(wb_storage, c("site", "year"), "wb_storage")

thermal_lowflow <- read_csv(
  file.path(eco_dir, "stream_thermal_lowflow_metrics_annual.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = if ("site" %in% names(.)) site else SITECODE,
    year = if ("year" %in% names(.)) year else wateryear
  ) %>%
  mutate(
    site = standardize_site_code(site),
    year = as.integer(year)
  ) %>%
  filter(year >= WY_START, year <= WY_END)

thermal_cols_required <- c(
  "T_7DMax", "Q_7Q5", "Pws", "precip_nov_may_mm"
)
thermal_cols_output <- c(
  "T_7DMax", "Q_7Q5", "Pws", "precip_nov_may_mm"
)
for (nm in thermal_cols_required) {
  if (!(nm %in% names(thermal_lowflow))) {
    thermal_lowflow[[nm]] <- NA_real_
  }
}

thermal_lowflow <- thermal_lowflow %>%
  select(
    site, year, all_of(thermal_cols_required)
  ) %>%
  mutate(
    Pws = dplyr::coalesce(Pws, precip_nov_may_mm),
    precip_nov_may_mm = dplyr::coalesce(precip_nov_may_mm, Pws)
  ) %>%
  select(
    site, year, all_of(thermal_cols_output)
  )
assert_unique_keys(thermal_lowflow, c("site", "year"), "thermal_lowflow")

HJA_annual <- rbi_recession %>%
  full_join(storage_fdc, by = c("site", "year")) %>%
  full_join(baseflow, by = c("site", "year")) %>%
  full_join(wb_storage, by = c("site", "year")) %>%
  full_join(thermal_lowflow, by = c("site", "year")) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site, year)
assert_unique_keys(HJA_annual, c("site", "year"), "HJA_annual")

write.csv(HJA_annual,
          file.path(master_dir, MASTER_ANNUAL_FILE),
          row.names = FALSE)

HJA_avg <- HJA_annual %>%
  group_by(site) %>%
  summarise(
    n_years = n(),
    across(
      where(is.numeric),
      list(mean = ~mean(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

isotope_site_path <- file.path(mobile_dir, "isotope_metrics_site.csv")
isotope_metrics <- read_csv(
  isotope_site_path,
  show_col_types = FALSE
) %>%
  mutate(
    site = if ("site" %in% names(.)) site else SITECODE
  ) %>%
  mutate(
    site = standardize_site_code(site),
    MTT = suppressWarnings(as.numeric(MTT)),
    Fyw = suppressWarnings(as.numeric(Fyw)),
    DR = suppressWarnings(as.numeric(DR))
  ) %>%
  select(site, MTT, Fyw, DR) %>%
  filter(!is.na(site), site != "", site %in% SITE_ORDER_HYDROMETRIC)
assert_unique_keys(isotope_metrics, c("site"), "isotope_metrics")

HJA_avg <- HJA_avg %>%
  left_join(isotope_metrics, by = "site")

catchment_chars <- read_csv(
  resolve_catchment_characteristics_file(),
  show_col_types = FALSE
) %>%
  rename(site = Site) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)
assert_unique_keys(catchment_chars, c("site"), "catchment_chars")

HJA_avg <- HJA_avg %>%
  left_join(catchment_chars, by = "site")

write.csv(HJA_avg,
          file.path(master_dir, MASTER_SITE_FILE),
          row.names = FALSE)
