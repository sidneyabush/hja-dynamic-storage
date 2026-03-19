# calculate mobile-storage metrics (chs + isotope metrics) by site and water year.
# inputs: discharge_dir/hf00402_v14.csv; ec_dir/cf01201_v4.txt; ec_dir/cf00201_v7.csv; isotope_dir/mtt_fyw.csv; isotope_dir/dampingratios_2025-07-07.csv; isotope_dir/mtt_fyw_2026-03-05.xlsx.
# author: keira johnson (original chs), sidney bush
# date: 2026-03-09

library(dplyr)
library(readr)
library(lubridate)
library(readxl)

rm(list = ls())

# load project config
source("config.R")

output_dir <- OUT_MET_MOBILE_DIR
discharge_dir <- DISCHARGE_DIR
ec_dir <- EC_DIR
isotope_dir <- ISOTOPE_DIR

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- part 1: chs from chemistry tracers and discharge ----

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

calc_chs_from_tracer <- function(tracer_daily, discharge_tbl, min_obs_per_wy) {
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

  if (nrow(tracer_q_kept) == 0) {
    annual <- tibble(
      SITECODE = character(),
      waterYear = integer(),
      CHS = numeric(),
      median_bf = numeric(),
      sd_bf = numeric(),
      n_obs = integer()
    )
    return(list(annual = annual, goodyears = goodyears))
  }

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
      CHS = mean(GW_prop, na.rm = TRUE),
      median_bf = median(GW_prop, na.rm = TRUE),
      sd_bf = sd(GW_prop, na.rm = TRUE),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(
      CHS = ifelse(is.nan(CHS), NA_real_, CHS),
      median_bf = ifelse(is.nan(median_bf), NA_real_, median_bf),
      sd_bf = ifelse(is.nan(sd_bf), NA_real_, sd_bf)
    )

  list(annual = annual, goodyears = goodyears)
}

ec_file <- file.path(ec_dir, "CF01201_v4.txt")
if (!file.exists(ec_file)) {
  stop("Missing required EC file: ", ec_file)
}
ec_inst <- read_csv(ec_file, show_col_types = FALSE)

ec_daily <- ec_inst %>%
  mutate(
    date = as.Date(DATE_TIME),
    SITECODE = standardize_site_code(SITECODE),
    tracer_value = suppressWarnings(as.numeric(EC_INST))
  ) %>%
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    is.finite(tracer_value)
  ) %>%
  group_by(SITECODE, date) %>%
  summarise(tracer_value = mean(tracer_value, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(tracer_value))

ec_chs <- calc_chs_from_tracer(ec_daily, discharge, CHS_MIN_DAYS_PER_WY)
annual_bf_prop <- ec_chs$annual %>%
  rename(n_days = n_obs)

write.csv(
  annual_bf_prop,
  file.path(output_dir, "annual_gw_prop.csv"),
  row.names = FALSE
)

write.csv(
  ec_chs$goodyears %>% arrange(SITECODE, waterYear),
  file.path(OUT_MET_SUPPORT_DIR, "chs_wy_day_counts_kept.csv"),
  row.names = FALSE
)

chem_file <- file.path(ec_dir, "CF00201_v7.csv")
if (!file.exists(chem_file)) {
  stop("Missing required chemistry file: ", chem_file)
}
chem <- read_csv(chem_file, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(DATE_TIME),
    SITECODE = standardize_site_code(SITECODE),
    COND = suppressWarnings(as.numeric(COND)),
    CA = suppressWarnings(as.numeric(CA)),
    waterYear = get_water_year(date)
  ) %>%
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    waterYear >= WY_START,
    waterYear <= WY_END
  )

chem_ec_daily <- chem %>%
  group_by(SITECODE, date) %>%
  summarise(tracer_value = mean(COND, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(tracer_value))

chem_ca_daily <- chem %>%
  group_by(SITECODE, date) %>%
  summarise(tracer_value = mean(CA, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(tracer_value))

ec_chem_chs <- calc_chs_from_tracer(chem_ec_daily, discharge, CHS_MIN_OBS_PER_WY_CHEM)
ca_chs <- calc_chs_from_tracer(chem_ca_daily, discharge, CHS_MIN_OBS_PER_WY_CHEM)

annual_bf_prop_ec_chem <- ec_chem_chs$annual
annual_bf_prop_ca <- ca_chs$annual

write.csv(
  annual_bf_prop_ec_chem,
  file.path(output_dir, "annual_gw_prop_ec_chem.csv"),
  row.names = FALSE
)

write.csv(
  annual_bf_prop_ca,
  file.path(output_dir, "annual_gw_prop_ca.csv"),
  row.names = FALSE
)

write.csv(
  ec_chem_chs$goodyears %>% arrange(SITECODE, waterYear),
  file.path(OUT_MET_SUPPORT_DIR, "chs_wy_obs_counts_ec_chem_kept.csv"),
  row.names = FALSE
)

write.csv(
  ca_chs$goodyears %>% arrange(SITECODE, waterYear),
  file.path(OUT_MET_SUPPORT_DIR, "chs_wy_obs_counts_ca_kept.csv"),
  row.names = FALSE
)

comparison <- annual_bf_prop %>%
  select(
    SITECODE,
    waterYear,
    CHS_EC = CHS,
    n_days_ec = n_days
  ) %>%
  full_join(
    annual_bf_prop_ca %>%
      select(
        SITECODE,
        waterYear,
        CHS_CA = CHS,
        n_obs_ca = n_obs
      ),
    by = c("SITECODE", "waterYear")
  ) %>%
  mutate(
    diff_ca_minus_ec = CHS_CA - CHS_EC
  ) %>%
  arrange(SITECODE, waterYear)

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  cor(x[keep], y[keep])
}

safe_cor_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  cor(x[keep], y[keep], method = "spearman")
}

safe_mean <- function(x) {
  if (!any(is.finite(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

safe_rmse <- function(x, y) {
  d <- x - y
  keep <- is.finite(d)
  if (sum(keep) == 0) {
    return(NA_real_)
  }
  sqrt(mean(d[keep]^2))
}

safe_lm_intercept <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  coef(lm(y[keep] ~ x[keep]))[1]
}

safe_lm_slope <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  coef(lm(y[keep] ~ x[keep]))[2]
}

site_summary <- comparison %>%
  group_by(SITECODE) %>%
  summarise(
    n_years_ec_ca = sum(is.finite(CHS_EC) & is.finite(CHS_CA)),
    corr_ec_vs_ca = safe_cor(CHS_EC, CHS_CA),
    mean_diff_ca_minus_ec = safe_mean(diff_ca_minus_ec),
    rmse_ca_vs_ec = safe_rmse(CHS_CA, CHS_EC),
    .groups = "drop"
  )

overall_n_years_ec_ca <- sum(is.finite(comparison$CHS_EC) & is.finite(comparison$CHS_CA))

overall_summary <- tibble(
  n_years_ec_ca = overall_n_years_ec_ca,
  corr_ec_vs_ca = safe_cor(comparison$CHS_EC, comparison$CHS_CA),
  mean_diff_ca_minus_ec = safe_mean(comparison$diff_ca_minus_ec),
  rmse_ca_vs_ec = safe_rmse(comparison$CHS_CA, comparison$CHS_EC)
)

site_year_pairs_ec_ca <- comparison %>%
  filter(is.finite(CHS_EC), is.finite(CHS_CA)) %>%
  transmute(
    SITECODE,
    waterYear,
    CHS_EC,
    CHS_CA,
    diff_ca_minus_ec = CHS_CA - CHS_EC,
    abs_diff_ca_minus_ec = abs(CHS_CA - CHS_EC),
    pct_diff_ca_vs_ec = ifelse(
      abs(CHS_EC) > 0,
      100 * (CHS_CA - CHS_EC) / CHS_EC,
      NA_real_
    )
  ) %>%
  arrange(SITECODE, waterYear)

site_by_site_stats <- site_year_pairs_ec_ca %>%
  group_by(SITECODE) %>%
  summarise(
    n_years_paired = n(),
    corr_pearson = safe_cor(CHS_EC, CHS_CA),
    corr_spearman = safe_cor_spearman(CHS_EC, CHS_CA),
    mean_diff_ca_minus_ec = safe_mean(diff_ca_minus_ec),
    median_diff_ca_minus_ec = median(diff_ca_minus_ec, na.rm = TRUE),
    rmse_ca_vs_ec = safe_rmse(CHS_CA, CHS_EC),
    lm_intercept_ca_on_ec = safe_lm_intercept(CHS_EC, CHS_CA),
    lm_slope_ca_on_ec = safe_lm_slope(CHS_EC, CHS_CA),
    .groups = "drop"
  ) %>%
  arrange(SITECODE)

write.csv(
  comparison,
  file.path(output_dir, "annual_gw_prop_ec_ca_comparison.csv"),
  row.names = FALSE
)

write.csv(
  site_summary,
  file.path(output_dir, "annual_gw_prop_ec_ca_site_summary.csv"),
  row.names = FALSE
)

write.csv(
  overall_summary,
  file.path(output_dir, "annual_gw_prop_ec_ca_overall_summary.csv"),
  row.names = FALSE
)

write.csv(
  site_year_pairs_ec_ca,
  file.path(output_dir, "annual_gw_prop_ec_ca_site_year_pairs.csv"),
  row.names = FALSE
)

write.csv(
  site_by_site_stats,
  file.path(output_dir, "annual_gw_prop_ec_ca_site_by_site_stats.csv"),
  row.names = FALSE
)

# ---- part 2: isotope-derived site metrics (mtt/fyw/dr) ----

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
    # collapse period-specific legacy labels into one MTT entity.
    MTT = rowMeans(cbind(MTT_early, MTT_late), na.rm = TRUE),
    MTT = ifelse(is.nan(MTT), NA_real_, MTT),
    Fyw = suppressWarnings(as.numeric(FYWM))
  ) %>%
  select(site, MTT, Fyw) %>%
  filter(!is.na(site), site != "")

damping <- read_csv(
  file.path(isotope_dir, "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  select(site, DR = DR_Overall, DR_err = DR__err) %>%
  filter(!is.na(site), site != "", site %in% SITE_ORDER_ALL)

site_template <- tibble(site = SITE_ORDER_ALL)

isotope_metrics_mean <- site_template %>%
  left_join(mtt_fyw, by = "site") %>%
  left_join(damping, by = "site") %>%
  mutate(site = factor(site, levels = SITE_ORDER_ALL)) %>%
  arrange(site)

write.csv(
  isotope_metrics_mean,
  file.path(output_dir, "isotope_metrics_site_mean.csv"),
  row.names = FALSE
)

# ---- part 3: annual isotope metrics from segura workbook (separate path) ----

annual_isotope_segura <- tibble(
  site = character(),
  year = integer(),
  MTT = numeric(),
  Fyw = numeric()
)

if (isTRUE(USE_SEGURA_ANNUAL_ISOTOPE_METRICS)) {
  segura_file <- file.path(isotope_dir, SEGURA_ANNUAL_ISOTOPE_FILE)
  if (!file.exists(segura_file)) {
    stop("Missing required annual isotope workbook: ", segura_file)
  }

  parse_leading_numeric <- function(x) {
    x_chr <- trimws(as.character(x))
    has_num <- grepl("[-+]?[0-9]*\\.?[0-9]+", x_chr)
    out <- suppressWarnings(as.numeric(sub("^\\s*([-+]?[0-9]*\\.?[0-9]+).*$", "\\1", x_chr)))
    out[!has_num] <- NA_real_
    out
  }

  safe_excel_col <- function(df, idx) {
    if (ncol(df) < idx) {
      return(rep(NA, nrow(df)))
    }
    df[[idx]]
  }

  normalize_segura_site <- function(x) {
    out <- toupper(trimws(as.character(x)))
    out <- ifelse(grepl("^WS[0-9]$", out), paste0("WS0", sub("^WS", "", out)), out)
    out <- dplyr::case_when(
      out == "MACK" ~ "Mack",
      out == "LOOK" ~ "Look",
      TRUE ~ out
    )
    standardize_site_code(out)
  }

  mtt_raw <- read_excel(
    segura_file,
    sheet = "MTT_Segura_2021",
    skip = 5,
    col_names = FALSE,
    .name_repair = "minimal"
  )

  mtt_wide <- tibble(
    site = normalize_segura_site(safe_excel_col(mtt_raw, 1)),
    MTT_2015 = parse_leading_numeric(safe_excel_col(mtt_raw, 2)),
    MTT_2016 = parse_leading_numeric(safe_excel_col(mtt_raw, 4)),
    MTT_2017 = parse_leading_numeric(safe_excel_col(mtt_raw, 6)),
    MTT_2018 = parse_leading_numeric(safe_excel_col(mtt_raw, 8))
  ) %>%
    filter(!is.na(site), site != "")

  mtt_long <- bind_rows(
    mtt_wide %>% transmute(site, year = 2015L, MTT = MTT_2015),
    mtt_wide %>% transmute(site, year = 2016L, MTT = MTT_2016),
    mtt_wide %>% transmute(site, year = 2017L, MTT = MTT_2017),
    mtt_wide %>% transmute(site, year = 2018L, MTT = MTT_2018)
  ) %>%
    distinct(site, year, .keep_all = TRUE)

  fyw_raw <- read_excel(
    segura_file,
    sheet = "YWF_Segura_2021",
    skip = 2,
    col_names = FALSE,
    .name_repair = "minimal"
  )

  fyw_wide <- tibble(
    site = normalize_segura_site(safe_excel_col(fyw_raw, 1)),
    Fyw_2015 = parse_leading_numeric(safe_excel_col(fyw_raw, 2)),
    Fyw_2016 = parse_leading_numeric(safe_excel_col(fyw_raw, 3)),
    Fyw_2017 = parse_leading_numeric(safe_excel_col(fyw_raw, 4)),
    Fyw_2018 = parse_leading_numeric(safe_excel_col(fyw_raw, 5))
  ) %>%
    filter(!is.na(site), site != "")

  fyw_long <- bind_rows(
    fyw_wide %>% transmute(site, year = 2015L, Fyw = Fyw_2015),
    fyw_wide %>% transmute(site, year = 2016L, Fyw = Fyw_2016),
    fyw_wide %>% transmute(site, year = 2017L, Fyw = Fyw_2017),
    fyw_wide %>% transmute(site, year = 2018L, Fyw = Fyw_2018)
  ) %>%
    distinct(site, year, .keep_all = TRUE)

  annual_isotope_segura <- mtt_long %>%
    full_join(fyw_long, by = c("site", "year")) %>%
    mutate(
      site = standardize_site_code(site),
      year = as.integer(year),
      MTT = suppressWarnings(as.numeric(MTT)),
      Fyw = suppressWarnings(as.numeric(Fyw))
    ) %>%
    filter(
      site %in% SITE_ORDER_HYDROMETRIC,
      year >= WY_START,
      year <= WY_END
    ) %>%
    arrange(site, year)

  write.csv(
    annual_isotope_segura,
    file.path(output_dir, SEGURA_ANNUAL_ISOTOPE_OUTPUT_FILE),
    row.names = FALSE
  )
}

annual_isotope_summary <- annual_isotope_segura %>%
  group_by(site) %>%
  summarise(
    n_MTT_years = sum(is.finite(MTT)),
    n_Fyw_years = sum(is.finite(Fyw)),
    n_isotope_years = sum(is.finite(MTT) | is.finite(Fyw)),
    MTT = ifelse(any(is.finite(MTT)), mean(MTT, na.rm = TRUE), NA_real_),
    Fyw = ifelse(any(is.finite(Fyw)), mean(Fyw, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  select(site, MTT, Fyw, n_MTT_years, n_Fyw_years, n_isotope_years)

isotope_metrics_annual <- site_template %>%
  left_join(annual_isotope_summary, by = "site") %>%
  left_join(damping, by = "site") %>%
  mutate(site = factor(site, levels = SITE_ORDER_ALL)) %>%
  arrange(site)

write.csv(
  isotope_metrics_annual,
  file.path(output_dir, "isotope_metrics_site_annual.csv"),
  row.names = FALSE
)

isotope_metrics <- isotope_metrics_mean

write.csv(
  isotope_metrics,
  file.path(output_dir, "isotope_metrics_site.csv"),
  row.names = FALSE
)
