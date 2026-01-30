# =============================================================================
# Chemical Hydrograph Separation - Baseflow Proportion Calculation
# =============================================================================
# Purpose: Use specific conductance (EC) to estimate annual mean baseflow
#          proportion via two-component hydrograph separation
#
# Method: Chemical hydrograph separation using specific conductance as tracer
#   - Groundwater endmember = 99th percentile SC (high SC = deep groundwater)
#   - Runoff endmember = 1st percentile SC (low SC = quick runoff)
#   - Daily baseflow proportion = (SC_obs - SC_runoff) / (SC_gw - SC_runoff)
#
# Timeline: Water Years with ≥365 days of concurrent SC and discharge data
#
# Inputs:
#   - CF01201_v3.txt: Continuous specific conductance (EC) data
#   - HF00402_v14.csv: Daily discharge data
#
# Outputs:
#   - Annual_GW_Prop.csv: Annual mean baseflow proportion by site
#   - QA plots: Continuous baseflow timeseries and annual summaries
#
# Author: Keira Johnson (original), Sidney Bush (adapted)
# =============================================================================

# Load libraries
library(dplyr)
library(ggplot2)
library(lubridate)

theme_set(theme_classic(base_size = 12))

# Clear environment
rm(list = ls())

# Source configuration (paths, site definitions, water year range)
# Get script directory (works with source() and Rscript)
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

config_path <- file.path(dirname(script_dir), "config.R")
if (!file.exists(config_path)) {
  config_path <- file.path(getwd(), "config.R")
}
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found. Please ensure config.R exists in the repo root.")
}

# =============================================================================
# 1. SETUP: Directories (from config.R)
# =============================================================================

base_dir <- BASE_DATA_DIR
output_dir <- OUTPUT_DIR

discharge_dir <- DISCHARGE_DIR
ec_dir <- EC_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# 2. LOAD DISCHARGE DATA
# =============================================================================

discharge <- read.csv(file.path(discharge_dir, "HF00402_v14.csv")) %>%
  mutate(date = as.Date(DATE, "%m/%d/%Y")) %>%
  # Filter to chemistry sites and water year range
  filter(
    SITECODE %in% SITE_ORDER_CHEMISTRY,
    WATERYEAR >= WY_START,
    WATERYEAR <= WY_END
  ) %>%
  select(SITECODE, date, MEAN_Q, WATERYEAR)

# =============================================================================
# 3. LOAD & PROCESS SPECIFIC CONDUCTANCE DATA
# =============================================================================

EC <- read.delim(file.path(ec_dir, "CF01201_v3.txt"), sep = ",")

# Aggregate to daily mean SC
EC_daily <- EC %>%
  mutate(date = as.Date(DATE_TIME, "%Y-%m-%d")) %>%
  # Fix site name mismatch (GSMACK → GSWSMC)
  mutate(
    SITECODE = case_when(
      SITECODE == "GSMACK" ~ "GSWSMC",
      .default = SITECODE
    )
  ) %>%
  # Filter to chemistry sites
  filter(SITECODE %in% SITE_ORDER_CHEMISTRY) %>%
  group_by(SITECODE, date) %>%
  summarise(daily_SC = mean(EC_INST, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(daily_SC))

# =============================================================================
# 4. MERGE SC AND DISCHARGE
# =============================================================================

EC_Q <- left_join(
  EC_daily,
  discharge %>% select(SITECODE, date, MEAN_Q, WATERYEAR),
  by = c("SITECODE", "date")
) %>%
  filter(!is.na(MEAN_Q), !is.na(daily_SC))

# =============================================================================
# 5. FILTER TO COMPLETE WATER YEARS
# =============================================================================

# Get water year for each record
EC_Q <- EC_Q %>%
  mutate(waterYear = get_water_year(date))

# Identify complete water years (≥365 days of data)
goodyears <- EC_Q %>%
  group_by(SITECODE, waterYear) %>%
  summarise(num_days = n_distinct(date), .groups = "drop") %>%
  filter(num_days >= 365)

# Filter to keep only complete water years
EC_Q <- EC_Q %>%
  semi_join(goodyears, by = c("SITECODE", "waterYear"))

# =============================================================================
# 6. CALCULATE BASEFLOW PROPORTION
# =============================================================================

# Define endmembers (per site, across all years)
EC_Q <- EC_Q %>%
  group_by(SITECODE) %>%
  mutate(
    SC_runoff = quantile(daily_SC, 0.01, na.rm = TRUE), # 1st percentile = quick runoff
    SC_groundwater = quantile(daily_SC, 0.99, na.rm = TRUE) # 99th percentile = groundwater
  ) %>%
  ungroup()

# Two-component mixing model
# Q_baseflow = Q_total × (SC_obs - SC_runoff) / (SC_gw - SC_runoff)
EC_Q <- EC_Q %>%
  mutate(
    Q_baseflow = MEAN_Q * (daily_SC - SC_runoff) / (SC_groundwater - SC_runoff),
    GW_prop = Q_baseflow / MEAN_Q
  )

# Clip to realistic range [0, 1]
EC_Q <- EC_Q %>%
  mutate(GW_prop = pmax(0, pmin(1, GW_prop)))

# =============================================================================
# 7. AGGREGATE TO ANNUAL MEAN BASEFLOW
# =============================================================================

annual_bf_prop <- EC_Q %>%
  group_by(SITECODE, waterYear) %>%
  summarise(
    mean_bf = mean(GW_prop, na.rm = TRUE),
    median_bf = median(GW_prop, na.rm = TRUE),
    sd_bf = sd(GW_prop, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  )

# =============================================================================
# 8. SAVE OUTPUTS
# =============================================================================

# Save annual baseflow proportions
output_file <- file.path(output_dir, "Annual_GW_Prop.csv")
write.csv(annual_bf_prop, output_file, row.names = FALSE)

# =============================================================================
# 9. QA PLOTS
# =============================================================================

# Plot 1: Continuous baseflow proportion timeseries
p1 <- ggplot(EC_Q, aes(x = date, y = GW_prop)) +
  geom_line(color = "steelblue", alpha = 0.7) +
  facet_wrap(~SITECODE, ncol = 2) +
  ylim(-0.1, 1.1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  labs(
    title = "Daily Baseflow Proportion (Chemical Hydrograph Separation)",
    subtitle = "Red dashed lines = physical bounds [0, 1]",
    x = "Date",
    y = "Baseflow Proportion"
  ) +
  theme_bw(base_size = 14)

ggsave(
  file.path(output_dir, "QA_Continuous_Baseflow_Prop.png"),
  p1,
  width = 12,
  height = 10,
  dpi = 300
)

# Plot 2: Annual mean baseflow by site and year
p2 <- ggplot(
  annual_bf_prop,
  aes(x = mean_bf, y = SITECODE, color = waterYear)
) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_c() +
  labs(
    title = "Annual Mean Baseflow Proportion by Site",
    x = "Mean Baseflow Proportion",
    y = "",
    color = "Water Year"
  ) +
  theme_bw(base_size = 14)

ggsave(
  file.path(output_dir, "QA_WY_Mean_Baseflow_Prop.png"),
  p2,
  width = 10,
  height = 6,
  dpi = 300
)

# Plot 3: Distribution of annual mean baseflow by site
p3 <- ggplot(annual_bf_prop, aes(x = SITECODE, y = mean_bf)) +
  geom_boxplot(fill = "lightblue", alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  labs(
    title = "Distribution of Annual Mean Baseflow by Site",
    x = "Site",
    y = "Mean Baseflow Proportion"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(output_dir, "QA_Baseflow_Distribution_by_Site.png"),
  p3,
  width = 10,
  height = 6,
  dpi = 300
)

# =============================================================================
# 10. SUMMARY STATISTICS
# =============================================================================

summary_stats <- annual_bf_prop %>%
  group_by(SITECODE) %>%
  summarise(
    n_years = n(),
    mean_bf_overall = mean(mean_bf, na.rm = TRUE),
    sd_bf_overall = sd(mean_bf, na.rm = TRUE),
    min_bf = min(mean_bf, na.rm = TRUE),
    max_bf = max(mean_bf, na.rm = TRUE),
    .groups = "drop"
  )
