# -----------------------------------------------------------------------------
# Load Isotope-Derived Mobile Storage Metrics (MTT, Fyw, DR)
# -----------------------------------------------------------------------------
# Purpose: Load and standardize isotope-derived metrics for mobile storage
#
# Metrics Loaded:
#   - MTT (Mean Transit Time): Average time water spends in catchment
#   - Fyw (Young Water Fraction): Proportion of water younger than ~3 months
#   - DR (Damping Ratio): Ratio of stream to precipitation isotope amplitude
#
# Inputs:
#   - MTT_FYW.csv: Mean transit times and young water fractions by site
#   - DampingRatios_2025-07-07.csv: Isotopic damping ratios by site
#
# Outputs:
#   - Isotope_Metrics_Site.csv: Site-level isotope metrics for aggregation
#
# Note: These metrics are site-averages (not annual) because isotope sampling
#       is typically too sparse for annual resolution
# -----------------------------------------------------------------------------

# Load libraries
library(dplyr)
library(readr)
library(tidyr)

# Clear environment
rm(list = ls())

# Source configuration (paths, site definitions, water year range)
# Get script directory (works with source() and Rscript)
script_dir <- tryCatch(
  {
    dirname(sys.frame(1)$ofile)
  },
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      getwd()
    }
  }
)
if (is.null(script_dir) || script_dir == "" || script_dir == ".") {
  script_dir <- getwd()
}

config_path <- file.path(script_dir, "config.R")
if (!file.exists(config_path)) {
  config_path <- file.path(dirname(script_dir), "config.R")
}
if (!file.exists(config_path)) {
  config_path <- file.path(getwd(), "config.R")
}
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found. Please ensure config.R exists in the repo root.")
}

# -----------------------------------------------------------------------------
# SETUP: Directories (from config.R)
# -----------------------------------------------------------------------------

output_dir <- OUT_MET_MOBILE_DIR
isotope_dir <- ISOTOPE_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# LOAD MTT AND FYW DATA
# -----------------------------------------------------------------------------

mtt_fyw <- read_csv(
  file.path(isotope_dir, "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  # Select key columns: use mean values (M suffix) for MTT and Fyw
  select(
    site,
    MTT = MTTM, # Mean Transit Time (mean of methods)
    Fyw = FYWM # Young Water Fraction (mean of methods)
  ) %>%
  # Remove empty rows
  filter(!is.na(site), site != "")

# -----------------------------------------------------------------------------
# LOAD DAMPING RATIOS
# -----------------------------------------------------------------------------

damping <- read_csv(
  file.path(isotope_dir, "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = standardize_site_code(site)
  ) %>%
  # Select key columns
  select(
    site,
    DR = DR_Overall, # Overall damping ratio (average of methods)
    DR_err = DR__err # Error estimate
  )

# -----------------------------------------------------------------------------
# MERGE ISOTOPE METRICS
# -----------------------------------------------------------------------------

isotope_metrics <- mtt_fyw %>%
  full_join(damping, by = "site") %>%
  # Set factor levels for consistent ordering
  mutate(site = factor(site, levels = SITE_ORDER_ALL)) %>%
  arrange(site)

# -----------------------------------------------------------------------------
# SUMMARY: DATA AVAILABILITY
# -----------------------------------------------------------------------------

availability <- isotope_metrics %>%
  mutate(
    has_MTT = !is.na(MTT),
    has_Fyw = !is.na(Fyw),
    has_DR = !is.na(DR),
    in_hydrometric = site %in% SITE_ORDER_HYDROMETRIC
  ) %>%
  select(site, in_hydrometric, has_MTT, has_Fyw, has_DR)

print(availability)


# -----------------------------------------------------------------------------
# SAVE OUTPUT
# -----------------------------------------------------------------------------

write.csv(
  isotope_metrics,
  file.path(output_dir, "isotope_metrics_site.csv"),
  row.names = FALSE
)

# End----
