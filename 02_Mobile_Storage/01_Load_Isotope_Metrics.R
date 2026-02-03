# =============================================================================
# Load Isotope-Derived Mobile Storage Metrics (MTT, Fyw, DR)
# =============================================================================
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
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(tidyr)

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

base_dir    <- BASE_DATA_DIR
output_dir  <- OUTPUT_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 2. LOAD MTT AND FYW DATA
# =============================================================================

mtt_fyw <- read_csv(
  file.path(base_dir, "Isotopes", "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  # Clean site names (remove trailing spaces, standardize)
  mutate(
    site = trimws(site),
    site = case_when(
      site == "MCRAEC"  ~ "MR",    # McRae Creek
      site == "GSLOOK " ~ "Look",  # Remove trailing space
      site == "GSLOOK"  ~ "Look",  # Lookout Creek
      TRUE ~ site
    )
  ) %>%
  # Select key columns: use mean values (M suffix) for MTT and Fyw
  select(
    site,
    MTT = MTTM,           # Mean Transit Time (mean of methods)
    Fyw = FYWM            # Young Water Fraction (mean of methods)
  ) %>%
  # Remove empty rows
  filter(!is.na(site), site != "")

# =============================================================================
# 3. LOAD DAMPING RATIOS
# =============================================================================

damping <- read_csv(
  file.path(base_dir, "Isotopes", "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  # Clean site names
  mutate(
    site = trimws(site),
    site = case_when(
      site == "GSMACK" ~ "GSWSMC",
      site == "MCRAEC" ~ "MCRAE",
      TRUE ~ site
    )
  ) %>%
  # Select key columns
  select(
    site,
    DR = DR_Overall,      # Overall damping ratio (average of methods)
    DR_err = DR__err      # Error estimate
  )

# =============================================================================
# 4. MERGE ISOTOPE METRICS
# =============================================================================

isotope_metrics <- mtt_fyw %>%
  full_join(damping, by = "site") %>%
  # Set factor levels for consistent ordering
  mutate(site = factor(site, levels = SITE_ORDER_ALL)) %>%
  arrange(site)

# =============================================================================
# 5. SUMMARY: DATA AVAILABILITY
# =============================================================================

cat("\n=== Isotope Metrics Data Availability ===\n\n")

availability <- isotope_metrics %>%
  mutate(
    has_MTT = !is.na(MTT),
    has_Fyw = !is.na(Fyw),
    has_DR = !is.na(DR),
    in_hydrometric = site %in% SITE_ORDER_HYDROMETRIC
  ) %>%
  select(site, in_hydrometric, has_MTT, has_Fyw, has_DR)

print(availability)

cat("\n\nSummary:\n")
cat("  Sites with MTT:", sum(availability$has_MTT), "\n")
cat("  Sites with Fyw:", sum(availability$has_Fyw), "\n")
cat("  Sites with DR:", sum(availability$has_DR), "\n")
cat("  Hydrometric sites:", sum(availability$in_hydrometric), "\n")

# =============================================================================
# 6. SAVE OUTPUT
# =============================================================================

write.csv(
  isotope_metrics,
  file.path(output_dir, "Isotope_Metrics_Site.csv"),
  row.names = FALSE
)

cat("\n\nSaved: Isotope_Metrics_Site.csv\n")

# =============================================================================
# 7. DISPLAY FINAL TABLE
# =============================================================================

cat("\n=== Isotope Metrics by Site ===\n\n")
print(isotope_metrics, n = 20)
