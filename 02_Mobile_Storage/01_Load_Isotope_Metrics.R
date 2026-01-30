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
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(tidyr)

# Clear environment
rm(list = ls())

# =============================================================================
# SITE ORDERING CONSTANTS
# =============================================================================
# Hydrometric sites (continuous streamflow - for Dynamic metrics: RBI, RCS, FDC, SD, WB)
SITE_ORDER_HYDROMETRIC <- c(
  "GSWS09",   # WS 09
  "GSWS10",   # WS 10
  "GSWS01",   # WS 01
  "GSLOOK",   # Lookout Creek
  "GSWS02",   # WS 02
  "GSWS03",   # WS 03
  "GSWS06",   # WS 06 (no chemistry data - exclude from CHS)
  "GSWS07",   # WS 07
  "GSWS08",   # WS 08
  "GSWSMC"    # Mack Creek
)

# Sites with chemistry data (for CHS/mean_bf)
SITE_ORDER_CHEMISTRY <- c(
  "GSWS09",   # WS 09
  "GSWS10",   # WS 10
  "GSWS01",   # WS 01
  "GSLOOK",   # Lookout Creek
  "GSWS02",   # WS 02
  "GSWS03",   # WS 03
  # GSWS06 excluded - no chemistry data
  "GSWS07",   # WS 07
  "GSWS08",   # WS 08
  "GSWSMC"    # Mack Creek
)

# All sites including isotope-only sites (for MTT, Fyw, DR)
# Site naming: GSWSMC=GSMACK=Mack, NC=Nostoc, MR=McRae, LC=Longer, LO1/LO2=Upper Lookout
SITE_ORDER_ALL <- c(
  "GSWS09",   # WS 09
  "GSWS10",   # WS 10
  "GSWS01",   # WS 01
  "GSLOOK",   # Lookout Creek
  "GSWS02",   # WS 02
  "GSWS03",   # WS 03
  "MR",       # McRae Creek (isotope only - no hydrometric)
  "GSWS06",   # WS 06 (hydrometric only - no chemistry/isotope)
  "GSWS07",   # WS 07
  "GSWS08",   # WS 08
  "NC",       # Nostoc Creek (isotope only - no hydrometric)
  "GSWSMC",   # Mack Creek (also GSMACK in some files)
  "LC",       # Longer Creek (isotope only - no hydrometric)
  "LO2",      # Upper Lookout Creek 2 (isotope only - no hydrometric)
  "CC",       # Cold Creek (isotope only - no hydrometric)
  "LO1"       # Upper Lookout Creek 1 (isotope only - no hydrometric)
)

# =============================================================================
# 1. SETUP: Directories
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir  <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"

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
      site == "MCRAEC" ~ "MR",        # McRae Creek
      site == "GSLOOK " ~ "GSLOOK",   # Remove trailing space
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
