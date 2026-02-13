# -----------------------------------------------------------------------------
# ANOVA & Tukey HSD Site Comparisons
# -----------------------------------------------------------------------------
# Purpose: Test for significant differences among sites in storage metrics
#          using one-way ANOVA followed by Tukey HSD post-hoc tests
#
# Workflow:
#   Load annual storage metrics (site-year observations)
#   For each storage metric, run one-way ANOVA (site as factor)
#   If ANOVA is significant (p < 0.05), run Tukey HSD post-hoc test
#   Visualize site comparisons with boxplots and Tukey groupings
#   Save results to CSV
#
# Inputs:
#   - master_annual.csv (from aggregate_metrics.R)
#
# Outputs:
#   - anova_results.csv: ANOVA F-statistics and p-values
#   - tukey_hsd_results.csv: Pairwise site comparisons with adjusted p-values
#
# Author: Sidney Bush
# Date: 2026-01-23
# -----------------------------------------------------------------------------

# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(multcompView)

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

theme_set(theme_pub(base_size = 12))

# Use configuration values
site_order <- SITE_ORDER_HYDROMETRIC
base_dir   <- BASE_DATA_DIR
output_dir <- OUTPUT_DIR
output_dir <- OUT_STATS_ANOVA_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# LOAD ANNUAL STORAGE METRICS
# -----------------------------------------------------------------------------

annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUT_MASTER_DIR, LEGACY_ANNUAL_FILE)
}

HJA_annual <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD) %>%
  mutate(site = factor(site, levels = site_order))

# -----------------------------------------------------------------------------
# SELECT STORAGE METRICS FOR ANOVA
# -----------------------------------------------------------------------------
# NOTE: Q5norm, CV_Q5norm are NOT storage metrics - they are response variables
# Storage metrics by type (using method abbreviations):
#   Dynamic: RBI, RCS, FDC, SD
#   Mobile: CHS - note: MTT, Fyw, DR are site-level only
#   Extended Dynamic: WB

storage_metrics <- c(
  "RCS",   # RCS - Recession Curve Slope - Dynamic
  "RBI",   # RBI - Richards-Baker Index - Dynamic
  "FDC",   # FDC - Flow Duration Curve - Dynamic
  "SD",    # SD  - Storage-Discharge - Dynamic
  "CHS",   # CHS - Chemical Hydrograph Separation - Mobile
  "WB"     # WB  - Water Balance - Extended Dynamic
)

# -----------------------------------------------------------------------------
# RUN ANOVA FOR EACH METRIC
# -----------------------------------------------------------------------------

anova_results <- data.frame()

for (metric in storage_metrics) {
  # Check if metric exists in data
  if (!metric %in% colnames(HJA_annual)) next

  # Prepare data for ANOVA
  anova_data <- HJA_annual %>%
    select(site, value = all_of(metric)) %>%
    filter(!is.na(value), !is.na(site))

  # Check if enough data
  if (nrow(anova_data) < 10) next

  # Run one-way ANOVA
  anova_fit <- aov(value ~ site, data = anova_data)
  anova_summary <- summary(anova_fit)

  # Extract F-statistic and p-value
  f_stat <- anova_summary[[1]]$`F value`[1]
  p_value <- anova_summary[[1]]$`Pr(>F)`[1]

  # Store results
  anova_results <- rbind(anova_results, data.frame(
    metric = metric,
    F_statistic = f_stat,
    p_value = p_value,
    df_between = anova_summary[[1]]$Df[1],
    df_within = anova_summary[[1]]$Df[2],
    significant = ifelse(p_value < 0.05, "Yes", "No")
  ))
}

# Save ANOVA results
anova_results <- anova_results %>%
  arrange(metric)

write.csv(anova_results,
          file.path(output_dir, "anova_results.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# RUN TUKEY HSD POST-HOC TESTS
# -----------------------------------------------------------------------------

tukey_results <- data.frame()
tukey_group_letters <- data.frame()

for (metric in storage_metrics) {
  # Check if metric exists and ANOVA was significant
  if (!metric %in% colnames(HJA_annual)) next

  anova_row <- anova_results %>% filter(metric == !!metric)
  if (nrow(anova_row) == 0 || anova_row$p_value >= 0.05) next

  # Prepare data
  anova_data <- HJA_annual %>%
    select(site, value = all_of(metric)) %>%
    filter(!is.na(value), !is.na(site))

  # Run ANOVA and Tukey HSD
  anova_fit <- aov(value ~ site, data = anova_data)
  tukey_fit <- TukeyHSD(anova_fit)

  # Extract Tukey results
  tukey_df <- as.data.frame(tukey_fit$site)
  tukey_df$comparison <- rownames(tukey_df)
  tukey_df$metric <- metric

  tukey_results <- rbind(tukey_results, tukey_df %>%
    select(metric, comparison, diff, lwr, upr, p_adj = `p adj`))

  # Grouping letters by site for compact significance display
  pvals <- setNames(tukey_df$`p adj`, tukey_df$comparison)
  letters <- multcompLetters(pvals)$Letters
  letters_df <- data.frame(
    metric = metric,
    site = names(letters),
    group_letter = as.character(letters),
    stringsAsFactors = FALSE
  )
  tukey_group_letters <- rbind(tukey_group_letters, letters_df)
}

# Save Tukey HSD results
if (nrow(tukey_results) > 0) {
  tukey_results <- tukey_results %>%
    arrange(metric, comparison)
  write.csv(tukey_results,
            file.path(output_dir, "tukey_hsd_results.csv"),
            row.names = FALSE)
}

if (nrow(tukey_group_letters) > 0) {
  tukey_group_letters <- tukey_group_letters %>%
    mutate(site = factor(site, levels = site_order)) %>%
    arrange(metric, site) %>%
    mutate(site = as.character(site))
  write.csv(tukey_group_letters,
            file.path(output_dir, "tukey_group_letters.csv"),
            row.names = FALSE)
}

# -----------------------------------------------------------------------------
# SUMMARY TABLE: SITE MEANS
# -----------------------------------------------------------------------------

site_means <- HJA_annual %>%
  group_by(site) %>%
  summarise(
    across(
      all_of(storage_metrics),
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(site = factor(site, levels = site_order)) %>%
  arrange(site) %>%
  mutate(site = as.character(site))

write.csv(site_means,
          file.path(output_dir, "site_means_storage_metrics.csv"),
          row.names = FALSE)
