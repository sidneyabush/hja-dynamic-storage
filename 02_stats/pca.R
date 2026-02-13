# -----------------------------------------------------------------------------
# Principal Component Analysis (PCA) of Storage Metrics
# -----------------------------------------------------------------------------
# Purpose: Perform PCA on annual storage metrics to identify dominant
#          patterns of covariation among hydrometric storage indicators
#
# Workflow:
#   Load annual storage metrics (from aggregate_metrics.R output)
#   Select storage metrics for PCA
#   Remove outliers (z-score > 3)
#   Normalize/scale features
#   Run PCA
#   Visualize:
#      - PC1 vs PC2 biplot with loadings
#      - Variance explained by each PC
#
# Inputs:
#   - master_annual.csv: Annual master data table
#
# Outputs:
#   - pca_loadings.csv: PC loadings by metric
#   - pca_variance_explained.csv: Variance explained by component
#   - pca_scores_pc1_pc2.csv: Site-year scores for PC1 and PC2
#
# Author: Pamela Sullivan (original), Sidney Bush (adapted)
# Date: 2026-01-23
# -----------------------------------------------------------------------------

# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)

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
output_dir <- OUT_STATS_PCA_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Remove deprecated catchment-characteristics PCA outputs.
unlink(file.path(output_dir, c(
  "watershed_char_storage_mlr_pca_loadings.csv",
  "watershed_char_storage_mlr_pca_variance_explained.csv",
  "watershed_char_storage_mlr_pca_scores_pc1_pc2.csv"
)))

# -----------------------------------------------------------------------------
# LOAD ANNUAL STORAGE METRICS
# -----------------------------------------------------------------------------
# This file was created by 06_Aggregate_All_Metrics.R and contains all annual
# storage metrics, temperature metrics, and catchment characteristics

annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUT_MASTER_DIR, LEGACY_ANNUAL_FILE)
}

HJA_Yr <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

# -----------------------------------------------------------------------------
# SELECT FEATURES FOR PCA
# -----------------------------------------------------------------------------

# Core hydrometric storage metrics (annual values)
# Using method abbreviations: RBI, RCS, FDC, SD, WB
# NOTE: Q5norm is a response variable, not a storage metric
features <- c(
  "RCS",
  "RBI",
  "FDC",
  "SD",
  "WB"
)

site_column <- "site"
year_column <- "year"

# Select features - keep all sites even with partial data
HJA_selected <- HJA_Yr %>%
  select(all_of(c(site_column, year_column)), all_of(features))

# -----------------------------------------------------------------------------
# OUTLIER REMOVAL (Z-SCORE > 3)
# -----------------------------------------------------------------------------
# Remove outliers only from non-missing values

HJA_clean <- HJA_selected %>%
  filter(if_all(all_of(features), ~ is.na(.) | abs((. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)) < 3))

# -----------------------------------------------------------------------------
# IMPUTE MISSING VALUES & NORMALIZE FEATURES
# -----------------------------------------------------------------------------
# Impute missing values with column means so all sites can be included

scaled_features <- HJA_clean %>%
  select(all_of(c(site_column, year_column)), all_of(features)) %>%
  mutate(across(all_of(features), ~ {
    if (all(is.na(.))) {
      .
    } else {
      ifelse(is.na(.), mean(., na.rm = TRUE), .)
    }
  })) %>%
  mutate(across(all_of(features), scale))

# -----------------------------------------------------------------------------
# RUN PCA
# -----------------------------------------------------------------------------

pca_result <- prcomp(
  scaled_features %>% select(-site, -year),
  center = TRUE,
  scale. = TRUE
)

# -----------------------------------------------------------------------------
# EXTRACT PCA SCORES AND LOADINGS
# -----------------------------------------------------------------------------

# PCA scores (PC1 and PC2)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$site <- factor(scaled_features$site, levels = site_order)
pca_df$year <- scaled_features$year

# Loadings (rotation matrix)
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$feature <- rownames(loadings)

# Scale loadings for visualization (multiply by 5 for visibility)
loadings_scaled <- loadings %>%
  mutate(
    PC1 = PC1 * 5,
    PC2 = PC2 * 5
  )


# -----------------------------------------------------------------------------
# VARIANCE EXPLAINED BY EACH PC
# -----------------------------------------------------------------------------

# Calculate variance explained
explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Create data frame (all PCs)
explained_df <- data.frame(
  PC = paste0("PC", 1:length(explained_var)),
  Variance_Explained = explained_var
)

write.csv(loadings,
          file.path(output_dir, "pca_loadings.csv"),
          row.names = FALSE)
write.csv(explained_df,
          file.path(output_dir, "pca_variance_explained.csv"),
          row.names = FALSE)
write.csv(pca_df %>% select(site, year, PC1, PC2),
          file.path(output_dir, "pca_scores_pc1_pc2.csv"),
          row.names = FALSE)

# Remove deprecated eco-MLR PCA outputs (eco models no longer use PCA screening).
unlink(file.path(output_dir, c(
  "storage_ecovar_mlr_pca_loadings.csv",
  "storage_ecovar_mlr_pca_variance_explained.csv",
  "storage_ecovar_mlr_pca_scores_pc1_pc2.csv"
)))

# Catchment-characteristics PCA removed by design.
# Catchment characteristics are screened directly in MLR using
# constrained predictor sets plus iterative VIF filtering.
