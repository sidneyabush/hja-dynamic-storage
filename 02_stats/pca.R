# Perform PCA on annual storage metrics to identify dominant.
# Inputs: No direct CSV file reads in this script.
# Author: Pamela Sullivan (original), Sidney Bush (adapted)
# Date: 2026-01-23

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)

# Clear environment
rm(list = ls())

# Source configuration (paths, site definitions, water year range)
# Get script directory (works with source() and Rscript)
# Load project config
source("config.R")


theme_set(theme_pub(base_size = 12))

# Use configuration values
site_order <- SITE_ORDER_HYDROMETRIC
base_dir   <- BASE_DATA_DIR
output_dir <- OUTPUT_DIR
output_dir <- OUT_STATS_PCA_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Remove older catchment-characteristics PCA outputs.
unlink(file.path(output_dir, c(
  "watershed_char_storage_mlr_pca_loadings.csv",
  "watershed_char_storage_mlr_pca_variance_explained.csv",
  "watershed_char_storage_mlr_pca_scores_pc1_pc2.csv"
)))

# LOAD ANNUAL STORAGE METRICS
# This file was created by 06_Aggregate_All_Metrics.R and contains all annual
# storage metrics, temperature metrics, and catchment characteristics

annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)

HJA_Yr <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

# SELECT FEATURES FOR PCA

# Core streamflow storage metrics (annual values)
# Using method abbreviations: RBI, RCS, FDC, SD, WB
# Q5norm is an ecological response, so we do not include it in storage-metric PCA.
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

# OUTLIER REMOVAL (Z-SCORE > 3)
# Remove outliers only from non-missing values

HJA_clean <- HJA_selected %>%
  filter(if_all(all_of(features), ~ is.na(.) | abs((. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)) < 3))

# IMPUTE MISSING VALUES & NORMALIZE FEATURES
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

# RUN PCA

pca_result <- prcomp(
  scaled_features %>% select(-site, -year),
  center = TRUE,
  scale. = TRUE
)

# EXTRACT PCA SCORES AND LOADINGS

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

# VARIANCE EXPLAINED BY EACH PC

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

# Remove older eco-MLR PCA outputs (eco models no longer use PCA screening).
unlink(file.path(output_dir, c(
  "storage_ecovar_mlr_pca_loadings.csv",
  "storage_ecovar_mlr_pca_variance_explained.csv",
  "storage_ecovar_mlr_pca_scores_pc1_pc2.csv"
)))

# Catchment-characteristics PCA removed by design.
# Catchment characteristics are screened directly in MLR using
# constrained predictor sets plus iterative VIF filtering.
