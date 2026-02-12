# -----------------------------------------------------------------------------
# Principal Component Analysis (PCA) of Storage Metrics
# -----------------------------------------------------------------------------
# Purpose: Perform PCA on annual storage metrics to identify dominant
#          patterns of covariation among hydrometric storage indicators
#
# Workflow:
#   1. Load annual storage metrics (from aggregate_metrics.R output)
#   2. Select storage metrics for PCA
#   3. Remove outliers (z-score > 3)
#   4. Normalize/scale features
#   5. Run PCA
#   6. Visualize:
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

theme_set(theme_minimal(base_size = 12))

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

# Use configuration values
site_order <- SITE_ORDER_HYDROMETRIC
base_dir   <- BASE_DATA_DIR
output_dir <- OUTPUT_DIR
output_dir <- OUT_STATS_PCA_DIR

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. LOAD ANNUAL STORAGE METRICS
# -----------------------------------------------------------------------------
# This file was created by 06_Aggregate_All_Metrics.R and contains all annual
# storage metrics, temperature metrics, and catchment characteristics

annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
}

HJA_Yr <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

# -----------------------------------------------------------------------------
# 3. SELECT FEATURES FOR PCA
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
# 4. OUTLIER REMOVAL (Z-SCORE > 3)
# -----------------------------------------------------------------------------
# Remove outliers only from non-missing values

HJA_clean <- HJA_selected %>%
  filter(if_all(all_of(features), ~ is.na(.) | abs((. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)) < 3))

# -----------------------------------------------------------------------------
# 5. IMPUTE MISSING VALUES & NORMALIZE FEATURES
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
# 6. RUN PCA
# -----------------------------------------------------------------------------

pca_result <- prcomp(
  scaled_features %>% select(-site, -year),
  center = TRUE,
  scale. = TRUE
)

# -----------------------------------------------------------------------------
# 7. EXTRACT PCA SCORES AND LOADINGS
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
# 9. VARIANCE EXPLAINED BY EACH PC
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

# -----------------------------------------------------------------------------
# 10. ECO-MLR PCA (PREDICTOR POOL FOR ECOLOGICAL MODELS)
# -----------------------------------------------------------------------------
# Separate PCA so ecological-model predictor structure is evaluated directly.

eco_predictor_pool <- c(
  "RCS", "RBI", "FDC", "SD", "WB", "CHS", "MTT", "Fyw", "DR",
  "Lava1_per", "Lava2_per", "Ash_Per", "Landslide_Total", "Landslide_Young"
)
eco_predictor_pool <- eco_predictor_pool[eco_predictor_pool %in% names(HJA_Yr)]

if (length(eco_predictor_pool) >= 2) {
  eco_selected <- HJA_Yr %>%
    select(all_of(c(site_column, year_column)), all_of(eco_predictor_pool))

  eco_clean <- eco_selected %>%
    filter(if_all(all_of(eco_predictor_pool), ~ is.na(.) | abs((. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)) < 3))

  eco_scaled <- eco_clean %>%
    select(all_of(c(site_column, year_column)), all_of(eco_predictor_pool)) %>%
    mutate(across(all_of(eco_predictor_pool), ~ {
      if (all(is.na(.))) {
        .
      } else {
        ifelse(is.na(.), mean(., na.rm = TRUE), .)
      }
    })) %>%
    mutate(across(all_of(eco_predictor_pool), scale))

  eco_pca <- prcomp(
    eco_scaled %>% select(-site, -year),
    center = TRUE,
    scale. = TRUE
  )

  eco_scores <- as.data.frame(eco_pca$x[, 1:2, drop = FALSE])
  eco_scores$site <- factor(eco_scaled$site, levels = site_order)
  eco_scores$year <- eco_scaled$year

  eco_loadings <- as.data.frame(eco_pca$rotation[, 1:2, drop = FALSE])
  eco_loadings$feature <- rownames(eco_loadings)

  eco_explained <- data.frame(
    PC = paste0("PC", seq_len(length(eco_pca$sdev))),
    Variance_Explained = (eco_pca$sdev^2) / sum(eco_pca$sdev^2)
  )

  write.csv(
    eco_loadings,
    file.path(output_dir, "eco_mlr_pca_loadings.csv"),
    row.names = FALSE
  )
  write.csv(
    eco_explained,
    file.path(output_dir, "eco_mlr_pca_variance_explained.csv"),
    row.names = FALSE
  )
  write.csv(
    eco_scores %>% select(site, year, PC1, PC2),
    file.path(output_dir, "eco_mlr_pca_scores_pc1_pc2.csv"),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 11. CATCHMENT-CHARACTERISTICS PCA (PREDICTOR POOL FOR STORAGE MODELS)
# -----------------------------------------------------------------------------

site_file <- file.path(OUT_MASTER_DIR, MASTER_SITE_FILE)
if (!file.exists(site_file)) {
  site_file <- file.path(OUTPUT_DIR, MASTER_SITE_FILE)
}

if (file.exists(site_file)) {
  site_df <- read_csv(site_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC)

  catch_chars_pool <- c(
    "Slope_mean", "Harvest", "Landslide_Total", "Landslide_Young",
    "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
  )
  catch_chars_pool <- catch_chars_pool[catch_chars_pool %in% names(site_df)]

  if (length(catch_chars_pool) >= 2) {
    catch_chars_df <- site_df %>%
      select(site, all_of(catch_chars_pool))

    catch_chars_scaled <- catch_chars_df %>%
      mutate(across(all_of(catch_chars_pool), ~ {
        if (all(is.na(.))) {
          .
        } else {
          ifelse(is.na(.), mean(., na.rm = TRUE), .)
        }
      })) %>%
      mutate(across(all_of(catch_chars_pool), scale))

    catch_chars_pca <- prcomp(
      catch_chars_scaled %>% select(-site),
      center = TRUE,
      scale. = TRUE
    )

    catch_chars_scores <- as.data.frame(catch_chars_pca$x[, 1:2, drop = FALSE])
    catch_chars_scores$site <- factor(catch_chars_scaled$site, levels = site_order)

    catch_chars_loadings <- as.data.frame(catch_chars_pca$rotation[, 1:2, drop = FALSE])
    catch_chars_loadings$feature <- rownames(catch_chars_loadings)

    catch_chars_explained <- data.frame(
      PC = paste0("PC", seq_len(length(catch_chars_pca$sdev))),
      Variance_Explained = (catch_chars_pca$sdev^2) / sum(catch_chars_pca$sdev^2)
    )

    write.csv(
      catch_chars_loadings,
      file.path(output_dir, "catch_chars_mlr_pca_loadings.csv"),
      row.names = FALSE
    )
    write.csv(
      catch_chars_explained,
      file.path(output_dir, "catch_chars_mlr_pca_variance_explained.csv"),
      row.names = FALSE
    )
    write.csv(
      catch_chars_scores %>% select(site, PC1, PC2),
      file.path(output_dir, "catch_chars_mlr_pca_scores_pc1_pc2.csv"),
      row.names = FALSE
    )
  }
}
