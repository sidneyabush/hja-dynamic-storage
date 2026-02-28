# perform pca on annual storage metrics to identify dominant.
# inputs: no direct csv file reads in this script.
# author: pamela sullivan (original), sidney bush (adapted)
# date: 2026-01-23

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)

# clear environment
rm(list = ls())

# source configuration (paths, site definitions, water year range)
# get script directory (works with source() and rscript)
# load project config
source("config.R")


theme_set(theme_pub(base_size = 12))

# use configuration values
site_order <- SITE_ORDER_HYDROMETRIC
base_dir   <- BASE_DATA_DIR
output_dir <- OUTPUT_DIR
output_dir <- OUT_STATS_PCA_DIR

# create output directory if needed
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# load annual storage metrics
# this file was created by 06_aggregate_all_metrics.r and contains all annual
# storage metrics, temperature metrics, and catchment characteristics

master_dir <- file.path(OUTPUT_DIR, "master")
annual_file <- file.path(master_dir, MASTER_ANNUAL_FILE)

HJA_Yr <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

# select features for pca

# core streamflow storage metrics (annual values)
# using method abbreviations: rbi, rcs, fdc, sd, wb
# q5norm is an ecological response, so we do not include it in storage-metric pca.
features <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD", "WB")
]

site_column <- "site"
year_column <- "year"

# select features - keep all sites even with partial data
HJA_selected <- HJA_Yr %>%
  select(all_of(c(site_column, year_column)), all_of(features))

# outlier removal (z-score > 3)
# remove outliers only from non-missing values

HJA_clean <- HJA_selected %>%
  filter(if_all(all_of(features), ~ {
    s <- sd(., na.rm = TRUE)
    if (is.na(s) || s == 0) {
      TRUE
    } else {
      is.na(.) | abs((. - mean(., na.rm = TRUE)) / s) < 3
    }
  }))

# impute missing values & normalize features
# impute missing values with column means so all sites can be included

scaled_features <- HJA_clean %>%
  select(all_of(c(site_column, year_column)), all_of(features)) %>%
  mutate(across(all_of(features), ~ {
    if (all(is.na(.))) {
      .
    } else {
      ifelse(is.na(.), mean(., na.rm = TRUE), .)
    }
  }))

# drop constant/all-na columns that cannot be scaled by pca.
feature_sds <- scaled_features %>%
  summarise(across(all_of(features), ~ sd(.x, na.rm = TRUE)))
feature_sds_vec <- as.numeric(feature_sds[1, ])
names(feature_sds_vec) <- names(feature_sds)
features_kept <- names(feature_sds_vec)[is.finite(feature_sds_vec) & feature_sds_vec > 0]

if (length(features_kept) == 0) {
  stop("PCA failed: all candidate features are constant or missing after cleaning.")
}

# run pca

pca_result <- prcomp(
  scaled_features %>% select(all_of(features_kept)),
  center = TRUE,
  scale. = TRUE
)

# extract pca scores and loadings

# pca scores (pc1 and pc2)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$site <- factor(scaled_features$site, levels = site_order)
pca_df$year <- scaled_features$year

# loadings (rotation matrix)
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$feature <- rownames(loadings)

# variance explained by each pc

# calculate variance explained
explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# create data frame (all pcs)
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

# catchment-characteristics pca removed by design.
# catchment characteristics are screened directly in mlr using
# constrained predictor sets plus iterative vif filtering.
