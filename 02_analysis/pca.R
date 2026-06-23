# perform PCA on annual storage metrics to identify dominant structure
# inputs: outputs/master/master_annual.csv
# outputs: outputs/models/pca/pca_loadings.csv, pca_scores_pc1_pc2.csv, pca_variance_explained.csv
# author: Pamela Sullivan (original), Sidney Bush (adapted)
# date: 2026-01-23

librarian::shelf(dplyr, readr, tidyr, ggplot2, scales, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")


theme_set(theme_pub(base_size = 12))

site_order <- SITE_ORDER_HYDROMETRIC
output_dir <- OUT_STATS_PCA_DIR

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

master_dir <- file.path(OUTPUT_DIR, "master")
annual_file <- file.path(master_dir, MASTER_ANNUAL_FILE)

HJA_Yr <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD)

# Q5norm is an ecological response, not a storage metric.
features <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD", "WB")
]

site_column <- "site"
year_column <- "year"

HJA_selected <- HJA_Yr %>%
  select(all_of(c(site_column, year_column)), all_of(features))

HJA_clean <- HJA_selected %>%
  filter(if_all(all_of(features), ~ {
    s <- sd(., na.rm = TRUE)
    if (is.na(s) || s == 0) {
      TRUE
    } else {
      is.na(.) | abs((. - mean(., na.rm = TRUE)) / s) < 3
    }
  }))

# Mean imputation keeps all site years in the storage metric PCA.
scaled_features <- HJA_clean %>%
  select(all_of(c(site_column, year_column)), all_of(features)) %>%
  mutate(across(all_of(features), ~ {
    if (all(is.na(.))) {
      .
    } else {
      ifelse(is.na(.), mean(., na.rm = TRUE), .)
    }
  }))

# drop constant/all NA columns that cannot be scaled by PCA
feature_sds <- scaled_features %>%
  summarise(across(all_of(features), ~ sd(.x, na.rm = TRUE)))
feature_sds_vec <- as.numeric(feature_sds[1, ])
names(feature_sds_vec) <- names(feature_sds)
features_kept <- names(feature_sds_vec)[is.finite(feature_sds_vec) & feature_sds_vec > 0]

if (length(features_kept) == 0) {
  stop("PCA failed: all candidate features are constant or missing after cleaning.")
}

pca_result <- prcomp(
  scaled_features %>% select(all_of(features_kept)),
  center = TRUE,
  scale. = TRUE
)

pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$site <- factor(scaled_features$site, levels = site_order)
pca_df$year <- scaled_features$year

loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$feature <- rownames(loadings)

explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)

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

# Catchment characteristics are screened directly in MLR using constrained
# predictor sets and iterative VIF filtering.
