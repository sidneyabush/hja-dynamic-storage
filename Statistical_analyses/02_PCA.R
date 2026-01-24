# =============================================================================
# Principal Component Analysis (PCA) of Storage Metrics
# =============================================================================
# Purpose: Perform PCA on annual storage metrics to identify dominant
#          patterns of covariation among hydrometric storage indicators
#
# Workflow:
#   1. Load annual storage metrics (from Correlations_Metrics.R output)
#   2. Select storage metrics for PCA
#   3. Remove outliers (z-score > 3)
#   4. Normalize/scale features
#   5. Run PCA
#   6. Visualize:
#      - PC1 vs PC2 biplot with loadings
#      - Variance explained by each PC
#
# Inputs:
#   - HJA_Stor_Temp_Yr.csv: Annual storage + temperature data
#
# Outputs:
#   - QA_PCA_biplot.png: PC1 vs PC2 with feature loadings
#   - QA_PCA_variance_explained.png: Scree plot
#
# Author: Pamela Sullivan (original), Sidney Bush (adapted)
# Date: 2026-01-23
# =============================================================================

# Load libraries
library(tidyverse)
library(ggplot2)
library(scales)
library(ggrepel)

theme_set(theme_minimal(base_size = 12))

# Clear environment
rm(list = ls())

# =============================================================================
# 1. SETUP: Directories
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir  <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 2. LOAD ANNUAL STORAGE METRICS
# =============================================================================

HJA_Yr <- read_csv(
  file.path(output_dir, "HJA_Stor_Temp_Yr.csv"),
  show_col_types = FALSE
)

# =============================================================================
# 3. SELECT FEATURES FOR PCA
# =============================================================================

# Core hydrometric storage metrics (annual values)
features <- c(
  "recession_curve_slope",
  "RBI",
  "Q5norm",
  "fdc_slope",
  "S_annual_mm",
  "DS_sum"
)

site_column <- "site"

# Select features and remove missing data
HJA_selected <- HJA_Yr %>%
  select(all_of(site_column), all_of(features)) %>%
  drop_na()

# =============================================================================
# 4. OUTLIER REMOVAL (Z-SCORE > 3)
# =============================================================================

HJA_clean <- HJA_selected %>%
  filter(if_all(all_of(features), ~ abs((. - mean(.)) / sd(.)) < 3))

# =============================================================================
# 5. NORMALIZE FEATURES
# =============================================================================

scaled_features <- HJA_clean %>%
  select(all_of(site_column), all_of(features)) %>%
  mutate(across(all_of(features), scale))

# =============================================================================
# 6. RUN PCA
# =============================================================================

pca_result <- prcomp(
  scaled_features %>% select(-site),
  center = TRUE,
  scale. = TRUE
)

# =============================================================================
# 7. EXTRACT PCA SCORES AND LOADINGS
# =============================================================================

# PCA scores (PC1 and PC2)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$site <- scaled_features$site

# Loadings (rotation matrix)
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$feature <- rownames(loadings)

# Scale loadings for visualization (multiply by 5 for visibility)
loadings_scaled <- loadings %>%
  mutate(
    PC1 = PC1 * 5,
    PC2 = PC2 * 5
  )

# =============================================================================
# 8. DEFINE SITE COLORS
# =============================================================================

site_colors <- c(
  "GSWS10" = "#AA4499",
  "GSWS09" = "#882255",
  "GSWS01" = "#CC6677",
  "GSLOOK" = "#DDCC77",
  "GSWS02" = "#999933",
  "GSWS03" = "#117733",
  "GSWSMC" = "#44AA99",
  "GSWS06" = "#88CCEE",
  "GSWS07" = "#6699CC",
  "GSWS08" = "#332288"
)

# Ensure site is a factor with the correct order
pca_df$site <- factor(pca_df$site, levels = names(site_colors))

# =============================================================================
# 9. PLOT: PC1 vs PC2 BIPLOT WITH LOADINGS
# =============================================================================

p_biplot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = site)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(name = "Site", values = site_colors) +

  # Add loading vectors
  geom_segment(
    data = loadings_scaled,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "gray30",
    inherit.aes = FALSE
  ) +

  # Add loading labels
  geom_text(
    data = loadings_scaled,
    aes(x = PC1, y = PC2, label = feature),
    color = "gray20",
    size = 3,
    vjust = 1.5,
    inherit.aes = FALSE
  ) +

  labs(
    title = "PCA Biplot: Storage Metrics",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "right"
  )

ggsave(
  file.path(output_dir, "QA_PCA_biplot.png"),
  p_biplot, width = 10, height = 8, dpi = 300
)

# =============================================================================
# 10. PLOT: VARIANCE EXPLAINED BY EACH PC
# =============================================================================

# Calculate variance explained
explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Create data frame (all PCs)
explained_df <- data.frame(
  PC = paste0("PC", 1:length(explained_var)),
  Variance_Explained = explained_var
)

# Scree plot
p_scree <- ggplot(explained_df, aes(x = PC, y = Variance_Explained)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(
    aes(label = scales::percent(Variance_Explained, accuracy = 0.1)),
    vjust = -0.5,
    size = 3.5
  ) +
  labs(
    title = "Variance Explained by Principal Components",
    x = "Principal Component",
    y = "Proportion of Variance"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(output_dir, "QA_PCA_variance_explained.png"),
  p_scree, width = 10, height = 6, dpi = 300
)
