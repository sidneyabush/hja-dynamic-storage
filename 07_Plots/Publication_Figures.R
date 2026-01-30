# =============================================================================
# Publication Figures for HJA Dynamic Storage Manuscript
# =============================================================================
# Purpose: Generate all publication-ready figures for the storage manuscript
#
# Figures Generated:
#   Fig 1: Study site map (placeholder - requires spatial data)
#   Fig 2: Storage metrics comparison across watersheds (grid with Tukey)
#   Fig 3: Correlation matrix of storage metrics
#   Fig 4: PCA biplot of storage metrics
#   Fig 5: Storage-thermal/low-flow relationships (hypothesis testing)
#   Fig 6: Storage metric time series by site
#
# Inputs:
#   - HJA_StorageMetrics_Annual_All.csv: Annual storage metrics
#   - HJA_Ave_StorageMetrics_CatCharacter.csv: Site-averaged metrics
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(colorspace)
library(patchwork)
library(ggcorrplot)
library(ggrepel)
library(scales)

# Clear environment
rm(list = ls())

# =============================================================================
# SOURCE CONFIGURATION
# =============================================================================

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
  script_dir <- file.path(getwd(), "07_Plots")
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
# SETUP
# =============================================================================

base_dir <- BASE_DATA_DIR
output_dir <- file.path(FIGURES_DIR, "Publication")

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Site order and colors
site_order <- SITE_ORDER_HYDROMETRIC
site_cols <- SITE_COLORS

# Metric labels for publication
metric_labels <- c(
  "RBI" = "Richards-Baker Index",
  "recession_curve_slope" = "Recession Slope",
  "fdc_slope" = "FDC Slope",
  "S_annual_mm" = "Storage-Discharge (mm)",
  "DS_sum" = "Drawdown (mm)",
  "mean_bf" = "Baseflow Fraction",
  "MTT" = "Mean Transit Time (yr)",
  "Fyw" = "Young Water Fraction",
  "DR" = "Damping Ratio",
  "Q5norm" = "Normalized Q5",
  "CV_Q5norm" = "CV of Normalized Q5"
)

# Theme for publication
theme_pub <- theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12)
  )

theme_set(theme_pub)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

# Annual storage metrics
annual_file <- file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUTPUT_DIR, "HJA_StorageMetrics_Annual_All.csv")
}

if (file.exists(annual_file)) {
  annual_data <- read_csv(annual_file, show_col_types = FALSE) %>%
    filter(!site %in% c("GSLOOK_FULL", "GSWSMA", "GSWSMF", "GSMACK")) %>%
    mutate(site = factor(site, levels = site_order))
  cat("  Loaded annual data:", nrow(annual_data), "rows\n")
} else {
  stop("Annual storage metrics file not found")
}

# Site-averaged metrics with catchment characteristics
avg_file <- file.path(OUTPUT_DIR, "HJA_Ave_StorageMetrics_CatCharacter.csv")
if (!file.exists(avg_file)) {
  avg_file <- file.path(base_dir, "DynamicStorage", "HJA_Ave_StorageMetrics_CatCharacter.csv")
}

if (file.exists(avg_file)) {
  avg_data <- read_csv(avg_file, show_col_types = FALSE) %>%
    filter(!site %in% c("GSLOOK_FULL", "GSWSMA", "GSWSMF", "GSMACK")) %>%
    mutate(site = factor(site, levels = site_order))
  cat("  Loaded site-averaged data:", nrow(avg_data), "rows\n")
} else {
  cat("  Warning: Site-averaged file not found\n")
  avg_data <- NULL
}

# =============================================================================
# FIGURE 2: STORAGE METRICS COMPARISON (GRID WITH MEANS)
# =============================================================================

cat("\nCreating Figure 2: Storage metrics comparison...\n")

# Select storage metrics for comparison
storage_vars <- c("RBI", "recession_curve_slope", "fdc_slope",
                  "S_annual_mm", "DS_sum", "mean_bf")

# Calculate site means and SDs
summary_data <- annual_data %>%
  select(site, year, all_of(intersect(storage_vars, names(annual_data)))) %>%
  pivot_longer(cols = -c(site, year), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(site, metric) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    metric_label = metric_labels[metric],
    metric = factor(metric, levels = storage_vars)
  )

# Create grid plot
fig2 <- ggplot(summary_data, aes(x = site, y = mean_val, color = site)) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.3, linewidth = 0.5
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2,
             labeller = labeller(metric = metric_labels)) +
  scale_color_manual(values = site_cols, guide = "none") +
  labs(x = NULL, y = "Mean ± 1 SD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(
  file.path(output_dir, "Fig2_Storage_Metrics_Comparison.png"),
  fig2, width = 8, height = 10, dpi = 300
)
ggsave(
  file.path(output_dir, "Fig2_Storage_Metrics_Comparison.pdf"),
  fig2, width = 8, height = 10
)

cat("  Saved Figure 2\n")

# =============================================================================
# FIGURE 3: CORRELATION MATRIX
# =============================================================================

cat("Creating Figure 3: Correlation matrix...\n")

if (!is.null(avg_data)) {
  # Select numeric columns for correlation
  cor_vars <- c("RBI_mean", "recession_curve_slope_mean", "fdc_slope_mean",
                "S_annual_mm_mean", "DS_sum_mean", "mean_bf_mean")
  cor_vars <- intersect(cor_vars, names(avg_data))

  if (length(cor_vars) >= 3) {
    cor_matrix <- cor(avg_data[cor_vars], use = "pairwise.complete.obs")

    # Clean up variable names for display
    rownames(cor_matrix) <- gsub("_mean", "", rownames(cor_matrix))
    colnames(cor_matrix) <- gsub("_mean", "", colnames(cor_matrix))

    fig3 <- ggcorrplot(
      cor_matrix,
      type = "lower",
      outline.color = "white",
      lab = TRUE,
      lab_size = 3.5,
      colors = c("#6D9EC1", "white", "#E46726"),
      title = "Storage Metrics Correlation Matrix"
    ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    ggsave(
      file.path(output_dir, "Fig3_Correlation_Matrix.png"),
      fig3, width = 8, height = 7, dpi = 300
    )
    ggsave(
      file.path(output_dir, "Fig3_Correlation_Matrix.pdf"),
      fig3, width = 8, height = 7
    )

    cat("  Saved Figure 3\n")
  } else {
    cat("  Skipping Figure 3: Insufficient correlation variables\n")
  }
} else {
  cat("  Skipping Figure 3: No site-averaged data\n")
}

# =============================================================================
# FIGURE 4: PCA BIPLOT
# =============================================================================

cat("Creating Figure 4: PCA biplot...\n")

# PCA on annual data
pca_vars <- c("RBI", "recession_curve_slope", "fdc_slope", "S_annual_mm", "DS_sum")
pca_vars <- intersect(pca_vars, names(annual_data))

if (length(pca_vars) >= 3) {
  pca_data <- annual_data %>%
    select(site, year, all_of(pca_vars)) %>%
    na.omit()

  if (nrow(pca_data) >= 20) {
    # Scale and run PCA
    pca_scaled <- pca_data %>%
      mutate(across(all_of(pca_vars), scale))

    pca_result <- prcomp(pca_scaled[pca_vars], center = TRUE, scale. = TRUE)

    # Extract scores
    pca_scores <- as.data.frame(pca_result$x[, 1:2])
    pca_scores$site <- pca_data$site

    # Extract loadings
    loadings <- as.data.frame(pca_result$rotation[, 1:2])
    loadings$feature <- rownames(loadings)
    loadings <- loadings %>%
      mutate(
        PC1_scaled = PC1 * 4,
        PC2_scaled = PC2 * 4,
        label = metric_labels[feature]
      )

    # Variance explained
    var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

    fig4 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = site)) +
      geom_point(size = 2, alpha = 0.7) +
      geom_segment(
        data = loadings,
        aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "gray30", linewidth = 0.6,
        inherit.aes = FALSE
      ) +
      geom_text_repel(
        data = loadings,
        aes(x = PC1_scaled, y = PC2_scaled, label = feature),
        color = "gray20", size = 3,
        inherit.aes = FALSE
      ) +
      scale_color_manual(values = site_cols, name = "Site") +
      labs(
        title = "PCA of Annual Storage Metrics",
        x = paste0("PC1 (", var_explained[1], "%)"),
        y = paste0("PC2 (", var_explained[2], "%)")
      ) +
      coord_fixed() +
      theme(legend.position = "right")

    ggsave(
      file.path(output_dir, "Fig4_PCA_Biplot.png"),
      fig4, width = 9, height = 7, dpi = 300
    )
    ggsave(
      file.path(output_dir, "Fig4_PCA_Biplot.pdf"),
      fig4, width = 9, height = 7
    )

    cat("  Saved Figure 4\n")
  } else {
    cat("  Skipping Figure 4: Insufficient data for PCA\n")
  }
} else {
  cat("  Skipping Figure 4: Insufficient PCA variables\n")
}

# =============================================================================
# FIGURE 5: STORAGE-RESPONSE RELATIONSHIPS (HYPOTHESIS TESTING)
# =============================================================================

cat("Creating Figure 5: Storage-response relationships...\n")

# Check if response variables exist
response_vars <- c("max_temp_7d_C", "min_Q_7d_mm_d", "temp_at_min_Q_7d_C")
response_present <- intersect(response_vars, names(annual_data))

if (length(response_present) >= 1) {
  # Create scatterplots for key relationships
  plot_list <- list()

  # Storage predictor to use
  predictor <- "DS_sum"  # Drawdown as main predictor
  if (!predictor %in% names(annual_data)) {
    predictor <- "S_annual_mm"
  }

  for (resp in response_present) {
    plot_data <- annual_data %>%
      select(site, all_of(c(predictor, resp))) %>%
      na.omit()

    if (nrow(plot_data) >= 10) {
      # Fit linear model
      lm_fit <- lm(as.formula(paste(resp, "~", predictor)), data = plot_data)
      r2 <- round(summary(lm_fit)$r.squared, 3)
      p_val <- summary(lm_fit)$coefficients[2, 4]
      p_text <- ifelse(p_val < 0.001, "p < 0.001",
                       ifelse(p_val < 0.01, paste0("p = ", round(p_val, 3)),
                              paste0("p = ", round(p_val, 2))))

      p <- ggplot(plot_data, aes(x = .data[[predictor]], y = .data[[resp]])) +
        geom_point(aes(color = site), size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
        scale_color_manual(values = site_cols, guide = "none") +
        labs(
          x = metric_labels[predictor],
          y = metric_labels[resp],
          subtitle = paste0("R² = ", r2, ", ", p_text)
        ) +
        theme(plot.subtitle = element_text(size = 9))

      plot_list[[resp]] <- p
    }
  }

  if (length(plot_list) >= 1) {
    fig5 <- wrap_plots(plot_list, ncol = min(length(plot_list), 2)) +
      plot_annotation(
        title = "Storage-Response Relationships",
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )

    ggsave(
      file.path(output_dir, "Fig5_Storage_Response_Relationships.png"),
      fig5, width = 10, height = 4 * ceiling(length(plot_list) / 2), dpi = 300
    )
    ggsave(
      file.path(output_dir, "Fig5_Storage_Response_Relationships.pdf"),
      fig5, width = 10, height = 4 * ceiling(length(plot_list) / 2)
    )

    cat("  Saved Figure 5\n")
  } else {
    cat("  Skipping Figure 5: Insufficient data for response plots\n")
  }
} else {
  cat("  Skipping Figure 5: No response variables found\n")
}

# =============================================================================
# FIGURE 6: TIME SERIES PANELS
# =============================================================================

cat("Creating Figure 6: Time series panels...\n")

# Key metrics for time series
ts_vars <- c("RBI", "DS_sum", "S_annual_mm")
ts_vars <- intersect(ts_vars, names(annual_data))

if (length(ts_vars) >= 1) {
  ts_data <- annual_data %>%
    select(site, year, all_of(ts_vars)) %>%
    pivot_longer(cols = all_of(ts_vars), names_to = "metric", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(metric_label = metric_labels[metric])

  fig6 <- ggplot(ts_data, aes(x = year, y = value, color = site)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    facet_grid(metric ~ site, scales = "free_y",
               labeller = labeller(metric = metric_labels)) +
    scale_color_manual(values = site_cols, guide = "none") +
    labs(x = "Water Year", y = "Value") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      strip.text.y = element_text(angle = 0),
      panel.spacing = unit(0.3, "lines")
    )

  ggsave(
    file.path(output_dir, "Fig6_Time_Series_Panels.png"),
    fig6, width = 14, height = 8, dpi = 300
  )
  ggsave(
    file.path(output_dir, "Fig6_Time_Series_Panels.pdf"),
    fig6, width = 14, height = 8
  )

  cat("  Saved Figure 6\n")
}

# =============================================================================
# SUPPLEMENTARY: CATCHMENT CHARACTERISTICS
# =============================================================================

cat("Creating supplementary figures...\n")

if (!is.null(avg_data)) {
  # Catchment characteristic variables
  catch_vars <- c("Area_km2", "Elevation_mean_m", "Slope_mean", "Harvest")
  catch_vars <- intersect(catch_vars, names(avg_data))

  if (length(catch_vars) >= 2) {
    catch_long <- avg_data %>%
      select(site, all_of(catch_vars)) %>%
      pivot_longer(cols = -site, names_to = "variable", values_to = "value") %>%
      filter(!is.na(value))

    fig_catch <- ggplot(catch_long, aes(x = site, y = value, fill = site)) +
      geom_col() +
      facet_wrap(~variable, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = site_cols, guide = "none") +
      labs(x = NULL, y = "Value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(
      file.path(output_dir, "FigS1_Catchment_Characteristics.png"),
      fig_catch, width = 10, height = 8, dpi = 300
    )

    cat("  Saved supplementary catchment figure\n")
  }
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== PUBLICATION FIGURES COMPLETE ===\n")
cat("Output directory:", output_dir, "\n")
cat("Files created:\n")
list.files(output_dir, pattern = "\\.(png|pdf)$") %>%
  paste0("  - ", .) %>%
  cat(sep = "\n")
