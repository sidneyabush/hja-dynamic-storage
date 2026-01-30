# =============================================================================
# Hydrometric Storage Metrics - Time Series & Summary Plots
# =============================================================================
# Purpose: Create time series and summary (mean ± SD with Tukey letters) plots
#          for all storage metrics across HJA watersheds
#
# Plots generated:
#   - storage_<metric>_by_site_wy.png: Time series by site for each metric
#   - storage_summary_<metric>.png: Mean ± SD with Tukey HSD letters
#   - grid_all_methods.png: Combined grid of all storage metrics
#
# Inputs:
#   - HJA_StorageMetrics_Annual.csv: Annual storage metrics
#   - storage_discharge_method_annual_mm_metrics_per_site_wateryear.csv
#   - DS_drawdown_annual.csv
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(colorspace)
library(tidyr)
library(multcompView)
library(lubridate)

theme_set(theme_classic(base_size = 16))

# Clear environment
rm(list = ls())

# =============================================================================
# SOURCE CONFIGURATION
# =============================================================================

# Get script directory (works with source() and Rscript)
script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  # For Rscript, use command line args
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
# DIRECTORIES AND SETTINGS
# =============================================================================

base_dir <- BASE_DATA_DIR
output_dir <- file.path(FIGURES_DIR, "Hydrometric")

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Site order & color palette (from config.R)
site_order <- c(
  "GSWS09", "GSWS10", "GSLOOK", "GSWS01", "GSWS02",
  "GSWS03", "Mack", "GSWS06", "GSWS07", "GSWS08"
)

# Use SITE_COLORS from config, adjusted for plot order
site_cols <- c(
  "GSWS09" = "#882255",
  "GSWS10" = "#AA4499",
  "GSLOOK" = "#DDCC77",
  "GSWS01" = "#CC6677",
  "GSWS02" = "#999933",
  "GSWS03" = "#117733",
  "Mack" = "#332288",
  "GSWS06" = "#44AA99",
  "GSWS07" = "#88CCEE",
  "GSWS08" = "#6699CC"
)
site_cols <- lighten(site_cols, amount = 0.1)

# Axis-label mappings
# Axis labels using method abbreviations
axis_labels <- c(
  FDC = "Flow Duration Curve Slope (FDC)",
  CHS = "Baseflow Fraction (CHS)",
  RBI = "Richards-Baker Index (RBI)",
  RCS = "Recession Curve Slope (RCS)",
  SD = "Storage-Discharge (SD, mm)",
  WB = "Water Balance (WB, mm)"
)

# =============================================================================
# 1. LOAD STORAGE METRICS
# =============================================================================

cat("Loading storage metrics...\n")

# Primary storage metrics file
hydrometric_file <- file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual.csv")
if (!file.exists(hydrometric_file)) {
  # Try output directory
  hydrometric_file <- file.path(OUTPUT_DIR, "HJA_StorageMetrics_Annual_All.csv")
}

hydrometric_storage <- read_csv(hydrometric_file, show_col_types = FALSE)

# Clean up columns
if ("...1" %in% names(hydrometric_storage)) {
  hydrometric_storage <- hydrometric_storage %>% select(-`...1`)
}
if ("rbfi" %in% names(hydrometric_storage)) {
  hydrometric_storage <- hydrometric_storage %>% rename(RBI = rbfi)
}

# Standardize site codes
hydrometric_storage <- hydrometric_storage %>%
  mutate(
    site = case_when(
      site == "Mack" ~ "Mack",
      site == "GSLOOK_FULL" ~ "GSLOOK",
      TRUE ~ site
    )
  ) %>%
  filter(site %in% site_order) %>%
  mutate(site = factor(site, levels = site_order))

hydrometric_storage_long <- hydrometric_storage %>%
  pivot_longer(
    cols = -c(site, year),
    names_to = "metric",
    values_to = "value"
  )

# =============================================================================
# 2. LOAD STORAGE-DISCHARGE DATA
# =============================================================================

ds_discharge_file <- file.path(
  base_dir, "DynamicStorage",
  "storage_discharge_method_annual_mm_metrics_per_site_wateryear.csv"
)

if (file.exists(ds_discharge_file)) {
  ds_discharge_long <- read_csv(ds_discharge_file, show_col_types = FALSE) %>%
    select(site = site, year = wateryear, SD = S_annual_mm) %>%
    mutate(
      site = case_when(
        site == "GSLOOK_FULL" ~ "GSLOOK",
        site == "Mack" ~ "Mack",
        TRUE ~ site
      )
    ) %>%
    pivot_longer(cols = SD, names_to = "metric", values_to = "value") %>%
    mutate(site = factor(site, levels = site_order))
} else {
  ds_discharge_long <- tibble()
  cat("Warning: Storage-discharge file not found\n")
}

# =============================================================================
# 3. LOAD DRAWDOWN DATA
# =============================================================================

ds_draw_file <- file.path(base_dir, "DynamicStorage", "DS_drawdown_annual.csv")

if (file.exists(ds_draw_file)) {
  ds_draw_long <- read_csv(ds_draw_file, show_col_types = FALSE)

  # Clean up columns
  if ("...1" %in% names(ds_draw_long)) {
    ds_draw_long <- ds_draw_long %>% select(-`...1`)
  }

  ds_draw_long <- ds_draw_long %>%
    rename(site = SITECODE, year = waterYear) %>%
    mutate(
      site = case_when(
        site == "Mack" ~ "Mack",
        site == "GSLOOK_FULL" ~ "GSLOOK",
        TRUE ~ site
      )
    ) %>%
    pivot_longer(cols = WB, names_to = "metric", values_to = "value") %>%
    mutate(site = factor(site, levels = site_order))
} else {
  ds_draw_long <- tibble()
  cat("Warning: Drawdown file not found\n")
}

# =============================================================================
# 4. COMBINE ALL METRICS
# =============================================================================

storage_long <- bind_rows(
  hydrometric_storage_long,
  ds_discharge_long,
  ds_draw_long
) %>%
  filter(!is.na(value), site %in% site_order)

cat("Metrics available:", paste(unique(storage_long$metric), collapse = ", "), "\n")

# =============================================================================
# 5. TIME SERIES & SUMMARY PLOTS (LOOP)
# =============================================================================

cat("Creating individual metric plots...\n")

for (m in unique(storage_long$metric)) {
  df <- filter(storage_long, metric == m)

  if (nrow(df) < 10) {
    cat("  Skipping", m, "- insufficient data\n")
    next
  }

  axis_label <- ifelse(m %in% names(axis_labels), axis_labels[m], m)

  # Time-series by site
  p_ts <- ggplot(df, aes(year, value, color = site, fill = site)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1) +
    facet_wrap(~site, ncol = 2) +
    scale_color_manual(values = site_cols, guide = "none") +
    scale_fill_manual(values = alpha(site_cols, 0.3), guide = "none") +
    labs(x = "Water Year", y = axis_label) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(color = "black", hjust = 0)
    )

  ggsave(
    paste0("storage_", m, "_by_site_wy.png"),
    p_ts, path = output_dir,
    width = 8, height = 9, units = "in", dpi = 300
  )

  # Summary: mean ± SD + Tukey letters
  sum_df <- df %>%
    group_by(site) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(mean_val))

  # Check if we have enough groups for ANOVA
  if (nrow(sum_df) >= 3 && length(unique(df$site)) >= 3) {
    aov_r <- aov(value ~ site, data = df)
    tuk <- TukeyHSD(aov_r, "site")$site
    pvals <- setNames(tuk[, "p adj"], rownames(tuk))
    lets <- tryCatch(
      multcompLetters(pvals)$Letters,
      error = function(e) setNames(rep("a", nrow(sum_df)), sum_df$site)
    )
    sum_df$group <- lets[as.character(sum_df$site)]
  } else {
    sum_df$group <- "a"
  }

  y_max <- max(sum_df$mean_val + sum_df$sd_val, na.rm = TRUE)
  y_min <- min(sum_df$mean_val - sum_df$sd_val, na.rm = TRUE)
  offset <- 0.05 * (y_max - y_min)
  sum_df <- sum_df %>% mutate(label_y = mean_val + sd_val + offset)

  p_sum <- ggplot(sum_df, aes(site, mean_val, color = site)) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
      width = 0.2
    ) +
    geom_text(aes(label = group, y = label_y), size = 4, vjust = 0) +
    scale_x_discrete(limits = site_order) +
    scale_color_manual(values = site_cols, guide = "none") +
    labs(x = NULL, y = paste0("Mean ±1 SD of ", axis_label)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank()
    )

  ggsave(
    paste0("storage_summary_", m, ".png"),
    p_sum, path = output_dir,
    width = 6, height = 4, units = "in", dpi = 300
  )

  cat("  Created plots for:", m, "\n")
}

# =============================================================================
# 6. GRID OF ALL STORAGE METRICS
# =============================================================================

cat("Creating combined grid plot...\n")

summary_all <- storage_long %>%
  group_by(metric, site) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_val)) %>%
  group_by(metric) %>%
  do({
    dat <- .
    dfm <- filter(storage_long, metric == dat$metric[1])

    if (length(unique(dfm$site)) >= 3 && nrow(dfm) >= 10) {
      aov_r <- aov(value ~ site, data = dfm)
      tuk_r <- TukeyHSD(aov_r, "site")$site
      pv <- setNames(tuk_r[, "p adj"], rownames(tuk_r))
      let <- tryCatch(
        multcompLetters(pv)$Letters,
        error = function(e) setNames(rep("a", nrow(dat)), dat$site)
      )
      dat$group <- let[as.character(dat$site)]
    } else {
      dat$group <- "a"
    }

    y_mx <- max(dat$mean_val + dat$sd_val, na.rm = TRUE)
    y_mn <- min(dat$mean_val - dat$sd_val, na.rm = TRUE)
    dat$label_y <- dat$mean_val + dat$sd_val + 0.05 * (y_mx - y_mn)
    dat
  }) %>%
  ungroup()

# Filter to storage metrics only (exclude response variables like Q5norm, CV_Q5norm)
# Using method abbreviations: RBI, RCS, FDC, SD, WB, CHS
storage_metrics_to_plot <- c("RBI", "RCS", "FDC", "SD", "WB", "CHS")
summary_all <- summary_all %>% filter(metric %in% storage_metrics_to_plot)

p_grid <- ggplot(summary_all, aes(site, mean_val, color = site)) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.2
  ) +
  geom_text(aes(label = group, y = label_y), size = 3.5, vjust = 0) +
  facet_wrap(
    ~metric,
    ncol = 2,
    scales = "free_y",
    labeller = as_labeller(axis_labels)
  ) +
  scale_x_discrete(limits = site_order) +
  scale_color_manual(values = site_cols, guide = "none") +
  labs(x = NULL, y = "Mean ± 1 SD") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(color = "black", hjust = 0)
  )

ggsave(
  "grid_all_methods.png",
  plot = p_grid, path = output_dir,
  width = 9, height = 11, units = "in", dpi = 300
)

cat("\n=== PLOTS COMPLETE ===\n")
cat("Output directory:", output_dir, "\n")
