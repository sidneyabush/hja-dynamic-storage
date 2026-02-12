# -----------------------------------------------------------------------------
# Main Text and Supplement Figures - HJA Dynamic Storage
# -----------------------------------------------------------------------------
#
# Main Text Figures:
#   Figure 3: Dynamic Storage (RBI, RCS, FDC, SD) - faceted
#   Figure 4: Mobile Isotope (MTT, Fyw, DR) - faceted
#   Figure 5: CHS (baseflow fraction) - single panel
#
# Supplement Figures:
#   - Faceted time series for Dynamic Storage metrics
#   - Faceted time series for Mobile Storage metrics
#
# Author: Sidney Bush
# Date: 2026-02-02
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(colorspace)

rm(list = ls())

# -----------------------------------------------------------------------------
# SOURCE CONFIGURATION
# -----------------------------------------------------------------------------

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
source(config_path)

# -----------------------------------------------------------------------------
# SETUP
# -----------------------------------------------------------------------------

base_dir <- BASE_DATA_DIR
main_dir <- file.path(FIGURES_DIR, "main")
supp_dir <- file.path(FIGURES_DIR, "supp")

for (d in c(main_dir, supp_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

site_order <- SITE_ORDER_HYDROMETRIC
site_colors <- SITE_COLORS

# Metric labels
metric_labels <- c(
  "RBI" = "RBI", "RCS" = "RCS", "FDC" = "FDC", "SD" = "SD (mm)",
  "WB" = "WB (mm)", "CHS" = "CHS", "MTT" = "MTT (yr)",
  "Fyw" = "Fyw", "DR" = "DR"
)

# Publication theme
theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, hjust = 0),
    legend.position = "none"
  )

theme_set(theme_pub)

# -----------------------------------------------------------------------------
# HELPER: Standardize site names to WS## format
# -----------------------------------------------------------------------------

standardize_sites <- function(df) {
  df %>%
    mutate(
      site = standardize_site_code(site)
    ) %>%
    filter(site %in% site_order) %>%
    mutate(site = factor(site, levels = site_order))
}

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

cat("Loading data...\n")

# Annual storage metrics
annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
}

annual_data <- read_csv(annual_file, show_col_types = FALSE) %>%
  rename_legacy_storage_metrics() %>%
  standardize_sites()

cat("  Loaded annual data:", nrow(annual_data), "rows\n")

# Isotope data for mobile storage (site-level metrics)
isotope_file <- file.path(OUT_MET_MOBILE_DIR, "isotope_metrics_site.csv")
if (!file.exists(isotope_file)) {
  isotope_file <- file.path(OUT_MET_MOBILE_DIR, "Isotope_Metrics_Site.csv")
}
if (!file.exists(isotope_file)) {
  isotope_file <- file.path(base_dir, "Isotopes", "Isotope_Metrics_Site.csv")
}

if (file.exists(isotope_file)) {
  isotope_data <- read_csv(isotope_file, show_col_types = FALSE) %>%
    standardize_sites()

  # Add MTT1/MTT2 columns from raw isotope table when available.
  mtt_file <- file.path(base_dir, "Isotopes", "MTT_FYW.csv")
  if (file.exists(mtt_file)) {
    mtt_versions <- read_csv(mtt_file, show_col_types = FALSE) %>%
      mutate(site = standardize_site_code(site)) %>%
      filter(site %in% site_order) %>%
      mutate(
        MTT2 = if ("MTT2M" %in% names(.)) {
          MTT2M
        } else if (all(c("MTT2L", "MTT2H") %in% names(.))) {
          rowMeans(cbind(MTT2L, MTT2H), na.rm = TRUE)
        } else {
          NA_real_
        }
      ) %>%
      select(site, MTT1, MTT2)

    isotope_data <- isotope_data %>%
      left_join(mtt_versions, by = "site")
  }

  if (!("MTT1" %in% names(isotope_data))) isotope_data$MTT1 <- NA_real_
  if (!("MTT2" %in% names(isotope_data))) {
    isotope_data$MTT2 <- if ("MTT" %in% names(isotope_data)) isotope_data$MTT else NA_real_
  }
  cat("  Loaded isotope data:", nrow(isotope_data), "rows\n")
} else {
  isotope_data <- NULL
  cat("  Warning: Isotope data not found\n")
}

# CHS data
chs_file <- file.path(base_dir, "EC", "CHS_annual_baseflow.csv")
if (!file.exists(chs_file)) {
  chs_file <- file.path(OUT_MET_MOBILE_DIR, "annual_gw_prop.csv")
}
if (!file.exists(chs_file)) {
  chs_file <- file.path(OUT_MET_MOBILE_DIR, "Annual_GW_Prop.csv")
}

if (file.exists(chs_file)) {
  chs_data <- read_csv(chs_file, show_col_types = FALSE) %>%
    mutate(
      site = if ("site" %in% names(.)) site else if ("SITECODE" %in% names(.)) SITECODE else NA_character_,
      year = if ("year" %in% names(.)) year else if ("waterYear" %in% names(.)) waterYear else NA_real_,
      CHS = if ("CHS" %in% names(.)) CHS else if ("mean_bf" %in% names(.)) mean_bf else CHS
    ) %>%
    select(site, year, CHS) %>%
    standardize_sites()
  cat("  Loaded CHS data:", nrow(chs_data), "rows\n")
} else {
  # Try to get CHS from annual data
  if ("CHS" %in% names(annual_data)) {
    chs_data <- annual_data %>%
      select(site, year, CHS) %>%
      filter(!is.na(CHS))
    cat("  Using CHS from annual data:", nrow(chs_data), "rows\n")
  } else {
    chs_data <- NULL
    cat("  Warning: CHS data not found\n")
  }
}

# -----------------------------------------------------------------------------
# FIGURE 3: DYNAMIC STORAGE (RBI, RCS, FDC, SD)
# -----------------------------------------------------------------------------

cat("\nCreating Figure 3: Dynamic Storage...\n")

dynamic_metrics <- c("RBI", "RCS", "FDC", "SD")

dynamic_long <- annual_data %>%
  select(site, year, any_of(dynamic_metrics)) %>%
  pivot_longer(cols = -c(site, year), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric, levels = dynamic_metrics))

# Calculate mean ± SD per site
dynamic_summary <- dynamic_long %>%
  group_by(site, metric) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    site = factor(site_order, levels = site_order),
    metric = factor(dynamic_metrics, levels = dynamic_metrics),
    fill = list(mean_val = NA_real_, sd_val = NA_real_)
  )

# Add mean labels for each facet
dynamic_labels <- dynamic_summary %>%
  group_by(metric) %>%
  summarise(
    x_pos = -Inf,
    y_pos = -Inf,
    .groups = "drop"
  )

fig3 <- ggplot(dynamic_summary, aes(x = site, y = mean_val, color = site)) +
  geom_point(size = 2.5, na.rm = TRUE) +
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.3, linewidth = 0.5, na.rm = TRUE
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2,
             labeller = labeller(metric = metric_labels)) +
  scale_color_manual(values = site_colors) +
  scale_x_discrete(limits = site_order, drop = FALSE) +
  labs(x = NULL, y = "Mean ± 1 SD")

ggsave(file.path(main_dir, "ds_summary.png"), fig3, width = 10, height = 8, dpi = 300)
ggsave(file.path(main_dir, "ds_summary.pdf"), fig3, width = 10, height = 8)
cat("  Saved Figure 3\n")

# -----------------------------------------------------------------------------
# FIGURE 4: MOBILE ISOTOPE (MTT, Fyw, DR)
# -----------------------------------------------------------------------------

cat("Creating Figure 4: Mobile Isotope Storage...\n")

if (!is.null(isotope_data) && nrow(isotope_data) > 0) {
  make_mobile_panel <- function(df, y_col, y_lab, err_col = NULL) {
    p <- ggplot(df, aes(x = site, y = .data[[y_col]], color = site)) +
      geom_point(size = 3, na.rm = TRUE) +
      scale_color_manual(values = site_colors) +
      scale_x_discrete(limits = site_order, drop = FALSE) +
      labs(x = NULL, y = y_lab)

    if (!is.null(err_col) && err_col %in% names(df)) {
      p <- p + geom_errorbar(
        aes(ymin = .data[[y_col]] - .data[[err_col]], ymax = .data[[y_col]] + .data[[err_col]]),
        width = 0.3, linewidth = 0.5, na.rm = TRUE
      )
    }
    p
  }

  p_mtt1 <- make_mobile_panel(isotope_data, "MTT1", "MTT1 (yr)")
  p_mtt2 <- make_mobile_panel(isotope_data, "MTT2", "MTT2 (yr)")
  p_fyw <- make_mobile_panel(isotope_data, "Fyw", "Fyw")
  p_dr <- make_mobile_panel(isotope_data, "DR", "DR", err_col = "DR_err")

  fig4 <- (p_mtt1 | p_mtt2) / (p_fyw | p_dr)
  ggsave(file.path(main_dir, "ms_isotope.png"), fig4, width = 12, height = 8, dpi = 300)
  ggsave(file.path(main_dir, "ms_isotope.pdf"), fig4, width = 12, height = 8)
  cat("  Saved Figure 4\n")
} else {
  cat("  Skipping Figure 4: No isotope data\n")
}

# -----------------------------------------------------------------------------
# FIGURE 5: CHS (BASEFLOW FRACTION)
# -----------------------------------------------------------------------------

cat("Creating Figure 5: CHS...\n")

if (!is.null(chs_data) && nrow(chs_data) > 0) {
  # Summary stats
  chs_summary <- chs_data %>%
    group_by(site) %>%
    summarise(
      mean_val = mean(CHS, na.rm = TRUE),
      sd_val = sd(CHS, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    complete(
      site = factor(site_order, levels = site_order),
      fill = list(mean_val = NA_real_, sd_val = NA_real_, n = 0)
    )

  fig5 <- ggplot(chs_summary, aes(x = site, y = mean_val, color = site)) +
    geom_point(size = 3, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
      width = 0.3, linewidth = 0.6, na.rm = TRUE
    ) +
    scale_color_manual(values = site_colors) +
    scale_x_discrete(limits = site_order, drop = FALSE) +
    labs(x = NULL, y = "Baseflow Fraction (CHS)")

  ggsave(file.path(main_dir, "ms_chs.png"), fig5, width = 8, height = 5, dpi = 300)
  ggsave(file.path(main_dir, "ms_chs.pdf"), fig5, width = 8, height = 5)

  fig5_box <- ggplot(chs_data, aes(x = site, y = CHS, color = site, fill = site)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2, na.rm = TRUE) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.4, alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = site_colors) +
    scale_fill_manual(values = site_colors) +
    scale_x_discrete(limits = site_order, drop = FALSE) +
    labs(x = NULL, y = "Baseflow Fraction (CHS)")

  ggsave(file.path(main_dir, "ms_chs_boxplot.png"), fig5_box, width = 8, height = 5, dpi = 300)
  ggsave(file.path(main_dir, "ms_chs_boxplot.pdf"), fig5_box, width = 8, height = 5)
  cat("  Saved Figure 5\n")
} else {
  cat("  Skipping Figure 5: No CHS data\n")
}

# -----------------------------------------------------------------------------
# SUPPLEMENT: FACETED TIME SERIES - DYNAMIC STORAGE
# -----------------------------------------------------------------------------

cat("\nCreating Supplement: Dynamic Storage Time Series...\n")

# Calculate site means for labels
site_means_dynamic <- dynamic_long %>%
  group_by(site, metric) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    x_pos = min(year, na.rm = TRUE),
    y_pos = max(value, na.rm = TRUE),
    .groups = "drop"
  )

dynamic_line_data <- dynamic_long %>%
  group_by(site, metric) %>%
  filter(sum(!is.na(value)) >= 2) %>%
  ungroup()

supp_dynamic_ts <- ggplot(dynamic_long, aes(x = year, y = value, color = site)) +
  geom_line(data = dynamic_line_data, linewidth = 0.5, na.rm = TRUE) +
  geom_point(size = 1, na.rm = TRUE) +
  geom_text(
    data = site_means_dynamic,
    aes(x = x_pos, y = y_pos, label = sprintf("%.2f", mean_val)),
    hjust = 0, vjust = -0.3, size = 2.5, color = "black", na.rm = TRUE
  ) +
  facet_grid(metric ~ site, scales = "free_y",
             labeller = labeller(metric = metric_labels), drop = FALSE) +
  scale_color_manual(values = site_colors) +
  labs(x = "Water Year", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text.y = element_text(angle = 0, size = 9),
    strip.text.x = element_text(size = 8)
  )

ggsave(file.path(supp_dir, "ds_annual_ts.png"), supp_dynamic_ts, width = 14, height = 10, dpi = 300)
ggsave(file.path(supp_dir, "ds_annual_ts.pdf"), supp_dynamic_ts, width = 14, height = 10)
cat("  Saved Dynamic Storage Time Series\n")

# -----------------------------------------------------------------------------
# SUPPLEMENT: FACETED TIME SERIES - CHS
# -----------------------------------------------------------------------------

if (!is.null(chs_data) && nrow(chs_data) > 0) {
  cat("Creating Supplement: CHS Time Series...\n")

  chs_means <- chs_data %>%
    group_by(site) %>%
    summarise(
      mean_val = mean(CHS, na.rm = TRUE),
      x_pos = min(year, na.rm = TRUE),
      y_pos = max(CHS, na.rm = TRUE),
      .groups = "drop"
    )
  chs_line_data <- chs_data %>%
    group_by(site) %>%
    filter(sum(!is.na(CHS)) >= 2) %>%
    ungroup()

  supp_chs_ts <- ggplot(chs_data, aes(x = year, y = CHS, color = site)) +
    geom_line(data = chs_line_data, linewidth = 0.6, na.rm = TRUE) +
    geom_point(size = 1.5, na.rm = TRUE) +
    geom_text(
      data = chs_means,
      aes(x = x_pos, y = y_pos, label = sprintf("%.2f", mean_val)),
      hjust = 0, vjust = -0.3, size = 3, color = "black", na.rm = TRUE
    ) +
    facet_wrap(~site, ncol = 2, drop = FALSE) +
    scale_color_manual(values = site_colors) +
    labs(x = "Water Year", y = "Baseflow Fraction (CHS)")

  ggsave(file.path(supp_dir, "ms_chs_annual_ts.png"), supp_chs_ts, width = 10, height = 10, dpi = 300)
  ggsave(file.path(supp_dir, "ms_chs_annual_ts.pdf"), supp_chs_ts, width = 10, height = 10)
  cat("  Saved CHS Time Series\n")
}

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------

cat("\n=== FIGURES COMPLETE ===\n")
cat("\nMain Text figures:", main_dir, "\n")
for (f in list.files(main_dir, pattern = "\\.(png|pdf)$")) {
  cat("  -", f, "\n")
}

cat("\nSupplement figures:", supp_dir, "\n")
for (f in list.files(supp_dir, pattern = "\\.(png|pdf)$")) {
  cat("  -", f, "\n")
}
