# =============================================================================
# Recession Curve Analysis Plots
# =============================================================================
# Purpose: Create recession curve diagnostic plots showing -dQ/dt vs Q
#          relationships on log-log scale
#
# These plots are used to:
#   1. Visualize the linear reservoir assumption
#   2. Compare recession behavior across watersheds
#   3. Validate recession curve slope calculations
#
# Inputs:
#   - Daily discharge data (from HF00402_v14.csv or water balance file)
#
# Outputs:
#   - recession_curves_all_sites.png: Combined log-log recession plots
#   - recession_curves_<site>.png: Individual site recession curves
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(colorspace)
library(patchwork)

theme_set(theme_classic(base_size = 12))

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
output_dir <- file.path(FIGURES_DIR, "Recession")

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Site order and colors
site_order <- SITE_ORDER_HYDROMETRIC
site_cols <- SITE_COLORS

# =============================================================================
# LOAD DISCHARGE DATA
# =============================================================================

cat("Loading discharge data...\n")

# Try water balance file first (has normalized Q)
wb_file <- file.path(base_dir, "DynamicStorage", "daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv")
q_file <- file.path(base_dir, "Q", "HF00402_v14.csv")

if (file.exists(wb_file)) {
  q_data <- read_csv(wb_file, show_col_types = FALSE) %>%
    mutate(
      DATE = as.Date(DATE, tryFormats = c("%Y-%m-%d", "%m/%d/%Y")),
      Q = Q_mm_d
    ) %>%
    select(DATE, SITECODE, Q) %>%
    rename(site = SITECODE, date = DATE)
  cat("  Using water balance file\n")
} else if (file.exists(q_file)) {
  # Load discharge and drainage areas
  da_file <- file.path(base_dir, "Q", "drainage_area.csv")
  da <- read_csv(da_file, show_col_types = FALSE)

  q_data <- read_csv(q_file, show_col_types = FALSE) %>%
    mutate(DATE = as.Date(DATE, format = "%m/%d/%Y")) %>%
    left_join(da, by = "SITECODE") %>%
    mutate(
      Q = MEAN_Q * 0.0283168 * 86400 / DA_M2  # Convert to mm/day
    ) %>%
    select(date = DATE, site = SITECODE, Q)
  cat("  Using raw discharge file\n")
} else {
  stop("No discharge data found")
}

# Standardize site codes
q_data <- q_data %>%
  mutate(
    site = case_when(
      site == "GSWSMC" ~ "GSMACK",
      site == "GSLOOK_FULL" ~ "GSLOOK",
      TRUE ~ site
    )
  ) %>%
  filter(site %in% site_order, !is.na(Q), Q > 0)

cat("  Loaded", nrow(q_data), "records for", length(unique(q_data$site)), "sites\n")

# =============================================================================
# CALCULATE RECESSION SEGMENTS
# =============================================================================

cat("Identifying recession segments...\n")

# Parameters
min_recession_days <- 3  # Minimum consecutive days of recession
max_gap <- 1             # Allow 1-day gaps in recession

recession_data <- q_data %>%
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(
    # Calculate -dQ/dt (rate of change)
    Q_lag = lag(Q),
    Q_lead = lead(Q),
    dQ_dt = -(Q - Q_lag),  # Negative because recession is decreasing
    Q_mid = (Q + Q_lag) / 2,

    # Identify recession days (Q decreasing)
    is_recession = Q < Q_lag & !is.na(Q_lag),

    # Group consecutive recession days
    recession_break = !is_recession | is.na(is_recession),
    recession_group = cumsum(recession_break)
  ) %>%
  filter(is_recession, dQ_dt > 0, Q_mid > 0) %>%  # Keep only valid recession points
  group_by(site, recession_group) %>%
  mutate(recession_length = n()) %>%
  filter(recession_length >= min_recession_days) %>%  # Filter short recessions
  ungroup()

cat("  Found", nrow(recession_data), "recession points\n")

# =============================================================================
# FIT RECESSION SLOPES (LOG-LOG)
# =============================================================================

cat("Fitting recession slopes...\n")

# Fit log-log linear models by site
recession_fits <- recession_data %>%
  group_by(site) %>%
  summarise(
    n_points = n(),
    # Fit: log(-dQ/dt) = a + b * log(Q)
    slope = coef(lm(log10(dQ_dt) ~ log10(Q_mid)))[2],
    intercept = coef(lm(log10(dQ_dt) ~ log10(Q_mid)))[1],
    r_squared = summary(lm(log10(dQ_dt) ~ log10(Q_mid)))$r.squared,
    .groups = "drop"
  ) %>%
  mutate(site = factor(site, levels = site_order))

cat("\nRecession Slope Fits:\n")
print(recession_fits, n = 12)

# =============================================================================
# FIGURE: COMBINED RECESSION CURVES
# =============================================================================

cat("\nCreating recession curve plots...\n")

# Add site factor for plotting
recession_data <- recession_data %>%
  mutate(site = factor(site, levels = site_order))

# Combined plot (all sites, faceted)
p_all <- ggplot(recession_data, aes(x = Q_mid, y = dQ_dt, color = site)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  facet_wrap(~site, ncol = 2, scales = "fixed") +
  scale_x_log10(
    labels = scales::label_number(accuracy = 0.01),
    limits = c(0.01, NA)
  ) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.001),
    limits = c(0.001, NA)
  ) +
  scale_color_manual(values = site_cols, guide = "none") +
  labs(
    title = "Recession Curves: -dQ/dt vs Q (Log-Log)",
    x = expression(Q~(mm~d^{-1})),
    y = expression(-dQ/dt~(mm~d^{-2}))
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0)
  )

ggsave(
  file.path(output_dir, "recession_curves_all_sites.png"),
  p_all, width = 10, height = 12, dpi = 300
)

# =============================================================================
# FIGURE: COMBINED ON SINGLE PANEL
# =============================================================================

# All sites on one plot
p_combined <- ggplot(recession_data, aes(x = Q_mid, y = dQ_dt, color = site)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10(
    labels = scales::label_number(accuracy = 0.01),
    limits = c(0.01, NA)
  ) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.001),
    limits = c(0.001, NA)
  ) +
  scale_color_manual(values = site_cols, name = "Site") +
  labs(
    title = "Recession Curves by Site",
    subtitle = "Slope of log-log relationship indicates storage behavior",
    x = expression(Q~(mm~d^{-1})),
    y = expression(-dQ/dt~(mm~d^{-2}))
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    legend.position = "right"
  )

ggsave(
  file.path(output_dir, "recession_curves_combined.png"),
  p_combined, width = 10, height = 7, dpi = 300
)

# =============================================================================
# FIGURE: RECESSION SLOPE COMPARISON
# =============================================================================

p_slope <- ggplot(recession_fits, aes(x = site, y = slope, fill = site)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = site_cols, guide = "none") +
  labs(
    title = "Recession Curve Slopes by Site",
    subtitle = "b = 1: linear reservoir; b = 2: quadratic reservoir",
    x = NULL,
    y = "Recession Slope (b)"
  ) +
  annotate("text", x = 0.5, y = 1.05, label = "b = 1", hjust = 0, size = 3, color = "gray40") +
  annotate("text", x = 0.5, y = 2.05, label = "b = 2", hjust = 0, size = 3, color = "gray40") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank()
  )

ggsave(
  file.path(output_dir, "recession_slopes_comparison.png"),
  p_slope, width = 8, height = 5, dpi = 300
)

# =============================================================================
# SAVE RECESSION FIT DATA
# =============================================================================

write.csv(
  recession_fits,
  file.path(output_dir, "recession_curve_fits.csv"),
  row.names = FALSE
)

cat("\n=== RECESSION PLOTS COMPLETE ===\n")
cat("Output directory:", output_dir, "\n")
cat("Files created:\n")
list.files(output_dir) %>%
  paste0("  - ", .) %>%
  cat(sep = "\n")
