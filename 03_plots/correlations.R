# Correlation and Scatter-Matrix Plots.
# Inputs: No direct CSV file reads in this script.
# Author: Sidney Bush
# Date: 2026-02-13

library(dplyr)
library(readr)
library(ggplot2)
library(ggcorrplot)

rm(list = ls())

# get script directory (works with source() and Rscript)
# Load project config
source("config.R")


theme_set(
  theme_pub() +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
)

output_dir <- OUT_MASTER_DIR
plot_dir <- file.path(FIGURES_DIR, "supp", "analysis", "correlations")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Clear prior correlation files so the folder reflects current workflow outputs.
old_corr_files <- list.files(plot_dir, pattern = "corrplot\\.png$", full.names = TRUE)
if (length(old_corr_files) > 0) {
  unlink(old_corr_files)
}

# load data

site_file <- file.path(output_dir, MASTER_SITE_FILE)
if (!file.exists(site_file)) {
  site_file <- file.path(OUTPUT_DIR, MASTER_SITE_FILE)
}

HJA_Ave <- read_csv(
  site_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

if (!("basin_slope" %in% names(HJA_Ave)) && ("Slope_mean" %in% names(HJA_Ave))) {
  HJA_Ave <- HJA_Ave %>% mutate(basin_slope = Slope_mean)
}

# catch_chars MLR family correlation matrix

watershed_predictors <- c(
  "basin_slope", "Harvest", "Landslide_Total", "Landslide_Young",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
)
watershed_predictors <- watershed_predictors[watershed_predictors %in% names(HJA_Ave)]

if (length(watershed_predictors) >= 2) {
  cor_catchment <- cor(HJA_Ave[watershed_predictors], use = "pairwise.complete.obs")
  corr_value_text_size <- FIG_TILE_TEXT_SIZE + 2
  corr_axis_text_size <- FIG_AXIS_TEXT_SIZE + 3

  p_catchment <- ggcorrplot(
    cor_catchment,
    hc.order = FALSE,
    type = "upper",
    outline.col = "white",
    lab = TRUE,
    lab_size = corr_value_text_size,
    tl.cex = corr_axis_text_size / ggplot2::.pt
  ) +
    labs(x = NULL, y = NULL) +
    theme(
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE + 2),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE + 2),
      legend.key.height = grid::unit(10, "mm"),
      legend.key.width = grid::unit(4, "mm"),
      axis.text.x = element_text(size = corr_axis_text_size),
      axis.text.y = element_text(size = corr_axis_text_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  ggsave(
    file.path(plot_dir, "watershed_char_storage_mlr_corr.png"),
    p_catchment,
    width = 9 * FIG_WIDTH_SCALE,
    height = 9 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
}

# storage and ecological-variable MLR correlation matrix

annual_file <- file.path(output_dir, MASTER_ANNUAL_FILE)

HJA_Yr <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

if (!("q5_7d_mm_d" %in% names(HJA_Yr)) && ("min_Q_7d_mm_d" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(q5_7d_mm_d = min_Q_7d_mm_d)
}
if (!("T_7DMax" %in% names(HJA_Yr)) && ("max_temp_7d_C" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(T_7DMax = max_temp_7d_C)
}
if (!("Q_7Q5" %in% names(HJA_Yr)) && ("q5_7d_mm_d" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(Q_7Q5 = q5_7d_mm_d)
}
if (!("temp_during_q5_7d_C" %in% names(HJA_Yr)) && ("temp_during_min_Q_7d_C" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(temp_during_q5_7d_C = temp_during_min_Q_7d_C)
}
if (!("temp_during_q5_7d_C" %in% names(HJA_Yr)) && ("temp_at_min_Q_7d_C" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(temp_during_q5_7d_C = temp_at_min_Q_7d_C)
}
if (!("T_Q7Q5" %in% names(HJA_Yr)) && ("temp_during_q5_7d_C" %in% names(HJA_Yr))) {
  HJA_Yr <- HJA_Yr %>%
    mutate(T_Q7Q5 = temp_during_q5_7d_C)
}

eco_response_vars <- c("T_7DMax", "Q_7Q5", "T_Q7Q5")
eco_response_vars <- eco_response_vars[eco_response_vars %in% names(HJA_Yr)]

storage_predictor_vars <- c("RCS", "RBI", "FDC", "SD", "WB", "CHS", "DR")
storage_predictor_vars <- storage_predictor_vars[storage_predictor_vars %in% names(HJA_Yr)]

eco_corr_vars <- unique(c(eco_response_vars, storage_predictor_vars))
eco_corr_vars <- eco_corr_vars[eco_corr_vars %in% names(HJA_Yr)]

if (length(eco_corr_vars) >= 2) {
  cor_eco <- cor(HJA_Yr[eco_corr_vars], use = "pairwise.complete.obs")
  corr_value_text_size <- FIG_TILE_TEXT_SIZE + 2
  corr_axis_text_size <- FIG_AXIS_TEXT_SIZE + 3

  p_eco <- ggcorrplot(
    cor_eco,
    hc.order = FALSE,
    type = "upper",
    outline.col = "white",
    lab = TRUE,
    lab_size = corr_value_text_size,
    tl.cex = corr_axis_text_size / ggplot2::.pt
  ) +
    labs(x = NULL, y = NULL) +
    theme(
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE + 2),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE + 2),
      legend.key.height = grid::unit(10, "mm"),
      legend.key.width = grid::unit(4, "mm"),
      axis.text.x = element_text(size = corr_axis_text_size),
      axis.text.y = element_text(size = corr_axis_text_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  ggsave(
    file.path(plot_dir, "storage_ecovar_mlr_corr.png"),
    p_eco,
    width = 11 * FIG_WIDTH_SCALE,
    height = 11 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
}
