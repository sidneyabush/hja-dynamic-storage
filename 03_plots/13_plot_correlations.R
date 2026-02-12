# -----------------------------------------------------------------------------
# Correlation and Scatter-Matrix Plots
# -----------------------------------------------------------------------------
# This script makes correlation heatmaps and all ggpairs matrices used in the
# legacy workflow.
#
# Inputs:
#   - master_site.csv
#   - master_annual.csv
#
# Outputs:
#   - QA_Storage_Metrics_CorrPlot.png
#   - QA_Catchment_Attributes_CorrPlot.png
#   - QA_Storage_Catchment_Combined_CorrPlot.png
#   - Avg_Corr_Plot.pdf
#   - Avg_Corr_Plot_Colored.pdf
#   - Yearly_Corr_Plot_Colored.pdf
#   - Yearly_Corr_Plot_bySite.pdf
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(ggcorrplot)
library(GGally)

theme_set(theme_classic(base_size = 12))

rm(list = ls())

# get script directory (works with source() and Rscript)
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

output_dir <- OUTPUT_DIR
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. load data
# -----------------------------------------------------------------------------

site_file <- file.path(output_dir, MASTER_SITE_FILE)
if (!file.exists(site_file)) {
  site_file <- file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv")
}

HJA_Ave <- read_csv(
  site_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

annual_file <- file.path(output_dir, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(BASE_DATA_DIR, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
}

HJA_Yr <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

site_colors <- SITE_COLORS[names(SITE_COLORS) %in% levels(HJA_Ave$site)]

no_se_smoother <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_smooth(data = data, mapping = mapping, se = FALSE, method = "lm", color = "black", ...) +
    geom_point()
}

# -----------------------------------------------------------------------------
# 2. correlation heatmaps
# -----------------------------------------------------------------------------

storage_vars <- c("RCS_mean", "RBI_mean", "CHS_mean", "FDC_mean", "SD_mean", "WB_mean")
storage_vars <- storage_vars[storage_vars %in% names(HJA_Ave)]

if (length(storage_vars) >= 2) {
  cor_storage <- cor(HJA_Ave[storage_vars], use = "pairwise.complete.obs")

  p_storage <- ggcorrplot(
    cor_storage,
    hc.order = FALSE,
    type = "lower",
    outline.col = "white",
    lab = TRUE
  ) + labs(title = "Storage Metrics Correlation Matrix (Site Averages)")

  ggsave(
    file.path(output_dir, "QA_Storage_Metrics_CorrPlot.png"),
    p_storage,
    width = 10,
    height = 10,
    dpi = 300
  )
}

vars_catchment <- c(
  "Area_km2", "Elevation_mean_m", "Slope_mean", "Aspec_Mean_deg",
  "Harvest", "Landslide_Young", "Landslide_Mod", "Landslide_Old",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
)
vars_catchment <- vars_catchment[vars_catchment %in% names(HJA_Ave)]

if (length(vars_catchment) >= 2) {
  cor_catchment <- cor(HJA_Ave[vars_catchment], use = "pairwise.complete.obs")

  p_catchment <- ggcorrplot(
    cor_catchment,
    hc.order = FALSE,
    type = "lower",
    outline.col = "white",
    lab = TRUE
  ) + labs(title = "Catchment Attributes Correlation Matrix")

  ggsave(
    file.path(output_dir, "QA_Catchment_Attributes_CorrPlot.png"),
    p_catchment,
    width = 10,
    height = 10,
    dpi = 300
  )
}

vars_combined <- c(storage_vars, vars_catchment)
vars_combined <- vars_combined[vars_combined %in% names(HJA_Ave)]

if (length(vars_combined) >= 2) {
  cor_combined <- cor(HJA_Ave[vars_combined], use = "pairwise.complete.obs")

  p_combined <- ggcorrplot(
    cor_combined,
    hc.order = FALSE,
    type = "lower",
    outline.col = "white",
    lab = FALSE
  ) + labs(title = "Storage Metrics and Catchment Attributes Correlation Matrix")

  ggsave(
    file.path(output_dir, "QA_Storage_Catchment_Combined_CorrPlot.png"),
    p_combined,
    width = 14,
    height = 14,
    dpi = 300
  )
}

# -----------------------------------------------------------------------------
# 3. full ggpairs matrices (legacy outputs)
# -----------------------------------------------------------------------------

avg_keep_cols <- c(
  "Area_km2", "Elevation_mean_m", "Slope_mean", "Aspec_Mean_deg",
  "Harvest", "Landslide_Young", "Landslide_Mod", "Landslide_Old",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per",
  "RCS_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean", "CHS_mean",
  "FDC_mean", "SD_mean", "WB_mean", "DR", "temp_at_min_Q_7d_C_mean"
)

avg_keep_cols <- avg_keep_cols[avg_keep_cols %in% names(HJA_Ave)]
avg_met_filt <- HJA_Ave %>% select(site, all_of(avg_keep_cols))

if (length(avg_keep_cols) >= 2) {
  pdf(file.path(output_dir, "Avg_Corr_Plot.pdf"), width = 25, height = 25)
  p_avg <- ggpairs(
    avg_met_filt,
    columns = 2:ncol(avg_met_filt),
    lower = list(continuous = no_se_smoother)
  ) + theme_bw()
  print(p_avg)
  dev.off()

  pdf(file.path(output_dir, "Avg_Corr_Plot_Colored.pdf"), width = 25, height = 25)
  p_avg_col <- ggpairs(
    avg_met_filt,
    columns = 2:ncol(avg_met_filt),
    aes(color = site),
    upper = list(continuous = "blank"),
    lower = list(continuous = no_se_smoother)
  ) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = site_colors)
  print(p_avg_col)
  dev.off()
}

yr_keep_cols <- c("RCS", "RBI", "Q5norm", "CHS", "FDC", "SD", "WB", "temp_at_min_Q_7d_C")
yr_keep_cols <- yr_keep_cols[yr_keep_cols %in% names(HJA_Yr)]
yr_met_filt <- HJA_Yr %>% select(site, all_of(yr_keep_cols))

if (length(yr_keep_cols) >= 2) {
  pdf(file.path(output_dir, "Yearly_Corr_Plot_Colored.pdf"), width = 15, height = 15)
  p_yr_col <- ggpairs(
    yr_met_filt,
    columns = 2:ncol(yr_met_filt),
    aes(color = site),
    upper = list(continuous = "blank"),
    lower = list(continuous = no_se_smoother)
  ) +
    theme_bw() +
    scale_color_manual(values = site_colors)
  print(p_yr_col)
  dev.off()
}

yr_site_cols <- c("RCS", "RBI", "Q5norm", "FDC", "SD", "WB", "temp_at_min_Q_7d_C")
yr_site_cols <- yr_site_cols[yr_site_cols %in% names(HJA_Yr)]
yr_met_filt2 <- HJA_Yr %>% select(site, all_of(yr_site_cols))
site_list <- unique(yr_met_filt2$site)

if (length(yr_site_cols) >= 2 && length(site_list) > 0) {
  pdf(file.path(output_dir, "Yearly_Corr_Plot_bySite.pdf"), width = 15, height = 15)
  for (i in seq_along(site_list)) {
    one_site <- yr_met_filt2 %>%
      filter(site == site_list[i]) %>%
      select(-site)

    valid_cols <- names(one_site)[
      sapply(one_site, function(x) {
        sum(!is.na(x)) > 1 && stats::sd(x, na.rm = TRUE) > 0
      })
    ]

    if (length(valid_cols) >= 2) {
      p_site <- ggpairs(
        one_site %>% select(all_of(valid_cols)),
        lower = list(continuous = no_se_smoother)
      ) +
        theme_bw() +
        ggtitle(as.character(site_list[i]))
      print(p_site)
    }
  }
  dev.off()
}
