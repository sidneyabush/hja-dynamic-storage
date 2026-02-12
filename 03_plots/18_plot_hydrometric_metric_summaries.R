# -----------------------------------------------------------------------------
# Hydrometric Metric Time Series and Site Summaries
# -----------------------------------------------------------------------------
# Legacy-style plotting for hydrometric storage metrics.
#
# Inputs:
#   - HJA_StorageMetrics_Annual_All.csv
#
# Outputs (FIGURES_DIR/Supplement/Hydrometric_Legacy):
#   - legacy_timeseries_<metric>_by_site_wy.png
#   - legacy_summary_<metric>_site_mean_sd.png
#   - legacy_selected_metric_summary_grid.png
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(multcompView)

rm(list = ls())

theme_set(theme_classic(base_size = 12))

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

plot_dir <- file.path(FIGURES_DIR, "Supplement", "Hydrometric_Legacy")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

annual_file <- file.path(BASE_DATA_DIR, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUTPUT_DIR, "HJA_StorageMetrics_Annual_All.csv")
}
if (!file.exists(annual_file)) {
  stop("HJA_StorageMetrics_Annual_All.csv not found.")
}

storage <- read_csv(annual_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

metric_order <- c("FDC", "Q5norm", "RBI", "RCS", "SD", "WB", "CHS")
metric_order <- metric_order[metric_order %in% names(storage)]

axis_labels <- c(
  FDC = "Flow Duration Curve Slope",
  Q5norm = "Normalized Q5",
  RBI = "Richards-Baker Index",
  RCS = "Recession Curve Slope",
  SD = "Storage-Discharge",
  WB = "Water Balance Drawdown",
  CHS = "Chemical Hydrograph Separation"
)

metric_labels_grid <- c(
  FDC = "A) Flow Duration Curve Slope",
  Q5norm = "B) Normalized Q5",
  RBI = "C) Richards-Baker Index",
  RCS = "D) Recession Curve Slope",
  SD = "E) Storage-Discharge",
  WB = "F) Water Balance Drawdown",
  CHS = "G) Chemical Hydrograph Separation"
)

site_cols <- SITE_COLORS

storage_long <- storage %>%
  select(site, year, all_of(metric_order)) %>%
  pivot_longer(cols = -c(site, year), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

get_tukey_letters <- function(df_metric) {
  if (n_distinct(df_metric$site) < 2) {
    return(setNames(rep("", n_distinct(df_metric$site)), unique(as.character(df_metric$site))))
  }

  fit <- tryCatch(aov(value ~ site, data = df_metric), error = function(e) NULL)
  if (is.null(fit)) {
    return(setNames(rep("", n_distinct(df_metric$site)), unique(as.character(df_metric$site))))
  }

  tuk <- tryCatch(TukeyHSD(fit, "site")$site, error = function(e) NULL)
  if (is.null(tuk) || nrow(tuk) == 0) {
    return(setNames(rep("", n_distinct(df_metric$site)), unique(as.character(df_metric$site))))
  }

  pvals <- setNames(tuk[, "p adj"], rownames(tuk))
  multcompLetters(pvals)$Letters
}

summary_all <- list()

for (m in metric_order) {
  df <- storage_long %>% filter(metric == m)
  if (nrow(df) == 0) {
    next
  }

  p_ts <- ggplot(df, aes(x = year, y = value, color = site, group = site)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1) +
    facet_wrap(~site, ncol = 2, scales = "free_y") +
    scale_color_manual(values = site_cols, guide = "none") +
    labs(x = "Water Year", y = axis_labels[[m]]) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0)
    )

  ggsave(
    file.path(plot_dir, paste0("legacy_timeseries_", m, "_by_site_wy.png")),
    p_ts,
    width = 8,
    height = 9,
    dpi = 300
  )

  sum_df <- df %>%
    group_by(site) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val = sd(value, na.rm = TRUE),
      .groups = "drop"
    )

  letters <- get_tukey_letters(df)
  sum_df$group <- letters[as.character(sum_df$site)]

  y_max <- max(sum_df$mean_val + sum_df$sd_val, na.rm = TRUE)
  y_min <- min(sum_df$mean_val - sum_df$sd_val, na.rm = TRUE)
  offset <- ifelse(is.finite(y_max - y_min), 0.05 * (y_max - y_min), 0)
  sum_df <- sum_df %>% mutate(label_y = mean_val + sd_val + offset, metric = m)

  p_sum <- ggplot(sum_df, aes(x = site, y = mean_val, color = site)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), width = 0.2) +
    geom_text(aes(label = group, y = label_y), size = 3, vjust = 0) +
    scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC) +
    scale_color_manual(values = site_cols, guide = "none") +
    labs(x = NULL, y = paste0("Mean +/- 1 SD: ", axis_labels[[m]])) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank()
    )

  ggsave(
    file.path(plot_dir, paste0("legacy_summary_", m, "_site_mean_sd.png")),
    p_sum,
    width = 6,
    height = 4,
    dpi = 300
  )

  summary_all[[m]] <- sum_df
}

summary_sel <- bind_rows(summary_all)
if (nrow(summary_sel) > 0) {
  summary_sel <- summary_sel %>%
    mutate(
      metric = factor(metric, levels = metric_order),
      metric_label = metric_labels_grid[as.character(metric)]
    )

  p_grid <- ggplot(summary_sel, aes(x = site, y = mean_val, color = site)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), width = 0.2) +
    geom_text(aes(label = group, y = label_y), size = 3, vjust = 0) +
    facet_wrap(~metric_label, ncol = 2, scales = "free_y") +
    scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC) +
    scale_color_manual(values = site_cols, guide = "none") +
    labs(x = NULL, y = "Mean +/- 1 SD") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0)
    )

  ggsave(
    file.path(plot_dir, "legacy_selected_metric_summary_grid.png"),
    p_grid,
    width = 8,
    height = 10,
    dpi = 300
  )
}
