# -----------------------------------------------------------------------------
# Hydrometric Metric Time Series and Site Summaries
# -----------------------------------------------------------------------------
# Hydrometric storage metric summaries for supplement figures.
#
# Inputs:
#   - master_annual.csv
#
# Outputs (FIGURES_DIR/supp/hydrometric):
#   - <prefix>_<metric>_annual_ts.png
#   - <prefix>_<metric>_ave.png
#   - storage_metric_ave_grid.png
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(multcompView)

rm(list = ls())

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

theme_set(theme_pub())

plot_dir <- file.path(FIGURES_DIR, "supp", "hydrometric")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Remove prior hydrometric summary outputs so this folder reflects current metrics.
old_files <- c(
  list.files(plot_dir, pattern = "_annual_ts\\.png$", full.names = TRUE),
  list.files(plot_dir, pattern = "_ave\\.png$", full.names = TRUE),
  file.path(plot_dir, "storage_metric_ave_grid.png")
)
old_files <- old_files[file.exists(old_files)]
if (length(old_files) > 0) {
  unlink(old_files)
}

annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUTPUT_DIR, MASTER_ANNUAL_FILE)
}
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUT_MASTER_DIR, LEGACY_ANNUAL_FILE)
}
if (!file.exists(annual_file)) {
  stop("master_annual.csv not found.")
}

storage <- read_csv(annual_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

metric_order <- c("FDC", "RBI", "RCS", "SD", "WB", "CHS")
metric_order <- metric_order[metric_order %in% names(storage)]
site_labels_panel <- make_panel_label_map(SITE_ORDER_HYDROMETRIC)

metric_prefix <- c(
  FDC = "ds",
  RBI = "ds",
  RCS = "ds",
  SD = "ds",
  WB = "eds",
  CHS = "ms"
)

axis_labels <- c(
  FDC = "Flow Duration Curve Slope",
  RBI = "Richards-Baker Index",
  RCS = "Recession Curve Slope",
  SD = "Storage-Discharge",
  WB = "Water Balance Drawdown",
  CHS = "Chemical Hydrograph Separation"
)

metric_labels_grid <- c(
  FDC = "A) Flow Duration Curve Slope",
  RBI = "B) Richards-Baker Index",
  RCS = "C) Recession Curve Slope",
  SD = "D) Storage-Discharge",
  WB = "E) Water Balance Drawdown",
  CHS = "F) Chemical Hydrograph Separation"
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

  df_line <- df %>%
    group_by(site) %>%
    filter(sum(!is.na(value)) >= 2) %>%
    ungroup()

  p_ts <- ggplot(df, aes(x = year, y = value, color = site, group = site)) +
    geom_line(data = df_line, linewidth = 0.5) +
    geom_point(size = FIG_POINT_SIZE_SMALL) +
    facet_wrap(
      ~site,
      ncol = 2,
      scales = "free_y",
      drop = FALSE,
      axes = "margins",
      axis.labels = "margins",
      labeller = labeller(site = site_labels_panel)
    ) +
    scale_color_manual(values = site_cols, guide = "none") +
    labs(x = "Water Year", y = axis_labels[[m]]) +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0, size = FIG_STRIP_TEXT_SIZE)
    )

  ggsave(
    file.path(plot_dir, paste0(metric_prefix[[m]], "_", tolower(m), "_annual_ts.png")),
    p_ts,
    width = 8 * FIG_WIDTH_SCALE,
    height = 9 * FIG_HEIGHT_SCALE,
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

  summary_all[[m]] <- sum_df
}

summary_sel <- bind_rows(summary_all)
if (nrow(summary_sel) > 0) {
  summary_sel <- summary_sel %>%
    filter(is.finite(label_y), !is.na(group))

  summary_sel <- summary_sel %>%
    mutate(
      metric = factor(metric, levels = metric_order),
      metric_label = metric_labels_grid[as.character(metric)]
    )
  metric_grid_panel <- make_panel_label_map(unique(as.character(summary_sel$metric_label)))

  p_grid <- ggplot(summary_sel, aes(x = site, y = mean_val, color = site)) +
    geom_point(size = FIG_POINT_SIZE_MED) +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), width = 0.2) +
    geom_text(aes(label = group, y = label_y), size = FIG_ANNOT_TEXT_SIZE, vjust = 0, check_overlap = FIG_LABEL_CHECK_OVERLAP) +
    facet_wrap(
      ~metric_label,
      ncol = 2,
      scales = "free_y",
      axes = "margins",
      axis.labels = "margins",
      labeller = labeller(metric_label = metric_grid_panel)
    ) +
    scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC) +
    scale_color_manual(values = site_cols, guide = "none") +
    labs(x = NULL, y = "Value") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0, size = FIG_STRIP_TEXT_SIZE),
      plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT)
    ) +
    coord_cartesian(clip = FIG_LABEL_CLIP)

  ggsave(
    file.path(plot_dir, "storage_metric_ave_grid.png"),
    p_grid,
    width = 8 * FIG_WIDTH_SCALE,
    height = 10 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
}
