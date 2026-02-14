# ANOVA/Tukey Plots.
# Inputs: ISOTOPE_DIR/MTT_FYW.csv; ISOTOPE_DIR/DampingRatios_2025-07-07.csv; OUT_STATS_ANOVA_DIR/anova_results.csv.
# Author: Sidney Bush
# Date: 2026-02-13

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)

rm(list = ls())

# Load project config
source("config.R")


main_dir <- file.path(FIGURES_DIR, "main")
supp_dir <- file.path(FIGURES_DIR, "supp", "analysis", "anova_tukey")
table_dir <- file.path(OUT_TABLES_DIR, "anova_tukey")
for (d in c(main_dir, supp_dir, table_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

safe_ggsave <- function(filename, plot, width, height, dpi = NULL) {
  tryCatch(
    {
      if (is.null(dpi)) {
        ggsave(filename, plot, width = width, height = height)
      } else {
        ggsave(filename, plot, width = width, height = height, dpi = dpi)
      }
      TRUE
    },
    error = function(e) {
      warning(sprintf("Could not write %s (%s)", filename, conditionMessage(e)))
      FALSE
    }
  )
}

annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
annual <- read_csv(annual_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

# Isotope predictors are site-level (not annual) and only used in the all-metrics panel.
iso_mtt_file <- file.path(ISOTOPE_DIR, "MTT_FYW.csv")
iso_dr_file <- file.path(ISOTOPE_DIR, "DampingRatios_2025-07-07.csv")

iso_mtt <- if (file.exists(iso_mtt_file)) {
  read_csv(iso_mtt_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    transmute(
      site,
      mtt1 = MTT1,
      mtt1_se = MTT1_SD,
      mtt2 = MTTM,
      mtt2_se = (MTT2L_SD + MTT2H_SD) / 2,
      fyw = FYWM,
      fyw_se = (FYWL_SD + FYWH_SD) / 2
    )
} else {
  tibble(
    site = character(),
    mtt1 = numeric(), mtt1_se = numeric(),
    mtt2 = numeric(), mtt2_se = numeric(),
    fyw = numeric(), fyw_se = numeric()
  )
}

iso_dr <- if (file.exists(iso_dr_file)) {
  read_csv(iso_dr_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    transmute(site, dr = DR_Overall, dr_se = DR__err)
} else {
  tibble(site = character(), dr = numeric(), dr_se = numeric())
}

iso_site <- full_join(iso_mtt, iso_dr, by = "site") %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

letters_file <- file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv")
letters_df <- if (file.exists(letters_file)) {
  read_csv(letters_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site))
} else {
  tibble(metric = character(), site = character(), group_letter = character())
}

anova_file <- file.path(OUT_STATS_ANOVA_DIR, "anova_results.csv")
anova_df <- if (file.exists(anova_file)) {
  read_csv(anova_file, show_col_types = FALSE)
} else {
  tibble(metric = character(), p_value = numeric())
}

metrics_main <- c("RBI", "RCS", "FDC", "SD")
metrics_all <- c("RBI", "RCS", "FDC", "SD", "CHS", "WB", "MTT1", "MTT2", "DR", "FYW")

build_faceted_plot <- function(df_long, metric_names, ncol_facets) {
  d <- df_long %>%
    filter(metric %in% metric_names) %>%
    mutate(metric = factor(metric, levels = metric_names))
  metric_labels_panel <- make_panel_label_map(metric_names)

  letters_m <- letters_df %>%
    filter(metric %in% metric_names) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
    left_join(
      d %>%
        group_by(metric, site) %>%
        summarise(y_site_max = suppressWarnings(max(value, na.rm = TRUE)), .groups = "drop") %>%
        left_join(
          d %>%
            group_by(metric) %>%
            summarise(
              y_span = {
                y_max <- suppressWarnings(max(value, na.rm = TRUE))
                y_min <- suppressWarnings(min(value, na.rm = TRUE))
                ifelse(is.finite(y_max - y_min) && (y_max - y_min) > 0, y_max - y_min, 1)
              },
              .groups = "drop"
            ),
          by = "metric"
        ),
      by = c("metric", "site")
    ) %>%
    mutate(y_label = ifelse(is.finite(y_site_max), y_site_max + 0.04 * y_span, NA_real_)) %>%
    mutate(metric = factor(metric, levels = metric_names))

  ggplot(d, aes(x = site, y = value, color = site, fill = site)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2, na.rm = TRUE) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = FIG_POINT_SIZE_SMALL, alpha = 0.7, na.rm = TRUE) +
    geom_text(
      data = letters_m,
      aes(x = site, y = y_label, label = group_letter),
      inherit.aes = FALSE,
      color = "black",
      size = FIG_ANNOT_TEXT_SIZE,
      check_overlap = FIG_LABEL_CHECK_OVERLAP,
      na.rm = TRUE
    ) +
    scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC, drop = FALSE) +
    scale_color_manual(values = SITE_COLORS) +
    scale_fill_manual(values = SITE_COLORS) +
    facet_wrap(
      ~metric,
      ncol = ncol_facets,
      scales = "free_y",
      labeller = labeller(metric = metric_labels_panel),
      axes = "margins",
      axis.labels = "margins"
    ) +
    labs(x = NULL, y = "Value") +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      strip.background = element_blank(),
      strip.text = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.line = element_blank(),
      plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT)
    ) +
    coord_cartesian(clip = FIG_LABEL_CLIP)
}

build_mean_se_plot <- function(df_long, metric_names, ncol_facets) {
  d <- df_long %>%
    filter(metric %in% metric_names) %>%
    mutate(metric = factor(metric, levels = metric_names))
  metric_labels_panel <- make_panel_label_map(metric_names)

  letters_m <- letters_df %>%
    mutate(metric = toupper(metric)) %>%
    filter(metric %in% metric_names) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
    left_join(
      d %>%
        group_by(metric, site) %>%
        summarise(y_site_max = suppressWarnings(max(mean_value + se_value, na.rm = TRUE)), .groups = "drop") %>%
        left_join(
          d %>%
            group_by(metric) %>%
            summarise(
              y_span = {
                y_max <- suppressWarnings(max(mean_value + se_value, na.rm = TRUE))
                y_min <- suppressWarnings(min(mean_value - se_value, na.rm = TRUE))
                ifelse(is.finite(y_max - y_min) && (y_max - y_min) > 0, y_max - y_min, 1)
              },
              .groups = "drop"
            ),
          by = "metric"
        ),
      by = c("metric", "site")
    ) %>%
    mutate(y_label = ifelse(is.finite(y_site_max), y_site_max + 0.04 * y_span, NA_real_)) %>%
    mutate(metric = factor(metric, levels = metric_names))

  ggplot(d, aes(x = site, y = mean_value, color = site)) +
    geom_point(size = FIG_POINT_SIZE_MED, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 0.3,
      linewidth = 0.5,
      na.rm = TRUE
    ) +
    geom_text(
      data = letters_m,
      aes(x = site, y = y_label, label = group_letter),
      inherit.aes = FALSE,
      color = "black",
      size = FIG_ANNOT_TEXT_SIZE,
      check_overlap = FIG_LABEL_CHECK_OVERLAP,
      na.rm = TRUE
    ) +
    scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC, drop = FALSE) +
    scale_color_manual(values = SITE_COLORS) +
    facet_wrap(
      ~metric,
      ncol = ncol_facets,
      scales = "free_y",
      labeller = labeller(metric = metric_labels_panel),
      axes = "margins",
      axis.labels = "margins"
    ) +
    labs(x = NULL, y = "Value") +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      strip.background = element_blank(),
      strip.text = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.line = element_blank(),
      plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT)
    ) +
    coord_cartesian(clip = FIG_LABEL_CLIP)
}

annual_long <- annual %>%
  select(site, all_of(metrics_main), CHS, WB) %>%
  pivot_longer(cols = -site, names_to = "metric", values_to = "value")

iso_long <- iso_site %>%
  select(site, mtt1, mtt2, dr, fyw) %>%
  pivot_longer(cols = c(mtt1, mtt2, dr, fyw), names_to = "metric", values_to = "value") %>%
  mutate(metric = toupper(metric))

annual_summary <- annual_long %>%
  group_by(metric, site) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    se_value = sd(value, na.rm = TRUE) / sqrt(sum(is.finite(value))),
    .groups = "drop"
  ) %>%
  mutate(metric = toupper(metric))

iso_se_long <- iso_site %>%
  select(site, mtt1_se, mtt2_se, dr_se, fyw_se) %>%
  pivot_longer(cols = -site, names_to = "metric_se", values_to = "se_value") %>%
  mutate(metric = toupper(gsub("_se$", "", metric_se))) %>%
  select(site, metric, se_value)

iso_summary <- iso_long %>%
  left_join(iso_se_long, by = c("site", "metric")) %>%
  transmute(metric, site, mean_value = value, se_value = se_value)

all_summary <- bind_rows(annual_summary, iso_summary)

fig_main <- build_faceted_plot(annual_long, metrics_main, ncol_facets = 2)
safe_ggsave(file.path(main_dir, "ds_anova_tukey.png"), fig_main, width = 11 * FIG_WIDTH_SCALE, height = 8 * FIG_HEIGHT_SCALE, dpi = 300)
safe_ggsave(file.path(main_dir, "ds_anova_tukey.pdf"), fig_main, width = 11 * FIG_WIDTH_SCALE, height = 8 * FIG_HEIGHT_SCALE)

all_long <- bind_rows(annual_long, iso_long)
fig_all <- build_mean_se_plot(all_summary, metrics_all, ncol_facets = 2)
safe_ggsave(file.path(supp_dir, "storage_anova_tukey_all.png"), fig_all, width = 11 * FIG_WIDTH_SCALE, height = 12 * FIG_HEIGHT_SCALE, dpi = 300)
safe_ggsave(file.path(supp_dir, "storage_anova_tukey_all.pdf"), fig_all, width = 11 * FIG_WIDTH_SCALE, height = 12 * FIG_HEIGHT_SCALE)

# Flat table for manuscript/supplement use
if (nrow(anova_df) > 0 && nrow(letters_df) > 0) {
  summary_table <- letters_df %>%
    left_join(anova_df %>% select(metric, F_statistic, p_value), by = "metric") %>%
    arrange(metric, site)
  write_csv(summary_table, file.path(table_dir, "anova_tukey_group_table.csv"))
}
