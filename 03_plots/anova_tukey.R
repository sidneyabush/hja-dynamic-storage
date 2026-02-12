# -----------------------------------------------------------------------------
# ANOVA/Tukey Plots
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)

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
if (!file.exists(config_path)) config_path <- file.path(dirname(script_dir), "config.R")
if (!file.exists(config_path)) config_path <- file.path(getwd(), "config.R")
source(config_path)

main_dir <- file.path(FIGURES_DIR, "main")
supp_dir <- file.path(FIGURES_DIR, "supp", "stats", "anova_tukey")
table_dir <- file.path(OUT_TABLES_DIR, "anova_tukey")
for (d in c(main_dir, supp_dir, table_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(BASE_DATA_DIR, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
}
annual <- read_csv(annual_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
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
metrics_all <- c("RBI", "RCS", "FDC", "SD", "CHS", "WB")

build_panel <- function(metric_name) {
  d <- annual %>%
    select(site, value = all_of(metric_name)) %>%
    mutate(metric = metric_name)

  y_max <- suppressWarnings(max(d$value, na.rm = TRUE))
  y_min <- suppressWarnings(min(d$value, na.rm = TRUE))
  y_span <- ifelse(is.finite(y_max - y_min) && (y_max - y_min) > 0, y_max - y_min, 1)

  letters_m <- letters_df %>%
    filter(metric == metric_name) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
    left_join(
      d %>%
        group_by(site) %>%
        summarise(y_site_max = suppressWarnings(max(value, na.rm = TRUE)), .groups = "drop"),
      by = "site"
    ) %>%
    mutate(y_label = ifelse(is.finite(y_site_max), y_site_max + 0.04 * y_span, NA_real_))

  ggplot(d, aes(x = site, y = value, color = site, fill = site)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2, na.rm = TRUE) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = FIG_POINT_SIZE_SMALL, alpha = 0.7, na.rm = TRUE) +
    geom_text(
      data = letters_m,
      aes(x = site, y = y_label, label = group_letter),
      inherit.aes = FALSE,
      color = "black",
      size = FIG_ANNOT_TEXT_SIZE,
      na.rm = TRUE
    ) +
    scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC, drop = FALSE) +
    scale_color_manual(values = SITE_COLORS) +
    scale_fill_manual(values = SITE_COLORS) +
    labs(x = NULL, y = "Value", title = metric_name) +
    theme_classic(base_size = FIG_BASE_SIZE) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      plot.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.line = element_blank()
    )
}

main_panels <- lapply(metrics_main, build_panel)
fig_main <- wrap_plots(main_panels, ncol = 2)
ggsave(file.path(main_dir, "ds_anova_tukey.png"), fig_main, width = 11, height = 8, dpi = 300)
ggsave(file.path(main_dir, "ds_anova_tukey.pdf"), fig_main, width = 11, height = 8)

all_panels <- lapply(metrics_all, build_panel)
fig_all <- wrap_plots(all_panels, ncol = 2)
ggsave(file.path(supp_dir, "storage_anova_tukey_all.png"), fig_all, width = 11, height = 12, dpi = 300)
ggsave(file.path(supp_dir, "storage_anova_tukey_all.pdf"), fig_all, width = 11, height = 12)

# Flat table for manuscript/supplement use
if (nrow(anova_df) > 0 && nrow(letters_df) > 0) {
  summary_table <- letters_df %>%
    left_join(anova_df %>% select(metric, F_statistic, p_value), by = "metric") %>%
    arrange(metric, site)
  write_csv(summary_table, file.path(table_dir, "anova_tukey_group_table.csv"))
}
