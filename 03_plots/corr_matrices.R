# Correlation and Scatter-Matrix Plots.
# Inputs: No direct CSV file reads in this script.
# Author: Sidney Bush
# Date: 2026-02-13

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggcorrplot)
library(RColorBrewer)

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

output_dir <- file.path(OUTPUT_DIR, "master")
supp_plot_dir <- file.path(FIGURES_DIR, "supp")
main_plot_dir <- file.path(FIGURES_DIR, "main")
for (d in c(supp_plot_dir, main_plot_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# Use a ColorBrewer diverging palette (red negative -> blue positive).
CORR_COLORS <- RColorBrewer::brewer.pal(11, "RdBu")[c(2, 6, 10)]

# load data

site_file <- file.path(output_dir, MASTER_SITE_FILE)

HJA_Ave <- read_csv(
  site_file,
  show_col_types = FALSE
) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

if (
  !("basin_slope" %in% names(HJA_Ave)) && ("Slope_mean" %in% names(HJA_Ave))
) {
  HJA_Ave <- HJA_Ave %>% mutate(basin_slope = Slope_mean)
}

# catch_chars MLR family correlation matrix

watershed_predictors <- c(
  "basin_slope",
  "Harvest",
  "Landslide_Total",
  "Landslide_Young",
  "Lava1_per",
  "Lava2_per",
  "Ash_Per",
  "Pyro_per"
)
watershed_predictors <- watershed_predictors[
  watershed_predictors %in% names(HJA_Ave)
]

if (length(watershed_predictors) >= 2) {
  catchment_df <- HJA_Ave[watershed_predictors]
  colnames(catchment_df) <- label_catchment_predictor(colnames(catchment_df))
  cor_catchment <- cor(
    catchment_df,
    use = "pairwise.complete.obs"
  )
  corr_value_text_size <- FIG_TILE_TEXT_SIZE + 2
  corr_axis_text_size <- FIG_AXIS_TEXT_SIZE + 3

  p_catchment <- ggcorrplot(
    cor_catchment,
    hc.order = FALSE,
    type = "upper",
    outline.col = "white",
    colors = CORR_COLORS,
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
    file.path(supp_plot_dir, "watershed_char_storage_mlr_corr.png"),
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

if (
  !("q5_7d_mm_d" %in% names(HJA_Yr)) && ("min_Q_7d_mm_d" %in% names(HJA_Yr))
) {
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
if (
  !("temp_during_q5_7d_C" %in% names(HJA_Yr)) &&
    ("temp_during_min_Q_7d_C" %in% names(HJA_Yr))
) {
  HJA_Yr <- HJA_Yr %>%
    mutate(temp_during_q5_7d_C = temp_during_min_Q_7d_C)
}
if (
  !("temp_during_q5_7d_C" %in% names(HJA_Yr)) &&
    ("temp_at_min_Q_7d_C" %in% names(HJA_Yr))
) {
  HJA_Yr <- HJA_Yr %>%
    mutate(temp_during_q5_7d_C = temp_at_min_Q_7d_C)
}
if (
  !("T_Q7Q5" %in% names(HJA_Yr)) && ("temp_during_q5_7d_C" %in% names(HJA_Yr))
) {
  HJA_Yr <- HJA_Yr %>%
    mutate(T_Q7Q5 = temp_during_q5_7d_C)
}
if (!("P_WetSeason" %in% names(HJA_Yr))) {
  HJA_Yr$P_WetSeason <- NA_real_
}
if ("precip_nov_may_mm" %in% names(HJA_Yr)) {
  HJA_Yr$P_WetSeason <- dplyr::coalesce(HJA_Yr$P_WetSeason, HJA_Yr$precip_nov_may_mm)
}
if ("P_NovJan" %in% names(HJA_Yr)) {
  HJA_Yr$P_WetSeason <- dplyr::coalesce(HJA_Yr$P_WetSeason, HJA_Yr$P_NovJan)
}
if ("precip_nov_jan_mm" %in% names(HJA_Yr)) {
  HJA_Yr$P_WetSeason <- dplyr::coalesce(HJA_Yr$P_WetSeason, HJA_Yr$precip_nov_jan_mm)
}

eco_response_vars <- c("Q_7Q5", "T_Q7Q5", "T_7DMax")
eco_response_vars <- eco_response_vars[eco_response_vars %in% names(HJA_Yr)]

eco_predictor_vars <- c("P_WetSeason", PLOT_ORDER_DYNAMIC_STORAGE, "CHS")
eco_predictor_vars <- eco_predictor_vars[
  eco_predictor_vars %in% names(HJA_Yr)
]

eco_corr_vars <- unique(c(eco_response_vars, eco_predictor_vars))
eco_corr_vars <- eco_corr_vars[eco_corr_vars %in% names(HJA_Yr)]

if (length(eco_corr_vars) >= 2) {
  cor_eco <- cor(HJA_Yr[eco_corr_vars], use = "pairwise.complete.obs")
  corr_name_clean <- function(x) {
    out <- gsub("_", "", x)
    out[x == "P_WetSeason"] <- "Pwetseason"
    out
  }
  corr_axis_labels <- function(x) {
    vals <- as.character(x)
    vals[vals == "Pwetseason"] <- "P[wetseason]"
    parse(text = vals)
  }
  dimnames(cor_eco) <- lapply(dimnames(cor_eco), corr_name_clean)
  corr_value_text_size <- FIG_TILE_TEXT_SIZE + 2
  corr_axis_text_size <- FIG_AXIS_TEXT_SIZE + 3

  p_eco <- ggcorrplot(
    cor_eco,
    hc.order = FALSE,
    type = "upper",
    outline.col = "white",
    colors = CORR_COLORS,
    lab = TRUE,
    lab_size = corr_value_text_size,
    tl.cex = corr_axis_text_size / ggplot2::.pt
  ) +
    scale_x_discrete(labels = corr_axis_labels) +
    scale_y_discrete(labels = corr_axis_labels) +
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
    file.path(supp_plot_dir, "storage_ecovar_mlr_corr.png"),
    p_eco,
    width = 11 * FIG_WIDTH_SCALE,
    height = 11 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
}

# dynamic + extended dynamic storage metrics correlation matrix (annual)

dynamic_extended_metrics <- c("RBI", "RCS", "FDC", "SD", "WB")
dynamic_extended_metrics <- dynamic_extended_metrics[
  dynamic_extended_metrics %in% names(HJA_Yr)
]

if (length(dynamic_extended_metrics) >= 2) {
  cor_dynamic_extended <- cor(
    HJA_Yr[, dynamic_extended_metrics, drop = FALSE],
    use = "pairwise.complete.obs"
  )

  corr_value_text_size <- FIG_TILE_TEXT_SIZE + 2
  corr_axis_text_size <- FIG_AXIS_TEXT_SIZE + 1

  p_dynamic_extended <- ggcorrplot(
    cor_dynamic_extended,
    hc.order = FALSE,
    type = "upper",
    outline.col = "white",
    colors = CORR_COLORS,
    lab = TRUE,
    lab_size = corr_value_text_size,
    tl.cex = corr_axis_text_size / ggplot2::.pt
  ) +
    labs(x = NULL, y = NULL) +
    theme(
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.text.x = element_text(angle = 45, hjust = 1, size = corr_axis_text_size),
      axis.text.y = element_text(size = corr_axis_text_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  ggsave(
    file.path(supp_plot_dir, "dynamic_extended_storage_corr.png"),
    p_dynamic_extended,
    width = 8 * FIG_WIDTH_SCALE,
    height = 8 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
}

# dynamic-vs-mobile storage correlation matrix (site-level)

dynamic_metrics_site <- paste0(PLOT_ORDER_DYNAMIC_STORAGE, "_mean")
dynamic_metrics_site <- dynamic_metrics_site[
  dynamic_metrics_site %in% names(HJA_Ave)
]

mobile_site_map <- PLOT_MOBILE_STORAGE_SITE_COLS[PLOT_ORDER_MOBILE_STORAGE]
mobile_site_map <- mobile_site_map[
  !is.na(mobile_site_map) & (unname(mobile_site_map) %in% names(HJA_Ave))
]
mobile_metrics_site <- unname(mobile_site_map)

dynamic_metric_levels <- gsub("_mean$", "", dynamic_metrics_site)
mobile_metric_levels <- names(mobile_site_map)
mobile_col_to_display <- setNames(mobile_metric_levels, mobile_metrics_site)

if ((length(dynamic_metrics_site) >= 1) && (length(mobile_metrics_site) >= 1)) {
  cross_cor <- cor(
    HJA_Ave[, dynamic_metrics_site, drop = FALSE],
    HJA_Ave[, mobile_metrics_site, drop = FALSE],
    use = "pairwise.complete.obs"
  )

  cross_long <- as.data.frame(as.table(cross_cor)) %>%
    rename(dynamic_metric = Var1, mobile_metric = Var2, corr = Freq) %>%
    mutate(
      dynamic_metric = gsub("_mean$", "", dynamic_metric),
      mobile_metric = dplyr::recode(mobile_metric, !!!mobile_col_to_display),
      dynamic_metric = factor(
        dynamic_metric,
        levels = dynamic_metric_levels
      ),
      mobile_metric = factor(
        mobile_metric,
        levels = mobile_metric_levels
      )
    )

  p_cross <- ggplot(
    cross_long,
    aes(x = mobile_metric, y = dynamic_metric, fill = corr)
  ) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(
      aes(label = sprintf("%.2f", corr)),
      size = FIG_TILE_TEXT_SIZE + 1
    ) +
    scale_fill_gradient2(
      low = CORR_COLORS[1],
      mid = "white",
      high = CORR_COLORS[3],
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson's r"
    ) +
    scale_y_discrete(limits = rev(dynamic_metric_levels)) +
    labs(x = "Mobile", y = "Dynamic") +
    theme_pub() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = FIG_AXIS_TEXT_SIZE + 1
      ),
      axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
      axis.title.x = element_text(size = FIG_AXIS_TITLE_SIZE + 2),
      axis.title.y = element_text(size = FIG_AXIS_TITLE_SIZE + 2),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  ggsave(
    file.path(main_plot_dir, "dynamic_mobile_storage_cross_corr.png"),
    p_cross,
    width = 9 * FIG_WIDTH_SCALE,
    height = 7 * FIG_HEIGHT_SCALE,
    dpi = 300
  )

  dynamic_long <- HJA_Ave %>%
    select(site, all_of(dynamic_metrics_site)) %>%
    pivot_longer(
      cols = -site,
      names_to = "dynamic_metric",
      values_to = "dynamic_value"
    )

  mobile_long <- HJA_Ave %>%
    select(site, all_of(mobile_metrics_site)) %>%
    pivot_longer(
      cols = -site,
      names_to = "mobile_metric",
      values_to = "mobile_value"
    )

  cross_points <- dynamic_long %>%
    inner_join(mobile_long, by = "site", relationship = "many-to-many") %>%
    mutate(
      dynamic_metric = gsub("_mean$", "", dynamic_metric),
      mobile_metric = dplyr::recode(mobile_metric, !!!mobile_col_to_display),
      dynamic_metric = factor(
        dynamic_metric,
        levels = dynamic_metric_levels
      ),
      mobile_metric = factor(
        mobile_metric,
        levels = mobile_metric_levels
      )
    )
  cross_points_stats <- cross_points %>%
    filter(is.finite(dynamic_value), is.finite(mobile_value)) %>%
    group_by(dynamic_metric, mobile_metric) %>%
    summarise(
      n = n(),
      r2 = if (n() >= 2) {
        tryCatch(summary(lm(dynamic_value ~ mobile_value))$r.squared, error = function(e) NA_real_)
      } else {
        NA_real_
      },
      beta1 = if (n() >= 2) {
        fit <- tryCatch(lm(dynamic_value ~ mobile_value), error = function(e) NULL)
        if (is.null(fit)) {
          NA_real_
        } else {
          unname(coef(fit)[["mobile_value"]])
        }
      } else {
        NA_real_
      },
      p_value = if (n() >= 3) {
        fit <- tryCatch(lm(dynamic_value ~ mobile_value), error = function(e) NULL)
        if (is.null(fit)) {
          NA_real_
        } else {
          coef_tbl <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
          if (is.null(coef_tbl) || !("mobile_value" %in% rownames(coef_tbl))) {
            NA_real_
          } else {
            suppressWarnings(as.numeric(coef_tbl["mobile_value", "Pr(>|t|)"]))
          }
        }
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_label = case_when(
        !is.finite(p_value) ~ "NA",
        p_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p_value)
      ),
      anno_label = paste0(
        "R2=", ifelse(is.finite(r2), sprintf("%.2f", r2), "NA"),
        " | beta1=", ifelse(is.finite(beta1), sprintf("%.2f", beta1), "NA"),
        " | p=", p_label
      )
    )

  top_dynamic_metric <- dynamic_metric_levels[1]
  top_n_by_mobile <- cross_points_stats %>%
    filter(dynamic_metric == top_dynamic_metric) %>%
    select(mobile_metric, n)

  top_n_lookup <- setNames(top_n_by_mobile$n, as.character(top_n_by_mobile$mobile_metric))
  mobile_labels_with_n <- setNames(
    ifelse(
      is.na(top_n_lookup[mobile_metric_levels]),
      paste0(mobile_metric_levels, " (n=NA)"),
      paste0(mobile_metric_levels, " (n=", top_n_lookup[mobile_metric_levels], ")")
    ),
    mobile_metric_levels
  )

  p_cross_points <- ggplot(
    cross_points,
    aes(x = mobile_value, y = dynamic_value)
  ) +
    geom_point(
      color = "black",
      alpha = 0.7,
      size = 1.6,
      na.rm = TRUE
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.5,
      color = CORR_COLORS[3],
      na.rm = TRUE
    ) +
    geom_text(
      data = cross_points_stats,
      aes(x = -Inf, y = Inf, label = anno_label),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.1,
      size = FIG_ANNOT_TEXT_SIZE
    ) +
    facet_grid(
      dynamic_metric ~ mobile_metric,
      scales = "free",
      labeller = labeller(mobile_metric = as_labeller(mobile_labels_with_n))
    ) +
    labs(x = "Mobile", y = "Dynamic") +
    theme_pub() +
    theme(
      axis.text.x = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title.x = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
      axis.title.y = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
      strip.text.x = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
      strip.text.y = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  ggsave(
    file.path(supp_plot_dir, "dynamic_mobile_storage_points_facet.png"),
    p_cross_points,
    width = 12 * FIG_WIDTH_SCALE,
    height = 10 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
}

# CHS vs other mobile storage metrics correlation matrix (site-level)

mobile_corr_vars <- mobile_metrics_site

if (length(mobile_corr_vars) >= 2) {
  cor_mobile <- cor(
    HJA_Ave[, mobile_corr_vars, drop = FALSE],
    use = "pairwise.complete.obs"
  )
  colnames(cor_mobile) <- unname(mobile_col_to_display[colnames(cor_mobile)])
  rownames(cor_mobile) <- unname(mobile_col_to_display[rownames(cor_mobile)])

  corr_value_text_size <- FIG_TILE_TEXT_SIZE + 2
  corr_axis_text_size <- FIG_AXIS_TEXT_SIZE + 1

  p_mobile <- ggcorrplot(
    cor_mobile,
    hc.order = FALSE,
    type = "upper",
    outline.col = "white",
    colors = CORR_COLORS,
    lab = TRUE,
    lab_size = corr_value_text_size,
    tl.cex = corr_axis_text_size / ggplot2::.pt
  ) +
    labs(x = NULL, y = NULL) +
    theme(
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.text.x = element_text(angle = 45, hjust = 1, size = corr_axis_text_size),
      axis.text.y = element_text(size = corr_axis_text_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  ggsave(
    file.path(supp_plot_dir, "chs_mobile_storage_corr.png"),
    p_mobile,
    width = 8 * FIG_WIDTH_SCALE,
    height = 8 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
}
