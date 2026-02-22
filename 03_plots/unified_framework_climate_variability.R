# Plot annual climate-variability diagnostics for unified framework.
# Inputs: OUT_STATS_DIR/unified_framework/unified_framework_annual_anomalies.csv;
#         OUT_STATS_DIR/unified_framework/unified_framework_site_sensitivity_slopes.csv.
# Author: Sidney Bush
# Date: 2026-02-21

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

rm(list = ls())

# Load project config
source("config.R")

find_write_root <- function() {
  roots <- c(OUTPUT_DIR, file.path(REPO_DIR, "outputs"))
  for (r in roots) {
    if (!dir.exists(r)) {
      try(dir.create(r, recursive = TRUE, showWarnings = FALSE), silent = TRUE)
    }
    probe <- file.path(r, ".write_probe")
    ok <- tryCatch(
      {
        file.create(probe)
        file.remove(probe)
        TRUE
      },
      warning = function(w) FALSE,
      error = function(e) FALSE
    )
    if (isTRUE(ok)) return(r)
  }
  stop("No writable output root found for unified-framework climate plots.")
}

find_model_file <- function(rel_path) {
  candidates <- c(
    file.path(OUTPUT_DIR, rel_path),
    file.path(REPO_DIR, "outputs", rel_path)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Missing required file: ", candidates[1])
  }
  hit[1]
}

write_root <- find_write_root()
main_fig_dir <- file.path(write_root, "figs", "main")
if (!dir.exists(main_fig_dir)) {
  dir.create(main_fig_dir, recursive = TRUE, showWarnings = FALSE)
}

anom_file <- find_model_file(file.path("models", "unified_framework", "unified_framework_annual_anomalies.csv"))
slope_file <- find_model_file(file.path("models", "unified_framework", "unified_framework_site_sensitivity_slopes.csv"))
pca_scores_file <- find_model_file(file.path("models", "pca", "pca_scores_pc1_pc2.csv"))
pca_loadings_file <- find_model_file(file.path("models", "pca", "pca_loadings.csv"))
pca_variance_file <- find_model_file(file.path("models", "pca", "pca_variance_explained.csv"))

annual_anom <- read_csv(anom_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site, year)

site_slopes <- read_csv(slope_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

main_pca_scores <- read_csv(pca_scores_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  select(site, year, dynamic_pc1 = PC1, dynamic_pc2 = PC2)

pca_loadings <- read_csv(pca_loadings_file, show_col_types = FALSE) %>%
  transmute(metric = feature, PC1, PC2)

pca_variance <- read_csv(pca_variance_file, show_col_types = FALSE)
pc1_pct <- 100 * pca_variance$Variance_Explained[pca_variance$PC == "PC1"][1]
pc2_pct <- 100 * pca_variance$Variance_Explained[pca_variance$PC == "PC2"][1]

# Figure 1: annual dynamic-storage trajectories in PCA space.
traj_base <- annual_anom %>%
  select(-any_of(c("dynamic_pc1", "dynamic_pc2"))) %>%
  left_join(main_pca_scores, by = c("site", "year")) %>%
  filter(is.finite(dynamic_pc1), is.finite(dynamic_pc2))

site_centroids <- traj_base %>%
  group_by(site) %>%
  summarise(
    dynamic_pc1 = mean(dynamic_pc1, na.rm = TRUE),
    dynamic_pc2 = mean(dynamic_pc2, na.rm = TRUE),
    .groups = "drop"
  )

site_label_offsets <- tibble::tribble(
  ~site,  ~dx,    ~dy,
  "WS10", -0.22,  0.12,
  "WS01",  0.20,  0.10,
  "WS09", -0.04, -0.12
)

site_centroid_labels <- site_centroids %>%
  left_join(site_label_offsets, by = "site") %>%
  mutate(
    dx = dplyr::coalesce(dx, 0),
    dy = dplyr::coalesce(dy, 0.08),
    x_lab = dynamic_pc1 + dx,
    y_lab = dynamic_pc2 + dy,
    site_color = dplyr::coalesce(unname(SITE_COLORS[as.character(site)]), "#333333")
  )

loading_scale <- {
  x_rng <- diff(range(traj_base$dynamic_pc1, na.rm = TRUE))
  y_rng <- diff(range(traj_base$dynamic_pc2, na.rm = TRUE))
  max_x <- max(abs(pca_loadings$PC1), na.rm = TRUE)
  max_y <- max(abs(pca_loadings$PC2), na.rm = TRUE)
  0.2 * min(x_rng / max_x, y_rng / max_y)
}

loading_label_map <- c(
  "RBI" = "RBI",
  "RCS" = "RCS",
  "FDC" = "FDC",
  "SD" = "SD",
  "WB" = "WB"
)

loading_label_offsets <- tibble::tribble(
  ~metric, ~ldx,  ~ldy,
  "SD",     0.14,  0.14,
  "RBI",    0.08, -0.02,
  "WB",     0.08,  0.06
)

loading_vectors <- pca_loadings %>%
  left_join(loading_label_offsets, by = "metric") %>%
  mutate(
    x = 0,
    y = 0,
    xend = PC1 * loading_scale,
    yend = PC2 * loading_scale,
    metric_label = dplyr::recode(metric, !!!loading_label_map, .default = metric),
    ldx = dplyr::coalesce(ldx, 0.03),
    ldy = dplyr::coalesce(ldy, 0.03)
  )

make_traj_plot <- function(data_in, color_var, color_label, loadings_df = NULL) {
  data_plot <- data_in %>%
    mutate(site_color = dplyr::coalesce(unname(SITE_COLORS[as.character(site)]), "#999999"))

  p <- ggplot(data_plot, aes(x = dynamic_pc1, y = dynamic_pc2)) +
    stat_ellipse(
      aes(group = site, fill = I(site_color), color = I(site_color)),
      geom = "polygon",
      type = "norm",
      level = 0.80,
      linewidth = 0.5,
      alpha = 0.10
    ) +
    geom_point(
      aes(color = .data[[color_var]]),
      size = FIG_POINT_SIZE_SMALL + 0.6,
      alpha = 0.9
    ) +
    geom_text(
      data = site_centroid_labels,
      aes(x = x_lab, y = y_lab, label = site, color = I(site_color)),
      inherit.aes = FALSE,
      size = FIG_ANNOT_TEXT_SIZE,
      check_overlap = FALSE
    )

  if (!is.null(loadings_df)) {
    p <- p +
      geom_segment(
        data = loadings_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        color = "grey35",
        linewidth = 0.42,
        alpha = 1,
        arrow = grid::arrow(length = grid::unit(0.13, "cm"), type = "closed")
      ) +
      geom_label(
        data = loadings_df,
        aes(x = xend, y = yend, label = metric_label),
        inherit.aes = FALSE,
        color = "grey35",
        fill = scales::alpha("white", 0.72),
        linewidth = 0.25,
        label.r = grid::unit(0.12, "lines"),
        label.padding = grid::unit(0.12, "lines"),
        size = FIG_ANNOT_TEXT_SIZE - 0.45,
        nudge_x = loadings_df$ldx,
        nudge_y = loadings_df$ldy
      )
  }

  p +
    scale_color_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      name = color_label
    ) +
    labs(
      x = paste0("PC1 (", scales::number(pc1_pct, accuracy = 0.1), "%)"),
      y = paste0("PC2\n(", scales::number(pc2_pct, accuracy = 0.1), "%)")
    ) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE - 1),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
    )
}

# Main trajectory: color by winter-precipitation anomaly.
traj_main <- traj_base %>%
  filter(is.finite(P_NovJan_anom))
p_traj_main <- make_traj_plot(
  traj_main,
  color_var = "P_NovJan_anom",
  color_label = "Annual Pwinter\nanomoly",
  loadings_df = loading_vectors
)

ggsave(
  file.path(main_fig_dir, "unified_framework_dynamic_trajectory_pwinter_anom.png"),
  p_traj_main,
  width = 8.5 * FIG_WIDTH_SCALE,
  height = 6.2 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_fig_dir, "unified_framework_dynamic_trajectory_pwinter_anom.pdf"),
  p_traj_main,
  width = 8.5 * FIG_WIDTH_SCALE,
  height = 6.2 * FIG_HEIGHT_SCALE
)

# Supplementary trajectory: color by WB drawdown anomaly.
traj_supp <- traj_base %>%
  filter(is.finite(WB_drawdown_anom))
p_traj_supp <- make_traj_plot(
  traj_supp,
  color_var = "WB_drawdown_anom",
  color_label = "Annual WB\nanomaly",
  loadings_df = loading_vectors
)

ggsave(
  file.path(main_fig_dir, "unified_framework_dynamic_trajectory_wb_anom.png"),
  p_traj_supp,
  width = 8.5 * FIG_WIDTH_SCALE,
  height = 6.2 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_fig_dir, "unified_framework_dynamic_trajectory_wb_anom.pdf"),
  p_traj_supp,
  width = 8.5 * FIG_WIDTH_SCALE,
  height = 6.2 * FIG_HEIGHT_SCALE
)

# Figure 2: site-level anomaly sensitivity slopes vs unified state index.
slope_plot_base <- site_slopes %>%
  filter(response %in% c("Q_7Q5_anom", "T_Q7Q5_anom", "T_7DMax_anom")) %>%
  filter(is.finite(slope), is.finite(unified_state_index))

response_labels <- c(
  "Q_7Q5_anom" = "Q7Q5 anomaly",
  "T_Q7Q5_anom" = "TQ7Q5 anomaly",
  "T_7DMax_anom" = "T7DMax anomaly"
)

make_slope_plot <- function(data_in, predictor_keep, predictor_labels, predictor_colors) {
  df <- data_in %>%
    filter(predictor %in% predictor_keep) %>%
    mutate(
      response = factor(response, levels = names(response_labels), labels = unname(response_labels)),
      predictor = factor(predictor, levels = names(predictor_labels), labels = unname(predictor_labels))
    ) %>%
    mutate(
      pred_idx = as.integer(predictor),
      x_label = unified_state_index + ifelse(pred_idx == 1L, -0.03, 0.03),
      site_label = as.character(site)
    )

  stat_df <- df %>%
    group_by(response, predictor) %>%
    group_modify(~{
      vv <- is.finite(.x$unified_state_index) & is.finite(.x$slope)
      if (sum(vv) < 3 || stats::sd(.x$unified_state_index[vv]) == 0) {
        return(tibble(r2 = NA_real_, p = NA_real_))
      }
      fit <- tryCatch(stats::lm(slope ~ unified_state_index, data = .x[vv, ]), error = function(e) NULL)
      if (is.null(fit)) {
        return(tibble(r2 = NA_real_, p = NA_real_))
      }
      s <- summary(fit)
      p_val <- NA_real_
      if ("unified_state_index" %in% rownames(s$coefficients)) {
        p_val <- as.numeric(s$coefficients["unified_state_index", "Pr(>|t|)"])
      }
      tibble(r2 = as.numeric(s$r.squared), p = p_val)
    }) %>%
    ungroup() %>%
    left_join(
      df %>%
        group_by(response) %>%
        summarise(
          x_min = min(unified_state_index, na.rm = TRUE),
          x_max = max(unified_state_index, na.rm = TRUE),
          y_min = min(slope, na.rm = TRUE),
          y_max = max(slope, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "response"
    ) %>%
    group_by(response) %>%
    mutate(
      pred_row = row_number(),
      x_span = pmax(x_max - x_min, 1e-9),
      y_span = pmax(y_max - y_min, 1e-9),
      x_pos = x_min + 0.03 * x_span,
      y_pos = y_max - (0.04 + 0.10 * (pred_row - 1)) * y_span,
      predictor_label = dplyr::case_when(
        grepl("WB", as.character(predictor), ignore.case = TRUE) ~ "WB",
        grepl("Pwinter", as.character(predictor), ignore.case = TRUE) ~ "P",
        grepl("Dynamic storage", as.character(predictor), ignore.case = TRUE) ~ "dS",
        TRUE ~ gsub("\n", " ", as.character(predictor), fixed = TRUE)
      ),
      p_label = dplyr::case_when(
        !is.finite(p) ~ "NA",
        p < 0.001 ~ "<.001",
        TRUE ~ sprintf("%.2f", p)
      ),
      label = paste0(
        predictor_label, " R2=",
        ifelse(is.finite(r2), sprintf("%.2f", r2), "NA"),
        " p=",
        p_label
      )
    ) %>%
    ungroup()

  ggplot(df, aes(x = unified_state_index, y = slope, color = predictor)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_point(size = FIG_POINT_SIZE_MED, alpha = 0.9) +
    geom_text(
      aes(x = x_label, label = site_label),
      color = "grey20",
      size = FIG_ANNOT_TEXT_SIZE - 0.8,
      show.legend = FALSE,
      check_overlap = TRUE
    ) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    geom_label(
      data = stat_df,
      aes(x = x_pos, y = y_pos, label = label, color = predictor),
      inherit.aes = FALSE,
      show.legend = FALSE,
      hjust = 0,
      vjust = 1,
      fill = scales::alpha("white", 0.74),
      linewidth = 0,
      label.padding = grid::unit(0.10, "lines"),
      size = FIG_ANNOT_TEXT_SIZE - 0.55
    ) +
    facet_wrap(~response, scales = "free_y") +
    scale_color_manual(values = predictor_colors) +
    labs(
      x = "Unified state index",
      y = "Sensitivity to year-to-year anomaly\n(per-site slope)",
      color = "Predictor"
    ) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      strip.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE - 1),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
    )
}

make_slope_plot_by_site <- function(data_in, predictor_name, predictor_display, y_axis_label) {
  df <- data_in %>%
    filter(predictor == predictor_name) %>%
    mutate(
      response = factor(response, levels = names(response_labels), labels = unname(response_labels)),
      site = factor(site, levels = SITE_ORDER_HYDROMETRIC)
    )

  label_df <- df %>%
    group_by(response) %>%
    mutate(
      y_span = pmax(max(slope, na.rm = TRUE) - min(slope, na.rm = TRUE), 1e-9),
      y_label = slope + 0.035 * y_span
    ) %>%
    ungroup()

  stat_df <- df %>%
    group_by(response) %>%
    group_modify(~{
      vv <- is.finite(.x$unified_state_index) & is.finite(.x$slope)
      if (sum(vv) < 3 || stats::sd(.x$unified_state_index[vv]) == 0) {
        return(tibble(r2 = NA_real_, p = NA_real_))
      }
      fit <- tryCatch(stats::lm(slope ~ unified_state_index, data = .x[vv, ]), error = function(e) NULL)
      if (is.null(fit)) {
        return(tibble(r2 = NA_real_, p = NA_real_))
      }
      s <- summary(fit)
      p_val <- NA_real_
      if ("unified_state_index" %in% rownames(s$coefficients)) {
        p_val <- as.numeric(s$coefficients["unified_state_index", "Pr(>|t|)"])
      }
      tibble(r2 = as.numeric(s$r.squared), p = p_val)
    }) %>%
    ungroup() %>%
    left_join(
      df %>%
        group_by(response) %>%
        summarise(
          x_min = min(unified_state_index, na.rm = TRUE),
          x_max = max(unified_state_index, na.rm = TRUE),
          y_min = min(slope, na.rm = TRUE),
          y_max = max(slope, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "response"
    ) %>%
    mutate(
      x_span = pmax(x_max - x_min, 1e-9),
      y_span = pmax(y_max - y_min, 1e-9),
      x_pos = x_min + 0.03 * x_span,
      y_pos = y_max - 0.04 * y_span,
      p_label = dplyr::case_when(
        !is.finite(p) ~ "NA",
        p < 0.001 ~ "<.001",
        TRUE ~ sprintf("%.2f", p)
      ),
      label = paste0(
        "R2=",
        ifelse(is.finite(r2), sprintf("%.2f", r2), "NA"),
        " p=",
        p_label
      )
    )

  ggplot(df, aes(x = unified_state_index, y = slope, color = site)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_point(
      size = FIG_POINT_SIZE_LARGE,
      alpha = 0.9
    ) +
    geom_text(
      data = label_df,
      aes(y = y_label, label = site),
      size = FIG_ANNOT_TEXT_SIZE - 0.2,
      show.legend = FALSE,
      check_overlap = FIG_LABEL_CHECK_OVERLAP
    ) +
    geom_smooth(
      data = df,
      aes(x = unified_state_index, y = slope, group = 1),
      inherit.aes = FALSE,
      method = "lm",
      se = FALSE,
      color = "grey30",
      linewidth = 0.55
    ) +
    geom_label(
      data = stat_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      show.legend = FALSE,
      hjust = -0.02,
      vjust = 1.02,
      fill = scales::alpha("white", 0.74),
      linewidth = 0,
      label.padding = grid::unit(0.10, "lines"),
      size = FIG_ANNOT_TEXT_SIZE - 0.45
    ) +
    facet_wrap(~response, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
    scale_color_manual(values = SITE_COLORS) +
    labs(
      x = "Unified state index",
      y = y_axis_label,
      color = "Site"
    ) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      strip.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
      legend.position = "none"
    )
}

# Remove legacy combined slope figures (P + WB in same facets).
unlink(file.path(main_fig_dir, c(
  "unified_framework_state_dependence_slopes.png",
  "unified_framework_state_dependence_slopes.pdf",
  "unified_framework_state_dependence_slopes_supp_wb.png",
  "unified_framework_state_dependence_slopes_supp_wb.pdf"
)))

# Supplementary split plots by predictor (site-colored points).
p_slope_supp_wb_site <- make_slope_plot_by_site(
  slope_plot_base,
  predictor_name = "WB_drawdown_anom",
  predictor_display = "WB drawdown anomaly",
  y_axis_label = "Sensitivity to year-to-year anomaly:\nd(response anomaly)/d(WB drawdown anomaly)"
)
ggsave(
  file.path(main_fig_dir, "unified_framework_state_dependence_slopes_supp_wb_by_site.png"),
  p_slope_supp_wb_site,
  width = 10.5 * FIG_WIDTH_SCALE,
  height = 4.8 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_fig_dir, "unified_framework_state_dependence_slopes_supp_wb_by_site.pdf"),
  p_slope_supp_wb_site,
  width = 10.5 * FIG_WIDTH_SCALE,
  height = 4.8 * FIG_HEIGHT_SCALE
)

p_slope_supp_p_site <- make_slope_plot_by_site(
  slope_plot_base,
  predictor_name = "P_NovJan_anom",
  predictor_display = "Pwinter anomaly",
  y_axis_label = "Sensitivity to year-to-year anomaly:\nd(response anomaly)/d(Pwinter anomaly)"
)
ggsave(
  file.path(main_fig_dir, "unified_framework_state_dependence_slopes_supp_pwinter_by_site.png"),
  p_slope_supp_p_site,
  width = 10.5 * FIG_WIDTH_SCALE,
  height = 4.8 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_fig_dir, "unified_framework_state_dependence_slopes_supp_pwinter_by_site.pdf"),
  p_slope_supp_p_site,
  width = 10.5 * FIG_WIDTH_SCALE,
  height = 4.8 * FIG_HEIGHT_SCALE
)

# Figure 3: site-mean ecological response vs unified structural state.
eco_site_means <- annual_anom %>%
  group_by(site) %>%
  summarise(
    unified_state_index = first(unified_state_index),
    Q_7Q5 = mean(Q_7Q5, na.rm = TRUE),
    T_Q7Q5 = mean(T_Q7Q5, na.rm = TRUE),
    T_7DMax = mean(T_7DMax, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(unified_state_index))

eco_long <- eco_site_means %>%
  pivot_longer(
    cols = c(Q_7Q5, T_Q7Q5, T_7DMax),
    names_to = "response",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  mutate(
    response = factor(
      response,
      levels = c("Q_7Q5", "T_Q7Q5", "T_7DMax"),
      labels = c(
        "Mean Q7Q5 (mm d-1)",
        "Mean TQ7Q5 (C)",
        "Mean T7DMax (C)"
      )
      )
    )

eco_stats <- eco_long %>%
  group_by(response) %>%
  group_modify(~{
    fit <- stats::lm(value ~ unified_state_index, data = .x)
    s <- summary(fit)
    tibble(
      r2 = as.numeric(s$r.squared),
      p = as.numeric(s$coefficients["unified_state_index", "Pr(>|t|)"]),
      label = paste0(
        "R2 = ", sprintf("%.2f", as.numeric(s$r.squared)),
        "\np = ", formatC(as.numeric(s$coefficients["unified_state_index", "Pr(>|t|)"]), format = "e", digits = 2)
      )
    )
  }) %>%
  ungroup()

eco_label_pos <- eco_long %>%
  group_by(response) %>%
  mutate(
    y_span = pmax(max(value, na.rm = TRUE) - min(value, na.rm = TRUE), 1e-9),
    y_label = value + 0.035 * y_span
  ) %>%
  ungroup()

p_eco_state <- ggplot(eco_long, aes(x = unified_state_index, y = value, color = site)) +
  geom_point(size = FIG_POINT_SIZE_LARGE, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "grey35", linewidth = 0.45) +
  geom_text(
    data = eco_label_pos,
    aes(y = y_label, label = site),
    size = FIG_ANNOT_TEXT_SIZE - 0.2,
    show.legend = FALSE,
    check_overlap = FIG_LABEL_CHECK_OVERLAP
  ) +
  geom_label(
    data = eco_stats,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.02,
    vjust = 1.02,
    fill = scales::alpha("white", 0.75),
    linewidth = 0,
    label.padding = grid::unit(0.12, "lines"),
    size = FIG_ANNOT_TEXT_SIZE - 0.2
  ) +
  facet_wrap(~response, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = "Unified state index",
    y = "Site-mean annual response\n(WY1997-2020)"
  ) +
  theme_pub() +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    strip.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
    legend.position = "none"
  )

ggsave(
  file.path(main_fig_dir, "unified_framework_eco_site_means_vs_state.png"),
  p_eco_state,
  width = 10.8 * FIG_WIDTH_SCALE,
  height = 4.9 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_fig_dir, "unified_framework_eco_site_means_vs_state.pdf"),
  p_eco_state,
  width = 10.8 * FIG_WIDTH_SCALE,
  height = 4.9 * FIG_HEIGHT_SCALE
)
