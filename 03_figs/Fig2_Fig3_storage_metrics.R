# Figures 2 and 3 storage metrics

# inputs:
# outputs/master/master_annual.csv
# outputs/master/master_site.csv
# outputs/models/anova_tukey/tukey_group_letters.csv
# outputs/models/pca/pca_scores_pc1_pc2.csv
# outputs/models/pca/pca_loadings.csv
# outputs/models/pca/pca_variance_explained.csv
# outputs/metrics/dynamic/fdc_slopes_wy.csv
# isotope_dir/MTT_FYW.csv
# isotope_dir/DampingRatios_2025-07-07.csv

# outputs:
# figs_tables_pub/main/Fig2_dynamic_storage_pca.*
# figs_tables_pub/main/Fig3_mobile_storage.*

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(dplyr, readr, tidyr, ggplot2, patchwork, cowplot, cran_repo = "https://cloud.r-project.org")

rm(list = ls())
source("config.R")

# set output folders for the main text figures
main_dir <- MS_FIG_MAIN_DIR
main_pdf_dir <- MS_FIG_MAIN_PDF_DIR
main_tiff_dir <- MS_FIG_MAIN_TIFF_DIR
for (d in c(main_dir, main_pdf_dir, main_tiff_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
letters_file <- file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv")
isotope_file <- file.path(ISOTOPE_DIR, "MTT_FYW.csv")
damping_file <- file.path(ISOTOPE_DIR, "DampingRatios_2025-07-07.csv")
fdc_wy_file <- file.path(OUT_MET_DYNAMIC_DIR, "fdc_slopes_wy.csv")
pca_scores_file <- file.path(OUT_STATS_PCA_DIR, "pca_scores_pc1_pc2.csv")
pca_loadings_file <- file.path(OUT_STATS_PCA_DIR, "pca_loadings.csv")
pca_variance_file <- file.path(OUT_STATS_PCA_DIR, "pca_variance_explained.csv")

# load annual values, site values, PCA results, and Tukey letters
annual <- read_csv(annual_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

# use annual site year FDC slopes in Figure 2
fdc_wy <- read_csv(fdc_wy_file, show_col_types = FALSE) %>%
  mutate(
    site = standardize_site_code(as.character(site)),
    year = suppressWarnings(as.numeric(
      if ("WaterYear" %in% names(.)) WaterYear else if ("year" %in% names(.)) year else NA_real_
    )),
    FDC_annual = suppressWarnings(as.numeric(
      if ("Slope" %in% names(.)) Slope else if ("fdc_slope" %in% names(.)) fdc_slope else if ("FDC" %in% names(.)) FDC else NA_real_
    ))
  ) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC, is.finite(year)) %>%
  select(site, year, FDC_annual) %>%
  distinct(site, year, .keep_all = TRUE)

annual <- annual %>%
  mutate(site_chr = as.character(site)) %>%
  left_join(fdc_wy, by = c("site_chr" = "site", "year" = "year")) %>%
  mutate(
    FDC = dplyr::coalesce(FDC_annual, FDC),
    site = factor(site_chr, levels = SITE_ORDER_HYDROMETRIC)
  ) %>%
  select(-site_chr, -FDC_annual)

site_summary <- read_csv(site_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

pca_scores <- read_csv(pca_scores_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

pca_loadings <- read_csv(pca_loadings_file, show_col_types = FALSE)
pca_variance <- read_csv(pca_variance_file, show_col_types = FALSE)

letters_df <- read_csv(letters_file, show_col_types = FALSE) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

# shared display options for the storage metric panels
BOX_FILL_ALPHA <- 0.22
POINT_ALPHA <- 0.55
PANEL_BG <- "white"
PANEL_BORDER <- "grey55"

# functions for uncertainty bars, axis limits, and repeated styling
first_finite <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x[[1]]
}

calc_ylim <- function(values, labels = numeric()) {
  values <- suppressWarnings(as.numeric(values))
  labels <- suppressWarnings(as.numeric(labels))
  all_vals <- c(values[is.finite(values)], labels[is.finite(labels)])
  if (length(all_vals) == 0) return(c(0, 1))
  y_min <- min(all_vals, na.rm = TRUE)
  y_max <- max(all_vals, na.rm = TRUE)
  span <- y_max - y_min
  if (!is.finite(span) || span <= 0) {
    span <- max(abs(y_max), 1) * 0.1
  }
  c(y_min - 0.05 * span, y_max + 0.12 * span)
}

theme_storage_panel <- function() {
  theme_pub() +
    theme(
      panel.background = element_rect(fill = PANEL_BG, color = NA),
      plot.background = element_rect(fill = PANEL_BG, color = NA),
      panel.border = element_rect(color = PANEL_BORDER, fill = NA, linewidth = 0.7),
      axis.line = element_blank(),
      plot.title = element_text(
        size = FIG_STRIP_TEXT_SIZE + 2,
        hjust = 0,
        margin = margin(t = 1, r = 0, b = 9, l = 0)
      ),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT_SIZE),
      axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.position = "none"
    )
}

equalize_plot_widths <- function(plot_list) {
  grobs <- lapply(plot_list, function(p) {
    if (inherits(p, "ggplot")) {
      ggplotGrob(p)
    } else {
      patchwork::patchworkGrob(p)
    }
  })
  max_widths <- Reduce(grid::unit.pmax, lapply(grobs, function(g) g$widths))

  lapply(grobs, function(g) {
    g$widths <- max_widths
    patchwork::wrap_elements(full = g)
  })
}

# Figure 2 annual dynamic storage boxplots
metric_map_fig2 <- c(
  RBI = "RBI",
  RCS = "RCS",
  FDC = "FDC",
  SD = "SD",
  WB = "WB"
)

fig2_long <- annual %>%
  select(site, year, any_of(unname(metric_map_fig2))) %>%
  pivot_longer(cols = all_of(unname(metric_map_fig2)), names_to = "metric", values_to = "value") %>%
  filter(is.finite(value)) %>%
  mutate(metric = factor(metric, levels = names(metric_map_fig2)))

metric_offsets <- fig2_long %>%
  group_by(metric) %>%
  summarise(
    span = diff(range(value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(offset = ifelse(is.finite(span) & span > 0, 0.06 * span, 0.06))

letters_fig2 <- letters_df %>%
  filter(metric %in% names(metric_map_fig2)) %>%
  left_join(
    fig2_long %>% group_by(metric, site) %>% summarise(ymax_site = max(value, na.rm = TRUE), .groups = "drop"),
    by = c("metric", "site")
  ) %>%
  left_join(metric_offsets, by = "metric") %>%
  mutate(y = ymax_site + offset)

fig2_titles <- c(
  RBI = "a) Richards-Baker Index (RBI)",
  RCS = "b) Recession-curve Slope (RCS)",
  FDC = "c) Flow-duration Curve Slope (FDC)",
  SD = "d) Storage-Discharge (SD)",
  WB = "e) Water-balance deficit (WB)"
)
fig2_y_labels <- c(
  RBI = "Unitless",
  RCS = "Unitless",
  FDC = "Unitless",
  SD = "Depth (mm)",
  WB = "Deficit (mm)"
)

build_fig2_panel <- function(metric_key) {
  # boxplots and Tukey letters share the same structure for all five metrics
  dat <- fig2_long %>% filter(metric == metric_key)
  let <- letters_fig2 %>% filter(metric == metric_key)
  y_lim <- calc_ylim(dat$value, let$y)
  show_x <- TRUE
  reserve_x_space <- FALSE
  panel_margin <- if (identical(metric_key, "FDC")) {
    margin(5.5, 5.5, -7, 5.5)
  } else if (identical(metric_key, "WB")) {
    # add top spacing so panel e does not crowd panel c
    margin(11, 5.5, 5.5, 5.5)
  } else {
    margin(5.5, 5.5, 5.5, 5.5)
  }

  ggplot(dat, aes(x = site, y = value, fill = site, color = site)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = BOX_FILL_ALPHA, linewidth = 0.8) +
    geom_jitter(width = 0.13, alpha = POINT_ALPHA, size = FIG_POINT_SIZE_SMALL + 0.2) +
    geom_text(
      data = let,
      aes(x = site, y = y, label = group_letter),
      inherit.aes = FALSE,
      size = FIG_ANNOT_TEXT_SIZE + 0.2,
      color = "grey35",
      family = "sans"
    ) +
    scale_fill_manual(values = SITE_COLORS, drop = FALSE) +
    scale_color_manual(values = SITE_COLORS, drop = FALSE) +
    scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC, drop = FALSE) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    labs(
      x = NULL,
      y = fig2_y_labels[[metric_key]],
      title = fig2_titles[[metric_key]]
    ) +
    theme_storage_panel() +
    theme(
      axis.text.x = if (show_x) {
        element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT_SIZE)
      } else if (reserve_x_space) {
        element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT_SIZE, colour = NA)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) {
        element_line(colour = "grey35", linewidth = 0.7)
      } else if (reserve_x_space) {
        element_line(colour = NA)
      } else {
        element_blank()
      },
      axis.ticks.length.x = grid::unit(4.5, "pt"),
      plot.margin = panel_margin
    )
}

build_fig2_pca_panel <- function() {
  # the PCA panel uses site year scores and metric loading arrows
  pc1_pct <- ifelse(nrow(pca_variance) >= 1, 100 * pca_variance$Variance_Explained[1], NA_real_)
  pc2_pct <- ifelse(nrow(pca_variance) >= 2, 100 * pca_variance$Variance_Explained[2], NA_real_)
  pca_point_alpha <- 0.55
  site_shapes <- setNames(c(21, 22, 23, 24, 25, 15, 16, 17, 18, 19), SITE_ORDER_HYDROMETRIC)
  site_point_sizes <- setNames(rep(FIG_POINT_SIZE_LARGE + 0.7, length(SITE_ORDER_HYDROMETRIC)), SITE_ORDER_HYDROMETRIC)
  diamond_sites <- names(site_shapes)[site_shapes %in% c(18, 23)]
  site_point_sizes[diamond_sites] <- site_point_sizes[diamond_sites] + 0.8

  score_lim <- max(abs(c(pca_scores$PC1, pca_scores$PC2)), na.rm = TRUE)
  load_lim <- max(abs(c(pca_loadings$PC1, pca_loadings$PC2)), na.rm = TRUE)
  arrow_scale <- ifelse(
    is.finite(score_lim) && is.finite(load_lim) && load_lim > 0,
    0.75 * score_lim / load_lim,
    1
  )

  loadings_plot <- pca_loadings %>%
    mutate(
      PC1s = PC1 * arrow_scale,
      PC2s = PC2 * arrow_scale,
      pretty_label = toupper(feature)
    )

  ggplot(pca_scores, aes(x = PC1, y = PC2, color = site, fill = site, shape = site)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_point(aes(size = site), alpha = pca_point_alpha, stroke = 0.8) +
    geom_segment(
      data = loadings_plot,
      aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
      inherit.aes = FALSE,
      arrow = arrow(length = grid::unit(0.16, "cm")),
      color = "black",
      linewidth = 0.5
    ) +
    geom_label(
      data = loadings_plot,
      aes(x = PC1s, y = PC2s, label = pretty_label),
      inherit.aes = FALSE,
      color = "black",
      fill = scales::alpha("white", 0.78),
      linewidth = 0,
      label.padding = grid::unit(0.10, "lines"),
      size = FIG_ANNOT_TEXT_SIZE + 0.2,
      nudge_x = 0.05,
      nudge_y = 0.05
    ) +
    scale_color_manual(values = SITE_COLORS, name = NULL, drop = FALSE) +
    scale_fill_manual(values = SITE_COLORS, guide = "none", drop = FALSE) +
    scale_shape_manual(values = site_shapes, name = NULL, drop = FALSE) +
    scale_size_manual(values = site_point_sizes, guide = "none", drop = FALSE) +
    guides(
      color = guide_legend(
        ncol = 1,
        override.aes = list(
          shape = unname(site_shapes[SITE_ORDER_HYDROMETRIC]),
          fill = unname(SITE_COLORS[SITE_ORDER_HYDROMETRIC]),
          alpha = pca_point_alpha,
          size = unname(site_point_sizes[SITE_ORDER_HYDROMETRIC]),
          stroke = 0.8
        )
      ),
      shape = "none",
      size = "none"
    ) +
    labs(
      x = paste0("PC1 (", scales::number(pc1_pct, accuracy = 0.1), "%)"),
      y = paste0("PC2 (", scales::number(pc2_pct, accuracy = 0.1), "%)"),
      title = "f) Principal Components Analysis (PCA)"
    ) +
    theme_storage_panel() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = FIG_AXIS_TEXT_SIZE),
      axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.position = "right",
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE - 2),
      legend.key.height = grid::unit(9, "pt"),
      legend.key.width = grid::unit(9, "pt"),
      legend.spacing.y = grid::unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(5.5, 3, 5.5, 5.5)
    ) +
    coord_cartesian(clip = FIG_LABEL_CLIP)
}

p_fig2_a <- build_fig2_panel("RBI")
p_fig2_b <- build_fig2_panel("RCS")
p_fig2_c <- build_fig2_panel("FDC")
p_fig2_d <- build_fig2_panel("SD")
p_fig2_e <- build_fig2_panel("WB")
p_fig2_f <- build_fig2_pca_panel()

# build the columns separately so panel d labels do not push panel e down
p_fig2_left <- p_fig2_a / p_fig2_c / p_fig2_e +
  plot_layout(heights = c(1, 1, 1))
p_fig2_right <- p_fig2_b / p_fig2_d / p_fig2_f +
  plot_layout(heights = c(1, 1, 1))

p_fig2 <- (p_fig2_left | p_fig2_right) +
  plot_layout(widths = c(1, 1))

# baseflow fraction data for Figure 3
fig4_df <- annual %>%
  select(site, year, BF) %>%
  filter(is.finite(BF))

letters_fig4 <- letters_df %>%
  filter(metric == "BF") %>%
  left_join(
    fig4_df %>% group_by(site) %>% summarise(ymax_site = max(BF, na.rm = TRUE), .groups = "drop"),
    by = "site"
  ) %>%
  mutate(y = ymax_site + 0.04)

chs_counts <- tibble(site = factor(SITE_ORDER_HYDROMETRIC, levels = SITE_ORDER_HYDROMETRIC)) %>%
  left_join(
    annual %>%
      group_by(site) %>%
      summarise(n = sum(is.finite(BF)), .groups = "drop"),
    by = "site"
  ) %>%
  mutate(n = ifelse(is.na(n), 0L, as.integer(n)))

fig4_x_labels <- setNames(as.character(chs_counts$site), as.character(chs_counts$site))

# Figure 3 mobile storage panels
fig4_mobile_map <- c(
  DR = "DR",
  Fyw = "Fyw",
  MTT = "MTT"
)

# use DR, Fyw, and MTT uncertainty columns from the isotope source tables
isotope_err <- tibble(site = SITE_ORDER_HYDROMETRIC)

mtt_err <- read_csv(isotope_file, show_col_types = FALSE) %>%
  mutate(
    site = standardize_site_code(as.character(site)),
    Fyw_err = rowMeans(
      cbind(
        suppressWarnings(as.numeric(FYWL_SD)),
        suppressWarnings(as.numeric(FYWH_SD))
      ),
      na.rm = TRUE
    ),
    MTT_early_window_err = suppressWarnings(as.numeric(MTT1_SD)),
    MTT_late_window_err = rowMeans(
      cbind(
        suppressWarnings(as.numeric(MTT2L_SD)),
        suppressWarnings(as.numeric(MTT2H_SD))
      ),
      na.rm = TRUE
    ),
    MTT_err = rowMeans(cbind(MTT_early_window_err, MTT_late_window_err), na.rm = TRUE),
    Fyw_err = ifelse(is.nan(Fyw_err), NA_real_, Fyw_err),
    MTT_late_window_err = ifelse(is.nan(MTT_late_window_err), NA_real_, MTT_late_window_err),
    MTT_err = ifelse(is.nan(MTT_err), NA_real_, MTT_err)
  ) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  group_by(site) %>%
  summarise(
    Fyw_err = first_finite(Fyw_err),
    MTT_err = first_finite(MTT_err),
    .groups = "drop"
  )

isotope_err <- isotope_err %>% left_join(mtt_err, by = "site")

dr_err <- read_csv(damping_file, show_col_types = FALSE) %>%
  mutate(
    site = standardize_site_code(as.character(site)),
    DR_err = suppressWarnings(as.numeric(DR__err))
  ) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  group_by(site) %>%
  summarise(DR_err = first_finite(DR_err), .groups = "drop")

isotope_err <- isotope_err %>% left_join(dr_err, by = "site")

fig5_err_long <- isotope_err %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  pivot_longer(
    cols = c(DR_err, Fyw_err, MTT_err),
    names_to = "metric",
    values_to = "err"
  ) %>%
  mutate(metric = sub("_err$", "", metric))

fig4_mobile_df <- site_summary %>%
  select(site, any_of(unname(fig4_mobile_map))) %>%
  pivot_longer(cols = all_of(unname(fig4_mobile_map)), names_to = "metric", values_to = "value") %>%
  left_join(fig5_err_long, by = c("site", "metric")) %>%
  filter(is.finite(value)) %>%
  mutate(metric = factor(metric, levels = names(fig4_mobile_map)))

fig4_titles <- list(
  BF = expression("a) Baseflow Fraction (BF)"),
  DR = expression("b) Damping Ratio (DR)"),
  Fyw = expression("c) Young Water Fraction (" * F[yw] * ")"),
  MTT = expression("d) Mean Transit Time (MTT)")
)
fig4_y_labels <- c(
  BF = "Unitless",
  DR = "Unitless",
  Fyw = "Unitless",
  MTT = "Years"
)

build_fig4_panel <- function(metric_key) {
  # BF is annual and plotted as boxplots, isotope metrics are site values with uncertainty
  if (identical(metric_key, "BF")) {
    return(
      ggplot(fig4_df, aes(x = site, y = BF, fill = site, color = site)) +
        geom_boxplot(width = 0.65, outlier.shape = NA, alpha = BOX_FILL_ALPHA, linewidth = 0.8) +
        geom_jitter(width = 0.13, alpha = POINT_ALPHA, size = FIG_POINT_SIZE_SMALL + 0.2) +
        geom_text(
          data = letters_fig4,
          aes(x = site, y = y, label = group_letter),
          inherit.aes = FALSE,
          size = FIG_ANNOT_TEXT_SIZE + 0.2,
          color = "grey35"
        ) +
        scale_fill_manual(values = SITE_COLORS, drop = FALSE) +
        scale_color_manual(values = SITE_COLORS, drop = FALSE) +
        scale_x_discrete(
          limits = SITE_ORDER_HYDROMETRIC,
          labels = fig4_x_labels,
          drop = FALSE
        ) +
        coord_cartesian(ylim = calc_ylim(fig4_df$BF, letters_fig4$y), clip = "off") +
        labs(x = NULL, y = fig4_y_labels[[metric_key]]) +
        ggtitle(fig4_titles[[metric_key]]) +
        theme_storage_panel() +
        theme(
          axis.text.x = element_text(
            angle = 45, hjust = 1, vjust = 1,
            size = FIG_AXIS_TEXT_SIZE
          ),
          axis.ticks.x = element_line(colour = "grey35", linewidth = 0.7),
          axis.ticks.length.x = grid::unit(4.5, "pt"),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5)
        )
    )
  }

  dat <- fig4_mobile_df %>% filter(metric == metric_key)
  # keep x axis labels only on the bottom row
  show_x <- TRUE
  reserve_x_space <- metric_key %in% c("DR")
  panel_margin <- margin(5.5, 5.5, 5.5, 5.5)
  y_lim <- calc_ylim(
    c(dat$value, dat$value - dat$err, dat$value + dat$err),
    numeric()
  )

  ggplot(dat, aes(x = site, y = value, color = site)) +
    geom_errorbar(
      aes(ymin = value - err, ymax = value + err),
      width = 0.14,
      alpha = 0.85,
      linewidth = 0.75,
      na.rm = TRUE
    ) +
    geom_point(size = FIG_POINT_SIZE_LARGE + 0.6, alpha = 0.95) +
    scale_color_manual(values = SITE_COLORS, drop = FALSE) +
    scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC, drop = FALSE) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    labs(
      x = NULL,
      y = fig4_y_labels[[metric_key]]
    ) +
    ggtitle(fig4_titles[[metric_key]]) +
    theme_storage_panel() +
    theme(
      axis.title.y = if (identical(metric_key, "MTT")) {
        element_text(size = FIG_AXIS_TITLE_SIZE, margin = margin(r = 22))
      } else {
        element_text(size = FIG_AXIS_TITLE_SIZE)
      },
      axis.text.x = if (show_x) {
        element_text(angle = 45, hjust = 1, vjust = 1, size = FIG_AXIS_TEXT_SIZE)
      } else if (reserve_x_space) {
        element_text(angle = 45, hjust = 1, vjust = 1, size = FIG_AXIS_TEXT_SIZE, colour = NA)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) {
        element_line(colour = "grey35", linewidth = 0.7)
      } else if (reserve_x_space) {
        element_line(colour = NA)
      } else {
        element_blank()
      },
      axis.ticks.length.x = grid::unit(4.5, "pt"),
      plot.margin = panel_margin
    )
}

fig4_order <- c("BF", names(fig4_mobile_map))
fig4_raw_plots <- lapply(fig4_order, build_fig4_panel)

left_aligned <- equalize_plot_widths(list(
  p_fig2_a, p_fig2_c, p_fig2_e
))
p_fig2_a <- left_aligned[[1]]
p_fig2_c <- left_aligned[[2]]
p_fig2_e <- left_aligned[[3]]

right_aligned <- equalize_plot_widths(list(
  p_fig2_b, p_fig2_d
))
p_fig2_b <- right_aligned[[1]]
p_fig2_d <- right_aligned[[2]]

p_fig4 <- patchwork::wrap_plots(
  A = fig4_raw_plots[[1]],
  B = fig4_raw_plots[[2]],
  C = fig4_raw_plots[[3]],
  D = fig4_raw_plots[[4]],
  design = "
  AB
  CD
  "
)

save_figure_set(
  p_fig2,
  "Fig2_dynamic_storage_pca",
  width = 10.4 * FIG_WIDTH_SCALE,
  height = 12.0 * FIG_HEIGHT_SCALE,
  png_dir = main_dir,
  pdf_dir = main_pdf_dir,
  tiff_dir = main_tiff_dir
)

save_figure_set(
  p_fig4,
  "Fig3_mobile_storage",
  width = 10.4 * FIG_WIDTH_SCALE,
  height = 8.0 * FIG_HEIGHT_SCALE,
  png_dir = main_dir,
  pdf_dir = main_pdf_dir,
  tiff_dir = main_tiff_dir
)
