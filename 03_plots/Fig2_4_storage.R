# figures 2 4 5 main manuscript plots

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

rm(list = ls())
source("config.R")

main_dir <- MS_FIG_MAIN_DIR
main_pdf_dir <- MS_FIG_MAIN_PDF_DIR
for (d in c(main_dir, main_pdf_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
letters_file <- file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv")
isotope_file <- file.path(ISOTOPE_DIR, "MTT_FYW.csv")
damping_file <- file.path(ISOTOPE_DIR, "DampingRatios_2025-07-07.csv")
fdc_wy_file <- file.path(OUT_MET_DYNAMIC_DIR, "fdc_slopes_wy.csv")

for (f in c(annual_file, site_file)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

annual <- read_csv(annual_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

# Ensure Fig2 FDC panel uses annual site-year slopes.
if (file.exists(fdc_wy_file)) {
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
}

site_summary <- read_csv(site_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

letters_df <- if (file.exists(letters_file)) {
  read_csv(letters_file, show_col_types = FALSE) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))
} else {
  tibble(metric = character(), site = factor(levels = SITE_ORDER_HYDROMETRIC), group_letter = character())
}

BOX_FILL_ALPHA <- 0.22
POINT_ALPHA <- 0.55
PANEL_BG <- "white"
PANEL_BORDER <- "grey55"

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

build_corr_triangle_panel <- function(
  data_df,
  metric_map,
  metric_labels,
  title,
  legend_mode = c("inset", "side", "side_inner", "bottom", "left"),
  axis_text_size = FIG_AXIS_TEXT_SIZE + 1,
  title_margin_bottom = 9,
  legend_scale = 1.32
) {
  legend_mode <- match.arg(legend_mode)
  corr_input <- data_df %>%
    select(any_of(unname(metric_map))) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

  if (ncol(corr_input) < 2) {
    return(patchwork::plot_spacer())
  }

  colnames(corr_input) <- names(metric_map)
  corr_mat <- suppressWarnings(cor(corr_input, use = "pairwise.complete.obs"))
  idx <- which(lower.tri(corr_mat), arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(patchwork::plot_spacer())
  }

  x_levels <- names(metric_map)[-length(metric_map)]
  y_levels <- names(metric_map)[-1]

  corr_pvals <- apply(idx, 1, function(rc) {
    x <- corr_input[[rc[2]]]
    y <- corr_input[[rc[1]]]
    keep <- is.finite(x) & is.finite(y)

    if (sum(keep) < 3) {
      return(NA_real_)
    }

    suppressWarnings(
      tryCatch(
        cor.test(x[keep], y[keep], method = "pearson")$p.value,
        error = function(e) NA_real_
      )
    )
  })

  corr_tri <- tibble(
    row_metric = rownames(corr_mat)[idx[, 1]],
    col_metric = colnames(corr_mat)[idx[, 2]],
    r = corr_mat[idx],
    p_value = as.numeric(corr_pvals)
  ) %>%
    mutate(
      row_metric = factor(row_metric, levels = y_levels),
      col_metric = factor(col_metric, levels = x_levels),
      r_label = ifelse(abs(r) < 0.005, 0, r),
      fontface = ifelse(is.finite(p_value) & p_value < 0.05, "bold", "plain"),
      label = ifelse(
        is.finite(r_label),
        sprintf("%.2f", r_label),
        ""
      )
    )

  p_corr <- ggplot(corr_tri, aes(x = col_metric, y = row_metric, fill = r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label, fontface = fontface), size = FIG_TILE_TEXT_SIZE * 1.15) +
    scale_fill_gradient2(
      low = "firebrick3",
      mid = "white",
      high = "dodgerblue3",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey85",
      name = "Pearson's r"
    ) +
    scale_x_discrete(
      labels = function(x) {
        parsed <- metric_labels[x]
        parsed[is.na(parsed)] <- x[is.na(parsed)]
        parse(text = unname(parsed))
      }
    ) +
    scale_y_discrete(
      labels = function(x) {
        parsed <- metric_labels[x]
        parsed[is.na(parsed)] <- x[is.na(parsed)]
        parse(text = unname(parsed))
      }
    ) +
    labs(x = NULL, y = NULL, title = title) +
    theme_storage_panel() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(
        size = FIG_STRIP_TEXT_SIZE + 2,
        hjust = 0,
        margin = margin(t = 1, r = 0, b = title_margin_bottom, l = 0)
      ),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      legend.title = element_text(size = FIG_AXIS_TEXT_SIZE),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE - 1),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.key.height = grid::unit(10, "pt"),
      legend.key.width = grid::unit(8, "pt"),
      plot.margin = margin(11, 0, 5.5, 5.5)
    )

  if (identical(legend_mode, "inset")) {
    return(
      p_corr +
        theme(
          legend.position = c(0.9, 0.35),
          legend.justification = c(1, 0.5),
          legend.background = element_rect(fill = scales::alpha("white", 0.85), color = NA)
        )
    )
  }

  if (identical(legend_mode, "bottom")) {
    return(
      p_corr +
        guides(
          fill = guide_colorbar(
            title.position = "top",
            direction = "horizontal",
            barwidth = grid::unit(48, "pt"),
            barheight = grid::unit(8, "pt")
          )
        ) +
        theme(
          legend.position = "bottom",
          legend.justification = c(0.5, 0.5),
          legend.direction = "horizontal",
          legend.background = element_blank()
        )
    )
  }

  if (identical(legend_mode, "left")) {
    corr_grob <- ggplotGrob(
      p_corr +
        theme(
          legend.position = "left",
          legend.justification = c(0.5, 0.5),
          legend.background = element_blank()
        )
    )
    guide_idx <- which(vapply(corr_grob$grobs, function(x) x$name, character(1)) == "guide-box")

    if (length(guide_idx) == 0) {
      return(p_corr + theme(legend.position = "none"))
    }

    legend_grob <- corr_grob$grobs[[guide_idx[1]]]
    p_corr_main <- p_corr + theme(legend.position = "none")

    return(
      (patchwork::wrap_elements(legend_grob) | p_corr_main) +
        plot_layout(widths = c(0.22, 1))
    )
  }

  corr_grob <- ggplotGrob(
    p_corr +
      theme(
        legend.position = "right",
        legend.justification = c(0.5, 0.5),
        legend.background = element_blank()
      )
  )
  guide_idx <- which(vapply(corr_grob$grobs, function(x) x$name, character(1)) == "guide-box")

  if (length(guide_idx) == 0) {
    return(p_corr + theme(legend.position = "none"))
  }

  legend_grob <- corr_grob$grobs[[guide_idx[1]]]
  legend_grob$widths <- legend_grob$widths * legend_scale
  legend_grob$heights <- legend_grob$heights * legend_scale
  p_corr_main <- p_corr + theme(legend.position = "none")

  if (identical(legend_mode, "side_inner")) {
    return(
      (p_corr_main + theme(plot.margin = margin(4, 0, 0, 5.5)) |
         patchwork::wrap_elements(legend_grob)) +
        plot_layout(widths = c(1, 0.07))
    )
  }

  (p_corr_main + theme(plot.margin = margin(4, 0, 0, 0)) |
     patchwork::wrap_elements(legend_grob)) +
    plot_layout(widths = c(1, 0.15))
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

# fig2 dynamic annual boxplots
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
  dat <- fig2_long %>% filter(metric == metric_key)
  let <- letters_fig2 %>% filter(metric == metric_key)
  y_lim <- calc_ylim(dat$value, let$y)
  show_x <- TRUE
  reserve_x_space <- FALSE
  panel_margin <- if (identical(metric_key, "FDC")) {
    margin(5.5, 5.5, -7, 5.5)
  } else if (identical(metric_key, "WB")) {
    # Add extra top spacing so panel e title does not crowd panel c above.
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
      color = "grey35"
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

p_fig2_a <- build_fig2_panel("RBI")
p_fig2_b <- build_fig2_panel("RCS")
p_fig2_c <- build_fig2_panel("FDC")
p_fig2_d <- build_fig2_panel("SD")
p_fig2_e <- build_fig2_panel("WB")
p_fig2_f <- build_corr_triangle_panel(
  annual,
  metric_map_fig2,
  metric_labels = c(
    RBI = "plain(RBI)",
    RCS = "plain(RCS)",
    FDC = "plain(FDC)",
    SD = "plain(SD)",
    WB = "plain(WB)"
  ),
  title = "f) Dynamic Storage Metric Correlations",
  legend_mode = "side",
  axis_text_size = FIG_AXIS_TEXT_SIZE,
  title_margin_bottom = 12,
  legend_scale = 1.18
)

# Build columns separately so panel d x-axis labels do not force extra blank
# space beneath panel c and push panel e downward.
p_fig2_left <- p_fig2_a / p_fig2_c / p_fig2_e +
  plot_layout(heights = c(1, 1, 1))
p_fig2_right <- p_fig2_b / p_fig2_d / p_fig2_f +
  plot_layout(heights = c(1, 1, 1))

p_fig2 <- (p_fig2_left | p_fig2_right) +
  plot_layout(widths = c(1, 1))

# legacy standalone Fig4 CHS boxplot content retained only to reuse the BF panel
# inside the consolidated mobile-storage figure (current Fig4).
fig4_df <- annual %>%
  select(site, year, CHS) %>%
  filter(is.finite(CHS))

letters_fig4 <- letters_df %>%
  filter(metric == "CHS") %>%
  left_join(
    fig4_df %>% group_by(site) %>% summarise(ymax_site = max(CHS, na.rm = TRUE), .groups = "drop"),
    by = "site"
  ) %>%
  mutate(y = ymax_site + 0.04)

chs_counts <- tibble(site = factor(SITE_ORDER_HYDROMETRIC, levels = SITE_ORDER_HYDROMETRIC)) %>%
  left_join(
    annual %>%
      group_by(site) %>%
      summarise(n = sum(is.finite(CHS)), .groups = "drop"),
    by = "site"
  ) %>%
  mutate(n = ifelse(is.na(n), 0L, as.integer(n)))

fig4_x_labels <- setNames(as.character(chs_counts$site), as.character(chs_counts$site))

p_fig4_legacy <- ggplot(fig4_df, aes(x = site, y = CHS, fill = site, color = site)) +
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
  coord_cartesian(ylim = calc_ylim(fig4_df$CHS, letters_fig4$y), clip = "off") +
  labs(x = NULL, y = "Baseflow Fraction (BF)") +
  theme_storage_panel() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = FIG_AXIS_TEXT_SIZE + 1),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1)
  )

# current Fig4 mobile storage panels
fig4_mobile_map <- c(
  DR = "DR",
  Fyw = "Fyw",
  MTT = "MTT"
)

isotope_err <- tibble(site = SITE_ORDER_HYDROMETRIC)

if (file.exists(isotope_file)) {
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
} else {
  isotope_err <- isotope_err %>%
    mutate(Fyw_err = NA_real_, MTT_err = NA_real_)
}

if (file.exists(damping_file)) {
  dr_err <- read_csv(damping_file, show_col_types = FALSE) %>%
    mutate(
      site = standardize_site_code(as.character(site)),
      DR_err = suppressWarnings(as.numeric(DR__err))
    ) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
    group_by(site) %>%
    summarise(DR_err = first_finite(DR_err), .groups = "drop")

  isotope_err <- isotope_err %>% left_join(dr_err, by = "site")
} else {
  isotope_err$DR_err <- NA_real_
}

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
  if (identical(metric_key, "BF")) {
    return(
      ggplot(fig4_df, aes(x = site, y = CHS, fill = site, color = site)) +
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
        coord_cartesian(ylim = calc_ylim(fig4_df$CHS, letters_fig4$y), clip = "off") +
        labs(x = NULL, y = fig4_y_labels[[metric_key]]) +
        ggtitle(fig4_titles[[metric_key]]) +
        theme_storage_panel() +
        theme(
          aspect.ratio = 0.6,
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
  # Match panel layout: keep x-axis labels only on bottom row.
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
      aspect.ratio = 0.6,
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
p_fig4_e <- build_corr_triangle_panel(
  site_summary,
  c(
    BF = "CHS_mean",
    DR = "DR",
    Fyw = "Fyw",
    MTT = "MTT"
  ),
  metric_labels = c(
    BF = "plain(BF)",
    DR = "plain(DR)",
    Fyw = "plain(F)[yw]",
    MTT = "plain(MTT)"
  ),
  title = "e) Mobile Storage Metric Correlations",
  legend_mode = "side",
  axis_text_size = FIG_AXIS_TEXT_SIZE + 2,
  title_margin_bottom = 9,
  legend_scale = 1.65
)

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

p_fig4 <- (
  (fig4_raw_plots[[1]] | fig4_raw_plots[[2]]) /
  (fig4_raw_plots[[3]] | fig4_raw_plots[[4]]) /
  (patchwork::wrap_elements(full = patchwork::patchworkGrob(p_fig4_e)) | patchwork::plot_spacer())
) +
  patchwork::plot_layout(
    widths = c(1, 1),
    heights = c(1, 1, 1.1)
  )

ggsave(
  file.path(main_dir, "Fig2_ds_annual_boxplots.png"),
  p_fig2,
  width = 10.4 * FIG_WIDTH_SCALE,
  height = 12.0 * FIG_HEIGHT_SCALE,
  bg = "white",
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "Fig2_ds_annual_boxplots.pdf"),
  p_fig2,
  width = 10.4 * FIG_WIDTH_SCALE,
  height = 12.0 * FIG_HEIGHT_SCALE,
  bg = "white"
)

ggsave(
  file.path(main_dir, "Fig4_mobile_storage.png"),
  p_fig4,
  width = 10.4 * FIG_WIDTH_SCALE,
  height = 11.8 * FIG_HEIGHT_SCALE,
  bg = "white",
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "Fig4_mobile_storage.pdf"),
  p_fig4,
  width = 10.4 * FIG_WIDTH_SCALE,
  height = 11.8 * FIG_HEIGHT_SCALE,
  bg = "white"
)
