# figures 2 4 5 main manuscript plots

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
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
  RBI = "a) RBI",
  RCS = "b) RCS",
  FDC = "c) FDC",
  SD = "d) SD",
  WB = "e) WB"
)
fig2_y_labels <- c(
  RBI = "Unitless",
  RCS = "Unitless",
  FDC = "Unitless",
  SD = "Depth (mm)",
  WB = "Depletion (mm)"
)

build_fig2_panel <- function(metric_key) {
  dat <- fig2_long %>% filter(metric == metric_key)
  let <- letters_fig2 %>% filter(metric == metric_key)
  y_lim <- calc_ylim(dat$value, let$y)
  show_x <- metric_key %in% c("SD", "WB")
  panel_margin <- if (identical(metric_key, "FDC")) {
    margin(5.5, 5.5, -7, 5.5)
  } else if (identical(metric_key, "WB")) {
    margin(0, 5.5, 5.5, 5.5)
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
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      plot.margin = panel_margin
    )
}

p_fig2_a <- build_fig2_panel("RBI")
p_fig2_b <- build_fig2_panel("RCS")
p_fig2_c <- build_fig2_panel("FDC")
p_fig2_d <- build_fig2_panel("SD")
p_fig2_e <- build_fig2_panel("WB")

# Build columns separately so panel d x-axis labels do not force extra blank
# space beneath panel c and push panel e downward.
p_fig2_left <- p_fig2_a / p_fig2_c / p_fig2_e +
  plot_layout(heights = c(1, 1, 1))
p_fig2_right <- p_fig2_b / p_fig2_d / patchwork::plot_spacer() +
  plot_layout(heights = c(1, 1, 1))

p_fig2 <- p_fig2_left | p_fig2_right

# fig4 chs boxplot
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

fig4_x_labels <- setNames(
  paste0(as.character(chs_counts$site), "\n(n=", chs_counts$n, ")"),
  as.character(chs_counts$site)
)

p_fig4 <- ggplot(fig4_df, aes(x = site, y = CHS, fill = site, color = site)) +
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
  labs(x = NULL, y = "Baseflow Fraction (CHS)") +
  theme_storage_panel() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT_SIZE + 1),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1)
  )

# fig5 isotope metrics
fig5_map <- c(
  DR = "DR",
  Fyw = "Fyw",
  MTT1 = "MTT1",
  MTT2 = "MTT2"
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
      MTT1_err = suppressWarnings(as.numeric(MTT1_SD)),
      MTT2_err = rowMeans(
        cbind(
          suppressWarnings(as.numeric(MTT2L_SD)),
          suppressWarnings(as.numeric(MTT2H_SD))
        ),
        na.rm = TRUE
      ),
      Fyw_err = ifelse(is.nan(Fyw_err), NA_real_, Fyw_err),
      MTT2_err = ifelse(is.nan(MTT2_err), NA_real_, MTT2_err)
    ) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
    group_by(site) %>%
    summarise(
      Fyw_err = first_finite(Fyw_err),
      MTT1_err = first_finite(MTT1_err),
      MTT2_err = first_finite(MTT2_err),
      .groups = "drop"
    )

  isotope_err <- isotope_err %>% left_join(mtt_err, by = "site")
} else {
  isotope_err <- isotope_err %>%
    mutate(Fyw_err = NA_real_, MTT1_err = NA_real_, MTT2_err = NA_real_)
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
    cols = c(DR_err, Fyw_err, MTT1_err, MTT2_err),
    names_to = "metric",
    values_to = "err"
  ) %>%
  mutate(metric = sub("_err$", "", metric))

fig5_df <- site_summary %>%
  select(site, any_of(unname(fig5_map))) %>%
  pivot_longer(cols = all_of(unname(fig5_map)), names_to = "metric", values_to = "value") %>%
  left_join(fig5_err_long, by = c("site", "metric")) %>%
  filter(is.finite(value)) %>%
  mutate(metric = factor(metric, levels = names(fig5_map)))

fig5_titles <- c(
  DR = "a) DR",
  Fyw = "b) Fyw",
  MTT1 = "c) MTT1",
  MTT2 = "d) MTT2"
)
fig5_y_labels <- c(
  DR = "Ratio",
  Fyw = "Fraction",
  MTT1 = "Years",
  MTT2 = "Years"
)

build_fig5_panel <- function(metric_key) {
  dat <- fig5_df %>% filter(metric == metric_key)
  show_x <- metric_key %in% c("MTT1", "MTT2")
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
      y = fig5_y_labels[[metric_key]],
      title = fig5_titles[[metric_key]]
    ) +
    theme_storage_panel() +
    theme(
      axis.text.x = if (show_x) {
        element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT_SIZE)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) element_line() else element_blank()
    )
}

p_fig5 <- wrap_plots(
  lapply(names(fig5_map), build_fig5_panel),
  ncol = 2,
  byrow = TRUE
)

ggsave(
  file.path(main_dir, "Fig2_ds_annual_boxplots.png"),
  p_fig2,
  width = 10.4 * FIG_WIDTH_SCALE,
  height = 10.2 * FIG_HEIGHT_SCALE,
  bg = "white",
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "Fig2_ds_annual_boxplots.pdf"),
  p_fig2,
  width = 10.4 * FIG_WIDTH_SCALE,
  height = 10.2 * FIG_HEIGHT_SCALE,
  bg = "white"
)

ggsave(
  file.path(main_dir, "Fig4_chs_boxplots.png"),
  p_fig4,
  width = 10.0 * FIG_WIDTH_SCALE,
  height = 6.2 * FIG_HEIGHT_SCALE,
  bg = "white",
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "Fig4_chs_boxplots.pdf"),
  p_fig4,
  width = 10.0 * FIG_WIDTH_SCALE,
  height = 6.2 * FIG_HEIGHT_SCALE,
  bg = "white"
)

ggsave(
  file.path(main_dir, "Fig5_iso_annual.png"),
  p_fig5,
  width = 10.8 * FIG_WIDTH_SCALE,
  height = 8.4 * FIG_HEIGHT_SCALE,
  bg = "white",
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "Fig5_iso_annual.pdf"),
  p_fig5,
  width = 10.8 * FIG_WIDTH_SCALE,
  height = 8.4 * FIG_HEIGHT_SCALE,
  bg = "white"
)
