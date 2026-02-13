# -----------------------------------------------------------------------------
# Recession Relationship Plots: log(Q) vs log(-dQ/dt)
# -----------------------------------------------------------------------------
# Faceted relationship plot by site for supplement figures.
#
# Inputs:
#   - HF00402_v14.csv
#   - drainage_area.csv
#
# Outputs (FIGURES_DIR/supp/hydrometric):
#   - rcs_curve_by_site.png
#   - rcs_curve_by_site.pdf
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)

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
site_labels_panel <- make_panel_label_map(SITE_ORDER_HYDROMETRIC)

da_path <- resolve_drainage_area_file()
q_path <- file.path(DISCHARGE_DIR, "HF00402_v14.csv")

if (!file.exists(da_path) || !file.exists(q_path)) {
  stop("Required discharge files not found in DISCHARGE_DIR.")
}

da_df <- read_csv(da_path, show_col_types = FALSE)

discharge <- read_csv(q_path, show_col_types = FALSE) %>%
  mutate(
    site = standardize_site_code(SITECODE),
    date = as.Date(DATE, format = "%m/%d/%Y")
  ) %>%
  filter(
    !is.na(site),
    site %in% SITE_ORDER_HYDROMETRIC,
    WATERYEAR >= WY_START
  ) %>%
  left_join(da_df, by = "SITECODE") %>%
  filter(!is.na(DA_M2)) %>%
  mutate(
    Q_cms = MEAN_Q * 0.02831683199881,
    Q_mm_day = Q_cms / DA_M2 * 86400 * 1000
  ) %>%
  arrange(site, date)

recession_clean <- discharge %>%
  group_by(site) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    dQ = Q_cms - lag(Q_cms),
    dt_days = as.numeric(date - lag(date)),
    dQ_dt = dQ / dt_days,
    recession_slope = -dQ_dt,
    slope_mm_day = recession_slope / DA_M2 * 86400 * 1000
  ) %>%
  ungroup() %>%
  filter(
    !is.na(Q_mm_day),
    !is.na(slope_mm_day),
    Q_mm_day > 0,
    slope_mm_day > 0
  ) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

slopes_df <- recession_clean %>%
  group_by(site) %>%
  summarise(
    slope_mm = abs(round(coef(lm(log(slope_mm_day) ~ log(Q_mm_day)))[2], 2)),
    x_pos = quantile(log(Q_mm_day), 0.85, na.rm = TRUE),
    y_pos = quantile(log(slope_mm_day), 0.85, na.rm = TRUE),
    .groups = "drop"
  )

p_curve <- ggplot(recession_clean, aes(x = log(Q_mm_day), y = log(slope_mm_day), color = site)) +
  geom_point(alpha = 0.5, size = FIG_POINT_SIZE_SMALL) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  geom_text(
    data = slopes_df,
    aes(x = x_pos, y = y_pos, label = paste0("slope = ", slope_mm)),
    inherit.aes = FALSE,
    color = "black",
    size = FIG_ANNOT_TEXT_SIZE,
    hjust = 0,
    check_overlap = FIG_LABEL_CHECK_OVERLAP
  ) +
  facet_wrap(
    ~site,
    ncol = 2,
    scales = "fixed",
    drop = FALSE,
    axes = "margins",
    axis.labels = "margins",
    labeller = labeller(site = site_labels_panel)
  ) +
  scale_color_manual(values = SITE_COLORS, guide = "none") +
  labs(x = "log Q (mm day-1)", y = "log -dQ/dt (mm day-1)") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT)
  ) +
  coord_cartesian(clip = FIG_LABEL_CLIP)

ggsave(
  file.path(plot_dir, "rcs_curve_by_site.png"),
  p_curve,
  width = 8 * FIG_WIDTH_SCALE,
  height = 9 * FIG_HEIGHT_SCALE,
  dpi = 300
)

ggsave(
  file.path(plot_dir, "rcs_curve_by_site.pdf"),
  p_curve,
  width = 8 * FIG_WIDTH_SCALE,
  height = 9 * FIG_HEIGHT_SCALE
)
