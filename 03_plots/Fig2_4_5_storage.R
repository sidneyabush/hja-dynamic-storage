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

for (f in c(annual_file, site_file)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

annual <- read_csv(annual_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

site_summary <- read_csv(site_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

letters_df <- if (file.exists(letters_file)) {
  read_csv(letters_file, show_col_types = FALSE) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))
} else {
  tibble(metric = character(), site = factor(levels = SITE_ORDER_HYDROMETRIC), group_letter = character())
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
  select(site, year, any_of(metric_map_fig2)) %>%
  pivot_longer(cols = all_of(metric_map_fig2), names_to = "metric", values_to = "value") %>%
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

p_fig2 <- ggplot(fig2_long, aes(x = site, y = value, fill = site)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(color = site), width = 0.15, alpha = 0.35, size = FIG_POINT_SIZE_SMALL) +
  geom_text(
    data = letters_fig2,
    aes(x = site, y = y, label = group_letter),
    inherit.aes = FALSE,
    size = FIG_ANNOT_TEXT_SIZE,
    color = "black"
  ) +
  scale_fill_manual(values = SITE_COLORS) +
  scale_color_manual(values = SITE_COLORS) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT_SIZE),
    axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
    strip.text = element_text(size = FIG_STRIP_TEXT_SIZE),
    legend.position = "none"
  )

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
  mutate(y = ymax_site + 0.03)

p_fig4 <- ggplot(fig4_df, aes(x = site, y = CHS, fill = site)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(color = site), width = 0.15, alpha = 0.35, size = FIG_POINT_SIZE_SMALL) +
  geom_text(
    data = letters_fig4,
    aes(x = site, y = y, label = group_letter),
    inherit.aes = FALSE,
    size = FIG_ANNOT_TEXT_SIZE,
    color = "black"
  ) +
  scale_fill_manual(values = SITE_COLORS) +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = NULL, y = "CHS") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT_SIZE),
    axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    legend.position = "none"
  )

# fig5 isotope metrics
fig5_map <- c(
  DR = "DR",
  Fyw = "Fyw",
  MTT1 = "MTT1",
  MTT2 = "MTT2"
)

fig5_df <- site_summary %>%
  select(site, any_of(fig5_map)) %>%
  pivot_longer(cols = all_of(fig5_map), names_to = "metric", values_to = "value") %>%
  filter(is.finite(value)) %>%
  mutate(metric = factor(metric, levels = names(fig5_map)))

p_fig5 <- ggplot(fig5_df, aes(x = site, y = value, color = site)) +
  geom_point(size = FIG_POINT_SIZE_LARGE + 1.2, alpha = 0.95) +
  scale_color_manual(values = SITE_COLORS) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT_SIZE),
    axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
    strip.text = element_text(size = FIG_STRIP_TEXT_SIZE),
    legend.position = "none"
  )

ggsave(
  file.path(main_dir, "Fig2_ds_annual_boxplots.png"),
  p_fig2,
  width = 13.2 * FIG_WIDTH_SCALE,
  height = 8.5 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "Fig2_ds_annual_boxplots.pdf"),
  p_fig2,
  width = 13.2 * FIG_WIDTH_SCALE,
  height = 8.5 * FIG_HEIGHT_SCALE
)

ggsave(
  file.path(main_dir, "Fig4_chs_boxplots.png"),
  p_fig4,
  width = 7.2 * FIG_WIDTH_SCALE,
  height = 5.6 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "Fig4_chs_boxplots.pdf"),
  p_fig4,
  width = 7.2 * FIG_WIDTH_SCALE,
  height = 5.6 * FIG_HEIGHT_SCALE
)

ggsave(
  file.path(main_dir, "Fig5_iso_annual.png"),
  p_fig5,
  width = 10.6 * FIG_WIDTH_SCALE,
  height = 7.0 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "Fig5_iso_annual.pdf"),
  p_fig5,
  width = 10.6 * FIG_WIDTH_SCALE,
  height = 7.0 * FIG_HEIGHT_SCALE
)
