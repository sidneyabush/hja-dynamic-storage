# figure 9 conceptual diagram

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

rm(list = ls())
source("config.R")

main_dir <- MS_FIG_MAIN_DIR
main_pdf_dir <- MS_FIG_MAIN_PDF_DIR
explore_dir <- CONCEPTUAL_DIAGRAM_DIR
for (d in c(main_dir, main_pdf_dir, explore_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

axes_file <- file.path(OUTPUT_DIR, "models", "unified_framework", "unified_framework_site_axes.csv")
loadings_file <- file.path(OUTPUT_DIR, "models", "unified_framework", "geology_composition_pca_loadings.csv")
variance_file <- file.path(OUTPUT_DIR, "models", "unified_framework", "geology_composition_pca_variance.csv")
master_site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)

for (f in c(axes_file, loadings_file, variance_file, master_site_file)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

axes <- read_csv(axes_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site)

loadings <- read_csv(loadings_file, show_col_types = FALSE)
variance <- read_csv(variance_file, show_col_types = FALSE)
master_site <- read_csv(master_site_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

if (!all(c("geology_pc1", "geology_pc2") %in% names(axes))) {
  stop("unified_framework_site_axes.csv is missing geology_pc1/geology_pc2")
}

if (!all(c("geology_pc1_loading", "geology_pc2_loading") %in% names(loadings))) {
  stop("geology_composition_pca_loadings.csv is missing pca loading columns")
}

pc1_pct <- ifelse(nrow(variance) >= 1, 100 * variance$variance_explained[1], NA_real_)
pc2_pct <- ifelse(nrow(variance) >= 2, 100 * variance$variance_explained[2], NA_real_)

pca_scores <- axes %>%
  transmute(
    site = as.character(site),
    geology_pc1 = as.numeric(geology_pc1),
    geology_pc2 = as.numeric(geology_pc2)
  ) %>%
  filter(is.finite(geology_pc1), is.finite(geology_pc2))

loading_df <- loadings %>%
  transmute(
    variable = as.character(variable),
    lx = as.numeric(geology_pc1_loading),
    ly = as.numeric(geology_pc2_loading)
  ) %>%
  filter(is.finite(lx), is.finite(ly))

pretty_loading <- function(x) {
  out <- as.character(x)
  out[out == "Lava1_per"] <- "Lava-1"
  out[out == "Lava2_per"] <- "Lava-2"
  out[out == "Ash_Per"] <- "Ash"
  out[out == "Pyro_per"] <- "Pyro"
  out[out == "Landslide_Young"] <- "Landslide Young"
  out[out == "Landslide_Mod"] <- "Landslide Mod"
  out[out == "Landslide_Old"] <- "Landslide Old"
  out[out == "Landslide_Total"] <- "Landslide Total"
  out
}

loading_df <- loading_df %>% mutate(label = pretty_loading(variable))

# scale loading arrows to pca score space
score_lim_x <- max(abs(pca_scores$geology_pc1), na.rm = TRUE)
score_lim_y <- max(abs(pca_scores$geology_pc2), na.rm = TRUE)
load_lim <- max(abs(c(loading_df$lx, loading_df$ly)), na.rm = TRUE)
arrow_scale <- ifelse(is.finite(load_lim) && load_lim > 0, 0.36 * max(score_lim_x, score_lim_y) / load_lim, 1)

loading_plot <- loading_df %>%
  mutate(
    xend = lx * arrow_scale,
    yend = ly * arrow_scale
  )

site_label_size <- FIG_ANNOT_TEXT_SIZE + 1.2
loading_label_size <- FIG_ANNOT_TEXT_SIZE + 0.3

p_geo <- ggplot(pca_scores, aes(x = geology_pc1, y = geology_pc2)) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey70") +
  geom_point(aes(color = site), size = FIG_POINT_SIZE_LARGE + 2.1, alpha = 0.98) +
  geom_segment(
    data = loading_plot,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.16, "cm")),
    color = "grey45",
    linewidth = 0.5
  ) +
  geom_label_repel(
    data = loading_plot,
    aes(x = xend, y = yend, label = label),
    inherit.aes = FALSE,
    color = "grey45",
    fill = alpha("white", 0.55),
    label.size = 0,
    size = loading_label_size,
    box.padding = 0.12,
    point.padding = 0.08,
    segment.color = alpha("grey60", 0.55),
    force = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  geom_text_repel(
    aes(label = site, color = site),
    size = site_label_size,
    box.padding = 0.32,
    point.padding = 0.28,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = paste0("PC1 (", number(pc1_pct, accuracy = 0.1), "%)"),
    y = paste0("PC2 (", number(pc2_pct, accuracy = 0.1), "%)")
  ) +
  theme_pub() +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1.5),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1.6),
    legend.position = "none"
  )

# complete-case conceptual panel
# require chs and isotope metrics so mobile axis is fully supported
complete_mobile_sites <- master_site %>%
  mutate(
    has_mtt = is.finite(MTT1) | is.finite(MTT2),
    complete_mobile = is.finite(CHS_mean) & is.finite(DR) & is.finite(Fyw) & has_mtt
  ) %>%
  transmute(site, complete_mobile)

plot_b <- axes %>%
  mutate(site = as.character(site)) %>%
  left_join(complete_mobile_sites, by = "site") %>%
  transmute(
    site,
    dynamic = as.numeric(dynamic_storage_strength_z),
    mobile = as.numeric(mobile_mixing_with_chs_z),
    pc1 = as.numeric(geology_pc1),
    pc2 = as.numeric(geology_pc2),
    complete_mobile = ifelse(is.na(complete_mobile), FALSE, complete_mobile)
  ) %>%
  filter(complete_mobile, is.finite(dynamic), is.finite(mobile), is.finite(pc1), is.finite(pc2))

if (nrow(plot_b) < 3) {
  stop("Not enough complete-case sites for panel b")
}

fit <- lm(mobile ~ dynamic, data = plot_b)
fit_sum <- summary(fit)
slope <- unname(coef(fit)[2])
r2 <- unname(fit_sum$r.squared)
p_val <- unname(coef(fit_sum)[2, 4])

ann_stats <- sprintf("R² = %.3f\np = %.3f\nb = %.3f", r2, p_val, slope)

x_rng <- range(plot_b$dynamic, na.rm = TRUE)
y_rng <- range(plot_b$mobile, na.rm = TRUE)

quad_df <- tibble(
  x = c(x_rng[1], x_rng[2], x_rng[1], x_rng[2]),
  y = c(y_rng[2], y_rng[2], y_rng[1], y_rng[1]),
  hjust = c(0, 1, 0, 1),
  vjust = c(1, 1, 0, 0),
  label = c(
    "low dynamic\nhigh mobile",
    "high dynamic\nhigh mobile",
    "low dynamic\nlow mobile",
    "high dynamic\nlow mobile"
  )
)

p_state <- ggplot(plot_b, aes(x = dynamic, y = mobile)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  geom_point(
    aes(fill = pc1, size = pc2),
    shape = 21,
    color = "grey35",
    stroke = 0.35,
    alpha = 0.96
  ) +
  geom_text_repel(
    aes(label = site),
    size = site_label_size,
    box.padding = 0.32,
    point.padding = 0.30,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  geom_label(
    data = quad_df,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE,
    size = FIG_ANNOT_TEXT_SIZE + 0.8,
    color = "grey25",
    fill = alpha("grey95", 0.92),
    label.size = 0,
    label.padding = unit(0.25, "lines")
  ) +
  annotate(
    "text",
    x = x_rng[2] - 0.04 * diff(x_rng),
    y = y_rng[2] - 0.10 * diff(y_rng),
    label = ann_stats,
    hjust = 1,
    vjust = 1,
    size = FIG_ANNOT_TEXT_SIZE + 1.1
  ) +
  scale_fill_gradient2(
    low = "#1b9e77",
    mid = "white",
    high = "#d95f02",
    midpoint = 0,
    name = "PC1"
  ) +
  scale_size_continuous(
    name = "PC2",
    breaks = c(-2, -1, 0, 1),
    range = c(3.5, 9)
  ) +
  labs(
    x = "Dynamic Storage",
    y = "Mobile Storage"
  ) +
  theme_pub() +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1.5),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1.6),
    legend.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1.3),
    legend.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1.1),
    legend.position = "right"
  ) +
  guides(fill = guide_colorbar(order = 1), size = guide_legend(order = 2))

fig9 <- p_geo / p_state +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "plain", size = FIG_AXIS_TITLE_SIZE + 1))

for (nm in c("Fig9_conceptual_diagram", "geo_landslide_pca_dynamic_mobile_storage_geology_pcsize_all_sites_stacked")) {
  ggsave(
    file.path(main_dir, paste0(nm, ".png")),
    fig9,
    width = 10.2 * FIG_WIDTH_SCALE,
    height = 13.2 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
  ggsave(
    file.path(main_pdf_dir, paste0(nm, ".pdf")),
    fig9,
    width = 10.2 * FIG_WIDTH_SCALE,
    height = 13.2 * FIG_HEIGHT_SCALE
  )
  ggsave(
    file.path(explore_dir, paste0(nm, ".png")),
    fig9,
    width = 10.2 * FIG_WIDTH_SCALE,
    height = 13.2 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
}
