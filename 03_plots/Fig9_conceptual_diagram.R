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

safe_ggsave <- function(filename, plot_obj, width, height, dpi = NULL) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ext <- tools::file_ext(filename)
  tmp_file <- tempfile(
    pattern = "fig9_",
    tmpdir = tempdir(),
    fileext = ifelse(nzchar(ext), paste0(".", ext), "")
  )
  tryCatch(
    {
      if (is.null(dpi)) {
        ggplot2::ggsave(tmp_file, plot_obj, width = width, height = height, bg = "white")
      } else {
        ggplot2::ggsave(tmp_file, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
      }
      ok <- file.copy(tmp_file, filename, overwrite = TRUE)
      unlink(tmp_file)
      if (!isTRUE(ok)) {
        stop("Failed to copy rendered file to destination")
      }
      TRUE
    },
    error = function(e) {
      unlink(tmp_file)
      warning("Failed to save plot: ", filename, " (", conditionMessage(e), ")")
      FALSE
    }
  )
}

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
  mutate(site = standardize_site_code(as.character(site))) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  group_by(site) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site)

loadings <- read_csv(loadings_file, show_col_types = FALSE)
variance <- read_csv(variance_file, show_col_types = FALSE)
master_site <- read_csv(master_site_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(as.character(site))) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

missing_biplot_sites <- setdiff(SITE_ORDER_HYDROMETRIC, as.character(axes$site))
if (length(missing_biplot_sites) > 0) {
  warning(
    "Sites missing from unified_framework_site_axes.csv for biplot: ",
    paste(missing_biplot_sites, collapse = ", ")
  )
}

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

FIG9_EXPORT_SCALE <- 1.20
FIG9_ELEMENT_SCALE <- 1.45
FIG9_TEXT_SCALE <- 1.45 * FIG9_ELEMENT_SCALE
FIG9_PCA_POINT_SCALE <- 1.85
site_label_size <- (FIG_ANNOT_TEXT_SIZE + 1.2) * FIG9_TEXT_SCALE
loading_label_size <- (FIG_ANNOT_TEXT_SIZE + 0.3) * FIG9_TEXT_SCALE
axis_text_size <- (FIG_AXIS_TEXT_SIZE + 1.5) * FIG9_TEXT_SCALE
axis_title_size <- (FIG_AXIS_TITLE_SIZE + 1.6) * FIG9_TEXT_SCALE
legend_title_size <- (FIG_AXIS_TITLE_SIZE + 1.3) * FIG9_TEXT_SCALE
legend_text_size <- (FIG_AXIS_TEXT_SIZE + 1.1) * FIG9_TEXT_SCALE
quad_label_size <- (FIG_ANNOT_TEXT_SIZE + 0.8) * FIG9_TEXT_SCALE
stats_label_size <- (FIG_ANNOT_TEXT_SIZE + 0.9) * FIG9_TEXT_SCALE * 0.96
tag_label_size <- (FIG_AXIS_TITLE_SIZE + 1) * FIG9_TEXT_SCALE
pca_point_size <- (FIG_POINT_SIZE_LARGE + 1.8) * FIG9_ELEMENT_SCALE * FIG9_PCA_POINT_SCALE
pca_point_stroke <- 0.35 * FIG9_ELEMENT_SCALE * 1.20
pca_label_box_padding <- unit(0.65, "lines")
pca_label_point_padding <- unit(1.20, "lines")
state_label_box_padding <- unit(0.36, "lines")
state_label_point_padding <- unit(0.72, "lines")
pca_label_force <- 2.20
pca_label_force_pull <- 0.55
pca_label_nudge_scale <- 0.32
ws9_ws10_nudge <- 0.34
state_label_force <- 0.60
state_label_force_pull <- 1.30
state_label_nudge_scale <- 0.22

pca_labels_general <- pca_scores %>%
  filter(!site %in% c("Mack", "Look", "WS09", "WS10")) %>%
  mutate(
    label_radius = sqrt(geology_pc1^2 + geology_pc2^2),
    nudge_x = if_else(
      label_radius > 0,
      pca_label_nudge_scale * geology_pc1 / label_radius,
      0
    ),
    nudge_y = if_else(
      label_radius > 0,
      pca_label_nudge_scale * geology_pc2 / label_radius,
      pca_label_nudge_scale
    )
  )

pca_labels_ws9_ws10 <- pca_scores %>%
  filter(site %in% c("WS09", "WS10")) %>%
  mutate(
    nudge_x = -ws9_ws10_nudge,
    nudge_y = if_else(site == "WS09", ws9_ws10_nudge * 0.95, -ws9_ws10_nudge * 0.95)
  )

p_geo <- ggplot(pca_scores, aes(x = geology_pc1, y = geology_pc2)) +
  geom_hline(yintercept = 0, linewidth = 0.5 * FIG9_ELEMENT_SCALE, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.5 * FIG9_ELEMENT_SCALE, color = "grey70") +
  geom_segment(
    data = loading_plot,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.16 * FIG9_ELEMENT_SCALE, "cm")),
    color = "grey45",
    linewidth = 0.5 * FIG9_ELEMENT_SCALE
  ) +
  geom_label_repel(
    data = loading_plot,
    aes(x = xend, y = yend, label = label),
    inherit.aes = FALSE,
    seed = 20260228,
    color = "grey45",
    fill = alpha("white", 0.45),
    label.size = 0,
    label.padding = unit(0.08, "lines"),
    label.r = unit(0.10, "lines"),
    size = loading_label_size,
    box.padding = 0.10,
    point.padding = 0.12,
    force = 0.9,
    force_pull = 0.95,
    min.segment.length = 0,
    segment.size = 0.23 * FIG9_ELEMENT_SCALE,
    segment.alpha = 0.75,
    segment.color = "grey55",
    max.overlaps = Inf
  ) +
  geom_point(
    aes(fill = site),
    shape = 21,
    color = "grey20",
    stroke = pca_point_stroke,
    size = pca_point_size,
    alpha = 0.98
  ) +
  geom_text_repel(
    data = pca_labels_general,
    aes(label = site, color = site),
    seed = 20260228,
    nudge_x = pca_labels_general$nudge_x,
    nudge_y = pca_labels_general$nudge_y,
    size = site_label_size,
    box.padding = pca_label_box_padding,
    point.padding = pca_label_point_padding,
    force = pca_label_force,
    force_pull = pca_label_force_pull,
    min.segment.length = Inf,
    segment.color = NA,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = pca_labels_ws9_ws10,
    aes(label = site, color = site),
    seed = 20260228,
    nudge_x = pca_labels_ws9_ws10$nudge_x,
    nudge_y = pca_labels_ws9_ws10$nudge_y,
    size = site_label_size,
    box.padding = pca_label_box_padding,
    point.padding = pca_label_point_padding,
    force = pca_label_force,
    force_pull = pca_label_force_pull,
    min.segment.length = Inf,
    segment.color = NA,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = pca_scores %>% filter(site == "Mack"),
    aes(label = site, color = site),
    seed = 20260228,
    size = site_label_size,
    nudge_y = -0.24,
    box.padding = pca_label_box_padding,
    point.padding = pca_label_point_padding,
    force = pca_label_force,
    force_pull = pca_label_force_pull,
    min.segment.length = Inf,
    segment.color = NA,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = pca_scores %>% filter(site == "Look"),
    aes(label = site, color = site),
    seed = 20260228,
    size = site_label_size,
    nudge_y = 0.24,
    box.padding = pca_label_box_padding,
    point.padding = pca_label_point_padding,
    force = pca_label_force,
    force_pull = pca_label_force_pull,
    min.segment.length = Inf,
    segment.color = NA,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = SITE_COLORS) +
  scale_fill_manual(values = SITE_COLORS) +
  labs(
    x = paste0("PC1 (", number(pc1_pct, accuracy = 0.1), "%)"),
    y = paste0("PC2 (", number(pc2_pct, accuracy = 0.1), "%)")
  ) +
  coord_cartesian(clip = "off") +
  theme_pub() +
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    legend.position = "none",
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  )

required_mobile_dynamic_cols <- c(
  "RBI_mean", "RCS_mean", "FDC_mean", "SD_mean", "WB_mean",
  "CHS_mean", "DR", "Fyw", "MTT1", "MTT2"
)
missing_required_cols <- setdiff(required_mobile_dynamic_cols, names(master_site))
if (length(missing_required_cols) > 0) {
  stop(
    "master_site is missing required columns for panel b filter: ",
    paste(missing_required_cols, collapse = ", ")
  )
}

eligible_sites_panel_b <- master_site %>%
  transmute(
    site = as.character(site),
    eligible = if_all(
      all_of(c("RBI_mean", "RCS_mean", "FDC_mean", "SD_mean", "WB_mean", "CHS_mean", "DR", "Fyw")),
      ~ is.finite(.x)
    ) & (is.finite(MTT1) | is.finite(MTT2))
  ) %>%
  filter(eligible) %>%
  pull(site)

plot_b <- axes %>%
  mutate(site = as.character(site)) %>%
  transmute(
    site,
    dynamic = as.numeric(dynamic_storage_strength_z),
    mobile = as.numeric(mobile_mixing_with_chs_z),
    pc1 = as.numeric(geology_pc1),
    pc2 = as.numeric(geology_pc2)
  ) %>%
  filter(
    site %in% eligible_sites_panel_b,
    is.finite(dynamic),
    is.finite(mobile),
    is.finite(pc1),
    is.finite(pc2)
  )

if (nrow(plot_b) < 3) {
  stop("Not enough complete-case sites for panel b")
}

fit <- lm(mobile ~ dynamic, data = plot_b)
fit_sum <- summary(fit)
slope <- unname(coef(fit)[2])
r2 <- unname(fit_sum$r.squared)
p_val <- unname(coef(fit_sum)[2, 4])

ann_stats <- sprintf("italic(R)^2 == %.3f~~italic(p) == %.3f~~italic(b) == %.3f", r2, p_val, slope)

x_rng <- range(plot_b$dynamic, na.rm = TRUE)
y_rng <- range(plot_b$mobile, na.rm = TRUE)
x_span <- diff(x_rng)
y_span <- diff(y_rng)
if (!is.finite(x_span) || x_span <= 0) x_span <- 1
if (!is.finite(y_span) || y_span <= 0) y_span <- 1
x_pad <- 0.40 * x_span
y_pad <- 0.40 * y_span
x_lim <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
y_lim <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)
stats_annot_x <- mean(x_lim)
stats_annot_y <- y_lim[1] - 0.34 * y_span

quad_df <- tibble(
  x = c(x_lim[1], x_lim[2], x_lim[1], x_lim[2]),
  y = c(y_lim[2], y_lim[2], y_lim[1], y_lim[1]),
  hjust = c(0, 1, 0, 1),
  vjust = c(1, 1, 0, 0),
  label = c(
    "low dynamic\nhigh mobile",
    "high dynamic\nhigh mobile",
    "low dynamic\nlow mobile",
    "high dynamic\nlow mobile"
  )
) %>%
  mutate(
    x = if_else(label == "high dynamic\nhigh mobile", x + 0.12 * x_span, x),
    y = if_else(label == "high dynamic\nhigh mobile", y + 0.08 * y_span, y)
  )

label_nudge_x_br <- max(0.15, 0.10 * x_span)
label_nudge_y_br <- -max(0.10, 0.08 * y_span)
label_nudge_x_ws07 <- -max(0.12, 0.08 * x_span)
label_nudge_y_ws07 <- max(0.05, 0.04 * y_span)

plot_b_labels <- plot_b %>%
  mutate(
    nudge_x = if_else(site == "WS07", label_nudge_x_ws07, label_nudge_x_br),
    nudge_y = if_else(site == "WS07", label_nudge_y_ws07, label_nudge_y_br)
  )

p_state <- ggplot(plot_b, aes(x = dynamic, y = mobile)) +
  geom_vline(
    xintercept = 0,
    color = "grey60",
    linewidth = 0.5 * FIG9_ELEMENT_SCALE,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey60",
    linewidth = 0.5 * FIG9_ELEMENT_SCALE,
    linetype = "dashed"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = alpha("grey30", 0.70),
    linetype = "dashed",
    linewidth = 0.75 * FIG9_ELEMENT_SCALE
  ) +
  geom_point(
    aes(fill = pc1, size = pc2),
    shape = 21,
    color = "grey35",
    stroke = 0.35 * FIG9_ELEMENT_SCALE,
    alpha = 0.96
  ) +
  geom_text_repel(
    data = plot_b_labels,
    aes(label = site),
    seed = 20260228,
    nudge_x = plot_b_labels$nudge_x,
    nudge_y = plot_b_labels$nudge_y,
    size = site_label_size,
    box.padding = state_label_box_padding,
    point.padding = state_label_point_padding,
    force = state_label_force,
    force_pull = state_label_force_pull,
    min.segment.length = Inf,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  geom_label(
    data = quad_df,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE,
    size = quad_label_size,
    color = "grey25",
    fill = alpha("grey95", 0.92),
    linewidth = 0,
    label.padding = unit(0.25, "lines")
  ) +
  annotate(
    "text",
    x = stats_annot_x,
    y = stats_annot_y,
    label = ann_stats,
    parse = TRUE,
    hjust = 0.5,
    vjust = 1,
    size = stats_label_size,
    color = "grey20"
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
    range = c(3.5, 9) * FIG9_ELEMENT_SCALE
  ) +
  labs(
    x = "Dynamic Storage",
    y = "Mobile Storage"
  ) +
  coord_cartesian(xlim = x_lim, ylim = y_lim, clip = "off") +
  theme_pub() +
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    axis.title.x = element_text(size = axis_title_size, margin = margin(t = 10, b = 10)),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size),
    legend.box = "vertical",
    legend.spacing.y = unit(1.10, "cm"),
    legend.box.spacing = unit(0.80, "cm"),
    legend.key.height = unit(0.80, "cm"),
    legend.key.width = unit(0.72, "cm"),
    legend.margin = margin(8, 12, 8, 12),
    legend.box.margin = margin(2, 2, 2, 2),
    legend.position = c(1.03, 0.63),
    legend.justification = c(0, 0.5),
    plot.margin = margin(24, 220, 88, 5.5)
  ) +
  guides(
    fill = guide_colorbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      title.theme = element_text(margin = margin(b = 14)),
      barheight = unit(7.6, "cm"),
      barwidth = unit(0.90, "cm"),
      ticks.colour = "grey35",
      frame.colour = "grey35"
    ),
    size = guide_legend(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      title.theme = element_text(margin = margin(b = 14)),
      byrow = TRUE,
      keyheight = unit(0.92, "cm"),
      keywidth = unit(0.92, "cm"),
      override.aes = list(fill = "grey82", color = "grey35", alpha = 1)
    )
  )

fig9 <- p_geo / p_state +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "plain", size = tag_label_size))

fig9_width <- 12.8 * FIG_WIDTH_SCALE * FIG9_EXPORT_SCALE
fig9_height <- 16.4 * FIG_HEIGHT_SCALE * FIG9_EXPORT_SCALE

nm <- "Fig9_conceptual_diagram"
invisible(safe_ggsave(
  file.path(main_dir, paste0(nm, ".png")),
  fig9,
  width = fig9_width,
  height = fig9_height,
  dpi = 300
))
invisible(safe_ggsave(
  file.path(main_pdf_dir, paste0(nm, ".pdf")),
  fig9,
  width = fig9_width,
  height = fig9_height
))
invisible(safe_ggsave(
  file.path(explore_dir, paste0(nm, ".png")),
  fig9,
  width = fig9_width,
  height = fig9_height,
  dpi = 300
))
