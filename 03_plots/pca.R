# PCA Plots.
# Inputs: OUT_STATS_PCA_DIR/pca_scores_pc1_pc2.csv; OUT_STATS_PCA_DIR/pca_loadings.csv; OUT_STATS_PCA_DIR/pca_variance_explained.csv.
# Author: Sidney Bush
# Date: 2026-02-13

library(dplyr)
library(readr)
library(ggplot2)
library(scales)

rm(list = ls())

# Load project config
source("config.R")


main_dir <- file.path(FIGURES_DIR, "main")
supp_dir <- file.path(FIGURES_DIR, "supp", "analysis", "pca")
table_dir <- file.path(OUT_TABLES_DIR, "pca")
for (d in c(main_dir, supp_dir, table_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Remove older catchment-characteristics PCA outputs.
unlink(file.path(supp_dir, c(
  "watershed_char_storage_mlr_pca_biplot.png",
  "watershed_char_storage_mlr_pca_biplot.pdf",
  "watershed_char_storage_mlr_pca_scree.png",
  "watershed_char_storage_mlr_pca_scree.pdf"
)))
unlink(file.path(supp_dir, c(
  "storage_ecovar_mlr_pca_biplot.png",
  "storage_ecovar_mlr_pca_biplot.pdf",
  "storage_ecovar_mlr_pca_scree.png",
  "storage_ecovar_mlr_pca_scree.pdf"
)))
unlink(file.path(table_dir, c(
  "watershed_char_storage_mlr_pca_variance_table.csv",
  "watershed_char_storage_mlr_pca_loading_table.csv"
)))
unlink(file.path(table_dir, c(
  "storage_ecovar_mlr_pca_variance_table.csv",
  "storage_ecovar_mlr_pca_loading_table.csv"
)))

scores_file <- file.path(OUT_STATS_PCA_DIR, "pca_scores_pc1_pc2.csv")
loads_file <- file.path(OUT_STATS_PCA_DIR, "pca_loadings.csv")
var_file <- file.path(OUT_STATS_PCA_DIR, "pca_variance_explained.csv")

scores <- read_csv(scores_file, show_col_types = FALSE) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

loads <- read_csv(loads_file, show_col_types = FALSE)
vexp <- read_csv(var_file, show_col_types = FALSE)

pc1_pct <- ifelse(nrow(vexp) >= 1, 100 * vexp$Variance_Explained[1], NA_real_)
pc2_pct <- ifelse(nrow(vexp) >= 2, 100 * vexp$Variance_Explained[2], NA_real_)

# Scale loading arrows to score space.
score_lim <- max(abs(c(scores$PC1, scores$PC2)), na.rm = TRUE)
load_lim <- max(abs(c(loads$PC1, loads$PC2)), na.rm = TRUE)
arrow_scale <- ifelse(is.finite(score_lim) && is.finite(load_lim) && load_lim > 0, 0.75 * score_lim / load_lim, 1)

loads_plot <- loads %>%
  mutate(
    PC1s = PC1 * arrow_scale,
    PC2s = PC2 * arrow_scale,
    pretty_label = toupper(feature)
  )

p_biplot <- ggplot(scores, aes(x = PC1, y = PC2, color = site)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_point(size = FIG_POINT_SIZE_LARGE + 1, alpha = 0.85) +
  geom_segment(
    data = loads_plot,
    aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.18, "cm")),
    color = "black",
    linewidth = 0.5
  ) +
  geom_label(
    data = loads_plot,
    aes(x = PC1s, y = PC2s, label = pretty_label),
    inherit.aes = FALSE,
    color = "black",
    fill = scales::alpha("white", 0.75),
    label.size = 0,
    label.padding = grid::unit(0.12, "lines"),
    size = FIG_ANNOT_TEXT_SIZE,
    nudge_x = 0.05,
    nudge_y = 0.05,
    check_overlap = FIG_LABEL_CHECK_OVERLAP
  ) +
  scale_color_manual(values = SITE_COLORS) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = FIG_POINT_SIZE_LARGE + 1, alpha = 1))) +
  labs(
    x = paste0("PC1 (", number(pc1_pct, accuracy = 0.1), "%)"),
    y = paste0("PC2 (", number(pc2_pct, accuracy = 0.1), "%)")
  ) +
  theme_pub() +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line = element_blank(),
    plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT)
  ) +
  coord_cartesian(clip = FIG_LABEL_CLIP)

ggsave(file.path(main_dir, "pca_biplot.png"), p_biplot, width = 8.5 * FIG_WIDTH_SCALE, height = 6.5 * FIG_HEIGHT_SCALE, dpi = 300)
ggsave(file.path(main_dir, "pca_biplot.pdf"), p_biplot, width = 8.5 * FIG_WIDTH_SCALE, height = 6.5 * FIG_HEIGHT_SCALE)

p_scree <- vexp %>%
  mutate(PC = factor(PC, levels = PC)) %>%
  ggplot(aes(x = PC, y = Variance_Explained)) +
  geom_col(fill = "grey55", color = "black", linewidth = 0.3) +
  geom_text(aes(label = percent(Variance_Explained, accuracy = 0.1)), vjust = -0.3, size = FIG_ANNOT_TEXT_SIZE) +
  labs(x = NULL, y = "Variance explained") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line = element_blank()
  )

ggsave(file.path(supp_dir, "pca_scree.png"), p_scree, width = 7 * FIG_WIDTH_SCALE, height = 4.5 * FIG_HEIGHT_SCALE, dpi = 300)
ggsave(file.path(supp_dir, "pca_scree.pdf"), p_scree, width = 7 * FIG_WIDTH_SCALE, height = 4.5 * FIG_HEIGHT_SCALE)

write_csv(vexp, file.path(table_dir, "pca_variance_table.csv"))
write_csv(loads, file.path(table_dir, "pca_loading_table.csv"))

# Catchment-characteristics and eco-model PCA companion plots removed by design.
