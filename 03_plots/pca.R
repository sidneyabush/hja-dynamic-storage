# -----------------------------------------------------------------------------
# PCA Plots
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(scales)

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
if (!file.exists(config_path)) config_path <- file.path(dirname(script_dir), "config.R")
if (!file.exists(config_path)) config_path <- file.path(getwd(), "config.R")
source(config_path)

main_dir <- file.path(FIGURES_DIR, "main")
supp_dir <- file.path(FIGURES_DIR, "supp", "stats", "pca")
table_dir <- file.path(OUT_TABLES_DIR, "pca")
for (d in c(main_dir, supp_dir, table_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

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
    pretty_label = dplyr::recode(
      feature,
      "RCS" = "Recession Slope",
      "RBI" = "Flashiness (RBI)",
      "FDC" = "FDC Slope",
      "SD" = "Storage-Discharge",
      "WB" = "Water Balance",
      .default = feature
    )
  )

p_biplot <- ggplot(scores, aes(x = PC1, y = PC2, color = site)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_point(size = FIG_POINT_SIZE_SMALL, alpha = 0.7) +
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
    nudge_y = 0.05
  ) +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = paste0("PC1 (", number(pc1_pct, accuracy = 0.1), "%)"),
    y = paste0("PC2 (", number(pc2_pct, accuracy = 0.1), "%)")
  ) +
  theme_classic(base_size = FIG_BASE_SIZE) +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line = element_blank()
  )

ggsave(file.path(main_dir, "pca_biplot.png"), p_biplot, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(main_dir, "pca_biplot.pdf"), p_biplot, width = 8.5, height = 6.5)

p_scree <- vexp %>%
  mutate(PC = factor(PC, levels = PC)) %>%
  ggplot(aes(x = PC, y = Variance_Explained)) +
  geom_col(fill = "grey55", color = "black", linewidth = 0.3) +
  geom_text(aes(label = percent(Variance_Explained, accuracy = 0.1)), vjust = -0.3, size = FIG_ANNOT_TEXT_SIZE) +
  labs(x = NULL, y = "Variance explained") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
  theme_classic(base_size = FIG_BASE_SIZE) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line = element_blank()
  )

ggsave(file.path(supp_dir, "pca_scree.png"), p_scree, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(supp_dir, "pca_scree.pdf"), p_scree, width = 7, height = 4.5)

write_csv(vexp, file.path(table_dir, "pca_variance_table.csv"))
write_csv(loads, file.path(table_dir, "pca_loading_table.csv"))

# Eco-model PCA plots (if eco PCA outputs exist)
eco_scores_file <- file.path(OUT_STATS_PCA_DIR, "eco_mlr_pca_scores_pc1_pc2.csv")
eco_loads_file <- file.path(OUT_STATS_PCA_DIR, "eco_mlr_pca_loadings.csv")
eco_var_file <- file.path(OUT_STATS_PCA_DIR, "eco_mlr_pca_variance_explained.csv")

if (file.exists(eco_scores_file) && file.exists(eco_loads_file) && file.exists(eco_var_file)) {
  eco_scores <- read_csv(eco_scores_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

  eco_loads <- read_csv(eco_loads_file, show_col_types = FALSE)
  eco_vexp <- read_csv(eco_var_file, show_col_types = FALSE)

  eco_pc1_pct <- ifelse(nrow(eco_vexp) >= 1, 100 * eco_vexp$Variance_Explained[1], NA_real_)
  eco_pc2_pct <- ifelse(nrow(eco_vexp) >= 2, 100 * eco_vexp$Variance_Explained[2], NA_real_)

  eco_score_lim <- max(abs(c(eco_scores$PC1, eco_scores$PC2)), na.rm = TRUE)
  eco_load_lim <- max(abs(c(eco_loads$PC1, eco_loads$PC2)), na.rm = TRUE)
  eco_arrow_scale <- ifelse(
    is.finite(eco_score_lim) && is.finite(eco_load_lim) && eco_load_lim > 0,
    0.75 * eco_score_lim / eco_load_lim,
    1
  )

  eco_loads_plot <- eco_loads %>%
    mutate(
      PC1s = PC1 * eco_arrow_scale,
      PC2s = PC2 * eco_arrow_scale,
      pretty_label = dplyr::recode(
        feature,
        "RCS" = "Recession Slope",
        "RBI" = "Flashiness (RBI)",
        "FDC" = "FDC Slope",
        "SD" = "Storage-Discharge",
        "WB" = "Water Balance",
        "CHS" = "CHS",
        "MTT" = "MTT",
        "Fyw" = "Fyw",
        "DR" = "DR",
        "Lava1_per" = "Lava-1",
        "Lava2_per" = "Lava-2",
        "Ash_Per" = "Ash",
        "Landslide_Total" = "LS Total",
        "Landslide_Young" = "LS Young",
        .default = feature
      )
    )

  p_eco_biplot <- ggplot(eco_scores, aes(x = PC1, y = PC2, color = site)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_point(size = FIG_POINT_SIZE_SMALL, alpha = 0.7) +
    geom_segment(
      data = eco_loads_plot,
      aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.18, "cm")),
      color = "black",
      linewidth = 0.5
    ) +
    geom_label(
      data = eco_loads_plot,
      aes(x = PC1s, y = PC2s, label = pretty_label),
      inherit.aes = FALSE,
      color = "black",
      fill = scales::alpha("white", 0.75),
      label.size = 0,
      label.padding = grid::unit(0.12, "lines"),
      size = FIG_ANNOT_TEXT_SIZE,
      nudge_x = 0.05,
      nudge_y = 0.05
    ) +
    scale_color_manual(values = SITE_COLORS) +
    labs(
      x = paste0("PC1 (", number(eco_pc1_pct, accuracy = 0.1), "%)"),
      y = paste0("PC2 (", number(eco_pc2_pct, accuracy = 0.1), "%)")
    ) +
    theme_classic(base_size = FIG_BASE_SIZE) +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.line = element_blank()
    )

  ggsave(file.path(supp_dir, "eco_mlr_pca_biplot.png"), p_eco_biplot, width = 8.5, height = 6.5, dpi = 300)
  ggsave(file.path(supp_dir, "eco_mlr_pca_biplot.pdf"), p_eco_biplot, width = 8.5, height = 6.5)

  p_eco_scree <- eco_vexp %>%
    mutate(PC = factor(PC, levels = PC)) %>%
    ggplot(aes(x = PC, y = Variance_Explained)) +
    geom_col(fill = "grey55", color = "black", linewidth = 0.3) +
    geom_text(aes(label = percent(Variance_Explained, accuracy = 0.1)), vjust = -0.3, size = FIG_ANNOT_TEXT_SIZE) +
    labs(x = NULL, y = "Variance explained") +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
    theme_classic(base_size = FIG_BASE_SIZE) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.line = element_blank()
    )

  ggsave(file.path(supp_dir, "eco_mlr_pca_scree.png"), p_eco_scree, width = 7, height = 4.5, dpi = 300)
  ggsave(file.path(supp_dir, "eco_mlr_pca_scree.pdf"), p_eco_scree, width = 7, height = 4.5)

  write_csv(eco_vexp, file.path(table_dir, "eco_mlr_pca_variance_table.csv"))
  write_csv(eco_loads, file.path(table_dir, "eco_mlr_pca_loading_table.csv"))
}

# Catchment-characteristics PCA plots (if catchment PCA outputs exist)
catch_scores_file <- file.path(OUT_STATS_PCA_DIR, "catch_chars_mlr_pca_scores_pc1_pc2.csv")
catch_loads_file <- file.path(OUT_STATS_PCA_DIR, "catch_chars_mlr_pca_loadings.csv")
catch_var_file <- file.path(OUT_STATS_PCA_DIR, "catch_chars_mlr_pca_variance_explained.csv")

if (file.exists(catch_scores_file) && file.exists(catch_loads_file) && file.exists(catch_var_file)) {
  catch_scores <- read_csv(catch_scores_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

  catch_loads <- read_csv(catch_loads_file, show_col_types = FALSE)
  catch_vexp <- read_csv(catch_var_file, show_col_types = FALSE)

  catch_pc1_pct <- ifelse(nrow(catch_vexp) >= 1, 100 * catch_vexp$Variance_Explained[1], NA_real_)
  catch_pc2_pct <- ifelse(nrow(catch_vexp) >= 2, 100 * catch_vexp$Variance_Explained[2], NA_real_)

  catch_score_lim <- max(abs(c(catch_scores$PC1, catch_scores$PC2)), na.rm = TRUE)
  catch_load_lim <- max(abs(c(catch_loads$PC1, catch_loads$PC2)), na.rm = TRUE)
  catch_arrow_scale <- ifelse(
    is.finite(catch_score_lim) && is.finite(catch_load_lim) && catch_load_lim > 0,
    0.75 * catch_score_lim / catch_load_lim,
    1
  )

  catch_loads_plot <- catch_loads %>%
    mutate(
      PC1s = PC1 * catch_arrow_scale,
      PC2s = PC2 * catch_arrow_scale,
      pretty_label = dplyr::recode(
        feature,
        "Slope_mean" = "Slope",
        "Landslide_Total" = "LS Total",
        "Landslide_Young" = "LS Young",
        "Lava1_per" = "Lava-1",
        "Lava2_per" = "Lava-2",
        "Ash_Per" = "Ash",
        "Pyro_per" = "Pyroclastic",
        .default = feature
      )
    )

  p_catch_biplot <- ggplot(catch_scores, aes(x = PC1, y = PC2, color = site)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_point(size = FIG_POINT_SIZE_SMALL, alpha = 0.7) +
    geom_segment(
      data = catch_loads_plot,
      aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.18, "cm")),
      color = "black",
      linewidth = 0.5
    ) +
    geom_label(
      data = catch_loads_plot,
      aes(x = PC1s, y = PC2s, label = pretty_label),
      inherit.aes = FALSE,
      color = "black",
      fill = scales::alpha("white", 0.75),
      label.size = 0,
      label.padding = grid::unit(0.12, "lines"),
      size = FIG_ANNOT_TEXT_SIZE,
      nudge_x = 0.05,
      nudge_y = 0.05
    ) +
    scale_color_manual(values = SITE_COLORS) +
    labs(
      x = paste0("PC1 (", number(catch_pc1_pct, accuracy = 0.1), "%)"),
      y = paste0("PC2 (", number(catch_pc2_pct, accuracy = 0.1), "%)")
    ) +
    theme_classic(base_size = FIG_BASE_SIZE) +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.line = element_blank()
    )

  ggsave(file.path(supp_dir, "catch_chars_mlr_pca_biplot.png"), p_catch_biplot, width = 8.5, height = 6.5, dpi = 300)
  ggsave(file.path(supp_dir, "catch_chars_mlr_pca_biplot.pdf"), p_catch_biplot, width = 8.5, height = 6.5)

  p_catch_scree <- catch_vexp %>%
    mutate(PC = factor(PC, levels = PC)) %>%
    ggplot(aes(x = PC, y = Variance_Explained)) +
    geom_col(fill = "grey55", color = "black", linewidth = 0.3) +
    geom_text(aes(label = percent(Variance_Explained, accuracy = 0.1)), vjust = -0.3, size = FIG_ANNOT_TEXT_SIZE) +
    labs(x = NULL, y = "Variance explained") +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
    theme_classic(base_size = FIG_BASE_SIZE) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.line = element_blank()
    )

  ggsave(file.path(supp_dir, "catch_chars_mlr_pca_scree.png"), p_catch_scree, width = 7, height = 4.5, dpi = 300)
  ggsave(file.path(supp_dir, "catch_chars_mlr_pca_scree.pdf"), p_catch_scree, width = 7, height = 4.5)

  write_csv(catch_vexp, file.path(table_dir, "catch_chars_mlr_pca_variance_table.csv"))
  write_csv(catch_loads, file.path(table_dir, "catch_chars_mlr_pca_loading_table.csv"))
}
