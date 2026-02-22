# Plot unified subsurface storage framework outputs.
# Inputs: OUT_STATS_DIR/unified_framework/unified_framework_site_axes.csv.
# Author: Sidney Bush
# Date: 2026-02-21

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
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
    if (isTRUE(ok)) {
      return(r)
    }
  }
  stop("No writable output root found for unified framework plots.")
}

find_axes_file <- function() {
  candidates <- c(
    file.path(OUTPUT_DIR, "models", "unified_framework", "unified_framework_site_axes.csv"),
    file.path(REPO_DIR, "outputs", "models", "unified_framework", "unified_framework_site_axes.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Missing required unified-framework file: ", candidates[1])
  }
  hit[1]
}

write_root <- find_write_root()
main_fig_dir <- file.path(write_root, "figs", "main")
main_pdf_dir <- file.path(main_fig_dir, "pdf")
for (d in c(main_fig_dir, main_pdf_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

axes_file <- find_axes_file()

axes <- read_csv(axes_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site)

scatter_df <- axes %>%
  filter(is.finite(dynamic_storage_strength_z), is.finite(mobile_mixing_z))

p_state <- ggplot(
  scatter_df,
  aes(
    x = dynamic_storage_strength_z,
    y = mobile_mixing_z,
    fill = flow_path_partitioning_z,
    label = site
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_point(
    shape = 21,
    size = FIG_POINT_SIZE_LARGE + 0.8,
    color = "grey20",
    stroke = 0.35
  ) +
  geom_text(
    nudge_y = 0.08,
    size = FIG_ANNOT_TEXT_SIZE,
    check_overlap = FIG_LABEL_CHECK_OVERLAP,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(
    low = "#b2182b",
    mid = "white",
    high = "#2166ac",
    midpoint = 0,
    na.value = "grey70",
    name = "Flow-path\npartitioning (z)"
  ) +
  labs(
    x = "Dynamic storage strength axis (z)",
    y = "Mobile turnover/mixing\naxis (z)"
  ) +
  theme_pub() +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    legend.title = element_text(size = FIG_AXIS_TITLE_SIZE - 1),
    legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
  )

ggsave(
  file.path(main_fig_dir, "unified_framework_state_space.png"),
  p_state,
  width = 8 * FIG_WIDTH_SCALE,
  height = 6 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "unified_framework_state_space.pdf"),
  p_state,
  width = 8 * FIG_WIDTH_SCALE,
  height = 6 * FIG_HEIGHT_SCALE
)

axis_long <- axes %>%
  select(
    site,
    dynamic_storage_strength_z,
    mobile_mixing_z,
    flow_path_partitioning_z,
    unified_state_index
  ) %>%
  pivot_longer(
    cols = -site,
    names_to = "axis",
    values_to = "score"
  ) %>%
  mutate(
    axis = factor(
      axis,
      levels = c(
        "dynamic_storage_strength_z",
        "mobile_mixing_z",
        "flow_path_partitioning_z",
        "unified_state_index"
      ),
      labels = c(
        "Dynamic\nStorage Strength",
        "Mobile Mixing",
        "Flow-Path Partitioning",
        "Unified State Index"
      )
    )
  )

p_heat <- ggplot(axis_long, aes(x = axis, y = site, fill = score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(is.finite(score), sprintf("%.2f", score), "NA")),
    size = FIG_TILE_TEXT_SIZE
  ) +
  scale_y_discrete(limits = rev(SITE_ORDER_HYDROMETRIC)) +
  scale_fill_gradient2(
    low = "#b2182b",
    mid = "white",
    high = "#2166ac",
    midpoint = 0,
    na.value = "grey85",
    name = "z-score"
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(
    axis.text.x = element_text(
      size = FIG_AXIS_TEXT_SIZE,
      angle = 0,
      hjust = 0.5,
      lineheight = 0.95
    ),
    axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
    legend.title = element_text(size = FIG_AXIS_TITLE_SIZE - 1),
    legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
  )

ggsave(
  file.path(main_fig_dir, "unified_framework_axes_heatmap.png"),
  p_heat,
  width = 9 * FIG_WIDTH_SCALE,
  height = 5.5 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(main_pdf_dir, "unified_framework_axes_heatmap.pdf"),
  p_heat,
  width = 9 * FIG_WIDTH_SCALE,
  height = 5.5 * FIG_HEIGHT_SCALE
)
