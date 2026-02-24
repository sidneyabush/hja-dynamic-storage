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

find_master_annual_file <- function() {
  candidates <- c(
    file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE),
    file.path(REPO_DIR, "outputs", "master", MASTER_ANNUAL_FILE)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Missing required master annual file: ", candidates[1])
  }
  hit[1]
}

find_pca_scores_file <- function() {
  candidates <- c(
    file.path(OUTPUT_DIR, "models", "pca", "pca_scores_pc1_pc2.csv"),
    file.path(REPO_DIR, "outputs", "models", "pca", "pca_scores_pc1_pc2.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Missing required PCA scores file: ", candidates[1])
  }
  hit[1]
}

write_root <- find_write_root()
model_out_dir <- file.path(write_root, "models", "unified_framework")
main_fig_dir <- file.path(write_root, "figs", "main")
main_pdf_dir <- file.path(main_fig_dir, "pdf")
dynamic_v_mobile_fig_dir <- file.path(main_fig_dir, "dynamic_v_mobile")
dynamic_v_mobile_pdf_dir <- file.path(main_pdf_dir, "dynamic_v_mobile")
supp_fig_dir <- file.path(write_root, "figs", "supp")
supp_pdf_dir <- file.path(supp_fig_dir, "pdf")
for (d in c(
  model_out_dir,
  main_fig_dir,
  main_pdf_dir,
  dynamic_v_mobile_fig_dir,
  dynamic_v_mobile_pdf_dir,
  supp_fig_dir,
  supp_pdf_dir
)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

axes_file <- find_axes_file()
annual_file <- find_master_annual_file()
pca_scores_file <- find_pca_scores_file()

axes <- read_csv(axes_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site)

if (!("geology_class" %in% names(axes))) {
  axes$geology_class <- NA_character_
}
if (!("geomorphology_class" %in% names(axes))) {
  axes$geomorphology_class <- NA_character_
}
if (!("mobile_mixing_with_chs_z" %in% names(axes))) {
  axes$mobile_mixing_with_chs_z <- suppressWarnings(as.numeric(axes$mobile_mixing_z))
}
if (!("mobile_mixing_no_chs_z" %in% names(axes))) {
  axes$mobile_mixing_no_chs_z <- NA_real_
}
if (!("CHS_mean" %in% names(axes))) {
  axes$CHS_mean <- NA_real_
}
if (!("Landslide_Total" %in% names(axes))) {
  axes$Landslide_Total <- NA_real_
}
if (!("Lava1_per" %in% names(axes))) {
  axes$Lava1_per <- NA_real_
}
if (!("Lava2_per" %in% names(axes))) {
  axes$Lava2_per <- NA_real_
}
if (!("Ash_Per" %in% names(axes))) {
  axes$Ash_Per <- NA_real_
}
if (!("Pyro_per" %in% names(axes))) {
  axes$Pyro_per <- NA_real_
}
if (!("geology_pc1" %in% names(axes))) {
  axes$geology_pc1 <- NA_real_
}
if (!("geology_pc2" %in% names(axes))) {
  axes$geology_pc2 <- NA_real_
}

axes <- axes %>%
  mutate(
    geology_class = ifelse(is.na(geology_class) | !nzchar(geology_class), "Unknown", geology_class),
    geomorphology_class = ifelse(
      is.na(geomorphology_class) | !nzchar(geomorphology_class),
      "Unknown",
      geomorphology_class
    ),
    geology_class = factor(
      geology_class,
      levels = c("Pyroclastic-dominant", "Lava-dominant", "Ash-dominant", "Mixed", "Unknown")
    ),
    geomorphology_class = factor(
      geomorphology_class,
      levels = c("Harvest-dominated", "Landslide-dominated", "Harvest+landslide", "Low-disturbance", "Unknown")
    )
  )

scatter_df_all <- axes %>%
  transmute(
    site = as.character(site),
    dynamic_storage_strength_z = as.numeric(dynamic_storage_strength_z),
    mobile_mixing_with_chs_z = as.numeric(dplyr::coalesce(mobile_mixing_with_chs_z, mobile_mixing_z)),
    mobile_mixing_no_chs_z = as.numeric(mobile_mixing_no_chs_z),
    CHS_mean = as.numeric(CHS_mean),
    Landslide_Total = as.numeric(Landslide_Total),
    Lava1_per = as.numeric(Lava1_per),
    Lava2_per = as.numeric(Lava2_per),
    Ash_Per = as.numeric(Ash_Per),
    Pyro_per = as.numeric(Pyro_per),
    geology_pc1 = as.numeric(geology_pc1),
    geology_pc2 = as.numeric(geology_pc2)
  )

format_p_value <- function(p_val) {
  if (!is.finite(p_val)) {
    return("NA")
  }
  if (p_val < 0.001) {
    return("<0.001")
  }
  sprintf("%.3f", p_val)
}

fit_stats <- function(data_in, x_col, y_col) {
  fit_df <- data_in %>%
    filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))

  if (nrow(fit_df) < 3) {
    return(list(
      slope = NA_real_,
      r2 = NA_real_,
      p_value = NA_real_,
      label = "R2 = NA\np = NA\nb = NA"
    ))
  }

  fit <- lm(reformulate(x_col, y_col), data = fit_df)
  fit_sum <- summary(fit)
  slope <- unname(coef(fit)[2])
  r2 <- unname(fit_sum$r.squared)
  p_val <- unname(coef(fit_sum)[2, 4])

  list(
    slope = slope,
    r2 = r2,
    p_value = p_val,
    label = sprintf(
      "R2 = %.3f\np = %s\nb = %.3f",
      r2,
      format_p_value(p_val),
      slope
    )
  )
}

build_quadrant_annotations <- function(x_rng, y_rng, include_chs = TRUE) {
  x_min <- min(c(x_rng[1], 0), na.rm = TRUE)
  x_max <- max(c(x_rng[2], 0), na.rm = TRUE)
  y_min <- min(c(y_rng[1], 0), na.rm = TRUE)
  y_max <- max(c(y_rng[2], 0), na.rm = TRUE)
  y_span <- ifelse(is.finite(y_max - y_min) && (y_max - y_min) > 0, y_max - y_min, 1)

  x_left <- (x_min + 0) / 2
  x_right <- (x_max + 0) / 2
  y_top <- y_max - 0.06 * y_span
  y_bottom <- y_min + 0.06 * y_span

  dyn_high <- "RBI-, RCS/SD+, FDC flatter, WBdep+"
  dyn_low <- "RBI+, RCS/SD-, FDC steeper, WBdep-"
  mobile_high <- if (isTRUE(include_chs)) {
    "MTT+, DR-, Fyw-, CHS+"
  } else {
    "MTT+, DR-, Fyw-"
  }
  mobile_low <- if (isTRUE(include_chs)) {
    "MTT-, DR+, Fyw+, CHS-"
  } else {
    "MTT-, DR+, Fyw+"
  }

  tibble(
    x = c(x_left, x_right, x_left, x_right),
    y = c(y_top, y_top, y_bottom, y_bottom),
    label = c(
      paste(dyn_low, mobile_high, sep = "\n"),
      paste(dyn_high, mobile_high, sep = "\n"),
      paste(dyn_low, mobile_low, sep = "\n"),
      paste(dyn_high, mobile_low, sep = "\n")
    )
  )
}

make_state_plot <- function(
  data_in,
  mobile_col,
  color_col,
  color_label,
  panel_title,
  palette_type,
  plot_name,
  mobile_axis_definition,
  require_chs = FALSE
) {
  plot_df <- data_in %>%
    filter(is.finite(dynamic_storage_strength_z), is.finite(.data[[mobile_col]]))
  if (isTRUE(require_chs)) {
    plot_df <- plot_df %>% filter(is.finite(CHS_mean))
  }

  if (nrow(plot_df) == 0) {
    return(list(
      plot = ggplot() + theme_void(),
      stats = tibble(
        plot_name = plot_name,
        mobile_axis_definition = mobile_axis_definition,
        color_variable = color_label,
        slope = NA_real_,
        r2 = NA_real_,
        p_value = NA_real_,
        n = 0L
      )
    ))
  }

  x_rng <- range(plot_df$dynamic_storage_strength_z, na.rm = TRUE)
  y_rng <- range(plot_df[[mobile_col]], na.rm = TRUE)
  x_span <- ifelse(is.finite(diff(x_rng)) && diff(x_rng) > 0, diff(x_rng), 1)
  y_span <- ifelse(is.finite(diff(y_rng)) && diff(y_rng) > 0, diff(y_rng), 1)
  ann_x <- x_rng[2] - 0.03 * x_span
  ann_y <- y_rng[2] - 0.22 * y_span
  fit_out <- fit_stats(plot_df, "dynamic_storage_strength_z", mobile_col)
  ann_txt <- fit_out$label
  quad_df <- build_quadrant_annotations(
    x_rng = x_rng,
    y_rng = y_rng,
    include_chs = grepl("^includes_CHS", mobile_axis_definition)
  )
  quad_label_size <- max(FIG_ANNOT_TEXT_SIZE - 1, 2.8)

  p <- ggplot(
    plot_df,
    aes(
      x = dynamic_storage_strength_z,
      y = .data[[mobile_col]],
      fill = .data[[color_col]]
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey75") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey75") +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.35,
      size = FIG_POINT_SIZE_LARGE + 0.6,
      alpha = 0.95
    ) +
    geom_text(
      aes(label = site),
      color = "black",
      nudge_y = 0.05,
      size = FIG_ANNOT_TEXT_SIZE,
      show.legend = FALSE
    ) +
    geom_smooth(
      data = plot_df,
      aes(x = dynamic_storage_strength_z, y = .data[[mobile_col]]),
      inherit.aes = FALSE,
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    geom_label(
      data = quad_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = quad_label_size,
      color = "grey20",
      lineheight = 0.9,
      fill = grDevices::adjustcolor("white", alpha.f = 0.72),
      linewidth = 0.15,
      label.padding = grid::unit(0.10, "lines"),
      label.r = grid::unit(0.08, "lines")
    ) +
    annotate(
      "label",
      x = ann_x,
      y = ann_y,
      label = ann_txt,
      hjust = 1,
      vjust = 1,
      size = FIG_ANNOT_TEXT_SIZE,
      linewidth = 0.2,
      fill = "white",
      alpha = 0.9
    ) +
    labs(
      title = panel_title,
      x = "Dynamic Storage",
      y = "Mobile Storage",
      fill = color_label
    ) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      plot.title = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE - 1),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
    )

  if (identical(palette_type, "chs")) {
    p <- p + scale_fill_gradient2(
      low = "#d73027",
      mid = "white",
      high = "#2166ac",
      midpoint = 0.5,
      na.value = "#bdbdbd"
    )
  } else if (identical(palette_type, "geology_pc")) {
    p <- p + scale_fill_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      na.value = "#bdbdbd"
    )
  } else if (identical(palette_type, "landslide")) {
    p <- p + scale_fill_gradient(
      low = "#fff7bc",
      high = "#b10026",
      na.value = "#bdbdbd"
    )
  }

  list(
    plot = p,
    stats = tibble(
      plot_name = plot_name,
      mobile_axis_definition = mobile_axis_definition,
      color_variable = color_label,
      slope = fit_out$slope,
      r2 = fit_out$r2,
      p_value = fit_out$p_value,
      n = nrow(plot_df)
    )
  )
}

make_geology_pc_bivariate_plot <- function(
  data_in,
  mobile_col,
  panel_title,
  plot_name,
  mobile_axis_definition,
  require_chs = FALSE
) {
  plot_df <- data_in %>%
    filter(
      is.finite(dynamic_storage_strength_z),
      is.finite(.data[[mobile_col]]),
      is.finite(geology_pc1),
      is.finite(geology_pc2)
    )
  if (isTRUE(require_chs)) {
    plot_df <- plot_df %>% filter(is.finite(CHS_mean))
  }

  if (nrow(plot_df) == 0) {
    return(list(
      plot = ggplot() + theme_void(),
      stats = tibble(
        plot_name = plot_name,
        mobile_axis_definition = mobile_axis_definition,
        color_variable = "Geology PC1 + Geology PC2",
        slope = NA_real_,
        r2 = NA_real_,
        p_value = NA_real_,
        n = 0L
      )
    ))
  }

  x_rng <- range(plot_df$dynamic_storage_strength_z, na.rm = TRUE)
  y_rng <- range(plot_df[[mobile_col]], na.rm = TRUE)
  x_span <- ifelse(is.finite(diff(x_rng)) && diff(x_rng) > 0, diff(x_rng), 1)
  y_span <- ifelse(is.finite(diff(y_rng)) && diff(y_rng) > 0, diff(y_rng), 1)
  ann_x <- x_rng[2] - 0.03 * x_span
  ann_y <- y_rng[2] - 0.22 * y_span
  fit_out <- fit_stats(plot_df, "dynamic_storage_strength_z", mobile_col)
  quad_df <- build_quadrant_annotations(
    x_rng = x_rng,
    y_rng = y_rng,
    include_chs = grepl("^includes_CHS", mobile_axis_definition)
  )
  quad_label_size <- max(FIG_ANNOT_TEXT_SIZE - 1, 2.8)

  p <- ggplot(
    plot_df,
    aes(x = dynamic_storage_strength_z, y = .data[[mobile_col]])
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey75") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey75") +
    geom_point(
      aes(fill = geology_pc1, color = geology_pc2),
      shape = 21,
      stroke = 0.9,
      size = FIG_POINT_SIZE_LARGE + 1.0,
      alpha = 0.95
    ) +
    geom_text(
      aes(label = site),
      color = "black",
      nudge_y = 0.05,
      size = FIG_ANNOT_TEXT_SIZE,
      show.legend = FALSE
    ) +
    geom_smooth(
      data = plot_df,
      aes(x = dynamic_storage_strength_z, y = .data[[mobile_col]]),
      inherit.aes = FALSE,
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    geom_label(
      data = quad_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = quad_label_size,
      color = "grey20",
      lineheight = 0.9,
      fill = grDevices::adjustcolor("white", alpha.f = 0.72),
      linewidth = 0.15,
      label.padding = grid::unit(0.10, "lines"),
      label.r = grid::unit(0.08, "lines")
    ) +
    annotate(
      "label",
      x = ann_x,
      y = ann_y,
      label = fit_out$label,
      hjust = 1,
      vjust = 1,
      size = FIG_ANNOT_TEXT_SIZE,
      linewidth = 0.2,
      fill = "white",
      alpha = 0.9
    ) +
    scale_fill_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      na.value = "#bdbdbd",
      name = "Geology PC1"
    ) +
    scale_color_gradient2(
      low = "#1b9e77",
      mid = "white",
      high = "#7570b3",
      midpoint = 0,
      na.value = "#636363",
      name = "Geology PC2"
    ) +
    labs(
      title = panel_title,
      x = "Dynamic Storage",
      y = "Mobile Storage"
    ) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      plot.title = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE - 1),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
    )

  list(
    plot = p,
    stats = tibble(
      plot_name = plot_name,
      mobile_axis_definition = mobile_axis_definition,
      color_variable = "Geology PC1 + Geology PC2",
      slope = fit_out$slope,
      r2 = fit_out$r2,
      p_value = fit_out$p_value,
      n = nrow(plot_df)
    )
  )
}

make_geology_pc_size_plot <- function(
  data_in,
  mobile_col,
  panel_title,
  plot_name,
  mobile_axis_definition,
  require_chs = FALSE
) {
  plot_df <- data_in %>%
    filter(
      is.finite(dynamic_storage_strength_z),
      is.finite(.data[[mobile_col]]),
      is.finite(geology_pc1),
      is.finite(geology_pc2)
    )
  if (isTRUE(require_chs)) {
    plot_df <- plot_df %>% filter(is.finite(CHS_mean))
  }

  if (nrow(plot_df) == 0) {
    return(list(
      plot = ggplot() + theme_void(),
      stats = tibble(
        plot_name = plot_name,
        mobile_axis_definition = mobile_axis_definition,
        color_variable = "Geology PC1 (fill) + Geology PC2 (size)",
        slope = NA_real_,
        r2 = NA_real_,
        p_value = NA_real_,
        n = 0L
      )
    ))
  }

  x_rng <- range(plot_df$dynamic_storage_strength_z, na.rm = TRUE)
  y_rng <- range(plot_df[[mobile_col]], na.rm = TRUE)
  x_span <- ifelse(is.finite(diff(x_rng)) && diff(x_rng) > 0, diff(x_rng), 1)
  y_span <- ifelse(is.finite(diff(y_rng)) && diff(y_rng) > 0, diff(y_rng), 1)
  ann_x <- x_rng[2] - 0.03 * x_span
  ann_y <- y_rng[2] - 0.22 * y_span
  fit_out <- fit_stats(plot_df, "dynamic_storage_strength_z", mobile_col)
  quad_df <- build_quadrant_annotations(
    x_rng = x_rng,
    y_rng = y_rng,
    include_chs = grepl("^includes_CHS", mobile_axis_definition)
  )
  quad_label_size <- max(FIG_ANNOT_TEXT_SIZE - 1, 2.8)

  p <- ggplot(
    plot_df,
    aes(x = dynamic_storage_strength_z, y = .data[[mobile_col]])
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey75") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey75") +
    geom_point(
      aes(fill = geology_pc1, size = geology_pc2),
      shape = 21,
      color = "black",
      stroke = 0.55,
      alpha = 0.95
    ) +
    geom_text(
      aes(label = site),
      color = "black",
      nudge_y = 0.05,
      size = FIG_ANNOT_TEXT_SIZE,
      show.legend = FALSE
    ) +
    geom_smooth(
      data = plot_df,
      aes(x = dynamic_storage_strength_z, y = .data[[mobile_col]]),
      inherit.aes = FALSE,
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    geom_label(
      data = quad_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = quad_label_size,
      color = "grey20",
      lineheight = 0.9,
      fill = grDevices::adjustcolor("white", alpha.f = 0.72),
      linewidth = 0.15,
      label.padding = grid::unit(0.10, "lines"),
      label.r = grid::unit(0.08, "lines")
    ) +
    annotate(
      "label",
      x = ann_x,
      y = ann_y,
      label = fit_out$label,
      hjust = 1,
      vjust = 1,
      size = FIG_ANNOT_TEXT_SIZE,
      linewidth = 0.2,
      fill = "white",
      alpha = 0.9
    ) +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "#F7F7F7",
      high = "#B40426",
      midpoint = 0,
      na.value = "#bdbdbd",
      name = "Geology PC1"
    ) +
    scale_size_continuous(
      range = c(2.2, 8.0),
      name = "Geology PC2"
    ) +
    labs(
      title = panel_title,
      x = "Dynamic Storage",
      y = "Mobile Storage"
    ) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      plot.title = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE - 1),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
    )

  list(
    plot = p,
    stats = tibble(
      plot_name = plot_name,
      mobile_axis_definition = mobile_axis_definition,
      color_variable = "Geology PC1 (fill) + Geology PC2 (size)",
      slope = fit_out$slope,
      r2 = fit_out$r2,
      p_value = fit_out$p_value,
      n = nrow(plot_df)
    )
  )
}

res_chs_with_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "CHS_mean",
  color_label = "CHS",
  panel_title = "Dynamic vs Mobile Storage (all sites)",
  palette_type = "chs",
  plot_name = "dynamic_mobile_storage_chs_color_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_chs_no_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "CHS_mean",
  color_label = "CHS",
  panel_title = "Dynamic vs Mobile Storage (CHS-complete cases)",
  palette_type = "chs",
  plot_name = "dynamic_mobile_storage_chs_color_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

res_landslide_with_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Landslide_Total",
  color_label = "Total Landslide",
  panel_title = "Dynamic vs Mobile Storage (all sites)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_total_landslide_color_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_landslide_no_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Landslide_Total",
  color_label = "Total Landslide",
  panel_title = "Dynamic vs Mobile Storage (CHS-complete cases)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_total_landslide_color_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

res_lava1_with_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Lava1_per",
  color_label = "Lava 1 (%)",
  panel_title = "Dynamic vs Mobile Storage (Lava 1; all sites)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_lava1_color_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_lava1_no_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Lava1_per",
  color_label = "Lava 1 (%)",
  panel_title = "Dynamic vs Mobile Storage (Lava 1; CHS-complete cases)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_lava1_color_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

res_lava2_with_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Lava2_per",
  color_label = "Lava 2 (%)",
  panel_title = "Dynamic vs Mobile Storage (Lava 2; all sites)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_lava2_color_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_lava2_no_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Lava2_per",
  color_label = "Lava 2 (%)",
  panel_title = "Dynamic vs Mobile Storage (Lava 2; CHS-complete cases)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_lava2_color_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

res_ash_with_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Ash_Per",
  color_label = "Ash (%)",
  panel_title = "Dynamic vs Mobile Storage (Ash; all sites)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_ash_color_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_ash_no_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Ash_Per",
  color_label = "Ash (%)",
  panel_title = "Dynamic vs Mobile Storage (Ash; CHS-complete cases)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_ash_color_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

res_pyro_with_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Pyro_per",
  color_label = "Pyro (%)",
  panel_title = "Dynamic vs Mobile Storage (Pyro; all sites)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_pyro_color_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_pyro_no_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "Pyro_per",
  color_label = "Pyro (%)",
  panel_title = "Dynamic vs Mobile Storage (Pyro; CHS-complete cases)",
  palette_type = "landslide",
  plot_name = "dynamic_mobile_storage_pyro_color_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

res_geology_pc1_with_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "geology_pc1",
  color_label = "Geology PC1",
  panel_title = "Dynamic vs Mobile Storage (Geology PC1; all sites)",
  palette_type = "geology_pc",
  plot_name = "dynamic_mobile_storage_geology_pc1_color_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_geology_pc1_no_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "geology_pc1",
  color_label = "Geology PC1",
  panel_title = "Dynamic vs Mobile Storage (Geology PC1; CHS-complete cases)",
  palette_type = "geology_pc",
  plot_name = "dynamic_mobile_storage_geology_pc1_color_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

res_geology_pc2_with_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "geology_pc2",
  color_label = "Geology PC2",
  panel_title = "Dynamic vs Mobile Storage (Geology PC2; all sites)",
  palette_type = "geology_pc",
  plot_name = "dynamic_mobile_storage_geology_pc2_color_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_geology_pc2_no_chs <- make_state_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  color_col = "geology_pc2",
  color_label = "Geology PC2",
  panel_title = "Dynamic vs Mobile Storage (Geology PC2; CHS-complete cases)",
  palette_type = "geology_pc",
  plot_name = "dynamic_mobile_storage_geology_pc2_color_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

res_geology_pc_bivar_with_chs <- make_geology_pc_bivariate_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  panel_title = "Dynamic vs Mobile Storage (Geology PC1 fill + PC2 outline; all sites)",
  plot_name = "dynamic_mobile_storage_geology_pc_bivariate_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_geology_pc_bivar_no_chs <- make_geology_pc_bivariate_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  panel_title = "Dynamic vs Mobile Storage (Geology PC1 fill + PC2 outline; CHS-complete cases)",
  plot_name = "dynamic_mobile_storage_geology_pc_bivariate_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

res_geology_pc_size_with_chs <- make_geology_pc_size_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  panel_title = "Dynamic vs Mobile Storage (Geology PC1 fill + PC2 size; all sites)",
  plot_name = "dynamic_mobile_storage_geology_pcsize_all_sites",
  mobile_axis_definition = "includes_CHS_all_sites",
  require_chs = FALSE
)

res_geology_pc_size_no_chs <- make_geology_pc_size_plot(
  scatter_df_all,
  mobile_col = "mobile_mixing_with_chs_z",
  panel_title = "Dynamic vs Mobile Storage (Geology PC1 fill + PC2 size; CHS-complete cases)",
  plot_name = "dynamic_mobile_storage_geology_pcsize_chs_complete_cases",
  mobile_axis_definition = "includes_CHS_complete_cases",
  require_chs = TRUE
)

# Remove older state-space color variants from current workflow outputs.
unlink(file.path(main_fig_dir, c(
  "unified_framework_state_space_geology_all_sites.png",
  "unified_framework_state_space_geomorphology_all_sites.png",
  "unified_framework_state_space_geology_no_ws09_look.png",
  "unified_framework_state_space_geomorphology_no_ws09_look.png",
  "unified_framework_state_space.png"
)))
unlink(file.path(main_pdf_dir, c(
  "unified_framework_state_space_geology_all_sites.pdf",
  "unified_framework_state_space_geomorphology_all_sites.pdf",
  "unified_framework_state_space_geology_no_ws09_look.pdf",
  "unified_framework_state_space_geomorphology_no_ws09_look.pdf",
  "unified_framework_state_space.pdf"
)))

plot_specs <- list(
  list(name = "dynamic_mobile_storage_chs_color_all_sites", plot = res_chs_with_chs$plot),
  list(name = "dynamic_mobile_storage_chs_color_chs_complete_cases", plot = res_chs_no_chs$plot),
  list(name = "dynamic_mobile_storage_total_landslide_color_all_sites", plot = res_landslide_with_chs$plot),
  list(name = "dynamic_mobile_storage_total_landslide_color_chs_complete_cases", plot = res_landslide_no_chs$plot),
  list(name = "dynamic_mobile_storage_lava1_color_all_sites", plot = res_lava1_with_chs$plot),
  list(name = "dynamic_mobile_storage_lava1_color_chs_complete_cases", plot = res_lava1_no_chs$plot),
  list(name = "dynamic_mobile_storage_lava2_color_all_sites", plot = res_lava2_with_chs$plot),
  list(name = "dynamic_mobile_storage_lava2_color_chs_complete_cases", plot = res_lava2_no_chs$plot),
  list(name = "dynamic_mobile_storage_ash_color_all_sites", plot = res_ash_with_chs$plot),
  list(name = "dynamic_mobile_storage_ash_color_chs_complete_cases", plot = res_ash_no_chs$plot),
  list(name = "dynamic_mobile_storage_pyro_color_all_sites", plot = res_pyro_with_chs$plot),
  list(name = "dynamic_mobile_storage_pyro_color_chs_complete_cases", plot = res_pyro_no_chs$plot),
  list(name = "dynamic_mobile_storage_geology_pc1_color_all_sites", plot = res_geology_pc1_with_chs$plot),
  list(name = "dynamic_mobile_storage_geology_pc1_color_chs_complete_cases", plot = res_geology_pc1_no_chs$plot),
  list(name = "dynamic_mobile_storage_geology_pc2_color_all_sites", plot = res_geology_pc2_with_chs$plot),
  list(name = "dynamic_mobile_storage_geology_pc2_color_chs_complete_cases", plot = res_geology_pc2_no_chs$plot),
  list(name = "dynamic_mobile_storage_geology_pc_bivariate_all_sites", plot = res_geology_pc_bivar_with_chs$plot),
  list(name = "dynamic_mobile_storage_geology_pc_bivariate_chs_complete_cases", plot = res_geology_pc_bivar_no_chs$plot),
  list(name = "dynamic_mobile_storage_geology_pcsize_all_sites", plot = res_geology_pc_size_with_chs$plot),
  list(name = "dynamic_mobile_storage_geology_pcsize_chs_complete_cases", plot = res_geology_pc_size_no_chs$plot)
)

state_fit_summary <- bind_rows(
  res_chs_with_chs$stats,
  res_chs_no_chs$stats,
  res_landslide_with_chs$stats,
  res_landslide_no_chs$stats,
  res_lava1_with_chs$stats,
  res_lava1_no_chs$stats,
  res_lava2_with_chs$stats,
  res_lava2_no_chs$stats,
  res_ash_with_chs$stats,
  res_ash_no_chs$stats,
  res_pyro_with_chs$stats,
  res_pyro_no_chs$stats,
  res_geology_pc1_with_chs$stats,
  res_geology_pc1_no_chs$stats,
  res_geology_pc2_with_chs$stats,
  res_geology_pc2_no_chs$stats,
  res_geology_pc_bivar_with_chs$stats,
  res_geology_pc_bivar_no_chs$stats,
  res_geology_pc_size_with_chs$stats,
  res_geology_pc_size_no_chs$stats
)
write_csv(
  state_fit_summary,
  file.path(model_out_dir, "dynamic_mobile_storage_regression_summary.csv")
)

# Remove previous dynamic-v-mobile files so naming updates do not leave stale plots.
unlink(file.path(dynamic_v_mobile_fig_dir, "dynamic_mobile_storage_*.png"))
unlink(file.path(dynamic_v_mobile_pdf_dir, "dynamic_mobile_storage_*.pdf"))

for (spec in plot_specs) {
  # Remove legacy copies from main root to keep these plots under dynamic_v_mobile.
  unlink(file.path(main_fig_dir, paste0(spec$name, ".png")))
  unlink(file.path(main_pdf_dir, paste0(spec$name, ".pdf")))

  ggsave(
    file.path(dynamic_v_mobile_fig_dir, paste0(spec$name, ".png")),
    spec$plot,
    width = 8 * FIG_WIDTH_SCALE,
    height = 6 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
  ggsave(
    file.path(dynamic_v_mobile_pdf_dir, paste0(spec$name, ".pdf")),
    spec$plot,
    width = 8 * FIG_WIDTH_SCALE,
    height = 6 * FIG_HEIGHT_SCALE
  )
}

geology_loadings_file <- file.path(model_out_dir, "geology_composition_pca_loadings.csv")
geology_variance_file <- file.path(model_out_dir, "geology_composition_pca_variance.csv")

geology_scores <- axes %>%
  transmute(
    site = as.character(site),
    geology_pc1 = as.numeric(geology_pc1),
    geology_pc2 = as.numeric(geology_pc2)
  ) %>%
  filter(is.finite(geology_pc1), is.finite(geology_pc2))

if (nrow(geology_scores) >= 2 &&
  file.exists(geology_loadings_file) &&
  file.exists(geology_variance_file)) {
  geology_loadings <- read_csv(geology_loadings_file, show_col_types = FALSE) %>%
    transmute(
      variable = dplyr::recode(
        variable,
        "Lava1_per" = "Lava 1",
        "Lava2_per" = "Lava 2",
        "Ash_Per" = "Ash",
        "Pyro_per" = "Pyro",
        .default = variable
      ),
      lx = as.numeric(geology_pc1_loading),
      ly = as.numeric(geology_pc2_loading)
    ) %>%
    filter(is.finite(lx), is.finite(ly))

  geology_var <- read_csv(geology_variance_file, show_col_types = FALSE)
  pc1_var <- geology_var %>%
    filter(pc == "PC1") %>%
    pull(variance_explained)
  pc2_var <- geology_var %>%
    filter(pc == "PC2") %>%
    pull(variance_explained)
  pc1_pct <- if (length(pc1_var) > 0 && is.finite(pc1_var[1])) {
    100 * pc1_var[1]
  } else {
    NA_real_
  }
  pc2_pct <- if (length(pc2_var) > 0 && is.finite(pc2_var[1])) {
    100 * pc2_var[1]
  } else {
    NA_real_
  }

  if (nrow(geology_loadings) > 0) {
    score_span_x <- diff(range(geology_scores$geology_pc1, na.rm = TRUE))
    score_span_y <- diff(range(geology_scores$geology_pc2, na.rm = TRUE))
    load_span_x <- diff(range(geology_loadings$lx, na.rm = TRUE))
    load_span_y <- diff(range(geology_loadings$ly, na.rm = TRUE))
    if (!is.finite(score_span_x) || score_span_x <= 0) score_span_x <- 1
    if (!is.finite(score_span_y) || score_span_y <= 0) score_span_y <- 1
    if (!is.finite(load_span_x) || load_span_x <= 0) load_span_x <- 1
    if (!is.finite(load_span_y) || load_span_y <= 0) load_span_y <- 1
    arrow_scale <- 0.35 * min(score_span_x / load_span_x, score_span_y / load_span_y)
    if (!is.finite(arrow_scale) || arrow_scale <= 0) arrow_scale <- 1
    geology_loadings <- geology_loadings %>%
      mutate(
        xend = lx * arrow_scale,
        yend = ly * arrow_scale,
        xlab = xend * 1.10,
        ylab = yend * 1.10
      )
  }

  x_lab <- if (is.finite(pc1_pct)) {
    sprintf("Geology PC1 (%.1f%%)", pc1_pct)
  } else {
    "Geology PC1"
  }
  y_lab <- if (is.finite(pc2_pct)) {
    sprintf("Geology PC2 (%.1f%%)", pc2_pct)
  } else {
    "Geology PC2"
  }

  p_geology_pca <- ggplot(
    geology_scores,
    aes(x = geology_pc1, y = geology_pc2)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
    geom_point(aes(color = site), size = FIG_POINT_SIZE_LARGE + 0.8) +
    geom_text(
      aes(label = site),
      color = "black",
      nudge_y = 0.06,
      size = FIG_ANNOT_TEXT_SIZE,
      show.legend = FALSE
    ) +
    scale_color_manual(values = SITE_COLORS, guide = "none")

  if (nrow(geology_loadings) > 0) {
    p_geology_pca <- p_geology_pca +
      geom_segment(
        data = geology_loadings,
        aes(x = 0, y = 0, xend = xend, yend = yend),
        inherit.aes = FALSE,
        arrow = arrow(length = grid::unit(0.015, "npc")),
        linewidth = 0.6,
        color = "black"
      ) +
      geom_text(
        data = geology_loadings,
        aes(x = xlab, y = ylab, label = variable),
        inherit.aes = FALSE,
        size = FIG_ANNOT_TEXT_SIZE,
        color = "black"
      )
  }

  p_geology_pca <- p_geology_pca +
    labs(
      title = "Geology Composition PCA",
      x = x_lab,
      y = y_lab
    ) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      plot.title = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0)
    )

  unlink(file.path(main_fig_dir, "geology_composition_pca_biplot.png"))
  unlink(file.path(main_pdf_dir, "geology_composition_pca_biplot.pdf"))

  ggsave(
    file.path(dynamic_v_mobile_fig_dir, "geology_composition_pca_biplot.png"),
    p_geology_pca,
    width = 8 * FIG_WIDTH_SCALE,
    height = 6 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
  ggsave(
    file.path(dynamic_v_mobile_pdf_dir, "geology_composition_pca_biplot.pdf"),
    p_geology_pca,
    width = 8 * FIG_WIDTH_SCALE,
    height = 6 * FIG_HEIGHT_SCALE
  )
}

annual_df <- read_csv(annual_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(year = as.integer(year))

if (!("P_WetSeason" %in% names(annual_df))) {
  annual_df$P_WetSeason <- NA_real_
}
if ("precip_nov_may_mm" %in% names(annual_df)) {
  annual_df$P_WetSeason <- dplyr::coalesce(annual_df$P_WetSeason, annual_df$precip_nov_may_mm)
}
if ("P_NovJan" %in% names(annual_df)) {
  annual_df$P_WetSeason <- dplyr::coalesce(annual_df$P_WetSeason, annual_df$P_NovJan)
}
if ("precip_nov_jan_mm" %in% names(annual_df)) {
  annual_df$P_WetSeason <- dplyr::coalesce(annual_df$P_WetSeason, annual_df$precip_nov_jan_mm)
}

pca_scores <- read_csv(pca_scores_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  transmute(
    site = as.character(site),
    year = as.integer(year),
    dynamic_pc1 = as.numeric(PC1),
    dynamic_pc2 = as.numeric(PC2)
  )

traj_df <- pca_scores %>%
  left_join(
    annual_df %>%
      transmute(site = as.character(site), year, P_WetSeason),
    by = c("site", "year")
  ) %>%
  group_by(site) %>%
  mutate(P_WetSeason_anom = P_WetSeason - mean(P_WetSeason, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(
    is.finite(dynamic_pc1),
    is.finite(dynamic_pc2),
    is.finite(P_WetSeason_anom)
  ) %>%
  mutate(
    site = factor(site, levels = SITE_ORDER_HYDROMETRIC),
    site_color = dplyr::coalesce(unname(SITE_COLORS[as.character(site)]), "#999999")
  )

site_label_offsets <- tibble::tribble(
  ~site,  ~dx,    ~dy,
  "WS10", -0.22,  0.12,
  "WS01",  0.20,  0.10,
  "WS09", -0.04, -0.12
)

traj_centroids <- traj_df %>%
  group_by(site) %>%
  summarise(
    dynamic_pc1 = mean(dynamic_pc1, na.rm = TRUE),
    dynamic_pc2 = mean(dynamic_pc2, na.rm = TRUE),
    site_color = dplyr::first(site_color),
    .groups = "drop"
  ) %>%
  left_join(site_label_offsets, by = "site") %>%
  mutate(
    dx = dplyr::coalesce(dx, 0),
    dy = dplyr::coalesce(dy, 0.08),
    x_lab = dynamic_pc1 + dx,
    y_lab = dynamic_pc2 + dy
  )

p_pca_wetseason <- ggplot(traj_df, aes(x = dynamic_pc1, y = dynamic_pc2)) +
  stat_ellipse(
    aes(group = site, fill = I(site_color), color = I(site_color)),
    geom = "polygon",
    type = "norm",
    level = 0.80,
    linewidth = 0.5,
    alpha = 0.10
  ) +
  geom_point(
    aes(color = P_WetSeason_anom),
    size = FIG_POINT_SIZE_SMALL + 0.6,
    alpha = 0.92
  ) +
  geom_text(
    data = traj_centroids,
    aes(x = x_lab, y = y_lab, label = site, color = I(site_color)),
    size = FIG_ANNOT_TEXT_SIZE,
    show.legend = FALSE,
    check_overlap = FALSE
  ) +
  scale_color_gradient2(
    low = "#b2182b",
    mid = "white",
    high = "#2166ac",
    midpoint = 0,
    name = "Annual\nPwetseason anomaly"
  ) +
  labs(
    title = "Precip anomaly PCA",
    x = "Annual Dynamic Storage PC1",
    y = "Annual Dynamic Storage PC2"
  ) +
  theme_pub() +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    plot.title = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0),
    legend.title = element_text(size = FIG_AXIS_TITLE_SIZE - 1),
    legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
  )

unlink(file.path(main_fig_dir, c(
  "PCA_precip_anomaly.png",
  "unified_framework_dynamic_trajectory_pwetseason_anom.png",
  "unified_framework_dynamic_trajectory_wb_anom.png",
  "unified_framework_dynamic_trajectory_pwinter_anom.png"
)))
unlink(file.path(main_pdf_dir, c(
  "PCA_precip_anomaly.pdf",
  "unified_framework_dynamic_trajectory_pwetseason_anom.pdf",
  "unified_framework_dynamic_trajectory_wb_anom.pdf",
  "unified_framework_dynamic_trajectory_pwinter_anom.pdf"
)))
unlink(file.path(supp_fig_dir, c(
  "unified_framework_dynamic_trajectory_pwetseason_anom.png",
  "PCA_precip_anomaly.png",
  "unified_framework_dynamic_trajectory_wb_anom.png",
  "unified_framework_dynamic_trajectory_pwinter_anom.png"
)))
unlink(file.path(supp_pdf_dir, c(
  "unified_framework_dynamic_trajectory_pwetseason_anom.pdf",
  "PCA_precip_anomaly.pdf",
  "unified_framework_dynamic_trajectory_wb_anom.pdf",
  "unified_framework_dynamic_trajectory_pwinter_anom.pdf"
)))

ggsave(
  file.path(supp_fig_dir, "PCA_precip_anomaly.png"),
  p_pca_wetseason,
  width = 8.5 * FIG_WIDTH_SCALE,
  height = 6.2 * FIG_HEIGHT_SCALE,
  dpi = 300
)
ggsave(
  file.path(supp_pdf_dir, "PCA_precip_anomaly.pdf"),
  p_pca_wetseason,
  width = 8.5 * FIG_WIDTH_SCALE,
  height = 6.2 * FIG_HEIGHT_SCALE
)
