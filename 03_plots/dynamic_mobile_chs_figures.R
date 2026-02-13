# -----------------------------------------------------------------------------
# Main Text and Supplement Figures - HJA Dynamic Storage
# -----------------------------------------------------------------------------
#
# Main Text Figures:
#   Figure 3: Dynamic Storage (RBI, RCS, FDC, SD) - faceted
#   Figure 4: Mobile Isotope (MTT, Fyw, DR) - faceted
#   Figure 5: CHS (baseflow fraction) - single panel
#
# Supplement Figures:
#   - Faceted time series for Dynamic Storage metrics
#   - Faceted time series for Mobile Storage metrics
#
# Author: Sidney Bush
# Date: 2026-02-02
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(colorspace)

rm(list = ls())

# -----------------------------------------------------------------------------
# SOURCE CONFIGURATION
# -----------------------------------------------------------------------------

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
source(config_path)

# -----------------------------------------------------------------------------
# SETUP
# -----------------------------------------------------------------------------

base_dir <- BASE_DATA_DIR
main_dir <- file.path(FIGURES_DIR, "main")
supp_dir <- file.path(FIGURES_DIR, "supp")

for (d in c(main_dir, supp_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

site_order <- SITE_ORDER_HYDROMETRIC
site_order_iso <- SITE_ORDER_ALL
site_colors <- SITE_COLORS
site_labels_panel <- make_panel_label_map(site_order)

# Metric labels
metric_labels <- c(
  "RBI" = "RBI", "RCS" = "RCS", "FDC" = "FDC", "SD" = "SD (mm)",
  "WB" = "WB (mm)", "CHS" = "CHS", "MTT" = "MTT (yr)",
  "Fyw" = "Fyw", "DR" = "DR"
)

# Publication theme
theme_pub_base <- theme_pub() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT_SIZE),
    axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    strip.background = element_blank(),
    strip.text = element_text(size = FIG_STRIP_TEXT_SIZE, hjust = 0),
    legend.position = "none"
  )

theme_set(theme_pub_base)

# -----------------------------------------------------------------------------
# HELPER: Standardize site names to WS## format
# -----------------------------------------------------------------------------

standardize_sites <- function(df, allowed_sites = site_order) {
  df %>%
    mutate(
      site = standardize_site_code(site)
    ) %>%
    filter(site %in% allowed_sites) %>%
    mutate(site = factor(site, levels = allowed_sites))
}

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------


# Annual storage metrics
annual_file <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUT_MASTER_DIR, LEGACY_ANNUAL_FILE)
}

annual_data <- read_csv(annual_file, show_col_types = FALSE) %>%
  rename_legacy_storage_metrics() %>%
  standardize_sites(allowed_sites = site_order)


# Isotope data for mobile storage (site-level metrics)
isotope_file <- file.path(OUT_MET_MOBILE_DIR, "isotope_metrics_site.csv")
if (!file.exists(isotope_file)) {
  isotope_file <- file.path(OUT_MET_MOBILE_DIR, "Isotope_Metrics_Site.csv")
}
if (!file.exists(isotope_file)) {
  isotope_file <- file.path(ISOTOPE_DIR, "Isotope_Metrics_Site.csv")
}

if (file.exists(isotope_file)) {
  isotope_data <- read_csv(isotope_file, show_col_types = FALSE) %>%
    standardize_sites(allowed_sites = site_order_iso)

  # Add MTT1/MTT2 columns from raw isotope table when available.
  mtt_file <- file.path(ISOTOPE_DIR, "MTT_FYW.csv")
  if (file.exists(mtt_file)) {
    mean_or_na <- function(x) {
      if (all(is.na(x))) {
        NA_real_
      } else {
        mean(x, na.rm = TRUE)
      }
    }

    mtt_raw <- read_csv(mtt_file, show_col_types = FALSE) %>%
      mutate(site = standardize_site_code(site)) %>%
      filter(site %in% site_order_iso)

    first_col <- function(df, candidates) {
      hit <- candidates[candidates %in% names(df)]
      if (length(hit) == 0) {
        return(rep(NA_real_, nrow(df)))
      }
      as.numeric(df[[hit[1]]])
    }

    mtt_raw$MTT1_val <- first_col(mtt_raw, c("MTT1"))
    mtt_raw$MTT1_sd <- first_col(mtt_raw, c("MTT1_SD"))
    mtt_raw$MTT1_low <- first_col(mtt_raw, c("MTT1L", "MTT1_low", "MTT1_MIN"))
    mtt_raw$MTT1_high <- first_col(mtt_raw, c("MTT1H", "MTT1_high", "MTT1_MAX"))
    mtt_raw$MTT1_err <- ifelse(
      !is.na(mtt_raw$MTT1_sd),
      abs(mtt_raw$MTT1_sd),
      abs((mtt_raw$MTT1_high - mtt_raw$MTT1_low) / 2)
    )

    mtt_raw$MTT2_val <- first_col(mtt_raw, c("MTT2M", "MTT2"))
    if (all(is.na(mtt_raw$MTT2_val))) {
      mtt2_low <- first_col(mtt_raw, c("MTT2L", "MTT2_low", "MTT2_MIN"))
      mtt2_high <- first_col(mtt_raw, c("MTT2H", "MTT2_high", "MTT2_MAX"))
      mtt_raw$MTT2_val <- rowMeans(cbind(mtt2_low, mtt2_high), na.rm = TRUE)
    }
    mtt_raw$MTT2_low <- first_col(mtt_raw, c("MTT2L", "MTT2_low", "MTT2_MIN"))
    mtt_raw$MTT2_high <- first_col(mtt_raw, c("MTT2H", "MTT2_high", "MTT2_MAX"))
    mtt_raw$MTT2_low_sd <- first_col(mtt_raw, c("MTT2L_SD", "MTT2_low_SD"))
    mtt_raw$MTT2_high_sd <- first_col(mtt_raw, c("MTT2H_SD", "MTT2_high_SD"))
    mtt_raw$MTT2_err <- ifelse(
      is.finite(mtt_raw$MTT2_low_sd) & is.finite(mtt_raw$MTT2_high_sd),
      abs((mtt_raw$MTT2_low_sd + mtt_raw$MTT2_high_sd) / 2),
      abs((mtt_raw$MTT2_high - mtt_raw$MTT2_low) / 2)
    )

    mtt_raw$Fyw_low <- first_col(mtt_raw, c("FywL", "FYWL", "Fyw_low", "FYWL"))
    mtt_raw$Fyw_high <- first_col(mtt_raw, c("FywH", "FYWH", "Fyw_high", "FYWH"))
    mtt_raw$Fyw_val <- first_col(mtt_raw, c("Fyw", "FYW", "fyw", "FYWM"))
    if (all(is.na(mtt_raw$Fyw_val))) {
      mtt_raw$Fyw_val <- rowMeans(cbind(mtt_raw$Fyw_low, mtt_raw$Fyw_high), na.rm = TRUE)
    }
    mtt_raw$Fyw_low_sd <- first_col(mtt_raw, c("FYWL_SD", "FywL_SD", "Fyw_low_SD"))
    mtt_raw$Fyw_high_sd <- first_col(mtt_raw, c("FYWH_SD", "FywH_SD", "Fyw_high_SD"))
    mtt_raw$Fyw_err <- ifelse(
      is.finite(mtt_raw$Fyw_low_sd) & is.finite(mtt_raw$Fyw_high_sd),
      abs((mtt_raw$Fyw_low_sd + mtt_raw$Fyw_high_sd) / 2),
      abs((mtt_raw$Fyw_high - mtt_raw$Fyw_low) / 2)
    )

    mtt_versions <- mtt_raw %>%
      mutate(
        Fyw_val = ifelse(is.na(Fyw_val), NA_real_, Fyw_val),
        Fyw_err = ifelse(is.na(Fyw_err), NA_real_, abs(Fyw_err))
      ) %>%
      group_by(site) %>%
      summarise(
        MTT1 = mean_or_na(MTT1_val),
        MTT1_err = mean_or_na(abs(MTT1_err)),
        MTT2 = mean_or_na(MTT2_val),
        MTT2_err = mean_or_na(abs(MTT2_err)),
        Fyw = mean_or_na(Fyw_val),
        Fyw_err = mean_or_na(Fyw_err),
        .groups = "drop"
      )

    isotope_data <- isotope_data %>%
      left_join(mtt_versions, by = "site")

    # Resolve duplicated columns after join so plotting uses the merged value.
    merge_col <- function(df, nm) {
      x <- paste0(nm, ".x")
      y <- paste0(nm, ".y")
      if (x %in% names(df) || y %in% names(df)) {
        xv <- if (x %in% names(df)) df[[x]] else NA_real_
        yv <- if (y %in% names(df)) df[[y]] else NA_real_
        df[[nm]] <- dplyr::coalesce(suppressWarnings(as.numeric(xv)), suppressWarnings(as.numeric(yv)))
        if (x %in% names(df)) df[[x]] <- NULL
        if (y %in% names(df)) df[[y]] <- NULL
      }
      df
    }
    for (nm in c("MTT1", "MTT1_err", "MTT2", "MTT2_err", "Fyw", "Fyw_err", "DR", "DR_err")) {
      isotope_data <- merge_col(isotope_data, nm)
    }
  }

  if (!("MTT1" %in% names(isotope_data))) isotope_data$MTT1 <- NA_real_
  if (!("MTT1_err" %in% names(isotope_data))) isotope_data$MTT1_err <- NA_real_
  if (!("MTT2" %in% names(isotope_data))) {
    isotope_data$MTT2 <- if ("MTT" %in% names(isotope_data)) isotope_data$MTT else NA_real_
  }
  if (!("MTT2_err" %in% names(isotope_data))) isotope_data$MTT2_err <- NA_real_
  if (!("Fyw" %in% names(isotope_data))) isotope_data$Fyw <- NA_real_
  if (!("Fyw_err" %in% names(isotope_data))) {
    fyw_low_name <- names(isotope_data)[tolower(names(isotope_data)) %in% c("fywl", "fyw_low")]
    fyw_high_name <- names(isotope_data)[tolower(names(isotope_data)) %in% c("fywh", "fyw_high")]
    if (length(fyw_low_name) > 0 && length(fyw_high_name) > 0) {
      isotope_data$Fyw_err <- abs((isotope_data[[fyw_high_name[1]]] - isotope_data[[fyw_low_name[1]]]) / 2)
    } else {
      isotope_data$Fyw_err <- NA_real_
    }
  }
  isotope_data <- isotope_data %>%
    complete(site = factor(site_order_iso, levels = site_order_iso))

  # Save the exact values used for isotope plotting so figure inputs are traceable
  isotope_plot_values <- isotope_data %>%
    transmute(
      site = as.character(site),
      MTT1,
      MTT1_err,
      MTT2,
      MTT2_err,
      Fyw,
      Fyw_err,
      DR,
      DR_err
    )

  write_csv(
    isotope_plot_values,
    file.path(OUT_MET_SUPPORT_DIR, "ms_isotope_plot_values.csv")
  )
  cat("  Loaded isotope data:", nrow(isotope_data), "rows\n")
} else {
  isotope_data <- NULL
  cat("  Warning: Isotope data not found\n")
}

# CHS data
chs_file <- file.path(EC_DIR, "CHS_annual_baseflow.csv")
if (!file.exists(chs_file)) {
  chs_file <- file.path(OUT_MET_MOBILE_DIR, "annual_gw_prop.csv")
}
if (!file.exists(chs_file)) {
  chs_file <- file.path(OUT_MET_MOBILE_DIR, "Annual_GW_Prop.csv")
}

if (file.exists(chs_file)) {
  chs_data <- read_csv(chs_file, show_col_types = FALSE) %>%
    mutate(
      site = if ("site" %in% names(.)) site else if ("SITECODE" %in% names(.)) SITECODE else NA_character_,
      year = if ("year" %in% names(.)) year else if ("waterYear" %in% names(.)) waterYear else NA_real_,
      CHS = if ("CHS" %in% names(.)) CHS else if ("mean_bf" %in% names(.)) mean_bf else CHS
    ) %>%
    select(site, year, CHS) %>%
    standardize_sites()
  cat("  Loaded CHS data:", nrow(chs_data), "rows\n")
} else {
  # Try to get CHS from annual data
  if ("CHS" %in% names(annual_data)) {
    chs_data <- annual_data %>%
      select(site, year, CHS) %>%
      filter(!is.na(CHS))
    cat("  Using CHS from annual data:", nrow(chs_data), "rows\n")
  } else {
    chs_data <- NULL
    cat("  Warning: CHS data not found\n")
  }
}

# -----------------------------------------------------------------------------
# FIGURE 3: DYNAMIC STORAGE (RBI, RCS, FDC, SD)
# -----------------------------------------------------------------------------


dynamic_metrics <- c("RBI", "RCS", "FDC", "SD")
metric_labels_panel <- make_panel_label_map(metric_labels[dynamic_metrics])

dynamic_long <- annual_data %>%
  select(site, year, any_of(dynamic_metrics)) %>%
  pivot_longer(cols = -c(site, year), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric, levels = dynamic_metrics))

# Calculate mean ± SD per site
dynamic_summary <- dynamic_long %>%
  group_by(site, metric) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    site = factor(site_order, levels = site_order),
    metric = factor(dynamic_metrics, levels = dynamic_metrics),
    fill = list(mean_val = NA_real_, sd_val = NA_real_)
  )

# Add mean labels for each facet
dynamic_labels <- dynamic_summary %>%
  group_by(metric) %>%
  summarise(
    x_pos = -Inf,
    y_pos = -Inf,
    .groups = "drop"
  )

fig3 <- ggplot(dynamic_summary, aes(x = site, y = mean_val, color = site)) +
  geom_point(size = FIG_POINT_SIZE_MED, na.rm = TRUE) +
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.3, linewidth = 0.5, na.rm = TRUE
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2,
             axes = "margins", axis.labels = "margins",
             labeller = labeller(metric = metric_labels_panel)) +
  scale_color_manual(values = site_colors) +
  scale_x_discrete(limits = site_order, drop = FALSE) +
  labs(x = NULL, y = "Value")

ggsave(file.path(main_dir, "ds_summary.png"), fig3, width = 10 * FIG_WIDTH_SCALE, height = 8 * FIG_HEIGHT_SCALE, dpi = 300)
ggsave(file.path(main_dir, "ds_summary.pdf"), fig3, width = 10 * FIG_WIDTH_SCALE, height = 8 * FIG_HEIGHT_SCALE)

# -----------------------------------------------------------------------------
# FIGURE 4: MOBILE ISOTOPE (MTT, Fyw, DR)
# -----------------------------------------------------------------------------


if (!is.null(isotope_data) && nrow(isotope_data) > 0) {
  make_mobile_panel <- function(df, y_col, y_lab, err_col = NULL) {
    p <- ggplot(df, aes(x = site, y = .data[[y_col]], color = site)) +
      geom_point(size = FIG_POINT_SIZE_LARGE, na.rm = TRUE) +
      scale_color_manual(values = site_colors) +
      scale_x_discrete(limits = site_order_iso, drop = FALSE) +
      labs(x = NULL, y = y_lab)

    if (!is.null(err_col) && err_col %in% names(df)) {
      p <- p + geom_errorbar(
        aes(ymin = .data[[y_col]] - .data[[err_col]], ymax = .data[[y_col]] + .data[[err_col]]),
        width = 0.3, linewidth = 0.5, na.rm = TRUE
      )
    }
    p
  }

  p_mtt1 <- make_mobile_panel(isotope_data, "MTT1", "MTT1 (yr)", err_col = "MTT1_err") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p_mtt2 <- make_mobile_panel(isotope_data, "MTT2", "MTT2 (yr)", err_col = "MTT2_err") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p_fyw <- make_mobile_panel(isotope_data, "Fyw", "Fyw", err_col = "Fyw_err")
  p_dr <- make_mobile_panel(isotope_data, "DR", "DR", err_col = "DR_err")

  fig4 <- (p_mtt1 | p_mtt2) / (p_fyw | p_dr)
  ggsave(file.path(main_dir, "ms_isotope.png"), fig4, width = 12 * FIG_WIDTH_SCALE, height = 8 * FIG_HEIGHT_SCALE, dpi = 300)
  ggsave(file.path(main_dir, "ms_isotope.pdf"), fig4, width = 12 * FIG_WIDTH_SCALE, height = 8 * FIG_HEIGHT_SCALE)
  cat("  Saved Figure 4\n")
} else {
  cat("  Skipping Figure 4: No isotope data\n")
}

# -----------------------------------------------------------------------------
# FIGURE 5: CHS (BASEFLOW FRACTION)
# -----------------------------------------------------------------------------


if (!is.null(chs_data) && nrow(chs_data) > 0) {
  # Summary stats
  chs_summary <- chs_data %>%
    group_by(site) %>%
    summarise(
      mean_val = mean(CHS, na.rm = TRUE),
      sd_val = sd(CHS, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    complete(
      site = factor(site_order, levels = site_order),
      fill = list(mean_val = NA_real_, sd_val = NA_real_, n = 0)
    )

  fig5 <- ggplot(chs_summary, aes(x = site, y = mean_val, color = site)) +
    geom_point(size = FIG_POINT_SIZE_LARGE, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
      width = 0.3, linewidth = 0.6, na.rm = TRUE
    ) +
    scale_color_manual(values = site_colors) +
    scale_x_discrete(limits = site_order, drop = FALSE) +
    labs(x = NULL, y = "Baseflow Fraction (CHS)")

  ggsave(file.path(main_dir, "ms_chs.png"), fig5, width = 8 * FIG_WIDTH_SCALE, height = 5 * FIG_HEIGHT_SCALE, dpi = 300)
  ggsave(file.path(main_dir, "ms_chs.pdf"), fig5, width = 8 * FIG_WIDTH_SCALE, height = 5 * FIG_HEIGHT_SCALE)

  fig5_box <- ggplot(chs_data, aes(x = site, y = CHS, color = site, fill = site)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2, na.rm = TRUE) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = FIG_POINT_SIZE_SMALL, alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = site_colors) +
    scale_fill_manual(values = site_colors) +
    scale_x_discrete(limits = site_order, drop = FALSE) +
    labs(x = NULL, y = "Baseflow Fraction (CHS)")

  ggsave(file.path(main_dir, "ms_chs_boxplot.png"), fig5_box, width = 8 * FIG_WIDTH_SCALE, height = 5 * FIG_HEIGHT_SCALE, dpi = 300)
  ggsave(file.path(main_dir, "ms_chs_boxplot.pdf"), fig5_box, width = 8 * FIG_WIDTH_SCALE, height = 5 * FIG_HEIGHT_SCALE)
  cat("  Saved Figure 5\n")
} else {
  cat("  Skipping Figure 5: No CHS data\n")
}

# -----------------------------------------------------------------------------
# SUPPLEMENT: FACETED TIME SERIES - DYNAMIC STORAGE
# -----------------------------------------------------------------------------


# Calculate site means for dashed reference lines
site_means_dynamic <- dynamic_long %>%
  group_by(site, metric) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    x_pos = max(year, na.rm = TRUE),
    .groups = "drop"
  )

dynamic_line_data <- dynamic_long %>%
  group_by(site, metric) %>%
  filter(sum(!is.na(value)) >= 2) %>%
  ungroup()

supp_dynamic_ts <- ggplot(dynamic_long, aes(x = year, y = value, color = site)) +
  geom_line(data = dynamic_line_data, linewidth = 0.5, na.rm = TRUE) +
  geom_point(size = FIG_POINT_SIZE_SMALL, na.rm = TRUE) +
  geom_hline(
    data = site_means_dynamic,
    aes(yintercept = mean_val),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.35,
    linetype = "dashed",
    na.rm = TRUE
  ) +
  geom_text(
    data = site_means_dynamic,
    aes(x = x_pos, y = mean_val, label = sprintf(paste0("%.", FIG_MEAN_LABEL_DIGITS, "f"), mean_val)),
    inherit.aes = FALSE,
    hjust = 1.02,
    vjust = -0.2,
    size = FIG_ANNOT_TEXT_SIZE,
    color = "black",
    check_overlap = FIG_LABEL_CHECK_OVERLAP,
    na.rm = TRUE
  ) +
  facet_grid(metric ~ site, scales = "free_y",
             labeller = labeller(metric = metric_labels_panel, site = site_labels_panel), drop = FALSE) +
  scale_color_manual(values = site_colors) +
  labs(x = "Water Year", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = FIG_AXIS_TEXT_SIZE),
    strip.text.y = element_text(angle = 0, size = FIG_STRIP_TEXT_SIZE),
    strip.text.x = element_text(size = FIG_STRIP_TEXT_SIZE),
    plot.margin = margin(
      FIG_LABEL_PLOT_MARGIN_PT,
      FIG_LABEL_PLOT_MARGIN_PT,
      FIG_LABEL_PLOT_MARGIN_PT,
      FIG_LABEL_PLOT_MARGIN_PT
    )
  ) +
  coord_cartesian(clip = FIG_LABEL_CLIP)

ggsave(file.path(supp_dir, "ds_annual_ts.png"), supp_dynamic_ts, width = 14 * FIG_WIDTH_SCALE, height = 10 * FIG_HEIGHT_SCALE, dpi = 300)
ggsave(file.path(supp_dir, "ds_annual_ts.pdf"), supp_dynamic_ts, width = 14 * FIG_WIDTH_SCALE, height = 10 * FIG_HEIGHT_SCALE)

# -----------------------------------------------------------------------------
# SUPPLEMENT: FACETED TIME SERIES - CHS
# -----------------------------------------------------------------------------

if (!is.null(chs_data) && nrow(chs_data) > 0) {
  cat("Creating Supplement: CHS Time Series...\n")

  chs_means <- chs_data %>%
    group_by(site) %>%
    summarise(
      mean_val = mean(CHS, na.rm = TRUE),
      x_pos = max(year, na.rm = TRUE),
      .groups = "drop"
    )
  chs_line_data <- chs_data %>%
    group_by(site) %>%
    filter(sum(!is.na(CHS)) >= 2) %>%
    ungroup()

  supp_chs_ts <- ggplot(chs_data, aes(x = year, y = CHS, color = site)) +
    geom_line(data = chs_line_data, linewidth = 0.6, na.rm = TRUE) +
    geom_point(size = FIG_POINT_SIZE_SMALL, na.rm = TRUE) +
    geom_hline(
      data = chs_means,
      aes(yintercept = mean_val, linetype = site),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.35,
      na.rm = TRUE
    ) +
    geom_text(
      data = chs_means,
      aes(x = x_pos, y = mean_val, label = sprintf(paste0("%.", FIG_MEAN_LABEL_DIGITS, "f"), mean_val)),
      inherit.aes = FALSE,
      hjust = 1.02,
      vjust = -0.2,
      size = FIG_ANNOT_TEXT_SIZE,
      color = "black",
      check_overlap = FIG_LABEL_CHECK_OVERLAP,
      na.rm = TRUE
    ) +
    facet_wrap(~site, ncol = 2, drop = FALSE, axes = "margins", axis.labels = "margins",
               labeller = labeller(site = site_labels_panel)) +
    scale_color_manual(values = site_colors) +
    scale_linetype_manual(values = FIG_MEAN_LINE_LINETYPES, guide = "none") +
    labs(x = "Water Year", y = "Baseflow Fraction (CHS)") +
    theme(
      plot.margin = margin(
        FIG_LABEL_PLOT_MARGIN_PT,
        FIG_LABEL_PLOT_MARGIN_PT,
        FIG_LABEL_PLOT_MARGIN_PT,
        FIG_LABEL_PLOT_MARGIN_PT
      )
    ) +
    coord_cartesian(clip = FIG_LABEL_CLIP)

  ggsave(file.path(supp_dir, "ms_chs_annual_ts.png"), supp_chs_ts, width = 10 * FIG_WIDTH_SCALE, height = 10 * FIG_HEIGHT_SCALE, dpi = 300)
  ggsave(file.path(supp_dir, "ms_chs_annual_ts.pdf"), supp_chs_ts, width = 10 * FIG_WIDTH_SCALE, height = 10 * FIG_HEIGHT_SCALE)
  cat("  Saved CHS Time Series\n")
}

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------

for (f in list.files(main_dir, pattern = "\\.(png|pdf)$")) {
  cat("  -", f, "\n")
}

for (f in list.files(supp_dir, pattern = "\\.(png|pdf)$")) {
  cat("  -", f, "\n")
}
