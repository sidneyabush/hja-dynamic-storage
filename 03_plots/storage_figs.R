# Combined storage figure generation (ANOVA/Tukey + dynamic/mobile/CHS figures).
# Inputs: OUTPUT_DIR/master/master_annual.csv; OUT_STATS_ANOVA_DIR/anova_results.csv; OUT_STATS_ANOVA_DIR/tukey_group_letters.csv; OUT_MET_MOBILE_DIR/isotope_metrics_site.csv; OUT_MET_MOBILE_DIR/annual_gw_prop.csv; ISOTOPE_DIR/MTT_FYW.csv; ISOTOPE_DIR/DampingRatios_2025-07-07.csv.
# Author: Sidney Bush
# Date: 2026-02-14

source("config.R")

# ---- Block 1: ANOVA/Tukey figures ----
local({
  # ANOVA/Tukey Plots.
  # Inputs: ISOTOPE_DIR/MTT_FYW.csv; ISOTOPE_DIR/DampingRatios_2025-07-07.csv; OUT_STATS_ANOVA_DIR/anova_results.csv.
  # Author: Sidney Bush
  # Date: 2026-02-13

  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)

  rm(list = ls())

  # Load project config

  main_dir <- file.path(FIGURES_DIR, "main")
  main_pdf_dir <- file.path(main_dir, "pdf")
  supp_dir <- file.path(FIGURES_DIR, "supp")
  supp_pdf_dir <- file.path(supp_dir, "pdf")
  dir_targets <- c(main_dir, main_pdf_dir, supp_dir, supp_pdf_dir)
  if (isTRUE(WRITE_TABLE_OUTPUTS)) {
    table_dir <- file.path(OUT_TABLES_DIR, "anova_tukey")
    dir_targets <- c(dir_targets, table_dir)
  }
  for (d in dir_targets) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  safe_ggsave <- function(filename, plot, width, height, dpi = NULL) {
    tryCatch(
      {
        if (is.null(dpi)) {
          ggsave(filename, plot, width = width, height = height)
        } else {
          ggsave(filename, plot, width = width, height = height, dpi = dpi)
        }
        TRUE
      },
      error = function(e) {
        warning(sprintf(
          "Could not write %s (%s)",
          filename,
          conditionMessage(e)
        ))
        FALSE
      }
    )
  }

  annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
  annual <- read_csv(annual_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

  # For ANOVA/Tukey visualization only, use annual FDC (site-year) values.
  fdc_wy_file <- file.path(OUT_MET_DYNAMIC_DIR, "fdc_slopes_wy.csv")
  if (file.exists(fdc_wy_file)) {
    fdc_wy <- read_csv(fdc_wy_file, show_col_types = FALSE)
    if ("WaterYear" %in% names(fdc_wy) && !("year" %in% names(fdc_wy))) {
      fdc_wy <- fdc_wy %>% rename(year = WaterYear)
    }
    if ("Slope" %in% names(fdc_wy) && !("FDC_anova" %in% names(fdc_wy))) {
      fdc_wy <- fdc_wy %>% rename(FDC_anova = Slope)
    }
    fdc_wy <- fdc_wy %>%
      transmute(
        site = standardize_site_code(site),
        year = as.integer(year),
        FDC_anova = as.numeric(FDC_anova)
      ) %>%
      filter(
        site %in% SITE_ORDER_HYDROMETRIC,
        year >= WY_START,
        year <= WY_END
      )

    annual <- annual %>%
      mutate(year = as.integer(year)) %>%
      left_join(fdc_wy, by = c("site", "year")) %>%
      mutate(FDC = dplyr::coalesce(FDC_anova, FDC)) %>%
      select(-FDC_anova)
  }

  # Isotope predictors are site-level (not annual) and only used in the all-metrics panel.
  iso_mtt_file <- file.path(ISOTOPE_DIR, "MTT_FYW.csv")
  iso_dr_file <- file.path(ISOTOPE_DIR, "DampingRatios_2025-07-07.csv")

  iso_mtt <- if (file.exists(iso_mtt_file)) {
    read_csv(iso_mtt_file, show_col_types = FALSE) %>%
      mutate(site = standardize_site_code(site)) %>%
      transmute(
        site,
        mtt1 = MTT1,
        mtt1_se = MTT1_SD,
        mtt2 = MTTM,
        mtt2_se = (MTT2L_SD + MTT2H_SD) / 2,
        fyw = FYWM,
        fyw_se = (FYWL_SD + FYWH_SD) / 2
      )
  } else {
    tibble(
      site = character(),
      mtt1 = numeric(),
      mtt1_se = numeric(),
      mtt2 = numeric(),
      mtt2_se = numeric(),
      fyw = numeric(),
      fyw_se = numeric()
    )
  }

  iso_dr <- if (file.exists(iso_dr_file)) {
    read_csv(iso_dr_file, show_col_types = FALSE) %>%
      mutate(site = standardize_site_code(site)) %>%
      transmute(site, dr = DR_Overall, dr_se = DR__err)
  } else {
    tibble(site = character(), dr = numeric(), dr_se = numeric())
  }

  iso_site <- full_join(iso_mtt, iso_dr, by = "site") %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

  letters_file <- file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv")
  letters_df <- if (file.exists(letters_file)) {
    read_csv(letters_file, show_col_types = FALSE) %>%
      mutate(site = standardize_site_code(site))
  } else {
    tibble(metric = character(), site = character(), group_letter = character())
  }

  anova_file <- file.path(OUT_STATS_ANOVA_DIR, "anova_results.csv")
  anova_df <- if (file.exists(anova_file)) {
    read_csv(anova_file, show_col_types = FALSE)
  } else {
    tibble(metric = character(), p_value = numeric())
  }

  metrics_main <- PLOT_ORDER_DYNAMIC_STORAGE
  metric_display_labels <- c(
    RBI = "RBI",
    RCS = "RCS",
    FDC = "FDC",
    SD = "SD",
    WB = "WB drawdown (mm)",
    CHS = "CHS",
    MTT1 = "MTT1",
    MTT2 = "MTT2",
    DR = "DR",
    FYW = "FYW"
  )

  build_faceted_plot <- function(df_long, metric_names, ncol_facets) {
    d <- df_long %>%
      filter(metric %in% metric_names) %>%
      mutate(value = ifelse(metric == "WB", -value, value)) %>%
      mutate(metric = factor(metric, levels = metric_names))
    metric_titles <- ifelse(
      metric_names %in% names(metric_display_labels),
      unname(metric_display_labels[metric_names]),
      metric_names
    )
    metric_labels_panel <- stats::setNames(
      paste0(letters[seq_along(metric_names)], ") ", metric_titles),
      metric_names
    )

    letters_m <- letters_df %>%
      filter(metric %in% metric_names) %>%
      mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
      left_join(
        d %>%
          group_by(metric, site) %>%
          summarise(
            y_site_max = suppressWarnings(max(value, na.rm = TRUE)),
            .groups = "drop"
          ) %>%
          left_join(
            d %>%
              group_by(metric) %>%
              summarise(
                y_span = {
                  y_max <- suppressWarnings(max(value, na.rm = TRUE))
                  y_min <- suppressWarnings(min(value, na.rm = TRUE))
                  ifelse(
                    is.finite(y_max - y_min) && (y_max - y_min) > 0,
                    y_max - y_min,
                    1
                  )
                },
                .groups = "drop"
              ),
            by = "metric"
          ),
        by = c("metric", "site")
      ) %>%
      mutate(
        y_label = ifelse(
          is.finite(y_site_max),
          y_site_max + 0.04 * y_span,
          NA_real_
        )
      ) %>%
      mutate(metric = factor(metric, levels = metric_names))

    ggplot(d, aes(x = site, y = value, color = site, fill = site)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.2, na.rm = TRUE) +
      geom_point(
        position = position_jitter(width = 0.12, height = 0),
        size = FIG_POINT_SIZE_SMALL,
        alpha = 0.7,
        na.rm = TRUE
      ) +
      geom_text(
        data = letters_m,
        aes(x = site, y = y_label, label = group_letter),
        inherit.aes = FALSE,
        color = "black",
        size = FIG_ANNOT_TEXT_SIZE,
        check_overlap = FIG_LABEL_CHECK_OVERLAP,
        na.rm = TRUE
      ) +
      scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC, drop = FALSE) +
      scale_color_manual(values = SITE_COLORS) +
      scale_fill_manual(values = SITE_COLORS) +
      facet_wrap(
        ~metric,
        ncol = ncol_facets,
        scales = "free_y",
        labeller = labeller(metric = metric_labels_panel),
        axes = "margins",
        axis.labels = "margins"
      ) +
      labs(x = NULL, y = "Metric value (WB as drawdown)") +
      theme_pub() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
        axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
        strip.background = element_blank(),
        strip.text = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0),
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.4
        ),
        axis.line = element_blank(),
        plot.margin = margin(
          FIG_LABEL_PLOT_MARGIN_PT,
          FIG_LABEL_PLOT_MARGIN_PT,
          FIG_LABEL_PLOT_MARGIN_PT,
          FIG_LABEL_PLOT_MARGIN_PT
        )
      ) +
      coord_cartesian(clip = FIG_LABEL_CLIP)
  }

  build_mean_se_plot <- function(df_long, metric_names, ncol_facets) {
    d <- df_long %>%
      filter(metric %in% metric_names) %>%
      mutate(metric = factor(metric, levels = metric_names))
    d <- d %>%
      mutate(mean_value = ifelse(metric == "WB", -mean_value, mean_value))
    metric_titles <- ifelse(
      metric_names %in% names(metric_display_labels),
      unname(metric_display_labels[metric_names]),
      metric_names
    )
    metric_labels_panel <- stats::setNames(
      paste0(letters[seq_along(metric_names)], ") ", metric_titles),
      metric_names
    )

    letters_m <- letters_df %>%
      mutate(metric = toupper(metric)) %>%
      filter(metric %in% metric_names) %>%
      mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
      left_join(
        d %>%
          group_by(metric, site) %>%
          summarise(
            y_site_max = suppressWarnings(max(
              mean_value + se_value,
              na.rm = TRUE
            )),
            .groups = "drop"
          ) %>%
          left_join(
            d %>%
              group_by(metric) %>%
              summarise(
                y_span = {
                  y_max <- suppressWarnings(max(
                    mean_value + se_value,
                    na.rm = TRUE
                  ))
                  y_min <- suppressWarnings(min(
                    mean_value - se_value,
                    na.rm = TRUE
                  ))
                  ifelse(
                    is.finite(y_max - y_min) && (y_max - y_min) > 0,
                    y_max - y_min,
                    1
                  )
                },
                .groups = "drop"
              ),
            by = "metric"
          ),
        by = c("metric", "site")
      ) %>%
      mutate(
        y_label = ifelse(
          is.finite(y_site_max),
          y_site_max + 0.04 * y_span,
          NA_real_
        )
      ) %>%
      mutate(metric = factor(metric, levels = metric_names))

    ggplot(d, aes(x = site, y = mean_value, color = site)) +
      geom_point(size = FIG_POINT_SIZE_MED, na.rm = TRUE) +
      geom_errorbar(
        aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
        width = 0.3,
        linewidth = 0.5,
        na.rm = TRUE
      ) +
      geom_text(
        data = letters_m,
        aes(x = site, y = y_label, label = group_letter),
        inherit.aes = FALSE,
        color = "black",
        size = FIG_ANNOT_TEXT_SIZE,
        check_overlap = FIG_LABEL_CHECK_OVERLAP,
        na.rm = TRUE
      ) +
      scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC, drop = FALSE) +
      scale_color_manual(values = SITE_COLORS) +
      facet_wrap(
        ~metric,
        ncol = ncol_facets,
        scales = "free_y",
        labeller = labeller(metric = metric_labels_panel),
        axes = "margins",
        axis.labels = "margins"
      ) +
      labs(x = NULL, y = "Value (WB as drawdown)") +
      theme_pub() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
        axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
        strip.background = element_blank(),
        strip.text = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0),
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.4
        ),
        axis.line = element_blank(),
        plot.margin = margin(
          FIG_LABEL_PLOT_MARGIN_PT,
          FIG_LABEL_PLOT_MARGIN_PT,
          FIG_LABEL_PLOT_MARGIN_PT,
          FIG_LABEL_PLOT_MARGIN_PT
        )
      ) +
      coord_cartesian(clip = FIG_LABEL_CLIP)
  }

  build_single_metric_plot <- function(df_long, metric_name, y_label, panel_title) {
    d <- df_long %>%
      filter(metric == metric_name) %>%
      mutate(metric = factor(metric, levels = metric_name))
    if (metric_name == "WB") {
      # Plot depletion magnitude as positive for readability.
      d <- d %>% mutate(value = -value)
    }

    letters_m <- letters_df %>%
      filter(metric == metric_name) %>%
      mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
      left_join(
        d %>%
          group_by(site) %>%
          summarise(
            y_site_max = suppressWarnings(max(value, na.rm = TRUE)),
            .groups = "drop"
          ),
        by = "site"
      )

    y_max <- suppressWarnings(max(d$value, na.rm = TRUE))
    y_min <- suppressWarnings(min(d$value, na.rm = TRUE))
    y_span <- ifelse(is.finite(y_max - y_min) && (y_max - y_min) > 0, y_max - y_min, 1)
    letters_m <- letters_m %>%
      mutate(y_label = ifelse(is.finite(y_site_max), y_site_max + 0.04 * y_span, NA_real_))

    ggplot(d, aes(x = site, y = value, color = site, fill = site)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.2, na.rm = TRUE) +
      geom_point(
        position = position_jitter(width = 0.12, height = 0),
        size = FIG_POINT_SIZE_SMALL,
        alpha = 0.7,
        na.rm = TRUE
      ) +
      geom_text(
        data = letters_m,
        aes(x = site, y = y_label, label = group_letter),
        inherit.aes = FALSE,
        color = "black",
        size = FIG_ANNOT_TEXT_SIZE,
        check_overlap = FIG_LABEL_CHECK_OVERLAP,
        na.rm = TRUE
      ) +
      scale_x_discrete(limits = SITE_ORDER_HYDROMETRIC, drop = FALSE) +
      scale_color_manual(values = SITE_COLORS) +
      scale_fill_manual(values = SITE_COLORS) +
      labs(x = NULL, y = y_label, title = panel_title) +
      theme_pub() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
        axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
        plot.title = element_text(
          size = FIG_AXIS_TITLE_SIZE,
          hjust = 0,
          margin = margin(b = 10)
        ),
        plot.title.position = "plot",
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.4
        ),
        axis.line = element_blank(),
        plot.margin = margin(
          FIG_LABEL_PLOT_MARGIN_PT,
          FIG_LABEL_PLOT_MARGIN_PT,
          FIG_LABEL_PLOT_MARGIN_PT,
          FIG_LABEL_PLOT_MARGIN_PT
        )
      ) +
      coord_cartesian(clip = FIG_LABEL_CLIP)
  }

  annual_long <- annual %>%
    select(site, all_of(unique(c(metrics_main, "CHS", "WB")))) %>%
    pivot_longer(cols = -site, names_to = "metric", values_to = "value")

  iso_long <- iso_site %>%
    select(site, mtt1, mtt2, dr, fyw) %>%
    pivot_longer(
      cols = c(mtt1, mtt2, dr, fyw),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(metric = toupper(metric))

  annual_summary <- annual_long %>%
    group_by(metric, site) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      se_value = sd(value, na.rm = TRUE) / sqrt(sum(is.finite(value))),
      .groups = "drop"
    ) %>%
    mutate(metric = toupper(metric))

  iso_se_long <- iso_site %>%
    select(site, mtt1_se, mtt2_se, dr_se, fyw_se) %>%
    pivot_longer(
      cols = -site,
      names_to = "metric_se",
      values_to = "se_value"
    ) %>%
    mutate(metric = toupper(gsub("_se$", "", metric_se))) %>%
    select(site, metric, se_value)

  iso_summary <- iso_long %>%
    left_join(iso_se_long, by = c("site", "metric")) %>%
    transmute(metric, site, mean_value = value, se_value = se_value)

  all_summary <- bind_rows(annual_summary, iso_summary)

  metric_units_main <- c(
    RBI = "Unitless",
    RCS = "Unitless",
    FDC = "Unitless",
    SD = "Depth (mm)",
    WB = "Drawdown (mm)"
  )
  metric_titles_main <- setNames(
    paste0(
      letters[seq_along(metrics_main)],
      ") ",
      c("RBI", "RCS", "FDC", "SD", "WB")
    ),
    metrics_main
  )
  fig_main_list <- lapply(metrics_main, function(m) {
    build_single_metric_plot(
      annual_long,
      metric_name = m,
      y_label = metric_units_main[[m]],
      panel_title = metric_titles_main[[m]]
    )
  })
  fig_main <- patchwork::wrap_plots(fig_main_list, ncol = 2)
  safe_ggsave(
    file.path(main_dir, "ds_anova_tukey.png"),
    fig_main,
    width = 11 * FIG_WIDTH_SCALE,
    height = 12 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
  safe_ggsave(
    file.path(main_pdf_dir, "ds_anova_tukey.pdf"),
    fig_main,
    width = 11 * FIG_WIDTH_SCALE,
    height = 12 * FIG_HEIGHT_SCALE
  )

  # Flat table for manuscript/supplement use (optional)
  if (isTRUE(WRITE_TABLE_OUTPUTS) && nrow(anova_df) > 0 && nrow(letters_df) > 0) {
    if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    summary_table <- letters_df %>%
      left_join(
        anova_df %>% select(metric, F_statistic, p_value),
        by = "metric"
      ) %>%
      arrange(metric, site)
    write_csv(
      summary_table,
      file.path(table_dir, "anova_tukey_group_table.csv")
    )
  }
})

# ---- Block 2: Dynamic/Mobile/CHS figures ----
local({
  # Main Text and Supplement Figures - HJA Dynamic Storage.
  # Inputs: OUT_MET_MOBILE_DIR/isotope_metrics_site.csv; ISOTOPE_DIR/MTT_FYW.csv; OUT_MET_MOBILE_DIR/Annual_GW_Prop.csv.
  # Author: Sidney Bush
  # Date: 2026-02-02

  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
  library(colorspace)

  rm(list = ls())

  # SOURCE CONFIGURATION

  # Load project config

  # SETUP

  base_dir <- BASE_DATA_DIR
  main_dir <- file.path(FIGURES_DIR, "main")
  main_pdf_dir <- file.path(main_dir, "pdf")
  supp_dir <- file.path(FIGURES_DIR, "supp")
  supp_pdf_dir <- file.path(supp_dir, "pdf")

  for (d in c(main_dir, main_pdf_dir, supp_dir, supp_pdf_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  # Remove legacy supplement figures that are no longer generated.
  unlink(file.path(supp_dir, c(
    "ds_summary.png",
    "ms_chs.png",
    "storage_anova_tukey_all.png"
  )))
  unlink(file.path(supp_pdf_dir, c(
    "ds_summary.pdf",
    "ms_chs.pdf",
    "storage_anova_tukey_all.pdf"
  )))

  site_order <- SITE_ORDER_HYDROMETRIC
  site_order_iso <- SITE_ORDER_ALL
  # EC/CHS plots: show all hydrometric sites; sites without CHS appear as blanks.
  CHS_SHOW_ALL_HYDRO_SITES <- TRUE
  # Keep these sites on axes, but do not plot CHS points/lines for them.
  CHS_EXCLUDE_DATA_SITES <- c("Look")
  site_order_chs <- if (CHS_SHOW_ALL_HYDRO_SITES) {
    SITE_ORDER_HYDROMETRIC
  } else {
    SITE_ORDER_CHEMISTRY
  }
  site_colors <- SITE_COLORS
  site_labels_panel <- make_panel_label_map(site_order)
  site_labels_panel_chs <- make_panel_label_map(site_order_chs)

  # Metric labels
  metric_labels <- c(
    "RBI" = "RBI",
    "RCS" = "RCS",
    "FDC" = "FDC",
    "SD" = "SD (mm)",
    "WB" = "WB drawdown (mm)",
    "CHS" = "CHS",
    "MTT1" = "MTT1 (yr)",
    "MTT2" = "MTT2 (yr)",
    "Fyw" = "Fyw",
    "DR" = "DR"
  )

  # Publication theme
  theme_pub_base <- theme_pub() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = FIG_AXIS_TEXT_SIZE
      ),
      axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      strip.background = element_blank(),
      strip.text = element_text(size = FIG_STRIP_TEXT_SIZE, hjust = 0),
      legend.position = "none"
    )

  theme_set(theme_pub_base)

  # HELPER: Standardize site names to WS## format

  standardize_sites <- function(df, allowed_sites = site_order) {
    df %>%
      mutate(
        site = standardize_site_code(site)
      ) %>%
      filter(site %in% allowed_sites) %>%
      mutate(site = factor(site, levels = allowed_sites))
  }

  # LOAD DATA

  # Annual storage metrics
  annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)

  annual_data <- read_csv(annual_file, show_col_types = FALSE) %>%
    standardize_sites(allowed_sites = site_order)

  letters_file <- file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv")
  letters_df <- if (file.exists(letters_file)) {
    read_csv(letters_file, show_col_types = FALSE) %>%
      mutate(site = standardize_site_code(site))
  } else {
    tibble(metric = character(), site = character(), group_letter = character())
  }

  # Isotope data for mobile storage (site-level metrics)
  isotope_file <- file.path(OUT_MET_MOBILE_DIR, "isotope_metrics_site.csv")

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
      mtt_raw$MTT1_high <- first_col(
        mtt_raw,
        c("MTT1H", "MTT1_high", "MTT1_MAX")
      )
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
      mtt_raw$MTT2_high <- first_col(
        mtt_raw,
        c("MTT2H", "MTT2_high", "MTT2_MAX")
      )
      mtt_raw$MTT2_low_sd <- first_col(mtt_raw, c("MTT2L_SD", "MTT2_low_SD"))
      mtt_raw$MTT2_high_sd <- first_col(mtt_raw, c("MTT2H_SD", "MTT2_high_SD"))
      mtt_raw$MTT2_err <- ifelse(
        is.finite(mtt_raw$MTT2_low_sd) & is.finite(mtt_raw$MTT2_high_sd),
        abs((mtt_raw$MTT2_low_sd + mtt_raw$MTT2_high_sd) / 2),
        abs((mtt_raw$MTT2_high - mtt_raw$MTT2_low) / 2)
      )

      mtt_raw$Fyw_low <- first_col(
        mtt_raw,
        c("FywL", "FYWL", "Fyw_low", "FYWL")
      )
      mtt_raw$Fyw_high <- first_col(
        mtt_raw,
        c("FywH", "FYWH", "Fyw_high", "FYWH")
      )
      mtt_raw$Fyw_val <- first_col(mtt_raw, c("Fyw", "FYW", "fyw", "FYWM"))
      if (all(is.na(mtt_raw$Fyw_val))) {
        mtt_raw$Fyw_val <- rowMeans(
          cbind(mtt_raw$Fyw_low, mtt_raw$Fyw_high),
          na.rm = TRUE
        )
      }
      mtt_raw$Fyw_low_sd <- first_col(
        mtt_raw,
        c("FYWL_SD", "FywL_SD", "Fyw_low_SD")
      )
      mtt_raw$Fyw_high_sd <- first_col(
        mtt_raw,
        c("FYWH_SD", "FywH_SD", "Fyw_high_SD")
      )
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
          df[[nm]] <- dplyr::coalesce(
            suppressWarnings(as.numeric(xv)),
            suppressWarnings(as.numeric(yv))
          )
          if (x %in% names(df)) {
            df[[x]] <- NULL
          }
          if (y %in% names(df)) df[[y]] <- NULL
        }
        df
      }
      for (nm in c(
        "MTT1",
        "MTT1_err",
        "MTT2",
        "MTT2_err",
        "Fyw",
        "Fyw_err",
        "DR",
        "DR_err"
      )) {
        isotope_data <- merge_col(isotope_data, nm)
      }
    }

    if (!("MTT1" %in% names(isotope_data))) {
      isotope_data$MTT1 <- NA_real_
    }
    if (!("MTT1_err" %in% names(isotope_data))) {
      isotope_data$MTT1_err <- NA_real_
    }
    if (!("MTT2" %in% names(isotope_data))) {
      isotope_data$MTT2 <- NA_real_
    }
    if (!("MTT2_err" %in% names(isotope_data))) {
      isotope_data$MTT2_err <- NA_real_
    }
    if (!("Fyw" %in% names(isotope_data))) {
      isotope_data$Fyw <- NA_real_
    }
    if (!("Fyw_err" %in% names(isotope_data))) {
      fyw_low_name <- names(isotope_data)[
        tolower(names(isotope_data)) %in% c("fywl", "fyw_low")
      ]
      fyw_high_name <- names(isotope_data)[
        tolower(names(isotope_data)) %in% c("fywh", "fyw_high")
      ]
      if (length(fyw_low_name) > 0 && length(fyw_high_name) > 0) {
        isotope_data$Fyw_err <- abs(
          (isotope_data[[fyw_high_name[1]]] - isotope_data[[fyw_low_name[1]]]) /
            2
        )
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

    if (isTRUE(WRITE_AUX_OUTPUTS)) {
      write_csv(
        isotope_plot_values,
        file.path(OUT_MET_SUPPORT_DIR, "ms_isotope_plot_values.csv")
      )
    }
    cat("  Loaded isotope data:", nrow(isotope_data), "rows\n")
  } else {
    isotope_data <- NULL
    cat("  Warning: Isotope data not found\n")
  }

  # CHS data
  chs_file <- file.path(OUT_MET_MOBILE_DIR, "annual_gw_prop.csv")

  if (file.exists(chs_file)) {
    chs_data <- read_csv(chs_file, show_col_types = FALSE) %>%
      mutate(
        site = if ("site" %in% names(.)) {
          site
        } else if ("SITECODE" %in% names(.)) {
          SITECODE
        } else {
          NA_character_
        },
        year = if ("year" %in% names(.)) {
          year
        } else if ("waterYear" %in% names(.)) {
          waterYear
        } else {
          NA_real_
        },
        CHS = if ("CHS" %in% names(.)) {
          CHS
        } else if ("mean_bf" %in% names(.)) {
          mean_bf
        } else {
          CHS
        }
      ) %>%
      select(site, year, CHS) %>%
      standardize_sites(allowed_sites = site_order_chs) %>%
      filter(!(as.character(site) %in% CHS_EXCLUDE_DATA_SITES))
    cat("  Loaded CHS data:", nrow(chs_data), "rows\n")
  } else {
    # Try to get CHS from annual data
    if ("CHS" %in% names(annual_data)) {
      chs_data <- annual_data %>%
        select(site, year, CHS) %>%
        filter(!is.na(CHS), as.character(site) %in% site_order_chs) %>%
        mutate(site = factor(as.character(site), levels = site_order_chs)) %>%
        filter(!(as.character(site) %in% CHS_EXCLUDE_DATA_SITES))
      cat("  Using CHS from annual data:", nrow(chs_data), "rows\n")
    } else {
      chs_data <- NULL
      cat("  Warning: CHS data not found\n")
    }
  }

  # Dynamic storage annual long table used by downstream supplement plots.

  dynamic_metrics_target <- c("RBI", "RCS", "FDC", "SD", "WB")
  dynamic_metrics <- dynamic_metrics_target[dynamic_metrics_target %in% names(annual_data)]
  missing_dynamic_metrics <- setdiff(dynamic_metrics_target, dynamic_metrics)
  if (length(missing_dynamic_metrics) > 0) {
    warning(
      "Missing expected dynamic metrics in annual data: ",
      paste(missing_dynamic_metrics, collapse = ", ")
    )
  }
  if (!("FDC" %in% dynamic_metrics)) {
    stop("FDC column is required for ds_annual_ts.png but was not found in annual data.")
  }
  dynamic_metric_titles <- c(
    RBI = "RBI",
    RCS = "RCS",
    FDC = "FDC",
    SD = "SD",
    WB = "WB drawdown"
  )
  metric_labels_panel <- make_panel_label_map(dynamic_metric_titles[
    dynamic_metrics
  ])

  # Keep WB on depletion-magnitude convention across all figures.
  dynamic_long <- annual_data %>%
    select(site, year, any_of(dynamic_metrics)) %>%
    pivot_longer(
      cols = -c(site, year),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(value = ifelse(metric == "WB", -value, value)) %>%
    mutate(metric = factor(metric, levels = dynamic_metrics_target)) %>%
    mutate(
      metric_label = factor(
        as.character(metric),
        levels = dynamic_metrics_target,
        labels = dynamic_metric_titles[dynamic_metrics_target]
      )
    )

  # FIGURE 4: MOBILE ISOTOPE (MTT, Fyw, DR)

  if (!is.null(isotope_data) && nrow(isotope_data) > 0) {
    metric_order_iso <- c(
      unname(PLOT_MOBILE_STORAGE_SITE_COLS[PLOT_ORDER_MOBILE_STORAGE]),
      "MTT2"
    )
    metric_order_iso <- unique(metric_order_iso)
    metric_order_iso <- metric_order_iso[
      metric_order_iso %in% c("DR", "Fyw", "MTT1", "MTT2")
    ]

    isotope_values_long <- isotope_data %>%
      select(site, all_of(metric_order_iso)) %>%
      pivot_longer(cols = -site, names_to = "metric", values_to = "value")

    isotope_err_long <- isotope_data %>%
      transmute(
        site,
        MTT1 = MTT1_err,
        MTT2 = MTT2_err,
        Fyw = Fyw_err,
        DR = DR_err
      ) %>%
      pivot_longer(cols = -site, names_to = "metric", values_to = "err")

    isotope_plot_long <- isotope_values_long %>%
      left_join(isotope_err_long, by = c("site", "metric")) %>%
      mutate(metric = factor(metric, levels = metric_order_iso))

    metric_titles_iso <- c(
      "DR" = "DR",
      "Fyw" = "Fyw",
      "MTT1" = "MTT1",
      "MTT2" = "MTT2"
    )
    metric_labels_iso <- setNames(
      paste0(letters[seq_along(metric_order_iso)], ") ", metric_titles_iso[metric_order_iso]),
      metric_order_iso
    )
    base_iso_plot <- ggplot(
      isotope_plot_long,
      aes(x = site, y = value, color = site)
    ) +
      geom_point(size = FIG_POINT_SIZE_LARGE, na.rm = TRUE) +
      geom_errorbar(
        aes(ymin = value - err, ymax = value + err),
        width = 0.3,
        linewidth = 0.5,
        na.rm = TRUE
      ) +
      scale_color_manual(values = site_colors, guide = "none") +
      scale_x_discrete(limits = site_order_iso, drop = FALSE) +
      labs(x = NULL, y = NULL)

    # Keep MTT1 and MTT2 on the same y-axis scale for direct comparison.
    mtt_long <- isotope_plot_long %>% filter(metric %in% c("MTT1", "MTT2"))
    mtt_ylim <- mtt_long %>%
      mutate(
        ymin = value - ifelse(is.na(err), 0, err),
        ymax = value + ifelse(is.na(err), 0, err)
      ) %>%
      summarise(
        lo = min(ymin, na.rm = TRUE),
        hi = max(ymax, na.rm = TRUE)
      )
    if (!is.finite(mtt_ylim$lo) || !is.finite(mtt_ylim$hi)) {
      mtt_ylim$lo <- NA_real_
      mtt_ylim$hi <- NA_real_
    }

    make_iso_panel <- function(metric_name, y_label, y_limits = NULL) {
      panel_data <- isotope_plot_long %>% filter(metric == metric_name)
      p <- (base_iso_plot %+% panel_data) +
        labs(
          y = y_label,
          title = metric_labels_iso[[metric_name]]
        ) +
        theme(
          plot.title = element_text(
            size = FIG_STRIP_TEXT_SIZE + 2,
            hjust = 0,
            margin = margin(b = 10)
          ),
          plot.title.position = "plot"
        )
      if (!is.null(y_limits) && all(is.finite(y_limits))) {
        p <- p + coord_cartesian(ylim = y_limits)
      }
      p
    }

    fig4_dr <- make_iso_panel("DR", "Ratio") +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    fig4_fyw <- make_iso_panel("Fyw", "Fraction") +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    mtt_limits <- c(mtt_ylim$lo, mtt_ylim$hi)
    fig4_mtt1 <- make_iso_panel("MTT1", "Years", mtt_limits)
    fig4_mtt2 <- make_iso_panel("MTT2", "Years", mtt_limits)

    fig4 <- (fig4_dr + fig4_fyw) / (fig4_mtt1 + fig4_mtt2) +
      patchwork::plot_layout(heights = c(1, 1))

    ggsave(
      file.path(main_dir, "ms_isotope.png"),
      fig4,
      width = 9 * FIG_WIDTH_SCALE,
      height = 8 * FIG_HEIGHT_SCALE,
      dpi = 300
    )
    ggsave(
      file.path(main_pdf_dir, "ms_isotope.pdf"),
      fig4,
      width = 9 * FIG_WIDTH_SCALE,
      height = 8 * FIG_HEIGHT_SCALE
    )
    cat("  Saved Figure 4\n")

    # FIGURE 4B: MTT1 vs MTT2 comparison (sites with both values)
    mtt_compare <- isotope_data %>%
      transmute(site, MTT1 = as.numeric(MTT1), MTT2 = as.numeric(MTT2)) %>%
      filter(is.finite(MTT1), is.finite(MTT2))

    if (nrow(mtt_compare) > 0) {
      lim_lo <- min(c(mtt_compare$MTT1, mtt_compare$MTT2), na.rm = TRUE)
      lim_hi <- max(c(mtt_compare$MTT1, mtt_compare$MTT2), na.rm = TRUE)
      lim_pad <- 0.08 * (lim_hi - lim_lo)
      if (!is.finite(lim_pad) || lim_pad <= 0) {
        lim_pad <- 0.2
      }
      n_overlap <- nrow(mtt_compare)
      mtt1_sd_overlap <- sd(mtt_compare$MTT1, na.rm = TRUE)
      mtt2_sd_overlap <- sd(mtt_compare$MTT2, na.rm = TRUE)

      fig4b <- ggplot(mtt_compare, aes(x = MTT1, y = MTT2, color = site)) +
        geom_abline(
          intercept = 0,
          slope = 1,
          linetype = "dashed",
          color = "gray40",
          linewidth = 0.6
        ) +
        geom_point(size = FIG_POINT_SIZE_LARGE) +
        geom_text(
          aes(label = site),
          nudge_x = 0.03,
          nudge_y = 0.03,
          size = FIG_ANNOT_TEXT_SIZE,
          check_overlap = TRUE,
          show.legend = FALSE
        ) +
        scale_color_manual(values = site_colors) +
        coord_equal(
          xlim = c(lim_lo - lim_pad, lim_hi + lim_pad),
          ylim = c(lim_lo - lim_pad, lim_hi + lim_pad),
          clip = "off"
        ) +
        labs(
          x = "MTT1 (yr)",
          y = "MTT2 (yr)",
          subtitle = sprintf(
            "Overlapping sites: n=%d | SD(MTT1)=%.2f | SD(MTT2)=%.2f",
            n_overlap,
            mtt1_sd_overlap,
            mtt2_sd_overlap
          )
        ) +
        theme(
          plot.margin = margin(
            FIG_LABEL_PLOT_MARGIN_PT,
            FIG_LABEL_PLOT_MARGIN_PT,
            FIG_LABEL_PLOT_MARGIN_PT,
            FIG_LABEL_PLOT_MARGIN_PT
          )
        )

      ggsave(
        file.path(supp_dir, "ms_mtt1_mtt2_compare.png"),
        fig4b,
        width = 7 * FIG_WIDTH_SCALE,
        height = 7 * FIG_HEIGHT_SCALE,
        dpi = 300
      )
      ggsave(
        file.path(supp_pdf_dir, "ms_mtt1_mtt2_compare.pdf"),
        fig4b,
        width = 7 * FIG_WIDTH_SCALE,
        height = 7 * FIG_HEIGHT_SCALE
      )
      cat("  Saved Figure 4B\n")
    } else {
      cat("  Skipping Figure 4B: No overlapping MTT1/MTT2 sites\n")
    }
  } else {
    cat("  Skipping Figure 4: No isotope data\n")
  }

  # FIGURE 5: CHS (BASEFLOW FRACTION)

  if (!is.null(chs_data) && nrow(chs_data) > 0) {
    chs_n_lookup <- chs_data %>%
      count(site, name = "n") %>%
      complete(
        site = factor(site_order_chs, levels = site_order_chs),
        fill = list(n = 0)
      ) %>%
      { setNames(.$n, as.character(.$site)) }
    chs_site_labels_n <- setNames(
      paste0(
        site_order_chs,
        "\n(n=",
        ifelse(
          is.na(chs_n_lookup[site_order_chs]),
          0,
          chs_n_lookup[site_order_chs]
        ),
        ")"
      ),
      site_order_chs
    )

    chs_data_box <- chs_data %>%
      filter(as.character(site) != "Look")

    fig5_box <- ggplot(
      chs_data_box,
      aes(x = site, y = CHS, color = site, fill = site)
    ) +
      geom_boxplot(outlier.shape = NA, alpha = 0.2, na.rm = TRUE) +
      geom_point(
        position = position_jitter(width = 0.12, height = 0),
        size = FIG_POINT_SIZE_SMALL,
        alpha = 0.7,
        na.rm = TRUE
      ) +
      scale_color_manual(values = site_colors) +
      scale_fill_manual(values = site_colors) +
      scale_x_discrete(
        limits = site_order_chs,
        drop = FALSE,
        labels = chs_site_labels_n
      ) +
      geom_text(
        data = {
          chs_box_letters <- letters_df %>%
            filter(metric == "CHS", site %in% site_order_chs) %>%
            mutate(site = factor(site, levels = site_order_chs)) %>%
            left_join(
              chs_data %>%
                group_by(site) %>%
                summarise(y_site_max = max(CHS, na.rm = TRUE), .groups = "drop"),
              by = "site"
            )
          y_max <- suppressWarnings(max(chs_data$CHS, na.rm = TRUE))
          y_min <- suppressWarnings(min(chs_data$CHS, na.rm = TRUE))
          y_span <- ifelse(is.finite(y_max - y_min) && (y_max - y_min) > 0, y_max - y_min, 1)
          chs_box_letters %>%
            mutate(y_label_box = ifelse(is.finite(y_site_max), y_site_max + 0.06 * y_span, NA_real_))
        },
        aes(x = site, y = y_label_box, label = group_letter),
        inherit.aes = FALSE,
        color = "black",
        size = FIG_ANNOT_TEXT_SIZE,
        vjust = 0
      ) +
      labs(x = NULL, y = "Baseflow Fraction (CHS)")

    ggsave(
      file.path(supp_dir, "ms_chs_boxplot.png"),
      fig5_box,
      width = 8 * FIG_WIDTH_SCALE,
      height = 5 * FIG_HEIGHT_SCALE,
      dpi = 300
    )
    ggsave(
      file.path(supp_pdf_dir, "ms_chs_boxplot.pdf"),
      fig5_box,
      width = 8 * FIG_WIDTH_SCALE,
      height = 5 * FIG_HEIGHT_SCALE
    )
    cat("  Saved Figure 5\n")
  } else {
    cat("  Skipping Figure 5: No CHS data\n")
  }

  # SUPPLEMENT: FACETED TIME SERIES - DYNAMIC STORAGE

  # Calculate site means for dashed reference lines
  site_means_dynamic <- dynamic_long %>%
    group_by(site, metric) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      x_pos = max(year, na.rm = TRUE),
      x_label = max(year, na.rm = TRUE) + 0.6,
      .groups = "drop"
    )

  dynamic_line_data <- dynamic_long %>%
    group_by(site, metric) %>%
    filter(sum(!is.na(value)) >= 2) %>%
    ungroup()

  supp_dynamic_ts <- ggplot(
    dynamic_long,
    aes(x = year, y = value, color = site)
  ) +
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
      aes(
        x = x_label,
        y = mean_val,
        label = sprintf(paste0("%.", FIG_MEAN_LABEL_DIGITS, "f"), mean_val)
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.35,
      size = FIG_ANNOT_TEXT_SIZE,
      color = "black",
      check_overlap = FIG_LABEL_CHECK_OVERLAP,
      na.rm = TRUE
    ) +
    facet_grid(
      metric ~ site,
      scales = "free_y",
      labeller = labeller(
        metric = metric_labels_panel,
        site = site_labels_panel
      ),
      drop = FALSE
    ) +
    scale_color_manual(values = site_colors) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.16))) +
    labs(x = "Water Year", y = NULL) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = FIG_AXIS_TEXT_SIZE
      ),
      strip.text.y = element_text(
        angle = 0,
        size = FIG_STRIP_TEXT_SIZE,
        hjust = 0
      ),
      strip.text.x = element_text(size = FIG_STRIP_TEXT_SIZE, hjust = 0),
      plot.margin = margin(
        FIG_LABEL_PLOT_MARGIN_PT,
        FIG_LABEL_PLOT_MARGIN_PT,
        FIG_LABEL_PLOT_MARGIN_PT,
        FIG_LABEL_PLOT_MARGIN_PT
      )
    ) +
    coord_cartesian(clip = FIG_LABEL_CLIP)

  ggsave(
    file.path(supp_dir, "ds_annual_ts.png"),
    supp_dynamic_ts,
    width = 16 * FIG_WIDTH_SCALE,
    height = 10 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
  ggsave(
    file.path(supp_pdf_dir, "ds_annual_ts.pdf"),
    supp_dynamic_ts,
    width = 16 * FIG_WIDTH_SCALE,
    height = 10 * FIG_HEIGHT_SCALE
  )

  # SUPPLEMENT: FACETED TIME SERIES - CHS

  if (!is.null(chs_data) && nrow(chs_data) > 0) {
    cat("Creating Supplement: CHS Time Series...\n")

    chs_means <- chs_data %>%
      group_by(site) %>%
      summarise(
        mean_val = mean(CHS, na.rm = TRUE),
        x_pos = max(year, na.rm = TRUE),
        x_label = max(year, na.rm = TRUE) + 0.6,
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
        aes(
          x = x_label,
          y = mean_val,
          label = sprintf(paste0("%.", FIG_MEAN_LABEL_DIGITS, "f"), mean_val)
        ),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.35,
        size = FIG_ANNOT_TEXT_SIZE,
        color = "black",
        check_overlap = FIG_LABEL_CHECK_OVERLAP,
        na.rm = TRUE
      ) +
      facet_wrap(
        ~site,
        ncol = 2,
        drop = FALSE,
        axes = "margins",
        axis.labels = "margins",
        labeller = labeller(site = site_labels_panel_chs)
      ) +
      scale_color_manual(values = site_colors) +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.16))) +
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

    ggsave(
      file.path(supp_dir, "ms_chs_annual_ts.png"),
      supp_chs_ts,
      width = 12 * FIG_WIDTH_SCALE,
      height = 10 * FIG_HEIGHT_SCALE,
      dpi = 300
    )
    ggsave(
      file.path(supp_pdf_dir, "ms_chs_annual_ts.pdf"),
      supp_chs_ts,
      width = 12 * FIG_WIDTH_SCALE,
      height = 10 * FIG_HEIGHT_SCALE
    )
    cat("  Saved CHS Time Series\n")
  }

  # SUMMARY

  for (f in list.files(main_dir, pattern = "\\.(png|pdf)$")) {
    cat("  -", f, "\n")
  }

  for (f in list.files(supp_dir, pattern = "\\.(png|pdf)$")) {
    cat("  -", f, "\n")
  }
})
