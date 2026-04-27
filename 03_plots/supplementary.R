# supplementary figures

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

rm(list = ls())
source("config.R")

safe_ggsave <- function(filename, plot_obj, width, height, dpi = NULL) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    {
      if (is.null(dpi)) {
        ggplot2::ggsave(filename, plot_obj, width = width, height = height, bg = "white")
      } else {
        ggplot2::ggsave(filename, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
      }
      TRUE
    },
    error = function(e) {
      warning("Failed to save plot: ", filename, " (", conditionMessage(e), ")")
      FALSE
    }
  )
}

build_corr_triangle_plot <- function(
  data_df,
  metric_map,
  legend_title = "Pearson's r",
  metric_labels = NULL,
  parse_metric_labels = FALSE
) {
  metric_map <- metric_map[metric_map %in% names(data_df)]
  if (length(metric_map) < 2) {
    return(NULL)
  }
  if (is.null(metric_labels)) {
    metric_labels <- stats::setNames(label_metric_abbrev(names(metric_map)), names(metric_map))
  }

  corr_input <- as.data.frame(lapply(data_df[, unname(metric_map), drop = FALSE], function(x) {
    suppressWarnings(as.numeric(x))
  }))
  colnames(corr_input) <- names(metric_map)

  corr_mat <- suppressWarnings(cor(corr_input, use = "pairwise.complete.obs"))
  if (!is.matrix(corr_mat) || ncol(corr_mat) < 2) {
    return(NULL)
  }

  idx <- which(lower.tri(corr_mat), arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(NULL)
  }

  x_levels <- names(metric_map)[-length(metric_map)]
  # Standard orientation: lower-triangle correlations shown as a top-left wedge.
  # This keeps the broadest row at the top (manuscript style).
  y_levels <- names(metric_map)[-1]

  corr_tri <- tibble(
    row_metric = rownames(corr_mat)[idx[, 1]],
    col_metric = colnames(corr_mat)[idx[, 2]],
    r = corr_mat[idx]
  ) %>%
    mutate(
      row_metric = factor(row_metric, levels = y_levels),
      col_metric = factor(col_metric, levels = x_levels),
      label = ifelse(
        is.finite(r),
        ifelse(abs(r) < 0.005, "|r| < 0.005", sprintf("%.2f", r)),
        ""
      )
    )

  p <- ggplot(corr_tri, aes(x = col_metric, y = row_metric, fill = r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label), size = FIG_TILE_TEXT_SIZE * 1.25) +
    scale_fill_gradient2(
      low = "firebrick3",
      mid = "white",
      high = "dodgerblue3",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey85",
      name = legend_title
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = FIG_AXIS_TEXT_SIZE + 1),
      axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1)
    )

  if (isTRUE(parse_metric_labels)) {
    p <- p +
      scale_x_discrete(
        labels = function(x) {
          parsed <- metric_labels[x]
          parsed[is.na(parsed)] <- x[is.na(parsed)]
          parse(text = unname(parsed))
        }
      ) +
      scale_y_discrete(
        labels = function(x) {
          parsed <- metric_labels[x]
          parsed[is.na(parsed)] <- x[is.na(parsed)]
          parse(text = unname(parsed))
        }
      )
  } else {
    p <- p +
      scale_x_discrete(labels = metric_labels) +
      scale_y_discrete(labels = metric_labels)
  }

  p
}

build_diagnostic_p_plot <- function(
  diag_df,
  response_col,
  response_order,
  response_label_fn = identity
) {
  if (nrow(diag_df) == 0) {
    return(NULL)
  }

  diag_df <- diag_df %>%
    mutate(
      response = as.character(.data[[response_col]]),
      n_residuals = suppressWarnings(as.numeric(n_residuals)),
      shapiro_p = suppressWarnings(as.numeric(shapiro_p)),
      ncv_p = suppressWarnings(as.numeric(ncv_p))
    )

  long_df <- bind_rows(
    diag_df %>%
      transmute(
        response,
        n_residuals,
        test = "Shapiro-Wilk normality test",
        p_value = shapiro_p
      ),
    diag_df %>%
      transmute(
        response,
        n_residuals,
        test = "Non-constant variance test",
        p_value = ncv_p
      )
  ) %>%
    filter(is.finite(p_value)) %>%
    mutate(
      test = factor(
        test,
        levels = c("Shapiro-Wilk normality test", "Non-constant variance test"),
        labels = c("a) Shapiro-Wilk normality test", "b) Non-constant variance test")
      ),
      response = factor(response, levels = response_order),
      pass_p05 = p_value >= 0.05
    )

  if (nrow(long_df) == 0) {
    return(NULL)
  }

  ggplot(long_df, aes(x = response, y = p_value)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "firebrick3", linewidth = 0.5) +
    geom_point(
      aes(fill = pass_p05),
      shape = 21,
      color = "grey20",
      alpha = 0.92,
      stroke = 0.25,
      size = 3.6
    ) +
    facet_wrap(~ test, ncol = 1) +
    scale_x_discrete(labels = function(x) response_label_fn(as.character(x))) +
    scale_y_continuous(
      limits = c(-0.04, 1.04),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::number_format(accuracy = 0.1),
      oob = scales::squish,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(
      values = c(`TRUE` = "dodgerblue3", `FALSE` = "firebrick3"),
      breaks = c(TRUE, FALSE),
      labels = c(expression(p >= 0.05), expression(p < 0.05)),
      name = "p-value"
    ) +
    labs(x = NULL, y = "Diagnostic p-value") +
    coord_cartesian(clip = "off") +
    theme_pub() +
    theme(
      axis.text.x = element_text(size = FIG_AXIS_TEXT_SIZE + 1, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
      strip.text = element_text(size = FIG_STRIP_TEXT_SIZE + 1),
      strip.background = element_rect(fill = "grey97", color = "grey55", linewidth = 0.4),
      panel.border = element_rect(fill = NA, color = "grey55", linewidth = 0.45),
      panel.spacing = grid::unit(0.7, "lines"),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
    )
}

supp_dir <- MS_FIG_SUPP_DIR
supp_pdf_dir <- MS_FIG_SUPP_PDF_DIR
main_dir <- MS_FIG_MAIN_DIR
main_pdf_dir <- MS_FIG_MAIN_PDF_DIR
supp_explore_dir <- SUPP_EXPLORATORY_DIR
supp_explore_pdf_dir <- SUPP_EXPLORATORY_PDF_DIR
for (d in c(supp_dir, supp_pdf_dir, main_dir, main_pdf_dir, supp_explore_dir, supp_explore_pdf_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# figS1 met context is maintained in its own script but dispatched from here
figs1_script <- file.path(REPO_DIR, "03_plots", "FigS1_met_context.R")
if (file.exists(figs1_script)) {
  status <- system2("Rscript", shQuote(figs1_script), stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    warning("FigS1_met_context.R failed when called from supplementary.R")
  }
}

# pooled eco predicted vs observed
pred_file <- file.path(
  OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
  "storage_eco_response_mlr_predicted_observed.csv"
)

if (!is.na(pred_file) && nzchar(pred_file) && file.exists(pred_file)) {
  pred_df <- read_csv(pred_file, show_col_types = FALSE) %>%
    mutate(
      Response = as.character(Response),
      Response = ifelse(Response == "Q_7Q5", "Q7Q5", Response),
      Response = ifelse(Response == "T_7DMax", "T7DMax", Response),
      site = standardize_site_code(as.character(Site)),
      site = factor(site, levels = SITE_ORDER_HYDROMETRIC)
    ) %>%
    filter(Response %in% c("Q7Q5", "T7DMax"), is.finite(Observed), is.finite(Predicted))

  response_labels <- c(
    "Q7Q5" = "a) Q7Q5",
    "T7DMax" = "b) T7DMax"
  )
  response_labels <- stats::setNames(
    wrap_plot_label(response_labels, width = 36),
    names(response_labels)
  )
  pred_point_alpha <- 0.65
  # Keep site symbols consistent with Figure 3.
  site_shapes <- setNames(c(21, 22, 23, 24, 25, 15, 16, 17, 18, 19), SITE_ORDER_HYDROMETRIC)
  site_point_sizes <- setNames(rep(FIG_POINT_SIZE_LARGE + 1, length(SITE_ORDER_HYDROMETRIC)), SITE_ORDER_HYDROMETRIC)
  diamond_sites <- names(site_shapes)[site_shapes %in% c(18, 23)]
  site_point_sizes[diamond_sites] <- site_point_sizes[diamond_sites] + 1.1
  legend_site_levels <- SITE_ORDER_HYDROMETRIC

  r2_df <- pred_df %>%
    group_by(Response) %>%
    summarise(
      r2 = suppressWarnings(cor(Observed, Predicted, use = "complete.obs")^2),
      .groups = "drop"
    )

  one_to_one_df <- pred_df %>%
    group_by(Response) %>%
    summarise(
      min_v = min(c(Observed, Predicted), na.rm = TRUE),
      max_v = max(c(Observed, Predicted), na.rm = TRUE),
      .groups = "drop"
    )

  r2_annot_df <- r2_df %>%
    left_join(one_to_one_df, by = "Response") %>%
    mutate(
      x = min_v + 0.015 * (max_v - min_v),
      y = max_v - 0.015 * (max_v - min_v),
      label_parse = paste0("R^2==", sprintf("%.2f", r2))
    )

  p_pred <- ggplot(pred_df, aes(x = Observed, y = Predicted)) +
    geom_point(
      aes(color = site, fill = site, shape = site, size = site),
      alpha = pred_point_alpha,
      stroke = 0.9
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_text(
      data = one_to_one_df,
      aes(x = max_v - 0.08 * (max_v - min_v), y = max_v - 0.08 * (max_v - min_v), label = "1:1"),
      inherit.aes = FALSE,
      color = "grey35",
      size = FIG_ANNOT_TEXT_SIZE + 1.2
    ) +
    geom_text(
      data = r2_annot_df,
      aes(x = x, y = y, label = label_parse),
      inherit.aes = FALSE,
      parse = TRUE,
      size = FIG_ANNOT_TEXT_SIZE + 1.2,
      hjust = 0,
      vjust = 1
    ) +
    facet_wrap(~ Response, scales = "free", ncol = 1, labeller = as_labeller(response_labels)) +
    scale_color_manual(values = SITE_COLORS, drop = FALSE) +
    scale_fill_manual(values = SITE_COLORS, drop = FALSE, guide = "none") +
    scale_shape_manual(values = site_shapes, drop = FALSE) +
    scale_size_manual(values = site_point_sizes, drop = FALSE) +
    guides(
      color = guide_legend(
        title = NULL,
        override.aes = list(
          shape = unname(site_shapes[legend_site_levels]),
          fill = unname(SITE_COLORS[legend_site_levels]),
          size = unname(site_point_sizes[legend_site_levels]),
          alpha = pred_point_alpha,
          stroke = 0.9
        )
      ),
      shape = "none",
      size = "none"
    ) +
    labs(x = "Observed", y = "Predicted", color = NULL) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
      strip.text = element_text(size = FIG_STRIP_TEXT_SIZE + 1),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
    )

  invisible(safe_ggsave(
    file.path(main_dir, "Fig9_eco_observed_v_predicted.png"),
    p_pred,
    width = 8.8 * FIG_WIDTH_SCALE,
    height = 8.4 * FIG_HEIGHT_SCALE,
    dpi = 300
  ))
  invisible(safe_ggsave(
    file.path(main_pdf_dir, "Fig9_eco_observed_v_predicted.pdf"),
    p_pred,
    width = 8.8 * FIG_WIDTH_SCALE,
    height = 8.4 * FIG_HEIGHT_SCALE
  ))
}

if (!isTRUE(WRITE_AUX_OUTPUTS)) {
  unlink(file.path(supp_dir, "FigSX_pca_Pws_anomaly.png"))
  unlink(file.path(supp_pdf_dir, "FigSX_pca_Pws_anomaly.pdf"))
}

if (isTRUE(WRITE_AUX_OUTPUTS)) {
  # pca precipitation anomaly supplement figure (auxiliary)
  pca_scores_file <- file.path(OUT_STATS_PCA_DIR, "pca_scores_pc1_pc2.csv")
  variance_file <- file.path(OUT_STATS_PCA_DIR, "pca_variance_explained.csv")
  annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
  if (file.exists(pca_scores_file) && file.exists(annual_file)) {
  pca_scores <- read_csv(pca_scores_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC)

  annual <- read_csv(annual_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC)

  pc1_pct <- NA_real_
  pc2_pct <- NA_real_
  if (file.exists(variance_file)) {
    pca_var <- read_csv(variance_file, show_col_types = FALSE) %>%
      transmute(
        PC = as.character(PC),
        Variance_Explained = suppressWarnings(as.numeric(Variance_Explained))
      )
    pc1_raw <- pca_var %>% filter(PC == "PC1") %>% pull(Variance_Explained)
    pc2_raw <- pca_var %>% filter(PC == "PC2") %>% pull(Variance_Explained)
    if (length(pc1_raw) > 0 && is.finite(pc1_raw[[1]])) {
      pc1_pct <- ifelse(pc1_raw[[1]] <= 1, 100 * pc1_raw[[1]], pc1_raw[[1]])
    }
    if (length(pc2_raw) > 0 && is.finite(pc2_raw[[1]])) {
      pc2_pct <- ifelse(pc2_raw[[1]] <= 1, 100 * pc2_raw[[1]], pc2_raw[[1]])
    }
  }
  x_lab <- if (is.finite(pc1_pct)) {
    paste0("Annual Dynamic Storage PC1 (", sprintf("%.1f", pc1_pct), "%)")
  } else {
    "Annual Dynamic Storage PC1"
  }
  y_lab <- if (is.finite(pc2_pct)) {
    paste0("Annual Dynamic Storage PC2 (", sprintf("%.1f", pc2_pct), "%)")
  } else {
    "Annual Dynamic Storage PC2"
  }

  if (!("Pws" %in% names(annual))) {
    annual$Pws <- NA_real_
  }
  if ("P_WetSeason" %in% names(annual)) {
    annual$Pws <- dplyr::coalesce(annual$Pws, annual$P_WetSeason)
  }
  if ("precip_nov_may_mm" %in% names(annual)) {
    annual$Pws <- dplyr::coalesce(annual$Pws, annual$precip_nov_may_mm)
  }

  pca_anno <- pca_scores %>%
    left_join(
      annual %>%
        group_by(site) %>%
        mutate(Pws_anom = Pws - mean(Pws, na.rm = TRUE)) %>%
        ungroup() %>%
        select(site, year, Pws_anom),
      by = c("site", "year")
    ) %>%
    filter(is.finite(PC1), is.finite(PC2), is.finite(Pws_anom))

  if (nrow(pca_anno) > 0) {
    centroids <- pca_anno %>%
      group_by(site) %>%
      summarise(
        x = median(PC1, na.rm = TRUE),
        y = median(PC2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(site_col = unname(SITE_COLORS[as.character(site)]))

    p_pca_anom <- ggplot(pca_anno, aes(x = PC1, y = PC2)) +
      stat_ellipse(
        aes(
          group = site,
          fill = I(unname(SITE_COLORS[as.character(site)])),
          color = I(unname(SITE_COLORS[as.character(site)]))
        ),
        geom = "polygon",
        linewidth = 0.7,
        alpha = 0.14,
        show.legend = FALSE
      ) +
      geom_point(
        aes(fill = Pws_anom),
        shape = 21,
        color = "grey30",
        stroke = 0.25,
        size = FIG_POINT_SIZE_MED,
        alpha = 0.9
      ) +
      geom_text(
        data = centroids,
        aes(x = x, y = y, label = site, color = I(site_col)),
        inherit.aes = FALSE,
        size = FIG_ANNOT_TEXT_SIZE + 1.1,
        show.legend = FALSE
      ) +
      scale_fill_gradient2(
        low = "#b2182b",
        mid = "white",
        high = "#2166ac",
        midpoint = 0,
        name = expression(atop("Annual", P[ws]~"anomaly"))
      ) +
      labs(x = x_lab, y = y_lab) +
      theme_pub() +
      theme(
        axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
        axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
        legend.title = element_text(size = FIG_AXIS_TITLE_SIZE),
        legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
      )

    invisible(safe_ggsave(
      file.path(supp_dir, "FigSX_pca_Pws_anomaly.png"),
      p_pca_anom,
      width = 8.5 * FIG_WIDTH_SCALE,
      height = 6.2 * FIG_HEIGHT_SCALE,
      dpi = 300
    ))
    invisible(safe_ggsave(
      file.path(supp_pdf_dir, "FigSX_pca_Pws_anomaly.pdf"),
      p_pca_anom,
      width = 8.5 * FIG_WIDTH_SCALE,
      height = 6.2 * FIG_HEIGHT_SCALE
    ))
    unlink(file.path(supp_dir, "PCA_precip_anomaly.png"))
    unlink(file.path(supp_pdf_dir, "PCA_precip_anomaly.pdf"))
  }
  }
}

# supplementary within-group storage correlation heatmaps
site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
annual_corr_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
if (file.exists(site_file) && file.exists(annual_corr_file)) {
  site_df <- read_csv(site_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC)

  annual_corr <- read_csv(annual_corr_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC)

  dynamic_map <- c(
    RBI = "RBI",
    RCS = "RCS",
    FDC = "FDC",
    SD = "SD",
    WB = "WB"
  )
  mobile_map <- c(
    CHS = "CHS_mean",
    DR = "DR",
    Fyw = "Fyw",
    MTT = "MTT"
  )

  unlink(file.path(supp_dir, "FigS2_dynamic_storage_corr.png"))
  unlink(file.path(supp_pdf_dir, "FigS2_dynamic_storage_corr.pdf"))
  unlink(file.path(supp_dir, "FigSX_dynamic_metrics_corr.png"))
  unlink(file.path(supp_pdf_dir, "FigSX_dynamic_metrics_corr.pdf"))

  unlink(file.path(supp_dir, "FigS3_mobile_storage_corr.png"))
  unlink(file.path(supp_pdf_dir, "FigS3_mobile_storage_corr.pdf"))
  unlink(file.path(supp_dir, "FigSX_mobile_metrics_corr.png"))
  unlink(file.path(supp_pdf_dir, "FigSX_mobile_metrics_corr.pdf"))
}

# supplementary diagnostics summary for reported residual tests (table only)
catch_diag_file <- file.path(
  OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR,
  "catchment_char_storage_mlr_diagnostics.csv"
)
catch_diag_tbl <- tibble()
if (file.exists(catch_diag_file)) {
  catch_diag <- read_csv(catch_diag_file, show_col_types = FALSE)
  catch_diag_tbl <- catch_diag %>%
    transmute(
      model_set = "Catchment characteristics",
      response_variable = label_metric_abbrev(gsub("_mean$", "", as.character(Outcome))),
      n_residuals = suppressWarnings(as.numeric(n_residuals)),
      shapiro_p = suppressWarnings(as.numeric(shapiro_p)),
      ncv_p = suppressWarnings(as.numeric(ncv_p)),
      shapiro_p = ifelse(
        is.finite(shapiro_p),
        paste0(format(signif(shapiro_p, 3), scientific = TRUE, trim = TRUE), ifelse(shapiro_p < 0.05, "*", "")),
        NA_character_
      ),
      ncv_p = ifelse(
        is.finite(ncv_p),
        paste0(format(signif(ncv_p, 3), scientific = TRUE, trim = TRUE), ifelse(ncv_p < 0.05, "*", "")),
        NA_character_
      )
    )
}

eco_diag_file <- file.path(
  OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
  "storage_eco_response_mlr_diagnostics.csv"
)
eco_diag_tbl <- tibble()
if (file.exists(eco_diag_file)) {
  eco_diag <- read_csv(eco_diag_file, show_col_types = FALSE)
  eco_diag_tbl <- eco_diag %>%
    transmute(
      model_set = "Ecological responses",
      response_variable = as.character(Response),
      n_residuals = suppressWarnings(as.numeric(n_residuals)),
      shapiro_p = suppressWarnings(as.numeric(shapiro_p)),
      ncv_p = suppressWarnings(as.numeric(ncv_p)),
      shapiro_p = ifelse(
        is.finite(shapiro_p),
        paste0(format(signif(shapiro_p, 3), scientific = TRUE, trim = TRUE), ifelse(shapiro_p < 0.05, "*", "")),
        NA_character_
      ),
      ncv_p = ifelse(
        is.finite(ncv_p),
        paste0(format(signif(ncv_p, 3), scientific = TRUE, trim = TRUE), ifelse(ncv_p < 0.05, "*", "")),
        NA_character_
      )
    )
}

# supplementary diagnostics table used for SI reporting.
diag_table <- bind_rows(catch_diag_tbl, eco_diag_tbl)
if (nrow(diag_table) > 0) {
  diag_table <- diag_table %>%
    arrange(model_set, response_variable)
  write_csv(diag_table, file.path(supp_dir, "TableS7_mlr_model_diagnostics.csv"))
  unlink(file.path(supp_dir, "Table_SX_mlr_model_diagnostics.csv"))
}

# remove legacy diagnostics figures from SI deliverables.
unlink(file.path(supp_dir, "FigS4_catchment_mlr_diagnostics.png"))
unlink(file.path(supp_pdf_dir, "FigS4_catchment_mlr_diagnostics.pdf"))
unlink(file.path(supp_dir, "FigS5_eco_mlr_diagnostics.png"))
unlink(file.path(supp_pdf_dir, "FigS5_eco_mlr_diagnostics.pdf"))

# supplementary dynamic-mobile scatter matrix supporting Figure 6 interpretation
if (file.exists(site_file)) {
  site_df_scatter <- read_csv(site_file, show_col_types = FALSE) %>%
    mutate(
      site = standardize_site_code(site),
      across(c(RBI_mean, RCS_mean, FDC_mean, SD_mean, WB_mean, CHS_mean, DR, Fyw, MTT), ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC)

  dynamic_map_scatter <- c(
    RBI = "RBI_mean",
    RCS = "RCS_mean",
    FDC = "FDC_mean",
    SD = "SD_mean",
    WB = "WB_mean"
  )
  mobile_map_scatter <- c(
    BF = "CHS_mean",
    DR = "DR",
    Fyw = "Fyw",
    MTT = "MTT"
  )

  scatter_df <- bind_rows(lapply(names(dynamic_map_scatter), function(dn) {
    bind_rows(lapply(names(mobile_map_scatter), function(mn) {
      tibble(
        site = site_df_scatter$site,
        Dynamic = dn,
        Mobile = mn,
        x = site_df_scatter[[dynamic_map_scatter[[dn]]]],
        y = site_df_scatter[[mobile_map_scatter[[mn]]]]
      )
    }))
  })) %>%
    mutate(
      Dynamic = factor(Dynamic, levels = names(dynamic_map_scatter)),
      Mobile = factor(Mobile, levels = names(mobile_map_scatter))
    )

  # keep site-shape mapping aligned with Figure 3 and other manuscript plots.
  site_shapes <- setNames(c(21, 22, 23, 24, 25, 15, 16, 17, 18, 19), SITE_ORDER_HYDROMETRIC)

  panel_stats <- scatter_df %>%
    group_by(Mobile, Dynamic) %>%
    summarise(
      n = sum(is.finite(x) & is.finite(y)),
      r = suppressWarnings(cor(x, y, use = "complete.obs")),
      x_min = suppressWarnings(min(x, na.rm = TRUE)),
      x_max = suppressWarnings(max(x, na.rm = TRUE)),
      y_min = suppressWarnings(min(y, na.rm = TRUE)),
      y_max = suppressWarnings(max(y, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      label = ifelse(
        is.finite(r),
        paste0("r=", sprintf("%.2f", r), "\n", "n=", n),
        paste0("r=NA", "\n", "n=", n)
      ),
      x_span = ifelse(is.finite(x_max - x_min) & (x_max - x_min) > 0, x_max - x_min, 1),
      y_span = ifelse(is.finite(y_max - y_min) & (y_max - y_min) > 0, y_max - y_min, 1),
      # Place annotation in padded top-right space within panel bounds.
      x = x_max + 0.12 * x_span,
      y = y_max + 0.12 * y_span
    )

  smooth_df <- scatter_df %>%
    group_by(Mobile, Dynamic) %>%
    mutate(n_pair = sum(is.finite(x) & is.finite(y))) %>%
    ungroup() %>%
    filter(n_pair >= 3)

  p_scatter <- ggplot(scatter_df, aes(x = x, y = y)) +
    geom_point(
      aes(color = site, fill = site, shape = site),
      size = FIG_POINT_SIZE_MED + 0.1,
      stroke = 0.7,
      alpha = 0.88,
      na.rm = TRUE
    ) +
    geom_smooth(
      data = smooth_df,
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.45,
      color = "grey25",
      na.rm = TRUE
    ) +
    geom_text(
      data = panel_stats,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = FIG_ANNOT_TEXT_SIZE,
      color = "grey35"
    ) +
    facet_grid(
      Mobile ~ Dynamic,
      scales = "free",
      labeller = labeller(
        Dynamic = as_labeller(
          c(RBI = "RBI", RCS = "RCS", FDC = "FDC", SD = "SD", WB = "WB"),
          label_parsed
        ),
        Mobile = as_labeller(
          c(BF = "BF", DR = "DR", Fyw = "F[yw]", MTT = "MTT"),
          label_parsed
        )
      )
    ) +
    scale_color_manual(values = SITE_COLORS, drop = FALSE) +
    scale_fill_manual(values = SITE_COLORS, drop = FALSE, guide = "none") +
    scale_shape_manual(values = site_shapes, drop = FALSE) +
    guides(
      color = guide_legend(
        title = NULL,
        override.aes = list(
          shape = unname(site_shapes[SITE_ORDER_HYDROMETRIC]),
          fill = unname(SITE_COLORS[SITE_ORDER_HYDROMETRIC]),
          size = rep(FIG_POINT_SIZE_MED + 0.6, length(SITE_ORDER_HYDROMETRIC)),
          alpha = 0.88,
          stroke = 0.7
        )
      ),
      shape = "none"
    ) +
    labs(
      x = "Dynamic storage metric value",
      y = "Mobile storage metric value",
      color = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.20))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.20))) +
    theme_pub() +
    theme(
      axis.text.x = element_text(size = FIG_AXIS_TEXT_SIZE, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "grey65", linewidth = 0.3),
      strip.text = element_text(size = FIG_STRIP_TEXT_SIZE),
      strip.background = element_rect(fill = "grey97", color = "grey55", linewidth = 0.4),
      panel.border = element_rect(fill = NA, color = "grey55", linewidth = 0.45),
      panel.spacing = grid::unit(0.65, "lines"),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      plot.margin = margin(8, 8, 8, 8)
    )

  invisible(safe_ggsave(
    file.path(supp_dir, "FigS4_dynamic_mobile_scatter_matrix.png"),
    p_scatter,
    width = 12.8 * FIG_WIDTH_SCALE,
    height = 8.8 * FIG_HEIGHT_SCALE,
    dpi = 300
  ))
  invisible(safe_ggsave(
    file.path(supp_pdf_dir, "FigS4_dynamic_mobile_scatter_matrix.pdf"),
    p_scatter,
    width = 12.8 * FIG_WIDTH_SCALE,
    height = 8.8 * FIG_HEIGHT_SCALE
  ))
  unlink(file.path(supp_dir, "FigS6_dynamic_mobile_scatter_matrix.png"))
  unlink(file.path(supp_pdf_dir, "FigS6_dynamic_mobile_scatter_matrix.pdf"))
}

# supplementary ec-vs-ca chs relationship plots (auxiliary)
ec_ca_pairs_file <- file.path(OUT_MET_MOBILE_DIR, "annual_gw_prop_ec_ca_site_year_pairs.csv")
if (isTRUE(WRITE_AUX_OUTPUTS) && file.exists(ec_ca_pairs_file)) {
  ec_ca_pairs <- read_csv(ec_ca_pairs_file, show_col_types = FALSE) %>%
    mutate(
      SITECODE = standardize_site_code(SITECODE),
      waterYear = suppressWarnings(as.integer(waterYear)),
      CHS_EC = suppressWarnings(as.numeric(CHS_EC)),
      CHS_CA = suppressWarnings(as.numeric(CHS_CA))
    ) %>%
    filter(
      SITECODE %in% SITE_ORDER_CHEMISTRY,
      SITECODE != "Look",
      is.finite(CHS_EC),
      is.finite(CHS_CA)
    ) %>%
    mutate(
      SITECODE = factor(SITECODE, levels = SITE_ORDER_CHEMISTRY)
    )

  if (nrow(ec_ca_pairs) > 0) {
    safe_r2 <- function(x, y) {
      keep <- is.finite(x) & is.finite(y)
      if (sum(keep) < 2) {
        return(NA_real_)
      }
      suppressWarnings(cor(x[keep], y[keep], use = "complete.obs")^2)
    }

    safe_p <- function(x, y) {
      keep <- is.finite(x) & is.finite(y)
      if (sum(keep) < 3) {
        return(NA_real_)
      }
      suppressWarnings(cor.test(x[keep], y[keep], method = "pearson")$p.value)
    }

    format_p <- function(p) {
      if (!is.finite(p)) {
        return("NA")
      }
      if (p < 0.001) {
        return("<0.001")
      }
      sprintf("%.3f", p)
    }

    format_p_label <- function(p) {
      p_txt <- format_p(p)
      if (identical(p_txt, "<0.001")) {
        return("<0.001")
      }
      paste0("p=", p_txt)
    }

    site_rel_stats <- ec_ca_pairs %>%
      group_by(SITECODE) %>%
      summarise(
        n = n(),
        r2 = safe_r2(CHS_EC, CHS_CA),
        p_value = safe_p(CHS_EC, CHS_CA),
        .groups = "drop"
      )

    facet_labels <- site_rel_stats %>%
      mutate(
        label = ifelse(
          is.finite(r2),
          paste0(
            as.character(SITECODE), "\n",
            "n=", n, ", R²=", sprintf("%.2f", r2), ", ", vapply(p_value, format_p_label, FUN.VALUE = character(1))
          ),
          paste0(as.character(SITECODE), "\n", "n=", n, ", R²=NA, p=NA")
        )
      ) %>%
      {stats::setNames(.$label, as.character(.$SITECODE))}

    p_ec_ca_site <- ggplot(ec_ca_pairs, aes(x = CHS_EC, y = CHS_CA)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
      geom_point(
        aes(color = SITECODE),
        size = FIG_POINT_SIZE_MED,
        alpha = 0.85,
        show.legend = FALSE
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        color = "black",
        linewidth = 0.55
      ) +
      facet_wrap(~ SITECODE, ncol = 3, labeller = as_labeller(facet_labels)) +
      scale_color_manual(values = SITE_COLORS, drop = TRUE) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
      coord_equal() +
      labs(
        x = "Baseflow Fraction from Electrical Conductivity (BF from EC)",
        y = "Baseflow Fraction from Calcium (BF from Ca)"
      ) +
      theme_pub() +
      theme(
        axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
        axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
        strip.text = element_text(size = FIG_STRIP_TEXT_SIZE)
      )

    overall_r2 <- safe_r2(ec_ca_pairs$CHS_EC, ec_ca_pairs$CHS_CA)
    overall_p <- safe_p(ec_ca_pairs$CHS_EC, ec_ca_pairs$CHS_CA)
    overall_n <- nrow(ec_ca_pairs)
    anno_label <- if (is.finite(overall_r2)) {
      paste0("n=", overall_n, ", R²=", sprintf("%.2f", overall_r2), ", ", format_p_label(overall_p))
    } else {
      paste0("n=", overall_n, ", R²=NA, p=NA")
    }

    p_ec_ca_overall <- ggplot(ec_ca_pairs, aes(x = CHS_EC, y = CHS_CA)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
      geom_point(
        aes(color = SITECODE),
        size = FIG_POINT_SIZE_MED + 0.3,
        alpha = 0.85
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        color = "black",
        linewidth = 0.65
      ) +
      annotate(
        "text",
        x = 0.02,
        y = 0.98,
        label = anno_label,
        hjust = 0,
        vjust = 1,
        size = FIG_ANNOT_TEXT_SIZE + 1.2
      ) +
      scale_color_manual(values = SITE_COLORS, drop = TRUE) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
      coord_equal() +
      labs(
        x = "Baseflow Fraction from Electrical Conductivity (BF from EC)",
        y = "Baseflow Fraction from Calcium (BF from Ca)",
        color = NULL
      ) +
      theme_pub() +
      theme(
        axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
        axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
        legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
      )

    invisible(safe_ggsave(
      file.path(supp_explore_dir, "FigSX_chs_ec_vs_ca_by_site.png"),
      p_ec_ca_site,
      width = 9.6 * FIG_WIDTH_SCALE,
      height = 8.2 * FIG_HEIGHT_SCALE,
      dpi = 300
    ))
    invisible(safe_ggsave(
      file.path(supp_explore_pdf_dir, "FigSX_chs_ec_vs_ca_by_site.pdf"),
      p_ec_ca_site,
      width = 9.6 * FIG_WIDTH_SCALE,
      height = 8.2 * FIG_HEIGHT_SCALE
    ))
    invisible(safe_ggsave(
      file.path(supp_explore_dir, "FigSX_chs_ec_vs_ca_overall.png"),
      p_ec_ca_overall,
      width = 7.2 * FIG_WIDTH_SCALE,
      height = 6.4 * FIG_HEIGHT_SCALE,
      dpi = 300
    ))
    invisible(safe_ggsave(
      file.path(supp_explore_pdf_dir, "FigSX_chs_ec_vs_ca_overall.pdf"),
      p_ec_ca_overall,
      width = 7.2 * FIG_WIDTH_SCALE,
      height = 6.4 * FIG_HEIGHT_SCALE
    ))

    # keep SI output directory clean by removing older copies from ms_materials/supp.
    unlink(file.path(supp_dir, "FigSX_chs_ec_vs_ca_by_site.png"))
    unlink(file.path(supp_pdf_dir, "FigSX_chs_ec_vs_ca_by_site.pdf"))
    unlink(file.path(supp_dir, "FigSX_chs_ec_vs_ca_overall.png"))
    unlink(file.path(supp_pdf_dir, "FigSX_chs_ec_vs_ca_overall.pdf"))
  }
}
