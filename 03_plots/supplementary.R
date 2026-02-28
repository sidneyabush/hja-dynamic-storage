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

build_corr_triangle_plot <- function(data_df, metric_map, legend_title = "Corr") {
  metric_map <- metric_map[metric_map %in% names(data_df)]
  if (length(metric_map) < 2) {
    return(NULL)
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
      label = ifelse(is.finite(r), sprintf("%.2f", r), "")
    )

  ggplot(corr_tri, aes(x = col_metric, y = row_metric, fill = r)) +
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
}

supp_dir <- MS_FIG_SUPP_DIR
supp_pdf_dir <- MS_FIG_SUPP_PDF_DIR
for (d in c(supp_dir, supp_pdf_dir)) {
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
# Prefer legacy ecovar predicted-vs-observed output because it matches the
# manuscript-facing pooled diagnostics used for the supplementary figure.
pred_file_candidates <- c(
  file.path(OUTPUT_DIR, "models", "storage_ecovar_mlr", "storage_ecovar_mlr_predicted_observed.csv"),
  file.path(OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR, "storage_eco_response_mlr_predicted_observed.csv")
)
pred_file <- pred_file_candidates[file.exists(pred_file_candidates)][1]

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

  response_labels <- c("Q7Q5" = "a) Q7Q5", "T7DMax" = "b) T7DMax")

  r2_df <- pred_df %>%
    group_by(Response) %>%
    summarise(
      r2 = suppressWarnings(cor(Observed, Predicted, use = "complete.obs")^2),
      x = min(Observed, na.rm = TRUE) + 0.05 * diff(range(Observed, na.rm = TRUE)),
      y = max(Predicted, na.rm = TRUE) - 0.08 * diff(range(Predicted, na.rm = TRUE)),
      .groups = "drop"
    )

  one_to_one_df <- pred_df %>%
    group_by(Response) %>%
    summarise(
      min_v = min(c(Observed, Predicted), na.rm = TRUE),
      max_v = max(c(Observed, Predicted), na.rm = TRUE),
      .groups = "drop"
    )

  p_pred <- ggplot(pred_df, aes(x = Observed, y = Predicted)) +
    geom_point(aes(color = site), size = FIG_POINT_SIZE_MED, alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_text(
      data = one_to_one_df,
      aes(x = max_v - 0.08 * (max_v - min_v), y = max_v - 0.08 * (max_v - min_v), label = "1:1"),
      inherit.aes = FALSE,
      color = "grey35",
      size = FIG_ANNOT_TEXT_SIZE + 1.2
    ) +
    geom_text(
      data = r2_df,
      aes(x = x, y = y, label = paste0("R² = ", sprintf("%.2f", r2))),
      inherit.aes = FALSE,
      size = FIG_ANNOT_TEXT_SIZE + 1.2,
      hjust = 0,
      vjust = 1
    ) +
    facet_wrap(~ Response, scales = "free", ncol = 1, labeller = as_labeller(response_labels)) +
    scale_color_manual(values = SITE_COLORS, drop = FALSE) +
    labs(x = "Observed", y = "Predicted", color = "Site") +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
      strip.text = element_text(size = FIG_STRIP_TEXT_SIZE + 1),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
    )

  invisible(safe_ggsave(
    file.path(supp_dir, "FigSX_eco_observed_v_predicted.png"),
    p_pred,
    width = 8.8 * FIG_WIDTH_SCALE,
    height = 8.4 * FIG_HEIGHT_SCALE,
    dpi = 300
  ))
  invisible(safe_ggsave(
    file.path(supp_pdf_dir, "FigSX_eco_observed_v_predicted.pdf"),
    p_pred,
    width = 8.8 * FIG_WIDTH_SCALE,
    height = 8.4 * FIG_HEIGHT_SCALE
  ))

  unlink(file.path(supp_dir, "storage_ecovar_mlr_predicted_vs_observed.png"))
  unlink(file.path(supp_pdf_dir, "storage_ecovar_mlr_predicted_vs_observed.pdf"))
}

# pca precipitation anomaly supplement figure
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
    MTT1 = "MTT1",
    MTT2 = "MTT2"
  )

  p_dynamic_corr <- build_corr_triangle_plot(annual_corr, dynamic_map, legend_title = "Corr")
  if (!is.null(p_dynamic_corr)) {
    invisible(safe_ggsave(
      file.path(supp_dir, "FigSx_dynamic_metrics_corr.png"),
      p_dynamic_corr,
      width = 7.8 * FIG_WIDTH_SCALE,
      height = 6.6 * FIG_HEIGHT_SCALE,
      dpi = 300
    ))
    invisible(safe_ggsave(
      file.path(supp_pdf_dir, "FigSx_dynamic_metrics_corr.pdf"),
      p_dynamic_corr,
      width = 7.8 * FIG_WIDTH_SCALE,
      height = 6.6 * FIG_HEIGHT_SCALE
    ))
  }

  p_mobile_corr <- build_corr_triangle_plot(site_df, mobile_map, legend_title = "Corr")
  if (!is.null(p_mobile_corr)) {
    invisible(safe_ggsave(
      file.path(supp_dir, "FigSx_mobile_metrics_corr.png"),
      p_mobile_corr,
      width = 7.8 * FIG_WIDTH_SCALE,
      height = 6.6 * FIG_HEIGHT_SCALE,
      dpi = 300
    ))
    invisible(safe_ggsave(
      file.path(supp_pdf_dir, "FigSx_mobile_metrics_corr.pdf"),
      p_mobile_corr,
      width = 7.8 * FIG_WIDTH_SCALE,
      height = 6.6 * FIG_HEIGHT_SCALE
    ))
  }
}
