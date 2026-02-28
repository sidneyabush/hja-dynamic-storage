# supplementary figures

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

rm(list = ls())
source("config.R")

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
pred_file <- file.path(OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR, "storage_eco_response_mlr_predicted_observed.csv")
if (file.exists(pred_file)) {
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

  ggsave(
    file.path(supp_dir, "storage_ecovar_mlr_predicted_vs_observed.png"),
    p_pred,
    width = 8.8 * FIG_WIDTH_SCALE,
    height = 8.4 * FIG_HEIGHT_SCALE,
    dpi = 300
  )
  ggsave(
    file.path(supp_pdf_dir, "storage_ecovar_mlr_predicted_vs_observed.pdf"),
    p_pred,
    width = 8.8 * FIG_WIDTH_SCALE,
    height = 8.4 * FIG_HEIGHT_SCALE
  )
}

# pca precipitation anomaly supplement figure
pca_scores_file <- file.path(OUT_STATS_PCA_DIR, "pca_scores_pc1_pc2.csv")
annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
if (file.exists(pca_scores_file) && file.exists(annual_file)) {
  pca_scores <- read_csv(pca_scores_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC)

  annual <- read_csv(annual_file, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site)) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC)

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
      )

    p_pca_anom <- ggplot(pca_anno, aes(x = PC1, y = PC2)) +
      stat_ellipse(aes(color = site), linewidth = 0.7, alpha = 0.8, show.legend = FALSE) +
      geom_point(aes(color = Pws_anom), size = FIG_POINT_SIZE_MED, alpha = 0.85) +
      geom_text(
        data = centroids,
        aes(x = x, y = y, label = site, color = I(unname(SITE_COLORS[site]))),
        inherit.aes = FALSE,
        size = FIG_ANNOT_TEXT_SIZE + 1.1,
        show.legend = FALSE
      ) +
      scale_color_gradient2(
        low = "#b2182b",
        mid = "white",
        high = "#2166ac",
        midpoint = 0,
        name = expression(atop("Annual", P[ws]~"anomaly"))
      ) +
      labs(x = "Annual Dynamic Storage PC1", y = "Annual Dynamic Storage PC2") +
      theme_pub() +
      theme(
        axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
        axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
        legend.title = element_text(size = FIG_AXIS_TITLE_SIZE),
        legend.text = element_text(size = FIG_AXIS_TEXT_SIZE)
      )

    ggsave(
      file.path(supp_dir, "PCA_precip_anomaly.png"),
      p_pca_anom,
      width = 8.5 * FIG_WIDTH_SCALE,
      height = 6.2 * FIG_HEIGHT_SCALE,
      dpi = 300
    )
    ggsave(
      file.path(supp_pdf_dir, "PCA_precip_anomaly.pdf"),
      p_pca_anom,
      width = 8.5 * FIG_WIDTH_SCALE,
      height = 6.2 * FIG_HEIGHT_SCALE
    )
  }
}
