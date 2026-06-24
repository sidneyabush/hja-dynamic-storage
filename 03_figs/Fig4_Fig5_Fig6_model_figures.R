# Figures 4, 5, and 6 model figures

# inputs:
# outputs/models/catchment_char_storage_mlr/*.csv
# outputs/models/storage_eco_response_mlr/*.csv

# outputs:
# figs_tables_pub/main/Fig4_catchment_controls.*
# figs_tables_pub/main/Fig5_ecological_response_models.*
# figs_tables_pub/main/Fig6_observed_predicted_ecological_responses.*

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(dplyr, readr, tidyr, ggplot2, ggtext, cran_repo = "https://cloud.r-project.org")

rm(list = ls())
source("config.R")

safe_ggsave <- function(filename, plot_obj, width, height, dpi = NULL, ...) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    {
      if (is.null(dpi)) {
        ggplot2::ggsave(filename, plot_obj, width = width, height = height, bg = "white", ...)
      } else {
        ggplot2::ggsave(filename, plot_obj, width = width, height = height, dpi = dpi, bg = "white", ...)
      }
      TRUE
    },
    error = function(e) {
      warning("Failed to save plot: ", filename, " (", conditionMessage(e), ")")
      FALSE
    }
  )
}

make_output_dirs <- function(...) {
  for (d in c(...)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

norm_response <- function(x) {
  out <- as.character(x)
  out[out == "Q_7Q5"] <- "Q7Q5"
  out[out == "T_7DMax"] <- "T7DMax"
  out
}

make_fig4_catchment_controls <- function() {
  output_dir <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
  plot_dir <- MS_FIG_MAIN_DIR
  plot_pdf_dir <- MS_FIG_MAIN_PDF_DIR
  plot_tiff_dir <- MS_FIG_MAIN_TIFF_DIR
  make_output_dirs(plot_dir, plot_pdf_dir, plot_tiff_dir)

  results_file <- file.path(output_dir, "catchment_char_storage_mlr_results.csv")
  summary_file <- file.path(output_dir, "catchment_char_storage_mlr_summary.csv")

  mlr_results <- read_csv(results_file, show_col_types = FALSE)
  mlr_summary <- read_csv(summary_file, show_col_types = FALSE)

  outcome_order <- STORAGE_METRIC_ORDER
  outcome_axis_labels_plotmath <- c(
    "RBI" = "plain(RBI)",
    "RCS" = "plain(RCS)",
    "FDC" = "plain(FDC)",
    "SD" = "plain(SD)",
    "WB" = "plain(WB)",
    "BF" = "plain(BF)",
    "DR" = "plain(DR)",
    "Fyw" = "plain(F)[yw]",
    "MTT" = "plain(MTT)"
  )

  format_model_p <- function(p) {
    ifelse(
      is.finite(p),
      format(signif(p, 3), scientific = TRUE, trim = TRUE),
      "NA"
    )
  }

  predictor_count <- mlr_results %>%
    group_by(Outcome) %>%
    summarise(k_predictors = sum(is.finite(Beta_Std)), .groups = "drop")

  mlr_summary <- mlr_summary %>%
    left_join(predictor_count, by = "Outcome")

  if (!("model_p_global" %in% names(mlr_summary))) {
    r2_num <- suppressWarnings(as.numeric(mlr_summary$R2))
    n_num <- suppressWarnings(as.numeric(mlr_summary$N))
    k_num <- suppressWarnings(as.numeric(mlr_summary$k_predictors))
    f_stat <- ifelse(
      is.finite(r2_num) & is.finite(n_num) & is.finite(k_num) &
        k_num > 0 & (n_num - k_num - 1) > 0 & r2_num < 1,
      (r2_num / k_num) / ((1 - r2_num) / (n_num - k_num - 1)),
      NA_real_
    )
    mlr_summary$model_p_global <- ifelse(
      is.finite(f_stat),
      pf(f_stat, k_num, n_num - k_num - 1, lower.tail = FALSE),
      NA_real_
    )
  }

  beta_plot_df <- mlr_results %>%
    filter(!is.na(Beta_Std), is.finite(Beta_Std)) %>%
    mutate(
      Outcome_clean = gsub("_mean$", "", Outcome),
      Predictor = label_catchment_predictor(Predictor),
      p_value = suppressWarnings(as.numeric(p_value)),
      beta_label = sprintf("%.2f", Beta_Std),
      sig_fontface = ifelse(is.finite(p_value) & p_value <= 0.05, "bold", "plain")
    )

  adj_r2_lookup <- mlr_summary %>%
    mutate(
      Outcome_clean = gsub("_mean$", "", Outcome),
      Outcome_label = factor(Outcome_clean, levels = outcome_order),
      Outcome_stats = paste0(
        label_storage_metric(Outcome_clean),
        ": adj R2 ",
        sprintf("%.2f", R2_adj),
        ", model p ",
        format_model_p(suppressWarnings(as.numeric(model_p_global)))
      )
    ) %>%
    arrange(Outcome_label) %>%
    mutate(Outcome_label = as.character(Outcome_label)) %>%
    select(Outcome, Outcome_label, Outcome_stats)

  beta_plot_df <- beta_plot_df %>%
    left_join(adj_r2_lookup, by = "Outcome")

  predictor_order <- c(
    "basin_slope", "Harvest", "Landslide_Total", "Landslide_Young",
    "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
  )
  predictor_order <- label_catchment_predictor(predictor_order)

  beta_plot_df <- beta_plot_df %>%
    mutate(
      Predictor = factor(Predictor, levels = rev(predictor_order)),
      Outcome_label = factor(Outcome_label, levels = outcome_order)
    )

  fig4_text_scale <- 1.0
  fig4_tile_text <- FIG_TILE_TEXT_SIZE * fig4_text_scale
  fig4_axis_text <- (FIG_AXIS_TEXT_SIZE + 1) * fig4_text_scale
  fig4_axis_title <- (FIG_AXIS_TITLE_SIZE + 1) * fig4_text_scale
  fig4_legend_title <- (FIG_AXIS_TITLE_SIZE + 1) * fig4_text_scale
  fig4_legend_text <- (FIG_AXIS_TEXT_SIZE + 1) * fig4_text_scale

  p_beta <- ggplot(beta_plot_df, aes(x = Outcome_label, y = Predictor, fill = Beta_Std)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = beta_label, fontface = sig_fontface), size = fig4_tile_text) +
    scale_x_discrete(
      labels = function(x) {
        parsed <- outcome_axis_labels_plotmath[x]
        parsed[is.na(parsed)] <- x[is.na(parsed)]
        parse(text = unname(parsed))
      }
    ) +
    scale_fill_gradient2(
      low = "firebrick3",
      mid = "white",
      high = "dodgerblue3",
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish,
      name = expression(italic(beta))
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text = element_text(size = fig4_axis_text),
      axis.title = element_text(size = fig4_axis_title),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      legend.title = element_text(size = fig4_legend_title),
      legend.text = element_text(size = fig4_legend_text),
      legend.position = "right",
      plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT)
    ) +
    coord_cartesian(clip = FIG_LABEL_CLIP)

  invisible(safe_ggsave(
    file.path(plot_dir, "Fig4_catchment_controls.png"),
    p_beta,
    width = 8.6 * FIG_WIDTH_SCALE,
    height = 5.8 * FIG_HEIGHT_SCALE,
    dpi = FIG_PREVIEW_DPI
  ))

  invisible(safe_ggsave(
    file.path(plot_pdf_dir, "Fig4_catchment_controls.pdf"),
    p_beta,
    width = 8.6 * FIG_WIDTH_SCALE,
    height = 5.8 * FIG_HEIGHT_SCALE
  ))

  invisible(safe_ggsave(
    file.path(plot_tiff_dir, "Fig4_catchment_controls.tiff"),
    p_beta,
    width = 8.6 * FIG_WIDTH_SCALE,
    height = 5.8 * FIG_HEIGHT_SCALE,
    dpi = FIG_PRODUCTION_DPI,
    compression = "lzw"
  ))
}

make_fig5_ecological_response_models <- function() {
  output_dir <- OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR
  main_dir <- MS_FIG_MAIN_DIR
  main_pdf_dir <- MS_FIG_MAIN_PDF_DIR
  main_tiff_dir <- MS_FIG_MAIN_TIFF_DIR
  make_output_dirs(main_dir, main_pdf_dir, main_tiff_dir)

  results_file <- file.path(output_dir, "storage_eco_response_mlr_results.csv")
  summary_file <- file.path(output_dir, "storage_eco_response_mlr_summary.csv")

  norm_predictor <- function(x) {
    out <- as.character(x)
    out[out %in% c("P_WetSeason", "Pwetseason", "Pwetseason")] <- "Pws"
    out
  }

  response_order <- c("Q7Q5", "T7DMax")
  predictor_order <- c("Pws", "RBI", "RCS", "FDC", "SD", "WB", "BF", "DR", "Fyw", "MTT")
  response_axis_labels <- stats::setNames(response_order, response_order)
  predictor_axis_labels_plotmath <- c(
    Pws = "plain(P)[ws]",
    RBI = "plain(RBI)",
    RCS = "plain(RCS)",
    FDC = "plain(FDC)",
    SD = "plain(SD)",
    WB = "plain(WB)",
    BF = "plain(BF)",
    DR = "plain(DR)",
    Fyw = "plain(F)[yw]",
    MTT = "plain(MTT)"
  )

  coef_raw <- read_csv(results_file, show_col_types = FALSE)

  if ("Candidate_Set" %in% names(coef_raw)) {
    selected_sets <- read_csv(summary_file, show_col_types = FALSE) %>%
      transmute(
        Response = norm_response(Response),
        Candidate_Set = suppressWarnings(as.numeric(Candidate_Set))
      ) %>%
      filter(is.finite(Candidate_Set)) %>%
      distinct(Response, Candidate_Set)

    coef_raw <- coef_raw %>%
      mutate(
        Response = norm_response(Response),
        Candidate_Set = suppressWarnings(as.numeric(Candidate_Set))
      ) %>%
      inner_join(selected_sets, by = c("Response", "Candidate_Set"))
  }

  coef_df <- coef_raw %>%
    mutate(
      Response = norm_response(Response),
      Predictor = norm_predictor(Predictor),
      Beta_Std = suppressWarnings(as.numeric(Beta_Std)),
      p_value = suppressWarnings(as.numeric(p_value))
    ) %>%
    filter(Response %in% response_order, Predictor %in% predictor_order)

  # stop if the selected ecological response models are not present
  if (nrow(coef_df) == 0) {
    stop("No rows found for Q7Q5/T7DMax in eco model results")
  }

  plot_df <- expand_grid(
    Response = response_order,
    Predictor = predictor_order
  ) %>%
    left_join(
      coef_df %>%
        group_by(Response, Predictor) %>%
        slice(1) %>%
        ungroup(),
      by = c("Response", "Predictor")
    ) %>%
    mutate(
      Response = factor(Response, levels = response_order),
      Predictor = factor(Predictor, levels = rev(predictor_order)),
      sig_fontface = ifelse(is.finite(p_value) & p_value <= 0.05, "bold", "plain"),
      label = ifelse(is.finite(Beta_Std), sprintf("%.2f", Beta_Std), "")
    )

  p <- ggplot(plot_df, aes(x = Response, y = Predictor)) +
    geom_tile(fill = "white", color = "white", linewidth = 0.3) +
    geom_tile(
      data = plot_df %>% filter(is.finite(Beta_Std)),
      aes(fill = Beta_Std),
      color = "white",
      linewidth = 0.3
    ) +
    geom_text(aes(label = label, fontface = sig_fontface), size = FIG_TILE_TEXT_SIZE * 0.62) +
    scale_fill_gradient2(
      low = "firebrick3",
      mid = "white",
      high = "dodgerblue3",
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish,
      name = expression(italic(beta))
    ) +
    scale_x_discrete(labels = response_axis_labels) +
    scale_y_discrete(
      labels = function(x) {
        parsed <- predictor_axis_labels_plotmath[x]
        parsed[is.na(parsed)] <- x[is.na(parsed)]
        parse(text = unname(parsed))
      }
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE - 5),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE - 3),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      legend.title = element_text(size = FIG_AXIS_TITLE_SIZE - 4),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE - 4)
    )

  invisible(safe_ggsave(
    file.path(main_dir, "Fig5_ecological_response_models.png"),
    p,
    width = 3.9 * FIG_WIDTH_SCALE,
    height = 4.1 * FIG_HEIGHT_SCALE,
    dpi = FIG_PREVIEW_DPI
  ))

  invisible(safe_ggsave(
    file.path(main_pdf_dir, "Fig5_ecological_response_models.pdf"),
    p,
    width = 3.9 * FIG_WIDTH_SCALE,
    height = 4.1 * FIG_HEIGHT_SCALE
  ))

  invisible(safe_ggsave(
    file.path(main_tiff_dir, "Fig5_ecological_response_models.tiff"),
    p,
    width = 3.9 * FIG_WIDTH_SCALE,
    height = 4.1 * FIG_HEIGHT_SCALE,
    dpi = FIG_PRODUCTION_DPI,
    compression = "lzw"
  ))
}

make_fig6_observed_predicted_ecological_responses <- function() {
  main_dir <- MS_FIG_MAIN_DIR
  main_pdf_dir <- MS_FIG_MAIN_PDF_DIR
  main_tiff_dir <- MS_FIG_MAIN_TIFF_DIR
  make_output_dirs(main_dir, main_pdf_dir, main_tiff_dir)

  pred_file <- file.path(
    OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
    "storage_eco_response_mlr_predicted_observed.csv"
  )
  summary_file <- file.path(
    OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
    "storage_eco_response_mlr_summary.csv"
  )

  pred_df <- read_csv(pred_file, show_col_types = FALSE) %>%
    mutate(
      Response = norm_response(Response),
      site = standardize_site_code(as.character(Site)),
      site = factor(site, levels = SITE_ORDER_HYDROMETRIC)
    ) %>%
    filter(Response %in% c("Q7Q5", "T7DMax"), is.finite(Observed), is.finite(Predicted))

  response_labels <- c(
    "Q7Q5" = "a) Q7Q5",
    "T7DMax" = "b) T7DMax"
  )
  pred_point_alpha <- 0.65
  site_shapes <- setNames(c(21, 22, 23, 24, 25, 15, 16, 17, 18, 19), SITE_ORDER_HYDROMETRIC)
  site_point_sizes <- setNames(rep(FIG_POINT_SIZE_LARGE + 1, length(SITE_ORDER_HYDROMETRIC)), SITE_ORDER_HYDROMETRIC)
  diamond_sites <- names(site_shapes)[site_shapes %in% c(18, 23)]
  site_point_sizes[diamond_sites] <- site_point_sizes[diamond_sites] + 1.1
  legend_site_levels <- SITE_ORDER_HYDROMETRIC

  format_predictor_label <- function(x) {
    labels <- c(
      Pws = "P<sub>ws</sub>",
      Fyw = "F<sub>yw</sub>"
    )
    out <- trimws(as.character(x))
    ifelse(out %in% names(labels), unname(labels[out]), out)
  }

  format_predictors <- function(x) {
    if (is.na(x) || !nzchar(trimws(as.character(x)))) {
      return("")
    }
    parts <- strsplit(gsub(";", ",", as.character(x), fixed = TRUE), ",", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    paste(format_predictor_label(parts), collapse = ", ")
  }

  model_stats <- read_csv(summary_file, show_col_types = FALSE) %>%
    transmute(
      Response = norm_response(Response),
      r2 = suppressWarnings(as.numeric(R2)),
      rmse = suppressWarnings(as.numeric(RMSE)),
      predictors = vapply(Predictors_Final, format_predictors, character(1))
    ) %>%
    filter(Response %in% c("Q7Q5", "T7DMax"))

  one_to_one_df <- pred_df %>%
    group_by(Response) %>%
    summarise(
      min_v = min(c(Observed, Predicted), na.rm = TRUE),
      max_v = max(c(Observed, Predicted), na.rm = TRUE),
      .groups = "drop"
    )

  annot_df <- model_stats %>%
    left_join(one_to_one_df, by = "Response") %>%
    mutate(
      x = min_v + 0.015 * (max_v - min_v),
      y = max_v - 0.015 * (max_v - min_v),
      label = paste0(
        "R<sup>2</sup> = ", sprintf("%.2f", r2),
        "<br>RMSE = ", signif(rmse, 3),
        "<br>Predictors: ", predictors
      )
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
    ggtext::geom_richtext(
      data = annot_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = FIG_ANNOT_TEXT_SIZE + 0.5,
      hjust = 0,
      vjust = 1,
      lineheight = 1.12,
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt")
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
    file.path(main_dir, "Fig6_observed_predicted_ecological_responses.png"),
    p_pred,
    width = 8.8 * FIG_WIDTH_SCALE,
    height = 8.4 * FIG_HEIGHT_SCALE,
    dpi = FIG_PREVIEW_DPI
  ))

  invisible(safe_ggsave(
    file.path(main_pdf_dir, "Fig6_observed_predicted_ecological_responses.pdf"),
    p_pred,
    width = 8.8 * FIG_WIDTH_SCALE,
    height = 8.4 * FIG_HEIGHT_SCALE
  ))

  invisible(safe_ggsave(
    file.path(main_tiff_dir, "Fig6_observed_predicted_ecological_responses.tiff"),
    p_pred,
    width = 8.8 * FIG_WIDTH_SCALE,
    height = 8.4 * FIG_HEIGHT_SCALE,
    dpi = FIG_PRODUCTION_DPI,
    compression = "lzw"
  ))
}

make_fig4_catchment_controls()
make_fig5_ecological_response_models()
make_fig6_observed_predicted_ecological_responses()
