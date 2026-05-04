# figure 8 pooled ecological response observed vs predicted

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

rm(list = ls())
source("config.R")

safe_ggsave <- function(filename, plot_obj, width, height, dpi = NULL) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ext <- tools::file_ext(filename)
  tmp_file <- tempfile(
    pattern = "fig8_",
    tmpdir = tempdir(),
    fileext = ifelse(nzchar(ext), paste0(".", ext), "")
  )
  tryCatch(
    {
      if (is.null(dpi)) {
        ggplot2::ggsave(tmp_file, plot_obj, width = width, height = height, bg = "white")
      } else {
        ggplot2::ggsave(tmp_file, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
      }
      ok <- file.copy(tmp_file, filename, overwrite = TRUE)
      unlink(tmp_file)
      if (!isTRUE(ok)) {
        stop("Failed to copy rendered file to destination")
      }
      TRUE
    },
    error = function(e) {
      unlink(tmp_file)
      warning("Failed to save plot: ", filename, " (", conditionMessage(e), ")")
      FALSE
    }
  )
}

main_dir <- MS_FIG_MAIN_DIR
main_pdf_dir <- MS_FIG_MAIN_PDF_DIR
for (d in c(main_dir, main_pdf_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

pred_file <- file.path(
  OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
  "storage_eco_response_mlr_predicted_observed.csv"
)
if (!file.exists(pred_file)) {
  stop("Missing required file: ", pred_file)
}

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
pred_point_alpha <- 0.65
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
  file.path(main_dir, "Fig8_eco_observed_v_predicted.png"),
  p_pred,
  width = 8.8 * FIG_WIDTH_SCALE,
  height = 8.4 * FIG_HEIGHT_SCALE,
  dpi = 300
))
invisible(safe_ggsave(
  file.path(main_pdf_dir, "Fig8_eco_observed_v_predicted.pdf"),
  p_pred,
  width = 8.8 * FIG_WIDTH_SCALE,
  height = 8.4 * FIG_HEIGHT_SCALE
))
