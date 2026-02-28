# catchment controls mlr plots.
# inputs: output_dir/catchment_char_storage_mlr_results.csv; output_dir/catchment_char_storage_mlr_summary.csv.
# author: sidney bush
# date: 2026-02-13

library(dplyr)
library(readr)
library(ggplot2)

rm(list = ls())

# load project config
source("config.R")

safe_ggsave <- function(filename, plot_obj, width, height, dpi = NULL) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ext <- tools::file_ext(filename)
  tmp_file <- tempfile(
    pattern = "fig7_",
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


output_dir <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
ALPHA <- 0.05
plot_dir <- MS_FIG_MAIN_DIR
plot_pdf_dir <- MS_FIG_MAIN_PDF_DIR
supp_plot_dir <- MS_FIG_SUPP_DIR
supp_plot_pdf_dir <- MS_FIG_SUPP_PDF_DIR
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_plot_pdf_dir, recursive = TRUE, showWarnings = FALSE)

results_file <- file.path(output_dir, "catchment_char_storage_mlr_results.csv")
summary_file <- file.path(output_dir, "catchment_char_storage_mlr_summary.csv")

if (!file.exists(results_file)) stop("Missing file: catchment_char_storage_mlr_results.csv")
if (!file.exists(summary_file)) stop("Missing file: catchment_char_storage_mlr_summary.csv")

mlr_results <- read_csv(results_file, show_col_types = FALSE)
mlr_summary <- read_csv(summary_file, show_col_types = FALSE)
outcome_order <- STORAGE_METRIC_ORDER

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
if ("Predictors_Final" %in% names(mlr_summary)) {
  mlr_summary <- mlr_summary %>%
    mutate(Predictors_Final = label_catchment_predictor_list(Predictors_Final))
}

perf_cols <- c("Outcome", "Predictors_Final", "R2_adj", "RMSE", "AICc", "N")
perf_cols <- perf_cols[perf_cols %in% names(mlr_summary)]
perf_df <- mlr_summary %>%
  select(all_of(perf_cols)) %>%
  arrange(Outcome)

coef_cols <- c("Outcome", "Predictor", "Beta_Std", "p_value", "VIF", "R2", "R2_adj", "RMSE", "AICc")
coef_cols <- coef_cols[coef_cols %in% names(mlr_results)]
coef_df <- mlr_results %>%
  mutate(Predictor = label_catchment_predictor(Predictor)) %>%
  select(all_of(coef_cols)) %>%
  arrange(Outcome, desc(abs(Beta_Std)))

beta_plot_df <- mlr_results %>%
  filter(!is.na(Beta_Std), is.finite(Beta_Std)) %>%
  mutate(
    Outcome_clean = gsub("_mean$", "", Outcome),
    Predictor = label_catchment_predictor(Predictor),
    p_value = suppressWarnings(as.numeric(p_value)),
    sig_label = ifelse(is.finite(p_value) & p_value <= ALPHA, "*", ""),
    beta_label = paste0(sprintf("%.2f", Beta_Std), sig_label)
  )

adj_r2_lookup <- mlr_summary %>%
  mutate(
    Outcome_clean = gsub("_mean$", "", Outcome),
    Outcome_label = factor(Outcome_clean, levels = outcome_order),
    Outcome_stats = paste0(
      Outcome_clean,
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

p_beta <- ggplot(beta_plot_df, aes(x = Outcome_label, y = Predictor, fill = Beta_Std)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = beta_label), size = FIG_TILE_TEXT_SIZE) +
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
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
    legend.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
    legend.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
    legend.position = "right",
    plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT)
  ) +
  coord_cartesian(clip = FIG_LABEL_CLIP)

invisible(safe_ggsave(
  file.path(plot_dir, "Fig7_catch_char_mlr_beta.png"),
  p_beta,
  width = 10.2 * FIG_WIDTH_SCALE,
  height = 6.8 * FIG_HEIGHT_SCALE,
  dpi = 300
))

invisible(safe_ggsave(
  file.path(plot_pdf_dir, "Fig7_catch_char_mlr_beta.pdf"),
  p_beta,
  width = 10.2 * FIG_WIDTH_SCALE,
  height = 6.8 * FIG_HEIGHT_SCALE
))
