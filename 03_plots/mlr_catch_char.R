# Watershed Controls MLR Plots.
# Inputs: output_dir/watershed_char_storage_mlr_results.csv; output_dir/watershed_char_storage_mlr_summary.csv.
# Author: Sidney Bush
# Date: 2026-02-13

library(dplyr)
library(readr)
library(ggplot2)

rm(list = ls())

# Load project config
source("config.R")


output_dir <- OUT_MODELS_WATERSHED_CHAR_STORAGE_MLR_DIR
plot_dir <- file.path(FIGURES_DIR, "main")
plot_pdf_dir <- file.path(plot_dir, "pdf")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(plot_pdf_dir)) dir.create(plot_pdf_dir, recursive = TRUE)

results_file <- file.path(output_dir, "watershed_char_storage_mlr_results.csv")
summary_file <- file.path(output_dir, "watershed_char_storage_mlr_summary.csv")

if (!file.exists(results_file)) stop("Missing file: watershed_char_storage_mlr_results.csv")
if (!file.exists(summary_file)) stop("Missing file: watershed_char_storage_mlr_summary.csv")

mlr_results <- read_csv(results_file, show_col_types = FALSE)
mlr_summary <- read_csv(summary_file, show_col_types = FALSE)
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
    beta_label = sprintf("%.2f", Beta_Std)
  )

adj_r2_lookup <- mlr_summary %>%
  mutate(
    Outcome_clean = gsub("_mean$", "", Outcome),
    Outcome_label = paste0(Outcome_clean, "\nAdj R² = ", sprintf("%.2f", R2_adj))
  ) %>%
  select(Outcome, Outcome_label)

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
    Outcome_label = factor(Outcome_label, levels = unique(adj_r2_lookup$Outcome_label))
  )

p_beta <- ggplot(beta_plot_df, aes(x = Outcome_label, y = Predictor, fill = Beta_Std)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = beta_label), size = FIG_TILE_TEXT_SIZE) +
  scale_fill_gradient2(
    low = "firebrick4",
    mid = "white",
    high = "dodgerblue3",
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish,
    name = "Beta"
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    legend.position = "right",
    plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT)
  ) +
  coord_cartesian(clip = FIG_LABEL_CLIP)

ggsave(
  file.path(plot_dir, "watershed_char_storage_mlr_beta.png"),
  p_beta,
  width = 11 * FIG_WIDTH_SCALE,
  height = 7 * FIG_HEIGHT_SCALE,
  dpi = 300
)

ggsave(
  file.path(plot_pdf_dir, "watershed_char_storage_mlr_beta.pdf"),
  p_beta,
  width = 11 * FIG_WIDTH_SCALE,
  height = 7 * FIG_HEIGHT_SCALE
)
