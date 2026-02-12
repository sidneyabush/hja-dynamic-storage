# -----------------------------------------------------------------------------
# Watershed Controls MLR Plots
# -----------------------------------------------------------------------------
# This script makes standalone plots for watershed-controls MLR results.
#
# Inputs:
#   - catch_chars_storage_mlr_model_avg_coef.csv
#   - catch_chars_storage_mlr_summary.csv
#
# Outputs:
#   - catch_chars_storage_mlr_beta.png
#   - catch_chars_storage_mlr_beta.pdf
#   - catch_chars_storage_mlr_model_perf.csv
#   - catch_chars_storage_mlr_coef.csv
#   - catch_chars_storage_mlr_table.csv
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)

rm(list = ls())

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
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found. Please ensure config.R exists in the repo root.")
}

output_dir <- OUT_STATS_MLR_CATCH_CHARS_DIR
plot_dir <- file.path(FIGURES_DIR, "supp", "stats", "mlr")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
table_dir <- OUT_TABLES_MLR_DIR
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

model_avg_file <- file.path(output_dir, "catch_chars_storage_mlr_model_avg_coef.csv")
summary_file <- file.path(output_dir, "catch_chars_storage_mlr_summary.csv")

if (!file.exists(model_avg_file)) stop("Missing file: catch_chars_storage_mlr_model_avg_coef.csv")
if (!file.exists(summary_file)) stop("Missing file: catch_chars_storage_mlr_summary.csv")

mlr_model_avg <- read_csv(model_avg_file, show_col_types = FALSE)
mlr_summary <- read_csv(summary_file, show_col_types = FALSE)

perf_cols <- c("Outcome", "Predictors_Final", "R2_adj", "RMSE", "AIC", "AICc", "N")
perf_cols <- perf_cols[perf_cols %in% names(mlr_summary)]
perf_df <- mlr_summary %>%
  select(all_of(perf_cols)) %>%
  arrange(Outcome)

coef_cols <- c("Outcome", "Predictor", "Beta_Std_Avg", "Inclusion_Weight", "Models_Present", "Models_DeltaAICc_LTE2")
coef_cols <- coef_cols[coef_cols %in% names(mlr_model_avg)]
coef_df <- mlr_model_avg %>%
  select(all_of(coef_cols)) %>%
  arrange(Outcome, desc(abs(Beta_Std_Avg)))

beta_plot_df <- mlr_model_avg %>%
  filter(!is.na(Beta_Std_Avg), is.finite(Beta_Std_Avg), Inclusion_Weight > 0) %>%
  mutate(
    Outcome_clean = gsub("_mean$", "", Outcome),
    beta_label = sprintf("%.2f", Beta_Std_Avg)
  )

adj_r2_lookup <- mlr_summary %>%
  mutate(
    Outcome_clean = gsub("_mean$", "", Outcome),
    Outcome_label = paste0(Outcome_clean, "\nAdj R2 = ", sprintf("%.2f", R2_adj))
  ) %>%
  select(Outcome, Outcome_label)

beta_plot_df <- beta_plot_df %>%
  left_join(adj_r2_lookup, by = "Outcome")

predictor_order <- c(
  "basin_slope", "Harvest", "Landslide_Total", "Landslide_Young",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
)

beta_plot_df <- beta_plot_df %>%
  mutate(
    Predictor = factor(Predictor, levels = rev(predictor_order)),
    Outcome_label = factor(Outcome_label, levels = unique(adj_r2_lookup$Outcome_label))
  )

p_beta <- ggplot(beta_plot_df, aes(x = Outcome_label, y = Predictor, fill = Beta_Std_Avg)) +
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
  theme_classic(base_size = FIG_BASE_SIZE) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "catch_chars_storage_mlr_beta.png"),
  p_beta,
  width = 11 * FIG_WIDTH_SCALE,
  height = 7 * FIG_HEIGHT_SCALE,
  dpi = 300
)

ggsave(
  file.path(plot_dir, "catch_chars_storage_mlr_beta.pdf"),
  p_beta,
  width = 11 * FIG_WIDTH_SCALE,
  height = 7 * FIG_HEIGHT_SCALE
)

write_csv(
  perf_df,
  file.path(table_dir, "catch_chars_storage_mlr_model_perf.csv")
)

write_csv(
  coef_df,
  file.path(table_dir, "catch_chars_storage_mlr_coef.csv")
)

results_table <- coef_df %>%
  left_join(perf_df, by = "Outcome") %>%
  arrange(Outcome, desc(abs(Beta_Std_Avg)))

write_csv(
  results_table,
  file.path(table_dir, "catch_chars_storage_mlr_table.csv")
)
