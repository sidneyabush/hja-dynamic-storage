# -----------------------------------------------------------------------------
# Catchment Controls MLR Plots
# -----------------------------------------------------------------------------
# This script makes standalone plots for catchment-controls MLR results.
#
# Inputs:
#   - MLR_Storage_Catchment_Results.csv
#   - MLR_Storage_Catchment_ModelSummary.csv
#
# Outputs:
#   - Catchment_Controls_MLR_Beta_Weights.png
#   - Catchment_Controls_MLR_Beta_Weights.pdf
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

output_dir <- OUTPUT_DIR
plot_dir <- file.path(FIGURES_DIR, "Statistical")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

results_file <- file.path(output_dir, "MLR_Storage_Catchment_Results.csv")
summary_file <- file.path(output_dir, "MLR_Storage_Catchment_ModelSummary.csv")

if (!file.exists(results_file)) {
  stop("Missing file: MLR_Storage_Catchment_Results.csv")
}
if (!file.exists(summary_file)) {
  stop("Missing file: MLR_Storage_Catchment_ModelSummary.csv")
}

mlr_results <- read_csv(results_file, show_col_types = FALSE)
mlr_summary <- read_csv(summary_file, show_col_types = FALSE)

beta_plot_df <- mlr_results %>%
  filter(!is.na(Beta_Std)) %>%
  mutate(
    Outcome_clean = gsub("_mean$", "", Outcome),
    sig_label = ifelse(!is.na(p_value) & p_value < 0.05, "*", ""),
    beta_label = paste0(sprintf("%.2f", Beta_Std), sig_label)
  )

adj_r2_lookup <- mlr_summary %>%
  mutate(
    Outcome_clean = gsub("_mean$", "", Outcome),
    Outcome_label = paste0(Outcome_clean, "\nAdjR2=", sprintf("%.2f", R2_adj))
  ) %>%
  select(Outcome, Outcome_label)

beta_plot_df <- beta_plot_df %>%
  left_join(adj_r2_lookup, by = "Outcome")

predictor_order <- c(
  "Slope_mean", "Harvest", "Landslide_Total", "Landslide_Young",
  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
)

beta_plot_df <- beta_plot_df %>%
  mutate(
    Predictor = factor(Predictor, levels = rev(predictor_order)),
    Outcome_label = factor(Outcome_label, levels = unique(adj_r2_lookup$Outcome_label))
  )

p_beta <- ggplot(beta_plot_df, aes(x = Outcome_label, y = Predictor, fill = Beta_Std)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = beta_label), size = 3) +
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
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "Catchment_Controls_MLR_Beta_Weights.png"),
  p_beta,
  width = 11,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(plot_dir, "Catchment_Controls_MLR_Beta_Weights.pdf"),
  p_beta,
  width = 11,
  height = 7
)
