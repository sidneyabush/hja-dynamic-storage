# -----------------------------------------------------------------------------
# Plots: Storage Predictors for Thermal and Low-Flow Responses
# -----------------------------------------------------------------------------
# Inputs:
#   - storage_eco_mlr_results.csv
#   - storage_eco_mlr_summary.csv
#
# Outputs:
#   - storage_eco_mlr_beta.png
#   - storage_eco_mlr_beta.pdf
#   - storage_eco_mlr_model_perf.csv
#   - storage_eco_mlr_coef.csv
#   - storage_eco_mlr_table.csv
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

output_dir <- OUT_STATS_MLR_ECO_DIR
plot_dir <- file.path(FIGURES_DIR, "supp", "stats", "mlr")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
table_dir <- OUT_TABLES_MLR_DIR
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

models_file <- file.path(output_dir, "storage_eco_mlr_results.csv")
summary_file <- file.path(output_dir, "storage_eco_mlr_summary.csv")

if (!file.exists(models_file)) {
  stop("Missing file: storage_eco_mlr_results.csv")
}
if (!file.exists(summary_file)) {
  stop("Missing file: storage_eco_mlr_summary.csv")
}

model_coefs <- read_csv(models_file, show_col_types = FALSE)
model_summary <- read_csv(summary_file, show_col_types = FALSE)

response_order <- c("T_7DMax", "Q_7Q5", "T_Q7Q5")

beta_plot_df <- model_coefs %>%
  filter(!is.na(Beta_Std)) %>%
  mutate(
    Response = factor(Response, levels = response_order),
    Predictor = factor(Predictor),
    beta_label = sprintf("%.2f", Beta_Std)
  )

perf_df <- model_summary %>%
  mutate(Response = factor(Response, levels = response_order)) %>%
  arrange(Response) %>%
  select(any_of(c("Response", "Predictors_Final", "R2_adj", "RMSE", "AIC", "AICc", "n")))

coef_df <- model_coefs %>%
  select(any_of(c("Response", "Predictor", "Beta_Std", "p_value", "VIF"))) %>%
  arrange(Response, desc(abs(Beta_Std)))

write_csv(
  perf_df,
  file.path(table_dir, "storage_eco_mlr_model_perf.csv")
)

write_csv(
  coef_df,
  file.path(table_dir, "storage_eco_mlr_coef.csv")
)

manuscript_table <- coef_df %>%
  left_join(perf_df, by = "Response") %>%
  arrange(Response, desc(abs(Beta_Std)))

write_csv(
  manuscript_table,
  file.path(table_dir, "storage_eco_mlr_table.csv")
)

# Preserve predictor ordering from strongest average absolute beta to weakest.
pred_order <- beta_plot_df %>%
  group_by(Predictor) %>%
  summarise(mean_abs_beta = mean(abs(Beta_Std), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_beta)) %>%
  pull(Predictor)

beta_plot_df <- beta_plot_df %>%
  mutate(Predictor = factor(Predictor, levels = rev(pred_order)))

p_beta <- ggplot(beta_plot_df, aes(x = Response, y = Predictor, fill = Beta_Std)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = beta_label), size = 3) +
  scale_fill_gradient2(
    low = "firebrick4",
    mid = "white",
    high = "dodgerblue3",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    oob = scales::squish,
    name = "Beta"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "storage_eco_mlr_beta.png"),
  p_beta,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(plot_dir, "storage_eco_mlr_beta.pdf"),
  p_beta,
  width = 9,
  height = 5
)
