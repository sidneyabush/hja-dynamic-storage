# -----------------------------------------------------------------------------
# Plots: Storage Predictors for Thermal and Low-Flow Responses
# -----------------------------------------------------------------------------
# Inputs:
#   - storage_ecovar_mlr_results.csv
#   - storage_ecovar_mlr_summary.csv
#
# Outputs:
#   - storage_ecovar_mlr_beta.png
#   - storage_ecovar_mlr_beta.pdf
#   - storage_ecovar_mlr_model_perf.csv
#   - storage_ecovar_mlr_coef.csv
#   - storage_ecovar_mlr_table.csv
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

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

output_dir <- OUT_MODELS_STORAGE_ECOVAR_MLR_DIR
plot_dir <- file.path(FIGURES_DIR, "main")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
table_dir <- OUT_TABLES_MLR_DIR
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

# Remove legacy duplicate table names from prior workflow versions.
legacy_tables <- c(
  "storage_eco_mlr_model_perf.csv",
  "storage_eco_mlr_coef.csv",
  "storage_eco_mlr_table.csv",
  "storage_eco_mlr_r2_heatmap_table.csv"
)
legacy_paths <- file.path(table_dir, legacy_tables)
legacy_paths <- legacy_paths[file.exists(legacy_paths)]
if (length(legacy_paths) > 0) unlink(legacy_paths)

models_file <- file.path(output_dir, "storage_ecovar_mlr_results.csv")
summary_file <- file.path(output_dir, "storage_ecovar_mlr_summary.csv")

if (!file.exists(models_file)) {
  stop("Missing file: storage_ecovar_mlr_results.csv")
}
if (!file.exists(summary_file)) {
  stop("Missing file: storage_ecovar_mlr_summary.csv")
}

model_coefs <- read_csv(models_file, show_col_types = FALSE)
model_summary <- read_csv(summary_file, show_col_types = FALSE)
coverage_file <- file.path(output_dir, "storage_ecovar_mlr_coverage.csv")
coverage_df <- if (file.exists(coverage_file)) {
  read_csv(coverage_file, show_col_types = FALSE) %>%
    transmute(
      Site = as.character(Site),
      Response = as.character(Response),
      model_status,
      reason_not_fit,
      n_response,
      usable_predictors
    )
} else {
  tibble(
    Site = character(),
    Response = character(),
    model_status = character(),
    reason_not_fit = character(),
    n_response = numeric(),
    usable_predictors = numeric()
  )
}

response_order <- c("T_7DMax", "Q_7Q5", "T_Q7Q5")
site_order <- SITE_ORDER_HYDROMETRIC
predictor_order <- c("RBI", "RCS", "FDC", "SD", "CHS", "WB")

# Predictor availability by site (distinguish true no-data from not-selected).
annual_path <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_path)) annual_path <- file.path(OUT_MASTER_DIR, LEGACY_ANNUAL_FILE)

annual_master <- if (file.exists(annual_path)) read_csv(annual_path, show_col_types = FALSE) else tibble()

annual_predictors <- predictor_order
annual_avail <- if (nrow(annual_master) > 0) {
  annual_master %>%
    filter(site %in% site_order) %>%
    group_by(site) %>%
    summarise(
      across(any_of(annual_predictors), ~ any(is.finite(.x)), .names = "{.col}"),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = -site, names_to = "Predictor", values_to = "has_predictor_data")
} else {
  tibble(site = character(), Predictor = character(), has_predictor_data = logical())
}
predictor_avail <- annual_avail %>%
  group_by(site, Predictor) %>%
  summarise(has_predictor_data = any(has_predictor_data), .groups = "drop")

perf_df <- model_summary %>%
  mutate(
    Response = factor(Response, levels = response_order),
    Site = factor(Site, levels = site_order)
  ) %>%
  arrange(Site, Response) %>%
  mutate(
    Predictors_Final = ifelse(is.na(Predictors_Final) | Predictors_Final == "", "-", Predictors_Final)
  ) %>%
  select(any_of(c(
    "Site", "Response", "Predictors_Final", "R2", "R2_adj", "RMSE", "AICc", "n",
    "RMSE_LOOCV", "RMSE_LOOCV_MEAN_RUNS", "R2_LOOCV",
    "delta_RMSE_LOOCV_minus_model", "delta_RMSE_LOOCV_mean_runs_minus_model"
  )))

coef_df <- model_coefs %>%
  filter(Predictor %in% predictor_order) %>%
  mutate(
    Response = factor(Response, levels = response_order),
    Site = factor(Site, levels = site_order),
    Predictor = factor(Predictor, levels = predictor_order),
    Beta_Std_Display = ifelse(is.na(Beta_Std), "-", sprintf("%.2f", Beta_Std))
  ) %>%
  select(any_of(c("Site", "Response", "Predictor", "Beta_Std", "Beta_Std_Display", "p_value", "VIF"))) %>%
  arrange(Response, Predictor, Site)

response_labels <- tibble(
  Response = response_order,
  Response_label = response_order
)
response_labels_panel <- make_panel_label_map(response_labels$Response_label)

r2_lookup <- model_summary %>%
  transmute(
    Site = as.character(Site),
    Response = as.character(Response),
    R2_model = R2
  )

beta_values <- model_coefs %>%
  filter(Predictor %in% predictor_order) %>%
  transmute(
    Site = as.character(Site),
    Response = as.character(Response),
    Predictor = as.character(Predictor),
    Beta_Std = Beta_Std
  )

full_grid <- expand_grid(
  Response = response_order,
  Predictor = predictor_order,
  Site = site_order
)

beta_plot_df <- full_grid %>%
  mutate(
    Response = as.character(Response),
    Predictor = as.character(Predictor),
    Site = as.character(Site)
  ) %>%
  left_join(r2_lookup, by = c("Site", "Response")) %>%
  left_join(beta_values, by = c("Site", "Response", "Predictor")) %>%
  left_join(predictor_avail %>% rename(Site = site), by = c("Site", "Predictor")) %>%
  left_join(coverage_df, by = c("Site", "Response")) %>%
  left_join(response_labels %>% mutate(Response = as.character(Response)), by = "Response") %>%
  mutate(
    has_predictor_data = ifelse(is.na(has_predictor_data), FALSE, has_predictor_data),
    beta_fill = case_when(
      !has_predictor_data ~ NA_real_,
      model_status == "not_fit" ~ NA_real_,
      TRUE ~ Beta_Std
    ),
    tile_label = case_when(
      !has_predictor_data ~ "-",
      model_status == "not_fit" ~ "-",
      is.finite(Beta_Std) ~ sprintf("%.2f", Beta_Std),
      TRUE ~ ""
    ),
    Site = factor(Site, levels = site_order),
    Predictor = factor(Predictor, levels = rev(predictor_order)),
    Response_label = factor(Response_label, levels = response_labels$Response_label)
  )

p_beta <- ggplot(beta_plot_df, aes(x = Site, y = Predictor, fill = beta_fill)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = tile_label), size = FIG_TILE_TEXT_SIZE, lineheight = 0.9) +
  facet_wrap(
    ~Response_label,
    ncol = 1,
    scales = "fixed",
    labeller = labeller(Response_label = response_labels_panel),
    axes = "margins",
    axis.labels = "margins"
  ) +
  scale_fill_gradient2(
    low = "firebrick4",
    mid = "white",
    high = "dodgerblue3",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    oob = scales::squish,
    na.value = "white",
    name = "Beta"
  ) +
  labs(x = "Site", y = "Storage predictor") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    strip.text = element_text(size = FIG_STRIP_TEXT_SIZE, hjust = 0),
    plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT),
    legend.position = "right"
  ) +
  coord_cartesian(clip = FIG_LABEL_CLIP)

save_plot_safe <- function(path, plot_obj, width, height, dpi = NULL) {
  tryCatch(
    {
      if (is.null(dpi)) {
        ggsave(path, plot_obj, width = width, height = height)
      } else {
        ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
      }
      TRUE
    },
    error = function(e) {
      warning(sprintf("Failed to save plot to %s: %s", path, e$message))
      FALSE
    }
  )
}

save_plot_safe(
  file.path(plot_dir, "storage_ecovar_mlr_beta.png"),
  p_beta,
  width = 11 * FIG_WIDTH_SCALE,
  height = 13 * FIG_HEIGHT_SCALE,
  dpi = 300
)

save_plot_safe(
  file.path(plot_dir, "storage_ecovar_mlr_beta.pdf"),
  p_beta,
  width = 11 * FIG_WIDTH_SCALE,
  height = 13 * FIG_HEIGHT_SCALE
)

# Companion performance heatmap: site x response adjusted R2
r2_heat_df <- expand_grid(
  Site = site_order,
  Response = response_order
) %>%
  mutate(
    Site = as.character(Site),
    Response = as.character(Response)
  ) %>%
  left_join(
    model_summary %>%
      transmute(Site = as.character(Site), Response = as.character(Response), adj_r2 = R2_adj),
    by = c("Site", "Response")
  ) %>%
  left_join(coverage_df, by = c("Site", "Response")) %>%
  mutate(
    Site = factor(Site, levels = site_order),
    Response = factor(Response, levels = response_order),
    r2_label = case_when(
      is.finite(adj_r2) ~ sprintf("%.2f", adj_r2),
      TRUE ~ "-"
    )
  )

p_r2 <- ggplot(r2_heat_df, aes(x = Site, y = Response, fill = adj_r2)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = r2_label), size = FIG_TILE_TEXT_SIZE) +
  scale_fill_gradient2(
    low = "#f7fbff",
    mid = "#6baed6",
    high = "#08306b",
    midpoint = 0.4,
    limits = c(0, 1),
    oob = scales::squish,
    na.value = "white",
    name = "Adj R2"
  ) +
  labs(x = "Site", y = "Response") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    strip.text = element_text(size = FIG_STRIP_TEXT_SIZE, hjust = 0)
  )

save_plot_safe(
  file.path(plot_dir, "storage_ecovar_mlr_r2_heatmap.png"),
  p_r2,
  width = 10 * FIG_WIDTH_SCALE,
  height = 4 * FIG_HEIGHT_SCALE,
  dpi = 300
)
save_plot_safe(
  file.path(plot_dir, "storage_ecovar_mlr_r2_heatmap.pdf"),
  p_r2,
  width = 10 * FIG_WIDTH_SCALE,
  height = 4 * FIG_HEIGHT_SCALE
)

save_csv_safe <- function(df, path) {
  tryCatch(
    {
      write_csv(df, path)
      TRUE
    },
    error = function(e) {
      warning(sprintf("Failed to write table to %s: %s", path, e$message))
      FALSE
    }
  )
}

pretty_text <- function(x) gsub("_", " ", x)

perf_export <- perf_df %>%
  mutate(
    Response = pretty_text(as.character(Response)),
    Predictors_Final = pretty_text(as.character(Predictors_Final))
  ) %>%
  rename(
    site = Site,
    response = Response,
    predictors_final = Predictors_Final,
    r2 = R2,
    r2_adj = R2_adj,
    rmse = RMSE,
    aicc = AICc,
    rmse_loocv = RMSE_LOOCV,
    rmse_loocv_mean_runs = RMSE_LOOCV_MEAN_RUNS,
    r2_loocv = R2_LOOCV,
    delta_rmse_loocv_minus_model = delta_RMSE_LOOCV_minus_model,
    delta_rmse_loocv_mean_runs_minus_model = delta_RMSE_LOOCV_mean_runs_minus_model
  ) %>%
  mutate(
    r2 = signif(r2, 3),
    r2_adj = signif(r2_adj, 3),
    rmse = signif(rmse, 3),
    aicc = signif(aicc, 3),
    rmse_loocv = signif(rmse_loocv, 3),
    rmse_loocv_mean_runs = signif(rmse_loocv_mean_runs, 3),
    r2_loocv = signif(r2_loocv, 3),
    delta_rmse_loocv_minus_model = signif(delta_rmse_loocv_minus_model, 3),
    delta_rmse_loocv_mean_runs_minus_model = signif(delta_rmse_loocv_mean_runs_minus_model, 3)
  )

save_csv_safe(
  perf_export,
  file.path(table_dir, "storage_ecovar_mlr_model_perf.csv")
)
save_csv_safe(
  r2_heat_df %>%
    transmute(
      site = as.character(Site),
      response = pretty_text(as.character(Response)),
      adj_r2 = signif(adj_r2, 3),
      model_status,
      reason_not_fit,
      n_response
    ),
  file.path(table_dir, "storage_ecovar_mlr_r2_heatmap_table.csv")
)

save_csv_safe(
  coef_df,
  file.path(table_dir, "storage_ecovar_mlr_coef.csv")
)

manuscript_table <- beta_plot_df %>%
  transmute(
    site = as.character(Site),
    response = pretty_text(as.character(Response)),
    predictor = as.character(Predictor),
    beta_std = signif(Beta_Std, 3),
    r2_model = signif(R2_model, 3),
    tile_label = tile_label,
    model_status,
    reason_not_fit,
    n_response,
    usable_predictors
  ) %>%
  left_join(
    perf_export %>%
      mutate(
        site = as.character(site),
        response = as.character(response)
      ),
    by = c("site", "response")
  ) %>%
  arrange(factor(site, levels = site_order), response, predictor)

save_csv_safe(
  manuscript_table,
  file.path(table_dir, "storage_ecovar_mlr_table.csv")
)
