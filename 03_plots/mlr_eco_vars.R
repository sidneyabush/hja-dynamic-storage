# plots for eco response predictors and model diagnostics
# inputs output_dir storage_ecovar_mlr_coverage csv
# inputs output_dir storage_ecovar_mlr_predicted_observed csv
# author sidney bush
# date 2026 02 13

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

rm(list = ls())

# load project config
source("config.R")


output_dir <- OUT_MODELS_STORAGE_ECOVAR_MLR_DIR
ALPHA <- 0.05
plot_dir <- file.path(FIGURES_DIR, "main")
plot_pdf_dir <- file.path(plot_dir, "pdf")
supp_plot_dir <- file.path(FIGURES_DIR, "supp")
supp_plot_pdf_dir <- file.path(supp_plot_dir, "pdf")

ensure_dirs <- function(paths) {
  for (p in paths) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
}

ensure_dirs(c(plot_dir, plot_pdf_dir, supp_plot_dir, supp_plot_pdf_dir))

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

save_plot_pair <- function(name, plot_obj, width, height, to_supp = FALSE) {
  png_dir <- if (to_supp) supp_plot_dir else plot_dir
  pdf_dir <- if (to_supp) supp_plot_pdf_dir else plot_pdf_dir
  png_ok <- save_plot_safe(
    file.path(png_dir, paste0(name, ".png")),
    plot_obj,
    width,
    height,
    dpi = 300
  )
  pdf_ok <- save_plot_safe(
    file.path(pdf_dir, paste0(name, ".pdf")),
    plot_obj,
    width,
    height
  )
  invisible(c(png = png_ok, pdf = pdf_ok))
}

models_file <- file.path(output_dir, "storage_ecovar_mlr_results.csv")
summary_file <- file.path(output_dir, "storage_ecovar_mlr_summary.csv")
pred_obs_file <- file.path(
  output_dir,
  "storage_ecovar_mlr_predicted_observed.csv"
)

if (!file.exists(models_file)) {
  stop("Missing file: storage_ecovar_mlr_results.csv")
}
if (!file.exists(summary_file)) {
  stop("Missing file: storage_ecovar_mlr_summary.csv")
}

to_internal_response <- function(x) {
  out <- as.character(x)
  out[out == "Q7Q5"] <- "Q_7Q5"
  out[out == "TQ7Q5"] <- "T_Q7Q5"
  out[out == "T7DMax"] <- "T_7DMax"
  out
}

to_internal_predictor <- function(x) {
  out <- as.character(x)
  out[out == "Pws"] <- "P_WetSeason"
  out
}

to_internal_predictor_list <- function(x) {
  vals <- as.character(x)
  vapply(
    vals,
    function(val) {
      if (is.na(val) || !nzchar(val) || identical(val, "-")) {
        return(val)
      }
      parts <- trimws(strsplit(val, ";", fixed = TRUE)[[1]])
      parts <- to_internal_predictor(parts)
      paste(parts, collapse = "; ")
    },
    character(1)
  )
}

model_coefs <- read_csv(models_file, show_col_types = FALSE) %>%
  mutate(
    Response = to_internal_response(Response),
    Predictor = to_internal_predictor(Predictor)
  )
model_summary <- read_csv(summary_file, show_col_types = FALSE) %>%
  mutate(
    Response = to_internal_response(Response),
    Predictors_Final = to_internal_predictor_list(Predictors_Final)
  )
coverage_file <- file.path(output_dir, "storage_ecovar_mlr_coverage.csv")
coverage_df <- if (file.exists(coverage_file)) {
  read_csv(coverage_file, show_col_types = FALSE) %>%
    transmute(
      Site = as.character(Site),
      Response = to_internal_response(as.character(Response)),
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

response_order <- c("Q_7Q5", "T_Q7Q5", "T_7DMax")
beta_response_order <- c("Q_7Q5", "T_7DMax")
pred_diag_response_order <- c("Q_7Q5", "T_7DMax")
site_order <- SITE_ORDER_HYDROMETRIC
predictor_order <- c(
  "P_WetSeason",
  STORAGE_METRIC_ORDER[
    STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD", "WB", "CHS")
  ]
)

response_display <- setNames(gsub("_", "", response_order), response_order)
predictor_display <- setNames(gsub("_", "", predictor_order), predictor_order)
predictor_axis_labels <- function(x) {
  keys <- as.character(x)
  mapped <- unname(predictor_display[keys])
  mapped[keys == "P_WetSeason"] <- "P[ws]"
  parse(text = mapped)
}

# predictor availability by site so we separate no data from not selected
master_dir <- file.path(OUTPUT_DIR, "master")
annual_path <- file.path(master_dir, MASTER_ANNUAL_FILE)

annual_master <- if (file.exists(annual_path)) {
  read_csv(annual_path, show_col_types = FALSE)
} else {
  tibble()
}
if (!("P_WetSeason" %in% names(annual_master))) {
  annual_master$P_WetSeason <- NA_real_
}
if ("precip_nov_may_mm" %in% names(annual_master)) {
  annual_master$P_WetSeason <- dplyr::coalesce(
    annual_master$P_WetSeason,
    annual_master$precip_nov_may_mm
  )
}
if ("P_NovJan" %in% names(annual_master)) {
  annual_master$P_WetSeason <- dplyr::coalesce(
    annual_master$P_WetSeason,
    annual_master$P_NovJan
  )
}
if ("precip_nov_jan_mm" %in% names(annual_master)) {
  annual_master$P_WetSeason <- dplyr::coalesce(
    annual_master$P_WetSeason,
    annual_master$precip_nov_jan_mm
  )
}

annual_predictors <- predictor_order
annual_avail <- if (nrow(annual_master) > 0) {
  annual_master %>%
    filter(site %in% site_order) %>%
    group_by(site) %>%
    summarise(
      across(
        any_of(annual_predictors),
        ~ any(is.finite(.x)),
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = -site,
      names_to = "Predictor",
      values_to = "has_predictor_data"
    )
} else {
  tibble(
    site = character(),
    Predictor = character(),
    has_predictor_data = logical()
  )
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
    Predictors_Final = ifelse(
      is.na(Predictors_Final) | Predictors_Final == "",
      "-",
      Predictors_Final
    )
  ) %>%
  select(any_of(c(
    "Site",
    "Response",
    "Predictors_Final",
    "R2",
    "R2_adj",
    "RMSE",
    "AICc",
    "n",
    "RMSE_LOOCV",
    "RMSE_LOOCV_MEAN_RUNS",
    "R2_LOOCV",
    "delta_RMSE_LOOCV_minus_model",
    "delta_RMSE_LOOCV_mean_runs_minus_model"
  )))

coef_df <- model_coefs %>%
  filter(Predictor %in% predictor_order) %>%
  mutate(
    Response = factor(Response, levels = response_order),
    Site = factor(Site, levels = site_order),
    Predictor = factor(Predictor, levels = predictor_order),
    Beta_Std_Display = ifelse(is.na(Beta_Std), "-", sprintf("%.2f", Beta_Std))
  ) %>%
  select(any_of(c(
    "Site",
    "Response",
    "Predictor",
    "Beta_Std",
    "Beta_Std_Display",
    "p_value",
    "VIF"
  ))) %>%
  arrange(Response, Predictor, Site)

response_labels_all <- tibble(
  Response = response_order,
  Response_label = unname(response_display[response_order])
)
response_labels_panel_all <- make_panel_label_map(
  response_labels_all$Response_label
)
beta_response_labels <- response_labels_all %>%
  filter(Response %in% beta_response_order)
beta_response_labels_panel <- make_panel_label_map(
  beta_response_labels$Response_label
)
pred_diag_response_labels <- response_labels_all %>%
  filter(Response %in% pred_diag_response_order)
pred_diag_response_labels_panel <- make_panel_label_map(
  pred_diag_response_labels$Response_label
)

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
    Beta_Std = suppressWarnings(as.numeric(Beta_Std)),
    p_value = suppressWarnings(as.numeric(p_value))
  )

full_grid <- expand_grid(
  Response = beta_response_order,
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
  left_join(
    predictor_avail %>% rename(Site = site),
    by = c("Site", "Predictor")
  ) %>%
  left_join(coverage_df, by = c("Site", "Response")) %>%
  left_join(
    beta_response_labels %>% mutate(Response = as.character(Response)),
    by = "Response"
  ) %>%
  mutate(
    has_predictor_data = ifelse(
      is.na(has_predictor_data),
      FALSE,
      has_predictor_data
    ),
    beta_fill = case_when(
      !has_predictor_data ~ NA_real_,
      model_status == "not_fit" ~ NA_real_,
      TRUE ~ Beta_Std
    ),
    tile_label = case_when(
      !has_predictor_data ~ "-",
      model_status == "not_fit" ~ "-",
      is.finite(Beta_Std) ~ paste0(
        sprintf("%.2f", Beta_Std),
        ifelse(is.finite(p_value) & p_value <= ALPHA, "*", "")
      ),
      TRUE ~ ""
    ),
    Site = factor(Site, levels = site_order),
    Predictor = factor(Predictor, levels = rev(predictor_order)),
    Response_label = factor(
      Response_label,
      levels = beta_response_labels$Response_label
    )
  )

all_sites_results_file <- file.path(
  output_dir,
  "storage_ecovar_mlr_all_sites_results.csv"
)
if (!file.exists(all_sites_results_file)) {
  stop("Missing file: storage_ecovar_mlr_all_sites_results.csv")
}

all_sites_coefs <- read_csv(all_sites_results_file, show_col_types = FALSE) %>%
  mutate(
    Response = to_internal_response(Response),
    Predictor = to_internal_predictor(Predictor),
    Beta_Std = suppressWarnings(as.numeric(Beta_Std)),
    p_value = suppressWarnings(as.numeric(p_value))
  ) %>%
  filter(
    Response %in% beta_response_order,
    Predictor %in% predictor_order,
    is.finite(Beta_Std)
  ) %>%
  left_join(
    beta_response_labels %>% mutate(Response = as.character(Response)),
    by = "Response"
  ) %>%
  mutate(
    Response_label = factor(
      Response_label,
      levels = beta_response_labels$Response_label
    ),
    Predictor = factor(Predictor, levels = rev(predictor_order))
  )

all_sites_beta_grid <- expand_grid(
  Response = beta_response_order,
  Predictor = predictor_order
) %>%
  mutate(
    Response = as.character(Response),
    Predictor = as.character(Predictor)
  ) %>%
  left_join(
    all_sites_coefs %>%
      transmute(
        Response = as.character(Response),
        Predictor = as.character(Predictor),
        Beta_Std,
        p_value
      ),
    by = c("Response", "Predictor")
  ) %>%
  mutate(
    Response = factor(Response, levels = beta_response_order),
    Predictor = factor(Predictor, levels = rev(predictor_order)),
    coef_label = case_when(
      is.finite(Beta_Std) ~ paste0(
        sprintf("%.2f", Beta_Std),
        ifelse(is.finite(p_value) & p_value <= ALPHA, "*", "")
      ),
      TRUE ~ ""
    )
  )

p_beta <- ggplot(
  all_sites_beta_grid,
  aes(x = Response, y = Predictor, fill = Beta_Std)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = coef_label),
    size = FIG_TILE_TEXT_SIZE * 1.10,
    lineheight = 0.9
  ) +
  scale_fill_gradient2(
    low = "firebrick4",
    mid = "white",
    high = "dodgerblue3",
    midpoint = 0,
    limits = c(-1.0, 1.0),
    oob = scales::squish,
    na.value = "white",
    name = expression(beta)
  ) +
  scale_x_discrete(labels = function(x) {
    unname(response_display[as.character(x)])
  }) +
  scale_y_discrete(labels = predictor_axis_labels) +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    plot.margin = margin(
      FIG_LABEL_PLOT_MARGIN_PT,
      FIG_LABEL_PLOT_MARGIN_PT,
      FIG_LABEL_PLOT_MARGIN_PT,
      FIG_LABEL_PLOT_MARGIN_PT
    ),
    legend.position = "right"
  )

save_plot_pair(
  "fig8_eco_mlr_beta",
  p_beta,
  width = 8.8 * FIG_WIDTH_SCALE,
  height = 6.2 * FIG_HEIGHT_SCALE
)
unlink(file.path(plot_dir, "storage_ecovar_mlr_beta.png"))
unlink(file.path(plot_pdf_dir, "storage_ecovar_mlr_beta.pdf"))

# remove old companion r2 heatmap outputs from prior workflow
unlink(file.path(plot_dir, "storage_ecovar_mlr_r2_heatmap.png"))
unlink(file.path(plot_pdf_dir, "storage_ecovar_mlr_r2_heatmap.pdf"))

# keep site-level r2 table data for optional exports
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
      transmute(
        Site = as.character(Site),
        Response = as.character(Response),
        adj_r2 = R2_adj
      ),
    by = c("Site", "Response")
  ) %>%
  left_join(coverage_df, by = c("Site", "Response")) %>%
  mutate(
    Site = factor(Site, levels = site_order),
    Response = factor(Response, levels = rev(response_order)),
    r2_label = case_when(
      is.finite(adj_r2) ~ sprintf("%.2f", adj_r2),
      TRUE ~ "-"
    )
  )

build_pred_obs_plot <- function(pred_obs_path) {
  if (!file.exists(pred_obs_path)) {
    return(NULL)
  }

  pred_obs_df <- read_csv(pred_obs_path, show_col_types = FALSE) %>%
    transmute(
      Site = factor(as.character(Site), levels = site_order),
      Response = to_internal_response(as.character(Response)),
      Observed = suppressWarnings(as.numeric(Observed)),
      Predicted = suppressWarnings(as.numeric(Predicted))
    ) %>%
    filter(
      Response %in% pred_diag_response_order,
      is.finite(Observed),
      is.finite(Predicted)
    ) %>%
    mutate(
      Response_label = factor(
        unname(response_display[Response]),
        levels = pred_diag_response_labels$Response_label
      )
    )

  if (nrow(pred_obs_df) == 0) {
    return(NULL)
  }

  facet_stats_df <- pred_obs_df %>%
    group_by(Response_label) %>%
    summarise(
      rng_min = min(c(Observed, Predicted), na.rm = TRUE),
      rng_max = max(c(Observed, Predicted), na.rm = TRUE),
      r2 = ifelse(
        n() >= 2,
        suppressWarnings(cor(Observed, Predicted, use = "complete.obs")^2),
        NA_real_
      ),
      .groups = "drop"
    )

  line_label_df <- facet_stats_df %>%
    mutate(
      x = rng_min + 0.9 * (rng_max - rng_min),
      y = x,
      line_label = "1:1"
    )

  r2_label_df <- facet_stats_df %>%
    mutate(
      r2_label = case_when(
        is.finite(r2) ~ sprintf("R^2 == %.2f", r2),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(r2_label)) %>%
    mutate(
      x = rng_min + 0.04 * (rng_max - rng_min),
      y = rng_max - 0.04 * (rng_max - rng_min)
    ) %>%
    dplyr::select(Response_label, x, y, r2_label)

  ggplot(
    pred_obs_df,
    aes(x = Observed, y = Predicted, color = Site)
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.6,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = line_label_df,
      aes(x = x, y = y, label = line_label),
      inherit.aes = FALSE,
      color = "grey30",
      hjust = 0,
      vjust = -0.35,
      size = FIG_AXIS_TEXT_SIZE / 2.7
    ) +
    geom_text(
      data = r2_label_df,
      aes(x = x, y = y, label = r2_label),
      inherit.aes = FALSE,
      parse = TRUE,
      color = "black",
      hjust = 0,
      vjust = 1,
      size = FIG_AXIS_TEXT_SIZE / 2.8
    ) +
    geom_point(
      size = FIG_POINT_SIZE_MED,
      alpha = 0.8
    ) +
    facet_wrap(
      ~Response_label,
      ncol = 1,
      scales = "free",
      labeller = labeller(Response_label = pred_diag_response_labels_panel),
      axes = "margins",
      axis.labels = "margins"
    ) +
    scale_color_manual(
      values = SITE_COLORS,
      drop = FALSE
    ) +
    labs(
      x = "Observed",
      y = "Predicted",
      color = "Site"
    ) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      strip.text = element_text(size = FIG_STRIP_TEXT_SIZE, hjust = 0),
      legend.position = "bottom",
      legend.title = element_text(size = FIG_AXIS_TEXT_SIZE),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE - 1)
    )
}

# predicted versus observed is treated as a standard companion figure
p_pred_obs <- build_pred_obs_plot(pred_obs_file)
if (!is.null(p_pred_obs)) {
  save_plot_pair(
    "FigSX_eco_observed_v_predicted",
    p_pred_obs,
    width = 9 * FIG_WIDTH_SCALE,
    height = 10 * FIG_HEIGHT_SCALE,
    to_supp = TRUE
  )
  unlink(file.path(supp_plot_dir, "storage_ecovar_mlr_predicted_vs_observed.png"))
  unlink(file.path(supp_plot_pdf_dir, "storage_ecovar_mlr_predicted_vs_observed.pdf"))
}

# pooled diagnostics heatmap with p values for residual checks
diagnostics_file <- file.path(
  output_dir,
  "storage_ecovar_mlr_all_sites_diagnostics.csv"
)
if (file.exists(diagnostics_file)) {
  diag_df <- read_csv(diagnostics_file, show_col_types = FALSE) %>%
    transmute(
      Response = to_internal_response(as.character(Response)),
      shapiro_p = suppressWarnings(as.numeric(shapiro_p)),
      ncv_p = suppressWarnings(as.numeric(ncv_p))
    ) %>%
    filter(Response %in% pred_diag_response_order)

  diag_long <- diag_df %>%
    pivot_longer(
      cols = c(shapiro_p, ncv_p),
      names_to = "Diagnostic",
      values_to = "p_value"
    ) %>%
    mutate(
      Diagnostic = recode(
        Diagnostic,
        shapiro_p = "Normality (Shapiro p)",
        ncv_p = "Homoscedasticity (NCV p)"
      ),
      Response = factor(Response, levels = pred_diag_response_order),
      p_label = case_when(
        is.finite(p_value) ~ sprintf("%.3f", p_value),
        TRUE ~ "-"
      ),
      fail = is.finite(p_value) & (p_value <= 0.05)
    )

  p_diag <- ggplot(
    diag_long,
    aes(x = Response, y = Diagnostic, fill = p_value)
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = p_label), size = FIG_TILE_TEXT_SIZE) +
    geom_point(
      data = diag_long %>% filter(fail),
      aes(x = Response, y = Diagnostic),
      inherit.aes = FALSE,
      shape = 4,
      size = 2.0,
      stroke = 0.75,
      color = "black"
    ) +
    scale_fill_gradientn(
      colors = c("firebrick4", "goldenrod2", "seagreen4"),
      values = scales::rescale(c(0, 0.05, 1)),
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "grey90",
      name = "p-value"
    ) +
    scale_x_discrete(labels = function(x) {
      unname(response_display[as.character(x)])
    }) +
    labs(
      x = "Response",
      y = "Diagnostic"
    ) +
    theme_pub() +
    theme(
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE)
    )

  save_plot_pair(
    "storage_ecovar_mlr_diagnostics",
    p_diag,
    width = 8 * FIG_WIDTH_SCALE,
    height = 4.5 * FIG_HEIGHT_SCALE,
    to_supp = TRUE
  )
}

if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  table_dir <- OUT_TABLES_MLR_DIR
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

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

  format_table_label <- function(x) {
    out <- as.character(x)
    out <- gsub("P_WetSeason", "Pws", out, fixed = TRUE)
    out <- gsub("_", "", out, fixed = TRUE)
    out
  }

  perf_export <- perf_df %>%
    mutate(
      Response = format_table_label(as.character(Response)),
      Predictors_Final = format_table_label(as.character(Predictors_Final))
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
      delta_rmse_loocv_mean_runs_minus_model = signif(
        delta_rmse_loocv_mean_runs_minus_model,
        3
      )
    )

  save_csv_safe(
    perf_export,
    file.path(table_dir, "storage_ecovar_mlr_model_perf.csv")
  )
  save_csv_safe(
    r2_heat_df %>%
      transmute(
        site = as.character(Site),
        response = format_table_label(as.character(Response)),
        adj_r2 = signif(adj_r2, 3),
        model_status,
        reason_not_fit,
        n_response
      ),
    file.path(table_dir, "storage_ecovar_mlr_r2_heatmap_table.csv")
  )

  coef_export <- coef_df %>%
    mutate(
      Response = format_table_label(as.character(Response)),
      Predictor = format_table_label(as.character(Predictor)),
      Beta_Std = signif(Beta_Std, 3),
      p_value = signif(p_value, 3),
      VIF = signif(VIF, 3),
      R2 = signif(R2, 3),
      R2_adj = signif(R2_adj, 3),
      RMSE = signif(RMSE, 3),
      AICc = signif(AICc, 3)
    )

  save_csv_safe(
    coef_export,
    file.path(table_dir, "storage_ecovar_mlr_coef.csv")
  )

  manuscript_table <- beta_plot_df %>%
    transmute(
      site = as.character(Site),
      response = format_table_label(as.character(Response)),
      predictor = format_table_label(as.character(Predictor)),
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
}
