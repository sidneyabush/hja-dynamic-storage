# plots: pooled eco-response predictor effects (all-sites model).
# inputs: output_dir/storage_ecovar_mlr_all_sites_results.csv.
# author: sidney bush
# date: 2026-02-24

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

rm(list = ls())

# load project config
source("config.R")

output_dir <- OUT_MODELS_STORAGE_ECOVAR_MLR_DIR
file_prefix <- "storage_ecovar_mlr_all_sites"
ALPHA <- 0.05

plot_dir <- MS_FIG_MAIN_DIR
plot_pdf_dir <- plot_dir
supp_plot_dir <- SUPP_LEGACY_DIR
supp_plot_pdf_dir <- SUPP_LEGACY_PDF_DIR

dir_targets <- c(plot_dir, plot_pdf_dir, supp_plot_dir, supp_plot_pdf_dir)
for (d in dir_targets) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

results_file <- file.path(output_dir, paste0(file_prefix, "_results.csv"))
summary_file <- file.path(output_dir, paste0(file_prefix, "_summary.csv"))

if (!file.exists(results_file)) {
  stop("Missing file: ", results_file)
}
if (!file.exists(summary_file)) {
  stop("Missing file: ", summary_file)
}

to_internal_response <- function(x) {
  out <- as.character(x)
  out[out == "Q7Q5"] <- "Q_7Q5"
  out[out == "T7DMax"] <- "T_7DMax"
  out
}

to_internal_predictor <- function(x) {
  out <- as.character(x)
  out[out == "Pws"] <- "P_WetSeason"
  out
}

response_order <- c("Q_7Q5", "T_7DMax")
predictor_order <- c(
  "P_WetSeason",
  STORAGE_METRIC_ORDER[STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD", "WB", "CHS")]
)
response_display <- c(
  "Q_7Q5" = "Q7Q5",
  "T_7DMax" = "T7DMax"
)
predictor_display <- c(
  "P_WetSeason" = "P[ws]",
  "RBI" = "RBI",
  "RCS" = "RCS",
  "FDC" = "FDC",
  "SD" = "SD",
  "WB" = "WB",
  "CHS" = "CHS"
)

format_model_p <- function(p) {
  ifelse(
    is.finite(p),
    format(signif(p, 3), scientific = TRUE, trim = TRUE),
    "NA"
  )
}

model_coefs <- read_csv(results_file, show_col_types = FALSE) %>%
  mutate(
    Response = to_internal_response(Response),
    Predictor = to_internal_predictor(Predictor)
  ) %>%
  filter(Response %in% response_order, Predictor %in% predictor_order)

model_summary <- read_csv(summary_file, show_col_types = FALSE) %>%
  mutate(Response = to_internal_response(Response)) %>%
  filter(Response %in% response_order)

if (!("model_p_global" %in% names(model_summary))) {
  n_pred <- vapply(
    as.character(model_summary$Predictors_Final),
    function(txt) {
      if (is.na(txt) || !nzchar(txt) || identical(txt, "-")) return(0L)
      parts <- trimws(strsplit(txt, ";", fixed = TRUE)[[1]])
      as.integer(sum(nzchar(parts)))
    },
    integer(1)
  )
  r2 <- suppressWarnings(as.numeric(model_summary$R2))
  n_obs <- suppressWarnings(as.numeric(model_summary$n))
  f_stat <- ifelse(
    n_pred > 0 & is.finite(r2) & is.finite(n_obs) & (n_obs - n_pred - 1) > 0,
    (r2 / n_pred) / ((1 - r2) / (n_obs - n_pred - 1)),
    NA_real_
  )
  model_summary$model_p_global <- ifelse(
    n_pred > 0 & is.finite(f_stat),
    pf(f_stat, n_pred, n_obs - n_pred - 1, lower.tail = FALSE),
    NA_real_
  )
}

response_meta <- model_summary %>%
  transmute(
    Response = as.character(Response),
    Response_label = unname(response_display[Response]),
    R2_adj = suppressWarnings(as.numeric(R2_adj)),
    model_p_global = suppressWarnings(as.numeric(model_p_global))
  ) %>%
  mutate(
    model_p_label = format_model_p(model_p_global),
    r2_label = ifelse(is.finite(R2_adj), sprintf("%.2f", R2_adj), "NA"),
    Response_facet = paste0(Response_label, " | adj R2 ", r2_label, " | model p ", model_p_label)
  )

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

unlink(file.path(plot_dir, paste0(file_prefix, "_beta_heatmap.png")))
unlink(file.path(plot_pdf_dir, paste0(file_prefix, "_beta_heatmap.pdf")))

beta_bar_df <- model_coefs %>%
  transmute(
    Response = as.character(Response),
    Predictor = as.character(Predictor),
    Beta_Std = suppressWarnings(as.numeric(Beta_Std)),
    p_value = suppressWarnings(as.numeric(p_value))
  ) %>%
  filter(is.finite(Beta_Std)) %>%
  left_join(
    response_meta %>% select(Response, Response_facet),
    by = "Response"
  ) %>%
  mutate(
    Response_label = factor(
      Response_facet,
      levels = unname(response_meta$Response_facet[match(response_order, response_meta$Response)])
    ),
    Predictor_label = unname(predictor_display[Predictor]),
    sig_label = ifelse(is.finite(p_value) & p_value <= ALPHA, "*", "")
  ) %>%
  group_by(Response_label) %>%
  arrange(Beta_Std, .by_group = TRUE) %>%
  mutate(Predictor_plot = factor(Predictor_label, levels = unique(Predictor_label))) %>%
  ungroup()

response_labels_panel <- make_panel_label_map(
  unname(response_meta$Response_facet[match(response_order, response_meta$Response)])
)

p_beta_bar <- ggplot(beta_bar_df, aes(x = Beta_Std, y = Predictor_plot, fill = Beta_Std)) +
  geom_col(width = 0.75, color = "white") +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_text(
    aes(label = sig_label, hjust = ifelse(Beta_Std >= 0, -0.3, 1.3)),
    size = FIG_TILE_TEXT_SIZE * 1.1
  ) +
  facet_wrap(
    ~Response_label,
    ncol = 1,
    scales = "free_y",
    labeller = labeller(Response_label = response_labels_panel)
  ) +
  scale_fill_gradient2(
    low = "firebrick4",
    mid = "white",
    high = "dodgerblue3",
    midpoint = 0,
    limits = c(-1.0, 1.0),
    oob = scales::squish,
    guide = "none"
  ) +
  labs(
    x = "Standardized coefficient (Beta)",
    y = "Predictor"
  ) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  theme_pub() +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
    strip.text = element_text(size = FIG_STRIP_TEXT_SIZE, hjust = 0)
  ) +
  coord_cartesian(clip = FIG_LABEL_CLIP)

save_plot_safe(
  file.path(supp_plot_dir, paste0(file_prefix, "_beta_bar.png")),
  p_beta_bar,
  width = 8 * FIG_WIDTH_SCALE,
  height = 8 * FIG_HEIGHT_SCALE,
  dpi = 300
)
save_plot_safe(
  file.path(supp_plot_pdf_dir, paste0(file_prefix, "_beta_bar.pdf")),
  p_beta_bar,
  width = 8 * FIG_WIDTH_SCALE,
  height = 8 * FIG_HEIGHT_SCALE
)
