# -----------------------------------------------------------------------------
# Supplemental Figure: Predictor vs Response Pairs for Both MLR Workflows
# Faceting is by site only (per user requirement).
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
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
if (!file.exists(config_path)) config_path <- file.path(dirname(script_dir), "config.R")
if (!file.exists(config_path)) config_path <- file.path(getwd(), "config.R")
if (!file.exists(config_path)) stop("config.R not found.")
source(config_path)

plot_dir <- file.path(FIGURES_DIR, "supp", "mlr")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Keep only the two site-faceted figures for clarity.
unlink(file.path(plot_dir, c(
  "mlr_predictor_response_pairs_combined.png",
  "mlr_predictor_response_pairs_combined.pdf"
)))

safe_scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || (rng[2] - rng[1]) == 0) return(rep(NA_real_, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# -----------------------------------------------------------------------------
# ECO MLR: annual points by site-year for fitted predictor-response pairs
# -----------------------------------------------------------------------------

eco_results_path <- file.path(OUT_MODELS_STORAGE_ECOVAR_MLR_DIR, "storage_ecovar_mlr_results.csv")
annual_path <- file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE)
if (!file.exists(annual_path)) annual_path <- file.path(OUT_MASTER_DIR, LEGACY_ANNUAL_FILE)

eco_plot <- NULL
if (file.exists(eco_results_path) && file.exists(annual_path)) {
  eco_results <- read_csv(eco_results_path, show_col_types = FALSE) %>%
    distinct(Response, Predictor)
  annual <- read_csv(annual_path, show_col_types = FALSE) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC)

  eco_long <- annual %>%
    select(any_of(c("site", "year", unique(eco_results$Response), unique(eco_results$Predictor)))) %>%
    pivot_longer(cols = any_of(unique(eco_results$Predictor)), names_to = "Predictor", values_to = "x") %>%
    pivot_longer(cols = any_of(unique(eco_results$Response)), names_to = "Response", values_to = "y") %>%
    semi_join(eco_results, by = c("Response", "Predictor")) %>%
    filter(is.finite(x), is.finite(y)) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
    mutate(combo = paste(Predictor, Response, sep = " -> ")) %>%
    group_by(Predictor, Response) %>%
    mutate(
      x_scaled = safe_scale01(x),
      y_scaled = safe_scale01(y)
    ) %>%
    ungroup() %>%
    filter(is.finite(x_scaled), is.finite(y_scaled))

  eco_line_data <- eco_long %>%
    count(site, combo, name = "n_obs") %>%
    filter(n_obs >= 3) %>%
    inner_join(eco_long, by = c("site", "combo"))

  eco_plot <- ggplot(eco_long, aes(x = x_scaled, y = y_scaled, color = combo)) +
    geom_point(size = 1.2, alpha = 0.7) +
    geom_smooth(
      data = eco_line_data,
      method = "lm",
      se = FALSE,
      linewidth = 0.4,
      alpha = 0.4,
      show.legend = FALSE
    ) +
    facet_wrap(~site, ncol = 2) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 2, alpha = 1))) +
    labs(x = "Predictor (scaled 0-1)", y = "Response (scaled 0-1)") +
    theme_pub() +
    theme(
      strip.text = element_text(size = FIG_STRIP_TEXT_SIZE),
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.position = "right",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(plot_dir, "storage_ecovar_mlr_predictor_response_pairs.png"),
    eco_plot,
    width = 13 * FIG_WIDTH_SCALE,
    height = 10 * FIG_HEIGHT_SCALE,
    dpi = 300,
    create.dir = TRUE
  )
  ggsave(
    file.path(plot_dir, "storage_ecovar_mlr_predictor_response_pairs.pdf"),
    eco_plot,
    width = 13 * FIG_WIDTH_SCALE,
    height = 10 * FIG_HEIGHT_SCALE,
    create.dir = TRUE
  )
}

# -----------------------------------------------------------------------------
# Catchment MLR: site-level points for fitted predictor-outcome pairs
# -----------------------------------------------------------------------------

catch_results_path <- file.path(OUT_MODELS_WATERSHED_CHAR_STORAGE_MLR_DIR, "watershed_char_storage_mlr_results_strict.csv")
site_path <- file.path(OUT_MASTER_DIR, MASTER_SITE_FILE)
if (!file.exists(site_path)) site_path <- file.path(OUT_MASTER_DIR, LEGACY_SITE_FILE)

catch_plot <- NULL
if (file.exists(catch_results_path) && file.exists(site_path)) {
  catch_results <- read_csv(catch_results_path, show_col_types = FALSE) %>%
    distinct(Outcome, Predictor)
  site_data <- read_csv(site_path, show_col_types = FALSE) %>%
    mutate(site = standardize_site_code(site))

  catch_long <- site_data %>%
    select(any_of(c("site", unique(catch_results$Outcome), unique(catch_results$Predictor)))) %>%
    pivot_longer(cols = any_of(unique(catch_results$Predictor)), names_to = "Predictor", values_to = "x") %>%
    pivot_longer(cols = any_of(unique(catch_results$Outcome)), names_to = "Outcome", values_to = "y") %>%
    semi_join(catch_results, by = c("Outcome", "Predictor")) %>%
    filter(is.finite(x), is.finite(y), site %in% SITE_ORDER_HYDROMETRIC) %>%
    mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
    mutate(combo = paste(Predictor, Outcome, sep = " -> ")) %>%
    group_by(Predictor, Outcome) %>%
    mutate(
      x_scaled = safe_scale01(x),
      y_scaled = safe_scale01(y)
    ) %>%
    ungroup() %>%
    filter(is.finite(x_scaled), is.finite(y_scaled))

  catch_plot <- ggplot(catch_long, aes(x = x_scaled, y = y_scaled, color = combo)) +
    geom_point(size = 2.0, alpha = 0.9) +
    facet_wrap(~site, ncol = 2) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.4, alpha = 1))) +
    labs(x = "Predictor (scaled 0-1)", y = "Outcome (scaled 0-1)") +
    theme_pub() +
    theme(
      strip.text = element_text(size = FIG_STRIP_TEXT_SIZE),
      axis.text = element_text(size = FIG_AXIS_TEXT_SIZE),
      axis.title = element_text(size = FIG_AXIS_TITLE_SIZE),
      legend.position = "right",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT, FIG_LABEL_PLOT_MARGIN_PT)
    ) +
    coord_cartesian(clip = FIG_LABEL_CLIP)

  ggsave(
    file.path(plot_dir, "watershed_char_storage_mlr_predictor_response_pairs.png"),
    catch_plot,
    width = 13 * FIG_WIDTH_SCALE,
    height = 10 * FIG_HEIGHT_SCALE,
    dpi = 300,
    create.dir = TRUE
  )
  ggsave(
    file.path(plot_dir, "watershed_char_storage_mlr_predictor_response_pairs.pdf"),
    catch_plot,
    width = 13 * FIG_WIDTH_SCALE,
    height = 10 * FIG_HEIGHT_SCALE,
    create.dir = TRUE
  )
}

# No combined mega-figure by design.
