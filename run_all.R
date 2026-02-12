#!/usr/bin/env Rscript

find_repo_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(10)) {
    has_config <- file.exists(file.path(cur, "config.R"))
    has_metrics <- dir.exists(file.path(cur, "01_metrics"))
    has_stats <- dir.exists(file.path(cur, "02_stats"))
    if (has_config && has_metrics && has_stats) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  return(normalizePath(start_dir, winslash = "/", mustWork = FALSE))
}

script_path <- tryCatch({
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE)
}, error = function(e) NA_character_)

if (is.na(script_path) || script_path == "") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
  }
}

start_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_root <- find_repo_root(start_dir)

run_script <- function(path) {
  full <- file.path(repo_root, path)
  if (!file.exists(full)) stop("Missing script: ", full)
  status <- system2("Rscript", shQuote(full), stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Script failed: ", path)
  }
}

# Preflight.
run_script("check_inputs.R")

# Metrics.
metric_scripts <- c(
  "01_metrics/00_create_master_hydrometric_dataset.R",
  "01_metrics/01_ds_rbi_rcs.R",
  "01_metrics/02_ds_sd_fdc.R",
  "01_metrics/03_ms_load_isotope.R",
  "01_metrics/04_ms_chs.R",
  "01_metrics/05_eds_wb_drawdown.R",
  "01_metrics/06_eco_vars.R",
  "01_metrics/07_aggregate_metrics.R",
  "01_metrics/08_data_availability.R"
)
for (s in metric_scripts) run_script(s)

# Stats.
stats_scripts <- c(
  "02_stats/anova_tukey.R",
  "02_stats/pca.R",
  "02_stats/catchment_storage_mlr.R",
  "02_stats/predict_thermal_lowflow.R"
)
for (s in stats_scripts) run_script(s)

# Plots.
plot_scripts <- c(
  "03_plots/anova_tukey.R",
  "03_plots/pca.R",
  "03_plots/correlations.R",
  "03_plots/catchment_controls_mlr.R",
  "03_plots/predict_thermal_lowflow.R",
  "03_plots/recession_q_dqdt.R",
  "03_plots/hydrometric_metric_summaries.R",
  "03_plots/met_context.R",
  "03_plots/main_figures.R",
  "03_plots/output_manifest.R"
)
for (s in plot_scripts) run_script(s)

# Post-run checks.
run_script("verify_outputs.R")

message("Run complete.")
