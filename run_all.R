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
    if (identical(parent, cur)) {
      break
    }
    cur <- parent
  }
  return(normalizePath(start_dir, winslash = "/", mustWork = FALSE))
}

# Prefer explicit override when launched from IDE with unusual cwd.
env_repo_root <- Sys.getenv("HJA_REPO_DIR", unset = "")
if (nzchar(env_repo_root)) {
  repo_root <- normalizePath(env_repo_root, winslash = "/", mustWork = FALSE)
} else {
  # Try to discover this script path robustly across Rscript/source/IDE run modes.
  script_path <- NA_character_

  # Look through call frames for an ofile ending in run_all.R
  frame_ofiles <- unlist(lapply(sys.frames(), function(fr) {
    tryCatch(fr$ofile, error = function(e) NA_character_)
  }))
  frame_ofiles <- frame_ofiles[
    is.character(frame_ofiles) & !is.na(frame_ofiles) & nzchar(frame_ofiles)
  ]
  if (length(frame_ofiles) > 0) {
    frame_hits <- frame_ofiles[basename(frame_ofiles) == "run_all.R"]
    if (length(frame_hits) > 0) {
      script_path <- normalizePath(
        frame_hits[1],
        winslash = "/",
        mustWork = FALSE
      )
    }
  }

  # Fall back to --file if running via Rscript
  if (is.na(script_path) || script_path == "") {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      script_path <- normalizePath(
        sub("^--file=", "", file_arg[1]),
        winslash = "/",
        mustWork = FALSE
      )
    }
  }

  start_dir <- if (!is.na(script_path) && nzchar(script_path)) {
    dirname(script_path)
  } else {
    getwd()
  }
  repo_root <- find_repo_root(start_dir)
}

if (!file.exists(file.path(repo_root, "check_inputs.R"))) {
  stop(
    paste0(
      "Could not locate repo root from current session.\n",
      "Set working directory to the repo, or set HJA_REPO_DIR to repo path.\n",
      "Current resolved repo_root: ",
      repo_root
    )
  )
}

run_script <- function(path) {
  full <- file.path(repo_root, path)
  if (!file.exists(full)) {
    stop("Missing script: ", full)
  }
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
for (s in metric_scripts) {
  run_script(s)
}

# Stats.
stats_scripts <- c(
  "02_stats/anova_tukey.R",
  "02_stats/pca.R",
  "02_stats/watershed_char_storage_mlr.R",
  "02_stats/storage_ecovar_mlr.R"
)
for (s in stats_scripts) {
  run_script(s)
}

# Plots: core manuscript set.
plot_scripts_core <- c(
  "03_plots/anova_tukey.R",
  "03_plots/pca.R",
  "03_plots/correlations.R",
  "03_plots/watershed_char_storage_mlr.R",
  "03_plots/storage_ecovar_mlr.R",
  "03_plots/met_context.R",
  "03_plots/dynamic_mobile_chs_figures.R"
)
for (s in plot_scripts_core) {
  run_script(s)
}

# Optional plot diagnostics/exploratory set (off by default).
include_optional_plots <- tolower(Sys.getenv(
  "HJA_INCLUDE_OPTIONAL_PLOTS",
  unset = "false"
)) ==
  "true"
if (include_optional_plots) {
  plot_scripts_optional <- c(
    "03_plots/recession_q_dqdt.R",
    "03_plots/hydrometric_metric_summaries.R",
    "03_plots/mlr_predictor_response_pairs.R"
  )
  for (s in plot_scripts_optional) {
    run_script(s)
  }
}

# Post-run checks.
run_script("verify_outputs.R")
