# inputs: no direct csv file reads in this script.
# author: sidney bush
# date: 2026-02-13

find_repo_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(10)) {
    has_config <- file.exists(file.path(cur, "config.R"))
    has_metrics <- dir.exists(file.path(cur, "01_storage_calcs"))
    has_stats <- dir.exists(file.path(cur, "02_analysis"))
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

# prefer explicit override when launched from ide with unusual cwd.
env_repo_root <- Sys.getenv("HJA_REPO_DIR", unset = "")
if (nzchar(env_repo_root)) {
  repo_root <- normalizePath(env_repo_root, winslash = "/", mustWork = FALSE)
} else {
  # try to discover this script path robustly across rscript/source/ide run modes.
  script_path <- NA_character_

  # look through call frames for an ofile ending in run_all.r
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

  # fall back to --file if running via rscript
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

if (!file.exists(file.path(repo_root, "helpers", "check_inputs.R"))) {
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

# preflight.
run_script("helpers/check_inputs.R")

# metrics preprocessing (daily met+q and et support tables).
metric_preprocess_scripts <- c(
  "00_data_preprocessing/create_hydromet_master.R",
  "00_data_preprocessing/prep_PT_methods.R",
  "00_data_preprocessing/prep_Hamon_methods.R"
)
for (s in metric_preprocess_scripts) {
  run_script(s)
}

# metrics.
metric_scripts <- c(
  "01_storage_calcs/calc_dynamic_storage.R",
  "01_storage_calcs/calc_mobile_storage.R",
  "01_storage_calcs/calc_eco_response.R",
  "01_storage_calcs/aggregate_all.R"
)
for (s in metric_scripts) {
  run_script(s)
}

# strict unit-consistency gate (fail fast before stats/plots).
run_script("helpers/check_units_consistency.R")

# stats.
# deprecated unified-framework climate/sensitivity scripts are kept under
# deprecated/02_analysis and are intentionally excluded from the core run.
stats_scripts <- c(
  "02_analysis/storage_sum_stats.R",
  "02_analysis/pca.R",
  "02_analysis/conceptual_diagram_calc.R",
  "02_analysis/mlr_catchment_char.R",
  "02_analysis/mlr_eco_response.R",
  "02_analysis/mlr_tables.R"
)
for (s in stats_scripts) {
  run_script(s)
}

# plots: core manuscript set.
# deprecated unified-framework climate-variability plotting is kept under
# deprecated/03_plots and is intentionally excluded from the core run.
plot_scripts_core <- c(
  "03_plots/Fig3_ds_pca_annual.R",
  "03_plots/Fig6_dynamic_mobile_corr.R",
  "03_plots/Fig9_conceptual_diagram.R",
  "03_plots/Fig7_catchment_mlr_beta.R",
  "03_plots/Fig8_eco_mlr_beta.R",
  "03_plots/Fig2_4_5_storage.R"
)
for (s in plot_scripts_core) {
  run_script(s)
}

# supplementary figures are managed in one script for easier triage.
run_script("03_plots/supplementary.R")

# post-run checks.
run_script("helpers/verify_outputs.R")
