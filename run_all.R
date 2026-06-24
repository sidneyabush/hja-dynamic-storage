# run the full analysis after running install_packages.R

# inputs:
# config.R
# figure_functions.R
# workflow_functions.R
# hydromet_functions.R
# model_functions.R
# workflow input files listed in each script

# outputs:
# outputs/*
# figs_tables_pub/*

# author: Sidney Bush
# date: 2026-02-13

# find the repository folder
repo_dir <- Sys.getenv("HJA_REPO_DIR", unset = "")
if (!nzchar(repo_dir)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  repo_dir <- if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
  } else {
    getwd()
  }
}

repo_dir <- normalizePath(repo_dir, winslash = "/", mustWork = FALSE)

setwd(repo_dir)

# require librarian because each script uses librarian::shelf
if (!requireNamespace("librarian", quietly = TRUE)) {
  stop(
    "Missing R package: librarian. Run Rscript install_packages.R once, ",
    "then rerun Rscript run_all.R.",
    call. = FALSE
  )
}

source(file.path(repo_dir, "config.R"))
make_output_dirs()

run_script <- function(path) {
  full <- file.path(repo_dir, path)
  status <- system2("Rscript", shQuote(full), stdout = "", stderr = "")

  # stop the workflow at the script that returns an error
  if (!identical(status, 0L)) {
    stop("Script failed: ", path)
  }
}

run_scripts <- function(paths) {
  for (path in paths) {
    run_script(path)
  }
}

# verify required input columns before running the workflow
check_inputs()

# build the daily catchment forcing and ET tables
metric_preprocess_scripts <- c(
  "00_data_preprocessing/create_hydromet_master.R",
  "00_data_preprocessing/prep_PT_methods.R",
  "00_data_preprocessing/prep_Hamon_methods.R"
)
run_scripts(metric_preprocess_scripts)

# calculate annual storage and ecological response metrics
metric_scripts <- c(
  "01_storage_calcs/calc_dynamic_storage.R",
  "01_storage_calcs/calc_mobile_storage.R",
  "01_storage_calcs/calc_eco_response.R",
  "01_storage_calcs/aggregate_all.R"
)
run_scripts(metric_scripts)

# stop here if the units do not line up
run_script("02_analysis/check_units_consistency.R")

# run summary stats, PCA, regression models, and MTT sensitivity
stats_scripts <- c(
  "02_analysis/storage_sum_stats.R",
  "02_analysis/pca.R",
  "02_analysis/unified_framework_calc.R",
  "02_analysis/mlr_catchment_eco.R",
  "02_analysis/mtt_sensitivity.R"
)
run_scripts(stats_scripts)

# make the main text figures
plot_scripts_core <- c(
  "03_figs/Fig2_Fig3_storage_metrics.R",
  "03_figs/Fig4_Fig5_Fig6_model_figures.R",
  "03_figs/Fig7_dynamic_mobile_framework.R"
)
run_scripts(plot_scripts_core)

# make the supporting information figures
plot_scripts_supp <- c(
  "03_figs/FigS1_met_context.R",
  "03_figs/FigS2_S3_storage_corr.R",
  "03_figs/FigS4_dynamic_mobile_corr.R"
)
run_scripts(plot_scripts_supp)

# write the supporting information tables
run_script("04_tables/SI_tables_S1_S6.R")
run_script("04_tables/SI_tables_S7_S12.R")

# verify required manuscript figures and supporting information tables
verify_outputs()

# close any open plots
try(grDevices::graphics.off(), silent = TRUE)

# remove accidental Rplots.pdf
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
