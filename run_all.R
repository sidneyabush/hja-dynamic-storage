# run the full analysis after running install_packages.R
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

if (!file.exists(file.path(repo_dir, "config.R"))) {
  stop(
    "Could not find the repository folder. Run this script from the repository folder ",
    "or set HJA_REPO_DIR to that folder."
  )
}

setwd(repo_dir)

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
  if (!file.exists(full)) {
    stop("Missing script: ", full)
  }
  status <- system2("Rscript", shQuote(full), stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Script failed: ", path)
  }
}

run_scripts <- function(paths) {
  for (path in paths) {
    run_script(path)
  }
}

# check input files before running the workflow
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
  "02_analysis/mlr_catchment_char.R",
  "02_analysis/mlr_eco_response.R",
  "02_analysis/mtt_sensitivity.R"
)
run_scripts(stats_scripts)

# make the main text figures
plot_scripts_core <- c(
  "03_figs/Fig2_Fig3_storage_metrics.R",
  "03_figs/Fig4_catchment_controls.R",
  "03_figs/Fig5_ecological_response_models.R",
  "03_figs/Fig6_observed_predicted_ecological_responses.R",
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

# remove old supporting information outputs
old_supp_files <- c(
  "FigS4_catchment_mlr_diagnostics.png",
  "FigS4_dynamic_mobile_scatter_matrix.png",
  "FigS5_eco_mlr_diagnostics.png",
  "FigS6_dynamic_mobile_scatter_matrix.png",
  "FigSX_dynamic_metrics_corr.png",
  "FigSX_mobile_metrics_corr.png",
  "FigSX_pca_Pws_anomaly.png",
  "FigSX_chs_ec_vs_ca_by_site.png",
  "FigSX_chs_ec_vs_ca_overall.png",
  "Table_SX_mlr_model_diagnostics.csv",
  "TableS9_mlr_model_diagnostics.csv"
)
old_supp_pdf_files <- unique(c(
  sub("\\.(png|csv)$", ".pdf", old_supp_files[grepl("\\.(png|csv)$", old_supp_files)]),
  "FigSX_pca_scree.pdf"
))
unlink(file.path(MS_FIG_SUPP_DIR, old_supp_files))
unlink(file.path(MS_FIG_SUPP_PDF_DIR, old_supp_pdf_files))

# write the supporting information tables
run_script("04_tables/SI_tables_S7_S12.R")

# check that the outputs were written
verify_outputs()

# close any open plots
try(grDevices::graphics.off(), silent = TRUE)

# remove any accidental Rplots.pdf
if (file.exists("Rplots.pdf")) {
  unlink("Rplots.pdf")
}
