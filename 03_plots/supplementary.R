# make the supplement figure used in the final SI
# inputs: figure scripts in this folder
# outputs: supplement figure files in ms_materials/supp/

rm(list = ls())
source("config.R")

supp_dir <- MS_FIG_SUPP_DIR
supp_pdf_dir <- MS_FIG_SUPP_PDF_DIR
for (d in c(supp_dir, supp_pdf_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

run_supp_script <- function(filename) {
  script_path <- file.path(REPO_DIR, "03_plots", filename)
  if (!file.exists(script_path)) {
    stop("Missing supplementary script: ", script_path)
  }
  status <- system2("Rscript", shQuote(script_path), stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Supplementary script failed: ", filename)
  }
}

# keep supplement outputs in separate scripts when that keeps things clear
run_supp_script("FigS1_met_context.R")
# remove older supplement files that are not part of the final SI
legacy_supp_files <- c(
  "FigS2_dynamic_storage_corr.png",
  "FigS3_mobile_storage_corr.png",
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
legacy_supp_pdf_files <- sub("\\.(png|csv)$", ".pdf", legacy_supp_files[grepl("\\.(png|csv)$", legacy_supp_files)])
legacy_supp_pdf_files <- unique(c(
  "FigS2_dynamic_storage_corr.pdf",
  "FigS3_mobile_storage_corr.pdf",
  "FigS4_catchment_mlr_diagnostics.pdf",
  "FigS4_dynamic_mobile_scatter_matrix.pdf",
  "FigS5_eco_mlr_diagnostics.pdf",
  "FigS6_dynamic_mobile_scatter_matrix.pdf",
  "FigSX_dynamic_metrics_corr.pdf",
  "FigSX_mobile_metrics_corr.pdf",
  "FigSX_pca_Pws_anomaly.pdf",
  "FigSX_pca_scree.pdf",
  "FigSX_chs_ec_vs_ca_by_site.pdf",
  "FigSX_chs_ec_vs_ca_overall.pdf"
))

unlink(file.path(supp_dir, legacy_supp_files))
unlink(file.path(supp_pdf_dir, legacy_supp_pdf_files))
