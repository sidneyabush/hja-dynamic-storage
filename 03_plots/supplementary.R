# supplementary exports used in the final SI only

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

rm(list = ls())
source("config.R")

supp_dir <- MS_FIG_SUPP_DIR
supp_pdf_dir <- MS_FIG_SUPP_PDF_DIR
for (d in c(supp_dir, supp_pdf_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Figure S1 ----
# kept in a dedicated script but dispatched from here so the SI workflow
# has a single entry point for supplementary figure/table generation.
figs1_script <- file.path(REPO_DIR, "03_plots", "FigS1_met_context.R")
if (file.exists(figs1_script)) {
  status <- system2("Rscript", shQuote(figs1_script), stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    warning("FigS1_met_context.R failed when called from supplementary.R")
  }
}

# ---- Table S8 ----
# residual diagnostics table used in the final SI.
catch_diag_file <- file.path(
  OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR,
  "catchment_char_storage_mlr_diagnostics.csv"
)
catch_diag_tbl <- tibble()
if (file.exists(catch_diag_file)) {
  catch_diag <- read_csv(catch_diag_file, show_col_types = FALSE)
  catch_diag_tbl <- catch_diag %>%
    transmute(
      `MLR model set` = "Catchment characteristics",
      `Response Variable` = label_metric_abbrev(gsub("_mean$", "", as.character(Outcome))),
      `Residual Sample Size (n)` = suppressWarnings(as.numeric(n_residuals)),
      `Shapiro-Wilk p-value` = suppressWarnings(as.numeric(shapiro_p)),
      `Non-constant Variance Test p-value` = suppressWarnings(as.numeric(ncv_p)),
      `Shapiro-Wilk p-value` = ifelse(
        is.finite(shapiro_p),
        paste0(format(signif(shapiro_p, 3), scientific = TRUE, trim = TRUE), ifelse(shapiro_p < 0.05, "*", "")),
        NA_character_
      ),
      `Non-constant Variance Test p-value` = ifelse(
        is.finite(ncv_p),
        paste0(format(signif(ncv_p, 3), scientific = TRUE, trim = TRUE), ifelse(ncv_p < 0.05, "*", "")),
        NA_character_
      )
    )
}

eco_diag_file <- file.path(
  OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
  "storage_eco_response_mlr_diagnostics.csv"
)
eco_diag_tbl <- tibble()
if (file.exists(eco_diag_file)) {
  eco_diag <- read_csv(eco_diag_file, show_col_types = FALSE)
  eco_diag_tbl <- eco_diag %>%
    transmute(
      `MLR model set` = "Ecological responses",
      `Response Variable` = as.character(Response),
      `Residual Sample Size (n)` = suppressWarnings(as.numeric(n_residuals)),
      `Shapiro-Wilk p-value` = suppressWarnings(as.numeric(shapiro_p)),
      `Non-constant Variance Test p-value` = suppressWarnings(as.numeric(ncv_p)),
      `Shapiro-Wilk p-value` = ifelse(
        is.finite(shapiro_p),
        paste0(format(signif(shapiro_p, 3), scientific = TRUE, trim = TRUE), ifelse(shapiro_p < 0.05, "*", "")),
        NA_character_
      ),
      `Non-constant Variance Test p-value` = ifelse(
        is.finite(ncv_p),
        paste0(format(signif(ncv_p, 3), scientific = TRUE, trim = TRUE), ifelse(ncv_p < 0.05, "*", "")),
        NA_character_
      )
    )
}

diag_table <- bind_rows(catch_diag_tbl, eco_diag_tbl)
if (nrow(diag_table) > 0) {
  diag_table <- diag_table %>%
    arrange(`MLR model set`, `Response Variable`)
  write_csv(diag_table, file.path(supp_dir, "TableS8_mlr_model_diagnostics.csv"))
}

# ---- cleanup of legacy supplement artifacts not used in the final SI ----
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
