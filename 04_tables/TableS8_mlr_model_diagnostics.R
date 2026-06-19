# write model diagnostics for the supplement
# inputs: catchment and ecological-response MLR diagnostics files
# outputs: ms_materials/supp/TableS8_mlr_model_diagnostics.csv

librarian::shelf(dplyr, readr, cran_repo = "https://cloud.r-project.org")

rm(list = ls())
source("config.R")

supp_dir <- MS_TABLES_SUPP_DIR
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

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
        paste0(
          format(signif(shapiro_p, 3), scientific = TRUE, trim = TRUE),
          ifelse(shapiro_p < 0.05, "*", "")
        ),
        NA_character_
      ),
      `Non-constant Variance Test p-value` = ifelse(
        is.finite(ncv_p),
        paste0(
          format(signif(ncv_p, 3), scientific = TRUE, trim = TRUE),
          ifelse(ncv_p < 0.05, "*", "")
        ),
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
        paste0(
          format(signif(shapiro_p, 3), scientific = TRUE, trim = TRUE),
          ifelse(shapiro_p < 0.05, "*", "")
        ),
        NA_character_
      ),
      `Non-constant Variance Test p-value` = ifelse(
        is.finite(ncv_p),
        paste0(
          format(signif(ncv_p, 3), scientific = TRUE, trim = TRUE),
          ifelse(ncv_p < 0.05, "*", "")
        ),
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
