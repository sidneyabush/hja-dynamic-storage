# write the ecological-response alternative-model table for the supplement
# inputs: outputs/models/storage_eco_response_mlr/storage_eco_response_mlr_aicc_lt2.csv
# outputs: ms_materials/supp/TableS7_eco_alt_models_unique_deltaAICc_le2_BF.csv

librarian::shelf(readr, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")
source("helpers/table_utils.R")

eco_alt <- load_eco_alt_models()
if (is.null(eco_alt)) {
  stop("Missing required ecological-response alternative-model input file")
}

dir.create(MS_TABLES_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(
  eco_alt,
  file.path(MS_TABLES_SUPP_DIR, "TableS7_eco_alt_models_unique_deltaAICc_le2_BF.csv")
)
