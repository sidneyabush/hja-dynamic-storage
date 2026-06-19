# write the catchment alternative-model table for the supplement
# inputs: outputs/models/catchment_char_storage_mlr/catchment_char_storage_mlr_aicc_lt2.csv
# outputs: ms_materials/supp/TableS6_catchment_alt_models_unique_deltaAICc_le2_BF.csv

librarian::shelf(readr, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")
source("helpers/table_utils.R")

catch_alt <- load_catchment_alt_models()
if (is.null(catch_alt)) {
  stop("Missing required catchment alternative-model input file")
}

dir.create(MS_TABLES_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(
  catch_alt,
  file.path(MS_TABLES_SUPP_DIR, "TableS6_catchment_alt_models_unique_deltaAICc_le2_BF.csv")
)
