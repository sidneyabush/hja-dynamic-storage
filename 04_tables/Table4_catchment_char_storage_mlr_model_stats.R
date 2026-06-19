# write Table 4 from catchment-characteristic storage MLR outputs
# inputs: outputs/models/catchment_char_storage_mlr/*.csv
# outputs: ms_materials/main/Table4_catchment_char_storage_mlr_model_stats.csv

librarian::shelf(dplyr, readr, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")
source("helpers/table_utils.R")

ws_models <- load_ws_models()

ws_summary_table <- ws_models %>%
  transmute(
    `Storage metric` = factor(
      label_metric_abbrev(response_variable),
      levels = label_metric_abbrev(WS_TABLE_ORDER)
    ),
    `Selected Predictor(s)` = clean_predictor_labels(selected_predictors),
    `R2` = r2,
    `Adj R2` = adj_r2_marked,
    `RMSE` = rmse,
    `LOOCV RMSE` = rmse_loocv,
    `AICc` = aicc
  ) %>%
  arrange(`Storage metric`) %>%
  mutate(`Storage metric` = as.character(`Storage metric`))

write_csv(
  ws_summary_table,
  file.path(
    OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR,
    "catchment_char_storage_mlr_model_stats_table.csv"
  )
)

dir.create(MS_TABLES_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(
  ws_summary_table,
  file.path(MS_TABLES_MAIN_DIR, "Table4_catchment_char_storage_mlr_model_stats.csv")
)
