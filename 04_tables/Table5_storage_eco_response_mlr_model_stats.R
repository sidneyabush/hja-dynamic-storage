# write Table 5 from storage-eco response MLR outputs
# inputs: outputs/models/storage_eco_response_mlr/*.csv
# outputs: ms_materials/main/Table5_storage_eco_response_mlr_model_stats.csv

librarian::shelf(dplyr, readr, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")
source("helpers/table_utils.R")

eco_models <- load_eco_models()
ws_models <- load_ws_models()

eco_summary_table <- eco_models %>%
  transmute(
    `Eco-response variable` = factor(response_variable, levels = ECO_ORDER),
    `Selected Predictor(s)` = clean_predictor_labels(selected_predictors),
    `R2` = r2,
    `Adj R2` = adj_r2_marked,
    `RMSE` = rmse,
    `LOOCV RMSE` = rmse_loocv,
    `AICc` = aicc
  ) %>%
  arrange(`Eco-response variable`) %>%
  mutate(`Eco-response variable` = as.character(`Eco-response variable`))

all_models <- bind_rows(eco_models, ws_models) %>%
  mutate(
    model_group = factor(model_group, levels = c("Catchment characteristics", "Eco response pooled")),
    row_rank = case_when(
      model_group == "Catchment characteristics" ~ match(response_variable, WS_ORDER),
      model_group == "Eco response pooled" ~ match(response_variable, ECO_ORDER),
      TRUE ~ 999L
    )
  ) %>%
  arrange(model_group, row_rank) %>%
  transmute(
    `Model Group` = as.character(model_group),
    `Response Variable` = ifelse(
      model_group == "Catchment characteristics",
      label_metric_abbrev(response_variable),
      response_variable
    ),
    `Selected Predictor(s)` = clean_predictor_labels(selected_predictors),
    `R2` = r2,
    `Adj R2 marked` = adj_r2_marked,
    `RMSE` = rmse,
    `LOOCV RMSE` = rmse_loocv,
    `AICc` = aicc,
    `Model p` = model_p_text,
    `Model significant (alpha=0.05)` = model_significant_alpha_0_05
  )

dir.create(OUT_STATS_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(all_models, file.path(OUT_STATS_DIR, "mlr_all_models_listed.csv"))

dir.create(MS_TABLES_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(
  eco_summary_table,
  file.path(MS_TABLES_MAIN_DIR, "Table5_storage_eco_response_mlr_model_stats.csv")
)
