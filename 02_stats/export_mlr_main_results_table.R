# Export a unified main MLR results table for manuscript reporting.
# Inputs: mlr_dir/storage_ecovar_mlr_model_perf.csv; mlr_dir/watershed_char_storage_mlr_model_perf.csv; val_dir/watershed_char_storage_mlr_loocv_validation.csv.
# Author: Sidney Bush
# Date: 2026-02-13

library(dplyr)
library(readr)

rm(list = ls())

# Load project config
source("config.R")


mlr_dir <- file.path(OUT_TABLES_DIR, "mlr")
val_dir <- file.path(OUT_TABLES_DIR, "validation")
if (!dir.exists(mlr_dir)) dir.create(mlr_dir, recursive = TRUE, showWarnings = FALSE)

eco_perf_file <- file.path(mlr_dir, "storage_ecovar_mlr_model_perf.csv")
ws_perf_file <- file.path(mlr_dir, "watershed_char_storage_mlr_model_perf.csv")
ws_val_file <- file.path(val_dir, "watershed_char_storage_mlr_loocv_validation.csv")

if (!file.exists(eco_perf_file)) stop("Missing file: ", eco_perf_file)
if (!file.exists(ws_perf_file)) stop("Missing file: ", ws_perf_file)

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)][1]
  if (is.na(hit)) return(rep(NA, nrow(df)))
  df[[hit]]
}

eco_raw <- read_csv(eco_perf_file, show_col_types = FALSE)
eco_perf <- read_csv(eco_perf_file, show_col_types = FALSE) %>%
  transmute(
    model_family = "storage_ecovar_mlr",
    site = pick_col(eco_raw, c("site", "Site")),
    response = pick_col(eco_raw, c("response", "Response")),
    n = pick_col(eco_raw, c("n", "N")),
    predictors_final = pick_col(eco_raw, c("predictors_final", "Predictors_Final")),
    r2 = pick_col(eco_raw, c("r2", "R2")),
    r2_adj = pick_col(eco_raw, c("r2_adj", "R2_adj")),
    rmse = pick_col(eco_raw, c("rmse", "RMSE")),
    rmse_loocv = pick_col(eco_raw, c("rmse_loocv", "RMSE_LOOCV")),
    aicc = pick_col(eco_raw, c("aicc", "AICc"))
  )

ws_perf <- read_csv(ws_perf_file, show_col_types = FALSE)

ws_r2 <- if (file.exists(ws_val_file)) {
  read_csv(ws_val_file, show_col_types = FALSE) %>%
    transmute(
      response = outcome,
      r2 = r2_model,
      rmse_loocv = rmse_loocv
    )
} else {
  tibble(response = character(), r2 = numeric(), rmse_loocv = numeric())
}

ws_main <- ws_perf %>%
  transmute(
    model_family = "watershed_char_storage_mlr",
    site = "all_sites",
    response = pick_col(ws_perf, c("Outcome", "response")),
    n = pick_col(ws_perf, c("N", "n")),
    predictors_final = pick_col(ws_perf, c("Predictors_Final", "predictors_final")),
    r2_adj = pick_col(ws_perf, c("R2_adj", "r2_adj")),
    rmse = pick_col(ws_perf, c("RMSE", "rmse")),
    aicc = pick_col(ws_perf, c("AICc", "aicc"))
  ) %>%
  left_join(ws_r2, by = "response") %>%
  dplyr::select(model_family, site, response, n, predictors_final, r2, r2_adj, rmse, rmse_loocv, aicc)

site_rank <- setNames(seq_along(SITE_ORDER_HYDROMETRIC), SITE_ORDER_HYDROMETRIC)
response_order_eco <- c("Q_7Q5", "T_7DMax", "T_Q7Q5")

main_table <- bind_rows(eco_perf, ws_main) %>%
  mutate(
    model_rank = if_else(model_family == "storage_ecovar_mlr", 1L, 2L),
    site_rank = if_else(site %in% names(site_rank), as.integer(site_rank[site]), 999L),
    response_rank = if_else(response %in% response_order_eco, match(response, response_order_eco), 999L)
  ) %>%
  arrange(model_rank, site_rank, response_rank, desc(r2_adj)) %>%
  dplyr::select(-model_rank, -site_rank, -response_rank) %>%
  mutate(
    response = gsub("_", " ", response),
    predictors_final = gsub("_", " ", predictors_final)
  ) %>%
  mutate(
    r2 = signif(r2, 3),
    r2_adj = signif(r2_adj, 3),
    rmse = signif(rmse, 3),
    rmse_loocv = signif(rmse_loocv, 3),
    aicc = signif(aicc, 3)
  )

out_file <- file.path(mlr_dir, "mlr_main_results_table.csv")
write_csv(main_table, out_file)
