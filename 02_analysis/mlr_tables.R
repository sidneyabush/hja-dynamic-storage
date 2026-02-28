# export model-stats tables for pooled eco-response and catchment mlr models.
# inputs:
# - out_models_storage_eco_response_mlr_dir/storage_eco_response_mlr_results.csv
# - out_models_storage_eco_response_mlr_dir/storage_eco_response_mlr_summary.csv
# - out_models_catchment_char_storage_mlr_dir/catchment_char_storage_mlr_results.csv
# - out_models_catchment_char_storage_mlr_dir/catchment_char_storage_mlr_summary.csv
# author: sidney bush
# date: 2026-02-24

library(dplyr)
library(readr)

rm(list = ls())

source("config.R")

ALPHA <- 0.05
ECO_ORDER <- c("Q7Q5", "T7DMax")
WS_ORDER <- STORAGE_METRIC_ORDER
WS_TABLE_ORDER <- STORAGE_METRIC_ORDER

format_p <- function(x) {
  ifelse(
    is.finite(x),
    format(signif(x, 3), scientific = TRUE, trim = TRUE),
    NA_character_
  )
}

format_adj_r2_marked <- function(adj_r2, model_sig) {
  base <- ifelse(is.finite(adj_r2), sprintf("%.3f", adj_r2), NA_character_)
  sig <- ifelse(is.na(model_sig), FALSE, as.logical(model_sig))
  ifelse(sig, paste0(base, "*"), base)
}

clean_predictor_labels <- function(x) {
  out <- as.character(x)
  out <- gsub("Pyroclastic \\(%\\)", "Pyroclastic", out)
  out <- gsub("Ash \\(%\\)", "Ash", out)
  out <- gsub("Lava 1 \\(%\\)", "Lava-1 (%)", out)
  out <- gsub("Lava 2 \\(%\\)", "Lava-2 (%)", out)
  out <- gsub("Total Landslide", "Landslide Total", out)
  out <- gsub("Young Landslide", "Landslide Young", out)
  out
}

compute_model_p_if_missing <- function(summary_df, id_col, n_col, r2_col, k_tbl) {
  if ("model_p_global" %in% names(summary_df)) return(summary_df)

  out <- summary_df %>%
    left_join(k_tbl, by = setNames("model_id", id_col)) %>%
    mutate(
      r2_num = suppressWarnings(as.numeric(.data[[r2_col]])),
      n_num = suppressWarnings(as.numeric(.data[[n_col]])),
      k_num = suppressWarnings(as.numeric(k_predictors)),
      f_stat = ifelse(
        is.finite(r2_num) & is.finite(n_num) & is.finite(k_num) &
          k_num > 0 & (n_num - k_num - 1) > 0 & r2_num < 1,
        (r2_num / k_num) / ((1 - r2_num) / (n_num - k_num - 1)),
        NA_real_
      ),
      model_p_global = ifelse(
        is.finite(f_stat),
        pf(f_stat, k_num, n_num - k_num - 1, lower.tail = FALSE),
        NA_real_
      )
    ) %>%
    select(-r2_num, -n_num, -k_num, -f_stat, -k_predictors)

  out
}

to_md_table <- function(df_in) {
  df <- as.data.frame(df_in, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    return(c("| (no rows) |", "| --- |"))
  }
  df[] <- lapply(df, function(x) {
    out <- as.character(x)
    out[is.na(out)] <- "-"
    out
  })
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(header, sep, rows)
}

# ---- inputs ----
eco_dir <- OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR
legacy_eco_dir <- file.path(OUT_STATS_DIR, "storage_ecovar_mlr")
eco_res_candidates <- c(
  file.path(legacy_eco_dir, "storage_ecovar_mlr_all_sites_results.csv"),
  file.path(eco_dir, "storage_eco_response_mlr_results.csv")
)
eco_sum_candidates <- c(
  file.path(legacy_eco_dir, "storage_ecovar_mlr_all_sites_summary.csv"),
  file.path(eco_dir, "storage_eco_response_mlr_summary.csv")
)
eco_res_file <- eco_res_candidates[file.exists(eco_res_candidates)][1]
eco_sum_file <- eco_sum_candidates[file.exists(eco_sum_candidates)][1]

ws_dir <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
ws_res_file <- file.path(ws_dir, "catchment_char_storage_mlr_results.csv")
ws_sum_file <- file.path(ws_dir, "catchment_char_storage_mlr_summary.csv")

if (is.na(eco_res_file) || is.na(eco_sum_file) || !file.exists(eco_res_file) || !file.exists(eco_sum_file)) {
  stop(
    "Missing eco model inputs; checked results: ",
    paste(eco_res_candidates, collapse = ", "),
    " ; checked summaries: ",
    paste(eco_sum_candidates, collapse = ", ")
  )
}
if (!file.exists(ws_res_file) || !file.exists(ws_sum_file)) {
  stop("Missing catchment model inputs: expected catchment_char_storage_mlr_results.csv and _summary.csv")
}

# ---- pooled eco-response models ----
eco_res <- read_csv(eco_res_file, show_col_types = FALSE) %>%
  transmute(
    model_id = gsub("_", "", as.character(Response), fixed = TRUE),
    Beta_Std = suppressWarnings(as.numeric(Beta_Std))
  )
eco_k <- eco_res %>%
  group_by(model_id) %>%
  summarise(k_predictors = sum(is.finite(Beta_Std)), .groups = "drop")

eco_sum <- read_csv(eco_sum_file, show_col_types = FALSE) %>%
  mutate(model_id = gsub("_", "", as.character(Response), fixed = TRUE))
eco_sum <- compute_model_p_if_missing(
  summary_df = eco_sum,
  id_col = "model_id",
  n_col = "n",
  r2_col = "R2",
  k_tbl = eco_k
)

eco_models <- eco_sum %>%
  transmute(
    model_group = "Eco response pooled",
    response_variable = as.character(model_id),
    selected_predictors = as.character(Predictors_Final),
    n = suppressWarnings(as.integer(n)),
    r2 = suppressWarnings(as.numeric(R2)),
    adj_r2 = suppressWarnings(as.numeric(R2_adj)),
    rmse = suppressWarnings(as.numeric(RMSE)),
    rmse_loocv = suppressWarnings(as.numeric(RMSE_LOOCV)),
    aicc = suppressWarnings(as.numeric(AICc)),
    model_p = suppressWarnings(as.numeric(model_p_global)),
    model_significant_alpha_0_05 = ifelse(is.finite(model_p), model_p <= ALPHA, NA)
  ) %>%
  mutate(
    response_variable = factor(response_variable, levels = ECO_ORDER),
    r2 = signif(r2, 3),
    adj_r2 = signif(adj_r2, 3),
    rmse = signif(rmse, 3),
    rmse_loocv = signif(rmse_loocv, 3),
    aicc = signif(aicc, 3),
    model_p_text = format_p(model_p),
    adj_r2_marked = format_adj_r2_marked(adj_r2, model_significant_alpha_0_05)
  ) %>%
  arrange(response_variable) %>%
  mutate(response_variable = as.character(response_variable))

# ---- catchment models ----
ws_res <- read_csv(ws_res_file, show_col_types = FALSE) %>%
    transmute(
      model_id = gsub("_mean$", "", as.character(Outcome)),
      Beta_Std = suppressWarnings(as.numeric(Beta_Std))
    )
ws_k <- ws_res %>%
  group_by(model_id) %>%
  summarise(k_predictors = sum(is.finite(Beta_Std)), .groups = "drop")

ws_sum <- read_csv(ws_sum_file, show_col_types = FALSE) %>%
  mutate(model_id = gsub("_mean$", "", as.character(Outcome)))
ws_sum <- compute_model_p_if_missing(
  summary_df = ws_sum,
  id_col = "model_id",
  n_col = "N",
  r2_col = "R2",
  k_tbl = ws_k
)

ws_models <- ws_sum %>%
  transmute(
    model_group = "Catchment characteristics",
    response_variable = as.character(model_id),
    selected_predictors = as.character(Predictors_Final),
    n = suppressWarnings(as.integer(N)),
    r2 = suppressWarnings(as.numeric(R2)),
    adj_r2 = suppressWarnings(as.numeric(R2_adj)),
    rmse = suppressWarnings(as.numeric(RMSE)),
    rmse_loocv = suppressWarnings(as.numeric(RMSE_LOOCV)),
    aicc = suppressWarnings(as.numeric(AICc)),
    model_p = suppressWarnings(as.numeric(model_p_global)),
    model_significant_alpha_0_05 = ifelse(is.finite(model_p), model_p <= ALPHA, NA)
  ) %>%
  mutate(
    response_variable = factor(response_variable, levels = WS_ORDER),
    r2 = signif(r2, 3),
    adj_r2 = signif(adj_r2, 3),
    rmse = signif(rmse, 3),
    rmse_loocv = signif(rmse_loocv, 3),
    aicc = signif(aicc, 3),
    model_p_text = format_p(model_p),
    adj_r2_marked = format_adj_r2_marked(adj_r2, model_significant_alpha_0_05)
  ) %>%
  arrange(response_variable) %>%
  mutate(response_variable = as.character(response_variable))

# ---- main-text catchment summary table ----
ws_summary_table <- ws_models %>%
  transmute(
    `Response Variable` = factor(response_variable, levels = WS_TABLE_ORDER),
    `Selected Predictor(s)` = clean_predictor_labels(selected_predictors),
    `R2` = r2,
    `Adj R2` = adj_r2_marked,
    `RMSE` = rmse,
    `LOOCV RMSE` = rmse_loocv,
    `AICc` = aicc
  ) %>%
  arrange(`Response Variable`) %>%
  mutate(`Response Variable` = as.character(`Response Variable`))

write_csv(ws_summary_table, file.path(ws_dir, "catchment_char_storage_mlr_model_stats_table.csv"))
dir.create(MS_TABLES_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(
  ws_summary_table,
  file.path(MS_TABLES_MAIN_DIR, "table4_catchment_char_storage_mlr_model_stats.csv")
)

# ---- single all-models csv (eco + catchment) ----
dir.create(OUT_STATS_DIR, recursive = TRUE, showWarnings = FALSE)

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
    `Response Variable` = response_variable,
    `Selected Predictor(s)` = clean_predictor_labels(selected_predictors),
    `R2` = r2,
    `Adj R2 marked` = adj_r2_marked,
    `RMSE` = rmse,
    `LOOCV RMSE` = rmse_loocv,
    `AICc` = aicc,
    `Model p` = model_p_text,
    `Model significant (alpha=0.05)` = model_significant_alpha_0_05
  )

write_csv(all_models, file.path(OUT_STATS_DIR, "mlr_all_models_listed.csv"))
dir.create(MS_TABLES_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(
  all_models,
  file.path(MS_TABLES_MAIN_DIR, "table5_mlr_all_models_listed.csv")
)

old_supp_all_models <- file.path(MS_TABLES_SUPP_DIR, "tableSX_mlr_all_models_listed.csv")
if (file.exists(old_supp_all_models)) file.remove(old_supp_all_models)
