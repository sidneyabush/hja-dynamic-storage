# write supporting information Tables S7 to S12

# inputs:
# outputs/MTT_sensitivity/TableS7_MTT_sensitivity.csv
# outputs/models/catchment_char_storage_mlr/*.csv
# outputs/models/storage_eco_response_mlr/*.csv

# outputs:
# figs_tables_pub/supp/TableS7_MTT_sensitivity.csv
# figs_tables_pub/supp/TableS8_catchment_char_storage_mlr_model_stats.csv
# figs_tables_pub/supp/TableS9_catchment_alt_models_unique_deltaAICc_le2_BF.csv
# figs_tables_pub/supp/TableS10_mlr_model_diagnostics.csv
# figs_tables_pub/supp/TableS11_storage_eco_response_mlr_model_stats.csv
# figs_tables_pub/supp/TableS12_eco_alt_models_unique_deltaAICc_le2_BF.csv

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(dplyr, readr, tibble, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

# constants used for formatting model tables
ALPHA <- 0.05
ECO_ORDER <- c("Q7Q5", "T7DMax")
WS_ORDER <- STORAGE_METRIC_ORDER
WS_TABLE_ORDER <- STORAGE_METRIC_ORDER

# formatting functions used by Tables S8 to S12
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
  out <- gsub("\\bCHS\\b", "BF", out)
  out <- gsub("Pyroclastic \\(%\\)", "Pyroclastic", out)
  out <- gsub("Ash \\(%\\)", "Ash", out)
  out <- gsub("Lava 1 \\(%\\)", "Lava-1 (%)", out)
  out <- gsub("Lava 2 \\(%\\)", "Lava-2 (%)", out)
  out <- gsub("Total Landslide", "Landslide Total", out)
  out <- gsub("Young Landslide", "Landslide Young", out)
  out
}

compute_model_p_if_missing <- function(summary_df, id_col, n_col, r2_col, k_tbl) {
  if ("model_p_global" %in% names(summary_df)) {
    return(summary_df)
  }

  summary_df %>%
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
}

load_eco_models <- function() {
  # selected ecological response models become Table S11
  eco_dir <- OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR
  eco_res_file <- file.path(eco_dir, "storage_eco_response_mlr_results.csv")
  eco_sum_file <- file.path(eco_dir, "storage_eco_response_mlr_summary.csv")

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

  eco_sum %>%
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
}

load_ws_models <- function() {
  # selected catchment characteristic models become Table S8
  ws_dir <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
  ws_res_file <- file.path(ws_dir, "catchment_char_storage_mlr_results.csv")
  ws_sum_file <- file.path(ws_dir, "catchment_char_storage_mlr_summary.csv")

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

  ws_sum %>%
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
}

load_catchment_alt_models <- function() {
  # supported catchment characteristic models become Table S9
  catch_alt_file <- file.path(
    OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR,
    "catchment_char_storage_mlr_aicc_lt2.csv"
  )
  read_csv(catch_alt_file, show_col_types = FALSE) %>%
    transmute(
      `Response Variable` = gsub("_mean$", "", as.character(Outcome)),
      `Selected Predictor(s)` = clean_predictor_labels(as.character(Predictors_Final)),
      AICc = suppressWarnings(as.numeric(AICc)),
      delta_AICc = suppressWarnings(as.numeric(delta_AICc)),
      R2 = suppressWarnings(as.numeric(R2)),
      Adj_R2 = suppressWarnings(as.numeric(R2_adj)),
      RMSE = suppressWarnings(as.numeric(RMSE)),
      LOOCV_RMSE = suppressWarnings(as.numeric(RMSE_LOOCV)),
      n = suppressWarnings(as.integer(N))
    ) %>%
    group_by(`Response Variable`, `Selected Predictor(s)`) %>%
    summarise(
      AICc = first(AICc),
      delta_AICc = first(delta_AICc),
      R2 = first(R2),
      Adj_R2 = first(Adj_R2),
      RMSE = first(RMSE),
      LOOCV_RMSE = first(LOOCV_RMSE),
      n = first(n),
      .groups = "drop"
    ) %>%
    mutate(
      `Response Variable` = factor(
        `Response Variable`,
        levels = c("RBI", "RCS", "FDC", "SD", "WB", "BF", "DR", "Fyw", "MTT")
      )
    ) %>%
    arrange(`Response Variable`, delta_AICc, `Selected Predictor(s)`) %>%
    mutate(`Response Variable` = as.character(`Response Variable`)) %>%
    transmute(
      `Response Variable`,
      `Selected Predictor(s)`,
      `R2` = R2,
      `Adj R2` = Adj_R2,
      `RMSE` = RMSE,
      `LOOCV RMSE` = LOOCV_RMSE,
      `AICc` = AICc
    )
}

load_eco_alt_models <- function() {
  # supported ecological response models become Table S12
  eco_alt_file <- file.path(
    OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
    "storage_eco_response_mlr_aicc_lt2.csv"
  )
  read_csv(eco_alt_file, show_col_types = FALSE) %>%
    transmute(
      `Response Variable` = as.character(Response),
      `Selected Predictor(s)` = clean_predictor_labels(as.character(Predictors_Final)),
      AICc = suppressWarnings(as.numeric(AICc)),
      delta_AICc = suppressWarnings(as.numeric(delta_AICc)),
      R2 = suppressWarnings(as.numeric(R2)),
      Adj_R2 = suppressWarnings(as.numeric(R2_adj)),
      RMSE = suppressWarnings(as.numeric(RMSE)),
      LOOCV_RMSE = suppressWarnings(as.numeric(RMSE_LOOCV)),
      n = suppressWarnings(as.integer(n))
    ) %>%
    group_by(`Response Variable`, `Selected Predictor(s)`) %>%
    summarise(
      AICc = first(AICc),
      delta_AICc = first(delta_AICc),
      R2 = first(R2),
      Adj_R2 = first(Adj_R2),
      RMSE = first(RMSE),
      LOOCV_RMSE = first(LOOCV_RMSE),
      n = first(n),
      .groups = "drop"
    ) %>%
    mutate(`Response Variable` = factor(`Response Variable`, levels = ECO_ORDER)) %>%
    arrange(`Response Variable`, delta_AICc, `Selected Predictor(s)`) %>%
    mutate(`Response Variable` = as.character(`Response Variable`)) %>%
    transmute(
      `Response Variable`,
      `Selected Predictor(s)`,
      `R2` = R2,
      `Adj R2` = Adj_R2,
      `RMSE` = RMSE,
      `LOOCV RMSE` = LOOCV_RMSE,
      `AICc` = AICc
    )
}

dir.create(MS_TABLES_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

write_table_s7 <- function() {
  # Table S7 is created by the MTT sensitivity script
  source_file <- file.path(OUTPUT_DIR, "MTT_sensitivity", "TableS7_MTT_sensitivity.csv")
  table_s7 <- read_csv(source_file, show_col_types = FALSE)
  write_csv(table_s7, file.path(MS_TABLES_SUPP_DIR, "TableS7_MTT_sensitivity.csv"))
}

write_table_s8 <- function() {
  # summarize the selected catchment characteristic models
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
    file.path(MS_TABLES_SUPP_DIR, "TableS8_catchment_char_storage_mlr_model_stats.csv")
  )
}

write_table_s9 <- function() {
  # list alternative catchment characteristic models with delta AICc <= 2
  catch_alt <- load_catchment_alt_models()

  write_csv(
    catch_alt,
    file.path(MS_TABLES_SUPP_DIR, "TableS9_catchment_alt_models_unique_deltaAICc_le2_BF.csv")
  )
}

write_table_s11 <- function() {
  # summarize the selected ecological response models
  eco_models <- load_eco_models()

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

  write_csv(
    eco_summary_table,
    file.path(MS_TABLES_SUPP_DIR, "TableS11_storage_eco_response_mlr_model_stats.csv")
  )
}

write_table_s12 <- function() {
  # list alternative ecological response models with delta AICc <= 2
  eco_alt <- load_eco_alt_models()

  write_csv(
    eco_alt,
    file.path(MS_TABLES_SUPP_DIR, "TableS12_eco_alt_models_unique_deltaAICc_le2_BF.csv")
  )
}

write_table_s10 <- function() {
  # combine residual diagnostics for both model sets
  catch_diag_file <- file.path(
    OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR,
    "catchment_char_storage_mlr_diagnostics.csv"
  )
  catch_diag <- read_csv(catch_diag_file, show_col_types = FALSE)
  catch_diag_tbl <- catch_diag %>%
    transmute(
      `MLR model set` = "Catchment characteristics",
      `Response Variable` = label_metric_abbrev(gsub("_mean$", "", as.character(Outcome))),
      `Residual Sample Size (n)` = suppressWarnings(as.numeric(n_residuals)),
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

  eco_diag_file <- file.path(
    OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
    "storage_eco_response_mlr_diagnostics.csv"
  )
  eco_diag <- read_csv(eco_diag_file, show_col_types = FALSE)
  eco_diag_tbl <- eco_diag %>%
    transmute(
      `MLR model set` = "Ecological responses",
      `Response Variable` = as.character(Response),
      `Residual Sample Size (n)` = suppressWarnings(as.numeric(n_residuals)),
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

  diag_table <- bind_rows(catch_diag_tbl, eco_diag_tbl) %>%
    arrange(`MLR model set`, `Response Variable`)

  write_csv(diag_table, file.path(MS_TABLES_SUPP_DIR, "TableS10_mlr_model_diagnostics.csv"))
}

write_table_s7()
write_table_s8()
write_table_s9()
write_table_s10()
write_table_s11()
write_table_s12()
