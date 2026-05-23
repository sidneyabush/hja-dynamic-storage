librarian::shelf(dplyr, readr, cran_repo = "https://cloud.r-project.org")

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

load_eco_models <- function() {
  eco_dir <- OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR
  eco_res_file <- file.path(eco_dir, "storage_eco_response_mlr_results.csv")
  eco_sum_file <- file.path(eco_dir, "storage_eco_response_mlr_summary.csv")

  if (!file.exists(eco_res_file) || !file.exists(eco_sum_file)) {
    stop(
      "Missing eco model inputs: expected storage_eco_response_mlr_results.csv and storage_eco_response_mlr_summary.csv in ",
      eco_dir
    )
  }

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
  ws_dir <- OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR
  ws_res_file <- file.path(ws_dir, "catchment_char_storage_mlr_results.csv")
  ws_sum_file <- file.path(ws_dir, "catchment_char_storage_mlr_summary.csv")

  if (!file.exists(ws_res_file) || !file.exists(ws_sum_file)) {
    stop(
      "Missing catchment model inputs: expected catchment_char_storage_mlr_results.csv and catchment_char_storage_mlr_summary.csv in ",
      ws_dir
    )
  }

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
  catch_alt_file <- file.path(
    OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR,
    "catchment_char_storage_mlr_aicc_lt2.csv"
  )
  if (!file.exists(catch_alt_file)) {
    return(NULL)
  }

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
    mutate(`Response Variable` = factor(`Response Variable`, levels = c("RBI", "RCS", "FDC", "SD", "WB", "BF", "DR", "Fyw", "MTT"))) %>%
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
  eco_alt_file <- file.path(
    OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR,
    "storage_eco_response_mlr_aicc_lt2.csv"
  )
  if (!file.exists(eco_alt_file)) {
    return(NULL)
  }

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
