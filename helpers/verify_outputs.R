# check that the expected output files were written
# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(readr, dplyr, cran_repo = "https://cloud.r-project.org")

find_repo_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(10)) {
    has_config <- file.exists(file.path(cur, "config.R"))
    has_metrics <- dir.exists(file.path(cur, "01_storage_calcs"))
    if (has_config && has_metrics) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(start_dir, winslash = "/", mustWork = FALSE)
}

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
if (is.na(script_path) || script_path == "") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
}
start_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_root <- find_repo_root(start_dir)

source(file.path(repo_root, "config.R"))
master_dir <- file.path(OUTPUT_DIR, "master")

required_outputs <- c(
  file.path(master_dir, MASTER_ANNUAL_FILE),
  file.path(master_dir, MASTER_SITE_FILE),
  file.path(OUT_STATS_ANOVA_DIR, "anova_results.csv"),
  file.path(OUT_STATS_ANOVA_DIR, "tukey_hsd_results.csv"),
  file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv"),
  file.path(OUT_STATS_ANOVA_DIR, "storage_metrics_summary_stats_by_site.csv"),
  file.path(OUT_STATS_PCA_DIR, "pca_loadings.csv"),
  file.path(OUT_STATS_PCA_DIR, "pca_variance_explained.csv"),
  file.path(OUT_MODELS_CATCHMENT_CHAR_STORAGE_MLR_DIR, "catchment_char_storage_mlr_summary.csv"),
  file.path(OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR, "storage_eco_response_mlr_summary.csv"),
  file.path(MS_FIG_MAIN_DIR, "Fig2_ds_annual_boxplots.png"),
  file.path(MS_FIG_MAIN_PDF_DIR, "Fig2_ds_annual_boxplots.pdf"),
  file.path(MS_FIG_MAIN_DIR, "Fig3_ds_pca_annual.png"),
  file.path(MS_FIG_MAIN_PDF_DIR, "Fig3_ds_pca_annual.pdf"),
  file.path(MS_FIG_MAIN_DIR, "Fig4_mobile_storage.png"),
  file.path(MS_FIG_MAIN_PDF_DIR, "Fig4_mobile_storage.pdf"),
  file.path(MS_FIG_MAIN_DIR, "Fig5_dynamic_mobile_corr.png"),
  file.path(MS_FIG_MAIN_PDF_DIR, "Fig5_dynamic_mobile_corr.pdf"),
  file.path(MS_FIG_MAIN_DIR, "Fig6_catch_char_mlr_beta.png"),
  file.path(MS_FIG_MAIN_PDF_DIR, "Fig6_catch_char_mlr_beta.pdf"),
  file.path(MS_FIG_MAIN_DIR, "Fig7_eco_mlr_beta.png"),
  file.path(MS_FIG_MAIN_PDF_DIR, "Fig7_eco_mlr_beta.pdf"),
  file.path(MS_FIG_MAIN_DIR, "Fig8_eco_observed_v_predicted.png"),
  file.path(MS_FIG_MAIN_PDF_DIR, "Fig8_eco_observed_v_predicted.pdf"),
  file.path(MS_FIG_MAIN_DIR, "Fig9_conceptual_diagram.png"),
  file.path(MS_FIG_MAIN_PDF_DIR, "Fig9_conceptual_diagram.pdf"),
  file.path(MS_TABLES_MAIN_DIR, "Table4_catchment_char_storage_mlr_model_stats.csv"),
  file.path(MS_TABLES_MAIN_DIR, "Table5_storage_eco_response_mlr_model_stats.csv"),
  file.path(MS_FIG_SUPP_DIR, "FigS1_met_context.png"),
  file.path(MS_FIG_SUPP_PDF_DIR, "FigS1_met_context.pdf"),
  file.path(MS_TABLES_SUPP_DIR, "TableS5_MTT_sensitivity.csv"),
  file.path(MS_TABLES_SUPP_DIR, "TableS6_catchment_alt_models_unique_deltaAICc_le2_BF.csv"),
  file.path(MS_TABLES_SUPP_DIR, "TableS7_eco_alt_models_unique_deltaAICc_le2_BF.csv"),
  file.path(MS_TABLES_SUPP_DIR, "TableS8_mlr_model_diagnostics.csv")
)

missing_outputs <- required_outputs[!file.exists(required_outputs)]
if (isTRUE(WRITE_AUX_OUTPUTS)) {
  aux_required <- c(
    file.path(master_dir, "master_site_metric_summary_stats.csv"),
    file.path(OUT_MET_SUPPORT_DIR, "site_metric_availability.csv")
  )
  missing_outputs <- c(
    missing_outputs,
    aux_required[!file.exists(aux_required)]
  )
}
if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  table_required <- c(file.path(OUT_TABLES_MLR_DIR, "mlr_main_results_table.csv"))
  missing_outputs <- c(
    missing_outputs,
    table_required[!file.exists(table_required)]
  )
}

if (length(missing_outputs) > 0) {
  msg <- paste0(
    "Missing required output file(s):\n- ",
    paste(missing_outputs, collapse = "\n- ")
  )
  stop(msg)
}

annual <- read_csv(file.path(master_dir, MASTER_ANNUAL_FILE), show_col_types = FALSE)
if (!all(c("site", "year") %in% names(annual))) {
  stop("master_annual is missing required columns: site/year")
}
if (nrow(annual) == 0) stop("master_annual has zero rows")
if (any(!is.finite(annual$year))) stop("master_annual has non-finite years")

# site checks
annual_sites <- unique(as.character(annual$site))
unknown_sites <- setdiff(annual_sites, SITE_ORDER_ALL)
if (length(unknown_sites) > 0) {
  stop(paste0("Unknown site code(s) in master_annual: ", paste(unknown_sites, collapse = ", ")))
}

# keep the configured site order in the Tukey group letters file
letters_path <- file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv")
letters_df <- read_csv(letters_path, show_col_types = FALSE)
if (!all(c("metric", "site") %in% names(letters_df))) {
  stop("tukey_group_letters is missing required columns: metric/site")
}

for (m in unique(letters_df$metric)) {
  s <- letters_df %>% filter(metric == m) %>% pull(site) %>% as.character()
  allowed <- SITE_ORDER_HYDROMETRIC[SITE_ORDER_HYDROMETRIC %in% s]
  if (!identical(s, allowed)) {
    stop(paste0("Site order mismatch in tukey_group_letters for metric: ", m))
  }
}
