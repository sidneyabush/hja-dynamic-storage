#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

find_repo_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(10)) {
    has_config <- file.exists(file.path(cur, "config.R"))
    has_metrics <- dir.exists(file.path(cur, "01_metrics"))
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

required_outputs <- c(
  file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE),
  file.path(OUT_MASTER_DIR, MASTER_SITE_FILE),
  file.path(OUT_STATS_ANOVA_DIR, "anova_results.csv"),
  file.path(OUT_STATS_ANOVA_DIR, "tukey_hsd_results.csv"),
  file.path(OUT_STATS_ANOVA_DIR, "tukey_group_letters.csv"),
  file.path(OUT_STATS_PCA_DIR, "pca_loadings.csv"),
  file.path(OUT_STATS_PCA_DIR, "pca_variance_explained.csv"),
  file.path(OUT_STATS_MLR_CATCH_CHARS_DIR, "catch_chars_storage_mlr_summary_strict.csv"),
  file.path(OUT_STATS_MLR_CATCH_CHARS_DIR, "catch_chars_storage_mlr_summary_non_strict.csv"),
  file.path(OUT_STATS_MLR_ECO_DIR, "storage_eco_mlr_summary_strict.csv"),
  file.path(OUT_STATS_MLR_ECO_DIR, "storage_eco_mlr_summary_non_strict.csv"),
  file.path(FIGURES_DIR, "main", "ds_summary.png"),
  file.path(FIGURES_DIR, "main", "pca_biplot.png")
)

missing_outputs <- required_outputs[!file.exists(required_outputs)]
if (length(missing_outputs) > 0) {
  stop(
    paste0(
      "Missing required output file(s):\n- ",
      paste(missing_outputs, collapse = "\n- ")
    )
  )
}

annual <- read_csv(file.path(OUT_MASTER_DIR, MASTER_ANNUAL_FILE), show_col_types = FALSE)
if (!all(c("site", "year") %in% names(annual))) {
  stop("master_annual is missing required columns: site/year")
}
if (nrow(annual) == 0) stop("master_annual has zero rows")
if (any(!is.finite(annual$year))) stop("master_annual has non-finite years")

# Site checks.
annual_sites <- unique(as.character(annual$site))
unknown_sites <- setdiff(annual_sites, SITE_ORDER_ALL)
if (length(unknown_sites) > 0) {
  stop(paste0("Unknown site code(s) in master_annual: ", paste(unknown_sites, collapse = ", ")))
}

# Enforce configured site order in Tukey group letters.
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

message("Output verification passed.")
