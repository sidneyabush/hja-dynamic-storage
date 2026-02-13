# -----------------------------------------------------------------------------
# Output Manifest and Coverage
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)

rm(list = ls())

script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    getwd()
  }
})
if (is.null(script_dir) || script_dir == "" || script_dir == ".") script_dir <- getwd()

config_path <- file.path(script_dir, "config.R")
if (!file.exists(config_path)) config_path <- file.path(dirname(script_dir), "config.R")
if (!file.exists(config_path)) config_path <- file.path(getwd(), "config.R")
source(config_path)

manifest_dir <- OUT_TABLES_DIR
if (!dir.exists(manifest_dir)) dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

collect_files <- function(root, kind) {
  if (!dir.exists(root)) return(tibble())
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  tibble(
    artifact_kind = kind,
    file_path = files,
    file_name = basename(files),
    rel_path = sub(paste0("^", normalizePath(OUTPUT_DIR, winslash = "/"), "/?"), "", normalizePath(files, winslash = "/")),
    bytes = as.numeric(file.info(files)$size)
  )
}

manifest <- bind_rows(
  collect_files(file.path(OUTPUT_DIR, "figs"), "figure"),
  collect_files(file.path(OUTPUT_DIR, "tables"), "table"),
  collect_files(file.path(OUTPUT_DIR, "models"), "model")
) %>%
  mutate(
    section = case_when(
      grepl("^figs/main/", rel_path) ~ "main",
      grepl("^figs/supp/", rel_path) ~ "supp",
      grepl("^tables/", rel_path) ~ "tables",
      grepl("^models/", rel_path) ~ "models",
      TRUE ~ "other"
    ),
    path_group = case_when(
      grepl("anova_tukey", rel_path) ~ "anova_tukey",
      grepl("pca", rel_path) ~ "pca",
      grepl("watershed_char_storage_mlr", rel_path) ~ "watershed_char_storage_mlr",
      grepl("storage_ecovar_mlr", rel_path) ~ "storage_ecovar_mlr",
      grepl("ms_isotope", rel_path) ~ "ms_isotope",
      grepl("ms_chs", rel_path) ~ "ms_chs",
      grepl("ds_summary|ds_annual_ts|hydrometric", rel_path) ~ "ds_hydrometric",
      TRUE ~ "other"
    )
  ) %>%
  arrange(section, path_group, rel_path)

write_csv(manifest, file.path(manifest_dir, "output_manifest.csv"))

has_file <- function(pattern) {
  any(grepl(pattern, manifest$rel_path))
}

coverage <- tribble(
  ~analysis, ~needs_main_figure, ~needs_supp_figure, ~needs_table, ~needs_stat,
  "site_differences_anova_tukey", TRUE, TRUE, TRUE, TRUE,
  "multivariate_pca", TRUE, TRUE, TRUE, TRUE,
  "watershed_char_storage_mlr", TRUE, FALSE, TRUE, TRUE,
  "storage_ecovar_mlr", TRUE, FALSE, TRUE, TRUE,
  "dynamic_storage_summary", TRUE, TRUE, FALSE, FALSE,
  "mobile_storage_isotope", TRUE, FALSE, FALSE, FALSE,
  "mobile_storage_chs", TRUE, TRUE, FALSE, FALSE
) %>%
  mutate(
    has_main_figure = case_when(
      analysis == "site_differences_anova_tukey" ~ has_file("figs/main/ds_anova_tukey\\.(png|pdf)$"),
      analysis == "multivariate_pca" ~ has_file("figs/main/pca_biplot\\.(png|pdf)$"),
      analysis == "watershed_char_storage_mlr" ~ has_file("figs/main/watershed_char_storage_mlr_beta\\.(png|pdf)$"),
      analysis == "storage_ecovar_mlr" ~ has_file("figs/main/storage_ecovar_mlr_beta\\.(png|pdf)$"),
      analysis == "dynamic_storage_summary" ~ has_file("figs/main/ds_summary\\.(png|pdf)$"),
      analysis == "mobile_storage_isotope" ~ has_file("figs/main/ms_isotope\\.(png|pdf)$"),
      analysis == "mobile_storage_chs" ~ has_file("figs/main/ms_chs\\.(png|pdf)$"),
      TRUE ~ FALSE
    ),
    has_supp_figure = case_when(
      analysis == "site_differences_anova_tukey" ~ has_file("figs/supp/analysis/anova_tukey/"),
      analysis == "multivariate_pca" ~ has_file("figs/supp/analysis/pca/pca_scree\\.(png|pdf)$"),
      analysis == "watershed_char_storage_mlr" ~ FALSE,
      analysis == "storage_ecovar_mlr" ~ FALSE,
      analysis == "dynamic_storage_summary" ~ has_file("figs/supp/ds_annual_ts\\.(png|pdf)$"),
      analysis == "mobile_storage_chs" ~ has_file("figs/supp/ms_chs_annual_ts\\.(png|pdf)$"),
      TRUE ~ FALSE
    ),
    has_table = case_when(
      analysis == "site_differences_anova_tukey" ~ has_file("tables/anova_tukey/anova_tukey_group_table\\.csv$"),
      analysis == "multivariate_pca" ~ has_file("tables/pca/pca_variance_table\\.csv$") & has_file("tables/pca/pca_loading_table\\.csv$"),
      analysis == "watershed_char_storage_mlr" ~ has_file("tables/mlr/watershed_char_storage_mlr_table\\.csv$"),
      analysis == "storage_ecovar_mlr" ~ has_file("tables/mlr/storage_ecovar_mlr_table\\.csv$"),
      TRUE ~ FALSE
    ),
    has_stat = case_when(
      analysis == "site_differences_anova_tukey" ~ has_file("models/anova_tukey/anova_results\\.csv$"),
      analysis == "multivariate_pca" ~ has_file("models/pca/pca_variance_explained\\.csv$"),
      analysis == "watershed_char_storage_mlr" ~ has_file("models/watershed_char_storage_mlr/watershed_char_storage_mlr_summary_strict\\.csv$"),
      analysis == "storage_ecovar_mlr" ~ has_file("models/storage_ecovar_mlr/storage_ecovar_mlr_summary_strict\\.csv$"),
      TRUE ~ FALSE
    ),
    is_complete = (!needs_main_figure | has_main_figure) &
      (!needs_supp_figure | has_supp_figure) &
      (!needs_table | has_table) &
      (!needs_stat | has_stat)
  )

write_csv(coverage, file.path(manifest_dir, "analysis_coverage.csv"))
