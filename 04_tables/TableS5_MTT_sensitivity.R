librarian::shelf(readr, cran_repo = "https://cloud.r-project.org")

rm(list = ls())

source("config.R")

source_file <- file.path(OUTPUT_DIR, "MTT_sensitivity", "TableS5_MTT_sensitivity.csv")
if (!file.exists(source_file)) {
  stop("Missing required MTT sensitivity table: ", source_file)
}

dir.create(MS_TABLES_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

table_s5 <- read_csv(source_file, show_col_types = FALSE)
write_csv(table_s5, file.path(MS_TABLES_SUPP_DIR, "TableS5_MTT_sensitivity.csv"))
