# Run all MLR workflows and export one manuscript-ready summary table.
# Inputs: No direct CSV file reads in this script.
# Author: Sidney Bush
# Date: 2026-02-14

run_subscript <- function(path) {
  status <- system2("Rscript", shQuote(path), stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Subscript failed: ", path)
  }
}

scripts <- c(
  "02_stats/watershed_char_storage_mlr.R",
  "02_stats/storage_ecovar_mlr.R",
  "02_stats/export_mlr_main_results_table.R"
)

for (s in scripts) {
  if (!file.exists(s)) {
    stop("Missing subscript: ", s)
  }
  run_subscript(s)
}
