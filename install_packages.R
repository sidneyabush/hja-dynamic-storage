# install required R packages for the workflow

# inputs:
# none

# outputs:
# required R packages installed locally

# author: Sidney Bush
# date: 2026-02-13

required_packages <- c(
  "car",
  "cowplot",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "ggtext",
  "librarian",
  "lubridate",
  "MASS",
  "multcompView",
  "patchwork",
  "purrr",
  "readr",
  "scales",
  "tibble",
  "tidyr",
  "zoo"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

message("Package setup complete. Next run: Rscript run_all.R")
