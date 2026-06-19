# install the R packages used by this code
# run this once before run_all.R if the required packages are not installed

required_packages <- c(
  "car",
  "cowplot",
  "dplyr",
  "ggplot2",
  "ggrepel",
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
