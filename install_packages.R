# install the R packages used by this code

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

if (length(missing_packages) == 0) {
  message("All required packages are already installed.")
} else {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
