# =============================================================================
# Install Required R Packages for HJA Dynamic Storage Analysis
# =============================================================================
# Run this script once to install all required packages:
#   source("install_packages.R")
# =============================================================================

# List of required packages
required_packages <- c(
  # Data manipulation
  "dplyr",
  "tidyr",
  "readr",
  "lubridate",
  "tibble",

  # Visualization
  "ggplot2",
  "patchwork",
  "colorspace",
  "scales",
  "GGally",
  "ggcorrplot",
  "ggrepel",

  # Time series and signal processing
  "zoo",
  "pracma",

  # Statistics
  "MASS",       # stepAIC
  "car",        # VIF
  "vegan"       # RDA, ordination
)

# Function to install missing packages
install_if_missing <- function(packages) {
  installed <- rownames(installed.packages())
  to_install <- packages[!(packages %in% installed)]

  if (length(to_install) > 0) {
    cat("Installing", length(to_install), "packages:\n")
    cat(paste(" -", to_install, collapse = "\n"), "\n\n")
    install.packages(to_install, dependencies = TRUE)
  } else {
    cat("All required packages are already installed.\n")
  }

  # Check for successful installation
  still_missing <- packages[!(packages %in% rownames(installed.packages()))]
  if (length(still_missing) > 0) {
    warning("The following packages failed to install:\n",
            paste(" -", still_missing, collapse = "\n"))
  } else {
    cat("\nAll packages installed successfully!\n")
  }
}

# Run installation
install_if_missing(required_packages)

# Print session info
cat("\n=== Session Info ===\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
