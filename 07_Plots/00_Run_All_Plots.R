# -----------------------------------------------------------------------------
# Run All Figures
# -----------------------------------------------------------------------------
# Purpose: Execute all plotting scripts to generate complete figure set
#
# Usage:
#   source("07_Plots/00_Run_All_Plots.R")
#   OR
#   Rscript 07_Plots/00_Run_All_Plots.R
#
# Output Structure:
#   Main_Text/     - Publication figures (Fig3, Fig4, Fig5)
#   Supplement/    - Supplementary time series figures
#
# Author: Sidney Bush
# Date: 2026-02-02
# -----------------------------------------------------------------------------

# Get script directory
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
if (is.null(script_dir) || script_dir == "" || script_dir == ".") {
  script_dir <- file.path(getwd(), "07_Plots")
}

cat("Script directory:", script_dir, "\n\n")

# -----------------------------------------------------------------------------
# RUN MAIN FIGURES SCRIPT
# -----------------------------------------------------------------------------

script_path <- file.path(script_dir, "Main_Figures.R")

if (file.exists(script_path)) {
  cat("================================================================\n")
  cat("Running: Main_Figures.R\n")
  cat("================================================================\n\n")

  tryCatch({
    source(script_path)
    cat("\n[OK] Main_Figures.R completed successfully\n\n")
  }, error = function(e) {
    cat("\n[ERROR] Main_Figures.R failed:", conditionMessage(e), "\n\n")
  })
} else {
  cat("[ERROR] Main_Figures.R not found at:", script_path, "\n")
}

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("  All Plots Complete\n")
cat("================================================================\n\n")

# Source config to get output directory
config_path <- file.path(dirname(script_dir), "config.R")
if (file.exists(config_path)) {
  source(config_path)

  cat("Output directories:\n")

  main_dir <- file.path(FIGURES_DIR, "Main_Text")
  supp_dir <- file.path(FIGURES_DIR, "Supplement")

  if (dir.exists(main_dir)) {
    n_files <- length(list.files(main_dir, pattern = "\\.(png|pdf)$"))
    cat(sprintf("  - Main_Text: %d files\n", n_files))
  }

  if (dir.exists(supp_dir)) {
    n_files <- length(list.files(supp_dir, pattern = "\\.(png|pdf)$"))
    cat(sprintf("  - Supplement: %d files\n", n_files))
  }
}

cat("\nDone.\n")
