# =============================================================================
# Run All Publication Plots
# =============================================================================
# Purpose: Execute all plotting scripts to generate complete figure set
#
# Usage:
#   source("07_Plots/00_Run_All_Plots.R")
#
# Scripts executed:
#   1. Hydrometric_Plots.R - Storage metrics time series and summaries
#   2. Recession_Curves.R - Recession curve log-log plots
#   3. Publication_Figures.R - Final publication-ready figures
#
# Prerequisites:
#   - Run all analysis scripts in 01-06 folders first
#   - Ensure config.R is properly configured
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  HJA Dynamic Storage - Generate All Plots\n")
cat("================================================================\n\n")

# Get script directory
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir) || script_dir == "") script_dir <- getwd()

# Verify we're in the right directory
if (!file.exists(file.path(script_dir, "Hydrometric_Plots.R"))) {
  script_dir <- file.path(getwd(), "07_Plots")
}

cat("Script directory:", script_dir, "\n\n")

# =============================================================================
# RUN PLOTTING SCRIPTS
# =============================================================================

scripts_to_run <- c(
  "Hydrometric_Plots.R",
  "Recession_Curves.R",
  "Publication_Figures.R"
)

for (script_name in scripts_to_run) {
  script_path <- file.path(script_dir, script_name)

  if (file.exists(script_path)) {
    cat("================================================================\n")
    cat("Running:", script_name, "\n")
    cat("================================================================\n\n")

    tryCatch({
      source(script_path)
      cat("\n[OK]", script_name, "completed successfully\n\n")
    }, error = function(e) {
      cat("\n[ERROR]", script_name, "failed:", conditionMessage(e), "\n\n")
    })
  } else {
    cat("[SKIP]", script_name, "not found\n\n")
  }
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  All Plots Complete\n")
cat("================================================================\n\n")

# Source config to get output directory
config_path <- file.path(dirname(script_dir), "config.R")
if (file.exists(config_path)) {
  source(config_path)

  cat("Output directories:\n")

  dirs <- c(
    "Hydrometric" = file.path(OUTPUT_DIR, "Hydrometric"),
    "Recession" = file.path(OUTPUT_DIR, "Recession"),
    "Publication" = file.path(OUTPUT_DIR, "Publication_Figures")
  )

  for (name in names(dirs)) {
    if (dir.exists(dirs[name])) {
      n_files <- length(list.files(dirs[name], pattern = "\\.(png|pdf)$"))
      cat(sprintf("  - %s: %d files\n", name, n_files))
    }
  }
}

cat("\nDone.\n")
