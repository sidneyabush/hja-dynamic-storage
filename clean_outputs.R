#!/usr/bin/env Rscript

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

allow_clean <- tolower(Sys.getenv("HJA_ALLOW_CLEAN", unset = "false")) == "true"
if (!allow_clean) {
  stop("Refusing to clean outputs. Re-run with HJA_ALLOW_CLEAN=true.")
}

if (!dir.exists(OUTPUT_DIR)) {
  message("Output directory does not exist: ", OUTPUT_DIR)
  quit(save = "no", status = 0)
}

# Safety guard: only clean the consolidated workflow output root.
if (!grepl("final_workflow$", OUTPUT_DIR)) {
  stop("Safety check failed. OUTPUT_DIR must end with 'final_workflow'.")
}

targets <- c("metrics", "master", "stats", "tables", "figs", "MET")
for (nm in targets) {
  p <- file.path(OUTPUT_DIR, nm)
  if (dir.exists(p)) unlink(p, recursive = TRUE, force = TRUE)
}

# Recreate expected tree via config side effects.
source(file.path(repo_root, "config.R"))

message("Cleaned output subfolders in: ", OUTPUT_DIR)
