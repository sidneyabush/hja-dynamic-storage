# Build unified subsurface storage framework axes from site-level metrics.
# Inputs: OUTPUT_DIR/master/master_site.csv.
# Author: Sidney Bush
# Date: 2026-02-21

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

rm(list = ls())

# Load project config
source("config.R")

find_read_root <- function() {
  primary <- OUTPUT_DIR
  fallback <- file.path(REPO_DIR, "outputs")
  if (file.exists(file.path(primary, "master", MASTER_SITE_FILE))) {
    return(primary)
  }
  fallback
}

is_writable_root <- function(root_dir) {
  if (!dir.exists(root_dir)) {
    ok_dir <- tryCatch(
      {
        dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
        TRUE
      },
      error = function(e) FALSE
    )
    if (!ok_dir) {
      return(FALSE)
    }
  }
  probe <- file.path(root_dir, ".write_probe")
  ok <- tryCatch(
    {
      file.create(probe)
      file.remove(probe)
      TRUE
    },
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
  isTRUE(ok)
}

find_write_root <- function() {
  primary <- OUTPUT_DIR
  fallback <- file.path(REPO_DIR, "outputs")
  if (is_writable_root(primary)) {
    return(primary)
  }
  if (is_writable_root(fallback)) {
    return(fallback)
  }
  stop("No writable output root found for unified framework.")
}

read_root <- find_read_root()
write_root <- find_write_root()

model_out_dir <- file.path(write_root, "models", "unified_framework")
table_out_dir <- file.path(write_root, "tables", "unified_framework")
for (d in c(model_out_dir, table_out_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

master_site_file <- file.path(read_root, "master", MASTER_SITE_FILE)
site_df <- read_csv(master_site_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site)

safe_z <- function(x) {
  x_num <- as.numeric(x)
  s <- suppressWarnings(sd(x_num, na.rm = TRUE))
  if (!is.finite(s) || s <= 0) {
    return(rep(NA_real_, length(x_num)))
  }
  (x_num - mean(x_num, na.rm = TRUE)) / s
}

row_mean_min <- function(df_like, min_non_na = 1L) {
  mat <- as.matrix(df_like)
  n_non_na <- rowSums(is.finite(mat))
  out <- rowMeans(mat, na.rm = TRUE)
  out[n_non_na < min_non_na] <- NA_real_
  list(mean = out, n = n_non_na)
}

# Orient metrics so larger values indicate stronger buffering, older/mobile-mixed
# behavior, or stronger baseflow partitioning.
metric_oriented <- site_df %>%
  transmute(
    site = as.character(site),
    RBI_inv = -RBI_mean,
    RCS = RCS_mean,
    FDC = FDC_mean,
    SD = SD_mean,
    WB_drawdown = -WB_mean,
    DR_inv = -DR,
    Fyw_inv = -Fyw,
    MTT1 = MTT1,
    MTT2 = MTT2,
    CHS = CHS_mean
  )

metric_z <- metric_oriented %>%
  mutate(across(-site, safe_z, .names = "{.col}_z"))

# Axis 1: Dynamic storage strength and release behavior.
dynamic_cols <- c("RBI_inv_z", "RCS_z", "FDC_z", "SD_z", "WB_drawdown_z")
dynamic_vals <- row_mean_min(metric_z[, dynamic_cols], min_non_na = 4L)

# Axis 2: Mobile-water turnover/mixing behavior.
# Combine MTT periods first so transit time is one component.
mtt_vals <- row_mean_min(metric_z[, c("MTT1_z", "MTT2_z")], min_non_na = 1L)
mobile_components <- tibble(
  MTT = mtt_vals$mean,
  DR = metric_z$DR_inv_z,
  Fyw = metric_z$Fyw_inv_z
)
mobile_vals <- row_mean_min(mobile_components, min_non_na = 2L)

# Axis 3: Flow-path partitioning (baseflow proportion).
partition_vals <- ifelse(is.finite(metric_z$CHS_z), metric_z$CHS_z, NA_real_)
partition_n <- ifelse(is.finite(partition_vals), 1L, 0L)

axes_raw <- tibble(
  site = metric_oriented$site,
  dynamic_storage_strength = dynamic_vals$mean,
  mobile_mixing = mobile_vals$mean,
  flow_path_partitioning = partition_vals,
  n_dynamic_components = dynamic_vals$n,
  n_mobile_components = mobile_vals$n,
  n_partition_components = partition_n
)

# Standardize each axis for direct comparison and optional averaging.
axes_site <- axes_raw %>%
  mutate(
    dynamic_storage_strength_z = safe_z(dynamic_storage_strength),
    mobile_mixing_z = safe_z(mobile_mixing),
    flow_path_partitioning_z = safe_z(flow_path_partitioning)
  ) %>%
  mutate(
    n_axes_available = rowSums(
      is.finite(
        cbind(dynamic_storage_strength_z, mobile_mixing_z, flow_path_partitioning_z)
      )
    ),
    unified_state_index = rowMeans(
      cbind(dynamic_storage_strength_z, mobile_mixing_z, flow_path_partitioning_z),
      na.rm = TRUE
    ),
    unified_state_index = ifelse(
      n_axes_available >= 2,
      unified_state_index,
      NA_real_
    )
  ) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site) %>%
  mutate(site = as.character(site))

axis_definitions <- tibble(
  axis = c(
    "dynamic_storage_strength",
    "dynamic_storage_strength",
    "dynamic_storage_strength",
    "dynamic_storage_strength",
    "dynamic_storage_strength",
    "mobile_mixing",
    "mobile_mixing",
    "mobile_mixing",
    "flow_path_partitioning"
  ),
  component = c(
    "RBI_mean",
    "RCS_mean",
    "FDC_mean",
    "SD_mean",
    "WB_mean",
    "MTT1/MTT2",
    "DR",
    "Fyw",
    "CHS_mean"
  ),
  direction_to_axis = c(
    "inverse",
    "direct",
    "direct",
    "direct",
    "inverse",
    "direct",
    "inverse",
    "inverse",
    "direct"
  ),
  note = c(
    "Lower flashiness contributes to higher dynamic storage strength axis score.",
    "Higher recession slope contributes to higher dynamic storage strength axis score.",
    "Flatter FDC (less negative) contributes to higher dynamic storage strength axis score.",
    "Higher storage-discharge depth contributes to higher dynamic storage strength axis score.",
    "Larger WB drawdown magnitude (-WB) contributes to higher dynamic storage strength axis score.",
    "Longer transit times contribute to higher mobile-mixing axis score.",
    "Lower damping ratio contributes to higher mobile-mixing axis score.",
    "Lower young-water fraction contributes to higher mobile-mixing axis score.",
    "Higher CHS contributes to higher flow-path partitioning axis score."
  )
)

axis_corr <- axes_site %>%
  select(
    dynamic_storage_strength_z,
    mobile_mixing_z,
    flow_path_partitioning_z,
    unified_state_index
  ) %>%
  cor(use = "pairwise.complete.obs")

axis_corr_long <- as.data.frame(as.table(axis_corr)) %>%
  transmute(axis_x = Var1, axis_y = Var2, pearson_r = Freq)

write_csv(
  axes_site,
  file.path(model_out_dir, "unified_framework_site_axes.csv")
)
write_csv(
  axis_definitions,
  file.path(model_out_dir, "unified_framework_axis_components.csv")
)
write_csv(
  axis_corr_long,
  file.path(model_out_dir, "unified_framework_axis_correlations.csv")
)

# Mirror core table into tables output for manuscript lookup.
write_csv(
  axes_site,
  file.path(table_out_dir, "unified_framework_site_axes.csv")
)
