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
dir_targets <- c(model_out_dir)
if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  dir_targets <- c(dir_targets, table_out_dir)
}
for (d in dir_targets) {
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

classify_geology <- function(lava1_per, lava2_per, ash_per, pyro_per, tie_tol = 5) {
  lava_total <- suppressWarnings(as.numeric(lava1_per)) + suppressWarnings(as.numeric(lava2_per))
  ash <- suppressWarnings(as.numeric(ash_per))
  pyro <- suppressWarnings(as.numeric(pyro_per))

  out <- rep(NA_character_, length(lava_total))
  for (i in seq_along(out)) {
    vals <- c(
      "Lava-dominant" = lava_total[i],
      "Ash-dominant" = ash[i],
      "Pyroclastic-dominant" = pyro[i]
    )
    if (!any(is.finite(vals))) next
    vals[!is.finite(vals)] <- -Inf
    ranked <- sort(vals, decreasing = TRUE)
    if (length(ranked) >= 2 && is.finite(ranked[2]) && (ranked[1] - ranked[2]) <= tie_tol) {
      out[i] <- "Mixed"
    } else {
      out[i] <- names(which.max(vals))
    }
  }
  out
}

classify_geomorphology <- function(harvest_pct, landslide_total_pct) {
  harvest <- suppressWarnings(as.numeric(harvest_pct))
  landslide <- suppressWarnings(as.numeric(landslide_total_pct))

  dplyr::case_when(
    is.finite(harvest) & is.finite(landslide) & harvest >= 50 & landslide >= 30 ~ "Harvest+landslide",
    is.finite(harvest) & harvest >= 50 ~ "Harvest-dominated",
    is.finite(landslide) & landslide >= 30 ~ "Landslide-dominated",
    is.finite(harvest) | is.finite(landslide) ~ "Low-disturbance",
    TRUE ~ NA_character_
  )
}

# Orient metrics so larger values indicate stronger buffering or stronger mobile-storage influence.
metric_oriented <- site_df %>%
  transmute(
    site = as.character(site),
    RBI_inv = -RBI_mean,
    RCS = RCS_mean,
    FDC = FDC_mean,
    SD = SD_mean,
    WB_depletion_mag = abs(WB_mean),
    DR_inv = -DR,
    Fyw_inv = -Fyw,
    MTT1 = MTT1,
    MTT2 = MTT2,
    CHS = CHS_mean
  )

metric_z <- metric_oriented %>%
  mutate(across(-site, safe_z, .names = "{.col}_z"))

site_classes <- site_df %>%
  transmute(
    site = as.character(site),
    geology_class = classify_geology(Lava1_per, Lava2_per, Ash_Per, Pyro_per),
    geomorphology_class = classify_geomorphology(Harvest, Landslide_Total),
    Landslide_Total = suppressWarnings(as.numeric(Landslide_Total)),
    Lava1_per = suppressWarnings(as.numeric(Lava1_per)),
    Lava2_per = suppressWarnings(as.numeric(Lava2_per)),
    Ash_Per = suppressWarnings(as.numeric(Ash_Per)),
    Pyro_per = suppressWarnings(as.numeric(Pyro_per))
  )

geology_cols <- c("Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")
geology_input <- site_classes %>%
  select(site, all_of(geology_cols))

geology_complete <- geology_input %>%
  filter(if_all(all_of(geology_cols), is.finite))

if (nrow(geology_complete) >= 3) {
  geology_pca <- prcomp(
    geology_complete[, geology_cols],
    center = TRUE,
    scale. = TRUE
  )

  geo_scores_mat <- geology_pca$x
  geology_scores <- tibble(
    site = geology_complete$site,
    geology_pc1 = as.numeric(geo_scores_mat[, 1]),
    geology_pc2 = if (ncol(geo_scores_mat) >= 2) {
      as.numeric(geo_scores_mat[, 2])
    } else {
      NA_real_
    }
  )

  geo_load_mat <- geology_pca$rotation
  geology_loadings <- tibble(
    variable = rownames(geo_load_mat),
    geology_pc1_loading = as.numeric(geo_load_mat[, 1]),
    geology_pc2_loading = if (ncol(geo_load_mat) >= 2) {
      as.numeric(geo_load_mat[, 2])
    } else {
      NA_real_
    }
  )

  geo_var <- geology_pca$sdev^2
  geo_var_explained <- geo_var / sum(geo_var)
  geology_variance <- tibble(
    pc = paste0("PC", seq_along(geo_var_explained)),
    variance_explained = as.numeric(geo_var_explained),
    cumulative_variance = as.numeric(cumsum(geo_var_explained))
  )
} else {
  geology_scores <- geology_input %>%
    transmute(site, geology_pc1 = NA_real_, geology_pc2 = NA_real_)
  geology_loadings <- tibble(
    variable = geology_cols,
    geology_pc1_loading = NA_real_,
    geology_pc2_loading = NA_real_
  )
  geology_variance <- tibble(
    pc = paste0("PC", seq_len(length(geology_cols))),
    variance_explained = NA_real_,
    cumulative_variance = NA_real_
  )
}

# Axis 1: Dynamic storage strength and release behavior.
dynamic_cols <- c("RBI_inv_z", "RCS_z", "FDC_z", "SD_z", "WB_depletion_mag_z")
dynamic_vals <- row_mean_min(metric_z[, dynamic_cols], min_non_na = 4L)

# Axis 2: Mobile-storage behavior.
# Combine MTT periods first so transit time is one component.
mtt_vals <- row_mean_min(metric_z[, c("MTT1_z", "MTT2_z")], min_non_na = 1L)
mobile_components_with_chs <- tibble(
  MTT = mtt_vals$mean,
  DR = metric_z$DR_inv_z,
  Fyw = metric_z$Fyw_inv_z,
  CHS = metric_z$CHS_z
)
mobile_components_no_chs <- tibble(
  MTT = mtt_vals$mean,
  DR = metric_z$DR_inv_z,
  Fyw = metric_z$Fyw_inv_z
)
mobile_vals_with_chs <- row_mean_min(mobile_components_with_chs, min_non_na = 2L)
mobile_vals_no_chs <- row_mean_min(mobile_components_no_chs, min_non_na = 2L)

axes_raw <- tibble(
  site = metric_oriented$site,
  dynamic_storage_strength = dynamic_vals$mean,
  mobile_mixing = mobile_vals_with_chs$mean,
  mobile_mixing_with_chs = mobile_vals_with_chs$mean,
  mobile_mixing_no_chs = mobile_vals_no_chs$mean,
  flow_path_partitioning = metric_oriented$CHS,
  CHS_mean = metric_oriented$CHS,
  n_dynamic_components = dynamic_vals$n,
  n_mobile_components = mobile_vals_with_chs$n,
  n_mobile_components_with_chs = mobile_vals_with_chs$n,
  n_mobile_components_no_chs = mobile_vals_no_chs$n,
  n_partition_components = ifelse(is.finite(metric_oriented$CHS), 1L, 0L)
)

# Standardize each axis for direct comparison and optional averaging.
axes_site <- axes_raw %>%
  mutate(
    dynamic_storage_strength_z = safe_z(dynamic_storage_strength),
    mobile_mixing_z = safe_z(mobile_mixing),
    mobile_mixing_with_chs_z = safe_z(mobile_mixing_with_chs),
    mobile_mixing_no_chs_z = safe_z(mobile_mixing_no_chs),
    flow_path_partitioning_z = safe_z(flow_path_partitioning)
  ) %>%
  mutate(
    n_axes_available = rowSums(
      is.finite(
        cbind(dynamic_storage_strength_z, mobile_mixing_with_chs_z)
      )
    ),
    unified_state_index = rowMeans(
      cbind(dynamic_storage_strength_z, mobile_mixing_with_chs_z),
      na.rm = TRUE
    ),
    unified_state_index = ifelse(
      n_axes_available >= 2,
      unified_state_index,
      NA_real_
    )
  ) %>%
  left_join(site_classes, by = "site") %>%
  left_join(geology_scores, by = "site") %>%
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
    "mobile_mixing"
  ),
  component = c(
    "RBI_mean",
    "RCS_mean",
    "FDC_mean",
    "SD_mean",
    "WB_depletion_mag",
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
    "direct",
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
    "Larger WB depletion magnitude (abs(WB)) contributes to higher dynamic storage strength axis score.",
    "Longer transit times contribute to higher mobile-storage axis score.",
    "Lower damping ratio contributes to higher mobile-storage axis score.",
    "Lower young-water fraction contributes to higher mobile-storage axis score.",
    "Higher CHS contributes to higher mobile-storage axis score."
  )
)

axis_corr <- axes_site %>%
  select(
    dynamic_storage_strength_z,
    mobile_mixing_z,
    mobile_mixing_no_chs_z,
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
write_csv(
  geology_loadings,
  file.path(model_out_dir, "geology_composition_pca_loadings.csv")
)
write_csv(
  geology_variance,
  file.path(model_out_dir, "geology_composition_pca_variance.csv")
)

if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  # Mirror core table into tables output for manuscript lookup.
  write_csv(
    axes_site,
    file.path(table_out_dir, "unified_framework_site_axes.csv")
  )
}
