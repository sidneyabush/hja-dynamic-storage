# sensitivity checks for the unified subsurface storage framework.
# inputs: output_dir/master/master_site.csv; unified_framework site axes output.
# author: sidney bush
# date: 2026-02-21

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

rm(list = ls())

# load project config
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
    if (!ok_dir) return(FALSE)
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
  if (is_writable_root(primary)) return(primary)
  if (is_writable_root(fallback)) return(fallback)
  stop("No writable output root found for unified-framework sensitivity.")
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
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

safe_z <- function(x) {
  x_num <- as.numeric(x)
  s <- suppressWarnings(sd(x_num, na.rm = TRUE))
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x_num)))
  (x_num - mean(x_num, na.rm = TRUE)) / s
}

z_from_ref <- function(x, ref) {
  mu <- mean(ref, na.rm = TRUE)
  s <- sd(ref, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mu) / s
}

row_mean_min <- function(df_like, min_non_na = 1L) {
  mat <- as.matrix(df_like)
  n_non_na <- rowSums(is.finite(mat))
  out <- rowMeans(mat, na.rm = TRUE)
  out[n_non_na < min_non_na] <- NA_real_
  list(mean = out, n = n_non_na)
}

row_median_min <- function(df_like, min_non_na = 1L) {
  mat <- as.matrix(df_like)
  n_non_na <- rowSums(is.finite(mat))
  out <- apply(mat, 1, function(r) {
    rr <- r[is.finite(r)]
    if (length(rr) < min_non_na) return(NA_real_)
    median(rr)
  })
  list(median = as.numeric(out), n = n_non_na)
}

scale_with_ref <- function(x, ref) {
  mu <- mean(ref, na.rm = TRUE)
  s <- sd(ref, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mu) / s
}

build_framework <- function(
  z_df,
  include_wb = TRUE,
  include_mtt = TRUE,
  include_partition = TRUE,
  axis_agg = c("mean", "median"),
  min_axes = 2L,
  dynamic_min_required = NULL,
  mobile_min_required = NULL,
  axis_scale_ref = NULL
) {
  axis_agg <- match.arg(axis_agg)

  dynamic_cols <- c("RBI_inv_z", "RCS_z", "FDC_z", "SD_z")
  if (isTRUE(include_wb)) {
    dynamic_cols <- c(dynamic_cols, "WB_drawdown_z")
  }
  if (is.null(dynamic_min_required)) {
    dynamic_min_required <- max(1L, length(dynamic_cols) - 1L)
  }
  dynamic_vals <- row_mean_min(z_df[, dynamic_cols, drop = FALSE], dynamic_min_required)

  mtt_vals <- row_mean_min(z_df[, c("MTT1_z", "MTT2_z"), drop = FALSE], 1L)
  if (isTRUE(include_mtt)) {
    mobile_components <- tibble(
      MTT = mtt_vals$mean,
      DR = z_df$DR_inv_z,
      Fyw = z_df$Fyw_inv_z
    )
    if (is.null(mobile_min_required)) mobile_min_required <- 2L
  } else {
    mobile_components <- tibble(
      DR = z_df$DR_inv_z,
      Fyw = z_df$Fyw_inv_z
    )
    if (is.null(mobile_min_required)) mobile_min_required <- 2L
  }
  mobile_vals <- row_mean_min(mobile_components, mobile_min_required)

  partition_vals <- if (isTRUE(include_partition)) z_df$CHS_z else rep(NA_real_, nrow(z_df))
  partition_n <- ifelse(is.finite(partition_vals), 1L, 0L)

  axes_raw <- tibble(
    dynamic_raw = dynamic_vals$mean,
    mobile_raw = mobile_vals$mean,
    partition_raw = partition_vals,
    n_dynamic_components = dynamic_vals$n,
    n_mobile_components = mobile_vals$n,
    n_partition_components = partition_n
  )

  if (is.null(axis_scale_ref)) {
    dynamic_z <- safe_z(axes_raw$dynamic_raw)
    mobile_z <- safe_z(axes_raw$mobile_raw)
    partition_z <- safe_z(axes_raw$partition_raw)
  } else {
    dynamic_z <- scale_with_ref(axes_raw$dynamic_raw, axis_scale_ref$dynamic_raw)
    mobile_z <- scale_with_ref(axes_raw$mobile_raw, axis_scale_ref$mobile_raw)
    partition_z <- scale_with_ref(axes_raw$partition_raw, axis_scale_ref$partition_raw)
  }

  axis_mat <- cbind(dynamic_z, mobile_z, partition_z)
  n_axes_available <- rowSums(is.finite(axis_mat))

  if (identical(axis_agg, "mean")) {
    unified <- rowMeans(axis_mat, na.rm = TRUE)
  } else {
    unified <- apply(axis_mat, 1, function(r) {
      rr <- r[is.finite(r)]
      if (length(rr) == 0L) return(NA_real_)
      median(rr)
    })
  }
  unified[n_axes_available < min_axes] <- NA_real_

  tibble(
    dynamic_storage_strength_z = dynamic_z,
    mobile_mixing_z = mobile_z,
    flow_path_partitioning_z = partition_z,
    n_dynamic_components = axes_raw$n_dynamic_components,
    n_mobile_components = axes_raw$n_mobile_components,
    n_partition_components = axes_raw$n_partition_components,
    n_axes_available = n_axes_available,
    unified_state_index = unified,
    dynamic_raw = axes_raw$dynamic_raw,
    mobile_raw = axes_raw$mobile_raw,
    partition_raw = axes_raw$partition_raw
  )
}

top_k_sites <- function(df, col, k = 3L) {
  df %>%
    filter(is.finite(.data[[col]])) %>%
    arrange(desc(.data[[col]])) %>%
    slice_head(n = k) %>%
    pull(site) %>%
    as.character()
}

bottom_k_sites <- function(df, col, k = 3L) {
  df %>%
    filter(is.finite(.data[[col]])) %>%
    arrange(.data[[col]]) %>%
    slice_head(n = k) %>%
    pull(site) %>%
    as.character()
}

compare_to_baseline <- function(df, baseline_col, scenario_col) {
  x <- df[[baseline_col]]
  y <- df[[scenario_col]]
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 3L) {
    return(tibble(
      n_sites = sum(valid),
      spearman_rho = NA_real_,
      kendall_tau = NA_real_,
      mean_abs_diff = NA_real_,
      mean_abs_rank_shift = NA_real_,
      top3_overlap = NA_real_,
      bottom3_overlap = NA_real_
    ))
  }

  base_df <- tibble(site = df$site, score = x) %>% filter(is.finite(score))
  scen_df <- tibble(site = df$site, score = y) %>% filter(is.finite(score))
  base_rank <- base_df %>% mutate(rank = rank(-score, ties.method = "average")) %>% select(site, rank)
  scen_rank <- scen_df %>% mutate(rank = rank(-score, ties.method = "average")) %>% select(site, rank)
  rank_diff <- base_rank %>%
    inner_join(scen_rank, by = "site", suffix = c("_base", "_scen")) %>%
    mutate(abs_shift = abs(rank_base - rank_scen))

  base_top <- top_k_sites(df, baseline_col, 3L)
  scen_top <- top_k_sites(df, scenario_col, 3L)
  base_bottom <- bottom_k_sites(df, baseline_col, 3L)
  scen_bottom <- bottom_k_sites(df, scenario_col, 3L)

  tibble(
    n_sites = sum(valid),
    spearman_rho = suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = "spearman")),
    kendall_tau = suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = "kendall")),
    mean_abs_diff = mean(abs(x[valid] - y[valid]), na.rm = TRUE),
    mean_abs_rank_shift = mean(rank_diff$abs_shift, na.rm = TRUE),
    top3_overlap = length(intersect(base_top, scen_top)),
    bottom3_overlap = length(intersect(base_bottom, scen_bottom))
  )
}

master_site_file <- file.path(read_root, "master", MASTER_SITE_FILE)
site_df <- read_csv(master_site_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC)) %>%
  arrange(site) %>%
  mutate(site = as.character(site))

ori <- site_df %>%
  transmute(
    site,
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

z_full <- ori %>%
  mutate(across(-site, safe_z, .names = "{.col}_z"))

baseline_axes <- build_framework(
  z_df = z_full,
  include_wb = TRUE,
  include_mtt = TRUE,
  include_partition = TRUE,
  axis_agg = "mean",
  min_axes = 2L
)

scoring_df <- tibble(site = ori$site) %>%
  bind_cols(
    baseline = baseline_axes$unified_state_index
  )

scenario_list <- list(
  no_wb = list(include_wb = FALSE, include_mtt = TRUE, include_partition = TRUE, axis_agg = "mean"),
  no_mtt = list(include_wb = TRUE, include_mtt = FALSE, include_partition = TRUE, axis_agg = "mean"),
  no_partition = list(include_wb = TRUE, include_mtt = TRUE, include_partition = FALSE, axis_agg = "mean"),
  median_axes = list(include_wb = TRUE, include_mtt = TRUE, include_partition = TRUE, axis_agg = "median")
)

for (nm in names(scenario_list)) {
  opts <- scenario_list[[nm]]
  sc <- build_framework(
    z_df = z_full,
    include_wb = opts$include_wb,
    include_mtt = opts$include_mtt,
    include_partition = opts$include_partition,
    axis_agg = opts$axis_agg,
    min_axes = 2L
  )
  scoring_df[[nm]] <- sc$unified_state_index
}

# alternative aggregation: equal-weight component score (no axis grouping).
mtt_combined <- row_mean_min(z_full[, c("MTT1_z", "MTT2_z"), drop = FALSE], 1L)$mean
component_score <- row_mean_min(
  tibble(
    RBI_inv = z_full$RBI_inv_z,
    RCS = z_full$RCS_z,
    FDC = z_full$FDC_z,
    SD = z_full$SD_z,
    WB_drawdown = z_full$WB_drawdown_z,
    MTT = mtt_combined,
    DR_inv = z_full$DR_inv_z,
    Fyw_inv = z_full$Fyw_inv_z,
    CHS = z_full$CHS_z
  ),
  min_non_na = 5L
)$mean
scoring_df$metric_equal_weight <- safe_z(component_score)

sensitivity_summary <- bind_rows(lapply(setdiff(names(scoring_df), c("site", "baseline")), function(nm) {
  compare_to_baseline(scoring_df, "baseline", nm) %>%
    mutate(scenario = nm, .before = 1)
}))

ranking_table <- scoring_df %>%
  pivot_longer(cols = -site, names_to = "scenario", values_to = "score") %>%
  group_by(scenario) %>%
  mutate(
    rank_desc = ifelse(is.finite(score), rank(-score, ties.method = "average"), NA_real_)
  ) %>%
  ungroup() %>%
  arrange(scenario, rank_desc)

# leave-one-site-out stability for baseline configuration.
loo_records <- lapply(ori$site, function(site_holdout) {
  train_raw <- ori %>% filter(site != site_holdout)
  test_raw <- ori %>% filter(site == site_holdout)

  train_z <- train_raw %>%
    mutate(across(-site, ~ z_from_ref(.x, .x), .names = "{.col}_z"))
  test_z <- test_raw %>%
    mutate(across(-site, ~ z_from_ref(.x, train_raw[[cur_column()]]), .names = "{.col}_z"))

  train_fw <- build_framework(
    z_df = train_z,
    include_wb = TRUE,
    include_mtt = TRUE,
    include_partition = TRUE,
    axis_agg = "mean",
    min_axes = 2L
  )
  axis_ref <- list(
    dynamic_raw = train_fw$dynamic_raw,
    mobile_raw = train_fw$mobile_raw,
    partition_raw = train_fw$partition_raw
  )
  test_fw <- build_framework(
    z_df = test_z,
    include_wb = TRUE,
    include_mtt = TRUE,
    include_partition = TRUE,
    axis_agg = "mean",
    min_axes = 2L,
    axis_scale_ref = axis_ref
  )

  tibble(
    site = site_holdout,
    loo_unified_state_index = test_fw$unified_state_index,
    loo_n_axes_available = test_fw$n_axes_available,
    loo_dynamic_components = test_fw$n_dynamic_components,
    loo_mobile_components = test_fw$n_mobile_components,
    loo_partition_components = test_fw$n_partition_components
  )
}) %>%
  bind_rows()

loo_results <- scoring_df %>%
  select(site, baseline) %>%
  left_join(loo_records, by = "site") %>%
  mutate(
    delta_loo_minus_baseline = loo_unified_state_index - baseline,
    abs_delta_loo_minus_baseline = abs(delta_loo_minus_baseline)
  )

loo_valid <- loo_results %>% filter(is.finite(baseline), is.finite(loo_unified_state_index))
loo_summary <- tibble(
  n_sites = nrow(loo_valid),
  spearman_rho = suppressWarnings(cor(loo_valid$baseline, loo_valid$loo_unified_state_index, method = "spearman")),
  kendall_tau = suppressWarnings(cor(loo_valid$baseline, loo_valid$loo_unified_state_index, method = "kendall")),
  mean_abs_delta = mean(loo_valid$abs_delta_loo_minus_baseline, na.rm = TRUE),
  max_abs_delta = max(loo_valid$abs_delta_loo_minus_baseline, na.rm = TRUE)
)

write_csv(
  scoring_df,
  file.path(model_out_dir, "unified_framework_sensitivity_site_scores.csv")
)
write_csv(
  sensitivity_summary,
  file.path(model_out_dir, "unified_framework_sensitivity_summary.csv")
)
write_csv(
  ranking_table,
  file.path(model_out_dir, "unified_framework_sensitivity_rankings.csv")
)
write_csv(
  loo_results,
  file.path(model_out_dir, "unified_framework_loo_site_scores.csv")
)
write_csv(
  loo_summary,
  file.path(model_out_dir, "unified_framework_loo_summary.csv")
)

if (isTRUE(WRITE_TABLE_OUTPUTS)) {
  write_csv(
    sensitivity_summary,
    file.path(table_out_dir, "unified_framework_sensitivity_summary.csv")
  )
  write_csv(
    loo_results,
    file.path(table_out_dir, "unified_framework_loo_site_scores.csv")
  )
  write_csv(
    loo_summary,
    file.path(table_out_dir, "unified_framework_loo_summary.csv")
  )
}
