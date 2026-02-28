# test for significant differences among sites in storage metrics.
# inputs: no direct csv file reads in this script.
# author: sidney bush
# date: 2026-01-23

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(multcompView)

# clear environment
rm(list = ls())

# source configuration (paths, site definitions, water year range)
# get script directory (works with source() and rscript)
# load project config
source("config.R")


theme_set(theme_pub(base_size = 12))

# use configuration values
site_order <- SITE_ORDER_HYDROMETRIC
base_dir   <- BASE_DATA_DIR
output_dir <- OUTPUT_DIR
output_dir <- OUT_STATS_ANOVA_DIR

# create output directory if needed
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# load annual storage metrics

master_dir <- file.path(OUTPUT_DIR, "master")
annual_file <- file.path(master_dir, MASTER_ANNUAL_FILE)

HJA_annual <- read_csv(
  annual_file,
  show_col_types = FALSE
) %>%
  mutate(site = standardize_site_code(site)) %>%
  filter(!site %in% SITE_EXCLUDE_STANDARD) %>%
  mutate(site = factor(site, levels = site_order))

# for anova/tukey only, use annual fdc slopes (site-year) from dynamic metrics output.
fdc_wy_file <- file.path(OUT_MET_DYNAMIC_DIR, "fdc_slopes_wy.csv")
fdc_annual <- if (file.exists(fdc_wy_file)) {
  fdc_raw <- read_csv(fdc_wy_file, show_col_types = FALSE)
  if ("WaterYear" %in% names(fdc_raw) && !("wateryear" %in% names(fdc_raw))) {
    fdc_raw <- fdc_raw %>% rename(wateryear = WaterYear)
  }
  if ("Slope" %in% names(fdc_raw) && !("value" %in% names(fdc_raw))) {
    fdc_raw <- fdc_raw %>% rename(value = Slope)
  }
  fdc_raw %>%
    transmute(
      site = standardize_site_code(site),
      wateryear = as.integer(wateryear),
      value = as.numeric(value)
    ) %>%
    filter(
      site %in% site_order,
      wateryear >= WY_START,
      wateryear <= WY_END,
      is.finite(value)
    ) %>%
    mutate(site = factor(site, levels = site_order))
} else {
  tibble(
    site = factor(levels = site_order),
    wateryear = integer(),
    value = numeric()
  )
}

# select storage metrics for anova
# q5norm and cv_q5norm are response variables, so they are excluded here.
# storage metrics by type (using method abbreviations):
#   dynamic: rbi, rcs, fdc, sd
#   mobile: chs - note: mtt, fyw, dr are site-level only
#   extended dynamic: wb

storage_metrics <- STORAGE_METRIC_ORDER[
  STORAGE_METRIC_ORDER %in% c("RBI", "RCS", "FDC", "SD", "WB", "CHS")
]

get_metric_anova_data <- function(metric_name) {
  if (identical(metric_name, "FDC")) {
    return(
      fdc_annual %>%
        select(site, value)
    )
  }

  if (!metric_name %in% colnames(HJA_annual)) {
    return(
      tibble(
        site = factor(levels = site_order),
        value = numeric()
      )
    )
  }

  HJA_annual %>%
    select(site, value = all_of(metric_name)) %>%
    filter(is.finite(value), !is.na(site))
}

metric_data_long <- bind_rows(lapply(storage_metrics, function(metric_name) {
  metric_df <- get_metric_anova_data(metric_name)
  if (nrow(metric_df) == 0) {
    return(NULL)
  }
  metric_df %>% mutate(metric = metric_name)
}))

if (nrow(metric_data_long) == 0) {
  stop("No valid metric data available for ANOVA/Tukey.")
}

# run anova for each metric

anova_results <- data.frame()

for (metric in storage_metrics) {
  # prepare data for anova
  anova_data <- metric_data_long %>%
    filter(metric == !!metric) %>%
    select(site, value)

  # check if enough data
  if (nrow(anova_data) < 10 || dplyr::n_distinct(anova_data$site) < 2) next

  # run one-way anova
  anova_fit <- aov(value ~ site, data = anova_data)
  anova_summary <- summary(anova_fit)

  # extract f-statistic and p-value
  f_stat <- anova_summary[[1]]$`F value`[1]
  p_value <- anova_summary[[1]]$`Pr(>F)`[1]

  # store results
  anova_results <- rbind(anova_results, data.frame(
    metric = metric,
    F_statistic = f_stat,
    p_value = p_value,
    df_between = anova_summary[[1]]$Df[1],
    df_within = anova_summary[[1]]$Df[2],
    significant = ifelse(p_value < 0.05, "Yes", "No")
  ))
}

# save anova results
anova_results <- anova_results %>%
  arrange(factor(metric, levels = storage_metrics))

write.csv(anova_results,
          file.path(output_dir, "anova_results.csv"),
          row.names = FALSE)

# run tukey hsd post-hoc tests

tukey_results <- data.frame()
tukey_group_letters <- data.frame()

for (metric in storage_metrics) {
  anova_row <- anova_results %>% filter(metric == !!metric)
  if (nrow(anova_row) == 0 || anova_row$p_value >= 0.05) next

  # prepare data
  anova_data <- metric_data_long %>%
    filter(metric == !!metric) %>%
    select(site, value)

  # run anova and tukey hsd
  anova_fit <- aov(value ~ site, data = anova_data)
  tukey_fit <- TukeyHSD(anova_fit)

  # extract tukey results
  tukey_df <- as.data.frame(tukey_fit$site)
  tukey_df$comparison <- rownames(tukey_df)
  tukey_df$metric <- metric

  tukey_results <- rbind(tukey_results, tukey_df %>%
    select(metric, comparison, diff, lwr, upr, p_adj = `p adj`))

  # grouping letters by site for compact significance display
  pvals <- setNames(tukey_df$`p adj`, tukey_df$comparison)
  letters <- multcompLetters(pvals)$Letters
  letters_df <- data.frame(
    metric = metric,
    site = names(letters),
    group_letter = as.character(letters),
    stringsAsFactors = FALSE
  )
  tukey_group_letters <- rbind(tukey_group_letters, letters_df)
}

# save tukey hsd results
if (nrow(tukey_results) > 0) {
  tukey_results <- tukey_results %>%
    arrange(factor(metric, levels = storage_metrics), comparison)
  write.csv(tukey_results,
            file.path(output_dir, "tukey_hsd_results.csv"),
            row.names = FALSE)
}

if (nrow(tukey_group_letters) > 0) {
  tukey_group_letters <- tukey_group_letters %>%
    mutate(site = factor(site, levels = site_order)) %>%
    arrange(factor(metric, levels = storage_metrics), site) %>%
    mutate(site = as.character(site))
  write.csv(tukey_group_letters,
            file.path(output_dir, "tukey_group_letters.csv"),
            row.names = FALSE)
}

# summary table: site means

site_means <- HJA_annual %>%
  mutate(site = factor(site, levels = site_order)) %>%
  select(site) %>%
  distinct() %>%
  tidyr::crossing(metric = storage_metrics) %>%
  left_join(
    metric_data_long %>%
      group_by(site, metric) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("site", "metric")
  ) %>%
  tidyr::pivot_wider(
    names_from = metric,
    values_from = c(mean, sd),
    names_glue = "{metric}_{.value}"
  ) %>%
  # keep legacy column order style: metric_mean then metric_sd.
  select(
    site,
    dplyr::all_of(unlist(lapply(storage_metrics, function(m) c(paste0(m, "_mean"), paste0(m, "_sd")))))
  ) %>%
  mutate(site = factor(site, levels = site_order)) %>%
  arrange(site) %>%
  mutate(site = as.character(site))

write.csv(site_means,
          file.path(output_dir, "site_means_storage_metrics.csv"),
          row.names = FALSE)

# summary stats table: long format per site and metric
site_summary_stats <- HJA_annual %>%
  mutate(site = factor(site, levels = site_order)) %>%
  select(site) %>%
  distinct() %>%
  tidyr::crossing(metric = storage_metrics) %>%
  left_join(
    metric_data_long %>%
      group_by(site, metric) %>%
      summarise(
        n = sum(is.finite(value)),
        mean = ifelse(n > 0, mean(value, na.rm = TRUE), NA_real_),
        sd = ifelse(n > 1, sd(value, na.rm = TRUE), NA_real_),
        min = ifelse(n > 0, min(value, na.rm = TRUE), NA_real_),
        median = ifelse(n > 0, median(value, na.rm = TRUE), NA_real_),
        max = ifelse(n > 0, max(value, na.rm = TRUE), NA_real_),
        .groups = "drop"
      ),
    by = c("site", "metric")
  ) %>%
  mutate(
    n = ifelse(is.na(n), 0L, as.integer(n))
  ) %>%
  mutate(site = factor(site, levels = site_order)) %>%
  arrange(site, factor(metric, levels = storage_metrics)) %>%
  mutate(site = as.character(site))

write.csv(
  site_summary_stats,
  file.path(output_dir, "storage_metrics_summary_stats_by_site.csv"),
  row.names = FALSE
)
