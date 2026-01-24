# =============================================================================
# ANOVA & Tukey HSD Site Comparisons
# =============================================================================
# Purpose: Test for significant differences among sites in storage metrics
#          using one-way ANOVA followed by Tukey HSD post-hoc tests
#
# Workflow:
#   1. Load annual storage metrics (site-year observations)
#   2. For each storage metric, run one-way ANOVA (site as factor)
#   3. If ANOVA is significant (p < 0.05), run Tukey HSD post-hoc test
#   4. Visualize site comparisons with boxplots and Tukey groupings
#   5. Save results to CSV
#
# Inputs:
#   - HJA_StorageMetrics_Annual_All.csv (from 04_Aggregate_All_Storage_Metrics.R)
#
# Outputs:
#   - ANOVA_results.csv: ANOVA F-statistics and p-values
#   - Tukey_HSD_results.csv: Pairwise site comparisons with adjusted p-values
#   - QA_ANOVA_boxplots.png: Boxplots of storage metrics by site
#
# Author: Sidney Bush
# Date: 2026-01-23
# =============================================================================

# Load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(broom)

theme_set(theme_bw(base_size = 12))

# Clear environment
rm(list = ls())

# =============================================================================
# 1. SETUP: Directories
# =============================================================================

base_dir    <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data"
output_dir  <- "/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs"

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 2. LOAD ANNUAL STORAGE METRICS
# =============================================================================

HJA_annual <- read_csv(
  file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv"),
  show_col_types = FALSE
)

# =============================================================================
# 3. SELECT STORAGE METRICS FOR ANOVA
# =============================================================================

storage_metrics <- c(
  "recession_curve_slope",
  "RBI",
  "Q5norm",
  "CV_Q5norm",
  "mean_bf",
  "fdc_slope",
  "S_annual_mm",
  "DS_sum"
)

# =============================================================================
# 4. RUN ANOVA FOR EACH METRIC
# =============================================================================

anova_results <- data.frame()

for (metric in storage_metrics) {
  # Check if metric exists in data
  if (!metric %in% colnames(HJA_annual)) next

  # Prepare data for ANOVA
  anova_data <- HJA_annual %>%
    select(site, value = all_of(metric)) %>%
    filter(!is.na(value), !is.na(site))

  # Check if enough data
  if (nrow(anova_data) < 10) next

  # Run one-way ANOVA
  anova_fit <- aov(value ~ site, data = anova_data)
  anova_summary <- summary(anova_fit)

  # Extract F-statistic and p-value
  f_stat <- anova_summary[[1]]$`F value`[1]
  p_value <- anova_summary[[1]]$`Pr(>F)`[1]

  # Store results
  anova_results <- rbind(anova_results, data.frame(
    metric = metric,
    F_statistic = f_stat,
    p_value = p_value,
    df_between = anova_summary[[1]]$Df[1],
    df_within = anova_summary[[1]]$Df[2],
    significant = ifelse(p_value < 0.05, "Yes", "No")
  ))
}

# Save ANOVA results
write.csv(anova_results,
          file.path(output_dir, "ANOVA_results.csv"),
          row.names = FALSE)

# =============================================================================
# 5. RUN TUKEY HSD POST-HOC TESTS
# =============================================================================

tukey_results <- data.frame()

for (metric in storage_metrics) {
  # Check if metric exists and ANOVA was significant
  if (!metric %in% colnames(HJA_annual)) next

  anova_row <- anova_results %>% filter(metric == !!metric)
  if (nrow(anova_row) == 0 || anova_row$p_value >= 0.05) next

  # Prepare data
  anova_data <- HJA_annual %>%
    select(site, value = all_of(metric)) %>%
    filter(!is.na(value), !is.na(site))

  # Run ANOVA and Tukey HSD
  anova_fit <- aov(value ~ site, data = anova_data)
  tukey_fit <- TukeyHSD(anova_fit)

  # Extract Tukey results
  tukey_df <- as.data.frame(tukey_fit$site)
  tukey_df$comparison <- rownames(tukey_df)
  tukey_df$metric <- metric

  tukey_results <- rbind(tukey_results, tukey_df %>%
    select(metric, comparison, diff, lwr, upr, p_adj = `p adj`))
}

# Save Tukey HSD results
if (nrow(tukey_results) > 0) {
  write.csv(tukey_results,
            file.path(output_dir, "Tukey_HSD_results.csv"),
            row.names = FALSE)
}

# =============================================================================
# 6. QA PLOTS: BOXPLOTS BY SITE
# =============================================================================

# Prepare data for plotting
plot_data <- HJA_annual %>%
  select(site, all_of(storage_metrics)) %>%
  pivot_longer(
    cols = all_of(storage_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

# Create faceted boxplot
p_boxplots <- ggplot(plot_data, aes(x = site, y = value, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Storage Metrics by Site (ANOVA Comparisons)",
    x = "Site",
    y = "Value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  file.path(output_dir, "QA_ANOVA_boxplots.png"),
  p_boxplots, width = 12, height = 14, dpi = 300
)

# =============================================================================
# 7. SUMMARY TABLE: SITE MEANS
# =============================================================================

site_means <- HJA_annual %>%
  group_by(site) %>%
  summarise(
    across(
      all_of(storage_metrics),
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write.csv(site_means,
          file.path(output_dir, "Site_Means_Storage_Metrics.csv"),
          row.names = FALSE)
