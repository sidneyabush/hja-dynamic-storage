# =============================================================================
# Main Text and Supplement Figures - HJA Dynamic Storage
# =============================================================================
#
# Main Text Figures:
#   Figure 3: Dynamic Storage (RBI, RCS, FDC, SD) - faceted
#   Figure 4: Mobile Isotope (MTT, Fyw, DR) - faceted
#   Figure 5: CHS (baseflow fraction) - single panel
#
# Supplement Figures:
#   - Faceted time series for Dynamic Storage metrics
#   - Faceted time series for Mobile Storage metrics
#
# Author: Sidney Bush
# Date: 2026-02-02
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(colorspace)

rm(list = ls())

# =============================================================================
# SOURCE CONFIGURATION
# =============================================================================

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

config_path <- file.path(dirname(script_dir), "config.R")
if (!file.exists(config_path)) {
  config_path <- file.path(getwd(), "config.R")
}
source(config_path)

# =============================================================================
# SETUP
# =============================================================================

base_dir <- BASE_DATA_DIR
main_dir <- file.path(FIGURES_DIR, "Main_Text")
supp_dir <- file.path(FIGURES_DIR, "Supplement")

for (d in c(main_dir, supp_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Site order using WS## format (not GSWS##)
site_order <- c("WS09", "WS10", "WS01", "Look", "WS02",
                "WS03", "WS06", "WS07", "WS08", "Mack")

# Site colors
site_colors <- c(
  "WS09" = "#882255", "WS10" = "#AA4499", "WS01" = "#CC6677",
  "Look" = "#DDCC77", "WS02" = "#999933", "WS03" = "#117733",
  "WS06" = "#44AA99", "WS07" = "#88CCEE", "WS08" = "#6699CC",
  "Mack" = "#332288"
)

# Metric labels
metric_labels <- c(
  "RBI" = "RBI", "RCS" = "RCS", "FDC" = "FDC", "SD" = "SD (mm)",
  "WB" = "WB (mm)", "CHS" = "CHS", "MTT" = "MTT (yr)",
  "Fyw" = "Fyw", "DR" = "DR"
)

# Publication theme
theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, hjust = 0),
    legend.position = "none"
  )

theme_set(theme_pub)

# =============================================================================
# HELPER: Standardize site names to WS## format
# =============================================================================

standardize_sites <- function(df) {
  df %>%
    mutate(
      site = case_when(
        site == "GSWS01" ~ "WS01",
        site == "GSWS02" ~ "WS02",
        site == "GSWS03" ~ "WS03",
        site == "GSWS06" ~ "WS06",
        site == "GSWS07" ~ "WS07",
        site == "GSWS08" ~ "WS08",
        site == "GSWS09" ~ "WS09",
        site == "GSWS10" ~ "WS10",
        site %in% c("GSWSMC", "GSMACK") ~ "Mack",
        site %in% c("GSLOOK", "GSLOOK_FULL", "Look_FULL") ~ "Look",
        TRUE ~ site
      )
    ) %>%
    filter(site %in% site_order) %>%
    mutate(site = factor(site, levels = site_order))
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

# Annual storage metrics
annual_file <- file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv")
if (!file.exists(annual_file)) {
  annual_file <- file.path(OUTPUT_DIR, "HJA_StorageMetrics_Annual_All.csv")
}

annual_data <- read_csv(annual_file, show_col_types = FALSE) %>%
  rename_with(~ case_when(
    .x == "recession_curve_slope" ~ "RCS",
    .x == "fdc_slope" ~ "FDC",
    .x == "S_annual_mm" ~ "SD",
    .x == "mean_bf" ~ "CHS",
    .x == "DS_sum" ~ "WB",
    TRUE ~ .x
  )) %>%
  standardize_sites()

cat("  Loaded annual data:", nrow(annual_data), "rows\n")

# Isotope data for mobile storage (site-level metrics)
isotope_file <- file.path(OUTPUT_DIR, "Isotope_Metrics_Site.csv")
if (!file.exists(isotope_file)) {
  isotope_file <- file.path(base_dir, "Isotopes", "Isotope_Metrics_Site.csv")
}

if (file.exists(isotope_file)) {
  isotope_data <- read_csv(isotope_file, show_col_types = FALSE) %>%
    standardize_sites()
  cat("  Loaded isotope data:", nrow(isotope_data), "rows\n")
} else {
  isotope_data <- NULL
  cat("  Warning: Isotope data not found\n")
}

# CHS data
chs_file <- file.path(base_dir, "EC", "CHS_annual_baseflow.csv")
if (!file.exists(chs_file)) {
  chs_file <- file.path(OUTPUT_DIR, "CHS_annual_baseflow.csv")
}

if (file.exists(chs_file)) {
  chs_data <- read_csv(chs_file, show_col_types = FALSE) %>%
    standardize_sites()
  cat("  Loaded CHS data:", nrow(chs_data), "rows\n")
} else {
  # Try to get CHS from annual data
  if ("CHS" %in% names(annual_data)) {
    chs_data <- annual_data %>%
      select(site, year, CHS) %>%
      filter(!is.na(CHS))
    cat("  Using CHS from annual data:", nrow(chs_data), "rows\n")
  } else {
    chs_data <- NULL
    cat("  Warning: CHS data not found\n")
  }
}

# =============================================================================
# FIGURE 3: DYNAMIC STORAGE (RBI, RCS, FDC, SD)
# =============================================================================

cat("\nCreating Figure 3: Dynamic Storage...\n")

dynamic_metrics <- c("RBI", "RCS", "FDC", "SD")

dynamic_long <- annual_data %>%
  select(site, year, any_of(dynamic_metrics)) %>%
  pivot_longer(cols = -c(site, year), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric, levels = dynamic_metrics))

# Calculate mean ± SD per site
dynamic_summary <- dynamic_long %>%
  group_by(site, metric) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

# Add mean labels for each facet
dynamic_labels <- dynamic_summary %>%
  group_by(metric) %>%
  summarise(
    x_pos = -Inf,
    y_pos = -Inf,
    .groups = "drop"
  )

fig3 <- ggplot(dynamic_summary, aes(x = site, y = mean_val, color = site)) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.3, linewidth = 0.5
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2,
             labeller = labeller(metric = metric_labels)) +
  scale_color_manual(values = site_colors) +
  labs(x = NULL, y = "Mean ± 1 SD")

ggsave(file.path(main_dir, "Fig3_Dynamic_Storage.png"), fig3, width = 10, height = 8, dpi = 300)
ggsave(file.path(main_dir, "Fig3_Dynamic_Storage.pdf"), fig3, width = 10, height = 8)
cat("  Saved Figure 3\n")

# =============================================================================
# FIGURE 4: MOBILE ISOTOPE (MTT, Fyw, DR)
# =============================================================================

cat("Creating Figure 4: Mobile Isotope Storage...\n")

if (!is.null(isotope_data) && nrow(isotope_data) > 0) {
  # Create individual panels and combine with patchwork
  panels <- list()

  # Panel A: MTT
  if ("MTT" %in% names(isotope_data)) {
    mtt_data <- isotope_data %>% filter(!is.na(MTT))
    if (nrow(mtt_data) > 0) {
      panels$MTT <- ggplot(mtt_data, aes(x = site, y = MTT, color = site)) +
        geom_point(size = 3) +
        scale_color_manual(values = site_colors) +
        labs(x = NULL, y = "MTT (yr)")
    }
  }

  # Panel B: Fyw
  if ("Fyw" %in% names(isotope_data)) {
    fyw_data <- isotope_data %>% filter(!is.na(Fyw))
    if (nrow(fyw_data) > 0) {
      panels$Fyw <- ggplot(fyw_data, aes(x = site, y = Fyw, color = site)) +
        geom_point(size = 3) +
        scale_color_manual(values = site_colors) +
        labs(x = NULL, y = "Fyw")
    }
  }

  # Panel C: DR (with error bars if available)
  if ("DR" %in% names(isotope_data)) {
    dr_data <- isotope_data %>% filter(!is.na(DR))
    if (nrow(dr_data) > 0) {
      p_dr <- ggplot(dr_data, aes(x = site, y = DR, color = site)) +
        geom_point(size = 3) +
        scale_color_manual(values = site_colors) +
        labs(x = NULL, y = "DR")

      if ("DR_err" %in% names(isotope_data)) {
        p_dr <- p_dr +
          geom_errorbar(
            aes(ymin = DR - DR_err, ymax = DR + DR_err),
            width = 0.3, linewidth = 0.5
          )
      }
      panels$DR <- p_dr
    }
  }

  if (length(panels) >= 1) {
    fig4 <- wrap_plots(panels, ncol = 1)
    ggsave(file.path(main_dir, "Fig4_Mobile_Isotope.png"), fig4, width = 8, height = 3 * length(panels), dpi = 300)
    ggsave(file.path(main_dir, "Fig4_Mobile_Isotope.pdf"), fig4, width = 8, height = 3 * length(panels))
    cat("  Saved Figure 4\n")
  } else {
    cat("  Skipping Figure 4: No mobile isotope metrics available\n")
  }
} else {
  cat("  Skipping Figure 4: No isotope data\n")
}

# =============================================================================
# FIGURE 5: CHS (BASEFLOW FRACTION)
# =============================================================================

cat("Creating Figure 5: CHS...\n")

if (!is.null(chs_data) && nrow(chs_data) > 0) {
  # Summary stats
  chs_summary <- chs_data %>%
    group_by(site) %>%
    summarise(
      mean_val = mean(CHS, na.rm = TRUE),
      sd_val = sd(CHS, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(mean_val))

  fig5 <- ggplot(chs_summary, aes(x = site, y = mean_val, color = site)) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
      width = 0.3, linewidth = 0.6
    ) +
    scale_color_manual(values = site_colors) +
    labs(x = NULL, y = "Baseflow Fraction (CHS)")

  ggsave(file.path(main_dir, "Fig5_CHS.png"), fig5, width = 8, height = 5, dpi = 300)
  ggsave(file.path(main_dir, "Fig5_CHS.pdf"), fig5, width = 8, height = 5)
  cat("  Saved Figure 5\n")
} else {
  cat("  Skipping Figure 5: No CHS data\n")
}

# =============================================================================
# SUPPLEMENT: FACETED TIME SERIES - DYNAMIC STORAGE
# =============================================================================

cat("\nCreating Supplement: Dynamic Storage Time Series...\n")

# Calculate site means for labels
site_means_dynamic <- dynamic_long %>%
  group_by(site, metric) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

supp_dynamic_ts <- ggplot(dynamic_long, aes(x = year, y = value, color = site)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_text(
    data = site_means_dynamic,
    aes(x = -Inf, y = -Inf, label = sprintf("%.2f", mean_val)),
    hjust = -0.1, vjust = -0.5, size = 2.5, color = "black"
  ) +
  facet_grid(metric ~ site, scales = "free_y",
             labeller = labeller(metric = metric_labels)) +
  scale_color_manual(values = site_colors) +
  labs(x = "Water Year", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text.y = element_text(angle = 0, size = 9),
    strip.text.x = element_text(size = 8)
  )

ggsave(file.path(supp_dir, "FigS_Dynamic_TimeSeries.png"), supp_dynamic_ts, width = 14, height = 10, dpi = 300)
ggsave(file.path(supp_dir, "FigS_Dynamic_TimeSeries.pdf"), supp_dynamic_ts, width = 14, height = 10)
cat("  Saved Dynamic Storage Time Series\n")

# =============================================================================
# SUPPLEMENT: FACETED TIME SERIES - CHS
# =============================================================================

if (!is.null(chs_data) && nrow(chs_data) > 0) {
  cat("Creating Supplement: CHS Time Series...\n")

  chs_means <- chs_data %>%
    group_by(site) %>%
    summarise(mean_val = mean(CHS, na.rm = TRUE), .groups = "drop")

  supp_chs_ts <- ggplot(chs_data, aes(x = year, y = CHS, color = site)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.5) +
    geom_text(
      data = chs_means,
      aes(x = -Inf, y = -Inf, label = sprintf("%.2f", mean_val)),
      hjust = -0.1, vjust = -0.5, size = 3, color = "black"
    ) +
    facet_wrap(~site, ncol = 2) +
    scale_color_manual(values = site_colors) +
    labs(x = "Water Year", y = "Baseflow Fraction (CHS)")

  ggsave(file.path(supp_dir, "FigS_CHS_TimeSeries.png"), supp_chs_ts, width = 10, height = 10, dpi = 300)
  ggsave(file.path(supp_dir, "FigS_CHS_TimeSeries.pdf"), supp_chs_ts, width = 10, height = 10)
  cat("  Saved CHS Time Series\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== FIGURES COMPLETE ===\n")
cat("\nMain Text figures:", main_dir, "\n")
for (f in list.files(main_dir, pattern = "\\.(png|pdf)$")) {
  cat("  -", f, "\n")
}

cat("\nSupplement figures:", supp_dir, "\n")
for (f in list.files(supp_dir, pattern = "\\.(png|pdf)$")) {
  cat("  -", f, "\n")
}
