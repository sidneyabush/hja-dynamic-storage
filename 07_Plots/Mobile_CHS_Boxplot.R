# =============================================================================
# Mobile Storage Metrics Figure - CHS (Chemical Hydrograph Separation)
# =============================================================================
# Purpose: Create box plot showing annual baseflow fraction (CHS) by site
#
# CHS = Chemical Hydrograph Separation method
# Measures mean baseflow contribution to streamflow using specific conductance
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)

# Clear environment
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
if (file.exists(config_path)) {
  source(config_path)
} else {
  stop("config.R not found.")
}

# =============================================================================
# SETUP
# =============================================================================

output_dir <- file.path(FIGURES_DIR, "Publication")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Site colors and order from config
site_colors <- SITE_COLORS
site_order <- SITE_ORDER_HYDROMETRIC

# Publication theme
theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

theme_set(theme_pub)

# =============================================================================
# LOAD CHS DATA
# =============================================================================

cat("Loading CHS data...\n")

# Annual baseflow proportion from chemical hydrograph separation
chs_data <- read_csv(
  file.path(OUTPUT_DIR, "Annual_GW_Prop.csv"),
  show_col_types = FALSE
) %>%
  rename(site = SITECODE, year = waterYear) %>%
  # Handle both old (mean_bf) and new (CHS) column names
  rename_with(~ ifelse(.x == "mean_bf", "CHS", .x)) %>%
  mutate(
    site = case_when(
      site == "GSWSMC" ~ "Mack",
      site == "GSMACK" ~ "Mack",
      TRUE ~ site
    )
  ) %>%
  filter(site %in% site_order) %>%
  mutate(site = factor(site, levels = site_order))

cat("Sites with CHS data:", paste(unique(chs_data$site), collapse = ", "), "\n")
cat("Years per site:\n")
print(chs_data %>% group_by(site) %>% summarise(n_years = n(), .groups = "drop"))

# =============================================================================
# CREATE BOX PLOT
# =============================================================================

cat("\nCreating CHS box plot...\n")

p_chs <- ggplot(chs_data, aes(x = site, y = CHS, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  geom_jitter(aes(color = site), width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_color_manual(values = site_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "Baseflow Fraction from Chemical Hydrograph Separation (CHS)",
    x = NULL,
    y = "Baseflow Fraction"
  ) +
  ylim(0, 1)

# Save figure
ggsave(
  file.path(output_dir, "Mobile_CHS_Boxplot.png"),
  p_chs, width = 8, height = 5, dpi = 300
)
ggsave(
  file.path(output_dir, "Mobile_CHS_Boxplot.pdf"),
  p_chs, width = 8, height = 5
)

cat("\nSaved: Mobile_CHS_Boxplot.png/pdf\n")
cat("Output directory:", output_dir, "\n")
