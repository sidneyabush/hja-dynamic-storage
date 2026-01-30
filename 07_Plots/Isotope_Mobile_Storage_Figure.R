# =============================================================================
# Mobile Storage Metrics Figure - Isotope-Based Metrics
# =============================================================================
# Purpose: Create 4-panel figure showing isotope-based mobile storage metrics:
#   Panel A: Damping Ratio (DR)
#   Panel B: Young Water Fraction (Fyw)
#   Panel C: Mean Transit Time - Method 1 (MTT1)
#   Panel D: Mean Transit Time - Method 2 (MTT2)
#
# Author: Sidney Bush
# Date: 2026-01-30
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)

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

base_dir <- BASE_DATA_DIR
output_dir <- file.path(FIGURES_DIR, "Publication")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Site colors from config
site_colors <- SITE_COLORS

# Extended site order for isotope sites (includes non-hydrometric sites)
isotope_site_order <- c(
  "GSWS09", "GSWS10", "GSWS01", "GSLOOK", "GSWS02",
  "GSWS03", "GSWS07", "GSWS08", "Mack",
  "MR", "NC", "LC", "LO2", "CC", "LO1"
)

# Extended color palette for additional isotope-only sites
extended_colors <- c(
  site_colors,
  "MR"  = "#666666",   # McRae Creek - gray
  "NC"  = "#999999",   # North Creek - light gray
  "LC"  = "#AAAAAA",   # Lookout Creek - lighter gray
  "LO2" = "#BBBBBB",   # LO2 Creek
  "CC"  = "#CCCCCC",   # CC Creek
  "LO1" = "#DDDDDD"    # LO1 Creek
)

# Publication theme
theme_pub <- theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

theme_set(theme_pub)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading isotope data...\n")

# Damping Ratios
dr_data <- read_csv(
  file.path(base_dir, "Isotopes", "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = case_when(
      site == "GSMACK" ~ "Mack",
      site == "MCRAEC" ~ "MR",
      site == "NCCREC" ~ "NC",
      site == "LCCREC" ~ "LC",
      site == "LO2CRE" ~ "LO2",
      site == "CCCREE" ~ "CC",
      site == "LO1CRE" ~ "LO1",
      TRUE ~ site
    )
  ) %>%
  filter(!is.na(DR_Overall)) %>%
  mutate(site = factor(site, levels = isotope_site_order))

# MTT and Fyw
mtt_fyw <- read_csv(
  file.path(base_dir, "Isotopes", "MTT_FYW.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    site = trimws(site),
    site = case_when(
      site == "GSWSMC" ~ "Mack",
      site == "MCRAEC" ~ "MR",
      site == "GSLOOK " ~ "GSLOOK",
      TRUE ~ site
    )
  ) %>%
  filter(!is.na(site), site != "") %>%
  mutate(site = factor(site, levels = isotope_site_order))

# =============================================================================
# PANEL A: DAMPING RATIO
# =============================================================================

cat("Creating Panel A: Damping Ratio...\n")

p_dr <- ggplot(dr_data, aes(x = site, y = DR_Overall, color = site)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = DR_Overall - DR__err, ymax = DR_Overall + DR__err),
    width = 0.3, linewidth = 0.6
  ) +
  scale_color_manual(values = extended_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "A) Damping Ratio (DR)",
    x = NULL,
    y = "Damping Ratio"
  ) +
  theme(plot.margin = margin(5, 10, 5, 5))

# =============================================================================
# PANEL B: YOUNG WATER FRACTION
# =============================================================================

cat("Creating Panel B: Young Water Fraction...\n")

fyw_data <- mtt_fyw %>%
  filter(!is.na(FYWM)) %>%
  mutate(
    Fyw_low = ifelse(is.na(FYWL), FYWM, FYWL),
    Fyw_high = ifelse(is.na(FYWH), FYWM, FYWH)
  )

p_fyw <- ggplot(fyw_data, aes(x = site, y = FYWM, color = site)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = Fyw_low, ymax = Fyw_high),
    width = 0.3, linewidth = 0.6
  ) +
  scale_color_manual(values = extended_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "B) Young Water Fraction (Fyw)",
    x = NULL,
    y = "Fraction"
  ) +
  theme(plot.margin = margin(5, 10, 5, 5))

# =============================================================================
# PANEL C: MTT1 (McGuire method)
# =============================================================================

cat("Creating Panel C: MTT1...\n")

mtt1_data <- mtt_fyw %>%
  filter(!is.na(MTT1))

p_mtt1 <- ggplot(mtt1_data, aes(x = site, y = MTT1, color = site)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = MTT1 - MTT1_SD, ymax = MTT1 + MTT1_SD),
    width = 0.3, linewidth = 0.6
  ) +
  scale_color_manual(values = extended_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "C) Mean Transit Time - Method 1 (MTT1)",
    x = NULL,
    y = "Years"
  ) +
  theme(plot.margin = margin(5, 10, 5, 5))

# =============================================================================
# PANEL D: MTT2 (Segura/Ortega method - show range)
# =============================================================================

cat("Creating Panel D: MTT2...\n")

mtt2_data <- mtt_fyw %>%
  filter(!is.na(MTT2L) | !is.na(MTT2H)) %>%
  mutate(
    MTT2_mid = (MTT2L + MTT2H) / 2,
    MTT2_low = MTT2L,
    MTT2_high = MTT2H
  )

p_mtt2 <- ggplot(mtt2_data, aes(x = site, y = MTT2_mid, color = site)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = MTT2_low, ymax = MTT2_high),
    width = 0.3, linewidth = 0.6
  ) +
  scale_color_manual(values = extended_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "D) Mean Transit Time - Method 2 (MTT2)",
    x = NULL,
    y = "Years"
  ) +
  theme(plot.margin = margin(5, 10, 5, 5))

# =============================================================================
# COMBINE PANELS
# =============================================================================

cat("Combining panels...\n")

fig_isotope <- (p_dr | p_fyw) / (p_mtt1 | p_mtt2) +
  plot_annotation(
    title = "Mobile Storage Metrics from Isotope Analysis",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

# Save figure
ggsave(
  file.path(output_dir, "Fig_Isotope_Mobile_Storage.png"),
  fig_isotope, width = 12, height = 10, dpi = 300
)
ggsave(
  file.path(output_dir, "Fig_Isotope_Mobile_Storage.pdf"),
  fig_isotope, width = 12, height = 10
)

cat("\nSaved: Fig_Isotope_Mobile_Storage.png/pdf\n")
cat("Output directory:", output_dir, "\n")
