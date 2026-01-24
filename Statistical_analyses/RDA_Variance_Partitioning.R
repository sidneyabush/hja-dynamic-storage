# =============================================================================
# Redundancy Analysis (RDA) & Variance Partitioning
# =============================================================================
# Purpose: Use RDA to partition variance in storage metrics explained by
#          catchment characteristics
#
# Workflow:
#   1. Load site-averaged storage metrics and catchment attributes
#   2. Center and scale all variables
#   3. Run RDA: storage metrics ~ catchment characteristics
#   4. Extract variance explained by RDA axes
#   5. Create RDA biplot showing:
#      - Site scores (points)
#      - Storage metric loadings (response arrows)
#      - Catchment attribute loadings (predictor arrows)
#
# Inputs:
#   - HJA_Ave_StorageMetrics_CatCharacter.csv (from Correlations_Metrics.R)
#
# Outputs:
#   - RDA_biplot.png: RDA1 vs RDA2 with loadings
#   - RDA_variance_explained.csv: Proportion of variance explained
#
# Author: Based on Keira Johnson SiSyn code, adapted by Sidney Bush
# Date: 2026-01-23
# =============================================================================

# Load libraries
library(vegan)       # for rda()
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(reshape2)
library(readr)

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
# 2. LOAD SITE-AVERAGED DATA
# =============================================================================

HJA_Ave <- read_csv(
  file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
  show_col_types = FALSE
)

# =============================================================================
# 3. SELECT RESPONSE VARIABLES (STORAGE METRICS)
# =============================================================================

response_vars <- c(
  "recession_curve_slope_mean",
  "RBI_mean",
  "Q5norm_mean",
  "fdc_slope_mean",
  "S_annual_mm_mean",
  "DS_sum_mean"
)

# =============================================================================
# 4. SELECT EXPLANATORY VARIABLES (CATCHMENT CHARACTERISTICS)
# =============================================================================

explanatory_vars <- c(
  "Area_km2",
  "Elevation_mean_m",
  "Slope_mean",
  "Harvest",
  "Landslide_Young",
  "Landslide_Total",
  "Lava1_per",
  "Lava2_per",
  "Ash_Per",
  "Pyro_per"
)

# =============================================================================
# 5. PREPARE DATA FOR RDA
# =============================================================================

# Remove rows with missing data in response or explanatory variables
rda_data <- HJA_Ave %>%
  select(site, all_of(c(response_vars, explanatory_vars))) %>%
  na.omit()

# Check if enough sites remain
if (nrow(rda_data) < 5) {
  stop("Insufficient data: fewer than 5 sites with complete data")
}

# Center and scale response variables
response_scaled <- rda_data %>%
  select(all_of(response_vars)) %>%
  mutate(across(everything(), scale))

# Center and scale explanatory variables
explanatory_scaled <- rda_data %>%
  select(all_of(explanatory_vars)) %>%
  mutate(across(everything(), scale))

# =============================================================================
# 6. RUN RDA
# =============================================================================

# RDA: storage metrics (response) ~ catchment characteristics (explanatory)
my_rda <- rda(response_scaled, explanatory_scaled)

# Summarize RDA
rda_sum <- summary(my_rda)

# Extract variance explained
rda_variance <- data.frame(
  Component = c("Total Variance", "Constrained (Explained)", "Unconstrained (Residual)"),
  Variance = c(
    sum(my_rda$CCA$eig) + sum(my_rda$CA$eig),
    sum(my_rda$CCA$eig),
    sum(my_rda$CA$eig)
  )
) %>%
  mutate(Proportion = Variance / sum(Variance))

# Save variance partitioning results
write.csv(rda_variance,
          file.path(output_dir, "RDA_variance_explained.csv"),
          row.names = FALSE)

# =============================================================================
# 7. EXTRACT RDA COMPONENTS FOR PLOTTING
# =============================================================================

# Site scores (points)
site_scores <- as.data.frame(rda_sum$sites[, 1:2])
site_scores$site <- rda_data$site

# Species scores (storage metrics - response variables)
species_scores <- as.data.frame(rda_sum$species[, 1:2]) * 2  # Scale for visibility
species_scores$variable <- rownames(species_scores)

# Biplot scores (catchment characteristics - explanatory variables)
biplot_scores <- as.data.frame(rda_sum$biplot[, 1:2])
biplot_scores$variable <- rownames(biplot_scores)

# =============================================================================
# 8. CREATE RDA BIPLOT
# =============================================================================

# Calculate axis labels with % variance explained
rda1_var <- round(100 * my_rda$CCA$eig[1] / sum(c(my_rda$CCA$eig, my_rda$CA$eig)), 1)
rda2_var <- round(100 * my_rda$CCA$eig[2] / sum(c(my_rda$CCA$eig, my_rda$CA$eig)), 1)

p_biplot <- ggplot() +
  # Site scores (points)
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2),
             size = 3, alpha = 0.7, color = "black") +
  geom_text_repel(data = site_scores, aes(x = RDA1, y = RDA2, label = site),
                  size = 3, color = "black") +

  # Storage metric arrows (response variables)
  geom_segment(data = species_scores,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 22.5, length = unit(0.3, "cm"), type = "closed"),
               linewidth = 0.6, color = "gray40") +
  geom_text_repel(data = species_scores,
                  aes(x = RDA1, y = RDA2, label = variable),
                  color = "gray40", size = 3) +

  # Catchment characteristic arrows (explanatory variables)
  geom_segment(data = biplot_scores,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 22.5, length = unit(0.3, "cm"), type = "closed"),
               linewidth = 0.6, color = "steelblue") +
  geom_text_repel(data = biplot_scores,
                  aes(x = RDA1, y = RDA2, label = variable),
                  color = "steelblue", size = 3) +

  labs(
    title = "RDA Biplot: Storage Metrics ~ Catchment Characteristics",
    subtitle = paste0("Constrained variance: ",
                      round(100 * rda_variance$Proportion[2], 1), "%"),
    x = paste0("RDA1 (", rda1_var, "%)"),
    y = paste0("RDA2 (", rda2_var, "%)")
  ) +

  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave(
  file.path(output_dir, "QA_RDA_biplot.png"),
  p_biplot, width = 12, height = 10, dpi = 300
)
