# PCA analysis for HJA Storage paper.
# Inputs: output_dir/HJA_Ave_StorageMetrics_CatCharacter.csv; output_dir/HJA_Stor_Temp_Yr.csv.
# Author: Sidney Bush
# Date: 2026-02-13

library(tidyverse)
library(ggplot2)
library(scales)

library(ggrepel)


base_dir   <-"/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/03_Data"

output_dir <-"/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/05_Outputs/Hydrometric"

HJA_Ave<-read_csv(
  file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
  show_col_types = FALSE
)

HJA_Yr<-read_csv(
  file.path(output_dir, "HJA_Stor_Temp_Yr.csv"),
  show_col_types = FALSE
)

#features<-c(  "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean", "mean_bf_mean","fdc_slope_mean", "S_annual_mm_mean", "DS_sum_mean","DR_Overall","JST_AT_mean", "Area_km2", "Elevation_mean_m","Slope_mean", "Aspec_Mean_deg", "Harvest", "Landslide_Young", "Landslide_Mod","Landslide_Old",  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")

#features<-c(  "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean","fdc_slope_mean", "S_annual_mm_mean", "DS_sum_mean", "Area_km2", "Elevation_mean_m","Slope_mean",  "Harvest", "Landslide_Young", "Landslide_Old",  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")


#features<-c( "recession_curve_slope", "RBI", "Q5norm","mean_bf", "fdc_slope", "S_annual_mm", "DS_sum","JST_AT")

features<-c( "recession_curve_slope", "RBI", "Q5norm", "fdc_slope", "S_annual_mm", "DS_sum")

# Add your site ID column name here (e.g., "Site_ID")
site_column <- "site"


HJA_selected <- HJA_Yr %>%
  select(all_of(site_column), all_of(features)) %>%
  drop_na()


# Remove outliers using z-score threshold
HJA_z <- HJA_selected %>%
  mutate(across(all_of(features), ~ (.-mean(.))/sd(.)))


HJA_clean <- HJA_selected %>%
  filter(if_all(all_of(features), ~ abs((. - mean(.)) / sd(.)) < 3))


# Normalize the features

scaled_features <- HJA_clean %>%
  select(all_of(site_column), all_of(features)) %>%
  mutate(across(all_of(features), scale))



# Run PCA
pca_result <- prcomp(scaled_features %>% select(-site), center = TRUE, scale. = TRUE)

summary(pca_result)

# Create a data frame with PCA scores and site names
pca_df <- as.data.frame(pca_result$x[, 1:2])  # PC1 and PC2
pca_df$site <- scaled_features$site     # Add site names back
#pca_df$Q5norm <- HJA_clean$Q5norm


# Extract loadings (rotation matrix)
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$feature <- rownames(loadings)

# Scale loadings for visualization
loadings_scaled <- loadings %>%
  mutate(PC1 = PC1 * 5,  # adjust multiplier for visibility
         PC2 = PC2 * 5)


#site_order <- c("GSWS10", "GSWS09", "GSWS01", "GSLOOK", "GSWS02","GSWS03", "GSWSMC", "GSWS06", "GSWS07", "GSWS08")


site_colors <- c(
  "GSWS10" = "#AA4499",
  "GSWS09" = "#882255",
  "GSWS01" = "#CC6677",
  "GSLOOK" = "#DDCC77",
  "GSWS02" = "#999933",
  "GSWS03" = "#117733",
  "GSWSMC" = "#44AA99",
  "GSWS06" = "#88CCEE",
  "GSWS07" = "#6699CC",
  "GSWS08" = "#332288"
)

# Ensure site is a factor with the correct order
pca_df$site <- factor(pca_df$site, levels = names(site_colors))
                                                  

#pca_df$site <- factor(pca_df$site, levels = site_order)

# Plot PC1 vs PC2 colored by percent lithology
ggplot(pca_df, aes(x = PC1, y = PC2, color = site)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(name = "Site", values = site_colors) +
  
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30") +
  geom_text(data = loadings_scaled,
            aes(x = PC1, y = PC2, label = feature),
            color = "gray20", size = 2, vjust = 1.5) +
  
  theme_minimal() +
  #geom_text_repel(data = pca_df, aes(label = site), size = 3, color = "black") +
  
  
  theme(
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),  # remove minor grid lines
    axis.line = element_line(color = "black"),  # keep axis lines
    axis.ticks = element_line(color = "black")  # keep tick marks
  ) +
  
  labs(
    x = "Principal Component 1",
    y = "Principal Component 2"
  )



# Calculate variance explained
explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Create a data frame
explained_df <- data.frame(
  PC = paste0("PC", 1:9),
  Variance_Explained = explained_var[1:9]
)

# Optional: Plot

ggplot(explained_df, aes(x = PC, y = Variance_Explained)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = scales::percent(Variance_Explained, accuracy = 0.1)),
            vjust = -0.5, size = 3.5) +
  theme_minimal() +
  labs(
    title = "Variance Explained by Principal Components",
    y = "Proportion of Variance",
    x = "Principal Component"
  )

