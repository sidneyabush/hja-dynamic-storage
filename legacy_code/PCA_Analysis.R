# pca analysis for hja storage paper.
# inputs: output_dir/hja_ave_storagemetrics_catcharacter.csv; output_dir/hja_stor_temp_yr.csv.
# author: sidney bush
# date: 2026-02-13

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

#features<-c(  "recession_curve_slope_mean", "rbi_mean", "q5norm_mean", "cv_q5norm_mean", "mean_bf_mean","fdc_slope_mean", "s_annual_mm_mean", "ds_sum_mean","dr_overall","jst_at_mean", "area_km2", "elevation_mean_m","slope_mean", "aspec_mean_deg", "harvest", "landslide_young", "landslide_mod","landslide_old",  "lava1_per", "lava2_per", "ash_per", "pyro_per")

#features<-c(  "recession_curve_slope_mean", "rbi_mean", "q5norm_mean", "cv_q5norm_mean","fdc_slope_mean", "s_annual_mm_mean", "ds_sum_mean", "area_km2", "elevation_mean_m","slope_mean",  "harvest", "landslide_young", "landslide_old",  "lava1_per", "lava2_per", "ash_per", "pyro_per")


#features<-c( "recession_curve_slope", "rbi", "q5norm","mean_bf", "fdc_slope", "s_annual_mm", "ds_sum","jst_at")

features<-c( "recession_curve_slope", "RBI", "Q5norm", "fdc_slope", "S_annual_mm", "DS_sum")

# add your site id column name here (e.g., "site_id")
site_column <- "site"


HJA_selected <- HJA_Yr %>%
  select(all_of(site_column), all_of(features)) %>%
  drop_na()


# remove outliers using z-score threshold
HJA_z <- HJA_selected %>%
  mutate(across(all_of(features), ~ (.-mean(.))/sd(.)))


HJA_clean <- HJA_selected %>%
  filter(if_all(all_of(features), ~ abs((. - mean(.)) / sd(.)) < 3))


# normalize the features

scaled_features <- HJA_clean %>%
  select(all_of(site_column), all_of(features)) %>%
  mutate(across(all_of(features), scale))



# run pca
pca_result <- prcomp(scaled_features %>% select(-site), center = TRUE, scale. = TRUE)

summary(pca_result)

# create a data frame with pca scores and site names
pca_df <- as.data.frame(pca_result$x[, 1:2])  # PC1 and PC2
pca_df$site <- scaled_features$site     # Add site names back
#pca_df$q5norm <- hja_clean$q5norm


# extract loadings (rotation matrix)
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$feature <- rownames(loadings)

# scale loadings for visualization
loadings_scaled <- loadings %>%
  mutate(PC1 = PC1 * 5,  # adjust multiplier for visibility
         PC2 = PC2 * 5)


#site_order <- c("gsws10", "gsws09", "gsws01", "gslook", "gsws02","gsws03", "gswsmc", "gsws06", "gsws07", "gsws08")


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

# ensure site is a factor with the correct order
pca_df$site <- factor(pca_df$site, levels = names(site_colors))
                                                  

#pca_df$site <- factor(pca_df$site, levels = site_order)

# plot pc1 vs pc2 colored by percent lithology
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



# calculate variance explained
explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# create a data frame
explained_df <- data.frame(
  PC = paste0("PC", 1:9),
  Variance_Explained = explained_var[1:9]
)

# optional: plot

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

