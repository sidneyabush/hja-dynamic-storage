# Goal is to develop correlations and relatinships among storage metrics and with catchement attributes.
# Inputs: base_dir/Catchment_Charc.csv; base_dir/mean_july.temp.csv; base_dir/mean_july_airtemp.csv; base_dir/HJA_StorageMetrics_Annual_All.csv; base_dir/DampingRatios_2025-07-07.csv.
# Author: Sidney Bush
# Date: 2026-02-13

library(dplyr)
library(readr)
library(ggplot2)
library(colorspace)
library(tidyr)
library(multcompView) 
library(ggcorrplot)

base_dir   <-"/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/03_Data"
output_dir <-"/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/05_Outputs/Hydrometric"

#reading in Catchment attributes from table that Zach made 
Catt<-read_csv(
  file.path(base_dir, "DynamicStorage", "Catchment_Charc.csv"),
  show_col_types = FALSE
)
Catt <- Catt %>%
  rename(
    site = Site
  )
#read in temperature data (They don't have years i the files)
#Mean stream water temp in July 
WtempM<-read_csv(
  file.path(base_dir, "Stream_T/Output", "mean_july.temp.csv"),
  show_col_types = FALSE
)
WtempM <- WtempM %>%
  rename(
    year = Year
  )
#reformatting 
WtempML<-pivot_longer(WtempM,
                      cols = -year,
                      names_to = "site",
                      values_to = "JulyM_ST")
WtempML$site[WtempML$site == "GSMACK"] <- "GSWSMC"

#Mean air temp in July 
AtempM<-read_csv(
  file.path(base_dir, "Stream_T/Output", "mean_july_airtemp.csv"),
  show_col_types = FALSE
)
AtempM <- AtempM %>%
  rename(
    year = Year
  )
#reformatting 
AtempML<-pivot_longer(AtempM,
                      cols = -year,
                      names_to = "site",
                      values_to = "JulyM_AT")
AtempML$site[AtempML$site == "GSMACK"] <- "GSWSMC"

HJA_Temp<-left_join(WtempML, AtempML, by = c("site", "year"))
HJA_Temp$JST_AT<-HJA_Temp$JulyM_ST/HJA_Temp$JulyM_AT


#reading in all storage metric data at the annual scale 
HJA_storage <- read_csv(
  file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual_All.csv"),
  show_col_types = FALSE,
  col_select = -1  
)

HJA_Stor_Temp_Yr<-left_join(HJA_storage, HJA_Temp, by = c("site", "year"))

#saving the yearly data for these sites 
write.csv(HJA_Stor_Temp_Yr,file.path(output_dir, "HJA_Stor_Temp_Yr.csv"))

#read the damping ratios are only available as a singular value, how we calculate this could be debated
Damp_R<-read_csv(
  file.path(base_dir, "Isotopes", "DampingRatios_2025-07-07.csv"),
  show_col_types = FALSE
)
Damp_R$site[Damp_R$site == "GSMACK"] <- "GSWSMC"

# Creating an average of the annual data 
cols <- c("recession_curve_slope", "RBI","Q5norm" , "CV_Q5norm", "mean_bf", "fdc_slope", "S_annual_mm", "DS_sum","JulyM_ST", "JulyM_AT","JST_AT")

HJA_storage_ave <- HJA_Stor_Temp_Yr %>%
  group_by(site) %>%
  summarise(across(all_of(cols), 
                   list(mean = ~mean(.x, na.rm = TRUE)
 ),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")

HJA_storage_ave<-left_join(HJA_storage_ave, Damp_R, by = "site")
HJA_storage_ave<- left_join(HJA_storage_ave, Catt, by = "site")

#Saving the all overall site based values 
write.csv(HJA_storage_ave,file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"))

HJA_storage_ave<-HJA_Ave

# List of attributes to include in the correlation
vars <- c(
  "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean", "mean_bf_mean",
  "fdc_slope_mean", "S_annual_mm_mean", "DS_sum_mean","DR_Overall","JST_AT_mean", "Area_km2", "Elevation_mean_m",
  "Slope_mean", "Aspec_Mean_deg", "Harvest", "Landslide_Young", "Landslide_Mod",
  "Landslide_Old",  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per"
)

vars1 <- c(
  "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean", "mean_bf_mean",
  "fdc_slope_mean", "S_annual_mm_mean", "DS_sum_mean",
  "DR_Overall", "JST_AT_mean")

vars2 <- c(
  "Area_km2", "Elevation_mean_m",
  "Slope_mean", "Aspec_Mean_deg", "Harvest", "Landslide_Young", "Landslide_Mod",
  "Landslide_Old", "Landslide_Total",  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per")

cor_storage1<-cor(HJA_storage_ave[vars2],use = "pairwise.complete.obs")


pdf("HJA_Storage_CorrPlot_CatchM.pdf", width = 8, height = 8)

ggcorrplot(cor_storage1, 
           hc.order = FALSE,
           type = "lower",
           outline.col = "white", lab = T)

dev.off()

cor_storage<-cor(HJA_storage_ave[vars],use = "pairwise.complete.obs")

pdf("HJA_Storage_CorrPlot.pdf", width = 12, height = 12)

ggcorrplot(cor_storage, 
           hc.order = FALSE,
           type = "lower",
           outline.col = "white", lab = T)

dev.off()

pdf("HJA_Storage_Metrics_Avg.pdf", width = 14, height = 14)

ggpairs(
  HJA_storage,
  lower = list(continuous = no_se_smoother),
  columns = 2:10
)+theme_bw()

dev.off()

