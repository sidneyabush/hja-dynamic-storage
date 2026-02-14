# Legacy analysis script retained for archival/reference use.
# Inputs: HJA_Ave_StorageMetrics_CatCharacter.csv; HJA_Stor_Temp_Yr.csv.
# Author: Legacy HJA storage team
# Date: 2026-02-13

require(GGally)

site_order <- c("GSWS10", "GSWS09", "GSWS01", "GSLOOK", "GSWS02","GSWS03", "GSWSMC", "GSWS06", "GSWS07", "GSWS08")
palette10 <- c(
  "#AA4499", "#882255", "#CC6677", "#DDCC77", "#999933",
  "#117733", "#44AA99", "#88CCEE", "#6699CC", "#332288"
)

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/05_Outputs/Hydrometric")

avg_met<-read.csv("HJA_Ave_StorageMetrics_CatCharacter.csv")

avg_met[3:38]<-lapply(avg_met[3:38], as.numeric)

keep_cols<-c("Area_km2", "Elevation_mean_m",
             "Slope_mean", "Aspec_Mean_deg", "Harvest", "Landslide_Young", "Landslide_Mod",
             "Landslide_Old",  "Lava1_per", "Lava2_per", "Ash_Per", "Pyro_per",
             "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean", "mean_bf_mean",
               "fdc_slope_mean", "S_annual_mm_mean", "DS_sum_mean","DR_Overall", "JST_AT_mean")

avg_met_filt<-avg_met %>%
  select(c(site, keep_cols))

avg_met_filt$site<-factor(avg_met_filt$site, levels = site_order)

yr_met<-read.csv("HJA_Stor_Temp_Yr.csv")

keep_cols<-c("recession_curve_slope", "RBI", "Q5norm", "mean_bf",
             "fdc_slope", "S_annual_mm", "DS_sum", "JST_AT")

keep_cols2<-c("recession_curve_slope", "RBI", "Q5norm",
             "fdc_slope", "S_annual_mm", "DS_sum", "JST_AT")

yr_met_filt<-yr_met %>%
  select(c(site,keep_cols))

yr_met_filt2<-yr_met %>%
  select(c(site, keep_cols2))

no_se_smoother <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_smooth(data=data, mapping=mapping, se = FALSE, method = "lm", color="black", ...) +
    geom_point()
}

pdf("Avg_Corr_Plot.pdf", width = 25, height = 25)

ggpairs(
  avg_met_filt,
  columns = 2:23,
  lower = list(continuous = no_se_smoother)
)+theme_bw()

dev.off()

pdf("Avg_Corr_Plot_Colored.pdf", width = 25, height = 25)

ggpairs(
  avg_met_filt,
  columns = 2:23,
  aes(col=site),
  upper = list(continuous = "blank"),
  lower = list(continuous = no_se_smoother)
  )+theme_bw()+theme(legend.position = "bottom")+
  scale_color_manual(values = palette10)

dev.off()

pdf("Yearly_Corr_Plot_Colored.pdf", width = 15, height = 15)

ggpairs(
  yr_met_filt,
  columns = 2:9,
  aes(col=site),
  upper = list(continuous = "blank"),
  lower = list(continuous = no_se_smoother))+theme_bw()+
  scale_color_manual(values = palette10)

dev.off()

site_list<-unique(yr_met_filt2$site)

pdf("Yearly_Corr_Plot_bySite.pdf", width = 15, height = 15)

for (i in 1:length(site_list)) {

  one_site<-subset(yr_met_filt2, yr_met_filt2$site==site_list[i])
  
  p1<-ggpairs(
    one_site,
    columns = 1:7,
    lower = list(continuous = no_se_smoother)
  )+theme_bw()+ggtitle(paste(site_list[i]))
    
  print(p1)
}

dev.off()


