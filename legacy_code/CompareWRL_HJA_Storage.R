# legacy analysis script retained for archival/reference use.
# inputs: /users/keirajohnson/box sync/hydrology_lab/projects/low flow/alldataharmonized_simple_05252025.csv; /users/keirajohnson/box sync/hydrology_lab/projects/low flow/rbfi.csv; /users/keirajohnson/box sync/05_storage_manuscript/03_data/storagemetrics/hja_storagemetrics_annual.csv.
# author: legacy hja storage team
# date: 2026-02-13

require(ggplot2)
require(dplyr)
require(ggpubr)

dat<-read.csv("/Users/keirajohnson/Box Sync/Hydrology_Lab/Projects/Low Flow/AllDataHarmonized_Simple_05252025.csv")

dat$Q5_norm<-(dat$Q5/(dat$drain_area_va*2.58999*1000000))*86400000

dat<-dat[complete.cases(dat),]

dat<-dat %>%
  dplyr::filter(!site_no %in% c(9153290, 13075983))

dat<-dat[,c(2,3,12,22,23)]

dat_avg<-dat %>%
  dplyr::group_by(site_no) %>%
  dplyr::summarise(bf=mean(mean_bf), recession=mean(recession_curve_slope), spcQ5=mean(Q5_norm),
            CV_spcQ5=sd(Q5_norm)/mean(Q5_norm))

rbfi<-read.csv("/Users/keirajohnson/Box Sync/Hydrology_Lab/Projects/Low Flow/RBFI.csv")

dat_avg<-left_join(dat_avg, rbfi[,c(2,3)])

HJA_storage<-read.csv("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/03_Data/StorageMetrics/HJA_StorageMetrics_Annual.csv")

HJA_storage$site<-factor(HJA_storage$site, levels = c("GSWS08","GSWS07", "GSWS06", "GSWSMC", "GSWS03",
                                                      "GSWS02", "GSWS01", "GSLOOK", "GSWS10", "GSWS09"))

HJA_storage_summary<-HJA_storage %>%
  group_by(site) %>%
  summarise(mean_rbfi=mean(rbfi), sd_rbfi=sd(rbfi), mean_Q5=mean(Q5norm), sd_Q5=sd(Q5norm),
            mean_mean_bf=mean(mean_bf, na.rm = T), sd_mean_bf=sd(mean_bf, na.rm=T), 
            mean_rcs=mean(recession_curve_slope),sd_rcs=sd(recession_curve_slope))

x_range <- range(dat_avg$rbfi, na.rm = TRUE)

p1<-ggplot(dat_avg, aes(x = rbfi)) +
  stat_density(geom = "raster", aes(y = 0, fill = ..density..), position = "identity") +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size = 20),
        legend.position = "left") +
  geom_rect(
    aes(xmin = x_range[1], xmax = x_range[2], ymin = -0.5, ymax = 0.5),
    fill = NA, color = "black", linewidth = 0.5
  ) +
  labs(x = "", fill = "Density")

p1

p2<-ggplot(HJA_storage_summary, aes(mean_rbfi, site))+
  geom_errorbar(aes(xmin=mean_rbfi-sd_rbfi, xmax=mean_rbfi+sd_rbfi))+
  geom_point()+theme_classic()+
  xlim(range(dat_avg$rbfi))+labs(x="RBFI", y="")+
  theme(text = element_text(size = 20))

p2

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/05_Outputs")

pdf("RBFI_compare.pdf", width = 10, height = 5.5)

ggarrange(p1, p2, nrow = 2, align = "hv", heights = c(0.3, 0.7))

dev.off()


x_range <- range(dat_avg$recession, na.rm = TRUE)

p1<-ggplot(dat_avg, aes(x = recession)) +
  stat_density(geom = "raster", aes(y = 0, fill = ..density..), position = "identity") +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size = 20),
        legend.position = "left") +
  geom_rect(
    aes(xmin = x_range[1], xmax = x_range[2], ymin = -0.5, ymax = 0.5),
    fill = NA, color = "black", linewidth = 0.5
  ) +
  labs(x = "Recession Curve Slope", fill = "Density")

p1

p2<-ggplot(HJA_storage_summary, aes(mean_rcs, site))+
  geom_errorbar(aes(xmin=mean_rcs-sd_rcs, xmax=mean_rcs+sd_rcs))+
  geom_point()+theme_classic()+
  xlim(range(dat_avg$recession))+labs(x="", y="")+
  theme(text = element_text(size = 20))

p2

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/05_Outputs")

pdf("RCS_compare.pdf", width = 10, height = 5.5)

ggarrange(p2, p1, nrow = 2, align = "hv")

dev.off()


x_range <- range(dat_avg$bf, na.rm = TRUE)

p1<-ggplot(dat_avg, aes(x = bf)) +
  stat_density(geom = "raster", aes(y = 0, fill = ..density..), position = "identity") +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size = 20),
        legend.position = "left") +
  geom_rect(
    aes(xmin = x_range[1], xmax = x_range[2], ymin = -0.5, ymax = 0.5),
    fill = NA, color = "black", linewidth = 0.5
  ) +
  labs(x = "Mean Baseflow Proportion", fill = "Density")

p1

p2<-ggplot(HJA_storage_summary, aes(mean_mean_bf, site))+
  geom_errorbar(aes(xmin=mean_mean_bf-sd_mean_bf, xmax=mean_mean_bf+sd_mean_bf))+
  geom_point()+theme_classic()+
  xlim(range(dat_avg$bf))+labs(x="", y="")+
  theme(text = element_text(size = 20))

p2

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/05_Outputs")

pdf("BF_prop_compare.pdf", width = 10, height = 5.5)

ggarrange(p2, p1, nrow = 2, align = "hv")

dev.off()


x_range <- range(dat_avg$spcQ5, na.rm = TRUE)

p1<-ggplot(dat_avg, aes(x = spcQ5)) +
  stat_density(geom = "raster", aes(y = 0, fill = ..density..), position = "identity") +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size = 20),
        legend.position = "left") +
  geom_rect(
    aes(xmin = x_range[1], xmax = x_range[2], ymin = -0.5, ymax = 0.5),
    fill = NA, color = "black", linewidth = 0.5
  ) +
  labs(x = "Mean Q5 (mm/day)", fill = "Density")

p1

p2<-ggplot(HJA_storage_summary, aes(mean_Q5, site))+
  geom_errorbar(aes(xmin=mean_Q5-sd_Q5, xmax=mean_Q5+sd_Q5))+
  geom_point()+theme_classic()+
  xlim(range(dat_avg$spcQ5))+labs(x="", y="")+
  theme(text = element_text(size = 20))

p2

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/05_Outputs")

pdf("Q5_compare.pdf", width = 10, height = 5.5)

ggarrange(p2, p1, nrow = 2, align = "hv")

dev.off()



