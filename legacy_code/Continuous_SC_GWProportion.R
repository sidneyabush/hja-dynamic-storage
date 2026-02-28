# legacy analysis script retained for archival/reference use.
# inputs: hf00402_v14.csv.
# author: legacy hja storage team
# date: 2026-02-13

require(EflowStats)

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/03_Data/Q")

discharge<-read.csv("HF00402_v14.csv")
discharge$date<-as.Date(discharge$DATE, "%m/%d/%Y")

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/03_Data/EC")

EC<-read.delim("CF01201_v3.txt", sep = ",")

EC_daily<-EC %>%
  dplyr::mutate(date=as.Date(DATE_TIME, "%Y-%m-%d")) %>%
  dplyr::mutate(SITECODE=case_when(
    SITECODE=="GSMACK"~"GSWSMC",
    .default = SITECODE
  )) %>%
  dplyr::group_by(SITECODE, date) %>%
  dplyr::summarise(daily_SC=mean(EC_INST, na.rm = T))

EC_Q<-left_join(EC_daily, discharge[,c(3,6,12)])

ggplot(EC_Q, aes(date))+geom_line(mapping = aes(y=MEAN_Q), col="blue")+
  geom_line(mapping=aes(y=daily_SC))+
  facet_wrap(~SITECODE, scales = "free")


goodyears<-EC_Q %>%
  dplyr::mutate(waterYear=get_waterYear(date)) %>%
  dplyr::group_by(SITECODE, waterYear) %>%
  dplyr::summarise(num_days=n_distinct(date), .groups = 'drop') %>%
  dplyr::filter(num_days >= 365)

EC_Q <- EC_Q %>%
  mutate(waterYear = get_waterYear(date)) %>%
  semi_join(goodyears, by = c("SITECODE", "waterYear"))

EC_Q <- EC_Q %>%
  group_by(SITECODE) %>%
  mutate(runoff = quantile(daily_SC, 0.01, na.rm = T), groundwater = quantile(daily_SC, 0.99, na.rm = T))

EC_Q <- EC_Q %>%
  dplyr::mutate(GW = MEAN_Q*(daily_SC-runoff)/(groundwater-runoff), GW_prop=GW/MEAN_Q)

pdf("Continuous_Baseflow_Prop.pdf", width = 10, height = 7)

ggplot(EC_Q, aes(date, GW_prop))+geom_line()+facet_wrap(~SITECODE)+
  ylim(-0.1,1.1)+geom_hline(yintercept = 0, col="red")+
  geom_hline(yintercept = 1, col="red")+
  theme_bw()+
  labs(y="Mean Groundwater Proportion", y="Date")+theme(text = element_text(size = 20))

dev.off()

pdf("WY_Mean_Baseflow_Prop.pdf", width = 10, height = 6)

EC_Q %>%
  dplyr::group_by(SITECODE, waterYear) %>%
  summarise(mean_bf=mean(GW_prop, na.rm = T)) %>%
  ggplot(aes(mean_bf, SITECODE, col=waterYear))+geom_point(size=2)+theme_bw()+
  labs(x="Mean Groundwater Proportion", y="")+theme(text = element_text(size = 20))

dev.off()  
  
annual_bf_prop<-EC_Q %>%
  dplyr::group_by(SITECODE, waterYear) %>%
  summarise(mean_bf=mean(GW_prop, na.rm = T))

write.csv(annual_bf_prop, "Annual_GW_Prop.csv")
