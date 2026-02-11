avg_merequire(dplyr)
require(lubridate)
require(GGally)

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/03_Data/Q")

da_df<-read.csv("drainage_area.csv")

discharge<-read.csv("HF00402_v14.csv")
discharge$Date<-as.Date(discharge$DATE, "%m/%d/%Y")

discharge<-discharge %>%
  dplyr::filter(WATERYEAR > 1997)

discharge<-left_join(discharge, da_df)

discharge<-discharge[complete.cases(discharge$DA_M2),]

sites<-unique(discharge$SITECODE)
slope_list<-list()
rbfi_list<-list()
q5_list<-list()
RBFI_list_annual<-list()
slope_list_annual<-list()

for (i in 1:length(sites)) {
  
  df<-subset(discharge, discharge$SITECODE==sites[i])
  df$Date<-as.Date(df$DATE, "%m/%d/%Y")
  df$Q<-df$MEAN_Q*0.02831683199881
  
  water_years<-unique(df$WATERYEAR)
  
  df <- df %>%
    dplyr::arrange(Date) %>% #order by date
    dplyr::mutate(dQ = Q - lag(Q),  # Change in discharge
           change_dQ = Q/lag(Q), #find change in discharge relative to previous day
           dQ_dt = dQ / as.numeric(Date - lag(Date))) %>%  # Daily rate of change
    dplyr::filter(!is.na(dQ_dt)) %>% # Remove NA values (first row)
    dplyr::filter(!change_dQ < 0.7) #remove anything where the relative difference is < 0.7, essentially where the change in discharge is small
  
  # Calculate the recession slope (-dQ/dt)
  recession_data <- df %>%
    filter(dQ < 0) %>%  # Keep only recession periods
    mutate(recession_slope = -dQ_dt)  # Make it positive for the slope
  
  # Fit a linear model to the recession data
  lm_model <- lm(log(recession_slope) ~ log(Q), data = recession_data)
  model_summary <- summary(lm_model)
  
  # Extract statistics
  slope <- coef(lm_model)[2]
  p_value <- model_summary$coefficients[2, 4]
  r_squared <- model_summary$r.squared
  
  slope_list[[i]]<-data.frame(slope, sites[i])
  
  annual_slope<-list()
  
  for (k in 1:length(water_years)) {
    
    recession_data_annual<-recession_data %>%
      filter(WATERYEAR==water_years[k])
    
    # Fit a linear model to the recession data
    lm_model <- lm(log(recession_slope) ~ log(Q), data = recession_data_annual)
    model_summary <- summary(lm_model)
    
    # Extract statistics
    slope <- coef(lm_model)[2]
    p_value <- model_summary$coefficients[2, 4]
    r_squared <- model_summary$r.squared
    
    annual_slope[[k]]<-data.frame(slope, water_years[k], sites[i])
    
  }
  
  slope_list_annual[[i]]<-bind_rows(annual_slope)
  
  # Calculate daily changes in discharge
  df <- df %>%
    arrange(Date) %>%
    mutate(dQ = Q - lag(Q),  # Daily change in discharge
           abs_dQ = abs(dQ)) %>%              # Absolute change in discharge
    filter(!is.na(abs_dQ))  # Remove NA values (first row)
  
  # Calculate the total discharge over the period
  total_discharge <- sum(df$Q)
  
  # Calculate the Richards-Baker Flashiness Index
  RBFI <- sum(df$abs_dQ) / total_discharge
  
  rbfi_list[[i]]<-data.frame(RBFI, sites[i])
  
  annual_RBFI<-list()
  
  for (k in 1:length(water_years)) {
    
    df_annual<-df %>%
      filter(WATERYEAR==water_years[k])
    
    # Calculate the total discharge over the period
    total_discharge <- sum(df_annual$Q)
    
    # Calculate the Richards-Baker Flashiness Index
    RBFI <- sum(df_annual$abs_dQ) / total_discharge
    
    annual_RBFI[[k]]<-data.frame(RBFI, water_years[k], sites[i])
    
  }
  
  RBFI_list_annual[[i]]<-bind_rows(annual_RBFI)
  
  q5<-df %>%
    filter(month(Date) < 11 & month(Date) > 7) %>%
    group_by(year(as.Date(Date))) %>%
    summarise(lowflow=quantile((Q/DA_M2)*86400000, 0.05))
  
  q5$CV_lowflow<-sd(q5$lowflow)/mean(q5$lowflow)
  
  q5_list[[i]]<-data.frame(q5, sites[i])
  
}


recession_curve_slope<-do.call(bind_rows, slope_list)
colnames(recession_curve_slope)<-c("recession_curve_slope", "site")

recession_curve_slope_annual<-do.call(bind_rows, slope_list_annual)
colnames(recession_curve_slope_annual)<-c("recession_curve_slope", "year", "site")

RBFI_df<-do.call(bind_rows, rbfi_list)
colnames(RBFI_df)<-c("rbfi", "site")

RBFI_df_annual<-do.call(bind_rows, RBFI_list_annual)
colnames(RBFI_df_annual)<-c("rbfi", "year", "site")

Q5_df<-do.call(bind_rows, q5_list)
colnames(Q5_df)<-c("year", "Q5norm", "CV_Q5norm", "site")

Q5_df_average<-Q5_df %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mean_Q5=mean(Q5norm), mean_CVQ5=mean(CV_Q5norm))

HJA_storage<-left_join(recession_curve_slope, RBFI_df)

HJA_storage<-left_join(Q5_df_average, HJA_storage)

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/03_Data/StorageMetrics")

gw_prop<-read.csv("Annual_GW_Prop.csv")
gw_prop<-gw_prop[,c(2:4)]
colnames(gw_prop)<-c("site", "year", "mean_bf")

gw_prop_avg<-gw_prop %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mean_bf=mean(mean_bf))

HJA_storage<-left_join(HJA_storage, gw_prop_avg)

aspect<-read.csv("HJA_aspect.csv")
slope<-read.csv("HJA_slope.csv")
elevation<-read.csv("HJA_elevation.csv")

topo<-bind_rows(aspect, elevation, slope)
topo_t<-data.frame(t(topo))
colnames(topo_t)<-c("aspect", "elevation", "slope")
topo_t$site<-c("GSWS01", "GSWS10", "GSWS02", "GSWS03", "GSWS06", "GSWS07", "GSWS08",
               "GSWS09", "GSLOOK", "GSWSMC")

HJA_storage<-left_join(HJA_storage, topo_t)

fdc_avg<-read.csv("FDC_slopes_overall.csv")
fdc_avg<-fdc_avg[,c(2:3)]
colnames(fdc_avg)<-c("fdc_slope", "site")

HJA_storage<-left_join(HJA_storage, fdc_avg)

write.csv(HJA_storage, "HJA_StorageMetrics_Average.csv")

pdf("HJA_Storage_Metrics_Avg.pdf", width = 14, height = 14)

ggpairs(
  HJA_storage,
  lower = list(continuous = no_se_smoother),
  columns = 2:10
)+theme_bw()

dev.off()

HJA_complete<-HJA_storage[complete.cases(HJA_storage),]

cor_storage<-cor(HJA_complete[2:10])

pdf("HJA_Storage_CorrPlot.pdf", width = 8, height = 8)

ggcorrplot(cor_storage, hc.order = TRUE,type = "lower",
           outline.col = "white", lab = T)

dev.off()

HJA_storage<-left_join(recession_curve_slope_annual, RBFI_df_annual)

HJA_storage<-left_join(HJA_storage, Q5_df)

HJA_storage<-left_join(HJA_storage, gw_prop)

fdc<-read.csv("FDC_slopes_WY.csv")
fdc<-fdc[,c(2:4)]
colnames(fdc)<-c("year","fdc_slope", "site")

HJA_storage<-left_join(HJA_storage, fdc)

write.csv(HJA_storage, "HJA_StorageMetrics_Annual.csv")

no_se_smoother <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_smooth(se = FALSE, method = "lm", color="black", ...) +
    geom_point()
}

pdf("HJA_Storage_Metrics_Annual.pdf", width = 10, height = 10)

for (i in 1:length(sites)) {
  
  HJA_one_site<-HJA_storage %>%
    dplyr::filter(site==sites[i]) 
  
  if(i %in% c(1,8)){
   
    p1<-ggpairs(HJA_one_site[,c(1,4,5,8)],
                lower = list(continuous = no_se_smoother))+
      ggtitle(sites[i])+theme_bw() 
    
  }else{
    
    p1<-ggpairs(HJA_one_site[,c(1,4,5,7,8)],
                lower = list(continuous = no_se_smoother))+ggtitle(sites[i])+theme_bw() 
    
  }
  
  print(p1)
  
}

dev.off()

