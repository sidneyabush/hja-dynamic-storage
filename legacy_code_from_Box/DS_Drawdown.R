require(pracma)
require(dplyr)
require(EflowStats)

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/03_Data/Q")

discharge<-read.csv("HF00402_v14.csv")
discharge$date<-as.Date(discharge$DATE, "%m/%d/%Y")

goodyears<-discharge %>%
  dplyr::filter(!(SITECODE %in% c("GSWSMC", "GSWSMF"))) %>%
  dplyr::mutate(waterYear=get_waterYear(date)) %>%
  dplyr::group_by(SITECODE, waterYear) %>%
  dplyr::summarise(num_days=n_distinct(date), .groups = 'drop') %>%
  dplyr::filter(num_days >= 365)

discharge <- discharge %>%
  mutate(waterYear = get_waterYear(date)) %>%
  semi_join(goodyears, by = c("SITECODE", "waterYear"))

unique(discharge$SITECODE)

discharge<-discharge %>%
  dplyr::group_by(SITECODE) %>%
  dplyr::mutate(Q_smoothed = rollmean(MEAN_Q, k = 7, fill = NA, align = "right"))

discharge<-discharge %>%
  dplyr::group_by(waterYear, SITECODE) %>%
  dplyr::mutate(hydro_flux=cumsum(MEAN_Q), hydro_flux_prop=hydro_flux/max(hydro_flux))

discharge <- discharge %>%
  dplyr::group_by(waterYear, SITECODE) %>%
  dplyr::filter(hydro_flux_prop > .30) %>%
  dplyr::mutate(tag = if_else(Q_smoothed == max(Q_smoothed), "peak", NA_character_))

discharge$wyd<-get_waterYearDay(discharge$date)

# Define percentage threshold (e.g., 20% = 0.2)
threshold_pct <- 0.08

find_last_peak <- function(data, threshold_pct) {
  # Extract the time series for the current group
  time_series <- data$MEAN_Q
  
  # Define the threshold (percentage of maximum peak discharge)
  max_peak_discharge <- max(time_series, na.rm = TRUE)
  threshold_value <- max_peak_discharge * threshold_pct
  
  # Use `findpeaks` to identify peaks
  peaks <- tryCatch({
    findpeaks(time_series)
  }, error = function(e) {
    cat("Error in findpeaks:", e$message, "\n")
    return(NULL)
  })
  
  # If peaks are found, process them
  if (!is.null(peaks)) {
    peaks_df <- as_tibble(peaks) %>%
      rename(peak_height = V1, peak_index = V2) %>%
      filter(peak_height >= threshold_value) %>%
      mutate(
        date = data$date[peak_index],
        wyd = data$wyd[peak_index]
      ) %>%
      filter(wyd < 300) %>%
      arrange(peak_index)  # Ensure peaks are sorted by index
    
    # Get the last peak that satisfies the threshold condition
    last_valid_peak <- peaks_df %>%
      slice_tail(n = 1)  # Get the last row, which is the last peak
    
    if (nrow(last_valid_peak) > 0) {
      last_peak_index <- last_valid_peak$peak_index
      last_peak_value <- time_series[last_peak_index]
      last_peak_date <- data$date[last_peak_index]  # Get the date corresponding to the peak
      return(tibble(last_peak_date = last_peak_date, last_peak_value = last_peak_value))
    } else {
      return(tibble(last_peak_date = NA, last_peak_value = NA))
    }
  } else {
    return(tibble(last_peak_date = NA, last_peak_value = NA))
  }
}

discharge<-discharge[complete.cases(discharge$Q_smoothed),]

last_peak<-discharge %>%
  group_by(SITECODE, waterYear) %>%
  do(find_last_peak(., 0.08)) %>%
  ungroup()

last_peak$wyd<-get_waterYearDay(last_peak$last_peak_date)

write.csv(last_peak, "DS_Drawdown_Date.csv")

pdf("LastPeak_DSDrawdown.pdf", width = 8, height = 6)

ggplot(last_peak, aes(wyd, SITECODE))+geom_boxplot()+
  theme_bw()+theme(text = element_text(size = 20))+
  labs(x="Day of Water Year", y="")

dev.off()

years<-unique(last_peak$waterYear)

years<-years[years > 1997 & years < 2020]

pdf("DS_Drawdown_LastDay.pdf", width = 12, height = 10)

for (i in 1:length(years)) {
  
  print(years[i])
  
  dis_df<-subset(discharge, discharge$waterYear==years[i])
  last_peak_df<-subset(last_peak, last_peak$waterYear==years[i])
  
  p1<-ggplot()+geom_line(dis_df, mapping=aes(get_waterYearDay(date), MEAN_Q, group=waterYear), alpha=0.4)+
    geom_point(last_peak_df, mapping=aes(wyd, last_peak_value), col="blue")+
    theme_bw()+
    facet_wrap(~SITECODE, scales = "free_y")+ggtitle(years[i])+
    labs(x="Water Year Day", y="Mean Discharge (cfs)")
  
  print(p1)
}

dev.off()

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/03_Data/DynamicStorage")

DS_dat<-read.csv("daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv")

unique(DS_dat$SITECODE)

DS_dat <- DS_dat %>%
  mutate(SITECODE=case_when(
    SITECODE=="GSLOOK_FULL"~"GSLOOK",
    SITECODE=="GSMACK"~"GSWSMA",
    .default = SITECODE
  ))

DS_dat$DATE<-as.Date(DS_dat$DATE, "%m/%d/%y")

DS_dat$waterYear<-get_waterYear(as.Date(DS_dat$DATE))

DS_dat <- DS_dat %>%
  filter(waterYear > 1997 & waterYear < 2020)

DS_dat<-left_join(DS_dat, last_peak)
DS_dat<-DS_dat[complete.cases(DS_dat$last_peak_date),]

DS_dat_cropped<-DS_dat %>%
  group_by(SITECODE, waterYear) %>%
  filter(as.Date(DATE) >= last_peak_date) %>%
  mutate(max_q = Q_mm_d[last_peak_date==DATE])

DS_compute<-DS_dat_cropped %>%
  group_by(SITECODE, waterYear) %>%
  mutate(DS_daily=P_mm_d-Q_mm_d-ET_mm_d, DS_sum=cumsum(DS_daily))

pdf("DS_Drawdown_Visualization.pdf", width = 12, height = 8)

ggplot(DS_compute, aes(get_waterYearDay(DATE), DS_sum))+geom_line(aes(col=waterYear, group=waterYear))+
  facet_wrap(~SITECODE)+labs(x="Day of Water Year", y="Dynamic Storage Drawdown (mm)", col="Water Year")+
  theme_classic()+geom_hline(yintercept = 0)+theme(text = element_text(size = 20))+
  scale_color_gradient(low = "black", high = "grey88")

dev.off()

DS_max<-DS_compute %>%
  group_by(SITECODE, waterYear) %>%
  slice_min(DS_sum)

ggplot(DS_max, aes(SITECODE, DS_sum))+geom_boxplot()

DS_site_annual<-DS_max[,c(2,7,13)]

write.csv(DS_site_annual, "DS_drawdown_annual.csv")
