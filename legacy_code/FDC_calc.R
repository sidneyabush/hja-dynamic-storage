# Legacy analysis script retained for archival/reference use.
# Inputs: HF00402_v14.csv.
# Author: Legacy HJA storage team
# Date: 2026-02-13

setwd("/Users/keirajohnson/Box Sync/05_Storage_Manuscript/03_Data/Q")

discharge<-read.csv("HF00402_v14.csv")

discharge<-discharge %>%
  dplyr::filter(WATERYEAR > 1997) %>%
  dplyr::filter(!c(SITECODE %in% c("GSWSMA", "GSWSMF")))

calculate_water_year_fdc_slopes <- function(
    dates, flows,
    prob_range = c(0,100),
    log_scale = TRUE,
    min_days = 360
) {
  library(dplyr)
  library(ggplot2)
  
  # Input check
  if (length(dates) != length(flows)) stop("dates and flows must be the same length")
  
  # Build main data frame
  df <- data.frame(Date = as.Date(dates), Flow = as.numeric(flows)) %>%
    filter(!is.na(Flow)) %>%
    mutate(
      Month = as.numeric(format(Date, "%m")),
      Year = as.numeric(format(Date, "%Y")),
      WaterYear = ifelse(Month >= 10, Year + 1, Year)
    )
  
  # Filter complete years
  complete_years <- df %>%
    group_by(WaterYear) %>%
    summarise(DaysAvailable = n(), .groups = "drop") %>%
    filter(DaysAvailable >= min_days) %>%
    pull(WaterYear)
  
  df <- df %>% filter(WaterYear %in% complete_years)
  
  # Compute FDCs and slopes
  fdc_df <- df %>%
    group_by(WaterYear) %>%
    arrange(desc(Flow), .by_group = TRUE) %>%
    mutate(
      Rank = row_number(),
      N = n(),
      ExceedanceProbability = 100 * Rank / (N + 1)
    ) %>%
    filter(ExceedanceProbability >= prob_range[1],
           ExceedanceProbability <= prob_range[2]) %>%
    ungroup()
  
  # Compute slope for each year using log(flow) ~ exceedance
  slopes_df <- fdc_df %>%
    group_by(WaterYear) %>%
    summarise(
      Slope = coef(lm(log10(Flow) ~ ExceedanceProbability))[2],
      .groups = "drop"
    )
  
  # Plot FDCs
  p <- ggplot(fdc_df, aes(x = ExceedanceProbability, y = Flow, color = factor(WaterYear))) +
    geom_line() +
    labs(
      title = paste0("Flow Duration Curves by Water Year for ", sites[i]),
      x = "Exceedance Probability (%)",
      y = "Flow",
      color = "Water Year"
    ) +
    theme_minimal()
  
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  
  print(p)
  
  return(slopes_df)
}

calculate_total_fdc_slopes <- function(
    dates, flows,
    prob_range = c(0,100),
    log_scale = TRUE,
    min_days = 360
) {
  library(dplyr)
  library(ggplot2)
  
  # Input check
  if (length(dates) != length(flows)) stop("dates and flows must be the same length")
  
  # Build main data frame
  df <- data.frame(Date = as.Date(dates), Flow = as.numeric(flows)) %>%
    filter(!is.na(Flow)) %>%
    mutate(
      Month = as.numeric(format(Date, "%m")),
      Year = as.numeric(format(Date, "%Y")),
      WaterYear = ifelse(Month >= 10, Year + 1, Year)
    )
  
  # Filter complete years
  complete_years <- df %>%
    group_by(WaterYear) %>%
    summarise(DaysAvailable = n(), .groups = "drop") %>%
    filter(DaysAvailable >= min_days) %>%
    pull(WaterYear)
  
  df <- df %>% filter(WaterYear %in% complete_years)
  
  # Compute FDCs and slopes
  fdc_df <- df %>%
    arrange(desc(Flow), .by_group = TRUE) %>%
    mutate(
      Rank = row_number(),
      N = n(),
      ExceedanceProbability = 100 * Rank / (N + 1)
    ) %>%
    filter(ExceedanceProbability >= prob_range[1],
           ExceedanceProbability <= prob_range[2]) %>%
    ungroup()
  
  # Compute slope for each year using log(flow) ~ exceedance
  slopes_df <- fdc_df %>%
    summarise(
      Slope = coef(lm(log10(Flow) ~ ExceedanceProbability))[2],
      .groups = "drop"
    )
  
  # Plot FDCs
  p <- ggplot(fdc_df, aes(x = ExceedanceProbability, y = Flow)) +
    geom_line() +
    labs(
      title = paste0("Flow Duration Curves for ", sites[i]),
      x = "Exceedance Probability (%)",
      y = "Flow"
    ) +
    theme_minimal()
  
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  
  print(p)
  
  return(slopes_df)
}

sites<-unique(discharge$SITECODE)

slopes_df_list<-list()
slopes_df_WY_list<-list()

pdf("FDC_HJA_Visualization.pdf", width = 8, height = 6)

for (i in 1:length(sites)) {
  
  df<-subset(discharge, discharge$SITECODE==sites[i])
  df$Date<-as.Date(df$DATE, "%m/%d/%Y")
  df$Q<-df$MEAN_Q*0.02831683199881
  
  df<-subset(df, df$Q > 0)
  
  slopes_df_WY<-calculate_water_year_fdc_slopes(df$Date, df$Q, prob_range = c(5,95))
  
  slopes_df<-calculate_total_fdc_slopes(df$Date, df$Q, prob_range = c(5,95))
  
  slopes_df_WY$site<-sites[i]
  slopes_df$site<-sites[i]
  
  slopes_df_WY_list[[i]]<-slopes_df_WY
  
  slopes_df_list[[i]]<-slopes_df
  
}

dev.off()

fdc_curves_overall<-bind_rows(slopes_df_list)

write.csv(fdc_curves_overall, "FDC_slopes_overall.csv")

fdc_curves_WY<-bind_rows(slopes_df_WY_list)

write.csv(fdc_curves_WY, "FDC_slopes_WY.csv")

