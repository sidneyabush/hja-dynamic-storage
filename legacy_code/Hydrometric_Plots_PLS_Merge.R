# legacy analysis script retained for archival/reference use.
# inputs: base_dir/ds_drawdown_annual.csv; base_dir/storage_discharge_method_annual_vol_per_site_wateryear.csv; base_dir/hja_storagemetrics_annual.csv; base_dir/drainage_area.csv; base_dir/hf00402_v14.csv.
# author: legacy hja storage team
# date: 2026-02-13

library(dplyr)
library(readr)
library(ggplot2)
library(colorspace)
library(tidyr)
library(multcompView)   

theme_set(theme_classic(base_size = 14))

rm(list = ls())

#base_dir   <- "/users/sidneybush/library/cloudstorage/box-box/05_storage_manuscript/03_data"
#output_dir <- "/users/sidneybush/library/cloudstorage/box-box/05_storage_manuscript/05_outputs/hydrometric"

base_dir   <-"/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/03_Data"
output_dir <-"/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/05_Outputs/Hydrometric"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

site_order <- c(
  "GSWS09","GSWS10","GSWS01","GSLOOK","GSWS02",
  "GSWS03","GSWSMC","GSWS06","GSWS07","GSWS08"
)

# axis label mappings
axis_labels <- c(
  fdc_slope            = "Flow Duration Curve Slope",
  Q5norm               = "Normalized Q5",
  RBI                  = "Richards–Baker Index",
  recession_curve_slope = "Recession Curve Slope",
  StorageDepth         = "Storage Depth",
  StorageVolume        = "Storage Volume",
  DS_sum  = "DS Drawdown",
  S_annual_mm = "DS Storage-Discharge",
  mean_bf ="Mean Baseflow"
)

DS_drawdown <- read_csv(
  file.path(base_dir, "DynamicStorage", "DS_drawdown_annual.csv"),
  show_col_types = FALSE
)

# sitecode needs to be replaced with site and wateryear to year 
DS_drawdown <- DS_drawdown %>%
  rename(
    site = SITECODE,
    year = waterYear
  )
#changing gswsma to gswsmc is labeled as 
DS_drawdown$site[DS_drawdown$site == "GSWSMA"] <- "GSWSMC"

DS_SQ <- read_csv(
  file.path(base_dir, "DynamicStorage", "storage_discharge_method_annual_vol_per_site_wateryear.csv"),
  show_col_types = FALSE
)

DS_SQ <- DS_SQ %>%
  rename(
    year = wateryear
  )
#changing gsmack to gswsmc is labeled as 
DS_SQ$site[DS_SQ$site == "GSMACK"] <- "GSWSMC"
DS_SQ$site[DS_SQ$site == "GSLOOK_FULL"] <- "GSLOOK"




# read in storage metrics 
storage <- read_csv(
  file.path(base_dir, "DynamicStorage", "HJA_StorageMetrics_Annual.csv"),
  show_col_types = FALSE
) %>%
  select(-`...1`) %>%      
  rename(RBI = rbfi) %>%
  filter(site %in% site_order) %>%
  mutate(site = factor(site, levels = site_order))

HJA_storage<- left_join(storage, DS_SQ %>% select(year, site, S_annual_mm), by = c("year", "site"))

HJA_storage<- left_join(HJA_storage, DS_drawdown %>% select(year, site, DS_sum), by = c("year", "site"))

write.csv(HJA_storage,file.path(output_dir, "HJA_StorageMetrics_Annual_All.csv"))


storage_long <- HJA_storage %>%
  pivot_longer(
    cols      = -c(site, year),
    names_to  = "metric",
    values_to = "value"
  )

palette10 <- c(
  "#AA4499", "#882255", "#CC6677", "#DDCC77", "#999933",
  "#117733", "#44AA99", "#88CCEE", "#6699CC", "#332288"
)

site_cols <- setNames(lighten(palette10, amount = 0.1), site_order)

for (m in unique(storage_long$metric)) {
  df <- storage_long %>% filter(metric == m)
  axis_label <- if (m %in% names(axis_labels)) axis_labels[m] else m
    p_ts <- ggplot(df, aes(x = year, y = value, color = site, fill = site)) +
    geom_line(size = 0.6) +
    geom_point(size = 1) +
    facet_wrap(~ site, ncol = 2) +
    scale_color_manual(values = site_cols, guide = FALSE) +
    scale_fill_manual(values = alpha(site_cols, 0.3), guide = FALSE) +
    labs(x = "Water Year", y = axis_label) +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, size = 1),
      axis.line        = element_blank(),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text       = element_text(color = "black", hjust = 0)
    )
  ggsave(
    paste0("storage_0803_", m, "_by_site_wy.png"),
    p_ts, path = output_dir, width = 8, height = 9, units = "in", dpi = 300
  )
  sum_df <- df %>%
    group_by(site) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val   = sd(value,   na.rm = TRUE),
      .groups  = "drop"
    )
  
  # anova/ tukey 
  aov_res <- aov(value ~ site, data = df)
  tuk      <- TukeyHSD(aov_res, "site")$site
  pvals    <- setNames(tuk[, "p adj"], rownames(tuk))
  letters  <- multcompLetters(pvals)$Letters
  sum_df$group <- letters[as.character(sum_df$site)]
  
  y_max <- max(sum_df$mean_val + sum_df$sd_val, na.rm = TRUE)
  y_min <- min(sum_df$mean_val - sum_df$sd_val, na.rm = TRUE)
  offset <- 0.05 * (y_max - y_min)
  sum_df <- sum_df %>% mutate(label_y = mean_val + sd_val + offset)
  
  p_sum <- ggplot(sum_df, aes(x = site, y = mean_val, color = site)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                  width = 0.2) +
    geom_text(aes(label = group, y = label_y), size = 4, vjust = 0) +
    scale_x_discrete(limits = site_order) +
    scale_color_manual(values = site_cols, guide = FALSE) +
    labs(x = NULL, y = paste0("Mean ±1 SD of ", axis_label)) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.border     = element_rect(color = "black", fill = NA, size = 1),
      axis.line        = element_blank(),
      strip.text       = element_text(color = "black", hjust = 0)
    )
  
  ggsave(
    paste0("storage_summary_0803_", m, ".png"),
    p_sum, path = output_dir, width = 6, height = 4, units = "in", dpi = 300
  )
}


selected_metrics <- c("fdc_slope", "Q5norm", "RBI", "recession_curve_slope","S_annual_mm","DS_sum","mean_bf")

desc_names <- c(
  "Flow Duration Curve Slope",
  "Normalized Q5",
  "Richards–Baker Index",
  "Recession Curve Slope",
  "DS Storage-Discharge",
  "DS Drawdown",
  "Mean Baseflow"
)

metric_labels <- setNames(
  paste0(LETTERS[1:length(selected_metrics)], ") ", desc_names),
  selected_metrics
)

# compute summary + letters by metric
summary_sel <- storage_long %>%
  filter(metric %in% selected_metrics) %>%
  group_by(metric, site) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value,   na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  group_by(metric) %>%
  do({
    dat   <- .
    dfm   <- storage_long %>% filter(metric == dat$metric[1])
    aov_r <- aov(value ~ site, data = dfm)
    tuk_r <- TukeyHSD(aov_r, "site")$site
    pv    <- setNames(tuk_r[, "p adj"], rownames(tuk_r))
    let   <- multcompLetters(pv)$Letters
    dat$group <- let[as.character(dat$site)]
    y_mx <- max(dat$mean_val + dat$sd_val, na.rm = TRUE)
    y_mn <- min(dat$mean_val - dat$sd_val, na.rm = TRUE)
    off   <- 0.05 * (y_mx - y_mn)
    dat$label_y <- dat$mean_val + dat$sd_val + off
    dat
  }) %>% ungroup()


summary_sel$metric <- factor(summary_sel$metric, levels = c("fdc_slope", "Q5norm" ,"RBI", "recession_curve_slope", "S_annual_mm","DS_sum","mean_bf" ))

# plot grid + letters + a–d labels
p_grid <- ggplot(summary_sel, aes(x = site, y = mean_val, color = site)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                width = 0.2) +
  geom_text(aes(label = group, y = label_y), size = 3, vjust = 0) +
  facet_wrap(
    ~ metric,
    ncol   = 2,
    scales = "free_y",
    labeller = labeller(metric = metric_labels)
  ) +
  scale_x_discrete(limits = site_order) +
  scale_color_manual(values = site_cols, guide = FALSE) +
  labs(x = NULL, y = "Mean ± 1 SD") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black", hjust = 0)
  )

ggsave(
  filename = "selected_summary_grid_0803.png",
  plot     = p_grid,
  path     = output_dir,
  width    = 8, height   = 10,
  units    = "in", dpi = 300
)

# create scatter plots of dq/dt vs q for supplement: 
da_df <- read_csv(
  file.path(base_dir, "Q", "drainage_area.csv"),
  show_col_types = FALSE
)

# read in discharge, filter, compute q in m3/s
discharge <- read_csv(
  file.path(base_dir, "Q", "HF00402_v14.csv"),
  show_col_types = FALSE
) %>%
  filter(WATERYEAR > 1997, SITECODE %in% site_order) %>%
  left_join(da_df, by = "SITECODE") %>%
  filter(!is.na(DA_M2)) %>%
  mutate(
    Date = as.Date(DATE, "%m/%d/%Y"),
    Q    = MEAN_Q * 0.02831683199881   # m3/s
  )

recession_clean <- discharge %>%
  arrange(SITECODE, Date) %>%
  mutate(
    dQ             = Q - lag(Q),
    dQ_dt          = dQ / as.numeric(Date - lag(Date)),        
    recession_slope = -dQ_dt,
    Q_mm_day       = Q / DA_M2 * 86400 * 1000,                  
    slope_mm_day   = recession_slope / DA_M2 * 86400 * 1000     
  ) %>%
  filter(!is.na(Q_mm_day), !is.na(slope_mm_day)) %>%
  transmute(
    site          = factor(SITECODE, levels = site_order),
    Q_mm_day,
    slope_mm_day
  )

slopes_df <- recession_clean %>%
  group_by(site) %>%
  summarise(
    slope_mm     = abs(round(coef(lm(slope_mm_day ~ Q_mm_day))[2], 2)),
    Q_pos_mm     = quantile(Q_mm_day,     0.9,  na.rm = TRUE),
    slope_pos_mm = quantile(slope_mm_day, 0.85, na.rm = TRUE),
    .groups      = "drop"
  )

palette10 <- c(
  "#AA4499", "#882255", "#CC6677", "#DDCC77", "#999933",
  "#117733", "#44AA99", "#88CCEE", "#6699CC", "#332288"
)
site_cols <- setNames(lighten(palette10, amount = 0.1), site_order)

p_curve <- ggplot(recession_clean, aes(x = log(Q_mm_day), y = log(abs(slope_mm_day)), color = site)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  geom_text(
    data        = slopes_df,
    inherit.aes = FALSE,
    aes(
      x     = Q_pos_mm,
      y     = slope_pos_mm,
      label = paste0("slope = ", slope_mm)
    ),
    hjust = 0.3,
    vjust = -9,
    size  = 4,
    color = "black"
  ) +
  facet_wrap(~ site, ncol = 2, scales = "fixed") +
  scale_color_manual(values = site_cols, guide = FALSE) +
  labs(
    x = "log Q (mm day⁻¹)",
    y = "log –dQ/dt (mm day⁻¹)"
  ) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.line        = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(color = "black", hjust = 0)
  )

ggsave(
  filename = "recession_curve_by_site_mmday.png",
  plot     = p_curve,
  path     = output_dir,
  width    = 8,
  height   = 9,
  units    = "in",
  dpi      = 300
)
