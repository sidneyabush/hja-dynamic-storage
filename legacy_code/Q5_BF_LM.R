# 2 emma baseflow correlation with peak swe.
# inputs: output_dir/hja_ave_storagemetrics_catcharacter.csv.
# author: keira johnson
# date: 2026-02-13 (imported)

require(ggpmisc)
require(EflowStats)
require(dplyr)
require(ggpubr)
require(corrplot)
require(GGally)
require(car)
require(QuantPsyc)
require(plot.matrix)
require(tibble)
require(raster)
#install.packages("quantpsyc")

#define wds
base_dir <- "/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/03_Data"

output_dir <- "/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/05_Outputs/Hydrometric"


#bf_files<-c("coalcreekbaseflow.csv", "hjabaseflow.csv", "sagehenbaseflow.csv")

#loadfiles
HJA_Ave <- read_csv(
  file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
  show_col_types = FALSE
)[, -1]

site_data <- HJA_Ave[, -c(0:20)]

outcome_vars <- c(
  "recession_curve_slope_mean",
  "RBI_mean",
  "Q5norm_mean",
  "CV_Q5norm_mean",
  "mean_bf_mean",
  "fdc_slope_mean",
  "S_annual_mm_mean",
  "DR_Overall",
  "JST_AT_mean"
)

# create a named list of data frames

#"slope_mean",  "harvest", "landslide_young",  "landslide_total", "lava1_per", "lava2_per","ash_per","pyro_per"

# loop through each metric and create a separate data frame
for (col in metric_columns) {
  df <- HJA_Ave[, c(
    col,
    "Slope_mean",
    "Harvest",
    "Landslide_Young",
    "Landslide_Total",
    "Lava1_per",
    "Lava2_per",
    "Ash_Per",
    "Pyro_per"
  )]
  names(df)[1] <- col # Ensure the second column retains its original name
  assign(col, df) # Create a data frame named after the column
}

site_data <- list(
  recession_curve_slope_mean,
  RBI_mean,
  Q5norm_mean,
  CV_Q5norm_mean,
  mean_bf_mean,
  fdc_slope_mean,
  S_annual_mm_mean,
  DR_Overall,
  JST_AT_mean
)


#open lists for coefficient and r2 values to be appended into
coef_df_list <- list()
r2_df_list <- list()


#cycle through outcome variables
for (k in 1:length(HJA_list)) {
  #create model forumla for aic
  lm_mod_vars <- formula(paste(outcome_vars[k], "~."))

  #put into model
  lm_mod <- lm(lm_mod_vars, na.omit(site_data[[k]]))

  #lm_ols<-stepaic(lm_mod, direction = "backward")

  #run stepwise aic
  lm_AIC <- stepAIC(lm_mod, direction = "backward")

  #extract aic coefficients
  AIC_coef <- as.data.frame(lm_AIC$coefficients)

  #turn the rows to columns
  AIC_coef <- rownames_to_column(AIC_coef)

  #remove the intercept
  remove_int <- c("(Intercept)")

  #remove intercept column, cannot be inclded to check vif
  AIC_coef <- AIC_coef[!AIC_coef$rowname %in% "(Intercept)", ]

  #aic_coef$lm_input<-ifelse(aic_coef$rowname=="(intercept)", paste0(aic_coef$`lm_aic$coefficients`),
  #paste0(aic_coef$rowname, "*", aic_coef$`lm_aic$coefficients`))

  #create new model imput from rownames of retained variables from aic
  new_lm_input <- paste(AIC_coef$rowname, collapse = "+")

  #create new model input
  new_lm_vars <- formula(paste(outcome_vars[k], "~", new_lm_input))

  #put inot model structure
  lm_post_AIC <- lm(new_lm_vars, site_data[[k]])

  #model summary
  lm_sum <- summary(lm_post_AIC)

  #get coefficients (mostly just interested in p value here)
  coefs <- as.data.frame(lm_sum$coefficients)

  #get vif of retained variables
  vif_df <- if (length(AIC_coef$rowname) > 1) {
    #get vif of retained variables
    as.data.frame(vif(lm_post_AIC))
  } else {
    as.data.frame(NA)
  }

  colnames(vif_df) <- "vif(lm_post_AIC)"

  #get beta coefficiencts for retained variables
  beta_df <- as.data.frame(lm.beta(lm_post_AIC))
  beta_df <- rownames_to_column(beta_df)
  names(beta_df)[1] <- "Row.names"

  #append coefficients and vif together
  df_master <- merge(coefs, vif_df, by = "row.names", all.x = TRUE)

  df_master <- merge(df_master, beta_df, by = "Row.names", all.x = TRUE)

  #turn rownames column back into rownames
  df_master <- df_master %>%
    remove_rownames %>%
    column_to_rownames(var = "Row.names")

  #extract p value and vif columns
  df_master <- df_master[, c(
    "Pr(>|t|)",
    "vif(lm_post_AIC)",
    "lm.beta(lm_post_AIC)"
  )]

  #rename columns based on site and outcome variable
  colnames(df_master) <- c(
    paste0(outcome_vars[k], "_pvalue"),
    paste0(outcome_vars[k], "_VIF"),
    paste0(outcome_vars[k], "_beta")
  )

  #get r2 and adjusted r2
  r2 <- lm_sum$r.squared
  adj_r2 <- lm_sum$adj.r.squared

  #bind together
  r2_tot <- rbind(r2, adj_r2)

  colnames(r2_tot) <- paste0(outcome_vars[k])

  coef_df_list[[k]] <- df_master

  r2_df_list[[k]] <- r2_tot
}


# initialize an empty data frame to store results
final_df <- data.frame(
  Prefix = character(),
  Variable = character(),
  Beta = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

# loop through each dataframe in the list
for (df in coef_df_list) {
  # find columns that end with "_vif"
  vif_cols <- grep("_VIF$", names(df), value = TRUE)

  for (vif_col in vif_cols) {
    # extract the prefix before "_vif"
    prefix <- sub("_VIF$", "", vif_col)

    # construct corresponding beta and pvalue column names
    beta_col <- paste0(prefix, "_beta")
    pval_col <- paste0(prefix, "_pvalue")

    # check if beta and pvalue columns exist
    if (beta_col %in% names(df) && pval_col %in% names(df)) {
      # filter rows where vif value is 10 or less
      valid_indices <- which(df[[vif_col]] <= 10)

      # extract row names and corresponding beta and pvalue values
      valid_rows <- rownames(df)[valid_indices]
      beta_vals <- df[[beta_col]][valid_indices]
      pval_vals <- df[[pval_col]][valid_indices]

      # create a temporary data frame
      temp_df <- data.frame(
        Prefix = prefix,
        Variable = valid_rows,
        Beta = beta_vals,
        pvalue = pval_vals,
        stringsAsFactors = FALSE
      )

      # append to the final data frame
      final_df <- rbind(final_df, temp_df)
    }
  }
}

master_data_frame_coefs < master_data_frame_coefs <- data.frame(matrix(
  ncol = 1,
  nrow = 7
))

colnames_site_data <- colnames(site_data)
colnames_site_data <- colnames_site_data[2:7]

rownames(master_data_frame_coefs) <- c("(Intercept)", colnames_site_data)

for (i in 1:length(coef_df_list_master)) {
  internal_list <- coef_df_list_master[[i]]

  for (k in 1:length(internal_list)) {
    df <- internal_list[[k]]

    master_data_frame_coefs <- merge(
      master_data_frame_coefs,
      df,
      by = "row.names",
      all.x = TRUE
    )
    master_data_frame_coefs <- master_data_frame_coefs %>%
      remove_rownames %>%
      column_to_rownames(var = "Row.names")
  }
}


master_data_frame_coefs <- master_data_frame_coefs[, -1]

beta_cols <- grep("beta", colnames(master_data_frame_coefs))

VIF_cols <- grep("VIF", colnames(master_data_frame_coefs))

betas <- master_data_frame_coefs[, c(beta_cols)]

vifs <- master_data_frame_coefs[, c(VIF_cols)]

betas <- betas[rowSums(is.na(betas)) != ncol(betas), ]

vifs <- vifs[rowSums(is.na(vifs)) != ncol(vifs), ]

##rename second "coal creek" column to hja

colnames(betas) <- c(
  colnames(betas)[1],
  "Coal_Q5_beta",
  colnames(betas)[3],
  "Lookout_Q5_beta",
  colnames(betas)[c(5, 6)]
)

colnames(vifs) <- c(
  colnames(vifs)[1],
  "Coal_Q5_VIF",
  colnames(vifs)[3],
  "Lookout_Q5_VIF",
  colnames(vifs)[c(5, 6)]
)

#plot
betas <- as.matrix(betas)

betas[, 4] <- NA

values <- seq(-3, 2, 0.1)
ii <- cut(
  values,
  breaks = seq(min(values), max(values), len = 45),
  include.lowest = TRUE
)
## use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c(
  "firebrick4",
  "indianred",
  "white",
  "skyblue2",
  "dodgerblue4"
))(44)[ii]

nHalf = ((3 + 2) * 100) / 2
Min = -3
Max = 2
Thresh = 0

## make vector of colors for values below threshold
rc1 = colorRampPalette(
  colors = c("firebrick4", "indianred", "white"),
  space = "Lab"
)(nHalf)
## make vector of colors for values above threshold
rc2 = colorRampPalette(
  colors = c("white", "skyblue2", "dodgerblue3"),
  space = "Lab"
)(nHalf)
rampcols = c(rc1, rc2)

## in your example, this line sets the color for values between 49 and 51.
rampcols[c(nHalf)] = rgb(t(col2rgb("white")), maxColorValue = 256)

rb1 = seq(Min, Thresh, length.out = nHalf)
rb2 = seq(Thresh, Max, length.out = nHalf)[-1]
rampbreaks = c(rb1, rb2)

breaks <- seq(-3, 2, 1)


setwd("/Users/keirajohnson/Box Sync/Keira_Johnson")

#plot
jpeg(
  "Beta_Weight_Plot_0913.jpeg",
  width = 8,
  height = 7,
  units = "in",
  res = 1000
)
#increase margin
par(mar = c(11, 6, 6, 6) + .1, family = "Times New Roman")

plot(
  betas,
  col = rampcols,
  breaks = rampbreaks,
  border = NA,
  las = 2,
  ann = FALSE,
  digits = 2,
  na.print = FALSE,
  polygon.key = NULL,
  fmt.key = "%.0f",
  cex = 1.5
)


dev.off()


##for predicting q5 using only ssgw
coef_df_list_BF <- list()
r2_df_list_BF <- list()

sites <- c("Coal", "Lookout", "Sagehen")

for (i in 1:length(df_list)) {
  #select dataframe for site
  site_data <- df_list[[i]]

  #put into model
  lm_mod <- lm(Q5 ~ BFPer, site_data)

  lm_coeff <- as.data.frame(summary(lm_mod)$coefficients)

  #get beta coefficiencts for retained variables
  beta_df <- as.data.frame(lm.beta(lm_mod))

  #append coefficients and vif together
  df_master <- merge(lm_coeff, beta_df, by = "row.names", all.x = TRUE)

  #turn rownames column back into rownames
  df_master <- df_master %>%
    remove_rownames %>%
    column_to_rownames(var = "Row.names")

  #extract p value and vif columns
  df_master <- df_master[, c("Pr(>|t|)", "lm.beta(lm_mod)")]

  #rename columns based on site and outcome variable
  colnames(df_master) <- c(
    paste0(sites[i], "_Q5_onlyBFPer", "_pvalue"),
    paste0(sites[i], "_Q5_onlyBFPer", "_beta")
  )

  #get r2 and adjusted r2
  r2 <- summary(lm_mod)$r.squared
  adj_r2 <- summary(lm_mod)$adj.r.squared

  #bind together
  r2_tot <- rbind(r2, adj_r2)

  colnames(r2_tot) <- paste0(sites[i], "_Q5_onlyBFPer")

  coef_df_list_BF[[i]] <- df_master

  r2_df_list_BF[[i]] <- r2_tot
}

master_data_frame_coefs_BF <- data.frame(matrix(ncol = 1, nrow = 2))

rownames(master_data_frame_coefs_BF) <- c("(Intercept)", "BFPer")

for (i in 1:length(coef_df_list_BF)) {
  df <- coef_df_list_BF[[i]]

  master_data_frame_coefs_BF <- merge(
    master_data_frame_coefs_BF,
    df,
    by = "row.names",
    all.x = TRUE
  )
  master_data_frame_coefs_BF <- master_data_frame_coefs_BF %>%
    remove_rownames %>%
    column_to_rownames(var = "Row.names")
}

master_data_frame_coefs_BF <- master_data_frame_coefs_BF[, -1]

beta_cols_BF <- grep("beta", colnames(master_data_frame_coefs_BF))

betas_BF <- master_data_frame_coefs_BF[, c(beta_cols_BF)]

betas_BF <- betas_BF[rowSums(is.na(betas_BF)) != ncol(betas_BF), ]


# ##for predicting q5 using bfper and snow parameters
# coef_df_list_q5<-list()
# r2_df_list_q5<-list()
#
# sites<-c("coal", "lookout", "sagehen")
#
# for (i in 1:length(df_list)) {
#
#   #select dataframe for site
#   site_data<-df_list[[i]]
#
#   site_data<-site_data[,-14]
#
#   retained_vars<-grep(paste(sites[i]), colnames(betas))
#
#   site_betas<-as.data.frame(betas[,retained_vars])
#   colnames(site_betas)[1]<-"bfper"
#
#   q5_beta<-as.numeric(grep("bfper", colnames(site_betas)))
#
#   q5_beta_df<-as.data.frame(site_betas[,q5_beta])
#   rownames(q5_beta_df)<-rownames(betas)
#
#   all_names<-rownames(q5_beta_df)
#
#   lag_vars<-grep("lag", all_names)
#
#    if(i==1){
#
#      all_names<-all_names[-c(lag_vars)]
#
#    }else{
#
#      all_names<-all_names
#
#    }
#
#   if(i==1){
#
#     site_data<-site_data
#
#   }else{
#
#     site_data<-site_data[complete.cases(site_data),]
#
#   }
#
#
#
#   q5_retained_vars<-which(!is.na(q5_beta_df$`site_betas[, q5_beta]`))
#
#   q5_names<-rownames(q5_beta_df)[q5_retained_vars]
#
#   use_these_vars<-setdiff(all_names, q5_names)
#
#   use_collapse<-paste(use_these_vars, collapse = "+")
#
#   use_final<-paste0("bfper", "+", use_collapse)
#
#   q5_lm_vars<-formula(paste("q5", "~", use_final))
#
#   #put into model
#   lm_mod<-lm(q5_lm_vars, site_data)
#
#   lm_aic<-stepaic(lm_mod, direction = "backward")
#
#   #run stepwise aic
#   #lm_aic<-stepaic(lm_mod, direction = "both")
#
#   #extract aic coefficients
#   aic_coef<-as.data.frame(lm_aic$coefficients)
#
#   #turn the rows to columns
#   aic_coef<-rownames_to_column(aic_coef)
#
#   #remove the intercept
#   remove_int<-c("(intercept)")
#
#   #remove intercept column, cannot be inclded to check vif
#   aic_coef<-aic_coef[!aic_coef$rowname %in% "(intercept)",]
#
#   #aic_coef$lm_input<-ifelse(aic_coef$rowname=="(intercept)", paste0(aic_coef$`lm_aic$coefficients`),
#   #paste0(aic_coef$rowname, "*", aic_coef$`lm_aic$coefficients`))
#
#   #create new model imput from rownames of retained variables from aic
#   new_lm_input<-paste(aic_coef$rowname, collapse = "+")
#
#   #create new model input
#   new_lm_vars<-formula(paste("q5", "~", new_lm_input))
#
#   #put inot model structure
#   lm_post_aic<-lm(new_lm_vars, site_data)
#
#   #model summary
#   lm_sum<-summary(lm_post_aic)
#
#   #get coefficients (mostly just interested in p value here)
#   coefs<-as.data.frame(lm_sum$coefficients)
#
#   vif_df<-if(length(aic_coef$rowname) > 1){
#
#     #get vif of retained variables
#     as.data.frame(vif(lm_post_aic))
#
#   }else{
#
#     as.data.frame(na)
#
#   }
#
#   colnames(vif_df)<-"vif(lm_post_aic)"
#
#   #get beta coefficiencts for retained variables
#   beta_df<-as.data.frame(lm.beta(lm_post_aic))
#   beta_df<-rownames_to_column(beta_df)
#   names(beta_df)[1]<-"row.names"
#
#   #append coefficients and vif together
#   df_master<-merge(coefs, vif_df, by="row.names", all.x = true)
#
#   df_master<-merge(df_master, beta_df, by="row.names", all.x = true)
#
#   #turn rownames column back into rownames
#   df_master<-df_master %>% remove_rownames %>% column_to_rownames(var="row.names")
#
#   #extract p value and vif columns
#   df_master<-df_master[,c("pr(>|t|)", "vif(lm_post_aic)", "lm.beta(lm_post_aic)")]
#
#   #rename columns based on site and outcome variable
#   colnames(df_master)<-c(paste0(sites[i], "_q5bfper", "_pvalue"),
#                          paste0(sites[i], "_q5bfper", "_vif"),
#                          paste0(sites[i], "_q5bfper", "_beta"))
#
#   #get r2 and adjusted r2
#   r2<-lm_sum$r.squared
#   adj_r2<-lm_sum$adj.r.squared
#
#   #bind together
#   r2_tot<-rbind(r2, adj_r2)
#
#   colnames(r2_tot)<-paste0(sites[i], "_q5bfper")
#
#   coef_df_list_q5[[i]]<-df_master
#
#   r2_df_list_q5[[i]]<-r2_tot
#
# }

master_data_frame_coefs_Q5 <- data.frame(matrix(
  ncol = 1,
  nrow = length(rownames(betas)) + 2
))

colnames_site_data <- rownames(betas)

rownames(master_data_frame_coefs_Q5) <- c(
  "(Intercept)",
  "BFPer",
  colnames_site_data
)

for (i in 1:length(coef_df_list_Q5)) {
  internal_list <- coef_df_list_Q5[[i]]

  df <- internal_list

  master_data_frame_coefs_Q5 <- merge(
    master_data_frame_coefs_Q5,
    df,
    by = "row.names",
    all.x = TRUE
  )
  master_data_frame_coefs_Q5 <- master_data_frame_coefs_Q5 %>%
    remove_rownames %>%
    column_to_rownames(var = "Row.names")
}

master_data_frame_coefs_Q5 <- master_data_frame_coefs_Q5[, -1]

coal_y <- grep(".y", colnames(master_data_frame_coefs))
colnames(master_data_frame_coefs)[coal_y] <- c(
  "Lookout_Q5_pvalue",
  "Lookout_Q5_VIF",
  "Lookout_Q5_beta"
)

master_data_frame_coefs[coal_y] <- NA

master_df <- merge(
  master_data_frame_coefs,
  master_data_frame_coefs_BF,
  by = "row.names",
  all = TRUE
)

master_df <- column_to_rownames(master_df, var = "Row.names")

master_df <- merge(
  master_df,
  master_data_frame_coefs_BF,
  by = "row.names",
  all.x = TRUE
)

write.csv(master_df, "MasterLM_Coefs_DF.csv")
getwd()

beta_cols <- grep("beta", colnames(master_df))

betas_all <- master_df[, c(beta_cols)]
rownames(betas_all) <- master_df$Row.names

betas_master <- betas_all

betas_master <- betas_master[
  rowSums(is.na(betas_master)) != ncol(betas_master),
]

BF_betas <- betas_master[grep("_BFPer_beta", colnames(betas_master))]
BF_betas <- BF_betas[rowSums(is.na(BF_betas)) != ncol(BF_betas), ]

Q5_betas <- betas_master[, setdiff(colnames(betas_master), colnames(BF_betas))]
Q5_betas <- Q5_betas[rowSums(is.na(Q5_betas)) != ncol(Q5_betas), ]

Q5_betas <- Q5_betas[, order(colnames(Q5_betas))]

betas_BF_matrix <- as.matrix(BF_betas)
betas_BF_matrix <- betas_BF_matrix[-1, ]

betas_Q5_matrix <- as.matrix(Q5_betas)

values <- seq(-3, 2, 0.2)
ii <- cut(
  values,
  breaks = seq(min(values), max(values), len = 45),
  include.lowest = TRUE
)
## use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c(
  "firebrick4",
  "indianred",
  "goldenrod",
  "white",
  "skyblue2",
  "dodgerblue4"
))(44)[ii]

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson")
#plot
pdf("BF_Beta_Weight_Plot.pdf", width = 8, height = 7)
#increase margin
par(mar = c(11, 6, 6, 6) + .1)

plot(
  betas_BF_matrix,
  col = colors,
  border = NA,
  las = 2,
  ann = FALSE,
  digits = 2,
  na.print = FALSE
)

dev.off()

pdf("Q5_Beta_Weight_Plot.pdf", width = 8, height = 7)
#increase margin
par(mar = c(11, 6, 6, 6) + .1)

plot(
  betas_Q5_matrix,
  col = colors,
  border = NA,
  las = 2,
  ann = FALSE,
  digits = 2,
  na.print = FALSE
)

dev.off()

# #sagehen
#
# shoulder_agg<-df_list[[3]]
#
# v1<-lm(formula = q5 ~ peakswe + meltrate + totprecip + liquidprecip + swelag +
#    snowfrac + q5lag, data = shoulder_agg)
#
# summary(lm(formula = q5 ~ peakswe + snowfrac+
#          q5lag, data = shoulder_agg))
#
# stepaic(v1, direction = "backward")
#
# v2<-lm(formula = bfperrominbf99 ~ peakswe + meltrate + totprecip + swelag +
#          peakswewyday + zerowyday + snowfrac + bflag + q5lag, data = shoulder_agg)
#
# v3<-lm(formula = bfperrominbf99 ~ peakswe + meltrate + totprecip + swelag +
#           zerowyday + snowfrac + bflag + q5lag, data = shoulder_agg)
#
# v4<-lm(formula = bfperrominbf99 ~ peakswe + meltrate + totprecip +
#          zerowyday + snowfrac + bflag + q5lag, data = shoulder_agg)
#
# v5<-lm(formula = bfper ~  totprecip +
#           snowfrac + q5lag, data = shoulder_agg)
#
# summary(v5)
#
# vif(v5)
#
# lm.beta(v5)
#
# shoulder_agg$lm<-v5$coefficients[1]+v5$coefficients[2]*shoulder_agg$totprecip+
#   v5$coefficients[3]*shoulder_agg$snowfrac+v5$coefficients[4]*shoulder_agg$q5lag
#
# bfq5<-lm(q5~bfperrominbf99, data = shoulder_agg)
#
# lm.beta(bfq5)
#
# summary(bfq5)
#
# k1<-lm(formula = q5 ~ meltrate +
#          zerowyday + snowfrac + q5lag, data = shoulder_agg)
#
# x1<-lm(formula = q5 ~ bfperrominbf99 +bflag,
#        data = shoulder_agg)
#
# summary(k1)
#
# vif(x1)
#
# lm.beta(x1)
#
# shoulder_agg$bfq5<-k1$coefficients[1]+k1$coefficients[2]*shoulder_agg$meltrate+
#   k1$coefficients[3]*shoulder_agg$zerowyday+k1$coefficients[4]*shoulder_agg$snowfrac+
#   k1$coefficients[5]*shoulder_agg$q5lag
#
# shoulder_agg$totlm<-x1$coefficients[1]+x1$coefficients[2]*shoulder_agg$bfperrominbf99+
#   x1$coefficients[3]*shoulder_agg$bflag
#
# #use these plots
# p1<-ggplot(shoulder_agg, aes(lm, bfperrominbf99))+geom_point(col="black")+
#   labs(y="shoulder season gw proportion", x="predicted gw proportion")+
#   #geom_text(x=0.13, y=0.73, label="model a \n r2adj=0.95", cex=6)+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
#
# p1
# p2<-ggplot(shoulder_agg, aes(bfperrominbf99, q5))+geom_point(col="black")+
#   labs(x="shoulder season gw proportion", y="5th percentile q volume")+
#   #geom_text(x=0.16, y=0.06, label="model b \n r2adj=0.44", cex=6)+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
#
#
# p3<-ggplot(shoulder_agg, aes(bfq5, q5))+geom_point(col="black")+
#   labs(x="modeled 5th percetile q volume - no gw prop.", y="")+
#   #geom_text(x=0.034, y=0.059, label="model c \n r2adj=0.73", cex=6)+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
#
#
# p4<-ggplot(shoulder_agg, aes(totlm, q5))+geom_point(col="black")+
#   #geom_text(x=0.033, y=0.059, label="model d \n r2adj=0.64", cex=6)+
#   labs(x="modeled 5th percetile q volume - including gw prop.", y="")+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# p4
#
#
# jpeg("sagelfmplot.jpeg", width = 530, height = 500)
#
# p1
#
# dev.off()
#
# jpeg("sagessgwplots.jpeg", width = 1500, height = 450)
#
# ggarrange(p2,p3,p4, nrow = 1)
#
# dev.off()
#
# #lookout
#
# shoulder_agg<-df_list[[2]]
#
# summary(lm(formula = bfper ~ peakswe + meltrate + totprecip + swelag +
#          peakswewyday + zerowyday + snowfrac + bflag +q5lag, data = shoulder_agg))
#
# v1<-lm(formula = bfper ~ peakswe+
#              q5lag, data = shoulder_agg)
#
# summary(v1)
# lm.beta(v1)
#
# stepaic(object = v1, direction = "both", trace = true)
#
# lm_best<-lm(formula = q5 ~ bfper + snowfrac,
#            data = shoulder_agg)
# #get stats
# summary(lm_best)
# #plot(lm_best)
# #get beta coefficients and make df
# beta_df<-as.data.frame(lm.beta(lm_best))
# beta_df$name<-rownames(beta_df)
#
# #get vif values and make df
# vif_df<-as.data.frame(vif(lm_best))
# vif_df$name<-rownames(vif_df)
#
# shoulder_agg$lm<-lm_best$coefficients[1]+lm_best$coefficients[2]*shoulder_agg$bfper+
#   lm_best$coefficients[3]*shoulder_agg$snowfrac
#
# ggplot(shoulder_agg, aes(lm, q5))+geom_point()+geom_smooth(method = "lm", se=false, col="black")+
#   theme_classic()+labs(x="ssgw + snow fraction", y="q5")+theme(text = element_text(size = 20))
#
#
# #
# # v1<-lm(formula = bfperro01bf99 ~ peakswe + meltrate + totprecip + liquidprecip + swelag +
# #           peakswewyday + zerowyday + snowfrac + bflag +q5lag, data = shoulder_agg)
# #
# # v2<-lm(formula = bfperro01bf99 ~ peakswe + meltrate + totprecip + swelag +
# #          peakswewyday + zerowyday + snowfrac + bflag +q5lag, data = shoulder_agg)
# #
# # v3<-lm(formula = bfperro01bf99 ~ peakswe + meltrate + swelag +
# #          peakswewyday + zerowyday + snowfrac + bflag +q5lag, data = shoulder_agg)
# #
# # v4<-lm(formula = bfperro01bf99 ~ peakswe + meltrate + swelag +
# #          peakswewyday + zerowyday + bflag +q5lag, data = shoulder_agg)
# #
# # v5<-lm(formula = bfperro01bf99 ~ peakswe+meltrate + swelag +
# #          peakswewyday+q5lag, data = shoulder_agg)
#
# shoulder_agg$lm<-v5$coefficients[1]+v5$coefficients[2]*shoulder_agg$peakswe+
#   v5$coefficients[3]*shoulder_agg$meltrate+v5$coefficients[4]*shoulder_agg$swelag+
#   v5$coefficients[5]*shoulder_agg$peakswewyday+v5$coefficients[6]*shoulder_agg$q5lag
#
# bfq5<-lm(q5~bfperro01bf99, data = shoulder_agg)
#
# summary(bfq5)
# lm.beta(bfq5)
#
# k1<-lm(formula = q5 ~ peakswe + meltrate + totprecip + liquidprecip + swelag +
#          peakswewyday + zerowyday + snowfrac + bflag +q5lag, data = shoulder_agg)
#
# stepaic(object = k1, direction = "both", trace = true)
#
# k1<-lm(formula = q5 ~ snowfrac+swelag+
#         q5lag, data = shoulder_agg)
#
# k_best<-lm(formula = q5 ~ meltrate + totprecip + swelag +
#              peakswewyday + zerowyday +q5lag, data = shoulder_agg)
#
# #get stats
# sum<-summary(k_best)
# #plot(k_best)
# #get beta coefficients and make df
# beta_df<-as.data.frame(lm.beta(k_best))
# beta_df$name<-rownames(beta_df)
#
# #get vif values and make df
# vif_df<-as.data.frame(vif(k_best))
# vif_df$name<-rownames(vif_df)
#
# #k1<-lm(formula = q5 ~ peakswe+liquidprecip + swelag +
#          #peakswewyday + zerowyday +q5lag, data = shoulder_agg)
#
# shoulder_agg$bfq5<-k1$coefficients[1]+k1$coefficients[2]*shoulder_agg$snowfrac+
#   k1$coefficients[3]*shoulder_agg$swelag+k1$coefficients[4]*shoulder_agg$q5lag
#
# x1<-lm(formula = q5 ~ bfperro01bf99 +
#          snowfrac, data = shoulder_agg)
#
# x1<-lm(formula = q5 ~ bfperro01bf99 + totprecip + liquidprecip +
#           zerowyday + snowfrac + bflag, data = shoulder_agg)
#
# stepaic(object=x1, direction = "both", trace = true)
#
# shoulder_agg$lmtot<-x1$coefficients[1]+x1$coefficients[2]*shoulder_agg$bfperro01bf99+
#   x1$coefficients[3]*shoulder_agg$snowfrac
#
#
# p1<-ggplot(shoulder_agg, aes(lm, bfperro01bf99))+geom_point(col="black")+
#   labs(y="shoulder season gw proportion", x="predicted gw proportion")+
#   #geom_text(x=0.37, y=0.59, label="model a \n r2adj=0.74", cex=6)+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
#
#
# p2<-ggplot(shoulder_agg, aes(bfperro01bf99, q5))+geom_point(col="black")+
#   labs(x="shoulder season gw proportion", y="5th percentile q volume")+
#   #geom_text(x=0.39, y=0.39, label="model b \n r2adj=0.70", cex=6)+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
#
#
# p3<-ggplot(shoulder_agg, aes(bfq5, q5))+geom_point(col="black")+
#   labs(x="modeled 5th percetile q volume - no gw prop.", y="")+
#   #geom_text(x=0.2, y=0.385, label="model c \n r2adj=0.64", cex=6)+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
#
#
# p4<-ggplot(shoulder_agg, aes(lmtot, q5))+geom_point(col="black")+
#   #geom_text(x=0.17, y=0.385, label="model d \n r2adj=0.76", cex=6)+
#   labs(x="modeled 5th percetile q volume - including gw prop.", y="")+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# p4
#
# jpeg("lookoutlfmplot.jpeg", width = 530, height = 500)
#
# p1
#
# dev.off()
#
# jpeg("lookoutlmplots.jpeg", width = 1500, height = 450)
#
# ggarrange(p2,p3,p4, nrow = 1)
#
# dev.off()
#
# #coal
#
# v1<-lm(formula = bfperro01bf99 ~ meltrate+liquidprecip,
#          data = shoulder_agg)
#
# x1<-lm(formula = q5 ~ peakswe+liquidprecip,
#        data = shoulder_agg)
#
# k1<-lm(formula = q5 ~ bfperro01bf99+peakswe,
#        data = shoulder_agg)
#
# q1<-lm(formula = q5 ~ bfperro01bf99,
#        data = shoulder_agg)
#
# summary(k1)
#
# lm.beta(k1)
#
# vif(k1)
#
# shoulder_agg$lm<-v1$coefficients[1]+v1$coefficients[2]*shoulder_agg$meltrate+
#   v1$coefficients[3]*shoulder_agg$liquidprecip
#
# shoulder_agg$bfq5<-x1$coefficients[1]+x1$coefficients[2]*shoulder_agg$peakswe+
#   x1$coefficients[3]*shoulder_agg$liquidprecip
#
# shoulder_agg$totlm<-k1$coefficients[1]+k1$coefficients[2]*shoulder_agg$bfperro01bf99+
#   k1$coefficients[3]*shoulder_agg$peakswe
#
#
# p1<-ggplot(shoulder_agg, aes(lm, bfperro01bf99))+geom_point(col="black")+
#   labs(y="shoulder season gw proportion", x="predicted gw proportion")+
#   #geom_text(x=0.25, y=0.63, label="model a \n r2adj=0.98", cex=6)+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# p1
#
# p2<-ggplot(shoulder_agg, aes(bfperro01bf99, q5))+geom_point(col="black")+
#   labs(x="shoulder season gw proportion", y="5th percentile q volume")+
#   #geom_text(x=0.25, y=0.06, label="model b \n r2adj=0.14", cex=6)+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
#
#
# p3<-ggplot(shoulder_agg, aes(bfq5, q5))+geom_point(col="black")+
#   labs(x="modeled 5th percetile q volume - no gw prop.", y="")+
#   #geom_text(x=0.027, y=0.059, label="model c \n r2adj=0.76", cex=6)+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
#
#
# p4<-ggplot(shoulder_agg, aes(totlm, q5))+geom_point(col="black")+
#   #geom_text(x=0.027, y=0.059, label="model d \n r2adj=0.82", cex=6)+
#   labs(x="modeled 5th percetile q volume - including gw prop.", y="")+
#   geom_smooth(method = "lm", se=f, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# p4
#
#
# jpeg("coallfmplot.jpeg", width = 530, height = 500)
#
# p1
#
# dev.off()
#
#
#
# jpeg("coalssgwplots.jpeg", width = 1500, height = 450)
# ggarrange(p2,p3, p4, nrow = 1)
# dev.off()
