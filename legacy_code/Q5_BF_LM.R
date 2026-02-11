#2 EMMA baseflow correlation with peak swe
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
#install.packages("QuantPsyc")

#define wds
base_dir   <-"/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/03_Data"

output_dir <-"/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/05_Outputs/Hydrometric"


#BF_files<-c("CoalCreekBaseflow.csv", "HJABaseflow.csv", "SagehenBaseflow.csv")

#Loadfiles 
HJA_Ave<-read_csv(
  file.path(output_dir, "HJA_Ave_StorageMetrics_CatCharacter.csv"),
  show_col_types = FALSE
)[,-1] 

site_data<-HJA_Ave[,-c(0:20)]

outcome_vars <- c(
  "recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", "CV_Q5norm_mean",
  "mean_bf_mean", "fdc_slope_mean", "S_annual_mm_mean", "DR_Overall", "JST_AT_mean"
)

# Create a named list of data frames

#"Slope_mean",  "Harvest", "Landslide_Young",  "Landslide_Total", "Lava1_per", "Lava2_per","Ash_Per","Pyro_per"

# Loop through each metric and create a separate data frame
for (col in metric_columns) {
  df <- HJA_Ave[, c(col,"Slope_mean",  "Harvest", "Landslide_Young",  "Landslide_Total", "Lava1_per", "Lava2_per","Ash_Per","Pyro_per")]
  names(df)[1] <- col  # Ensure the second column retains its original name
  assign(col, df)      # Create a data frame named after the column
}

site_data<-list(recession_curve_slope_mean, RBI_mean, Q5norm_mean, CV_Q5norm_mean,mean_bf_mean, fdc_slope_mean, S_annual_mm_mean, DR_Overall, JST_AT_mean)



#open lists for coefficient and r2 values to be appended into
coef_df_list<-list()
r2_df_list<-list()


  #cycle through outcome variables
  for (k in 1:length(HJA_list)) {
    
    #create model forumla for AIC
    lm_mod_vars<-formula(paste(outcome_vars[k], "~."))
    
    #put into model
    lm_mod<-lm(lm_mod_vars, na.omit(site_data[[k]]))
    
    #lm_ols<-stepAIC(lm_mod, direction = "backward")
    
    #run stepwise AIC
    lm_AIC<-stepAIC(lm_mod, direction = "backward")
    
    #extract AIC coefficients
    AIC_coef<-as.data.frame(lm_AIC$coefficients)
    
    #turn the rows to columns
    AIC_coef<-rownames_to_column(AIC_coef)
    
    #remove the intercept
    remove_int<-c("(Intercept)")
    
    #remove intercept column, cannot be inclded to check VIF
    AIC_coef<-AIC_coef[!AIC_coef$rowname %in% "(Intercept)",]
    
    #AIC_coef$lm_input<-ifelse(AIC_coef$rowname=="(Intercept)", paste0(AIC_coef$`lm_AIC$coefficients`),
                              #paste0(AIC_coef$rowname, "*", AIC_coef$`lm_AIC$coefficients`))
    
    #create new model imput from rownames of retained variables from AIC
    new_lm_input<-paste(AIC_coef$rowname, collapse = "+")
    
    #create new model input
    new_lm_vars<-formula(paste(outcome_vars[k], "~", new_lm_input))
    
    #put inot model structure
    lm_post_AIC<-lm(new_lm_vars, site_data[[k]])
    
    #model summary
    lm_sum<-summary(lm_post_AIC)
    
    #get coefficients (mostly just interested in p value here)
    coefs<-as.data.frame(lm_sum$coefficients)
    
    #get VIF of retained variables
    vif_df<-if(length(AIC_coef$rowname) > 1){
      
      #get VIF of retained variables
      as.data.frame(vif(lm_post_AIC))
      
    }else{
      
      as.data.frame(NA)
      
    }
    
    colnames(vif_df)<-"vif(lm_post_AIC)"
    
    #get beta coefficiencts for retained variables
    beta_df<-as.data.frame(lm.beta(lm_post_AIC))
    beta_df<-rownames_to_column(beta_df)
    names(beta_df)[1]<-"Row.names"
    
    #append coefficients and VIF together
    df_master<-merge(coefs, vif_df, by="row.names", all.x = TRUE)
    
    df_master<-merge(df_master, beta_df, by="Row.names", all.x = TRUE)
    
    #turn rownames column back into rownames
    df_master<-df_master %>% remove_rownames %>% column_to_rownames(var="Row.names")
    
    #extract p value and VIF columns
    df_master<-df_master[,c("Pr(>|t|)", "vif(lm_post_AIC)", "lm.beta(lm_post_AIC)")]
    
    #rename columns based on site and outcome variable
    colnames(df_master)<-c(paste0(outcome_vars[k], "_pvalue"), 
                           paste0(outcome_vars[k], "_VIF"),
                                  paste0(outcome_vars[k], "_beta"))
    
    #get r2 and adjusted r2
    r2<-lm_sum$r.squared
    adj_r2<-lm_sum$adj.r.squared
    
    #bind together
    r2_tot<-rbind(r2, adj_r2)
    
    colnames(r2_tot)<-paste0(outcome_vars[k])
    
    coef_df_list[[k]]<-df_master
    
    r2_df_list[[k]]<-r2_tot

  }
  


# Initialize an empty data frame to store results
final_df <- data.frame(Prefix = character(), Variable = character(), Beta = numeric(), pvalue = numeric(), stringsAsFactors = FALSE)

# Loop through each dataframe in the list
for (df in coef_df_list) {
  # Find columns that end with "_VIF"
  vif_cols <- grep("_VIF$", names(df), value = TRUE)
  
  for (vif_col in vif_cols) {
    # Extract the prefix before "_VIF"
    prefix <- sub("_VIF$", "", vif_col)
    
    # Construct corresponding beta and pvalue column names
    beta_col <- paste0(prefix, "_beta")
    pval_col <- paste0(prefix, "_pvalue")
    
    # Check if beta and pvalue columns exist
    if (beta_col %in% names(df) && pval_col %in% names(df)) {
      # Filter rows where VIF value is 10 or less
      valid_indices <- which(df[[vif_col]] <= 10)
      
      # Extract row names and corresponding beta and pvalue values
      valid_rows <- rownames(df)[valid_indices]
      beta_vals <- df[[beta_col]][valid_indices]
      pval_vals <- df[[pval_col]][valid_indices]
      
      # Create a temporary data frame
      temp_df <- data.frame(
        Prefix = prefix,
        Variable = valid_rows,
        Beta = beta_vals,
        pvalue = pval_vals,
        stringsAsFactors = FALSE
      )
      
      # Append to the final data frame
      final_df <- rbind(final_df, temp_df)
    }
  }
}

master_data_frame_coefs<

master_data_frame_coefs<-data.frame(matrix(ncol=1, nrow=7))

colnames_site_data<-colnames(site_data)
colnames_site_data<-colnames_site_data[2:7]

rownames(master_data_frame_coefs)<-c("(Intercept)", colnames_site_data)

for (i in 1:length(coef_df_list_master)) {
  
  internal_list<-coef_df_list_master[[i]]
  
  for (k in 1:length(internal_list)) {
    
    df<-internal_list[[k]]
    
    master_data_frame_coefs<-merge(master_data_frame_coefs, df, by="row.names", all.x = TRUE)
    master_data_frame_coefs<-master_data_frame_coefs %>% remove_rownames %>% 
      column_to_rownames(var="Row.names")
    
  }
  
}


master_data_frame_coefs<-master_data_frame_coefs[,-1]

beta_cols<-grep("beta", colnames(master_data_frame_coefs))

VIF_cols<-grep("VIF", colnames(master_data_frame_coefs))

betas<-master_data_frame_coefs[, c(beta_cols)]

vifs<-master_data_frame_coefs[, c(VIF_cols)]

betas<-betas[rowSums(is.na(betas)) != ncol(betas), ]

vifs<-vifs[rowSums(is.na(vifs)) != ncol(vifs), ]

##rename second "Coal Creek" column to HJA

colnames(betas)<-c(colnames(betas)[1], "Coal_Q5_beta", colnames(betas)[3], "Lookout_Q5_beta",
                   colnames(betas)[c(5,6)])

colnames(vifs)<-c(colnames(vifs)[1], "Coal_Q5_VIF", colnames(vifs)[3], "Lookout_Q5_VIF",
                   colnames(vifs)[c(5,6)])

#plot
betas<-as.matrix(betas)

betas[,4]<-NA

values<-seq(-3,2,0.1)
ii <- cut(values, breaks = seq(min(values), max(values), len = 45), 
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("firebrick4", "indianred","white", "skyblue2", "dodgerblue4"))(44)[ii]

nHalf = ((3+2)*100)/2
Min = -3
Max = 2
Thresh = 0

## Make vector of colors for values below threshold
rc1 = colorRampPalette(colors = c("firebrick4", "indianred", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("white","skyblue2", "dodgerblue3"), space="Lab")(nHalf)
rampcols = c(rc1, rc2)

## In your example, this line sets the color for values between 49 and 51. 
rampcols[c(nHalf)] = rgb(t(col2rgb("white")), maxColorValue=256) 

rb1 = seq(Min, Thresh, length.out=nHalf)
rb2 = seq(Thresh, Max, length.out=nHalf)[-1]
rampbreaks = c(rb1, rb2)

breaks<-seq(-3,2,1)


setwd("/Users/keirajohnson/Box Sync/Keira_Johnson")

#plot
jpeg("Beta_Weight_Plot_0913.jpeg", width = 8, height = 7, units = "in", res = 1000)
#increase margin
par(mar=c(11,6,6,6)+.1, family="Times New Roman")

plot(betas, col = rampcols, breaks = rampbreaks, border=NA, las=2, ann=FALSE,  digits = 2, 
     na.print = FALSE, polygon.key = NULL, fmt.key="%.0f", cex=1.5)



dev.off()




##for predicting Q5 using only SSGW
coef_df_list_BF<-list()
r2_df_list_BF<-list()

sites<-c("Coal", "Lookout", "Sagehen")

for (i in 1:length(df_list)) {
  
  #select dataframe for site
  site_data<-df_list[[i]]
  
  #put into model
  lm_mod<-lm(Q5~BFPer, site_data)
  
  lm_coeff<-as.data.frame(summary(lm_mod)$coefficients)
  
  #get beta coefficiencts for retained variables
  beta_df<-as.data.frame(lm.beta(lm_mod))
  
  #append coefficients and VIF together
  df_master<-merge(lm_coeff, beta_df, by="row.names", all.x = TRUE)
  
  #turn rownames column back into rownames
  df_master<-df_master %>% remove_rownames %>% column_to_rownames(var="Row.names")
  
  #extract p value and VIF columns
  df_master<-df_master[,c("Pr(>|t|)","lm.beta(lm_mod)")]
  
  #rename columns based on site and outcome variable
  colnames(df_master)<-c(paste0(sites[i], "_Q5_onlyBFPer", "_pvalue"), 
                         paste0(sites[i], "_Q5_onlyBFPer", "_beta"))
  
  #get r2 and adjusted r2
  r2<-summary(lm_mod)$r.squared
  adj_r2<-summary(lm_mod)$adj.r.squared
  
  #bind together
  r2_tot<-rbind(r2, adj_r2)
  
  colnames(r2_tot)<-paste0(sites[i], "_Q5_onlyBFPer")
  
  coef_df_list_BF[[i]]<-df_master
  
  r2_df_list_BF[[i]]<-r2_tot
  
}

master_data_frame_coefs_BF<-data.frame(matrix(ncol=1, nrow=2))

rownames(master_data_frame_coefs_BF)<-c("(Intercept)", "BFPer")

for (i in 1:length(coef_df_list_BF)) {
  
  df<-coef_df_list_BF[[i]]
    
  master_data_frame_coefs_BF<-merge(master_data_frame_coefs_BF, df, by="row.names", all.x = TRUE)
  master_data_frame_coefs_BF<-master_data_frame_coefs_BF %>% remove_rownames %>% 
      column_to_rownames(var="Row.names")
  
}

master_data_frame_coefs_BF<-master_data_frame_coefs_BF[,-1]

beta_cols_BF<-grep("beta", colnames(master_data_frame_coefs_BF))

betas_BF<-master_data_frame_coefs_BF[, c(beta_cols_BF)]

betas_BF<-betas_BF[rowSums(is.na(betas_BF)) != ncol(betas_BF), ]




# ##for predicting Q5 using BFPer and snow parameters
# coef_df_list_Q5<-list()
# r2_df_list_Q5<-list()
# 
# sites<-c("Coal", "Lookout", "Sagehen")
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
#   colnames(site_betas)[1]<-"BFPer"
#   
#   Q5_beta<-as.numeric(grep("BFPer", colnames(site_betas)))
#   
#   Q5_beta_df<-as.data.frame(site_betas[,Q5_beta])
#   rownames(Q5_beta_df)<-rownames(betas)
#   
#   all_names<-rownames(Q5_beta_df)
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
#   Q5_retained_vars<-which(!is.na(Q5_beta_df$`site_betas[, Q5_beta]`))
#   
#   Q5_names<-rownames(Q5_beta_df)[Q5_retained_vars]
#   
#   use_these_vars<-setdiff(all_names, Q5_names)
#   
#   use_collapse<-paste(use_these_vars, collapse = "+")
#   
#   use_final<-paste0("BFPer", "+", use_collapse)
#   
#   Q5_lm_vars<-formula(paste("Q5", "~", use_final))
#   
#   #put into model
#   lm_mod<-lm(Q5_lm_vars, site_data)
#   
#   lm_AIC<-stepAIC(lm_mod, direction = "backward")
#   
#   #run stepwise AIC
#   #lm_AIC<-stepAIC(lm_mod, direction = "both")
#   
#   #extract AIC coefficients
#   AIC_coef<-as.data.frame(lm_AIC$coefficients)
#   
#   #turn the rows to columns
#   AIC_coef<-rownames_to_column(AIC_coef)
#   
#   #remove the intercept
#   remove_int<-c("(Intercept)")
#   
#   #remove intercept column, cannot be inclded to check VIF
#   AIC_coef<-AIC_coef[!AIC_coef$rowname %in% "(Intercept)",]
#   
#   #AIC_coef$lm_input<-ifelse(AIC_coef$rowname=="(Intercept)", paste0(AIC_coef$`lm_AIC$coefficients`),
#   #paste0(AIC_coef$rowname, "*", AIC_coef$`lm_AIC$coefficients`))
#   
#   #create new model imput from rownames of retained variables from AIC
#   new_lm_input<-paste(AIC_coef$rowname, collapse = "+")
#   
#   #create new model input
#   new_lm_vars<-formula(paste("Q5", "~", new_lm_input))
#   
#   #put inot model structure
#   lm_post_AIC<-lm(new_lm_vars, site_data)
#   
#   #model summary
#   lm_sum<-summary(lm_post_AIC)
#   
#   #get coefficients (mostly just interested in p value here)
#   coefs<-as.data.frame(lm_sum$coefficients)
#   
#   vif_df<-if(length(AIC_coef$rowname) > 1){
# 
#     #get VIF of retained variables
#     as.data.frame(vif(lm_post_AIC))
#     
#   }else{
# 
#     as.data.frame(NA)
#     
#   }
#   
#   colnames(vif_df)<-"vif(lm_post_AIC)"
#   
#   #get beta coefficiencts for retained variables
#   beta_df<-as.data.frame(lm.beta(lm_post_AIC))
#   beta_df<-rownames_to_column(beta_df)
#   names(beta_df)[1]<-"Row.names"
#   
#   #append coefficients and VIF together
#   df_master<-merge(coefs, vif_df, by="row.names", all.x = TRUE)
#   
#   df_master<-merge(df_master, beta_df, by="Row.names", all.x = TRUE)
#   
#   #turn rownames column back into rownames
#   df_master<-df_master %>% remove_rownames %>% column_to_rownames(var="Row.names")
#   
#   #extract p value and VIF columns
#   df_master<-df_master[,c("Pr(>|t|)", "vif(lm_post_AIC)", "lm.beta(lm_post_AIC)")]
#   
#   #rename columns based on site and outcome variable
#   colnames(df_master)<-c(paste0(sites[i], "_Q5BFPer", "_pvalue"), 
#                          paste0(sites[i], "_Q5BFPer", "_VIF"),
#                          paste0(sites[i], "_Q5BFPer", "_beta"))
#   
#   #get r2 and adjusted r2
#   r2<-lm_sum$r.squared
#   adj_r2<-lm_sum$adj.r.squared
#   
#   #bind together
#   r2_tot<-rbind(r2, adj_r2)
#   
#   colnames(r2_tot)<-paste0(sites[i], "_Q5BFPer")
#   
#   coef_df_list_Q5[[i]]<-df_master
#   
#   r2_df_list_Q5[[i]]<-r2_tot
#   
# }

master_data_frame_coefs_Q5<-data.frame(matrix(ncol=1, nrow=length(rownames(betas))+2))

colnames_site_data<-rownames(betas)

rownames(master_data_frame_coefs_Q5)<-c("(Intercept)", "BFPer", colnames_site_data)

for (i in 1:length(coef_df_list_Q5)) {
  
  internal_list<-coef_df_list_Q5[[i]]
    
    df<-internal_list
    
    master_data_frame_coefs_Q5<-merge(master_data_frame_coefs_Q5, df, by="row.names", all.x = TRUE)
    master_data_frame_coefs_Q5<-master_data_frame_coefs_Q5 %>% remove_rownames %>% 
      column_to_rownames(var="Row.names")
  
}

master_data_frame_coefs_Q5<-master_data_frame_coefs_Q5[,-1]

coal_y<-grep(".y", colnames(master_data_frame_coefs))
colnames(master_data_frame_coefs)[coal_y]<-c("Lookout_Q5_pvalue", "Lookout_Q5_VIF",
                                             "Lookout_Q5_beta")

master_data_frame_coefs[coal_y]<-NA

master_df<-merge(master_data_frame_coefs, master_data_frame_coefs_BF, by="row.names", all=TRUE)

master_df<-column_to_rownames(master_df, var="Row.names")

master_df<-merge(master_df, master_data_frame_coefs_BF, by="row.names", all.x=TRUE)

write.csv(master_df, "MasterLM_Coefs_DF.csv")
getwd()

beta_cols<-grep("beta", colnames(master_df))

betas_all<-master_df[, c(beta_cols)]
rownames(betas_all)<-master_df$Row.names

betas_master<-betas_all

betas_master<-betas_master[rowSums(is.na(betas_master)) != ncol(betas_master), ]

BF_betas<-betas_master[grep("_BFPer_beta", colnames(betas_master))]
BF_betas<-BF_betas[rowSums(is.na(BF_betas)) != ncol(BF_betas), ]

Q5_betas<-betas_master[,setdiff(colnames(betas_master), colnames(BF_betas))]
Q5_betas<-Q5_betas[rowSums(is.na(Q5_betas)) != ncol(Q5_betas), ]

Q5_betas<-Q5_betas[,order(colnames(Q5_betas))]

betas_BF_matrix<-as.matrix(BF_betas)
betas_BF_matrix<-betas_BF_matrix[-1,]

betas_Q5_matrix<-as.matrix(Q5_betas)

values<-seq(-3,2,0.2)
ii <- cut(values, breaks = seq(min(values), max(values), len = 45), 
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("firebrick4", "indianred", "goldenrod","white", "skyblue2", "dodgerblue4"))(44)[ii]

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson")
#plot
pdf("BF_Beta_Weight_Plot.pdf", width = 8, height = 7)
#increase margin
par(mar=c(11,6,6,6)+.1)

plot(betas_BF_matrix, col = colors, border=NA, las=2, ann=FALSE,  digits = 2, na.print = FALSE)

dev.off()

pdf("Q5_Beta_Weight_Plot.pdf", width = 8, height = 7)
#increase margin
par(mar=c(11,6,6,6)+.1)

plot(betas_Q5_matrix, col = colors, border=NA, las=2, ann=FALSE, digits = 2, na.print = FALSE)

dev.off()

# #Sagehen
# 
# shoulder_agg<-df_list[[3]]
# 
# v1<-lm(formula = Q5 ~ PeakSWE + MeltRate + TotPrecip + LiquidPrecip + SWElag +
#    SnowFrac + Q5lag, data = shoulder_agg)
# 
# summary(lm(formula = Q5 ~ PeakSWE + SnowFrac+
#          Q5lag, data = shoulder_agg))
# 
# stepAIC(v1, direction = "backward")
# 
# v2<-lm(formula = BFPerROminBF99 ~ PeakSWE + MeltRate + TotPrecip + SWElag +
#          PeakSWEWYDay + ZeroWYDay + SnowFrac + BFlag + Q5lag, data = shoulder_agg)
# 
# v3<-lm(formula = BFPerROminBF99 ~ PeakSWE + MeltRate + TotPrecip + SWElag +
#           ZeroWYDay + SnowFrac + BFlag + Q5lag, data = shoulder_agg)
# 
# v4<-lm(formula = BFPerROminBF99 ~ PeakSWE + MeltRate + TotPrecip +
#          ZeroWYDay + SnowFrac + BFlag + Q5lag, data = shoulder_agg)
# 
# v5<-lm(formula = BFPer ~  TotPrecip +
#           SnowFrac + Q5lag, data = shoulder_agg)
# 
# summary(v5)
# 
# vif(v5)
# 
# lm.beta(v5)
# 
# shoulder_agg$lm<-v5$coefficients[1]+v5$coefficients[2]*shoulder_agg$TotPrecip+
#   v5$coefficients[3]*shoulder_agg$SnowFrac+v5$coefficients[4]*shoulder_agg$Q5lag
# 
# BFQ5<-lm(Q5~BFPerROminBF99, data = shoulder_agg)
# 
# lm.beta(BFQ5)
# 
# summary(BFQ5)
# 
# k1<-lm(formula = Q5 ~ MeltRate +
#          ZeroWYDay + SnowFrac + Q5lag, data = shoulder_agg)
# 
# x1<-lm(formula = Q5 ~ BFPerROminBF99 +BFlag,
#        data = shoulder_agg)
# 
# summary(k1)
# 
# vif(x1)
# 
# lm.beta(x1)
# 
# shoulder_agg$BFQ5<-k1$coefficients[1]+k1$coefficients[2]*shoulder_agg$MeltRate+
#   k1$coefficients[3]*shoulder_agg$ZeroWYDay+k1$coefficients[4]*shoulder_agg$SnowFrac+
#   k1$coefficients[5]*shoulder_agg$Q5lag
# 
# shoulder_agg$totlm<-x1$coefficients[1]+x1$coefficients[2]*shoulder_agg$BFPerROminBF99+
#   x1$coefficients[3]*shoulder_agg$BFlag
# 
# #use these plots
# p1<-ggplot(shoulder_agg, aes(lm, BFPerROminBF99))+geom_point(col="black")+
#   labs(y="Shoulder Season GW Proportion", x="Predicted GW Proportion")+
#   #geom_text(x=0.13, y=0.73, label="Model A \n R2adj=0.95", cex=6)+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# 
# p1
# p2<-ggplot(shoulder_agg, aes(BFPerROminBF99, Q5))+geom_point(col="black")+
#   labs(x="Shoulder Season GW Proportion", y="5th Percentile Q Volume")+
#   #geom_text(x=0.16, y=0.06, label="Model B \n R2adj=0.44", cex=6)+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# 
# 
# p3<-ggplot(shoulder_agg, aes(BFQ5, Q5))+geom_point(col="black")+
#   labs(x="Modeled 5th Percetile Q Volume - No GW Prop.", y="")+
#   #geom_text(x=0.034, y=0.059, label="Model C \n R2adj=0.73", cex=6)+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# 
# 
# p4<-ggplot(shoulder_agg, aes(totlm, Q5))+geom_point(col="black")+
#   #geom_text(x=0.033, y=0.059, label="Model D \n R2adj=0.64", cex=6)+
#   labs(x="Modeled 5th Percetile Q Volume - Including GW Prop.", y="")+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# p4
# 
# 
# jpeg("SageLFMPlot.jpeg", width = 530, height = 500)
# 
# p1
# 
# dev.off()
# 
# jpeg("SageSSGWPlots.jpeg", width = 1500, height = 450)
# 
# ggarrange(p2,p3,p4, nrow = 1)
# 
# dev.off()
# 
# #Lookout
# 
# shoulder_agg<-df_list[[2]]
# 
# summary(lm(formula = BFPer ~ PeakSWE + MeltRate + TotPrecip + SWElag +
#          PeakSWEWYDay + ZeroWYDay + SnowFrac + BFlag +Q5lag, data = shoulder_agg))
# 
# v1<-lm(formula = BFPer ~ PeakSWE+
#              Q5lag, data = shoulder_agg)
# 
# summary(v1)
# lm.beta(v1)
# 
# stepAIC(object = v1, direction = "both", trace = TRUE)
# 
# lm_best<-lm(formula = Q5 ~ BFPer + SnowFrac,
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
# shoulder_agg$lm<-lm_best$coefficients[1]+lm_best$coefficients[2]*shoulder_agg$BFPer+
#   lm_best$coefficients[3]*shoulder_agg$SnowFrac
# 
# ggplot(shoulder_agg, aes(lm, Q5))+geom_point()+geom_smooth(method = "lm", se=FALSE, col="black")+
#   theme_classic()+labs(x="SSGW + Snow Fraction", y="Q5")+theme(text = element_text(size = 20))
# 
# 
# # 
# # v1<-lm(formula = BFPerRO01BF99 ~ PeakSWE + MeltRate + TotPrecip + LiquidPrecip + SWElag +
# #           PeakSWEWYDay + ZeroWYDay + SnowFrac + BFlag +Q5lag, data = shoulder_agg)
# # 
# # v2<-lm(formula = BFPerRO01BF99 ~ PeakSWE + MeltRate + TotPrecip + SWElag +
# #          PeakSWEWYDay + ZeroWYDay + SnowFrac + BFlag +Q5lag, data = shoulder_agg)
# # 
# # v3<-lm(formula = BFPerRO01BF99 ~ PeakSWE + MeltRate + SWElag +
# #          PeakSWEWYDay + ZeroWYDay + SnowFrac + BFlag +Q5lag, data = shoulder_agg)
# # 
# # v4<-lm(formula = BFPerRO01BF99 ~ PeakSWE + MeltRate + SWElag +
# #          PeakSWEWYDay + ZeroWYDay + BFlag +Q5lag, data = shoulder_agg)
# # 
# # v5<-lm(formula = BFPerRO01BF99 ~ PeakSWE+MeltRate + SWElag +
# #          PeakSWEWYDay+Q5lag, data = shoulder_agg)
# 
# shoulder_agg$lm<-v5$coefficients[1]+v5$coefficients[2]*shoulder_agg$PeakSWE+
#   v5$coefficients[3]*shoulder_agg$MeltRate+v5$coefficients[4]*shoulder_agg$SWElag+
#   v5$coefficients[5]*shoulder_agg$PeakSWEWYDay+v5$coefficients[6]*shoulder_agg$Q5lag
# 
# BFq5<-lm(Q5~BFPerRO01BF99, data = shoulder_agg)
# 
# summary(BFq5)
# lm.beta(BFq5)
# 
# k1<-lm(formula = Q5 ~ PeakSWE + MeltRate + TotPrecip + LiquidPrecip + SWElag +
#          PeakSWEWYDay + ZeroWYDay + SnowFrac + BFlag +Q5lag, data = shoulder_agg)
# 
# stepAIC(object = k1, direction = "both", trace = TRUE)
# 
# k1<-lm(formula = Q5 ~ SnowFrac+SWElag+
#         Q5lag, data = shoulder_agg)
# 
# k_best<-lm(formula = Q5 ~ MeltRate + TotPrecip + SWElag +
#              PeakSWEWYDay + ZeroWYDay +Q5lag, data = shoulder_agg)
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
# #k1<-lm(formula = Q5 ~ PeakSWE+LiquidPrecip + SWElag +
#          #PeakSWEWYDay + ZeroWYDay +Q5lag, data = shoulder_agg)
# 
# shoulder_agg$BFQ5<-k1$coefficients[1]+k1$coefficients[2]*shoulder_agg$SnowFrac+
#   k1$coefficients[3]*shoulder_agg$SWElag+k1$coefficients[4]*shoulder_agg$Q5lag
# 
# x1<-lm(formula = Q5 ~ BFPerRO01BF99 + 
#          SnowFrac, data = shoulder_agg)
# 
# x1<-lm(formula = Q5 ~ BFPerRO01BF99 + TotPrecip + LiquidPrecip +
#           ZeroWYDay + SnowFrac + BFlag, data = shoulder_agg)
# 
# stepAIC(object=x1, direction = "both", trace = TRUE)
# 
# shoulder_agg$lmtot<-x1$coefficients[1]+x1$coefficients[2]*shoulder_agg$BFPerRO01BF99+
#   x1$coefficients[3]*shoulder_agg$SnowFrac
# 
# 
# p1<-ggplot(shoulder_agg, aes(lm, BFPerRO01BF99))+geom_point(col="black")+
#   labs(y="Shoulder Season GW Proportion", x="Predicted GW Proportion")+
#   #geom_text(x=0.37, y=0.59, label="Model A \n R2adj=0.74", cex=6)+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# 
# 
# p2<-ggplot(shoulder_agg, aes(BFPerRO01BF99, Q5))+geom_point(col="black")+
#   labs(x="Shoulder Season GW Proportion", y="5th Percentile Q Volume")+
#   #geom_text(x=0.39, y=0.39, label="Model B \n R2adj=0.70", cex=6)+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# 
# 
# p3<-ggplot(shoulder_agg, aes(BFQ5, Q5))+geom_point(col="black")+
#   labs(x="Modeled 5th Percetile Q Volume - No GW Prop.", y="")+
#   #geom_text(x=0.2, y=0.385, label="Model C \n R2adj=0.64", cex=6)+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# 
# 
# p4<-ggplot(shoulder_agg, aes(lmtot, Q5))+geom_point(col="black")+
#   #geom_text(x=0.17, y=0.385, label="Model D \n R2adj=0.76", cex=6)+
#   labs(x="Modeled 5th Percetile Q Volume - Including GW Prop.", y="")+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# p4
# 
# jpeg("LookoutLFMPlot.jpeg", width = 530, height = 500)
# 
# p1
# 
# dev.off()
# 
# jpeg("LookoutLMPlots.jpeg", width = 1500, height = 450)
# 
# ggarrange(p2,p3,p4, nrow = 1)
# 
# dev.off()
# 
# #Coal
# 
# v1<-lm(formula = BFPerRO01BF99 ~ MeltRate+LiquidPrecip, 
#          data = shoulder_agg)
# 
# x1<-lm(formula = Q5 ~ PeakSWE+LiquidPrecip, 
#        data = shoulder_agg)
# 
# k1<-lm(formula = Q5 ~ BFPerRO01BF99+PeakSWE, 
#        data = shoulder_agg)
# 
# q1<-lm(formula = Q5 ~ BFPerRO01BF99, 
#        data = shoulder_agg)
# 
# summary(k1)
# 
# lm.beta(k1)
# 
# vif(k1)
# 
# shoulder_agg$lm<-v1$coefficients[1]+v1$coefficients[2]*shoulder_agg$MeltRate+
#   v1$coefficients[3]*shoulder_agg$LiquidPrecip
# 
# shoulder_agg$BFQ5<-x1$coefficients[1]+x1$coefficients[2]*shoulder_agg$PeakSWE+
#   x1$coefficients[3]*shoulder_agg$LiquidPrecip
# 
# shoulder_agg$totlm<-k1$coefficients[1]+k1$coefficients[2]*shoulder_agg$BFPerRO01BF99+
#   k1$coefficients[3]*shoulder_agg$PeakSWE
# 
# 
# p1<-ggplot(shoulder_agg, aes(lm, BFPerRO01BF99))+geom_point(col="black")+
#   labs(y="Shoulder Season GW Proportion", x="Predicted GW Proportion")+
#   #geom_text(x=0.25, y=0.63, label="Model A \n R2adj=0.98", cex=6)+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# p1
# 
# p2<-ggplot(shoulder_agg, aes(BFPerRO01BF99, Q5))+geom_point(col="black")+
#   labs(x="Shoulder Season GW Proportion", y="5th Percentile Q Volume")+
#   #geom_text(x=0.25, y=0.06, label="Model B \n R2adj=0.14", cex=6)+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# 
# 
# p3<-ggplot(shoulder_agg, aes(BFQ5, Q5))+geom_point(col="black")+
#   labs(x="Modeled 5th Percetile Q Volume - No GW Prop.", y="")+
#   #geom_text(x=0.027, y=0.059, label="Model C \n R2adj=0.76", cex=6)+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# 
# 
# p4<-ggplot(shoulder_agg, aes(totlm, Q5))+geom_point(col="black")+
#   #geom_text(x=0.027, y=0.059, label="Model D \n R2adj=0.82", cex=6)+
#   labs(x="Modeled 5th Percetile Q Volume - Including GW Prop.", y="")+
#   geom_smooth(method = "lm", se=F, col="black")+
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                    text = element_text(size=18))
# p4
# 
# 
# jpeg("CoalLFMPlot.jpeg", width = 530, height = 500)
# 
# p1
# 
# dev.off()
# 
# 
# 
# jpeg("CoalSSGWPLots.jpeg", width = 1500, height = 450)
# ggarrange(p2,p3, p4, nrow = 1)
# dev.off()
