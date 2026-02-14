# 2 EMMA baseflow correlation with peak swe.
# Inputs: output_dir/HJA_Ave_StorageMetrics_CatCharacter.csv.
# Author: Sidney Bush
# Date: 2026-02-13

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
  for (k in 1:length(site_data)) {
    
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

### This is creating a matrix of the Beta value so that we can create a plot from them. 
# Define row and column names
row_names <- c("Ash_Per", "Harvest", "Landslide_Total", "Landslide_Young", 
               "Lava1_per", "Lava2_per", "Slope_mean")

col_names <- c("recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", 
               "CV_Q5norm_mean", "mean_bf_mean", "fdc_slope_mean", 
               "S_annual_mm_mean", "DR_Overall", "JST_AT_mean")

# Initialize an empty matrix
beta_matrix <- matrix(NA, nrow = length(row_names), ncol = length(col_names),
                      dimnames = list(row_names, col_names))

# Fill the matrix using df_final
for (i in seq_along(row_names)) {
  for (j in seq_along(col_names)) {
    # Filter df_final for matching Variable and Prefix
    match_row <- final_df[final_df$Variable == row_names[i] & final_df$Prefix == col_names[j], ]
    
    # If a match is found, assign the Beta value to the matrix
    if (nrow(match_row) == 1) {
      beta_matrix[row_names[i], col_names[j]] <- match_row$Beta
    }
  }
}

# View the matrix


#### At this point we are only using the variables that had a VIF less than 10 to re-run the models 


# Open lists for coefficient and r2 values to be appended into
coef_df_list <- list()
r2_df_list <- list()

# Cycle through outcome variables to come to a final model
for (k in 1:length(site_data)) {
  
  # Get the current prefix (outcome variable name)
  current_prefix <- outcome_vars[k]
  
  # Filter final_df to get variables associated with this prefix
  relevant_vars <- final_df$Variable[final_df$Prefix == current_prefix]
  
  # Subset site_data[[k]] to include only relevant variables
  site_subset <- site_data[[k]][, relevant_vars, drop = FALSE]
  
  # Add the outcome variable (Prefix) column to site_subset
  site_subset[[current_prefix]] <- site_data[[k]][[current_prefix]]
  
  # Create model formula
  lm_mod_vars <- formula(paste(current_prefix, "~", paste(relevant_vars, collapse = "+")))
  
  # Fit model with NA removed
  lm_mod <- lm(lm_mod_vars, data = na.omit(site_subset))
  lm_sum <- summary(lm_mod)
  
  # Extract coefficients
  coefs <- as.data.frame(lm_sum$coefficients)
  
  # Get beta coefficients
  beta_df <- as.data.frame(lm.beta(lm_mod))
  beta_df <- rownames_to_column(beta_df)
  names(beta_df)[1] <- "Row.names"
  
  
  # Merge results
  df_master <- merge(coefs, beta_df, by.x = "row.names", by.y = "Row.names", all.x = TRUE)
  
  # Rename the merged row name column to a consistent name
  names(df_master)[1] <- "Variable"
  
  # Set row names and clean up
  df_master <- df_master %>% column_to_rownames(var = "Variable")
  df_master <- df_master[, c("Pr(>|t|)", "lm.beta(lm_mod)")]
  
  # Rename columns
  colnames(df_master) <- c(paste0(current_prefix, "_pvalue"),
                           paste0(current_prefix, "_beta"))
  
  
  # Get R² and adjusted R²
  r2 <- lm_sum$r.squared
  adj_r2 <- lm_sum$adj.r.squared
  r2_tot <- rbind(r2, adj_r2)
  colnames(r2_tot) <- paste0(current_prefix)
  
  # Append results
  coef_df_list[[k]] <- df_master
  r2_df_list[[k]] <- r2_tot
}


#### Puts results back into a dataframe


# Initialize an empty data frame to store results
final_df <- data.frame(Prefix = character(), Variable = character(), Beta = numeric(), pvalue = numeric(), stringsAsFactors = FALSE)

# Loop through each dataframe in the list of coefficents an put it into a daframe, at this point we have only kept vairables with a VIF that was less than = 10. 
for (df in coef_df_list) {
  # Identify beta and pvalue columns
  beta_cols <- grep("_beta$", names(df), value = TRUE)
  pval_cols <- grep("_pvalue$", names(df), value = TRUE)
  
  # Loop through each prefix found in beta columns
  for (beta_col in beta_cols) {
    prefix <- sub("_beta$", "", beta_col)
    pval_col <- paste0(prefix, "_pvalue")
    
    # Check if corresponding pvalue column exists
    if (pval_col %in% names(df)) {
      # Extract row names and corresponding beta and pvalue values
      variables <- rownames(df)
      beta_vals <- df[[beta_col]]
      pval_vals <- df[[pval_col]]
      
      # Create a temporary data frame
      temp_df <- data.frame(
        Prefix = prefix,
        Variable = variables,
        Beta = beta_vals,
        pvalue = pval_vals,
        stringsAsFactors = FALSE
      )
      
      # Append to the final data frame
      final_df <- rbind(final_df, temp_df)
    }
  }
}

#Removing intercept 
final_df <- final_df[final_df$Variable != "(Intercept)", ]



### This is creating a matrix of the Beta value so that we can create a plot from them. 
# Define row and column names
row_names <- c("Ash_Per", "Harvest", "Landslide_Total", "Landslide_Young", 
               "Lava1_per", "Lava2_per", "Slope_mean")

col_names <- c("recession_curve_slope_mean", "RBI_mean", "Q5norm_mean", 
               "CV_Q5norm_mean", "mean_bf_mean", "fdc_slope_mean", 
               "S_annual_mm_mean", "DR_Overall", "JST_AT_mean")

# Initialize an empty matrix
beta_matrix <- matrix(NA, nrow = length(row_names), ncol = length(col_names),
                      dimnames = list(row_names, col_names))

# Fill the matrix using df_final
for (i in seq_along(row_names)) {
  for (j in seq_along(col_names)) {
    # Filter df_final for matching Variable and Prefix
    match_row <- final_df[final_df$Variable == row_names[i] & final_df$Prefix == col_names[j], ]
    
    # If a match is found, assign the Beta value to the matrix
    if (nrow(match_row) == 1) {
      beta_matrix[row_names[i], col_names[j]] <- match_row$Beta
    }
  }
}

# View the matrix
#plot

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


setwd("/Users/pamelasullivan/Box Sync/2024_NSF-Wildfire_WaterCycle/05_Storage_Manuscript/05_Outputs/Hydrometric")

pdf("Beta_Weight_Plot_HJA2.pdf", width = 9, height = 7,family = "Times")
#increase margin
par(mar=c(11,10,3,6)+.1)

# Remove "_mean" from column names
colnames(beta_matrix) <- sub("_mean$", "", colnames(beta_matrix))

adj_r2_values <- sapply(r2_df_list, function(df) round(df[2, 1], 2))

plot(beta_matrix, col = rampcols, breaks = rampbreaks, border=NA, las=2, ann=FALSE,  digits = 2, 
     na.print = FALSE, polygon.key = NULL, fmt.key="%.0f", cex=1.0)


text(x = 1:ncol(beta_matrix), y = nrow(beta_matrix) + 0.35, labels = adj_r2_values , cex = 0.8,font = 2)

dev.off()
