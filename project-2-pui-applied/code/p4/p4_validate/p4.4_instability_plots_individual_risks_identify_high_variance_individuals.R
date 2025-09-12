###
### Progream to explore which individuals have lots of variability in risk
###

###
### Read in bootstrapped models and create instability plot
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)
library(ggplot2)
library(survival)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- as.numeric(180)

### Set seed
set.seed(101)

### Assign inputs
knots4 = FALSE
n_pat = 3000
n_model = 500
  
### function to create and save data
get_instability_data <- function(gender){
  
  ### Define gender_char
  gender_char <- c("male", "female")[gender]
  
  ### Extract model fits
  if (knots4 == FALSE){
    
    ### Read in models, offsets and bhaz's
    fit_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_stability_cox_", gender, "_samp", x, ".rds", sep = ""))})
    bhaz_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_stability_cox_bhaz_", gender, "_samp", x, ".rds", sep = ""))})
    offset_means_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_stability_offset_means_", gender, "_samp", x, ".rds", sep = ""))})
    
    ### Read in prototype3
    fit_prototype3 <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz_prototype3 <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
  } else if (knots4 == TRUE) {
    
    ### Read in models, offsets and bhaz's
    fit_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_4knots_stability_cox_", gender, "_samp", x, ".rds", sep = ""))})
    bhaz_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_4knots_stability_cox_bhaz_", gender, "_samp", x, ".rds", sep = ""))})
    offset_means_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_4knots_stability_offset_means_", gender, "_samp", x, ".rds", sep = ""))})
    
    ### Read in prototype3
    fit_prototype3 <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz_prototype3 <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
  }
  
  ### Read in validation dataset
  data_valid <- readRDS(paste("data/p4/dfs_valid_", gender_char, ".rds", sep = ""))[[1]]
  
  ### Pick n_pat patients at random
  set.seed(101)
  data_valid_subset <- data_valid[sample(1:nrow(data_valid), n_pat, replace = FALSE), ]
  
  ### Add values for offset variables/treatment strategy
  data_valid_subset <- 
    dplyr::mutate(data_valid_subset,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  
  ### Estimate risks according to the non-bootstrapped model of prototype3
  data_valid_subset$pred_10y <- 1 - est_surv_offset_prototype3(newdata = data_valid_subset, 
                                                               fit = fit_prototype3, 
                                                               bhaz = bhaz_prototype3, 
                                                               time = 10*365.25)
  
  ### Write function to estimate risks
  estimate_risks <- function(t_eval){
    
    pred_list <- lapply(1:n_model, function(x){
      
      ### Get the lp
      lp <- predict(fit_list[[x]], newdata = data_valid_subset, reference = "sample") - sum(offset_means_list[[x]])
      
      ### Get the survival probability
      surv <- as.numeric(exp(-exp(lp)*bhaz_list[[x]]$hazard[max(which(bhaz_list[[x]]$time <= t_eval))]))
      pred <- 1-surv
      
      return(pred)
    })
    
    pred_list <- do.call("cbind", pred_list)
    colnames(pred_list) <- paste("pred", 1:n_model, sep = "")
    
    return(pred_list) 
    
  }
  
  ### Estiamte risks
  pred_10y <- estimate_risks(10*365.25)
  # pred_5y <- estimate_risks(5*365.25)
  
  ### Add to data frame
  data_valid_subset <- cbind(data_valid_subset, pred_10y)
  
  ### Create data frame from the predicted probabilities
  data_valid_subset_long <- data_valid_subset |>
    dplyr::select(patid, pred_10y, paste("pred",1:n_model, sep = "")) |>
    tidyr::pivot_longer(cols = paste("pred",1:n_model, sep = ""), names_to = "sample", values_to = "pred")
  
  ### Save this for working with
  saveRDS(data_valid_subset_long, paste("data/p4/instability_plots_data", gender, ".rds", sep = ""))
  
}
  
### Get the instability data
get_instability_data(1)
get_instability_data(2)

###
### Look at male data
###

### Read in male first, as this had more variation
data_valid_subset_long <- readRDS(paste("data/p4/instability_plots_data", 1, ".rds", sep = ""))

### Get the range for each individual
data_valid_range <- data_valid_subset_long |>
  dplyr::group_by(patid) |>
  subset(pred_10y < 0.1) |>
  dplyr::summarise(range = max(pred) - min(pred),
                   prototype3_risk = mean(pred_10y)) |>
  dplyr::arrange(dplyr::desc(range))

### Get the top five individuals
data_valid <- readRDS(paste("data/p4/dfs_valid_", "male", ".rds", sep = ""))[[1]]
data_valid_interest <- subset(data_valid, patid == 6603329221578)
data_valid_interest <- subset(data_valid, patid %in% data_valid_range$patid[1:5])


###
### Look at female data
###

### Read in male first, as this had more variation
data_valid_subset_long <- readRDS(paste("data/p4/instability_plots_data", 2, ".rds", sep = ""))

### Get the range for each individual
data_valid_range <- data_valid_subset_long |>
  dplyr::group_by(patid) |>
  subset(pred_10y < 0.3) |>
  dplyr::summarise(range = max(pred) - min(pred),
                   prototype3_risk = mean(pred_10y)) |>
  dplyr::arrange(dplyr::desc(range))

### Get the top five individuals
data_valid <- readRDS(paste("data/p4/dfs_valid_", "female", ".rds", sep = ""))[[1]]
data_valid_interest <- subset(data_valid, patid == 2158480020003)
data_valid_interest <- subset(data_valid, patid %in% data_valid_range$patid[1:5])
