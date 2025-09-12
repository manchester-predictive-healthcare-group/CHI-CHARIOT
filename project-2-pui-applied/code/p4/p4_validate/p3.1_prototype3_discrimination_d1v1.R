###
### Program to estimate calibration in 1 development and validation dataset (split sample)
###

### Calibration curve produced seperately for each development/validation dataset combination
### This is primarily because for each development dataset, we get a different set of counter-factual survival times in the validation dataset

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Load survival package
library(survival)
library(foreach)
library(doParallel)
library(doFuture)

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

###
### Write a function that will estimate discrimination for a given development and validation dataset
###
est_discrim_d1v1 <- function(gender, 
                             imps.valid, # Having this as separate argument to speed up computationally
                             age_knots,
                             caltime){
  
  ### Read in model and counter factual survival times
  if (age_knots == 3){
    ### Read in model
    fit <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
    ### Read in counterfactual survival times
    cf.surv.times <- readRDS(paste("data/p4/prototype3_cf_surv_times_", gender,
                                   "_devel", 1,
                                   "_valid", 1, ".rds", sep = ""))
  } else if (age_knots == 4){
    
    if (caltime == TRUE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_caltime_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     ".rds", sep = ""))
    }
    
    if (caltime == FALSE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     ".rds", sep = ""))
    }
  }
  
  ### Pick the first validation dataset
  data.valid <- imps.valid[[1]]
  
  ### Replace cvd_time with counterfactual survival times
  data.valid <- dplyr::select(data.valid, -cvd_time)
  data.valid <- merge(data.valid, cf.surv.times, by.x = "patid", by.y = "patid") |>
    dplyr::rename(cvd_time = cvd_time_cf)
  #   print("data with counterfactual survival times")
  
  ### Create variable with forms that can be used for prediction, and agrees with 'prescribed treatment strategy'
  ### which is what the individual is doing
  data.valid.pred <- 
    dplyr::mutate(data.valid,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  #############################
  ### Assess discrimination ###
  #############################
  
  ### Estimate risks
  ### Note that changing the time, which alter the risks, but not the order of the risks, and therefore
  ### will not effect the estimation of Harrels C
  surv <- as.numeric(est_surv_offset_prototype3(newdata = data.valid.pred, fit = fit, bhaz = bhaz, time = 10*365.25))
  
  ### Estimate C-statistics
  Cstat.object <- intsurv::cIndex(time = data.valid.pred$cvd_time, 
                                  event = data.valid.pred$cvd_indicator, 
                                  risk_score = 1 - surv)
  
  ### Calculate confidence interval and add to object
  Cstat <- Cstat.object[["index"]]
  
  ## Get logit
  logit_Cstat <- log(Cstat/(1-Cstat))
  ## Use delta method to get standard error (sd) of the logit of the C-statistic (which is a proportion)
  logit_se <- sqrt(1/(Cstat*(1-Cstat)*Cstat.object[["comparable"]]))
  ## Get confidence interval by transformation logit_Cstat and upper and lower bounds back onto proper scale
  C_lower <- 1/(1 + exp(-(logit_Cstat - qnorm(0.975)*logit_se)))
  C_upper <- 1/(1 + exp(-(logit_Cstat + qnorm(0.975)*logit_se)))
  
  ### Add to object
  Cstat.object[["index_lower"]] <- C_lower
  Cstat.object[["index_upper"]] <- C_upper
  
  ### Return output object
  return(Cstat.object)
  
}

###
### Write a function that will estimate discrimination for a given development and validation dataset, and a given by variable (must be a factor)
###
est_discrim_d1v1_byvar <- function(gender, 
                                   imps.valid, # Having this as separate argument to speed up computationally
                                   age_knots,
                                   caltime,
                                   byvar){
  
  # gender <- 2
  # byvar <- "age_cat"
  
  ### Read in model and counter factual survival times
  if (age_knots == 3){
    ### Read in model
    fit <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
    ### Read in counterfactual survival times
    cf.surv.times <- readRDS(paste("data/p4/prototype3_cf_surv_times_", gender,
                                   "_devel", 1,
                                   "_valid", 1, ".rds", sep = ""))
  } else if (age_knots == 4){
    
    if (caltime == TRUE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_caltime_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     ".rds", sep = ""))
    }
    
    if (caltime == FALSE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     ".rds", sep = ""))
    }
    
  }
  
  ### Pick the first validation dataset
  data.valid <- imps.valid[[1]]
  
  ### Remove people with missing region variable
  if (byvar == "region"){
    data.valid <- subset(data.valid, !is.na(region))
  }
  
  ### Replace cvd_time with counterfactual survival times
  data.valid <- dplyr::select(data.valid, -cvd_time)
  data.valid <- merge(data.valid, cf.surv.times, by.x = "patid", by.y = "patid") |>
    dplyr::rename(cvd_time = cvd_time_cf)
  
  ### Create variable with forms that can be used for prediction, and agrees with 'prescribed treatment strategy'
  ### which is what the individual is doing
  data.valid.pred <- 
    dplyr::mutate(data.valid,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  ### Subset within subgroups of interest
  bylevels <- levels(data.valid.pred[[byvar]])
  data.valid.pred.list <- lapply(bylevels, function(x) {data.valid.pred[data.valid.pred[,byvar] == x, ]})
  
  #############################
  ### Assess discrimination ###
  #############################
  
  discrim_out <- lapply(1:length(data.valid.pred.list), function(x) {
    ### Estimate risks
    ### Note that changing the time, which alter the risks, but not the order of the risks, and therefore
    ### will not effect the estimation of Harrels C
    surv <- as.numeric(est_surv_offset_prototype3(newdata = data.valid.pred.list[[x]], fit = fit, bhaz = bhaz, time = 10*365.25))
    
    ### Estimate C-statistics
    Cstat.object <- intsurv::cIndex(time = data.valid.pred.list[[x]]$cvd_time, 
                                    event = data.valid.pred.list[[x]]$cvd_indicator, 
                                    risk_score = 1 - surv)
    
    ### Calculate confidence interval and add to object
    Cstat <- Cstat.object[["index"]]
    
    ## Get logit
    logit_Cstat <- log(Cstat/(1-Cstat))
    ## Use delta method to get standard error (sd) of the logit of the C-statistic (which is a proportion)
    logit_se <- sqrt(1/(Cstat*(1-Cstat)*Cstat.object[["comparable"]]))
    ## Get confidence interval by transformation logit_Cstat and upper and lower bounds back onto proper scale
    C_lower <- 1/(1 + exp(-(logit_Cstat - qnorm(0.975)*logit_se)))
    C_upper <- 1/(1 + exp(-(logit_Cstat + qnorm(0.975)*logit_se)))
    
    ### Add to object
    Cstat.object[["index_lower"]] <- C_lower
    Cstat.object[["index_upper"]] <- C_upper
    
    ### Return output object
    return(Cstat.object)
    
  })
  
  ####################################################
  ### Function to extract list of ICI, E50 and E90 ###
  ####################################################
  print("create tables")
  
  ### Extract data
  Cstat <- unlist(lapply(1:length(data.valid.pred.list), function(x){discrim_out[[x]][["index"]]}))
  Cstat_lower <- unlist(lapply(1:length(data.valid.pred.list), function(x){discrim_out[[x]][["index_lower"]]}))
  Cstat_upper <- unlist(lapply(1:length(data.valid.pred.list), function(x){discrim_out[[x]][["index_upper"]]}))
  concordant <- unlist(lapply(1:length(data.valid.pred.list), function(x){discrim_out[[x]][["concordant"]]}))
  comparable <- unlist(lapply(1:length(data.valid.pred.list), function(x){discrim_out[[x]][["comparable"]]}))
  tied_risk <- unlist(lapply(1:length(data.valid.pred.list), function(x){discrim_out[[x]][["tied_risk"]]}))
  
  ### Get N
  N_out <- sapply(data.valid.pred.list, nrow)
  
  ### Create table for ICI, E50 and E90
  discrim_data_table <- data.frame("N" = N_out,
                                   "Cstat" = Cstat,
                                   "Cstat_lower" = Cstat_lower,
                                   "Cstat_upper" = Cstat_upper,
                                   "concordant" = concordant,
                                   "comparable" = comparable,
                                   "tied_risk" = tied_risk)
  rownames(discrim_data_table) <- bylevels
  
  ### Return output object
  return(discrim_data_table)
  
}

###
### Write a function to produce the discrimination table for a given gender, age_knots and caltime, which will run the above two functions, and
### save the output appropriately
###
create_and_save_output <- function(gender, age_knots, caltime){
  
  ### Define gender character version
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  print(paste("age_knots = ", age_knots))
  print(paste("caltime = ", caltime))
  
  ### Read in validation datasets
  imps.valid <- readRDS(paste("data/p4/dfs_valid_", gender_char, ".rds", sep = ""))
  
  ### Create age variable
  imps.valid <- lapply(imps.valid, function(x) {
    x$age_cat <- cut(x$age, breaks = c(17,30,40,50,60,70,80,Inf))
    return(x)
  })
  
  ###
  ### Load the region information
  region_df <- readRDS("../../Aurum_Jun2021_extract/data/extraction/cohort_exclu3.rds") |>
    dplyr::select(patid, region)
  
  ### Read in lookup
  region_lookup <- read.table("../../Aurum_Jun2021_extract/data/unzip/zLookups/region.txt", header = TRUE, sep = "\t") |>
    dplyr::rename(region = regionid, region_desc = Description)
  
  region_df <- dplyr::left_join(region_df, region_lookup, by = dplyr::join_by(region)) |>
    dplyr::select(-region) |>
    dplyr::rename(region = region_desc) |>
    dplyr::mutate(region = as.factor(region))
  
  ### Create region variable
  imps.valid <- lapply(imps.valid, function(x) {
    x <- dplyr::left_join(x, region_df, by = dplyr::join_by(patid))
    return(x)
  })
  
  #####################
  ### Entire cohort ###
  #####################
  
  ### Estimate calibration curves for development=1 and validation=1 datasets
  discrim_d1v1 <- est_discrim_d1v1(gender, imps.valid = imps.valid, age_knots = age_knots, caltime)
  
  saveRDS(discrim_d1v1, paste("data/p4/prototype3_discrim_d1v1_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  print(discrim_d1v1)
  
  ##############
  ### region ###
  ##############
  print("region")
  
  ### Get data for calibration plots
  discrim_data <- est_discrim_d1v1_byvar(gender = gender, 
                                         imps.valid = imps.valid, 
                                         age_knots = age_knots,
                                         caltime = caltime,
                                         byvar = "region")
  
  ### Save
  saveRDS(discrim_data, paste("data/p4/prototype3_discrim_table_d1v1_region_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  print(discrim_data)
  
  #################
  ### Ethnicity ###
  #################
  print("Ethnicity")
  
  ### Get data for calibration plots
  discrim_data <- est_discrim_d1v1_byvar(gender = gender, 
                                         imps.valid = imps.valid, 
                                         age_knots = age_knots,
                                         caltime = caltime,
                                         byvar = "ethnicity")
  
  ### Save
  saveRDS(discrim_data, paste("data/p4/prototype3_discrim_table_d1v1_ethnicity_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  print(discrim_data)
  
  ###########
  ### Age ###
  ###########
  print("AGE")
  
  ### Get data for calibration plots
  discrim_data <- est_discrim_d1v1_byvar(gender = gender, 
                                         imps.valid = imps.valid, 
                                         age_knots = age_knots,
                                         caltime = caltime,
                                         byvar = "age_cat")
  
  ### Save
  saveRDS(discrim_data, paste("data/p4/prototype3_discrim_table_d1v1_age_cat_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  print(discrim_data)
  
}

### Run this function
for (gender_in in c(1,2)){
  
  create_and_save_output(gender = gender_in, age_knots = 3, caltime = FALSE) 
  create_and_save_output(gender = gender_in, age_knots = 4, caltime = FALSE) 
  create_and_save_output(gender = gender_in, age_knots = 4, caltime = TRUE) 
  
}

print(paste("FINISHED", Sys.time()))


