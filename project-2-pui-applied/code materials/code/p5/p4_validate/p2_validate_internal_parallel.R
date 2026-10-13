###
### Internal validation. Written to run in parallel over imputation chains.
### Calibration and discrimination are assessed, in entire cohort, and within subgroups defined
### by protected characteristics.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Load packages
library(survival)
library(doParallel)
library(foreach)

### Define time to evaluate performance at (10-years)
t_eval <- round(10*365.25)

### Write function to run validation on a single imputation run
est_calibration_one_imputation <- function(gender, imp_num){
  
  ### Prelim
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  print(paste("imp_num = ", imp_num, Sys.time()))
  
  ### Pick the imputed dataset
  data_valid <- readRDS(paste("data/p5/dfs_all_", gender_char, imp_num, ".rds", sep = ""))
  
  ### Read in model fit
  fit <- readRDS(paste("data/p5/model_cox_", gender, "_imp", imp_num, ".rds", sep = ""))
  bhaz <- readRDS(paste("data/p5/model_cox_bhaz_", gender, "_imp", imp_num, ".rds", sep = ""))
  
  ### Read in counterfactual survival times
  cf_surv_times <- readRDS(paste("data/p5/cf_surv_times_", gender, "_imp", imp_num, ".rds", sep = ""))
  stopifnot(!anyDuplicated(cf_surv_times$patid))
  
  ### Replace cvd_time and cvd_indicator with counterfactual survival times
  n_before <- nrow(data_valid)
  data_valid <- dplyr::select(data_valid, -cvd_time, -cvd_indicator) |>
    dplyr::left_join(cf_surv_times, by = "patid") |>
    dplyr::rename(cvd_time = cvd_time_cf, cvd_indicator = cvd_indicator_cf)
  stopifnot(nrow(data_valid) == n_before)
  stopifnot(!anyNA(data_valid$cvd_time))
  stopifnot(!anyNA(data_valid$cvd_indicator))
  
  ### Create variable with forms that can be used for prediction, and agrees with 'prescribed treatment strategy'
  ### which is what the individual is doing
  data_valid_pred <- 
    dplyr::mutate(data_valid,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  ### Ensure data_valid_pred is sorted in the correct order
  data_valid_pred <- dplyr::arrange(data_valid_pred, patid)
  
  ### Estimate survival probabilities
  data_valid_pred$surv <- as.numeric(est_surv_offset_prototype3(newdata = data_valid_pred, fit = fit, bhaz = bhaz, time = t_eval))
  
  ####################################
  ### Calibration in entire cohort ###
  ####################################
  
  ### PH regression approach
  print(paste("PH", Sys.time()))
  calib_ph <- est_calib_ph(data = data_valid_pred, 
                           fit = fit, 
                           bhaz = bhaz, 
                           time = t_eval,
                           surv = data_valid_pred$surv)
  
  ### KM grouped approach
  print(paste("KM group", Sys.time()))
  calib_km_group <- est_calib_plot_group(data = data_valid_pred, 
                                         fit = fit, 
                                         bhaz = bhaz,  
                                         time = t_eval,
                                         n.groups = 100,
                                         surv = data_valid_pred$surv,
                                         single.CI = TRUE)
  
  ### Save
  saveRDS(calib_ph, paste("data/p5/calib_ph_", gender, "_", imp_num, ".rds", sep = ""))
  saveRDS(calib_km_group, paste("data/p5/calib_km_group_", gender, "_", imp_num, ".rds", sep = ""))
  
  ####################################
  ### Calibration within subgroups ###
  ####################################
  
  ### 
  ### Create variables for assessing performance within subgroups
  ###
  
  ### Create age variable categorised
  data_valid_pred$age_cat <- cut(data_valid_pred$age, breaks = c(17,30,40,50,60,70,80,Inf))
  
  ###
  ### Load the region information
  region_df <- readRDS("../../Aurum_Jun2021_extract/data/extraction/cohort_exclu3.rds") |>
    dplyr::select(patid, region)
  stopifnot(!anyDuplicated(region_df$patid))
  
  ### Read in lookup
  region_lookup <- read.table("../../Aurum_Jun2021_extract/data/unzip/zLookups/region.txt", header = TRUE, sep = "\t") |>
    dplyr::rename(region = regionid, region_desc = Description)
  
  region_df <- dplyr::left_join(region_df, region_lookup, by = dplyr::join_by(region)) |>
    dplyr::select(-region) |>
    dplyr::rename(region = region_desc) |>
    dplyr::mutate(region = as.factor(region))
  
  ### Create region variable
  data_valid_pred <- dplyr::left_join(data_valid_pred, region_df, by = dplyr::join_by(patid))
  
  ###
  ### Write a function that will estimate calibration within subgroups defined by a variable (must be a factor)
  ###
  est_calib_ph_wrapper_byvar <- function(data_valid_pred, byvar){
    
    print(paste("gender = ", gender, byvar, Sys.time()))
    
    ### Remove people with missing region variable
    if (byvar == "region"){
      data_valid_pred <- subset(data_valid_pred, !is.na(region))
    }
    
    ### Subset within subgroups of interest
    bylevels <- levels(data_valid_pred[[byvar]])
    data_valid_pred_list <- lapply(bylevels, function(x) {data_valid_pred[data_valid_pred[,byvar] == x, ]})
    
    ##########################
    ### Assess calibration ###
    ##########################
    
    ### PH regression approach
    calib_ph <- lapply(1:length(data_valid_pred_list), function(x) {
      est_calib_ph(data = data_valid_pred_list[[x]], 
                   fit = fit, 
                   bhaz = bhaz, 
                   time = t_eval,
                   surv = data_valid_pred_list[[x]]$surv,
                   plot = FALSE)
    })
    
    calib_km_group <- lapply(1:length(data_valid_pred_list), function(x) {
      est_calib_plot_group(data = data_valid_pred_list[[x]], 
                           fit = fit, 
                           bhaz = bhaz,  
                           time = t_eval,
                           n.groups = 10,
                           surv = data_valid_pred_list[[x]]$surv)
    })
    
    ### Save
    saveRDS(calib_ph, paste0("data/p5/calib_ph_", gender, "_", imp_num, "_", byvar, ".rds", sep = ""))
    saveRDS(calib_km_group, paste0("data/p5/calib_km_group_", gender, "_", imp_num, "_", byvar, ".rds", sep = ""))
    
    NULL
    
  }
  
  ### Run this function
  est_calib_ph_wrapper_byvar(data_valid_pred, byvar = "region")
  est_calib_ph_wrapper_byvar(data_valid_pred, byvar = "age_cat")
  est_calib_ph_wrapper_byvar(data_valid_pred, byvar = "ethnicity")
  print(paste("FINISHED CALIB"))
  
  ######################
  ### Discrimination ###
  ######################
  
  ########################
  ### In entire cohort ###
  ########################
  
  ### Estimate C-statistics
  Cstat_object <- survival::concordance(survival::Surv(cvd_time, cvd_indicator) ~ surv, 
                                        data = data_valid_pred,
                                        timewt = "n",
                                        influence = 0)
  
  ### Save
  saveRDS(Cstat_object, paste0("data/p5/discrim_", gender, "_", imp_num, ".rds"))
  
  ########################
  ### Within subgroups ###
  ########################
  
  ###
  ### Write a function that will estimate discrimination within subgroups defined by a variable (must be a factor)
  ### Variable type defines whether validation occurs in bootstrapped dataset, or original dataset
  ###
  est_discrim_wrapper_byvar <- function(data_valid_pred,
                                        byvar){
    
    ### Remove people with missing region variable
    if (byvar == "region"){
      data_valid_pred <- subset(data_valid_pred, !is.na(region))
    }
    
    ### Subset within subgroups of interest
    bylevels <- levels(data_valid_pred[[byvar]])
    data_valid_pred_list <- lapply(bylevels, function(x) {data_valid_pred[data_valid_pred[,byvar] == x, ]})
    
    #############################
    ### Assess discrimination ###
    #############################
    discrim_out <- lapply(1:length(data_valid_pred_list), function(x) {
      
      ### Estimate C-statistics
      Cstat.object <- survival::concordance(survival::Surv(cvd_time, cvd_indicator) ~ surv, 
                                            data = data_valid_pred_list[[x]],
                                            timewt = "n",
                                            influence = 0)
      
      ### Extract C-stat
      Cstat <- unname(Cstat.object$concordance)
      
      ### Return
      return(Cstat)
      
    })
    
    ### Extract data in vector
    Cstat <- unlist(discrim_out, use.names = FALSE)
    
    ### Get N
    N_out <- sapply(data_valid_pred_list, nrow)
    
    ### Create table
    discrim_data_table <- data.frame("N" = N_out,
                                     "Cstat" = Cstat,
                                     levels = bylevels)
    
    ### Save
    saveRDS(discrim_data_table, paste0("data/p5/discrim_", gender, "_", imp_num, "_", byvar, ".rds"))
    
    NULL
    
  }
  
  ### Run this function
  est_discrim_wrapper_byvar(data_valid_pred, byvar = "region")
  est_discrim_wrapper_byvar(data_valid_pred, byvar = "age_cat")
  est_discrim_wrapper_byvar(data_valid_pred, byvar = "ethnicity")
  print(paste("FINISHED DISCRIM"))
  
}

###
### Run the validation for a given gender
###
gender_wrapper <- function(gender){
  
  ### Set up parallel backend
  n_cores <- min(10L, parallel::detectCores())
  print(paste("Parallel workers =", n_cores))
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  ### Run validation in parallel over imputations
  foreach::foreach(
    imp_num = 1:10,
    .packages = c("dplyr", "survival", "rms", "ggplot2", "ggExtra"),
    .export = c(
      "est_calibration_one_imputation",
      "est_surv",
      "est_surv_offset_prototype3",
      "est_calib_ph",
      "est_calib_plot_group",
      "gender",
      "t_eval"
    )
  ) %dopar% {
    est_calibration_one_imputation(gender, imp_num)
  }
  
  ### Stop parallel workers
  parallel::stopCluster(cl)
  
  print(paste("FINISHED", c("male", "female")[gender], Sys.time()))
  
}

### Run for both genders
gender_wrapper(2)
gender_wrapper(1)
print(paste("FINISHED ALL", Sys.time()))
