###
### Program to estimate calibration.
### The focus will be a temporal validation of the initial risk estimation layer.
### i.e. we are evaluating performance of the model when the index date is defined as 1/2/3/4/5 years after
### the start of follow-up.
###
### We have not imputed these cohorts separately. We have re-extracted the data, but missing data is imputed using the
### data imputed at baseline.
###
### Predicted risks will be under the treatment strategy of 'continue on current treatment strategy'. Counterfactual survival
### times have also been calculated under this treatment strategy.
### 

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
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
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_validate_wrapper <- function(gender, 
                                 t_fup_int, # number of years defining the follow-up index date
                                 t_eval_int, # number of years prediction to evaluate performance at
                                 plot.range = NULL){
  
  print(paste("gender =", gender, "t_fup_int =", t_fup_int, Sys.time()))
  
  ### Define t_eval in days
  t_eval <- t_eval_int*365.25
  
  ### Read in model
  fit <- readRDS(paste("data/p5/model_cox_", gender, "_imp", 1, ".rds", sep = ""))
  bhaz <- readRDS(paste("data/p5/model_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
  
  ### Read in counterfactual survival times
  cf_surv_times <- readRDS(paste0("data/p5/temporalv_cf_surv_times_", gender, "_t", t_fup_int, ".rds"))
  
  ### Read in validation dataset
  data_valid <- readRDS(paste0("data/p5/df_valid_temporalv_", gender, "_", t_fup_int, ".rds"))
  
  ### Replace cvd_time with counterfactual survival times
  data_valid <- dplyr::select(data_valid, -cvd_time, -cvd_indicator)
  data_valid <- merge(data_valid, cf_surv_times, by.x = "patid", by.y = "patid") |>
    dplyr::rename(cvd_time = cvd_time_cf, cvd_indicator = cvd_indicator_cf)
  #   print("data with counterfactual survival times")
  
  ### Create variable with forms that can be used for prediction, and agrees with 'prescribed treatment strategy'
  ### which is what the individual is doing
  data_valid_pred <- 
    dplyr::mutate(data_valid,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  ### Estimate survival probabilities
  data_valid_pred$surv <- as.numeric(est_surv_offset_prototype3(newdata = data_valid_pred, fit = fit, bhaz = bhaz, time = t_eval))
  
  ##########################
  ### Assess calibration ###
  ##########################
  
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
  saveRDS(calib_ph, paste0("data/p5/temporalv_calib_ph_", gender, "_t", t_fup_int, "_teval", t_eval_int,".rds"))
  saveRDS(calib_km_group, paste0("data/p5/temporalv_calib_km_group_", gender, "_t", t_fup_int, "_teval", t_eval_int,".rds"))
  
  #############################
  ### Assess discrimination ###
  #############################
  
  ### Estimate C-statistics
  Cstat_object <- survival::concordance(survival::Surv(cvd_time, cvd_indicator) ~ surv, 
                                        data = data_valid_pred,
                                        timewt = "n",
                                        influence = 0)
  
  ### Save
  saveRDS(Cstat_object, paste0("data/p5/temporalv_discrim_", gender, "_t", t_fup_int, "_teval", t_eval_int,".rds"))
  
  NULL
}

### Run function
lapply(c(1,2), function(x){
  lapply(1:5, function(y){
    est_validate_wrapper(gender = x, t_fup_int = y, t_eval_int = 5)
    est_validate_wrapper(gender = x, t_fup_int = y, t_eval_int = 10)
  })
})

print(paste("FINISHED", Sys.time()))