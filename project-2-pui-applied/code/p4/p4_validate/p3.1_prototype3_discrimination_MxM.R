###
### Program to estimate calibration
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
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_discrim_MxM <- function(gender, 
                            age_knots,
                            caltime,
                            devel.num,
                            valid.num,
                            imps.valid){
  
  
  ### Read in model and counter factual survival times
  if (age_knots == 3){
    ### Read in model
    fit <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", devel.num, ".rds", sep = ""))
    bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", devel.num, ".rds", sep = ""))
    
    ### Read in counterfactual survival times
    cf.surv.times <- readRDS(paste("data/p4/prototype3_cf_surv_times_", gender,
                                   "_devel", devel.num,
                                   "_valid", valid.num, ".rds", sep = ""))
  } else if (age_knots == 4){
    
    if (caltime == TRUE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", devel.num, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", devel.num, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_caltime_cf_surv_times_", gender, 
                                     "_devel", devel.num, 
                                     "_valid", valid.num, 
                                     ".rds", sep = ""))
    }
    
    if (caltime == FALSE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", devel.num, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", devel.num, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_cf_surv_times_", gender, 
                                     "_devel", devel.num, 
                                     "_valid", valid.num, 
                                     ".rds", sep = ""))
    }
  }
  
  ### Pick one validation dataset
  data.valid <- imps.valid[[valid.num]]
  
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
  ### Note that changing the time, which alter the risks, but not hte order of the risks, and therefore
  ### will not effect the estimation of Harrels C
  surv <- as.numeric(est_surv_offset_prototype3(newdata = data.valid.pred, fit = fit, bhaz = bhaz, time = 10*365.25))
  
  ### Estimate C-statistics
  Cstat <- intsurv::cIndex(time = data.valid.pred$cvd_time, 
                           event = data.valid.pred$cvd_indicator, 
                           risk_score = 1 - surv)
  
  ### Return output object
  return(Cstat)
  
}

###
### Write a function to produce the calibration plots for a given gender, age_knots and caltime, which will run the above two functions, and
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
  
  ### Apply this functions over 10 development and validation datasets
  cl <- parallel::makeCluster(11)
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers()
  discrim_MxM <- (foreach(x = 1:10, .combine = list, .multicombine = TRUE, .packages = c("survival"), .export = c("est_discrim_MxM", "est_calib_plot", "est_calib_plot_group", "est_surv_offset_prototype3", "est_surv")) %dopar% {
    lapply(1:10, function(y) {est_discrim_MxM(gender, age_knots, caltime, x, y, imps.valid = imps.valid)})
  })
  stopCluster(cl)
  
  print(paste("DATA OBTAINED", Sys.time()))
  ### Extract the data
  discrim_comb <- unlist(discrim_MxM, recursive = FALSE) |>
    lapply(function(x) {x[["index"]]}) |>
    unlist()
  
  ### Create a table with mean, median, inter quartile range
  dp <- 2
  discrim_table <- data.frame("mean (sd)" = paste(round(mean(discrim_comb), dp), 
                                                  " (", 
                                                  round(sd(discrim_comb), dp),
                                                  ")", sep = ""),
                              "median (p025 p975)" = paste(round(median(discrim_comb), dp), 
                                                           " (", 
                                                           round(quantile(discrim_comb, p = 0.025), dp),
                                                           ",",
                                                           round(quantile(discrim_comb, p = 0.975), dp),
                                                           ")", sep = ""),
                              min = round(min(discrim_comb), dp),
                              max = round(max(discrim_comb), dp))
  
  saveRDS(discrim_table, paste("data/p4/prototype3_discrim_table_MxM", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  saveRDS(discrim_comb, paste("data/p4/prototype3_discrim_list_MxM", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  print(paste("FINISHED", gender_char, Sys.time()))
  
}

### Run this function
for (gender_in in c(1,2)){

  create_and_save_output(gender = gender_in, age_knots = 4, caltime = TRUE) 
  
}

print(paste("FINISHED", Sys.time()))

