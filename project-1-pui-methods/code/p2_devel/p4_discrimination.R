###
### Program to estimate discrimination in the validation cohort
###

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

### Define HR for offsets
lnHR_statins <- readRDS("data/offsets_total_lnHR_statins.rds")
lnHR_ah <- readRDS("data/offsets_total_lnHR_ah.rds")
lnHR_sbp <- readRDS("data/offsets_direct_lnHR_sbp.rds")
lnHR_bmi <- readRDS("data/offsets_direct_lnHR_bmi.rds")
lnHR_nonhdl <- readRDS("data/offsets_direct_lnHR_nonhdl.rds")

###
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_discrimination <- function(gender, model){
  
  # gender <- 2
  # plot.range = NULL
  
  ### Read in model
  fit <- readRDS(paste("data/fit_", gender, "_model", model, ".rds", sep = ""))
  bhaz <- readRDS(paste("data/bhaz_", gender, "_model", model, ".rds", sep = ""))
  
  ### Read in validation data
  df_valid <- readRDS(paste("data/df_imp_valid_", gender, sep = ""))
  
  if (model != 0){
    
    ### Read in counterfactual survival times
    cf_surv_times <- readRDS(paste("data/cf_surv_times_", gender, "_model", model, ".rds", sep = ""))
    
    ### Replace cvd_time with counterfactual survival times
    df_valid <- dplyr::select(df_valid, -cvd_time)
    df_valid <- merge(df_valid, cf_surv_times, by.x = "patid", by.y = "patid") |>
      dplyr::rename(cvd_time = cvd_time_cf)
    #   print("data with counterfactual survival times")
  }
  
  ###
  ### For model model 3, apply an offset for differences in sbp, nonhdl and BMI.
  ### We create adjusted versions of SBP, nonhdl and BMI to be relative to some baseline, to apply the offsets
  ### These are also capped:
  ### Assume no benefit for SBP below 120, so it is capped there
  ### Assume no benefit for BMI below 25, so it is capped there
  ###
  
  ### SBP relative to 120 (no benefit for being lower than 120)
  ### BMI relative to 25 (no benefit for being lower than 25)
  ### nonhdl relative to 4
  df_valid <- dplyr::mutate(df_valid,
                            sbp_adj = dplyr::case_when(sbp < 120 ~ 0,
                                                       TRUE ~ sbp - 120),
                            bmi_adj = dplyr::case_when(bmi < 25 ~ 0,
                                                       TRUE ~ bmi - 25),
                            nonhdl_adj = dplyr::case_when(nonhdl < 2.6 ~ 0,
                                                          TRUE ~ nonhdl - 2.6)
  )
  
  ### Create offset variable with appropriate names
  df_valid_pred <- 
    dplyr::mutate(df_valid,
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_bmi_lnHR = lnHR_bmi*bmi_adj, 
                  offset_nonhdl_lnHR = lnHR_nonhdl*nonhdl_adj,
                  offset_sbp_lnHR = lnHR_sbp*sbp_adj)
  
  #############################
  ### Assess discrimination ###
  #############################
  
  ### Estimate risks
  ### Note that changing the time, which alter the risks, but not the order of the risks, and therefore
  ### will not effect the estimation of Harrels C
  surv <- as.numeric(est_surv_offset(newdata = df_valid_pred, fit = fit, bhaz = bhaz, time = 10*365.25))
  
  ### Estimate C-statistics
  Cstat_object <- intsurv::cIndex(time = df_valid_pred$cvd_time, 
                                  event = df_valid_pred$cvd_indicator, 
                                  risk_score = 1 - surv)
  
  ### Calculate confidence interval and add to object
  Cstat <- Cstat_object[["index"]]
  
  ## Get logit
  logit_Cstat <- log(Cstat/(1-Cstat))
  ## Use delta method to get standard error (sd) of the logit of the C-statistic (which is a proportion)
  logit_se <- sqrt(1/(Cstat*(1-Cstat)*Cstat_object[["comparable"]]))
  ## Get confidence interval by transformation logit_Cstat and upper and lower bounds back onto proper scale
  C_lower <- 1/(1 + exp(-(logit_Cstat - qnorm(0.975)*logit_se)))
  C_upper <- 1/(1 + exp(-(logit_Cstat + qnorm(0.975)*logit_se)))
  
  ### Add to object
  Cstat_object[["index_lower"]] <- C_lower
  Cstat_object[["index_upper"]] <- C_upper
  
  ### Save output
  saveRDS(Cstat_object, paste("data/discrim_", gender, "_model", model, ".rds", sep = ""))
  
}

# ### Run this function
# for (gender_in in c(1,2)){
#   
#   print(paste("gender = ", gender_in))
#   
#   for (model_in in c(1,2,3,4,5,6,7)){
#     
#     print(paste("model = ", model_in))
#     est_discrimination(gender = gender_in, model = model_in) 
#     
#   }
# }

### Run this function
model_in <- 0
for (gender_in in c(1,2)){
  
  print(paste("gender = ", gender_in))
  
  print(paste("model = ", model_in))
  est_discrimination(gender = gender_in, model = model_in) 
  
}

print(paste("FINISHED", Sys.time()))

