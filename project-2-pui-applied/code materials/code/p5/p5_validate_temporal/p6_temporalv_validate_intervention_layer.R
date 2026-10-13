###
### Program to estimate calibration and discrimination for temporal validation.
### Risks are estimated in an interventional manner (i.e. how the model would be used in practice)
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

### Read in the entire unimputed cohort at baseline (this is both development and validation datasets)
### This will be used in subsequent functions
entire_cohort_baseline <- readRDS("data/p5/cohort_prototype3.rds")

### These are currently up in the air, so going to read in from Bowens table for now
direct_HR_smoking_initiation <- exp(readRDS("data/p5/offsets_direct_lnHR_smok_init.rds"))
direct_HR_smoking_cessation <- exp(readRDS("data/p5/offsets_direct_lnHR_smok_cess.rds"))
direct_HR_nonhdl <- exp(readRDS("data/p5/offsets_direct_lnHR_nonhdl.rds"))
direct_HR_bmi <- exp(readRDS("data/p5/offsets_direct_lnHR_bmi.rds"))
direct_HR_sbp <- exp(readRDS("data/p5/offsets_direct_lnHR_sbp.rds"))

###
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_validate_interventional_wrapper <- function(gender, 
                                                t_fup_int, # number of years defining the follow-up index date
                                                t_eval_int, # number of years prediction to evaluate performance at
                                                mfr,
                                                plot.range = NULL){
  
  print(paste("gender =", gender, "t_fup_int =", t_fup_int, "mfr" = mfr, Sys.time()))
  
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
  ### which is what the individual is doing, and name as temporal dataset
  df_valid_temporal <- 
    dplyr::mutate(data_valid,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  ### Reduce entire_cohort_baseline to those in the validation dataset through a merge
  df_valid_baseline <- dplyr::left_join(dplyr::select(df_valid_temporal, patid), 
                                        entire_cohort_baseline, 
                                        by = dplyr::join_by(patid))
  
  ### Now reduce to those who don't have missing data at baseline on the variable of interest
  df_valid_baseline <- df_valid_baseline[!is.na(df_valid_baseline[[mfr]]), ]
  
  ### Reduce temporal validation cohort to those who had non-missing data at baseline
  df_valid_temporal <- df_valid_temporal[!is.na(fastmatch::fmatch(df_valid_temporal$patid, df_valid_baseline$patid)), ]
  
  ### Create a new variable with altered name for the mfr at baseline
  df_valid_baseline[,"baseline_mfr"] <- df_valid_baseline[, mfr]
  
  ### Merge this into df_valid_temporal
  df_valid_temporal <- dplyr::left_join(df_valid_temporal, 
                                        df_valid_baseline[, c("patid", "baseline_mfr")], 
                                        by = dplyr::join_by(patid))
  
  ### Create a new variable with altered name for the mfr at baseline
  df_valid_temporal[,"followup_mfr"] <- df_valid_temporal[, mfr]
  
  ### For estimate the change, we are only interested in changes above some threshold, so create new variables which have this cut-off
  if (mfr == "sbp"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal, 
                                       baseline_mfr_cutoff = dplyr::case_when(baseline_mfr > 120 ~ baseline_mfr,
                                                                              TRUE ~ 120),
                                       followup_mfr_cutoff = dplyr::case_when(followup_mfr > 120 ~ followup_mfr,
                                                                              TRUE ~ 120))
  } else if (mfr == "nonhdl"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal, 
                                       baseline_mfr_cutoff = dplyr::case_when(baseline_mfr > 2.6 ~ baseline_mfr,
                                                                              TRUE ~ 2.6),
                                       followup_mfr_cutoff = dplyr::case_when(followup_mfr > 2.6 ~ followup_mfr,
                                                                              TRUE ~ 2.6))
  } else if (mfr == "bmi"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal, 
                                       baseline_mfr_cutoff = dplyr::case_when(baseline_mfr > 24 ~ baseline_mfr,
                                                                              TRUE ~ 24),
                                       followup_mfr_cutoff = dplyr::case_when(followup_mfr > 24 ~ followup_mfr,
                                                                              TRUE ~ 24))
  }
  
  ### Calculate the change in the modifiable risk factor between baseline and temporal validation cohorts
  ### This will be done differently depending on whether its a continuous mfr (sbp, bmi, nonhdl) or factor (smoking status)
  if (mfr != "smoking"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal,
                                       change_mfr = followup_mfr_cutoff - baseline_mfr_cutoff)
  } else if (mfr == "smoking"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal,
                                       change_smoking_initiation = dplyr::case_when(baseline_mfr == followup_mfr ~ 0,
                                                                                    baseline_mfr == "Ex-smoker" & followup_mfr == "Current" ~ 0,
                                                                                    baseline_mfr == "Current" & followup_mfr == "Ex-smoker" ~ 0,
                                                                                    baseline_mfr == "Non-smoker" & followup_mfr == "Current" ~ 1,
                                                                                    baseline_mfr == "Non-smoker" & followup_mfr == "Ex-smoker" ~ 1,
                                                                                    TRUE ~ 0),###NB: THIS SHOULD BE ~ NA, AND THERE SHOULD BE NO NA's, I AS GETTING ABOUT 800, FIXING BUG BUT LEAVING LIKE THIS FOR NOW
                                       change_smoking_cessation = dplyr::case_when(baseline_mfr == followup_mfr ~ 0,
                                                                                   baseline_mfr == "Ex-smoker" & followup_mfr == "Current" ~ -1,
                                                                                   baseline_mfr == "Current" & followup_mfr == "Ex-smoker" ~ 1,
                                                                                   baseline_mfr == "Non-smoker" & followup_mfr == "Current" ~ 0,
                                                                                   baseline_mfr == "Non-smoker" & followup_mfr == "Ex-smoker" ~ 1,
                                                                                   TRUE ~ 0))
  }
  
  ### Cap the reduction at 40 if variable is SBP
  if (mfr == "sbp"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal,
                                       change_mfr = dplyr::case_when(change_mfr > 40 ~ 40,
                                                                     change_mfr < -40 ~ -40,
                                                                     TRUE ~ change_mfr))
  }
  
  ###
  ### Estimate risks from initial risk estimation layer
  ###
  
  ### To do this, for all variables, we want to use the most recent value recorded, except the modifiable risk factor,
  ### where we will use the value from baseline.
  
  ### Set the modifiable risk factor to be the value from baseline
  df_valid_temporal[, mfr] <- df_valid_temporal[, "baseline_mfr"]
  
  ### Estimate risk from the initial risk estimation layer 
  df_valid_temporal$surv <- as.numeric(est_surv_offset_prototype3(newdata = df_valid_temporal, 
                                                                  fit = fit, 
                                                                  bhaz = bhaz, 
                                                                  time = t_eval))
  
  ###
  ### Estimate risks from intervention layer
  ###
  
  ### Estimate and apply the relative risk reduction
  if (mfr == "sbp"){
    direct_HR <- direct_HR_sbp
  } else if (mfr == "bmi"){
    direct_HR <- direct_HR_bmi
  } else if (mfr == "nonhdl"){
    direct_HR <- direct_HR_nonhdl
  } else if (mfr == "smoking"){
    direct_HR_initiation <- direct_HR_smoking_initiation
    direct_HR_cessation <- direct_HR_smoking_cessation
  }
  
  ### Apply the hazard ratio
  if (mfr != "smoking"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal,
                                       HR = direct_HR^change_mfr,
                                       apply_HR_scale = HR*(log(surv)),
                                       surv_intervention = exp(apply_HR_scale),
                                       pred_intervention = 1 - surv_intervention)
  } else if (mfr == "smoking"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal,
                                       HR = (direct_HR_initiation^change_smoking_initiation)*(direct_HR_cessation^change_smoking_cessation),
                                       apply_HR_scale = HR*(log(surv)),
                                       surv_intervention = exp(apply_HR_scale),
                                       pred_intervention = 1 - surv_intervention)
  }
  
  
  ##########################
  ### Assess calibration ###
  ##########################
  
  ### PH regression approach
  print(paste("PH", Sys.time()))
  calib_ph <- est_calib_ph(data = df_valid_temporal, 
                           fit = fit, 
                           bhaz = bhaz, 
                           time = t_eval,
                           surv = df_valid_temporal$surv_intervention)
  
  ### KM grouped approach
  print(paste("KM group", Sys.time()))
  calib_km_group <- est_calib_plot_group(data = df_valid_temporal, 
                                         fit = fit, 
                                         bhaz = bhaz,  
                                         time = t_eval,
                                         n.groups = 100,
                                         surv = df_valid_temporal$surv_intervention,
                                         single.CI = TRUE)
  
  ### Save
  saveRDS(calib_ph, paste("data/p5/temporalv_interventional_calib_ph_", gender, "_t", t_fup_int, "_", mfr, "_teval", t_eval_int,".rds", sep = ""))
  saveRDS(calib_km_group, paste("data/p5/temporalv_interventional_calib_km_group_", gender, "_t", t_fup_int, "_", mfr, "_teval", t_eval_int,".rds", sep = ""))
  
  #############################
  ### Assess discrimination ###
  #############################
  
  ### Estimate C-statistics
  Cstat_object <- survival::concordance(survival::Surv(cvd_time, cvd_indicator) ~ surv_intervention, 
                                        data = df_valid_temporal,
                                        timewt = "n",
                                        influence = 0)
  
  ### Save
  saveRDS(Cstat_object, paste0("data/p5/temporalv_interventional_discrim_", gender, "_t", t_fup_int, "_", mfr, "_teval", t_eval_int,".rds"))
  
  NULL
}

### Run function
lapply(c(1,2), function(x){
  lapply(1:5, function(y){
    lapply(c("bmi", "nonhdl", "sbp", "smoking"), function(z){
      est_validate_interventional_wrapper(gender = x, t_fup_int = y, t_eval_int = 5, mfr = z)
      est_validate_interventional_wrapper(gender = x, t_fup_int = y, t_eval_int = 10, mfr = z)
    })
  })
})

print(paste("FINISHED", Sys.time()))
