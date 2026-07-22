###
### Program to estimate calibration in 1 development and validation dataset (split sample)
### The focus will be a temporal validation of the intervention layer
### i.e. we are evaluating performance of the model when the index date is defined as 1/2/3/4/5 years after
### the start of follow-up, but we apply the relative risk increase/reduction for changes in one of the modifiable risk factors
### from baseline (i.e. how the model will be used in practice at follow-up visits).
### This program will assess each modifiable risk factor on its own.
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
setwd("")
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

### Read in the appropriate effect estimtae (relative risk)
# direct_RR_smoking_initiation <- readRDS("data/p4/direct_RR_smoking_initiation")
# direct_RR_smoking_cessation <- readRDS("data/p4/direct_RR_smoking_cessation")
# direct_RR_nonhdl <- readRDS("data/p4/direct_RR_nonhdl")
# direct_RR_bmi <- readRDS("data/p4/direct_RR_bmi")
# direct_RR_sbp <- readRDS("data/p4/direct_RR_sbp")

### These are currently up in the air, so going to read in from Bowens table for now
direct_OR_smoking_initiation <- 1.49995
direct_OR_smoking_cessation <- 0.8511231
direct_OR_nonhdl <- 1/0.8120427
direct_OR_bmi <- 1/0.9768588
direct_OR_sbp <- 1/0.9741304

### Read in the entire unimputed cohort at baseline (this is both development and validation datasets)
### This will be used in subsequent functions
entire_cohort_baseline <- readRDS("data/p4/cohort_prototype3.rds")

###
### Write a function to prepare data and generate risks, before assessing calibration or discrimination
### Note this function has no inputs, but in reliant on inputs that must be specified. This means it will only 
### work when being run inside the subsequent functions for estimating calibration (est_calib_ph_d1v1) and discrimination (est_discrim_d1v1)
###
prepare_data <- function(){
  
  ### Read in temporal validation dataset
  df_valid_temporal <- readRDS(paste("data/p4/df_valid_temporalv_", gender, "_", t_fup, ".rds", sep = ""))
  
  ### Reduce entire_cohort_baseline to those in the validation dataset through a merge
  df_valid_baseline <- dplyr::left_join(dplyr::select(df_valid_temporal, patid), 
                                        entire_cohort_baseline, 
                                        by = dplyr::join_by(patid))
  
  ### Now reduce to those who don't have missing data at baseline on the variable of interest
  df_valid_baseline <- df_valid_baseline[!is.na(df_valid_baseline[[mfr]]), ]
  
  ### Reduce temporal validation cohort to those who had non-missing data at baseline
  df_valid_temporal <- df_valid_temporal[!is.na(fastmatch::fmatch(df_valid_temporal$patid, df_valid_baseline$patid)), ]
  
  ### Read in model and counter factual survival times
  if (age_knots == 3){
    ### Read in model
    fit <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
    ### Read in counterfactual survival times
    cf_surv_times <- readRDS(paste("data/p4/prototype3_cf_surv_times_", gender, 
                                   "_devel", 1, 
                                   "_valid", 1, 
                                   "_t", t_fup, ".rds", sep = ""))
  } else if (age_knots == 4){
    
    if (caltime_in == TRUE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf_surv_times <- readRDS(paste("data/p4/prototype3_4knots_caltime_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     "_t", t_fup, ".rds", sep = ""))
    }
    
    if (caltime_in == FALSE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf_surv_times <- readRDS(paste("data/p4/prototype3_4knots_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     "_t", t_fup, ".rds", sep = ""))
    }
  }
  
  ### Replace cvd_time with counterfactual survival times
  df_valid_temporal <- dplyr::select(df_valid_temporal, -cvd_time)
  df_valid_temporal <- merge(df_valid_temporal, cf_surv_times, by.x = "patid", by.y = "patid") |>
    dplyr::rename(cvd_time = cvd_time_cf)
  #   print("data with counterfactual survival times")
  
  ### Create variable with forms that can be used for prediction, and agrees with 'prescribed treatment strategy'
  ### which is what the individual is doing
  df_valid_temporal <- 
    dplyr::mutate(df_valid_temporal,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  #############################################################################
  ### Estimate predicted risks under changes in the modifiable risk factors ###
  #############################################################################
  
  ###
  ### Preliminary calculations and data manipulation
  ###
  
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
  df_valid_temporal$pred <- 1 - as.numeric(est_surv_offset_prototype3(newdata = df_valid_temporal, 
                                                                      fit = fit, 
                                                                      bhaz = bhaz, 
                                                                      time = round(t_eval*365.25)))
  
  ###
  ### Estimate risks from intervention layer
  ###
  
  ### Estimate and apply the relative risk reduction
  if (mfr == "sbp"){
    direct_OR <- direct_OR_sbp
  } else if (mfr == "bmi"){
    direct_OR <- direct_OR_bmi
  } else if (mfr == "nonhdl"){
    direct_OR <- direct_OR_nonhdl
  } else if (mfr == "smoking"){
    direct_OR_initiation <- direct_OR_smoking_initiation
    direct_OR_cessation <- direct_OR_smoking_cessation
  }
  
  ### Apply the odds ratio
  if (mfr != "smoking"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal,
                                       OR = direct_OR^change_mfr,
                                       pred_intervention = pred*OR/(1 - pred + pred*OR),
                                       surv_intervention = 1 - pred_intervention)
  } else if (mfr == "smoking"){
    df_valid_temporal <- dplyr::mutate(df_valid_temporal,
                                       OR = (direct_OR_initiation^change_smoking_initiation)*(direct_OR_cessation^change_smoking_cessation),
                                       pred_intervention = pred*OR/(1 - pred + pred*OR),
                                       surv_intervention = 1 - pred_intervention)
  }
  
  return(df_valid_temporal)
  
}
  
###
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_calib_ph_d1v1 <- function(mfr,
                              gender,
                              t_fup, # number of days defining the follow-up index date
                              t_eval, # number of years prediction to evaluate performance at
                              age_knots,
                              caltime_in,
                              plot.range = NULL){
  
  # mfr <- "smoking"
  # gender <- 2
  # t_fup <- round(365.25*as.numeric(2))
  # t_eval <- 5
  # age_knots <- 4
  # caltime_in <- TRUE
  # plot.range = NULL
  
  ### Prepare data
  df_valid_temporal <- prepare_data()
 
  ##########################
  ### Assess calibration ###
  ##########################
  
  ##############################
  ### PH regression approach ###
  ##############################
  print(paste("PH", Sys.time()))
  calib.ph <- est_calib_plot(data = df_valid_temporal, 
                             fit = fit, 
                             bhaz = bhaz, 
                             time = round(t_eval*365.25),
                             fit.type = "prototype3",
                             surv = df_valid_temporal$surv_intervention,
                             type = "single")
  
  ###########################
  ### KM grouped approach ###
  ###########################
  print(paste("KM group", Sys.time()))
  calib.km.group <- est_calib_plot_group(data = df_valid_temporal, 
                                         fit = fit, 
                                         bhaz = bhaz,  
                                         time = round(t_eval*365.25),
                                         n.groups = 100,
                                         fit.type = "prototype3",
                                         surv = df_valid_temporal$surv_intervention,
                                         type = "single",
                                         single.CI = TRUE)
  
  ### Extract plot data
  df.calib.smooth <- calib.ph[["plotdata"]]
  df.calib.grouped <- calib.km.group[["plot"]]$data
  
  ### Assign label for mfr
  if (mfr == "sbp"){mfr_label <- "systolic blood pressure"}
  if (mfr == "nonhdl"){mfr_label <- "non-HDL cholesterol"}
  if (mfr == "smoking"){mfr_label <- "smoking status"}
  if (mfr == "bmi"){mfr_label <- "body mass index"}
  ### Create a ggplot out of this
  ggplot.comb <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), 
                       data = df.calib.smooth) + 
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
    ggplot2::geom_point(data = df.calib.smooth,  
                        ggplot2::aes(x = pred, y = pred.obs), col = grDevices::rgb(0, 0, 0, alpha = 0)) +
    ggplot2::geom_point(data = df.calib.grouped, 
                        ggplot2::aes(x = pred, y = obs, col = grDevices::rgb(0, 1, 0, alpha = 1))) + 
    ggplot2::xlim(c(0,0.3)) + ggplot2::ylim(c(0, 0.3)) +
    ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(paste("Modifiable risk factor = ", mfr_label, ", \nfollow-up visit ", round(t_fup/365.25), ", t = ", t_eval, " years", sep = ""))
  
  ### Create output.object
  output.object <- list("plot" = ggplot.comb,
                        "plotdata.smooth" = calib.ph[["plotdata"]],
                        "plotdata.range" = calib.ph[["plotdata.range"]],
                        "plotdata.grouped" = calib.km.group[["plot"]]$data,
                        "ICI" = as.numeric(calib.ph[["ICI"]]),
                        "E50" = as.numeric(calib.ph[["E50"]]),
                        "E90" =  as.numeric(calib.ph[["E90"]]))
  
  ### Return output object
  return(output.object)
  
}


###
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_discrim_d1v1 <- function(mfr,
                             gender,
                             t_fup, # number of days defining the follow-up index date
                             t_eval, # number of years prediction to evaluate performance at
                             age_knots,
                             caltime_in){
  
  # mfr <- "sbp"
  # gender <- 2
  # t_fup <- round(365.25*as.numeric(2))
  # t_eval <- 5
  # age_knots <- 4
  
  ### Prepare data
  df_valid_temporal <- prepare_data()
  
  #############################
  ### Assess discrimination ###
  #############################
  
  ### Estimate C-statistics
  Cstat_object <- intsurv::cIndex(time = df_valid_temporal$cvd_time, 
                                  event = df_valid_temporal$cvd_indicator, 
                                  risk_score = 1 - df_valid_temporal$surv_intervention)
  
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
  
  ### Return output object
  return(Cstat_object)
  
}


###
### Get calibration plots and discrimination
###

### Define gender
for (gender in c(1,2)){
  
  print(paste("gender = ", c("male", "female")[gender]))
  
  ### Define age knots
  for (age_knots in c(4)){
    
    ### Define whether calendar time included in model
    for (caltime_in in c(TRUE, FALSE)){
      
      print(paste("caltime = ", caltime_in))
      
      ### Define t_fup
      for (t in 1:5){
        
        t_fup <- round(t*365.25)
        print(paste("t_fup = ", t_fup))
        
        ### Define t_eval
        for (t_eval in c(5,10)){
          
          print(paste("t_eval = ", t_eval))
          
          ### Define modifiable risk factor
          for (mfr in c("smoking", "sbp", "bmi", "nonhdl")){
            
            print(paste("mfr = ", mfr))
            
            ###################
            ### Calibration ###
            ###################
            
            ### Estimate calibration curves for development=1 and validation=1 datasets
            calib_MxM <- est_calib_ph_d1v1(mfr = mfr, gender = gender, t_fup = t_fup, t_eval = t_eval, age_knots = age_knots, caltime_in = caltime_in)
            
            ### Save output
            # png(paste("figures/p4/temporalv_intervention_prototype3_calib_ph_d1v1_", gender, 
            #           "_nk", age_knots, 
            #           "_caltime", as.numeric(caltime),
            #           "_tfup", t_fup, 
            #           "_teval", t_eval, 
            #           "_", mfr, ".png", sep = ""), width = 6, height = 6, unit = "in", res = 600)
            # plot(calib_MxM[["plot"]])
            # dev.off()
            ragg::agg_png(paste("figures/p4/temporalv_intervention_prototype3_calib_ph_d1v1_", gender, 
                                "_nk", age_knots, 
                                "_caltime", as.numeric(caltime_in),
                                "_tfup", t_fup, 
                                "_teval", t_eval, 
                                "_", mfr, ".png", sep = ""), width = 1, height = 1, scaling = 1/6, unit = "in", res = 600)
            plot(calib_MxM[["plot"]])
            dev.off()
            
            saveRDS(calib_MxM[["plotdata.smooth"]], paste("data/p4/temporalv_intervention_prototype3_calib_ph_d1v1_plotdata_smooth_", gender, 
                                                          "_nk", age_knots, 
                                                          "_caltime", as.numeric(caltime_in),
                                                          "_tfup", t_fup, 
                                                          "_teval", t_eval, 
                                                          "_", mfr, ".rds", sep = ""))
            saveRDS(calib_MxM[["plotdata.grouped"]], paste("data/p4/temporalv_intervention_prototype3_calib_ph_d1v1_plotdata_grouped_", gender, 
                                                           "_nk", age_knots, 
                                                           "_caltime", as.numeric(caltime_in),
                                                           "_tfup", t_fup, 
                                                           "_teval", t_eval, 
                                                           "_", mfr, ".rds", sep = ""))
            
            ### print ICI, E50, E90
            print(paste("ICI = ", calib_MxM[["ICI"]]))
            print(paste("E50 = ", calib_MxM[["E50"]]))
            print(paste("E90 = ", calib_MxM[["E90"]]))
            saveRDS(calib_MxM[["ICI"]], paste("data/p4/temporalv_intervention_prototype3_calib_ph_ICI_d1v1_", gender, 
                                              "_nk", age_knots, 
                                              "_caltime", as.numeric(caltime_in),
                                              "_tfup", t_fup, 
                                              "_teval", t_eval, 
                                              "_", mfr, ".rds", sep = ""))
            saveRDS(calib_MxM[["E50"]], paste("data/p4/temporalv_intervention_prototype3_calib_ph_E50_d1v1_", gender, 
                                              "_nk", age_knots, 
                                              "_caltime", as.numeric(caltime_in),
                                              "_tfup", t_fup, 
                                              "_teval", t_eval, 
                                              "_", mfr, ".rds", sep = ""))
            saveRDS(calib_MxM[["E90"]], paste("data/p4/temporalv_intervention_prototype3_calib_ph_E90_d1v1_", gender, 
                                              "_nk", age_knots, 
                                              "_caltime", as.numeric(caltime_in),
                                              "_tfup", t_fup, 
                                              "_teval", t_eval, 
                                              "_", mfr, ".rds", sep = ""))
            
            ######################
            ### Discrimination ###
            ######################
            
            ### Estimate discrimination curves for development=1 and validation=1 datasets
            discrim_d1v1 <- est_discrim_d1v1(mfr = mfr, gender = gender, t_fup = t_fup, t_eval = t_eval, age_knots = age_knots, caltime_in = caltime_in)
            
            ### Save output
            saveRDS(discrim_d1v1, paste("data/p4/temporalv_intervention_prototype3_discrim_d1v1_", gender,
                                        "_nk", age_knots,
                                        "_caltime", as.numeric(caltime_in),
                                        "_tfup", t_fup,
                                        "_teval", t_eval,
                                        "_", mfr, ".rds", sep = ""))
            print(discrim_d1v1)
            
          }
        }
      }
    }
  }
}

print(paste("FINISHED", Sys.time()))


