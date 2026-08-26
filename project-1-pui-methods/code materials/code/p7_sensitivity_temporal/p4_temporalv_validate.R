###
### Program to estimate calibration and discrimination.
### Temporal validation of the initial risk estimation layer.
### We are evaluating performance of the model when the index date is defined as 1/2/3/4/5 years after
### the start of follow-up.
###
### Predicted risks are under the treatment strategy of 'continue on current treatment strategy'. Counterfactual survival
### times have been calculated under this treatment strategy.
### 

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- 180

### Define model = 4 (two component model) for these analyses
model <- 4

### Define number of years for evaluating performance
t_eval <- 5

### Load packages
library(survival)

### Extract gender from command line
args <- commandArgs(trailingOnly = T)
gender <- as.numeric(args[1])
gender_char <- c("male", "female")[gender]
print(paste("gender = ", gender_char))

### Extract followup time from command line
t_fup_integer <-as.numeric(args[2])
print(paste("t_fup_integer = ", t_fup_integer))

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

###
### Write a function that will estimate calibration curve for a given t_fup
###
est_performance <- function(gender, t_fup_integer, assess_calibration = TRUE, assess_discrim = TRUE){
  
  ### Define t_fup (want t_fup_integer for ggtitle)
  t_fup <- round(365.25*t_fup_integer)
  
  ### Read in model
  fit <- readRDS(paste("data/fit_", gender, "_model", model, ".rds", sep = ""))
  bhaz <- readRDS(paste("data/bhaz_", gender, "_model", model, ".rds", sep = ""))
  
  ### Read in validation data
  df_valid <- readRDS(paste("data/sens_temporalv_df_valid_", gender, "_", t_fup, ".rds", sep = ""))
  
    ### Read in counterfactual survival times
    cf_surv_times <- readRDS(paste("data/sens_temporalv__cf_surv_times_", gender, 
                                   "_t", t_fup, 
                                   ".rds", sep = ""))
 
    ### Replace cvd_time with counterfactual survival times
    df_valid <- dplyr::select(df_valid, -cvd_time)
    df_valid <- merge(df_valid, cf_surv_times, by.x = "patid", by.y = "patid") |>
      dplyr::rename(cvd_time = cvd_time_cf)
 
  ### Create offset variable with appropriate names
  df_valid_pred <- 
    dplyr::mutate(df_valid,
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0)
  
  ### Assign time and status variables for compatibility with est_calib_ph (used in boot_func_for_calib_metrics_CI)
  df_valid_pred <- dplyr::mutate(df_valid_pred, time = cvd_time, status = cvd_indicator)
  
  ##########################
  ### Assess calibration ###
  ##########################
  
  if (assess_calibration == TRUE){
    
    ##############################
    ### PH regression approach ###
    ##############################
    print(paste("PH", Sys.time()))
    calib_ph <- est_calib_plot(data = df_valid_pred, 
                               fit = fit, 
                               bhaz = bhaz, 
                               time = round(t_eval*365.25))
    
    #####################################
    ### Get CI's for ICI, E50 and E90 ###
    #####################################
    print(paste("PH boot", Sys.time()))
    calib_ph_boot <- boot::boot(data = df_valid_pred, 
                                statistic = boot_func_for_calib_metrics_CI, 
                                R = 500, fit = fit, bhaz = bhaz, time = round(t_eval*365.25), nk = 4)
    
    ###########################
    ### KM grouped approach ###
    ###########################
    print(paste("KM group", Sys.time()))
    calib_km_group <- est_calib_plot_group(data = df_valid_pred, 
                                           fit = fit, 
                                           bhaz = bhaz,  
                                           time = round(t_eval*365.25),
                                           n.groups = 25,
                                           CI = TRUE)
    
    ### Extract plot data and save
    df_calib_smooth <- calib_ph[["plotdata"]]
    df_calib_grouped <- calib_km_group[["plot"]]$data
    saveRDS(df_calib_smooth, paste("data/sens_temporalv_calib_ph_df_smooth_", gender, "_tfup", t_fup, ".rds", sep = ""))
    saveRDS(df_calib_grouped, paste("data/sens_temporalv_calib_ph_df_grouped_", gender, "_tfup", t_fup, ".rds", sep = ""))
    
    ### Create a ggplot out of this
    ggplot_comb <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), 
                         data = df_calib_smooth) + 
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
      ggplot2::geom_point(data = df_calib_smooth,  
                          ggplot2::aes(x = pred, y = pred.obs), col = grDevices::rgb(0, 0, 0, alpha = 0)) +
      ggplot2::geom_point(data = df_calib_grouped, 
                          ggplot2::aes(x = pred, y = obs, col = grDevices::rgb(0, 1, 0, alpha = 1))) + 
      ggplot2::theme(legend.position = "none") + 
      ggplot2::xlim(c(0,0.6)) + ggplot2::ylim(c(0,0.6))
    
    ### Add title
    ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(
      paste("Temporal validation; Calibration of two-component model at ", 
            t_fup_integer, " year", ifelse(t_fup_integer == 1, "", "s"), " post index date", sep = ""), width = 40)) 
    
    ### Save plot
    ragg::agg_png(paste("figures/sens_temporalv_calib_ph_", gender, "_tfup", t_fup, ".png", sep = ""), 
                  width = 1, height = 1, scaling = 1/5, unit = "in", res = 600)
    plot(ggplot_comb)
    dev.off()
    
    ### Print and save ICI, E50, E90 and the CI
    print(paste("ICI = ", calib_ph[["ICI"]]))
    print(paste("E50 = ", calib_ph[["E50"]]))
    print(paste("E90 = ", calib_ph[["E90"]]))
    saveRDS(calib_ph[["ICI"]], paste("data/sens_temporalv_calib_ph_ICI_", gender, "_tfup", t_fup, ".rds", sep = ""))
    saveRDS(calib_ph[["E50"]], paste("data/sens_temporalv_calib_ph_E50_", gender, "_tfup", t_fup, ".rds", sep = ""))
    saveRDS(calib_ph[["E90"]], paste("data/sens_temporalv_calib_ph_E90_", gender, "_tfup", t_fup, ".rds", sep = ""))
    saveRDS(calib_ph_boot$t, paste("data/sens_temporalv_calib_ph_CI_ICI_E50_E90_", gender, "_tfup", t_fup, ".rds", sep = ""))
    
    print(paste("CALIB FINISHED", Sys.time()))
    
  }
  
  #############################
  ### Assess discrimination ###
  #############################
  
  if (assess_discrim == TRUE){
    
    ### Estimate risks
    ### Note that changing the time, which alter the risks, but not hte order of the risks, and therefore
    ### will not effect the estimation of Harrels C
    surv <- as.numeric(est_surv_offset(newdata = df_valid_pred, fit = fit, bhaz = bhaz, time = t_eval*365.25))
    
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
    print(Cstat_object)
    
    ### Return output object
    saveRDS(Cstat_object, paste("data/sens_temporalv_discrim_", gender, "_tfup", t_fup, ".rds", sep = ""))
    print(paste("DISCRIM FINISHED", Sys.time()))
    
  }
  
}

### Run this function
est_performance(gender = gender, t_fup_integer = t_fup_integer)
print(paste("FINISHED", Sys.time()))