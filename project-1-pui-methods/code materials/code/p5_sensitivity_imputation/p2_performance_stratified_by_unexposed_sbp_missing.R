###
### Program to estimate calibration and discrimination stratified by whether unexposed SBP was imputed
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
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

### Define model to be unexposed mediator model
model <- 2

###
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_performance <- function(gender){
  
  ### Read in model
  fit <- readRDS(paste("data/fit_", gender, "_model", model, ".rds", sep = ""))
  bhaz <- readRDS(paste("data/bhaz_", gender, "_model", model, ".rds", sep = ""))
  
  ### Read in validation data
  df_valid <- readRDS(paste("data/df_imp_valid_", gender, sep = ""))
  
  ### Read in counterfactual survival times
  cf_surv_times <- readRDS(paste("data/cf_surv_times_", gender, "_model", model, ".rds", sep = ""))
  
  ### Replace cvd_time with counterfactual survival times
  df_valid <- dplyr::select(df_valid, -cvd_time)
  df_valid <- merge(df_valid, cf_surv_times, by.x = "patid", by.y = "patid") |>
    dplyr::rename(cvd_time = cvd_time_cf)
  #   print("data with counterfactual survival times")
  
  ### Create offset variable with appropriate names
  df_valid_pred <- 
    dplyr::mutate(df_valid,
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0)
  
  ### Read in original data, to identify those with missing sbp_unexposed
  df_unimputed <- readRDS(paste("data/cohort_", c("male", "female")[gender], "_pui.rds", sep = "")) |>
    dplyr::mutate(sbp_unexposed_miss_ind = as.numeric(is.na(sbp_unexposed))) |>
    dplyr::select(patid, sbp_unexposed_miss_ind)
  
  ### Merge with validation data
  df_valid_pred <- merge(df_valid_pred, df_unimputed,
                         by.x = "patid",
                         by.y = "patid",
                         all.x = TRUE)
  
  ### Assign time and status variables for compatibility with est_calib_ph (used in boot_func_for_calib_metrics_CI)
  df_valid_pred <- dplyr::mutate(df_valid_pred, time = cvd_time, status = cvd_indicator)
  
  ### Run loop to validation for missing and non-missing unexposed sbp
  for (miss_ind in c(0, 1)){
    
    print(paste("missing indicator = ", miss_ind, Sys.time()))
    
    ### Subset data
    df_valid_pred_subset <- dplyr::filter(df_valid_pred, sbp_unexposed_miss_ind == miss_ind)
    
    ##########################
    ### Assess calibration ###
    ##########################
    
    ##############################
    ### PH regression approach ###
    ##############################
    print(paste("PH", Sys.time()))
    calib_ph <- est_calib_plot(data = df_valid_pred_subset, 
                               fit = fit, 
                               bhaz = bhaz, 
                               time = round(10*365.25))
    
    #####################################
    ### Get CI's for ICI, E50 and E90 ###
    #####################################
    print(paste("PH boot", Sys.time()))
    calib_ph_boot <- boot::boot(data = df_valid_pred_subset, 
                                statistic = boot_func_for_calib_metrics_CI, 
                                R = 500, fit = fit, bhaz = bhaz, time = round(10*365.25), nk = 4)
    
    ###########################
    ### KM grouped approach ###
    ###########################
    print(paste("KM group", Sys.time()))
    calib_km_group <- est_calib_plot_group(data = df_valid_pred_subset, 
                                           fit = fit, 
                                           bhaz = bhaz,  
                                           time = round(10*365.25),
                                           n.groups = 50,
                                           CI = TRUE)
    
    ### Extract plot data and save
    df_calib_smooth <- calib_ph[["plotdata"]]
    df_calib_grouped <- calib_km_group[["plot"]]$data
    saveRDS(df_calib_smooth, paste("data/sens_imputation_stratified_calib_ph_df_smooth_", gender, "_model", model, "_missind", miss_ind, ".rds", sep = ""))
    saveRDS(df_calib_grouped, paste("data/sens_imputation_stratified_calib_ph_df_grouped_", gender, "_model", model, "_missind", miss_ind, ".rds", sep = ""))
    
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
    if (miss_ind == 0){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste(c("Male", "Female")[gender], " model: Recorded unexposed SBP", sep = "")) 
    } else if (miss_ind == 1){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste(c("Male", "Female")[gender], " model: Missing unexposed SBP", sep = "")) 
    } 
    
    ### Save plot
    ragg::agg_png(paste("figures/sens_imputation_stratified_calib_ph_", gender, "_model", model, 
                        "_missind", miss_ind, ".png", sep = ""), 
                  width = 1, height = 1, scaling = 1/5, unit = "in", res = 600)
    plot(ggplot_comb)
    dev.off()
    
    ### Print and save ICI, E50, E90 and the CI
    print(paste("ICI = ", calib_ph[["ICI"]]))
    print(paste("E50 = ", calib_ph[["E50"]]))
    print(paste("E90 = ", calib_ph[["E90"]]))
    saveRDS(calib_ph[["ICI"]], paste("data/sens_imputation_stratified_calib_ph_ICI_", gender, "_model", model, 
                                     "_missind", miss_ind, ".rds", sep = ""))
    saveRDS(calib_ph[["E50"]], paste("data/sens_imputation_stratified_calib_ph_E50_", gender, "_model", model, 
                                     "_missind", miss_ind, ".rds", sep = ""))
    saveRDS(calib_ph[["E90"]], paste("data/sens_imputation_stratified_calib_ph_E90_", gender, "_model", model, 
                                     "_missind", miss_ind, ".rds", sep = ""))
    saveRDS(calib_ph_boot$t, paste("data/sens_imputation_stratified_calib_ph_CI_ICI_E50_E90_", gender, "_model", model, 
                                   "_missind", miss_ind, ".rds", sep = ""))
    
    #############################
    ### Assess discrimination ###
    #############################
    
    ### Estimate risks
    ### Note that changing the time, which alter the risks, but not the order of the risks, and therefore
    ### will not effect the estimation of Harrels C
    surv <- as.numeric(est_surv_offset(newdata = df_valid_pred_subset, fit = fit, bhaz = bhaz, time = 10*365.25))
    
    ### Estimate C-statistics
    Cstat_object <- intsurv::cIndex(time = df_valid_pred_subset$cvd_time, 
                                    event = df_valid_pred_subset$cvd_indicator, 
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
    
    ### Save output
    saveRDS(Cstat_object, paste("data/sens_imputation_stratified_discrim_", gender, "_model", model, 
                                "_missind", miss_ind, ".rds", sep = ""))
    
  }
  
}

### Run this function
for (gender_in in c(2,1)){
  
  print(paste("gender = ", gender_in))
  
  est_performance(gender = gender_in)
  
}

print(paste("FINISHED", Sys.time()))

warnings()