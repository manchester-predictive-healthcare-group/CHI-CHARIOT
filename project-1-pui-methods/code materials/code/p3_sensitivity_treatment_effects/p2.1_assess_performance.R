### Assess performance for modified treatment effect in validation cohort

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load survival package
library(survival)

### Extract arguments from command line
args <- commandArgs(trailingOnly = T)
gender_in <- as.numeric(args[1])
scenario <- as.numeric(args[2])
print(paste("gender = ", gender_in))
print(paste("scenario = ", scenario))

### Define model
model <- 4

### Create table for amounts the hazard-ratios will be adjusted by in the validation dataset
scenario_adjustments <- expand.grid(statins = c(0.8, 0.9, 1, 1.1, 1.2),
                                    ah = c(0.8, 0.9, 1, 1.1, 1.2)) |>
  dplyr::filter(!(statins == 1 & ah == 1))

### Get row numbers where either effect is unchanged
rows_1_unchanged <- which(scenario_adjustments$statins == 1 | scenario_adjustments$ah == 1)

### Reduce to this
scenario_adjustments <- scenario_adjustments[rows_1_unchanged, ]

###
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_performance <- function(gender, adjustment_statins, adjustment_ah){
  
  ### Read in model
  fit <- readRDS(paste("data/fit_", gender, "_model", model, ".rds", sep = ""))
  bhaz <- readRDS(paste("data/bhaz_", gender, "_model", model, ".rds", sep = ""))
  
  ### Read in validation data
  df_valid <- readRDS(paste("data/df_imp_valid_", gender, sep = ""))
  
  ### Read in counterfactual survival times
  cf_surv_times <- readRDS(paste("data/cf_surv_times_treatment_effect_modified_", gender, "_model", model, 
                                 "_statins", adjustment_statins, "_ah", adjustment_ah, ".rds", sep = ""))
  
  ### Replace cvd_time with counterfactual survival times
  df_valid <- dplyr::select(df_valid, -cvd_time)
  df_valid <- merge(df_valid, cf_surv_times, by.x = "patid", by.y = "patid") |>
    dplyr::rename(cvd_time = cvd_time_cf)
  #   print("data with counterfactual survival times")
  
  ### Assign time and status variables for compatibility with est_calib_ph (used in boot_func_for_calib_metrics_CI)
  df_valid <- dplyr::mutate(df_valid, time = cvd_time, status = cvd_indicator)
  
  ### Create offset variable with appropriate names
  df_valid_pred <- 
    dplyr::mutate(df_valid,
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0)

  ##########################
  ### Assess calibration ###
  ##########################
  
  ##############################
  ### PH regression approach ###
  ##############################
  print(paste("PH", Sys.time()))
  calib_ph <- est_calib_ph(data = df_valid_pred, 
                             fit = fit, 
                             bhaz = bhaz, 
                             time = round(10*365.25))
  
  ### Define pred.plot.range for CI's, cut-off top 1%
  pred_plot_range <- seq(0.001, quantile(as.numeric(calib_ph[["calib_data"]]$pred), p = 0.99), length.out = 250)
  
  #####################################
  ### Get CI's for ICI, E50 and E90 ###
  #####################################
  print(paste("PH boot", Sys.time()))
  calib_ph_boot <- boot::boot(data = df_valid_pred, 
                              statistic = boot_func_for_calib_metrics_CI, 
                              R = 1000, fit = fit, bhaz = bhaz, time = round(10*365.25), nk = 4, 
                              pred.plot.range = pred_plot_range)
  
  ###########################
  ### KM grouped approach ###
  ###########################
  print(paste("KM group", Sys.time()))
  calib_km_group <- est_calib_plot_group(data = df_valid_pred, 
                                         fit = fit, 
                                         bhaz = bhaz,  
                                         time = round(10*365.25),
                                         n.groups = 50,
                                         CI = TRUE)
  
  ### Extract plot data and save
  df_calib_smooth <- calib_ph[["calib_data"]]
  df_calib_grouped <- calib_km_group[["plot"]]$data
  saveRDS(df_calib_smooth, paste("data/sens_treatment_effect_modified_calib_ph_df_smooth_", gender,  
                                 "_statins", adjustment_statins, "_ah", adjustment_ah, ".rds", sep = ""))
  saveRDS(df_calib_grouped, paste("data/sens_treatment_effect_modified_calib_ph_df_grouped_", gender,  
                                  "_statins", adjustment_statins, "_ah", adjustment_ah, ".rds", sep = ""))
  
  ### Create a ggplot out of this
  ggplot_comb <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), 
                       data = df_calib_smooth) + 
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
    ggplot2::geom_point(data = df_calib_smooth,  
                        ggplot2::aes(x = pred, y = pred.obs), 
                        colour = "black",
                        alpha = 0) +
    ggplot2::geom_point(data = df_calib_grouped, 
                        ggplot2::aes(x = pred, y = obs),
                        colour = "green3") + 
    ggplot2::theme(legend.position = "none") + 
    ggplot2::xlim(c(0,0.6)) + ggplot2::ylim(c(0,0.6))
  
  ### Add title
  ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Statins effect multiplier: ", adjustment_statins, "; ",
                                                      "Antihypertensives effect multiplier: ", adjustment_ah, 
                                                      sep = "")) 
  
  ### Save plot
  ragg::agg_png(paste("figures/sens_treatment_effect_modified_calib_ph_", gender, 
                      "_statins", adjustment_statins, "_ah", adjustment_ah, ".png", sep = ""), 
                width = 1, height = 1, scaling = 1/5, unit = "in", res = 600)
  plot(ggplot_comb)
  dev.off()
  
  ### Print and save ICI, E50, E90 and the CI
  print(paste("ICI = ", calib_ph[["ICI"]]))
  print(paste("E50 = ", calib_ph[["E50"]]))
  print(paste("E90 = ", calib_ph[["E90"]]))
  saveRDS(calib_ph[["ICI"]], paste("data/sens_treatment_effect_modified_calib_ph_ICI_", gender,  
                                   "_statins", adjustment_statins, "_ah", adjustment_ah, ".rds", sep = ""))
  saveRDS(calib_ph[["E50"]], paste("data/sens_treatment_effect_modified_calib_ph_E50_", gender,  
                                   "_statins", adjustment_statins, "_ah", adjustment_ah, ".rds", sep = ""))
  saveRDS(calib_ph[["E90"]], paste("data/sens_treatment_effect_modified_calib_ph_E90_", gender, 
                                   "_statins", adjustment_statins, "_ah", adjustment_ah, ".rds", sep = ""))
  saveRDS(calib_ph_boot$t, paste("data/sens_treatment_effect_modified_calib_ph_CI_ICI_E50_E90_", gender,  
                                 "_statins", adjustment_statins, "_ah", adjustment_ah, ".rds", sep = ""))
  saveRDS(pred_plot_range, paste("data/sens_treatment_effect_modified_calib_ph_pred_plot_range_", gender,  
                                 "_statins", adjustment_statins, "_ah", adjustment_ah, ".rds", sep = ""))
  
  #############################
  ### Assess discrimination ###
  #############################
  
  ### Estimate risks
  ### Note that changing the time, which alter the risks, but not the order of the risks, and therefore
  ### will not effect the estimation of Harrels C
  surv <- as.numeric(est_surv_offset(newdata = df_valid_pred, fit = fit, bhaz = bhaz, time = 10*365.25))
  
  ### Estimate C-statistics
  Cstat_object <- intsurv::cIndex(time = df_valid_pred$time, 
                                  event = df_valid_pred$status, 
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
  saveRDS(Cstat_object, paste("data/sens_treatment_effect_modified_discrim_", gender, 
                              "_statins", adjustment_statins, "_ah", adjustment_ah, ".rds", sep = ""))
  
}

### Estimate calibration
est_performance(gender = gender_in, 
                adjustment_statins = scenario_adjustments[scenario,1], 
                adjustment_ah = scenario_adjustments[scenario,2])

print(paste("FINISHED", Sys.time()))
    