###
### Program to estimate calibration in 1 development and validation dataset (split sample)
### The focus will be a temporal validation of the initial risk estimation layer
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
est_calib_ph_d1v1 <- function(gender, 
                              t_fup, # number of days defining the follow-up index date
                              t_eval, # number of years prediction to evaluate performance at
                              age_knots,
                              caltime_in,
                              plot.range = NULL){
  
  # gender <- 2
  # t_fup <- round(365.25*as.numeric(2))
  # age_knots = 4
  # plot.range = NULL
  # str(data.valid)
  # 
  ### Read in model and counter factual survival times
  if (age_knots == 3){
    ### Read in model
    fit <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
    ### Read in counterfactual survival times
    cf.surv.times <- readRDS(paste("data/p4/prototype3_cf_surv_times_", gender, 
                                   "_devel", 1, 
                                   "_valid", 1, 
                                   "_t", t_fup, ".rds", sep = ""))
  } else if (age_knots == 4){
    
    if (caltime_in == TRUE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_caltime_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     "_t", t_fup, ".rds", sep = ""))
    }
    
    if (caltime_in == FALSE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     "_t", t_fup, ".rds", sep = ""))
    }
    
  }
  
  ### Read in validation dataset
  data.valid <- readRDS(paste("data/p4/df_valid_temporalv_", gender, "_", t_fup, ".rds", sep = ""))
  colnames(data.valid)
  
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
  
  ##########################
  ### Assess calibration ###
  ##########################
  
  ##############################
  ### PH regression approach ###
  ##############################
  print(paste("PH", Sys.time()))
  calib.ph <- est_calib_plot(data = data.valid.pred, 
                             fit = fit, 
                             bhaz = bhaz, 
                             time = round(t_eval*365.25),
                             fit.type = "prototype3",
                             type = "single")
  
  ###########################
  ### KM grouped approach ###
  ###########################
  print(paste("KM group", Sys.time()))
  calib.km.group <- est_calib_plot_group(data = data.valid.pred, 
                                         fit = fit, 
                                         bhaz = bhaz,  
                                         time = round(t_eval*365.25),
                                         n.groups = 100,
                                         fit.type = "prototype3",
                                         type = "single",
                                         single.CI = TRUE
  )
  
  ### Extract plot data
  df.calib.smooth <- calib.ph[["plotdata"]]
  df.calib.grouped <- calib.km.group[["plot"]]$data
  
  ### Create a ggplot out of this
  ggplot.comb <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), 
                       data = df.calib.smooth) + 
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
    ggplot2::geom_point(data = df.calib.smooth,  
                        ggplot2::aes(x = pred, y = pred.obs), col = grDevices::rgb(0, 0, 0, alpha = 0)) +
    ggplot2::geom_point(data = df.calib.grouped, 
                        ggplot2::aes(x = pred, y = obs, col = grDevices::rgb(0, 1, 0, alpha = 1)))  + 
    ggplot2::xlim(c(0,0.3)) + ggplot2::ylim(c(0, 0.3)) +
    ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(paste("Follow-up visit ", round(t_fup/365.25), ", t = ", t_eval, " years", sep = ""))
  
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
est_discrim_d1v1 <- function(gender, 
                             t_fup, # number of days defining the follow-up index date
                             t_eval, # number of years prediction to evaluate performance at
                             age_knots,
                             caltime_in){
  
  # gender <- 2
  # t_fup <- round(365.25*as.numeric(2))
  # age_knots = 4
  # plot.range = NULL
  # str(data.valid)
  # 
  ### Read in model and counter factual survival times
  if (age_knots == 3){
    ### Read in model
    fit <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
    ### Read in counterfactual survival times
    cf.surv.times <- readRDS(paste("data/p4/prototype3_cf_surv_times_", gender, 
                                   "_devel", 1, 
                                   "_valid", 1, 
                                   "_t", t_fup, ".rds", sep = ""))
  } else if (age_knots == 4){
    
    if (caltime_in == TRUE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_caltime_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     "_t", t_fup, ".rds", sep = ""))
    }
    
    if (caltime_in == FALSE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     "_t", t_fup, ".rds", sep = ""))
    }
  }
  
  ### Read in validation dataset
  data.valid <- readRDS(paste("data/p4/df_valid_temporalv_", gender, "_", t_fup, ".rds", sep = ""))
  colnames(data.valid)
  
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
  surv <- as.numeric(est_surv_offset_prototype3(newdata = data.valid.pred, fit = fit, bhaz = bhaz, time = t_eval*365.25))
  
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
### Get calibration plots and discrimination
###

### Define gender
for (gender in c(1,2)){
  
  print(paste("gender = ", c("male", "female")[gender]))
  
  ### Define age knots
  for (age_knots in c(4)){
    
    ### Define whether calendar time included in model
    for (caltime_in in c(TRUE, FALSE)){
      
      print(paste("caltime_in = ", caltime_in))
      
      ### Define t_fup
      for (t in 1:5){
        
        t_fup <- round(t*365.25)
        print(paste("t_fup = ", t_fup))
        
        ### Define t_eval
        for (t_eval in c(5,10)){
          
          print(paste("t_eval = ", t_eval))
          
          ###################
          ### Calibration ###
          ###################
          
          ### Estimate calibration curves for development=1 and validation=1 datasets
          calib_MxM <- est_calib_ph_d1v1(gender = gender, t_fup = t_fup, t_eval = t_eval, age_knots = age_knots, caltime_in = caltime_in)
          
          ### Save plot
          # png(paste("figures/p4/temporalv_initial_prototype3_calib_ph_d1v1_", gender, 
          #           "_nk", age_knots, 
          #           "_caltime", as.numeric(caltime),
          #           "_tfup", t_fup, 
          #           "_teval", t_eval, ".png", sep = ""), width = 6, height = 6, unit = "in", res = 600)
          # plot(calib_MxM[["plot"]])
          # dev.off()
          ragg::agg_png(paste("figures/p4/temporalv_initial_prototype3_calib_ph_d1v1_", gender, 
                              "_nk", age_knots, 
                              "_caltime", as.numeric(caltime_in),
                              "_tfup", t_fup, 
                              "_teval", t_eval, ".png", sep = ""), width = 1, height = 1, scaling = 1/6, unit = "in", res = 600)
          plot(calib_MxM[["plot"]])
          dev.off()
          
          saveRDS(calib_MxM[["plotdata.smooth"]], paste("data/p4/temporalv_initial_prototype3_calib_ph_d1v1_plotdata_smooth_", gender, 
                                                        "_nk", age_knots, 
                                                        "_caltime", as.numeric(caltime_in),
                                                        "_tfup", t_fup, 
                                                        "_teval", t_eval, ".rds", sep = ""))
          saveRDS(calib_MxM[["plotdata.grouped"]], paste("data/p4/temporalv_initial_prototype3_calib_ph_d1v1_plotdata_grouped_", gender, 
                                                         "_nk", age_knots, 
                                                         "_caltime", as.numeric(caltime_in),
                                                         "_tfup", t_fup, 
                                                         "_teval", t_eval, ".rds", sep = ""))
          
          ### print ICI, E50, E90
          print(paste("ICI = ", calib_MxM[["ICI"]]))
          print(paste("E50 = ", calib_MxM[["E50"]]))
          print(paste("E90 = ", calib_MxM[["E90"]]))
          saveRDS(calib_MxM[["ICI"]], paste("data/p4/temporalv_initial_prototype3_calib_ph_ICI_d1v1_", gender, 
                                            "_nk", age_knots, 
                                            "_caltime", as.numeric(caltime_in),
                                            "_tfup", t_fup, 
                                            "_teval", t_eval, ".rds", sep = ""))
          saveRDS(calib_MxM[["E50"]], paste("data/p4/temporalv_initial_prototype3_calib_ph_E50_d1v1_", gender, 
                                            "_nk", age_knots, 
                                            "_caltime", as.numeric(caltime_in),
                                            "_tfup", t_fup, 
                                            "_teval", t_eval, ".rds", sep = ""))
          saveRDS(calib_MxM[["E90"]], paste("data/p4/temporalv_initial_prototype3_calib_ph_E90_d1v1_", gender, 
                                            "_nk", age_knots, 
                                            "_caltime", as.numeric(caltime_in),
                                            "_tfup", t_fup, 
                                            "_teval", t_eval, ".rds", sep = ""))
          
          ######################
          ### Discrimination ###
          ######################
          
          ### Estimate discrimination curves for development=1 and validation=1 datasets
          discrim_d1v1 <- est_discrim_d1v1(gender = gender, t_fup = t_fup, t_eval = t_eval, age_knots = age_knots, caltime_in = caltime_in)
          saveRDS(discrim_d1v1, paste("data/p4/temporalv_initial_prototype3_discrim_d1v1_", gender, 
                                      "_nk", age_knots, 
                                      "_caltime", as.numeric(caltime_in),
                                      "_tfup", t_fup, 
                                      "_teval", t_eval, ".rds", sep = ""))
          print(discrim_d1v1)
          
        }
      }
    }
  }
}

print(paste("FINISHED", Sys.time()))