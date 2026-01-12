###
### Program to estimate calibration in 1 development and validation dataset (split sample)
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

### Define HR for offsets
lnHR_statins <- readRDS("data/offsets_total_lnHR_statins.rds")
lnHR_ah <- readRDS("data/offsets_total_lnHR_ah.rds")
lnHR_sbp <- readRDS("data/offsets_direct_lnHR_sbp.rds")
lnHR_bmi <- readRDS("data/offsets_direct_lnHR_bmi.rds")
lnHR_nonhdl <- readRDS("data/offsets_direct_lnHR_nonhdl.rds")

###
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_calib_ph <- function(gender, model){
  
  # gender <- 1
  # model <- 1
  # # gender <- 2
  # # plot.range = NULL
  # test <- readRDS(paste("data/df_imp_devel_", gender, sep = ""))
  
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
  
  ##########################
  ### Assess calibration ###
  ##########################
  
  ##############################
  ### PH regression approach ###
  ##############################
  print(paste("PH", Sys.time()))
  calib_ph <- est_calib_plot(data = df_valid_pred, 
                             fit = fit, 
                             bhaz = bhaz, 
                             time = round(10*365.25))
  
  #####################################
  ### Get CI's for ICI, E50 and E90 ###
  #####################################
  print(paste("PH boot", Sys.time()))
  calib_ph_boot <- boot::boot(data = df_valid_pred, 
                         statistic = boot_func_for_calib_metrics_CI, 
                         R = 500, fit = fit, bhaz = bhaz, time = round(10*365.25), nk = 4)
  
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
  df_calib_smooth <- calib_ph[["plotdata"]]
  df_calib_grouped <- calib_km_group[["plot"]]$data
  saveRDS(df_calib_smooth, paste("data/calib_ph_df_smooth_", gender, "_model", model, ".rds", sep = ""))
  saveRDS(df_calib_grouped, paste("data/calib_ph_df_grouped_", gender, "_model", model, ".rds", sep = ""))
  
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
  if (model == 0){
    ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Non-causal model", sep = "")) 
  } else if (model == 1){
    ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Treatment offset model", sep = "")) 
  } else if (model == 2){
    ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Unexposed mediator model", sep = "")) 
  } else if (model == 3){
    ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Modifiable risk factor model", sep = "")) 
  } else if (model == 4){
    ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Two-component model", sep = "")) 
  } else if (model == 5){
    ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Modifiable risk factor model (SBP only)", sep = "")) 
  } else if (model == 6){
    ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Modifiable risk factor model (BMI only)", sep = "")) 
  } else if (model == 7){
    ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Modifiable risk factor model (non-HDL cholesterol only)", sep = "")) 
  }
  
  ### Save plot
  ragg::agg_png(paste("figures/calib_ph_", gender, "_model", model, ".png", sep = ""), 
                width = 1, height = 1, scaling = 1/5, unit = "in", res = 600)
  plot(ggplot_comb)
  dev.off()
  
  ### Save higher res for manuscript figures
  if (gender == 2){
    if (model == 0){
      ragg::agg_png(paste("figures/Figure1.png", sep = ""), 
                    width = 5, height = 5, scaling = 1, unit = "in", res = 600)
      plot(ggplot_comb)
      dev.off()
    } else if (model == 1){
      ragg::agg_png(paste("figures/Figure2.png", sep = ""), 
                    width = 5, height = 5, scaling = 1, unit = "in", res = 600)
      plot(ggplot_comb)
      dev.off()
    } else if (model == 2){
      ragg::agg_png(paste("figures/Figure3.png", sep = ""), 
                    width = 5, height = 5, scaling = 1, unit = "in", res = 600)
      plot(ggplot_comb)
      dev.off()
    } else if (model == 3){
      ragg::agg_png(paste("figures/Figure4.png", sep = ""), 
                    width = 5, height = 5, scaling = 1, unit = "in", res = 600)
      plot(ggplot_comb)
      dev.off()
    } else if (model == 4){
      ragg::agg_png(paste("figures/Figure5.png", sep = ""), 
                    width = 5, height = 5, scaling = 1, unit = "in", res = 600)
      plot(ggplot_comb)
      dev.off()
    } else if (model == 5){
      ragg::agg_png(paste("figures/Figure6.png", sep = ""), 
                    width = 5, height = 5, scaling = 1, unit = "in", res = 600)
      plot(ggplot_comb)
      dev.off()
    } else if (model == 6){
      ragg::agg_png(paste("figures/Figure7.png", sep = ""), 
                    width = 5, height = 5, scaling = 1, unit = "in", res = 600)
      plot(ggplot_comb)
      dev.off()
    } else if (model == 7){
      ragg::agg_png(paste("figures/Figure8.png", sep = ""), 
                    width = 5, height = 5, scaling = 1, unit = "in", res = 600)
      plot(ggplot_comb)
      dev.off()
    }
  }
  
  ### Print and save ICI, E50, E90 and the CI
  print(paste("ICI = ", calib_ph[["ICI"]]))
  print(paste("E50 = ", calib_ph[["E50"]]))
  print(paste("E90 = ", calib_ph[["E90"]]))
  saveRDS(calib_ph[["ICI"]], paste("data/calib_ph_ICI_", gender, "_model", model, ".rds", sep = ""))
  saveRDS(calib_ph[["E50"]], paste("data/calib_ph_E50_", gender, "_model", model, ".rds", sep = ""))
  saveRDS(calib_ph[["E90"]], paste("data/calib_ph_E90_", gender, "_model", model, ".rds", sep = ""))
  saveRDS(calib_ph_boot$t, paste("data/calib_ph_CI_ICI_E50_E90_", gender, "_model", model, ".rds", sep = ""))
  
  
}

### Run this function
for (gender_in in c(2,1)){
  
  print(paste("gender = ", gender_in))
  
  for (model_in in c(0,1,2,3,4,5,6,7)){
    
    print(paste("model = ", model_in))
    est_calib_ph(gender = gender_in, model = model_in)
    
  }
}

warnings()

### Create a combined plot for models 0 to 4
create_combined_calib_plot <- function(gender){
  
  ### Create a list to store the plots in
  output_plot_list <- vector("list", 5)
  
  ### For each model, create a plot and store it
  for (model in 0:4){
    
    ### Read in data
    df_calib_smooth <- readRDS(paste("data/calib_ph_df_smooth_", gender, "_model", model, ".rds", sep = ""))
    df_calib_grouped <- readRDS(paste("data/calib_ph_df_grouped_", gender, "_model", model, ".rds", sep = ""))
    
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
    if (model == 0){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Non-causal model", sep = "")) 
    } else if (model == 1){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Treatment offset model", sep = "")) 
    } else if (model == 2){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Unexposed mediator model", sep = "")) 
    } else if (model == 3){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Modifiable risk factor model", sep = "")) 
    } else if (model == 4){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(paste("Model ", model, ": Two-component model", sep = "")) 
    }
    
    ### Assign to output
    output_plot_list[[(model+1)]] <- ggplot_comb
  }
  
  ### Combine the plots into a single ggplot
  plot_out <- ggpubr::ggarrange(plotlist = output_plot_list, nrow = 3, ncol = 2)
  
  ### Save
  ragg::agg_png(paste("figures/calib_plot_grid", gender, ".png", sep = ""), 
                width = 10, height = 15, scaling = 1, unit = "in", res = 600)
  plot(plot_out)
  dev.off()
}

### Run this function
for (gender_in in c(2,1)){
  
  print(paste("gender = ", gender_in))
  create_combined_calib_plot(gender = gender_in)
  
}

print(paste("FINISHED", Sys.time()))

warnings()