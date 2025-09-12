###
### Program to estimate calibration in 1 development and validation dataset (split sample)
###

### Calibration curve produced seperately for each development/validation dataset combination
### This is primarily because for each development dataset, we get a different set of counter-factual survival times in the validation dataset

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

###
### Write a function that will produce a calibration curve for a given development and validation dataset
###
est_calib_ph_d1v1 <- function(gender, 
                              imps.valid, # Having this as separate argument to speed up computationally
                              age_knots,
                              caltime,
                              plot.range = NULL){
  
  # gender <- 2
  # plot.range = NULL
  
  ### Read in model and counter factual survival times
  if (age_knots == 3){
    ### Read in model
    fit <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
    ### Read in counterfactual survival times
    cf.surv.times <- readRDS(paste("data/p4/prototype3_cf_surv_times_", gender,
                                   "_devel", 1,
                                   "_valid", 1, ".rds", sep = ""))
  } else if (age_knots == 4){
    
    if (caltime == TRUE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_caltime_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     ".rds", sep = ""))
    }
    
    if (caltime == FALSE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     ".rds", sep = ""))
    }
  }
  
  ### Pick the first validation dataset
  data.valid <- imps.valid[[1]]
  
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
                             time = round(10*365.25),
                             fit.type = "prototype3",
                             type = "single")
  
  ###########################
  ### KM grouped approach ###
  ###########################
  print(paste("KM group", Sys.time()))
  calib.km.group <- est_calib_plot_group(data = data.valid.pred, 
                                         fit = fit, 
                                         bhaz = bhaz,  
                                         time = round(10*365.25),
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
                        ggplot2::aes(x = pred, y = obs, col = grDevices::rgb(0, 1, 0, alpha = 1))) + 
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle("Calibration of prototype 3") + ggplot2::xlim(c(0,0.6)) + ggplot2::ylim(c(0,0.6))
  
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
### Write a function that will produce a calibration curve for a given development and validation dataset, within subgroups defined by the by variable
###
est_calib_ph_d1v1_byvar <- function(gender, 
                                    imps.valid, # Having this as separate argument to speed up computationally
                                    plot.range = NULL,
                                    age_knots,
                                    caltime,
                                    byvar){
  
  # gender <- 2
  # plot.range = NULL
  # byvar <- "ethnicity"
  
  ### Read in model and counter factual survival times
  if (age_knots == 3){
    ### Read in model
    fit <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
    ### Read in counterfactual survival times
    cf.surv.times <- readRDS(paste("data/p4/prototype3_cf_surv_times_", gender,
                                   "_devel", 1,
                                   "_valid", 1, ".rds", sep = ""))
  } else if (age_knots == 4){
    
    if (caltime == TRUE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_caltime_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     ".rds", sep = ""))
    }
    
    if (caltime == FALSE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
      
      ### Read in counterfactual survival times
      cf.surv.times <- readRDS(paste("data/p4/prototype3_4knots_cf_surv_times_", gender, 
                                     "_devel", 1, 
                                     "_valid", 1, 
                                     ".rds", sep = ""))
    }
  }
  
  ### Pick the first validation dataset
  data.valid <- imps.valid[[1]]
  
  ### Remove people with missing region variable
  if (byvar == "region"){
    data.valid <- subset(data.valid, !is.na(region))
  }
  
  ### Replace cvd_time with counterfactual survival times
  data.valid <- dplyr::select(data.valid, -cvd_time)
  data.valid <- merge(data.valid, cf.surv.times, by.x = "patid", by.y = "patid") |>
    dplyr::rename(cvd_time = cvd_time_cf)
  
  ### Create variable with forms that can be used for prediction, and agrees with 'prescribed treatment strategy'
  ### which is what the individual is doing
  data.valid.pred <- 
    dplyr::mutate(data.valid,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  ### Subset within subgroups of interest
  bylevels <- levels(data.valid.pred[[byvar]])
  data.valid.pred.list <- lapply(bylevels, function(x) {data.valid.pred[data.valid.pred[,byvar] == x, ]})
  
  ##########################
  ### Assess calibration ###
  ##########################
  
  ##############################
  ### PH regression approach ###
  ##############################
  print(paste("PH", Sys.time()))
  calib.ph <- lapply(1:length(data.valid.pred.list), function(x) {
    est_calib_plot(data = data.valid.pred.list[[x]], 
                   fit = fit, 
                   bhaz = bhaz, 
                   time = round(10*365.25),
                   fit.type = "prototype3",
                   type = "single")
  })
  
  ###########################
  ### KM grouped approach ###
  ###########################
  print(paste("KM group", Sys.time()))
  calib.km.group <- lapply(1:length(data.valid.pred.list), function(x) {
    est_calib_plot_group(data = data.valid.pred.list[[x]], 
                         fit = fit, 
                         bhaz = bhaz,  
                         time = round(10*365.25),
                         n.groups = 10,
                         fit.type = "prototype3",
                         type = "single",
                         single.CI = TRUE)
  })
  
  ##############################################
  ### Function to combine these into ggplots ###
  ##############################################
  print("combine into single plot")
  ggplot.comb <- lapply(1:length(data.valid.pred.list), function(x){
    
    ### Extract plot data
    df.calib.smooth <- calib.ph[[x]][["plotdata"]]
    df.calib.grouped <- calib.km.group[[x]][["plot"]]$data
    
    ### Create a ggplot out of this
    ggplot.out <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), 
                         data = df.calib.smooth) + 
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
      ggplot2::geom_point(data = df.calib.smooth,  
                          ggplot2::aes(x = pred, y = pred.obs), col = grDevices::rgb(0, 0, 0, alpha = 0)) +
      ggplot2::geom_point(data = df.calib.grouped, 
                          ggplot2::aes(x = pred, y = obs, col = grDevices::rgb(0, 1, 0, alpha = 1))) + 
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(paste0(toupper(substr(bylevels[[x]], 1, 1)), substring(bylevels[[x]], 2)))
    
    if (byvar == "ethnicity"){
      ggplot.out <- ggplot.out + ggplot2::xlim(c(0,0.4)) + ggplot2::ylim(c(0,0.4))
    }
    
    if (byvar == "age_cat"){
      p99 <- as.numeric(quantile(df.calib.smooth$pred, p = 0.99))
      ggplot.out <- ggplot.out + ggplot2::xlim(c(0,p99)) + ggplot2::ylim(c(0,p99))
    }
    
    return(ggplot.out)
    
  })
  
  ####################################################
  ### Function to extract list of ICI, E50 and E90 ###
  ####################################################
  print("create tables")
  
  ### Get tables
  ICI_out <- lapply(1:length(data.valid.pred.list), function(x){calib.ph[[x]][["ICI"]]})
  E50_out <- lapply(1:length(data.valid.pred.list), function(x){calib.ph[[x]][["E50"]]})
  E90_out <- lapply(1:length(data.valid.pred.list), function(x){calib.ph[[x]][["E90"]]})
  names(ICI_out) <- paste0(toupper(substr(bylevels, 1, 1)), substring(bylevels, 2))
  names(E50_out) <- paste0(toupper(substr(bylevels, 1, 1)), substring(bylevels, 2))
  names(E90_out) <- paste0(toupper(substr(bylevels, 1, 1)), substring(bylevels, 2))
  
  ## Get N
  N_out <- sapply(data.valid.pred.list, nrow)
  
  ### Create output.object
  output.object <- list("plot" = ggplot.comb,
                        "ICI" = ICI_out,
                        "E50" = E50_out,
                        "E90" =  E90_out,
                        "N" = N_out)
  
  ### Return output object
  return(output.object)
  
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
  
  ### Create age variable
  imps.valid <- lapply(imps.valid, function(x) {
    x$age_cat <- cut(x$age, breaks = c(17,30,40,50,60,70,80,Inf))
    return(x)
  })
  
  ###
  ### Load the region information
  region_df <- readRDS("../../Aurum_Jun2021_extract/data/extraction/cohort_exclu3.rds") |>
    dplyr::select(patid, region)
  
  ### Read in lookup
  region_lookup <- read.table("../../Aurum_Jun2021_extract/data/unzip/zLookups/region.txt", header = TRUE, sep = "\t") |>
    dplyr::rename(region = regionid, region_desc = Description)
  
  region_df <- dplyr::left_join(region_df, region_lookup, by = dplyr::join_by(region)) |>
    dplyr::select(-region) |>
    dplyr::rename(region = region_desc) |>
    dplyr::mutate(region = as.factor(region))
  
  ### Create region variable
  imps.valid <- lapply(imps.valid, function(x) {
    x <- dplyr::left_join(x, region_df, by = dplyr::join_by(patid))
    return(x)
  })
  
  #####################
  ### Entire cohort ###
  #####################
  
  ### Estimate calibration curves for development=1 and validation=1 datasets
  calib_MxM <- est_calib_ph_d1v1(gender, imps.valid = imps.valid, age_knots = age_knots, caltime = caltime)
  
  ### Create
  ### Save plot
  ragg::agg_png(paste("figures/p4/prototype3_calib_ph_d1v1_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".png", sep = ""), 
                width = 1, height = 1, scaling = 1/5, unit = "in", res = 600)
  plot(calib_MxM[["plot"]])
  dev.off()
  
  ### print ICI, E50, E90
  print(paste("ICI = ", calib_MxM[["ICI"]]))
  print(paste("E50 = ", calib_MxM[["E50"]]))
  print(paste("E90 = ", calib_MxM[["E90"]]))
  saveRDS(calib_MxM[["ICI"]], paste("data/p4/prototype3_calib_ph_ICI_d1v1_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  saveRDS(calib_MxM[["E50"]], paste("data/p4/prototype3_calib_ph_E50_d1v1_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  saveRDS(calib_MxM[["E90"]], paste("data/p4/prototype3_calib_ph_E90_d1v1_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  
  ##############
  ### Region ###
  ##############
  
  ### Get data for calibration plots
  calib_data <- est_calib_ph_d1v1_byvar(gender = gender, 
                                        imps.valid = imps.valid, 
                                        age_knots = age_knots,
                                        caltime = caltime,
                                        byvar = "region")
  
  ### Combine into single image
  plots_combined <- ggpubr::ggarrange(plotlist = calib_data[["plot"]], nrow = 3, ncol = 4)
  
  ### Save plot
  ragg::agg_png(paste("figures/p4/prototype3_calib_ph_d1v1_region_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".png", sep = ""), 
                width = 12/9, height = 1, scaling = 1/9, unit = "in", res = 600)
  plot(plots_combined)
  dev.off()
  
  ### Create table for ICI, E50 and E90
  calib_data_table <- data.frame("N" = calib_data["N"],
                                 "ICI" = unlist(calib_data[["ICI"]]),
                                 "E50" = unlist(calib_data[["E50"]]),
                                 "E90" = unlist(calib_data[["E90"]]))
  
  ### Save
  saveRDS(calib_data_table, paste("data/p4/prototype3_calib_ph_table_d1v1_region_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  
  #################
  ### Ethnicity ###
  #################
  
  ### Get data for calibration plots
  calib_data <- est_calib_ph_d1v1_byvar(gender = gender, 
                                        imps.valid = imps.valid, 
                                        age_knots = age_knots,
                                        caltime = caltime,
                                        byvar = "ethnicity")
  
  ### Combine into single image
  plots_combined <- ggpubr::ggarrange(plotlist = calib_data[["plot"]], nrow = 3, ncol = 3)
  
  ### Save plot
  ragg::agg_png(paste("figures/p4/prototype3_calib_ph_d1v1_ethnicity_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".png", sep = ""), 
                width = 1, height = 1, scaling = 1/9, unit = "in", res = 600)
  plot(plots_combined)
  dev.off()
  
  ### Create table for ICI, E50 and E90
  calib_data_table <- data.frame("N" = calib_data["N"],
                                 "ICI" = unlist(calib_data[["ICI"]]),
                                 "E50" = unlist(calib_data[["E50"]]),
                                 "E90" = unlist(calib_data[["E90"]]))
  
  ### Save
  saveRDS(calib_data_table, paste("data/p4/prototype3_calib_ph_table_d1v1_ethnicity_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  
  ###########
  ### Age ###
  ###########
  
  ### Get data for calibration plots
  calib_data <- est_calib_ph_d1v1_byvar(gender = gender, 
                                        imps.valid = imps.valid, 
                                        age_knots = age_knots,
                                        caltime = caltime,
                                        byvar = "age_cat")
  
  ### Combine into single image
  plots_combined <- ggpubr::ggarrange(plotlist = calib_data[["plot"]], nrow = 3, ncol = 3)
  
  ### Save plot
  ragg::agg_png(paste("figures/p4/prototype3_calib_ph_d1v1_age_cat_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".png", sep = ""), 
                width = 1, height = 1, scaling = 1/9, unit = "in", res = 600)
  plot(plots_combined)
  dev.off()
  
  ### Create table for ICI, E50 and E90
  calib_data_table <- data.frame("N" = calib_data["N"],
                                 "ICI" = unlist(calib_data[["ICI"]]),
                                 "E50" = unlist(calib_data[["E50"]]),
                                 "E90" = unlist(calib_data[["E90"]]))
  
  ### Save
  saveRDS(calib_data_table, paste("data/p4/prototype3_calib_ph_table_d1v1_age_cat_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  
  
}

### Run this function
for (gender_in in c(1,2)){
  
  create_and_save_output(gender = gender_in, age_knots = 3, caltime = FALSE) 
  create_and_save_output(gender = gender_in, age_knots = 4, caltime = FALSE) 
  create_and_save_output(gender = gender_in, age_knots = 4, caltime = TRUE) 
  
}

print(paste("FINISHED", Sys.time()))

