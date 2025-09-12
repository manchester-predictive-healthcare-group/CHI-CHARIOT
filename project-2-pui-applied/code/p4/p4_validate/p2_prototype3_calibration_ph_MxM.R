###
### Program to estimate calibration
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
est_calib_MxM <- function(gender, 
                          age_knots,
                          caltime,
                          devel.num,
                          valid.num,
                          imps.valid, # Having this as separate argument to speed up computationally
                          plot.range = NULL){
  
  # fitest1 <- readRDS(paste("data/p4/prototype3_cox_", 1, "_imp", 1, ".rds", sep = ""))
  # # fitest2 <- readRDS(paste("data/p4/prototype3_cox_", 2, "_imp", 1, ".rds", sep = ""))
  # coefficients(fitest1)
  # fitest1$n
  # length(fitest1[["model"]]$patid)
  # # coefficients(fitest2)
  #   treatment_strategy <- "healthy"
  #   devel.num <- 1
  #   valid.num <- 1
  #   gender <- 1
  #   gender_char <- c("male", "female")[gender]
  #   print(paste("gender = ", gender_char))
  #   imps.valid <- readRDS(paste("data/dfs_valid_", gender_char, ".rds", sep = ""))
  
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
  
  ##########################
  ### Assess calibration ###
  ##########################
  
  ##############################
  ### PH regression approach ###
  ##############################
  print(paste(devel.num, valid.num, "PH", Sys.time()))
  calib.ph <- est_calib_plot(data = data.valid.pred, 
                             fit = fit, 
                             bhaz = bhaz, 
                             time = round(10*365.25),
                             fit.type = "prototype3",
                             type = "single",
                             plot.range = plot.range)
  
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
  
  ### Create output.object
  output.object <- list("plotdata.smooth" = calib.ph[["plotdata"]],
                        "plotdata.range" = calib.ph[["plotdata.range"]],
                        "plotdata.grouped" = calib.km.group[["plot"]]$data,
                        "ICI" = as.numeric(calib.ph[["ICI"]]),
                        "E50" = as.numeric(calib.ph[["E50"]]),
                        "E90" =  as.numeric(calib.ph[["E90"]]))
  
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
  
  ###
  ### Get the 1st and 99th percentile of predicted risks so we can plot calibration curves over the same set of points, and obtain a median curve
  ### We will do this using the first validation dataset
  ###
  
  ### Read in 1 model to obtain these
  if (age_knots == 3){
    ### Read in model
    fit <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
  } else if (age_knots == 4){
    if (caltime == TRUE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    }
    
    if (caltime == FALSE){
      ### Read in model
      fit <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    }
  }
  
  ### Pick first one
  data.valid <- imps.valid[[1]]
  
  ### Create dataset for prediction which agrees with 'prescribed treatment strategy', which is what the individual is doing (i.e. no changes)
  data.valid.pred <- 
    dplyr::mutate(data.valid,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  ### Calculate predicted risks for validation dataset 1, according to model developed in development dataset 1
  predrisk <- as.numeric(est_surv_offset_prototype3(newdata = data.valid.pred, fit = fit, bhaz = bhaz, time = round(10*365.25)))
  
  ### Use these to calculate limit for the plot
  pred.quantiles <- as.numeric(quantile(1 - predrisk, c(0.001,0.999)))
  rm(fit, bhaz, data.valid, predrisk)
  
  ###
  ### Estimate calibration curves for each combination of development and validation datasets
  ###
  
  ### Apply this functions over 10 development and validation datasets
  cl <- parallel::makeCluster(11)
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers()
  calib.data.MxM <- (foreach(x = 1:10, .combine = list, .multicombine = TRUE, .packages = c("survival"), .export = c("est_calib_MxM", "est_calib_plot", "est_calib_plot_group", "est_surv_offset_prototype3", "est_surv")) %dopar% {
    lapply(1:10, function(y) {est_calib_MxM(gender, age_knots, caltime, x, y, 
                                            imps.valid = imps.valid, 
                                            plot.range = seq(pred.quantiles[1], pred.quantiles[2], length.out = 1000))})
  })
  stopCluster(cl)
  
  print(paste("DATA OBTAINED", Sys.time()))
  saveRDS(calib.data.MxM, paste("data/p4/prototype3_calib_ph_data_MxM_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  
  ######################################################################
  ### Create a calibration plot with all these calibration curves on ###
  ### and with a median curve
  ######################################################################
  
  ###
  ### Extract and combine the smoothed plot data
  ###
  calib.curves.smooth.comb <- lapply(unlist(calib.data.MxM, recursive = FALSE), function(x) {x[["plotdata.smooth"]]})
  
  ### Add an indicator for each plot
  calib.curves.smooth.comb <- 
    lapply(1:length(calib.curves.smooth.comb), 
           function(x) {
             calib.curves.smooth.comb[[x]]$curve.num <- x
             return(calib.curves.smooth.comb[[x]])
           })
  
  ### Concatenate into a single dataset
  calib.curves.smooth.comb <- do.call("rbind", calib.curves.smooth.comb)
  
  ###
  ### Extract and combine the grouped plot data
  ###
  calib.curves.grouped.comb <- lapply(unlist(calib.data.MxM, recursive = FALSE), function(x) {x[["plotdata.grouped"]]})
  
  ### Extract pred.obs for each
  calib.curves.grouped.obs <-  lapply(calib.curves.grouped.comb, function(x) {x$obs})
  calib.curves.grouped.pred <-  lapply(calib.curves.grouped.comb, function(x) {x$pred})
  
  ### Get average of pred.obs
  obs.mean <- as.numeric(rowMeans(do.call("cbind", calib.curves.grouped.obs)))
  obs.median <- as.numeric(apply(do.call("cbind", calib.curves.grouped.obs), 1, function(x) {quantile(x, p = 0.5)}))
  pred.mean <- as.numeric(rowMeans(do.call("cbind", calib.curves.grouped.pred)))
  pred.median <- as.numeric(apply(do.call("cbind", calib.curves.grouped.pred), 1, function(x) {quantile(x, p = 0.5)}))
  
  ### Create a single dataset
  mean.grouped <- data.frame("obs" = obs.mean, 
                             "pred" = pred.mean)
  median.grouped <- data.frame("obs" = obs.median, 
                               "pred" = pred.median)
  
  ###
  ### Create median calibration curve
  ###
  calib.curves.smooth.range.comb <- lapply(unlist(calib.data.MxM, recursive = FALSE), function(x) {x[["plotdata.range"]]})
  
  ### Extract pred.obs for each
  pred.obs.range <-  lapply(calib.curves.smooth.range.comb, function(x) {x$pred.obs})
  
  ### Get average of pred.obs
  pred.obs.mean <- as.numeric(rowMeans(do.call("cbind", pred.obs.range)))
  pred.obs.median <- as.numeric(apply(do.call("cbind", pred.obs.range), 1, function(x) {quantile(x, p = 0.5)}))
  
  ### Get predicted risks
  pred.range <- calib.curves.smooth.range.comb[[1]]$pred
  
  ### Concatenate into a single dataset
  mean.curve <- data.frame("pred.obs" = pred.obs.mean, 
                           "pred" = pred.range)
  median.curve <- data.frame("pred.obs" = pred.obs.median, 
                             "pred" = pred.range)
  
  ### Add this to dataset with all the other calibration curves
  ### Get number of curves
  n.curves <- max(calib.curves.smooth.comb$curve.num)
  calib.curves.smooth.comb <- rbind(calib.curves.smooth.comb,
                                    dplyr::mutate(median.curve, curve.num = n.curves + 1))
  
  ### Create an indicator for whether its median or a normal curve
  calib.curves.smooth.comb <- dplyr::mutate(calib.curves.smooth.comb,
                                            "curve.ind" = dplyr::case_when(curve.num <= n.curves ~ "Calibration",
                                                                           curve.num > n.curves ~ "Median"))
  
  ### Define object for group colours
  curve.colours <- c("Calibration" = "red", "Median" = "blue")
  
  
  ### Create a ggplot out of this
  ggplot.comb <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, group = curve.num, color = curve.ind, alpha = curve.ind, linewidth = curve.ind), 
                       data = calib.curves.smooth.comb) + 
    ggplot2::scale_colour_manual(values = curve.colours) +
    ggplot2::scale_alpha_manual(values = c(0.25, 1), guide = "none") +
    ggplot2::scale_linewidth_manual(values = c(1, 0.5)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::xlim(c(0, max(mean.grouped$pred))) + ggplot2::ylim(c(0, max(mean.grouped$pred))) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()) +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
    ggplot2::geom_point(data = subset(calib.curves.smooth.comb, curve.num <= n.curves), 
                        ggplot2::aes(x = pred, y = pred.obs), col = grDevices::rgb(0, 0, 0, alpha = 0)) +
    ggplot2::geom_point(data = mean.grouped, 
                        ggplot2::aes(x = pred, y = obs, col = grDevices::rgb(0, 1, 0, alpha = 1))) +
    ggplot2::ggtitle(c("Male cohort", "Female cohort")[gender])
  
  # ## Add ggMarginal
  ggplot.comb <- ggExtra::ggMarginal(ggplot.comb,
                                     margins = "x",
                                     x = pred,
                                     size = 8,
                                     colour = "red")
  
  # ### Create a new plot, with no transparency, purely for getting a legend
  # ggplot.temp <- ggplot2::ggplot() +
  #   ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, group = curve.num, color = curve.ind, alpha = curve.ind), 
  #                      data = calib.curves.smooth.comb) + 
  #   ggplot2::scale_colour_manual(values = curve.colours) +
  #   ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
  #   ggplot2::geom_point(data = calib.curves.comb, ggplot2::aes(x = pred, y = pred.obs), col = grDevices::rgb(0, 0, 0, alpha = 0))
  # legend.save <- ggpubr::get_legend(ggplot.temp, position = "bottom")
  # 
  # ### Combine plot with legend
  # gglot.comb <- gridExtra::arrangeGrob(ggplot.comb, legend.save, nrow = 2, heights = c(15,1))
  
  ### Save the plot
  ragg::agg_png(paste("figures/p4/prototype3_calib_ph_MxM_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".png", sep = ""), 
                width = 1, height = 1, scaling = 1/6, unit = "in", res = 600)
  plot(ggplot.comb)
  dev.off()
  
  ###
  ### Create a table with mean, median, inter quartile range of ICI, E50 and E90
  ###
  
  ### Write a function to add a row
  create_row <- function(varname, dp = 2){
    ### Extract the plot data
    values <- as.numeric(lapply(unlist(calib.data.MxM, recursive = FALSE), function(x) {x[[varname]]}))
    
    output <- c(paste(round(mean(values), dp), 
                      " (", 
                      round(sd(values), dp),
                      ")", sep = ""),
                paste(round(median(values), dp), 
                      " (", 
                      round(quantile(values, p = 0.025), dp),
                      ",",
                      round(quantile(values, p = 0.975), dp),
                      ")", sep = ""),
                round(min(values), dp),
                round(max(values), dp))
    return(output)
  }
  
  
  ### Create an empty table
  calib.table <- data.frame(character(),
                            character(),
                            numeric(),
                            numeric())
  
  ### Add rows
  for (var in c("ICI", "E50", "E90")){
    calib.table <- rbind(calib.table, create_row(var))
  }
  rownames(calib.table) <- c("ICI", "E50", "E90")
  colnames(calib.table) <- c("Mean (sd)", "Median (p025, 0.975)", "min", "max")
  
  ### Save the table
  saveRDS(calib.table, paste("data/p4/prototype3_calib_ph_table_MxM_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  print(paste("FINISHED", gender_char, Sys.time()))
  
}

### Run this function
for (gender_in in c(1,2)){
  
  create_and_save_output(gender = gender_in, age_knots = 4, caltime = TRUE) 
  
}

print(paste("FINISHED", Sys.time()))
