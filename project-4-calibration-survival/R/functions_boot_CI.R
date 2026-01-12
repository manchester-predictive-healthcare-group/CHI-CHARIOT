#########################################
### Functions for confidence interval ###
#########################################

###
### Each of the functions beginning "boot_func_..." estimate a ccalibration curve in a bootstrapped
### dataset using one of the four methods (ph, pv, pv-ipcw or ipcw). They call on other functions to estimate the calibration curves (see functions.R),
### but are written to be compatible with the boot::boot function for bootstrapping, 
### meaning they can be used to estimate a confidence interval for the calibration curves
###
### They output a vector of predicted-observed values. It is important to
### specify pred.plot.range, so the vector of predicted-observed values, corresponds to the same vector
### of predicted values for each bootstrap iteration, meaning we can calculate 95% limits over the bootstrapped
### curves.
###

####################################
### Propotional Hazards approach ###
####################################
boot_func_for_CI_ph <- function(data, indices, fit, bhaz, t, nk, pred.plot.range = NULL){
  
  ### Create bootstrapped data
  data.b <- data[indices, ]
  
  ### Assess calibration of model in original dataset using existing function
  calib_object <- est_calib_ph(data = data.b, 
                               fit = fit, 
                               bhaz = bhaz, 
                               t = t, 
                               nk = nk, 
                               pred.plot.range = pred.plot.range)
  
  ### Extract the plot data which is of interest (predicted observed values)
  pred.obs <- calib_object[["plot"]]$data$pred.obs
  
  ### Return
  return(pred.obs)
  
}


#############################
### Pseudo-value approach ###
#############################
boot_func_for_CI_pv <- function(data, indices, fit, bhaz, t, nk, split.n.groups, pred.plot.range = NULL){
  
  ### Create bootstrapped data
  data.b <- data[indices, ]
  
  ### Assess calibration of model in original dataset using existing function
  calib_object <- est_calib_pv(data = data.b, 
                               fit  = fit,
                               bhaz = bhaz, 
                               t = t, 
                               nk = nk,
                               group.by = "predrisk",
                               split.n.groups = split.n.groups, 
                               pred.plot.range = pred.plot.range)
  
  ### Extract the plot data which is of interest (predicted observed values)
  pred.obs <- calib_object[["plot"]]$data$pred.obs
  
  ### Return
  return(pred.obs)
  
}

###################################################################
### Pseudo-value approach (with IPC weighting, binder approach) ###
###################################################################
boot_func_for_CI_pv_ipcw <- function(data, indices, fit, bhaz, t, nk, ipcw.formula, pred.plot.range = NULL){
  
  ### Create bootstrapped data
  data.b <- data[indices, ]
  
  ### Assess calibration of model in original dataset using existing function
  calib_object <- est_calib_pv_ipcw(data = data.b, 
                                    fit  = fit, 
                                    bhaz = bhaz, 
                                    t = t, 
                                    nk = nk, 
                                    ipcw.formula = ipcw.formula,
                                    pred.plot.range = pred.plot.range,
                                    eventglm.pv.method = "binder")
  
  ### Extract the plot data which is of interest (predicted observed values)
  pred.obs <- calib_object[["plot"]]$data$pred.obs
  
  ### Return
  return(pred.obs)
  
}

#####################
### IPCW approach ###
#####################
boot_func_for_CI_ipcw <- function(data, indices, fit, bhaz, t, nk, ipcw.formula, pred.plot.range){
  
  ### Create bootstrapped data
  data.b <- data[indices, ]
  
  ### Assess calibration of model in original dataset using existing function
  calib_object <- est_calib_ipcw(data = data.b,
                                 fit  = fit, 
                                 bhaz = bhaz, 
                                 t = t, 
                                 nk = nk,
                                 ipcw.formula = ipcw.formula,
                                 pred.plot.range = pred.plot.range)
  
  ### Extract the plot data which is of interest (predicted observed values)
  pred.obs <- calib_object[["plot"]]$data$pred.obs
  
  ### Return
  return(pred.obs)
  
}


################################################################################
### Wrapper function to estimate a confidence interval using bootstrapping,  ###
### for any of the above methods                                             ###
################################################################################
est_calib_CI_boot_wrapper <- function(calib.method = c("ph", "pv", "pv_ipcw", "ipcw"),  
                                      data, 
                                      fit, 
                                      bhaz, 
                                      t, 
                                      nk = 4, 
                                      pred.plot.range, 
                                      pv.split.n.groups = 1,
                                      ipcw.formula, 
                                      CI, 
                                      R.boot){
  
  # calib.method = "ph" 
  # data = dat.valid
  # fit = fit1
  # bhaz = bhaz1
  # t = t.eval
  # nk = 3 
  # pred.plot.range = pred.plot.range
  # CI = 95
  # R.boot = 10
  
  # ### Get the survival probabilities
  # if (is.null(surv)){
  #   data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  # } else {
  #   data$surv <- as.numeric(surv)
  # }
  # 
  # ### Add complementary log-log of predicted survival probabilities to data
  # data$pred <- 1 - data$surv
  # data$cloglog <- log(-log(data$surv))
  # ### Add log of follow-up time to data
  # data$fup_py <- log(data$time)
  
  ### Match.arg
  calib.method <- match.arg(calib.method)
  
  ### Run bootstrapping
  if (calib.method == "ph"){
    boot.pred.obs <- boot::boot(data = data, 
                                statistic = boot_func_for_CI_ph, 
                                R = R.boot, fit = fit, bhaz = bhaz, t = t, nk = nk, 
                                pred.plot.range = pred.plot.range)
  } else if (calib.method == "pv"){
    boot.pred.obs <- boot::boot(data = data, 
                                statistic = boot_func_for_CI_pv, 
                                R = R.boot, fit = fit, bhaz = bhaz, t = t, nk = nk,
                                split.n.groups = pv.split.n.groups,
                                pred.plot.range = pred.plot.range)
  } else if (calib.method == "pv_ipcw"){
    boot.pred.obs <- boot::boot(data = data, 
                                statistic = boot_func_for_CI_pv_ipcw, 
                                R = R.boot, fit = fit, bhaz = bhaz, t = t, nk = nk,
                                ipcw.formula = ipcw.formula,
                                pred.plot.range = pred.plot.range)
  } else if (calib.method == "ipcw"){
    boot.pred.obs <- boot::boot(data = data, 
                                statistic = boot_func_for_CI_ipcw, 
                                R = R.boot, fit = fit, bhaz = bhaz, t = t, nk = nk,
                                ipcw.formula = ipcw.formula,
                                pred.plot.range = pred.plot.range)
  } 
  
  ### Get the p.025 and p975
  alpha <- (1-CI/100)/2
  pred.obs.lower <- apply(boot.pred.obs$t, 2, quantile, probs = alpha, na.rm = TRUE)
  pred.obs.upper <- apply(boot.pred.obs$t, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
  pred.obs.median <- apply(boot.pred.obs$t, 2, quantile, probs = 0.5, na.rm = TRUE)
  
  ### Get the mean
  pred.obs.mean <- apply(boot.pred.obs$t, 2, mean, na.rm = TRUE)
  
  ### Create output object
  plot.data <- data.frame("pred" = pred.plot.range, 
                          "pred.obs.median" =  pred.obs.median,
                          "pred.obs.mean" =  pred.obs.mean,
                          "pred.obs.lower" = pred.obs.lower, 
                          "pred.obs.upper" = pred.obs.upper)
  
  # ### Create plot
  # plot <- ggplot2::ggplot(data = plot.data) +
  #   ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
  #   ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") + 
  #   ggplot2::ggtitle(calib.method) + 
  #   ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  
  # ### Create output.object
  # output.object <- list("plot" = plot,
  #                       "plotdata" = plot.data)
  
  return(plot.data)
  
}