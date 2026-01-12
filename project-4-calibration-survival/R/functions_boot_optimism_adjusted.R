#######################################
### Functions for optimism adjusted ###
#######################################

###
### Each of the functions beginning "boot_func_..." estimate a calibration curve in the original dataset, 
### for a model fitted in a bootstrapped dataset. The calibration curves are estimated using one of the four 
### methods (ph, pv, pv-ipcw or ipcw). They call on other functions to estimate the calibration curves (see functions.R),
### but are written to be compatible with the boot::boot function for bootstrapping, 
### meaning they can be used to estimate optimism adjusted survival curves using the approach of Steyerberg.
### See "A Practical Approach to Development, Validation and Updating. 2019" 
### Section 5.3.3 "Bootstrapping for Prediction: Optimism Correction
###
### They output a vector of predicted-observed values. It is important to
### specify pred.plot.range, so the vector of predicted-observed values, corresponds to the same vector
### of predicted values for each bootstrap iteration, meaning we can calculate the mean over the bootstrapped
### curves.
###

####################################
### Proportional Hazards approach ###
####################################
boot_func_for_optimism_adjusted_ph <- function(data, indices, fit, t, nk, pred.plot.range = NULL){
  
  ### Create bootstrapped data
  data.b <- data[indices, ]
  
  ### Fit model in bootstrapped dataset
  fit.b <- coxph(fit$formula, data = data.b, model = TRUE)
  ### Estimate cumulative baseline hazard
  bhaz.b <- basehaz(fit.b, centered = TRUE)
  
  ### Assess calibration of model in original dataset using existing function
  calib_object <- est_calib_ph(data = data, 
                               fit = fit.b, 
                               bhaz = bhaz.b, 
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
boot_func_for_optimism_adjusted_pv <- function(data, indices, fit, t, nk, pred.plot.range, split.n.groups = 1){
  
  ### Create bootstrapped data
  data.b <- data[indices, ]
  
  ### Fit model in bootstrapped dataset
  fit.b <- coxph(fit$formula, data = data.b, model = TRUE)
  ### Estimate cumulative baseline hazard
  bhaz.b <- basehaz(fit.b, centered = TRUE)
  
  ### Assess calibration of model in original dataset using existing function
  calib_object <- est_calib_pv(data = data, 
                               fit = fit.b, 
                               bhaz = bhaz.b, 
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
boot_func_for_optimism_adjusted_pv_ipcw <- function(data, indices, fit, t, nk, ipcw.formula, pred.plot.range = NULL){
  
  ### Create bootstrapped data
  data.b <- data[indices, ]
  
  ### Fit model in bootstrapped dataset
  fit.b <- coxph(fit$formula, data = data.b, model = TRUE)
  ### Estimate cumulative baseline hazard
  bhaz.b <- basehaz(fit.b, centered = TRUE)
  
  ### Assess calibration of model in original dataset using existing function
  calib_object <- est_calib_pv_ipcw(data = data,
                                    fit  = fit.b, 
                                    bhaz = bhaz.b, 
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
boot_func_for_optimism_adjusted_ipcw <- function(data, indices, fit, t, nk, ipcw.formula, pred.plot.range){
  
  ### Create bootstrapped data
  data.b <- data[indices, ]
  
  ### Fit model in bootstrapped dataset
  fit.b <- coxph(fit$formula, data = data.b, model = TRUE)
  ### Estimate cumulative baseline hazard
  bhaz.b <- basehaz(fit.b, centered = TRUE)
  
  ### Assess calibration of model in original dataset using existing function
  calib_object <- est_calib_ipcw(data = data,
                                 fit  = fit.b, 
                                 bhaz = bhaz.b, 
                                 t = t, 
                                 nk = nk,
                                 ipcw.formula = ipcw.formula,
                                 pred.plot.range = pred.plot.range)
  
  ### Extract the plot data which is of interest (predicted observed values)
  pred.obs <- calib_object[["plot"]]$data$pred.obs
  
  ### Return
  return(pred.obs)
  
}

#############################################################################################
### Wrapper function to estimate optimism adjusted calibration curve using boostrapping,  ###
### for any of the above methods                                                          ###
#############################################################################################
est_calib_optimism_adjusted_boot_wrapper <- function(calib.method = c("ph", "pv", "pv_ipcw", "ipcw"), 
                                                     data, 
                                                     fit, 
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
  # t = t.eval
  # nk = 3
  # pred.plot.range = pred.plot.range
  # CI = 95
  # R.boot = 100
  
  ### Match.arg
  calib.method <- match.arg(calib.method)
  
  ### Run bootstrapping
  if (calib.method == "ph"){
    boot.pred.obs <- boot::boot(data = data, 
                                statistic = boot_func_for_optimism_adjusted_ph, 
                                R = R.boot, fit = fit, t = t, nk = nk, 
                                pred.plot.range = pred.plot.range)
  } else if (calib.method == "pv"){
    boot.pred.obs <- boot::boot(data = data, 
                                statistic = boot_func_for_optimism_adjusted_pv, 
                                R = R.boot, fit = fit, t = t, nk = nk,
                                split.n.groups = pv.split.n.groups,
                                pred.plot.range = pred.plot.range)
  } else if (calib.method == "pv_ipcw"){
    boot.pred.obs <- boot::boot(data = data, 
                                statistic = boot_func_for_optimism_adjusted_pv_ipcw, 
                                R = R.boot, fit = fit, t = t, nk = nk,
                                ipcw.formula = ipcw.formula,
                                pred.plot.range = pred.plot.range)
  } else if (calib.method == "ipcw"){
    boot.pred.obs <- boot::boot(data = data, 
                                statistic = boot_func_for_optimism_adjusted_ipcw, 
                                R = R.boot, fit = fit, t = t, nk = nk,
                                ipcw.formula = ipcw.formula,
                                pred.plot.range = pred.plot.range)
  } 
  
  ### Get the p.025 and p975
  # alpha <- (1-CI/100)/2
  # lower <- apply(boot.pred.obs$t, 2, quantile, probs = alpha, na.rm = TRUE)
  # upper <- apply(boot.pred.obs$t, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
  
  ### Get the mean
  pred.obs.mean <- apply(boot.pred.obs$t, 2, mean, na.rm = TRUE)
  
  ### Create output object
  plot.data <- data.frame("pred" = pred.plot.range, "pred.obs" =  pred.obs.mean
                          #"pred.obs.lower" = lower, "pred.obs.upper" = upper)
  )
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") + 
    ggplot2::ggtitle(calib.method) + 
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "plotdata" = plot.data)
  
  return(output.object)
  
}


###
###
### OLD

###
###
###

# ###
# ### A wrapper, than can be used to estimate optimism adjusted performance metrics
# ###
# est_calib_optimism <- function(calib.method, data, fit, bhaz, t, surv = NULL, nk = 4, 
#                                pred.plot.range = NULL, CI, R.boot, 
#                                cens.max.follow = NULL, cens.formula, split.n.groups = 1){
#   # library(survival)
#   # calib.method = "ph"
#   # data = dat.devel
#   # fit = fit.qu
#   # t = t.eval
#   # nk = 4
#   # R.boot = 40
#   # pred.plot.range = NULL
#   
#   # ### Get the survival probabilities
#   # if (is.null(surv)){
#   #   data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
#   # } else {
#   #   data$surv <- as.numeric(surv)
#   # }
#   # 
#   # ### Add complementary log-log of predicted survival probabilities to data
#   # data$pred <- 1 - data$surv
#   # data$cloglog <- log(-log(data$surv))
#   # ### Add log of follow-up time to data
#   # data$fup_py <- log(data$time)
#   # 
#   # calib.method = "pv"
#   # data = dat.devel
#   # fit = fit.ph
#   # t = t.eval
#   # nk = 4
#   # R.boot = 50
#   # pred.plot.range <- NULL
#   ### Get the survival probabilities to estimate pred.plot.range
#   if (is.null(pred.plot.range)){
#     if (class(fit) == "coxph"){
#       pred <- 1 - as.numeric(est_surv(newdata = data, fit = fit, t = t))
#       pred.plot.range <- seq(quantile(pred, p = 0.01), quantile(pred, p = 0.99), length = 1000)
#     } else if (class(fit) == "flexsurvreg"){
#       pred <- dplyr::pull(predict(fit, newdata = data, times = t, type = "survival"), .pred_survival)
#       pred.plot.range <- seq(quantile(pred, p = 0.01), quantile(pred, p = 0.99), length = 1000)
#     }
#   }
#   
#   ### Run bootstrapping
#   if (calib.method == "ph"){
#     boot.pred.obs <- boot::boot(data = data, 
#                                 statistic = est_calib_optimism_ph, 
#                                 R = R.boot, fit = fit, t = t, nk = nk, pred.plot.range = pred.plot.range)
#   } else if (calib.method == "pois"){
#     # boot.pred.obs <- boot::boot(data = data, 
#     #                             statistic = est_calib_pois_boot, 
#     #                             R = R.boot, fit = fit, bhaz = bhaz, t = t, nk = nk, pred.plot.range = pred.plot.range)
#   } else if (calib.method == "ipcw"){
#     # boot.pred.obs <- boot::boot(data = data, 
#     #                             statistic = est_calib_ipcw_boot, 
#     #                             R = R.boot, fit = fit, bhaz = bhaz, t = t, nk = nk, pred.plot.range = pred.plot.range, 
#     #                             cens.max.follow = cens.max.follow, cens.formula = cens.formula)
#   } else if (calib.method == "pv"){
#     boot.pred.obs <- boot::boot(data = data, 
#                                 statistic = est_calib_optimism_pv, 
#                                 R = R.boot, fit = fit, t = t, nk = nk, pred.plot.range = pred.plot.range, 
#                                 split.n.groups = split.n.groups)
#   }
#   
#   ### Get the mean
#   mean <- apply(boot.pred.obs$t, 2, mean, na.rm = TRUE)
#   
#   ### Create plotdata
#   plotdata <- data.frame("pred" = pred.plot.range, "pred.obs" = mean)
#   
#   ### Create data for marginal density plot
#   margdata <- data.frame("pred" = pred, "y" = pred)
#   
#   ### Create plot
#   plot <- ggplot2::ggplot(data = plotdata) +
#     ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
#     ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
#     ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
#     ggplot2::geom_point(ggplot2::aes(x = pred, y = y), col = grDevices::rgb(0,0,0,alpha=0), data = margdata)
#   
#   ### Add marginal density
#   plot <- ggExtra::ggMarginal(plot, type = "density", margins = "x", size = 6)
#   
#   ### Create output object
#   output.object <- list("pred" = pred,
#                         "plotdata" = plotdata,
#                         "plot" = plot)
#   
#   return(output.object)
#   
# }

# ###
# ### Proportional-hazards approach
# ### 
# est_calib_optimism_ph_OLD <- function(data, indices, fit, t, nk, pred.plot.range){
#   
#   # data <- dat.devel
#   # indices <- sample(1:nrow(dat.devel), nrow(dat.devel), replace = TRUE)
#   # fit = fit
#   # t = t
#   # nk = nk
#   # pred.plot.range = pred.plot.range
#   # str(fit)
#   
#   ### Create bootstrapped data
#   data.b <- data[indices, ]
#   
#   ### Fit model in bootstrapped dataset, and get survival probabilities in original dataset
#   if (class(fit) == "coxph"){
#     ### Fit model in bootstrapped dataset
#     fit.b <- coxph(fit$formula, data = data.b, model = TRUE)
#     ### Estimate cumulative baseline hazard
#     bhaz.b <- basehaz(fit.b, centered = TRUE)
#     ### Calculate survival probabilities in original dataset
#     data$surv <- as.numeric(est_surv(newdata = data, fit = fit.b, bhaz = bhaz.b, t = t))
#   } else if (class(fit) == "flexsurvreg"){
#     ## Fit model in bootstrapped dataset
#     fit.b <- flexsurv::flexsurvreg(fit$all.formulae$scale, data = data.b, dist = "weibull")
#     ### Calculate survival probabilities in original dataset
#     data$surv <- dplyr::pull(predict(fit.b, newdata = data, times = t, type = "survival"), .pred_survival)
#   }
#   
#   ###
#   ### Assess calibration of model in original dataset
#   ###
#   
#   ### Add cloglog
#   data$cloglog = log(-log(data$surv))
#   
#   ### Fit calibration model
#   fit.calib <- survival::coxph(survival::Surv(time, status) ~ rms::rcs(cloglog, nk), data = data)
#   bhaz.calib <- survival::basehaz(fit.calib, centered = TRUE)
#   
#   ### Generate predicted observed values over the range of values pred.plot.range
#   ### Create temporary data frame
#   tmp.data <- data.frame("pred" = pred.plot.range, cloglog = log(-log(1 - pred.plot.range)))
#   
#   ### Calculate predicted observed values over pred.plot.range
#   pred.obs <- 1 - est_surv(newdata = tmp.data, fit = fit.calib, bhaz = bhaz.calib, t = t)
#   
#   ### Return
#   return(pred.obs)
#   
# }


# est_calib_optimism_pv_OLD <- function(data, indices, fit, bhaz, t, nk, pred.plot.range, split.n.groups = 1){
#   
#   #   data <- dat.devel
#   #   indices <- 1:nrow(data)
#   #   # indices <- sample(1:nrow(dat.devel), nrow(dat.devel), replace = TRUE)
#   #   fit = fit
#   #   t = t
#   #   nk = nk
#   #   pred.plot.range = pred.plot.range
#   # split.n.groups <- 1
#   # cens.max.follow = NULL
#   # cens.formula = as.formula("Surv(cens_time, cens_indicator) ~ x1+x2+x3+x4+x5")
#   # split.n.groups = 3
#   # indices <- sample(1:nrow(data), nrow(data), replace = TRUE)
#   # 
#   # ### Get the survival probabilities
#   # if (is.null(surv)){
#   #   data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
#   # } else {
#   #   data$surv <- as.numeric(surv)
#   # }
#   # 
#   # ### Add complementary log-log of predicted survival probabilities to data
#   # data$pred <- 1 - data$surv
#   # data$cloglog <- log(-log(data$surv))
#   # ### Add log of follow-up time to data
#   # data$fup_py <- log(data$time)
#   
#   ### Create bootstrapped data
#   data.b <- data[indices, ]
#   
#   ### Fit model in bootstrapped dataset, and get survival probabilities in original dataset
#   if (class(fit) == "coxph"){
#     ### Fit model in bootstrapped dataset
#     fit.b <- coxph(fit$formula, data = data.b, model = TRUE)
#     ### Estimate cumulative baseline hazard
#     bhaz.b <- basehaz(fit.b, centered = TRUE)
#     ### Calculate survival probabilities in original dataset
#     data$surv <- as.numeric(est_surv(newdata = data, fit = fit.b, bhaz = bhaz.b, t = t))
#   } else if (class(fit) == "flexsurvreg"){
#     ## Fit model in bootstrapped dataset
#     fit.b <- flexsurv::flexsurvreg(fit$all.formulae$scale, data = data.b, dist = "weibull")
#     ### Calculate survival probabilities in original dataset
#     data$surv <- dplyr::pull(predict(fit.b, newdata = data, times = t, type = "survival"), .pred_survival)
#   }
#   
#   ###
#   ### Assess calibration of model in original dataset
#   ###
#   
#   ### Calculate pred
#   data$pred <- 1 - data$surv
#   
#   ### Create new id variable
#   data$id <- 1:nrow(data)
#   
#   ### Split data by surv
#   data.pv.split <- split(data,
#                          cut(data$surv,
#                              c(-Inf, quantile(data$surv, probs = 1:split.n.groups/split.n.groups))))
#   
#   ### Calculate pseudo-values for each data split (note est_pv calculates pseudo-value for survival prob, so want to take 1 - pv)
#   pv <- lapply(data.pv.split, est_pv, t = t)
#   pvs <- 1 - do.call("c", pv)
#   ### Also get the corresponding ids
#   ids <- lapply(data.pv.split, function(x) {x$id})
#   ids <- do.call("c", ids)
#   
#   ### Combine
#   pvs <- data.frame("pv" = pvs, "id" = ids)
#   pvs <- dplyr::arrange(pvs, id)
#   
#   ### Add the pseudo-values to validation dataset
#   data$pv <- pvs$pv
#   
#   ### Transform est.surv onto logit scale
#   data$pred.logit <- log(data$pred/(1-data$pred))
#   
#   ### Fit the model using logit link function
#   calib.model.pv <- stats::glm(pv ~ rms::rcs(pred.logit, nk), 
#                                data = data, 
#                                family = stats::gaussian(link = "logit"),
#                                start = c(0,0,0,0))
#   
#   ### Generate predicted observed values over the range of values pred.plot.range
#   
#   ### Create temporary data frame
#   tmp.data <- data.frame("pred" = pred.plot.range, "pred.logit" = log(pred.plot.range/(1-pred.plot.range)))
#   
#   ### Estimate predict observed values
#   pred.obs <- predict(calib.model.pv, newdata = tmp.data, type = "response")
#   
#   ### Return
#   return(pred.obs)
#   
# }