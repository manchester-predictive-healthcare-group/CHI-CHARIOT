#################
### Functions ###
#################

###############################################################
### Functions to get counterfactual adjsuted survival times ###
###############################################################

###
### Function written so it can be applied in conjunction with group_modify
###
### Going to write a seperate function for models 1,2, vs 3,4,5,6,7, rather than having an if statement, for 
### computational reaons (this function will be applied 1 millon times)
###
### For models 1 and 2, we only adjust for antihypertensive use, but at baseline and changes during follow-up
###
### For models 3,4,5,6,7, we adjust for change in either statins or antihypertensives during foloow-up.
###
### Utilised in p2_calculate_cf_surv_times.R, and sensitivity analyses
###

### Write a function to extract bhaz for a specific time, used in subsequent functions
get_bhaz <- function(time, bhaz){
  if (time == 0){
    out <- 0
  } else {
    out <- as.numeric(bhaz$hazard[max(which(bhaz$time <= time))])
  } 
  return(out)
}

### Function for models 1,2
get_survtimes_adj_model12 <- function(data, id){
  
  ### Calculate cumulative hazard at transition times
  data$cumhaz_tstart <- unlist(lapply(data$tstart, get_bhaz, bhaz = bhaz))
  data$cumhaz_tstop <- unlist(lapply(data$cvd_time, get_bhaz, bhaz = bhaz))
  
  ### Intervals which are "on treatment" get multiplied by exp(B)
  data <- dplyr::mutate(data,
                        cumhaz_tstart = cumhaz_tstart*
                          exp(offset_ah_timevar_lnHR),
                        cumhaz_tstop = cumhaz_tstop*
                          exp(offset_ah_timevar_lnHR))
  
  ### Then calculate total hazard
  cumhaz_adj <- sum(data$cumhaz_tstop - data$cumhaz_tstart)
  
  ### Now see what value of t corresponds to a hazard of this
  ### Want minimum time, for which this amount of hazard is reached
  t_out <- as.numeric(bhaz$time[min(which(bhaz$hazard >= cumhaz_adj))])
  
  ### Get cvd_indicator
  cvd_indicator <- max(data$cvd_indicator)
  
  ### GET NA values for when we increase hazard over the max (this happens to people with long survival times,
  ### and increased risk factors. We can set these to the max followup time)
  if (is.na(t_out)){
    t_out <- max(bhaz$time)
    cvd_indicator <- 0}
  
  ### Create output tibble
  out <- data.frame("cvd_time_cf" = t_out, "cvd_indicator_cf" = cvd_indicator)
  
  return(tibble::tibble(out))
  
}

### Function for models 3,4,5,6,7
get_survtimes_adj_model34567 <- function(data, id){
  
  ### Calculate cumulative hazard at transition times
  data$cumhaz_tstart <- unlist(lapply(data$tstart, get_bhaz, bhaz = bhaz))
  data$cumhaz_tstop <- unlist(lapply(data$cvd_time, get_bhaz, bhaz = bhaz))
  
  ### Intervals which are "on treatment" get multiplied by exp(B)
  data <- dplyr::mutate(data,
                        cumhaz_tstart = cumhaz_tstart*
                          exp(offset_statins_timevar_lnHR)*
                          exp(offset_ah_timevar_lnHR),
                        cumhaz_tstop = cumhaz_tstop*
                          exp(offset_statins_timevar_lnHR)*
                          exp(offset_ah_timevar_lnHR))
  
  ### Then calculate total hazard
  cumhaz_adj <- sum(data$cumhaz_tstop - data$cumhaz_tstart)
  
  ### Now see what value of t corresponds to a hazard of this
  ### Want minimum time, for which this amount of hazard is reached
  t_out <- as.numeric(bhaz$time[min(which(bhaz$hazard >= cumhaz_adj))])
  
  ### Get cvd_indicator
  cvd_indicator <- max(data$cvd_indicator)
  
  ### GET NA values for when we increase hazard over the max (this happens to people with long survival times,
  ### and increased risk factors. We can set these to the max followup time)
  if (is.na(t_out)){
    t_out <- max(bhaz$time)
    cvd_indicator <- 0}
  
  ### Create output tibble
  out <- data.frame("cvd_time_cf" = t_out, "cvd_indicator_cf" = cvd_indicator)
  
  return(tibble::tibble(out))
  
}

######################################################
### Misc functions; used in other bigger functions ###
######################################################

###
### Function to estimate survival probability for a given model and fitted baseline hazard
### NB: Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE)
est_surv <- function(newdata, fit, bhaz, time){
  
  ### Get the lp
  lp <- predict(fit, newdata = newdata, reference = "sample")
  
  ### Calculate survival probability
  surv <- as.numeric(exp(-exp(lp)*bhaz$hazard[max(which(bhaz$time <= time))]))
  
  return(surv)
  
}

###
### Function to estimate survival time, at minimum of a variable (eventtime) or fixed numeric time (t)
### Used when estimating IP weights
est_surv_eventtime_or_t <- function(newdata, bhaz, fit, eventtime, t, type = "coxph"){
  
  if (type == "coxph"){
    
    ### Get the lp
    lp <- predict(fit, newdata = newdata, type = "lp", reference = "sample")
    
    ### Extract baseline hazard at minimum of eventtime, or t,
    # Vectorised baseline hazard lookup using findInterval for efficiency
    pat_times <- pmin(as.numeric(newdata[[eventtime]]), t)
    idx       <- pmax(findInterval(pat_times, bhaz$time), 1L)
    bhaz.pat  <- bhaz$hazard[idx]
    
    ### Calculate survival probability
    surv <- as.numeric(exp(-exp(lp) * bhaz.pat))
    
  } else if (type == "flexsurv"){
    
    ### Calculate survival probability
    surv <- unlist(lapply(1:nrow(newdata),
                          function(x) {dplyr::pull(predict(fit,
                                                           newdata = newdata[x, ],
                                                           times = min(as.numeric(newdata[x, eventtime]), t),
                                                           type = "survival"),
                                                   .pred_survival)}
    ))
    
  }
  
  return(surv)
  
}


###
### Function to estimate survival probability for a given model and fitted baseline ,
### when model has offset terms
### NB: Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE)
est_surv_offset <- function(newdata, fit, bhaz, time){
  
  ### Get the lp
  lp <- predict(fit, newdata = newdata, reference = "sample")
  
  ### Adjust for offsets, which are not deducting from the lp during predict.coxph
  ### See program: pX_verify_functions_est_surv.R
  if (grepl("offset_statins_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_statins_lnHR)"])
  }
  if (grepl("offset_ah_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_ah_lnHR)"])
  }
  if (grepl("offset_statins_timevar_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_statins_timevar_lnHR)"])
  }
  if (grepl("offset_ah_timevar_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_ah_timevar_lnHR)"])
  }
  if (grepl("offset_sbp_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_sbp_lnHR)"])
  }
  if (grepl("offset_bmi_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_bmi_lnHR)"])
  }
  if (grepl("offset_nonhdl_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_nonhdl_lnHR)"])
  }
  
  ### Calculate survival probability
  surv <- as.numeric(exp(-exp(lp)*bhaz$hazard[max(which(bhaz$time <= time))]))
  
  return(surv)
  
}

#########################################################
### Function to estimate calibration curves and plots ###
#########################################################

###
### Write a function to estimate calibration plots for a dataset using PH regression based approach
### (graphical calibration curves by Austin et al).
### fit is a fitted cox model for which calibration will be assessed 
### bhaz is corresponding baseline hazard (centered)
### Weighting this model is normally not appropriate given time-varying weights required.
###
est_calib_plot <- function(data, 
                           fit, 
                           bhaz, 
                           time, 
                           surv = NULL, 
                           nk = 5,
                           weights = NULL,
                           plot.range = NULL,
                           plot.subset = TRUE,
                           remove.low.risk = FALSE){
  

  ### Calculate survival probability
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv_offset(newdata = data, fit = fit, bhaz = bhaz, time = time))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Add predicted survival probabilities to data
  data$cloglog <- log(-log(data$surv))
  
  ### Fit calibration model
  fit.calib <- survival::coxph(survival::Surv(cvd_time, cvd_indicator) ~ rms::rcs(cloglog, nk), data = data, weights = weights)
  bhaz.calib <- survival::basehaz(fit.calib, centered = TRUE)
  
  ###
  ### Geneated predicted observed values
  
  ### First we do this over range of predicted probbilities in the validation dataset
  pred.obs <- est_surv(newdata = data, fit = fit.calib, bhaz = bhaz.calib, time = time)
  
  ### Create plot data
  plot.data <- data.frame("pred.obs" = 1 - pred.obs, "pred" = 1 - data$surv)
  
  ### Remove low risk if specified (this is purely to focus on calibration of individuals with risks over 5%, for manuscript)
  if (remove.low.risk == TRUE){
    plot.data <- subset(plot.data, pred > 0.05)
  }
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(plot.data$pred.obs - plot.data$pred))
  E50 <- median(abs(plot.data$pred.obs - plot.data$pred))
  E90 <- as.numeric(quantile(abs(plot.data$pred.obs - plot.data$pred), probs = .9))
  
  ## If plot.range is defined...
  if (!is.null(plot.range)){
    
    ## Create dataset over range of specified points, which is what the plot will be made over
    range.data <- data.frame("pred" = plot.range, "surv" = 1 - plot.range, "cloglog" = log(-log(1 - plot.range)))
    
    ## create predicted observed over this range of predicted values
    pred.obs.range <- est_surv(newdata = range.data, fit = fit.calib, bhaz = bhaz.calib, time = time)
    
    ### Create plot data
    plot.data.range <- data.frame("pred.obs" = 1 - pred.obs.range, "pred" = range.data$pred)
    
  }
  
  ### Take a random subset of 5000 values for the plot
  if (plot.subset == TRUE){
    plot.data <- plot.data[sample(1:nrow(plot.data), min(5000, nrow(data)), replace = FALSE), ] |>
      dplyr::arrange(pred)
  }
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") 
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "plotdata" = plot.data,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90)
  if (!is.null(plot.range)){
    output.object[["plotdata.range"]] <- plot.data.range
  }
  
  return(output.object)
  
}

###
### Write a function to estimate calibration plots for a dataset using Kaplan-Meier vs mean risk within subgroups.
### data is the data in which calibration will be assessed
### fit is fitted Cox-models for which calibration will be assessed
### bhaz is corresponding baseline hazard (centered)
###
est_calib_plot_group <- function(data, 
                                 fit, 
                                 bhaz, 
                                 time, 
                                 n.groups, 
                                 surv = NULL,
                                 CI = FALSE){
  
  ### Calculate survival probability
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv_offset(newdata = data, fit = fit, bhaz = bhaz, time = time))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Add predicted survival probabilities to data
  data$cloglog <- log(-log(data$surv))
  
  ### Split data by cloglog
  data.split <- split(data, cut(data$cloglog, c(-Inf, quantile(data$cloglog, probs = 1:n.groups/n.groups))))
  
  ### Create vectors to store estimates
  obs <- vector("numeric", n.groups)
  obs_upper <- vector("numeric", n.groups)
  obs_lower <- vector("numeric", n.groups)
  pred <- vector("numeric", n.groups)
  
  ### Run through
  for (i in 1:n.groups){
    
    ### Get Kaplan-Meier estimate
    survobj <- survival::survfit(survival::Surv(cvd_time, cvd_indicator) ~ 1, data = data.split[[i]])
    obs[i] <- as.numeric(survobj$surv[max(which(survobj$time <= time))])
    if (CI == TRUE){
      obs_upper[i] <- as.numeric(survobj$upper[max(which(survobj$time <= time))])
      obs_lower[i] <- as.numeric(survobj$lower[max(which(survobj$time <= time))])
    }
    
    ### Get predicted
    #     pred[i] <- 1 - as.numeric(exp(-exp(mean(data.split[[i]]$cloglog))))
    pred[i] <- 1 - (mean(data.split[[i]]$surv))
    
  }
  
  ### Create plot data
  plot.data <- data.frame("obs" = 1 - obs, "pred" = pred, "obs_lower" = 1 - obs_lower, "obs_upper" = 1 - obs_upper)
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") 
  
  if (CI == TRUE){
    plot <- plot +
      ggplot2::geom_errorbar(ggplot2::aes(pred, ymin = obs_lower, ymax = obs_upper), width = .005)
    
  }
  
  ### Create output.object
  output.object <- list("plot" = plot)
  
  return(output.object)
  
}

###
### Function to estimate calibration using proportional hazards approach 
### (Austin et al., graphical calibration curves)
###
### This function does the same thing as 'est_calib_plot', but is more flexible and organised in terms of inputs/outputs
### It has been written to implement the sensitivity analyses, in 
###
est_calib_ph <- function(data, 
                         fit, 
                         bhaz, 
                         time, 
                         surv = NULL, 
                         nk = 5,
                         weights = NULL,
                         pred.plot.range = NULL,
                         pred.cap.quantile = .99,
                         plot = TRUE,
                         remove.low.risk = FALSE){
  
  #   data = imps.valid[[1]]
  #   length(data)
  #   fit = fit.list[[1]]
  #   bhaz = bhaz.list[[1]]
  #   time = round(10*365.25)
  #   str(fit)
  #   coefficients(fit)
  #   str(data)
  
  ### Calculate survival probability
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv_offset(newdata = data, fit = fit, bhaz = bhaz, time = time))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### 2. Add predicted risk
  data$pred <- 1 - data$surv
  # Take predictions away from 0/1
  data$pred[data$pred == 1] <- 0.9999
  data$pred[data$pred == 0] <- 0.0001
  # Put onto cloglog scale
  data$cloglog <- log(-log(data$surv))
  
  ### Fit calibration model
  fit.calib <- survival::coxph(survival::Surv(time, status) ~ rms::rcs(cloglog, nk), data = data, weights = weights)
  bhaz.calib <- survival::basehaz(fit.calib, centered = TRUE)
  
  ###
  ### Geneated predicted observed values
  
  ### First we do this over range of predicted probbilities in the validation dataset
  data$pred.obs <- 1 - est_surv(newdata = data, fit = fit.calib, bhaz = bhaz.calib, time = time)
  
  ### 6. Calculate ICI, E50 and E90
  ICI <- mean(abs(data$pred.obs - data$pred))
  E50 <- median(abs(data$pred.obs - data$pred))
  E90 <- as.numeric(quantile(abs(data$pred.obs - data$pred), probs = .9, na.rm = TRUE))
  
  ### 7. Create output data
  if ("id" %in% colnames(data)){
    output.data <- data.frame("id" = data$id, "pred.obs" = data$pred.obs, "pred" = data$pred)
  } else {
    output.data <- data.frame("pred.obs" = data$pred.obs, "pred" = data$pred)
  }
  
  ### 8. Create Plotting Data
  ### Cap plot range at xth percentile of predicted risks
  pred_cap <- as.numeric(stats::quantile(data$pred, probs = pred.cap.quantile, na.rm = TRUE))
  
  if (!is.null(pred.plot.range)) {
    plot_df <- data.frame(
      pred       = pred.plot.range,
      cloglog = log(-log(1-pred.plot.range))
    )
    plot_df$pred.obs <- 
      1 - as.numeric(est_surv(newdata = plot_df, fit = fit.calib, bhaz = bhaz.calib, time = time))
  } else {
    plot_df <- data.frame(
      pred = seq(min(data$pred), pred_cap, length.out = 100)
    ) |>
      dplyr::mutate(cloglog = log(-log(1 - pred)))
    plot_df$pred.obs <- 
      1 - as.numeric(est_surv(newdata = plot_df, fit = fit.calib, bhaz = bhaz.calib, time = time))
  }
  
  ### 9. Generate Plot
  if (plot == TRUE) {
    
    ### Define axis max
    if (!is.null(pred.plot.range)){
      axis_max <- max(pred.plot.range) * 1.02
    } else {
      axis_max <- pred_cap * 1.02
    }
    
    ### Subsample raw data for the marginal histogram (max 10,000 individuals)
    data_hist <- data[sample(nrow(data), min(10000L, nrow(data))), ]
    
    ### Create plot object
    plot_obj_base <- ggplot2::ggplot() +
      ### Invisible geom_point from subsampled raw data — required by ggMarginal;
      ### drives the marginal histogram distribution
      ggplot2::geom_point(
        data    = data_hist,
        mapping = ggplot2::aes(x = pred, y = pred),
        alpha   = 0
      ) +
      ggplot2::geom_line(
        data      = plot_df,
        mapping   = ggplot2::aes(x = pred, y = pred.obs),
        color     = "red",
        linewidth = 1
      ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
      ggplot2::scale_x_continuous(limits = c(0, axis_max)) +
      ggplot2::scale_y_continuous(limits = c(0, axis_max)) +
      ggplot2::ggtitle("Proportional Hazards") +
      ggplot2::xlab("Predicted risk") +
      ggplot2::ylab("Predicted-observed risk") +
      ggplot2::theme_minimal()
    
    ### Add marginal histogram along x-axis
    plot_obj <- ggExtra::ggMarginal(
      p       = plot_obj_base,
      type    = "histogram",
      margins = "x",
      size    = 5,
      bins    = 50,
      fill    = "grey70",
      colour  = "white"
    )
    
  } else {
    plot_obj <- NULL
  }
  
  ### 10. Create and return output object
  output.object <- list("plot"       = plot_obj,
                        "ICI"        = ICI,
                        "E50"        = E50,
                        "E90"        = E90,
                        "calib_data" = output.data)
  
  ### Only return separate calibration curve data if pred.plot.range is specified
  if (!is.null(pred.plot.range)) {
    output.object$calib_data_pred_plot_range <- plot_df
  }
  
  return(output.object)
  
}

###################################################################
### Function to be used when bootstrapping confidence intervals ###
###################################################################

### NB: This function now calls on est_calib_ph, rather than est_calib_plot.
### These functions do the same thing, but est_calib_ph is superseeding est_calib_plot and will be used
### moving forwards.
boot_func_for_calib_metrics_CI <- function(data, indices, fit, bhaz, time, nk, pred.plot.range = NULL){
  
  ### Create bootstrapped data
  data_b <- data[indices, ]
  
  ### Assess calibration of model in bootstrapped dataset using existing function
  calib_object <- est_calib_ph(data = data_b, 
                               fit = fit, 
                               bhaz = bhaz, 
                               time = time, 
                               nk = nk,
                               pred.plot.range = pred.plot.range,
                               plot = FALSE)
  
  ### Extract the plot data which is of interest (predicted observed values)
  ICI <- calib_object[["ICI"]]
  E50 <- calib_object[["E50"]]
  E90 <- calib_object[["E90"]]
  
  ### Extract predicted-observed values for calibration curve
  pred_obs_range <- calib_object$calib_data_pred_plot_range$pred.obs
  
  ### output object
  output_object <- c(ICI, E50, E90, pred_obs_range)
  
  ### Return
  return(output_object)
  
}
