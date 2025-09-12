#################
### Functions ###
#################

###
### Function to estimate survival probability for a given model and fitted baseline hazard
### NB: Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE)
est_surv <- function(newdata, fit, bhaz, time){
  
  ### Get the lp
  lp <- predict(fit, newdata = newdata, reference = "sample")
  
  ### Get the linear predictor for ne wdata
  surv <- as.numeric(exp(-exp(lp)*bhaz$hazard[max(which(bhaz$time <= time))]))
  
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
  
  ### Get the linear predictor for new data
  surv <- as.numeric(exp(-exp(lp)*bhaz$hazard[max(which(bhaz$time <= time))]))
  
  return(surv)
  
}

###
### Write a function to estimate calibration plots for a dataset using PH regression based approach 
### (graphical calibration curves by Austin et al).
### Calibration is assessed in one validation dataset, but the survival probabilities can be estimated by averaging the predicted
### risks from each developed model (multiple imputation) using Rubins rules (type = 'imp' or 'single')
### data is the data in which calibration will be assessed
### fit is a fitted cox model, or list of fitted Cox-models, for which calibration will be assessed (each from a different imputed dataset)
### bhaz is corrresponding baseline hazard or list of baseline hazards (centered)
est_calib_plot <- function(data, 
                           fit, 
                           bhaz, 
                           time, 
                           surv = NULL, 
                           nk = 5,
                           plot.range = NULL,
                           plot.subset = TRUE,
                           remove.low.risk = FALSE){
  
  #   data = imps.valid[[1]]
  #   length(data)
  #   fit = fit.list[[1]]
  #   bhaz = bhaz.list[[1]]
  #   time = round(10*365.25)
  #   str(fit)
  #   coefficients(fit)
  #   str(data)
  
  # data = df_valid_pred
  # fit = fit
  # bhaz = bhaz
  # time = round(10*365.25)
  # surv = NULL
  # CI = FALSE
  # 
  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv_offset(newdata = data, fit = fit, bhaz = bhaz, time = time))
  } else {
    data$surv <- as.numeric(surv)
  }
  # max(data$surv)
  # min(data$surv)
  # head(data$surv)
  # sum(is.na(data$surv))
  # temp <- data[is.na(data$surv), ]
  # temp <- temp[1:100,]
  ### Add predicted survival probabilities to data
  data$cloglog <- log(-log(data$surv))
  
  ### Fit calibration model
  fit.calib <- survival::coxph(survival::Surv(cvd_time, cvd_indicator) ~ rms::rcs(cloglog, nk), data = data)
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
  
  ### Take a random subset of 1000 values for the plot
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
### Calibration is assessed in one validation dataset, but the survival probabilities are estimated by averaging the predicted
### risks from each developed model (multiple imputation) using Rubins rules.
### data is the data in which calibration will be assessed
### fit.list is list of fitted Cox-models for which calibration will be assessed (each from a different imputed dataset)
### bhaz.list is corrresponding list of baseline hazards (centered)
est_calib_plot_group <- function(data, 
                                 fit, 
                                 bhaz, 
                                 time, 
                                 n.groups, 
                                 surv = NULL,
                                 CI = FALSE){
  
  ### Get the survival probabilities
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