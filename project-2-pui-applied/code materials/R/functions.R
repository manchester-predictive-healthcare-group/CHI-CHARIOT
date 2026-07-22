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
### when model has offset term for statins and antihypertensives
### NB: Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE)
est_surv_offset <- function(newdata, fit, bhaz, time){
  
  ### Get the lp
  lp <- predict(fit, newdata = newdata, reference = "sample")
  
  ### Adjust for offsets
  if (grepl("med_status_statins", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(med_status_statins)"])
  }
  if (grepl("med_status_ah", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(med_status_ah)"])
  }
  
  ### Get the linear predictor for new data
  surv <- as.numeric(exp(-exp(lp)*bhaz$hazard[max(which(bhaz$time <= time))]))
  
  return(surv)
  
}

###
### Function to estimate survival probability for a given set of models and fitted baseline hazards,
### where each model was developed using a different imputed dataset.
### Survival probabilities are estimated using each model seperately, then combined on the cloglog scale
### NB: Baseline hazards must have been fitted using basehaz(surv.obj, centered = TRUE)
est_surv_mi <- function(newdata, fit.list, bhaz.list, time){
  
  #   newdata <- imps.valid[[1]][1:100, ]
  #   bhaz.list[[1]]$time
  #   time <- 10*365.25
  
  ### Get the lp
  lp <- lapply(fit.list, function(x) {predict(x, newdata = newdata, reference = "sample")})
  
  ### Get the linear predictor for ne wdata
  surv <- lapply(1:length(fit.list), function(x) {
    bhaz <- bhaz.list[[x]]
    as.numeric(exp(-exp(lp[[x]])*bhaz$hazard[max(which(bhaz$time <= time))]))
  })
  
  ### Now get average on cloglog scale
  surv.average <- exp(-exp(rowMeans(do.call("cbind", lapply(surv, function(x) {log(-log(x))})))))
  
  ### Return
  return(surv.average)
  
}


###
### Function to estimate survival probability for a given set of models and fitted baseline hazards,
### where each model was developed using a different imputed dataset.
### Survival probabilities are estimated using each model seperately, then combined on the cloglog scale
### NB: Baseline hazards must have been fitted using basehaz(surv.obj, centered = TRUE)
est_surv_offset_mi <- function(newdata, fit.list, bhaz.list, time){
  
  #   newdata <- imps.valid[[1]][1:100, ]
  #   bhaz.list[[1]]$time
  #   time <- 10*365.25
  
  ### Get the lp
  lp <- lapply(fit.list, function(x) {
    ## Get lp
    lp.out <- predict(x, newdata = newdata, reference = "sample")
    
    ## Adjust for offsets
    lp.out <- lp.out - mean(x$model[,"offset(med_status_statins)"]) - mean(x$model[,"offset(med_status_ah)"])
    return(lp.out)
  })
  
  ### Get the survival probabilities
  surv <- lapply(1:length(fit.list), function(x) {
    bhaz <- bhaz.list[[x]]
    as.numeric(exp(-exp(lp[[x]])*bhaz$hazard[max(which(bhaz$time <= time))]))
  })
  
  ### Now get average on cloglog scale
  surv.average <- exp(-exp(rowMeans(do.call("cbind", lapply(surv, function(x) {log(-log(x))})))))
  
  ### Return
  return(surv.average)
  
}

###
### Function to estimate survival probability for a given model and fitted baseline hazard, but survival probability
### is probability at the time they had their event.
### Variable eventtime is a character string, of the name of the variable denoting the event time
### NB: Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE)
est_surv_eventtime <- function(newdata, bhaz, fit, eventtime){
  
  ### Get lp
  lp <- predict(fit, newdata = newdata, type = "lp", reference = "sample")
  
  ### Get bhaz specific for each patient at the time they had their event
  bhaz.pat <- unlist(lapply(1:nrow(newdata), function(x) {bhaz$hazard[max(which(bhaz$time <= newdata[x, eventtime]))]}))
  
  ### Get expected
  surv <- as.numeric(exp(-exp(lp)*bhaz.pat))
  
  return(surv)
  
}

###
### Function to estimate expected number of events
### NB: Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE)
est_expected <- function(newdata, bhaz, fit){
  
  ### Get lp
  lp <- predict(fit, newdata = newdata, type = "lp", reference = "sample")
  
  ### Get bhaz specific for each patient at the time they had their event
  bhaz.pat <- unlist(lapply(1:nrow(newdata), function(x) {bhaz$hazard[max(which(bhaz$time <= newdata$cvd_time[x]))]}))
  
  ### Get expected
  expected <- log(as.numeric(exp(lp)*bhaz.pat))
  
  return(expected)
  
}


###
### Function to estimate expected number of events for a given set of models and fitted baseline hazards,
### where each model was developed using a different imputed dataset.
### Survival probabilities are estimated using each model seperately, then combined on the cloglog scale
### NB: Baseline hazards must have been fitted using basehaz(surv.obj, centered = TRUE)
est_expected_mi <- function(newdata, fit.list, bhaz.list){
  
  ### Get the linear predictor for ne wdata
  expected <- lapply(1:length(fit.list), function(x) {
    
    ### Get lp
    lp <- predict(fit.list[[x]], newdata = newdata, type = "lp", reference = "sample")
    
    ### Get bhaz specific for each patient at the time they had their event
    bhaz.pat <- unlist(lapply(1:nrow(newdata), function(y) {bhaz.list[[x]]$hazard[max(which(bhaz.list[[x]]$time <= newdata$cvd_time[y]))]}))
    
    ### Get expected
    expected <- log(as.numeric(exp(lp)*bhaz.pat))
    
    return(expected)
    
  })
  
  ### Now get average on cloglog scale
  expected.average <- rowMeans(do.call("cbind", expected))
  
  ### Return
  return(expected.average)
  
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
                           type = c("imp", "single"), 
                           fit.type = c("tn", "standard", "prototype3"),
                           surv = NULL, 
                           nk = 5,
                           plot.range = NULL,
                           plot.subset = TRUE){
  
  #   data = imps.valid[[1]]
  #   length(data)
  #   fit = fit.list[[1]]
  #   bhaz = bhaz.list[[1]]
  #   time = round(10*365.25)
  #   str(fit)
  #   coefficients(fit)
  #   str(data)
  
  # data = df_valid_temporal
  # time = round(t_eval*365.25)
  # fit.type = "prototype3"
  # surv = df_valid_temporal$surv_intervention
  # type = "single"
  
  ### Set args
  type <- match.arg(type)
  fit.type <- match.arg(fit.type)
  
  ### Get the survival probabilities
  if (is.null(surv)){
    if (fit.type == "tn"){
      if (type == "imp"){
        data$surv <- as.numeric(est_surv_offset_mi(newdata = data, fit.list = fit, bhaz.list = bhaz, time = time))
      } else if (type == "single"){
        data$surv <- as.numeric(est_surv_offset(newdata = data, fit = fit, bhaz = bhaz, time = time))
      }
    } else if (fit.type == "standard"){
      if (type == "imp"){
        data$surv <- as.numeric(est_surv_mi(newdata = data, fit.list = fit, bhaz.list = bhaz, time = time))
      } else if (type == "single"){
        data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, time = time))
      }
    } else if (fit.type == "prototype3"){
      if (type == "imp"){
        stop("CODE NOT WRITTEN FOR THIS YET")
      } else if (type == "single"){
        data$surv <- as.numeric(est_surv_offset_prototype3(newdata = data, fit = fit, bhaz = bhaz, time = time))
      }
    }
  } else {
    data$surv <- as.numeric(surv)
  }
  
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
### Function to estimate calibration using pseudo-value appraoch, where pseudo-values estimated using
### the kaplan-meier estimate of survival
###

### First write a function which will estimate pseudo-values for a given dataset at a time t
est_pv <- function(df, t){
  
  ### Create prodlim object
  prodlim.obj <- prodlim::prodlim(prodlim::Hist(cvd_time, cvd_indicator) ~ 1, data = df)
  
  ### Calculate psuedo-values
  pv <- prodlim::jackknife(prodlim.obj, times = t)
  
  return(as.numeric(pv))
}

###
### Write a function to estimate calibration plots for a dataset using pseudo-values 
### Calibration is assessed in one validation dataset, but the survival probabilities can be estimated by averaging the predicted
### risks from each developed model (multiple imputation) using Rubins rules (type = 'imp' or 'single')
### data is the data in which calibration will be assessed
### fit is a fitted cox model, or list of fitted Cox-models, for which calibration will be assessed (each from a different imputed dataset)
### bhaz is corrresponding baseline hazard or list of baseline hazards (centered)
est_calib_plot_pv <- function(data, 
                              fit, 
                              bhaz, 
                              time, 
                              type = c("imp", "single"), 
                              fit.type = c("tn", "standard"),
                              surv = NULL, 
                              plot.range = NULL,
                              plot.subset = TRUE,
                              nk = 5, 
                              split.n.groups = 1, 
                              split.var = NULL){
  
  # data = data.valid.off[1:10000, ]
  # fit = fit
  # bhaz = bhaz
  # time = round(10*365.25)
  # type = "single"
  # fit.type = "tn"
  # plot.range = NULL
  # nk <- 5
  # split.n.groups = 1
  # surv = NULL
  
  ### Set args
  type <- match.arg(type)
  fit.type <- match.arg(fit.type)
  
  ### Get the survival probabilities
  if (is.null(surv)){
    if (fit.type == "tn"){
      if (type == "imp"){
        data$surv <- as.numeric(est_surv_offset_mi(newdata = data, fit.list = fit, bhaz.list = bhaz, time = time))
      } else if (type == "single"){
        data$surv <- as.numeric(est_surv_offset(newdata = data, fit = fit, bhaz = bhaz, time = time))
      }
    } else if (fit.type == "standard"){
      if (type == "imp"){
        data$surv <- as.numeric(est_surv_mi(newdata = data, fit.list = fit, bhaz.list = bhaz, time = time))
      } else if (type == "single"){
        data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, time = time))
      }
    }
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Add predicted risk
  data$pred <- 1 - data$surv
  
  ### Split data by surv
  data.pv.split <- split(data, 
                         cut(data$surv, 
                             c(-Inf, quantile(data$surv, probs = 1:split.n.groups/split.n.groups))))
  
  ### Calculate pseudo-values for each data split (note est_pv calculates pseudo-value for survival prob, so want to take 1 - pv)
  pv <- lapply(data.pv.split, est_pv, t = time)
  pvs <- 1 - do.call("c", pv)
  ### Also get the corresponding ids
  ids <- lapply(data.pv.split, function(x) {x$patid})
  ids <- do.call("c", ids)
  
  ### Combine
  pvs <- data.frame("pv" = pvs, "patid" = ids)
  pvs <- dplyr::arrange(pvs, patid)
  
  ### Add the pseudo-values to validation dataset
  data$pv <- pvs$pv
  
  ### Transform est.surv onto logit scale
  data$pred.logit <- log(data$pred/(1-data$pred))
  
  ### Fit the model using logit link function
  calib.model.pv <- stats::glm(pv ~ rms::rcs(pred.logit, nk), 
                               data = data, 
                               family = stats::gaussian(link = "logit"),
                               start = rep(0, nk))
  
  ### Estimate predict observed values
  pred.obs <- predict(calib.model.pv, newdata = data, type = "response")
  
  ### Create plot data
  plot.data <- data.frame("pred.obs" = pred.obs, 
                          "pred" = data$pred)
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle(("Pseudo-value")) + 
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(plot.data$pred.obs - plot.data$pred))
  E50 <- median(abs(plot.data$pred.obs - plot.data$pred))
  E90 <- as.numeric(quantile(abs(plot.data$pred.obs - plot.data$pred), probs = .9, na.rm = TRUE))
  
  ## If plot.range is defined...
  if (!is.null(plot.range)){
    
    ## Create dataset over range of specified points, which is what the plot will be made over
    range.data <- data.frame("pred" = 1 - plot.range, 
                             "pred.logit" = log(plot.range/(1-plot.range)))
    
    ## create predicted observed over this range of predicted values
    pred.obs.range <- predict(calib.model.pv, newdata = range.data, type = "response")
    
    ### Create plot data
    plot.data.range <- data.frame("pred.obs" = pred.obs.range, "pred" = range.data$pred)
    
  }
  
  ### Take a random subset of 1000 values for the plot
  if (plot.subset == TRUE){
    plot.data <- plot.data[sample(1:nrow(plot.data), 1000, replace = FALSE), ] |>
      dplyr::arrange(pred)
  }
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "plotdata" = plot.data,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90)
  if (exists("plot.data.range")){
    output.object[["plotdata.range"]] <- plot.data.range
  }
  
  return(output.object)
  
}


###
### Write a function to estimate calibration plots for a dataset using IPCWs 
### Calibration is assessed in one validation dataset, but the survival probabilities can be estimated by averaging the predicted
### risks from each developed model (multiple imputation) using Rubins rules (type = 'imp' or 'single')
### data is the data in which calibration will be assessed
### fit is a fitted cox model, or list of fitted Cox-models, for which calibration will be assessed (each from a different imputed dataset)
### bhaz is corrresponding baseline hazard or list of baseline hazards (centered)
est_calib_plot_IPCW <- function(data, fit, bhaz, time, type = c("imp", "single"), surv = NULL, max.follow = NULL, calib.type, cens.formula){
  
  #   data = imps.valid[[1]][1:500000, ]
  #   fit = fit.list[[1]]
  #   bhaz = bhaz.list[[1]]
  #   time = round(10*365.25)
  #   type = "single"
  #   calib.type = "rcs"
  #   cens.formula = cens.formula
  #   max.follow = NULL
  #   surv = NULL
  
  ### Defeine data.valid
  data.valid <- data
  
  ### Get the survival probabilities
  if (is.null(surv)){
    if (type == "imp"){
      data.valid$surv <- as.numeric(est_surv_mi(newdata = data.valid, fit.list = fit, bhaz.list = bhaz, time = time))
    } else if (type == "single"){
      data.valid$surv <- as.numeric(est_surv(newdata = data.valid, fit = fit, bhaz = bhaz, time = time))
    }
  } else {
    data.valid$surv <- as.numeric(surv)
  }
  
  ### Assign a variable for whether individuals have had an event by time
  data.valid <- dplyr::mutate(data.valid, 
                              cvd_indicator_time = dplyr::case_when(cvd_time <= time ~ cvd_indicator,
                                                                    cvd_time > time ~ 0))
  
  ### Lets create a variable for censoring or not censored
  data.valid$dtcens <- data.valid$cvd_time
  data.valid$dtcens.s <- 1 - data.valid$cvd_indicator
  
  ### Don't want to model the effect of everybody becoming censored at the final follow up time, 
  ### so set these individuals to be uncensored at time point of interst
  if (!is.null(max.follow)){
    data.valid <- dplyr::mutate(data.valid, 
                                dtcens.s = dplyr::case_when(dtcens < max.follow + 2 ~ dtcens.s, 
                                                            dtcens >= max.follow + 2 ~ 0), 
                                dtcens = dplyr::case_when(dtcens < max.follow + 2 ~ dtcens, 
                                                          dtcens >= max.follow + 2 ~ max.follow + 2))
  }
  
  ### Create censoring model
  cens.model <- survival::coxph(cens.formula, data = data.valid, model = TRUE)
  cens.model.int <- survival::coxph(survival::Surv(dtcens, dtcens.s) ~ 1, data = data.valid)
  
  ### Create survival predictions from this model, probability of being uncensored at time t
  cens.bhaz <- survival::basehaz(cens.model, centered = TRUE)
  data.valid$surv.cens <- est_surv(newdata = data.valid, fit = cens.model, bhaz = cens.bhaz, time = time)
  
  ### For individual who have an event prior to time t (and are therefore uncesored at time t), 
  ### we want probability of being uncensored at t_event. For people cenosred before time t,
  ### they will not be included in the analysis
  obs.event.prior <- which(data.valid$dtcens <= time & data.valid$dtcens.s == 0)
  obs.censored.prior <- which(data.valid$dtcens <= time & data.valid$dtcens.s == 1)
  
  ### Calculate weights for individuals who have events prior to time t
  surv.event.prior <- est_surv_eventtime(newdata = data.valid[obs.event.prior, ], bhaz = cens.bhaz, fit = cens.model, eventtime = "dtcens")
  
  ### Replace these values in dataset
  data.valid$surv.cens[obs.event.prior] <- surv.event.prior
  
  ### Create weights
  data.valid$ipcw <- 1/data.valid$surv.cens
  
  ### Stabilise
  cens.bhaz.int <- survival::basehaz(cens.model.int, centered = TRUE)
  data.valid$ipcw.numer <- est_surv(newdata = data.valid, fit = cens.model.int, bhaz = cens.bhaz.int, time = time)
  data.valid$ipcw.stab <- data.valid$ipcw.numer*data.valid$ipcw
  
  ### Assign NA's to individuals censored prior to time t
  data.valid$ipcw[obs.censored.prior] <- NA
  data.valid$ipcw.stab[obs.censored.prior] <- NA
  
  ### Reduce data.valid to individuals uncesored at time t
  data.valid <- subset(data.valid, !is.na(ipcw))
  
  ### Convert predict risk onto logit scale
  data.valid$pred <- 1 - data.valid$surv
  data.valid$pred.logit <- log(data.valid$pred/(1 - data.valid$pred))
  
  ### Fit weighted calibration model and generate predicted observed
  if (calib.type == "rcs"){
    ## Fit model
    rcs.model.stab <- rms::lrm(cvd_indicator_time ~ rms::rcs(pred.logit, 5), 
                               data = data.valid, 
                               weights = data.valid[, "ipcw.stab"])
    ## Generate predicted
    data.valid$pred.obs <- predict(rcs.model.stab, newdata = data.valid, type = "fitted")
    
  } else if (calib.type == "loess"){
    ## Fit model
    loess.model <- stats::loess(cvd_indicator_time ~ pred, 
                                data = data.valid, 
                                weights = data.valid[, "ipcw.stab"])
    
    ## Generate predicted
    data.valid$pred.obs <- predict(loess.model, newdata = data.valid)
  }
  
  ### Take a random subset of 1000 values for the plot 
  ### (I should really be producing observed values over a fixed range of predicted risks, or over
  ### a range of predicted risks, equally spaced between max and min for the validation dataset)
  ### Leave for now, and until checked Angela wood paper on this topic, which references another paper.
  plot.data.small <- data.valid[sample(1:nrow(data.valid), 1000, replace = FALSE), ] |>
    dplyr::arrange(pred)
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data.small) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") 
  plot
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(data.valid$pred.obs - data.valid$pred))
  E50 <- median(abs(data.valid$pred.obs - data.valid$pred))
  E90 <- as.numeric(quantile(abs(data.valid$pred.obs - data.valid$pred), probs = .9))
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90)
  
  ### Return
  return(output.object)
  
}

###
### Write a function to estimate calibration plots using predRupdate
### NOTE THIS WILL ONLY WORK FOR A SPECIFIC MODEL, BECAUSE THE COEFFICIENTS FROM THE FIT ARE USED AS ARGUMENTS FOR PREDRUPDATE
### i.e. functionality of tyep = c("imp", "single") is not there
### Functionality for "surv" is also not there, as predicted risks are etimate from the model as part of predRupdate
est_calib_plot_predRupdate <- function(data, fit, time, type = c("imp", "single"), surv = NULL, formula){
  
  #   data = imps.valid[[1]]
  #   length(data)
  #   fit = fit.list[[1]]
  #   bhaz = bhaz.list[[1]]
  #   time = round(10*365.25)
  #   str(fit)
  #   coefficients(fit)
  #   str(data)
  
  ### Assign bhaz
  bhaz <- basehaz(fit.list[[1]], centered = FALSE)
  
  ### Get predictor matirx
  data.valid <- model.matrix(formula, data = data)
  data.valid <- data.valid[,-1]
  data.valid <- data.frame(data.valid)
  colnames(data.valid) <- paste("var", 1:length(coefficients(fit)), sep = "")
  
  ### Add outcome data to data.valid
  data.valid$cvd_time <- imps.valid[[1]][, "cvd_time"]
  data.valid$cvd_indicator <- imps.valid[[1]][, "cvd_indicator"]
  
  ### Create coefs_Table
  coefs_table <- data.frame(matrix(NA, nrow = 1, ncol = length(coefficients(fit))))
  colnames(coefs_table) <- paste("var", 1:length(coefficients(fit)), sep = "")
  coefs_table[1, ] <- as.numeric(coefficients(fit))
  
  ### Create cumhaz table
  bhaz_table <- data.frame("time" = time, "hazard" = bhaz$hazard[max(which(bhaz$time <= time))])
  
  ### Pass into pred_input_info()
  existing_tte_model <- predRupdate::pred_input_info(model_type = "survival",
                                                     model_info = coefs_table,
                                                     cum_hazard = bhaz_table)
  
  ### Pass into pred_validate
  validation_results <- predRupdate::pred_validate(x = existing_tte_model,
                                                   new_data = data.valid,
                                                   survival_time = "cvd_time",
                                                   event_indicator = "cvd_indicator",
                                                   time_horizon = time)
  
  plot <- validation_results$flex_calibrationplot
  
  ### Create output.object
  output.object <- list("plot" = plot)
  
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
                                 type = c("imp", "single"), 
                                 fit.type = c("tn", "standard","prototype3"),
                                 surv = NULL,
                                 single.CI = FALSE){
  
  #     data = imps.valid[[1]][1:100000, ]
  #     time = round(10*365.25)
  #     n.groups <- 10
  
  ### Set args
  type <- match.arg(type)
  fit.type <- match.arg(fit.type)
  
  ### Get the survival probabilities
  if (is.null(surv)){
    if (fit.type == "tn"){
      if (type == "imp"){
        data$surv <- as.numeric(est_surv_offset_mi(newdata = data, fit.list = fit, bhaz.list = bhaz, time = time))
      } else if (type == "single"){
        data$surv <- as.numeric(est_surv_offset(newdata = data, fit = fit, bhaz = bhaz, time = time))
      }
    } else if (fit.type == "standard"){
      if (type == "imp"){
        data$surv <- as.numeric(est_surv_mi(newdata = data, fit.list = fit, bhaz.list = bhaz, time = time))
      } else if (type == "single"){
        data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, time = time))
      }
    } else if (fit.type == "prototype3"){
      if (type == "imp"){
        stop("CODE NOT WRITTEN FOR THIS YET")
      } else if (type == "single"){
        data$surv <- as.numeric(est_surv_offset_prototype3(newdata = data, fit = fit, bhaz = bhaz, time = time))
      }
    }
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
    if (single.CI == TRUE){
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
  
  if (single.CI == TRUE){
    plot <- plot +
      ggplot2::geom_errorbar(ggplot2::aes(pred, ymin = obs_lower, ymax = obs_upper), width = .005)
    
  }
  
  ### Create output.object
  output.object <- list("plot" = plot)
  
  return(output.object)
  
}


###
### Write a function to estimate calibration plots for a dataset using Poisson approach (Crowson et al)
### Calibration is assessed in one validation dataset, but the survival probabilities are estimated by averaging the predicted
### risks from each developed model (multiple imputation) using Rubins rules.
### data is the data in which calibration will be assessed
### fit.list is list of fitted Cox-models for which calibration will be assessed (each from a different imputed dataset)
### bhaz.list is corrresponding list of baseline hazards (centered)
est_calib_plot_poisson <- function(data, fit, bhaz, time, type = c("imp", "single"), surv = NULL){
  
  #     data = imps.valid[[1]]
  #     time = round(10*365.25)
  #     str(fit)
  #     coefficients(fit)
  #     str(data)
  
  ### Get the survival probabilities
  if (is.null(surv)){
    if (type == "imp"){
      data$surv <- as.numeric(est_surv_mi(newdata = data, fit.list = fit, bhaz.list = bhaz, time = time))
    } else if (type == "single"){
      data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, time = time))
    }
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Add predicted survival probabilities to data
  data$cloglog <- log(-log(data$surv))
  data$fup_py <- log(data$cvd_time)
  
  ### Fit the calibration model
  fit.calib <- glm(cvd_indicator ~ rms::rcs(cloglog, 5) + offset(fup_py), family = poisson, data = data)
  
  ### Create dummy dataset for predictions with fixed follow up
  data.dummy <- data
  data.dummy$fup_py <- rep(log(time), nrow(data.dummy))
  data.dummy$pred.obs <- predict(fit.calib, newdata = data.dummy, type = "link")
  data.dummy$pred.obs <- exp(-exp(data.dummy$pred.obs))
  
  ### Create plot data
  plot.data <- data.frame("pred.obs" = 1 - data.dummy$pred.obs, "pred" = 1 - data.dummy$surv)
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(plot.data$pred.obs - plot.data$pred))
  E50 <- median(abs(plot.data$pred.obs - plot.data$pred))
  E90 <- as.numeric(quantile(abs(plot.data$pred.obs - plot.data$pred), probs = .9))
  
  ### Take a random subset of 2000 values for the plot 
  ### (I should really be producing observed values over a fixed range of predicted risks, or over
  ### a range of predicted risks, equally spaced between max and min for the validation dataset)
  ### Leave for now, and until checked Angela wood paper on this topic, which references another paper.
  plot.data.small <- plot.data[sample(1:nrow(plot.data), 2000, replace = FALSE), ] |>
    dplyr::arrange(pred)
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data.small) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") 
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90)
  
  return(output.object)
  
}


###
### Write a function to estimate calibration plots (grouped) for a dataset using Poisson approach (Crowson et al)
### Calibration is assessed in one validation dataset, but the expected probabilities are estimated by averaging the predicted
### risks from each developed model (multiple imputation) using Rubins rules.
### data is the data in which calibration will be assessed
### fit.list is list of fitted Cox-models for which calibration will be assessed (each from a different imputed dataset)
### bhaz.list is corrresponding list of baseline hazards (centered)
### NB: XXXX Need to get est_expected working properly before we can use this approach
est_calib_plot_poisson_group <- function(data, fit, bhaz, time, n.groups, type = c("imp", "single")){
  
  #   data <- imps.valid[[1]][1:10000, ]
  #   fit <- fit.list[[1]]
  #   bhaz <- bhaz.list[[1]]
  #   time <- 10*365.25
  #   n.groups <- 10
  #   type <- "single"
  
  ### Get the survival probabilities
  if (type == "imp"){
    surv <- est_surv_mi(newdata = data, fit.list = fit, bhaz.list = bhaz, time = time)
  } else if (type == "single"){
    surv <- est_surv(newdata = data, fit = fit, bhaz = bhaz, time = time)
  }
  
  ### Add predicted survival probabilities to data
  data$surv <- as.numeric(surv)
  data$cloglog <- log(-log(data$surv))
  
  ### Get the expectd number of events
  if (type == "imp"){
    data$p <- est_expected_mi(newdata = data, fit.list = fit, bhaz.list = bhaz)
  } else if (type == "single"){
    data$p <- est_expected(newdata = data, fit = fit, bhaz = bhaz)
  }
  
  ### When creating group, combine categories that are small
  data$group <- cut(data$cloglog,c(-Inf, quantile(data$cloglog, (1:(n.groups - 1))/n.groups), Inf))
  
  ### Fit calibration model
  fit.calib <- glm(cvd_indicator ~ -1 + group + offset(p), family = poisson, data = data)
  
  ###
  ### Produce data for plot
  ###
  
  ### Dummy dataset
  dummy.data <- data.frame(matrix(vector(), nrow = n.groups))
  
  ### Add mean risk for each group
  dummy.data$mean.risk <- unlist(lapply(names(table(data$group)), function(x) {
    mean(1-data$surv[data$group == x])}))
  
  ### Create expected number of events, from mean risk
  dummy.data$p <- log(dummy.data$mean.risk)
  dummy.data$group <- names(table(data$group))
  
  ### Create predicted-observed values from the model
  dummy.data$pred.obs <- predict(fit.calib, newdata = dummy.data, type = "response")
  
  ### Create plot data
  plot.data <- data.frame("pred.obs" = dummy.data$pred.obs, "pred" = dummy.data$mean.risk)
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_point(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") 
  
  ###
  ### Produce data for observed/expected
  ###
  
  ### Get coefs
  coefs <- exp(coefficients(fit.calib))
  names(coefs) <- names(table(data$group))
  
  ### Observed/expected
  OE <- unlist(lapply(names(table(data$group)), function(x) {sum(data$cvd_indicator[data$group == x])/sum(exp(data$p[data$group == x]))}))
  names(OE) <- names(table(data$group))
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "coefs" = coefs,
                        "OE" = OE)
  
  return(output.object)
  
}



# est_calib_plot_poisson_groupv2 <- function(data, fit.list, bhaz.list, time, n.groups){
#   
#   ### Get the survival probabilities
#   surv <- est_surv_mi(newdata = data, fit.list = fit.list, bhaz.list = bhaz.list, time = time)
#   
#   ### Add predicted survival probabilities to data
#   data$surv <- as.numeric(surv)
#   data$cloglog <- log(-log(data$surv))
#   
#   ### Split data by cloglog
#   data.split <- split(data, cut(data$cloglog, c(-Inf, quantile(data$cloglog, probs = 1:n.groups/n.groups))))
#   
#   ### Create vectors to store estimates
#   obs <- vector("numeric", n.groups)
#   pred <- vector("numeric", n.groups)
#   
#   ### Run through
#   for (i in 1:n.groups){
#     ### Get Kaplan-Meier estimate
#     total.fup <- sum(data.split[[i]]$cvd_time)
#     total.ind <- sum(data.split[[i]]$cvd_indicator)
#     rate <- total.ind/total.fup
#     obs[i] <- exp(-rate*10*365.25)
# 
#     ### Get predicted
#     pred[i] <- 1 - as.numeric(exp(-exp(mean(data.split[[i]]$cloglog))))
#     
#   }
#   
#   ### Create plot data
#   plot.data <- data.frame("obs" = 1 - obs, "pred" = pred)
#   
#   ### Create plot
#   plot <- ggplot2::ggplot(data = plot.data) +
#     ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), color = "red") +
#     ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") 
#   
#   ### Create output.object
#   output.object <- list("plot" = plot)
#   
#   return(output.object)
#   
# }


###
### Alternative approaches
###
est_calib_plot_poisson_B <- function(data, fit, bhaz, time, type = c("imp", "single")){
  
  #       data = imps.valid[[1]]
  #       time = round(10*365.25)
  #       fit <- fit.list[[1]]
  #       bhaz <- bhaz.list[[1]]
  #       str(fit)
  #       coefficients(fit)
  #       str(data)
  #       type <- "single"
  
  ### Get the survival probabilities
  if (type == "imp"){
    surv <- est_surv_mi(newdata = data, fit.list = fit, bhaz.list = bhaz, time = time)
  } else if (type == "single"){
    surv <- est_surv(newdata = data, fit = fit, bhaz = bhaz, time = time)
  }
  
  ### Get other required quantities
  data$surv <- as.numeric(surv)
  data$cloglog <- log(-log(data$surv))
  data$fup_py <- log(data$cvd_time)
  
  ### Fit the calibration model
  fit.calib <- glm(cvd_indicator ~ rms::rcs(cloglog, 5) + offset(fup_py), family = poisson, data = data)
  
  ### Create dummy dataset for predictions with fixed follow up
  data.dummy <- data
  data.dummy$fup_py <- rep(log(time), nrow(data.dummy))
  data.dummy$pred.obs <- predict(fit.calib, newdata = data.dummy, type = "response")
  #data.dummy$pred.obs <- exp(-data.dummy$pred.obs)
  head(data.dummy)
  
  ### Get the predicted number of expected events
  data.dummy$pred.exp <- -log(data.dummy$surv)
  #   data.dummy2 <- data
  #   data.dummy2$fup_py <- rep(log(1000), nrow(data.dummy2))
  #   data.dummy2$pred.obs <- predict(fit.calib, newdata = data.dummy2, type = "link")
  #   data.dummy2$exppred.obs <- exp(-data.dummy2$pred.obs)
  #   head(data.dummy2)
  ### Create plot data
  #   plot.data <- data.frame("pred.obs" = 1 - data.dummy$pred.obs, "pred" = 1 - data.dummy$surv)
  plot.data <- data.frame("pred.obs" = data.dummy$pred.obs, "pred" = data.dummy$pred.exp)
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(plot.data$pred.obs - plot.data$pred))
  E50 <- median(abs(plot.data$pred.obs - plot.data$pred))
  E90 <- as.numeric(quantile(abs(plot.data$pred.obs - plot.data$pred), probs = .9))
  
  ### Take a random subset of 2000 values for the plot 
  ### (I should really be producing observed values over a fixed range of predicted risks, or over
  ### a range of predicted risks, equally spaced between max and min for the validation dataset)
  ### Leave for now, and until checked Angela wood paper on this topic, which references another paper.
  plot.data.small <- plot.data[sample(1:nrow(plot.data), 2000, replace = FALSE), ] |>
    dplyr::arrange(pred)
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data.small) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") 
  plot
  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90)
  
  return(output.object)
  
}

###
### Function to take average of test observations on the same day
### Used in conjunction with group_modify.
### Data must be grouped by patid, obsdate contains date of observation, value contains value of observation
###
take_average <- function(data, id){
  
  ### Sort by issuedate
  data <- dplyr::arrange(data, obsdate)
  
  ### Get means
  means <- group_by(data, obsdate) |> summarise(value = mean(value)) |> pull(value)
  
  ### Remove duplicates
  data <- data[!duplicated(data$obsdate),]
  
  ### Add means
  data$value <- means
  
  return(data)
  
}

###
### Functions for summarising data in p4/p2.3
###

### Get individual trajectories
get_plot_trajectory <- function(df, title, step = FALSE){
  
  ### Get first 10 individuals
  patids_temp <- df$patid[!duplicated(df$patid)] |> magrittr::extract(1:10)
  
  ### Make plot
  if (step == FALSE){
    plot_trajectory <- ggplot(data = subset(df, patid %in% patids_temp)) + 
      geom_line(aes(x = obsdate_fup, y = value, group = patid, color = patid)) +
      theme(legend.position = "none") +
      ggtitle(title) + xlab("Years")
  } else {
    plot_trajectory <- ggplot(data = subset(df, patid %in% patids_temp)) + 
      geom_step(aes(x = obsdate_fup, y = value, group = patid, color = patid)) +
      theme(legend.position = "none") +
      ggtitle(title) + xlab("Years")
  }
  
  
  return(plot_trajectory)
  
}

### Get a smooth average of trajectories
get_plot_smooth_trajectory <- function(df, title, ah = NULL, statins = NULL){
  
  ### Filter df if ah is not NULL
  if (!is.null(ah)){
    df <- subset(df, ah_record_post == ah)
  }
  
  ### Filter df if statins is not NULL
  if (!is.null(statins)){
    df <- subset(df, statins_record_post == statins)
  }
  
  ### Create a smoothed curve for one of these
  model_slope <- lm(value ~ obsdate_fup, data = df)
  slope <- as.numeric(model_slope$coefficients["obsdate_fup"])
  
  ### Create a smoothed curve for one of these
  model_smooth <- lm(value ~ rcs(obsdate_fup, 5), data = df)
  
  ### Create data frame for plot
  df_plot_smooth <- df[sample(1:nrow(df), 2000, replace = FALSE), ]
  df_plot_smooth$pred_smooth <- predict(model_smooth, newdata = df_plot_smooth)
  
  ### Create a plot by group
  plot_smooth <- ggplot(data = arrange(df_plot_smooth, obsdate_fup)) +
    geom_point(aes(x = obsdate_fup, y = value)) + 
    geom_line(aes(x = obsdate_fup, y = pred_smooth), color = "red") +
    ggtitle(title) + xlab("Years")
  
  ### Create output object
  output_object <- list("plot_smooth" = plot_smooth, "slope" = slope)
  
  return(output_object)
  
}



#################################
### FUNCTIONS FOR PROTOTYPE 3 ###
#################################

###
### Function to estimate survival probability for a given model and fitted baseline ,
### when model has offset term for statins and antihypertensives
### NB: Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE)
est_surv_offset_prototype3 <- function(newdata, fit, bhaz, time){
  
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
  if (grepl("offset_smoking_dummy1_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_smoking_dummy1_lnHR)"])
  }
  if (grepl("offset_smoking_dummy2_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_smoking_dummy2_lnHR)"])
  }
  if (grepl("offset_smoking_timevar_dummy1_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_smoking_timevar_dummy1_lnHR)"])
  }
  if (grepl("offset_smoking_timevar_dummy2_lnHR", paste(as.character(fit$formula), collapse = " "))){
    lp <- lp - mean(fit$model[,"offset(offset_smoking_timevar_dummy2_lnHR)"])
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
### Function to estimate survival probability for a given set of models and fitted baseline hazards,
### where each model was developed using a different imputed dataset.
### Survival probabilities are estimated using each model seperately, then combined on the cloglog scale
### 
### NB: Baseline hazards must have been fitted using basehaz(surv.obj, centered = TRUE)
###
### NB: Unlike the other est_surv functions, the offsets are inputted manually, rather than estimating them
### from the fitted object. This is because to estimate them from the fitted object, the development data frame
### must be provided in the fit object, which means we could not share it for others to use. Instead, we provide
### the means off the offsets to be deducted off the linear predictor seperately.
###
### Prototype 3 model
est_surv_offset_prototype3_mi_manual_offset <- function(newdata, fit_list, bhaz_list, 
                                          means_statins, means_ah, means_smoking1, means_smoking2, 
                                          time){
  
  #   newdata <- imps.valid[[1]][1:100, ]
  #   bhaz_list[[1]]$time
  #   time <- 10*365.25
  
  # newdata = pat.prototype3.baseline
  # fit_list = fit.prototype3
  # bhaz_list = bhaz.prototype3
  # means_statins = means_statins
  # means_ah = means_ah
  # means_smoking1 = means_smoking1
  # means_smoking2 = means_smoking2
  
  ### Get the lp
  lp <- lapply(1:length(fit_list), function(x) {
    
    ## Get lp
    lp.out <- predict(fit_list[[x]], newdata = newdata, reference = "sample")
    
    ## Adjust for offsets
    lp.out <- lp.out - means_statins[[x]] - means_ah[[x]] - means_smoking1[[x]] - means_smoking2[[x]]
    
    return(lp.out)
  })
  
  ### Get the survival probabilities
  surv <- lapply(1:length(fit_list), function(x) {
    bhaz <- bhaz_list[[x]]
    as.numeric(exp(-exp(lp[[x]])*bhaz))
  })
  
  ### Now get average on cloglog scale
  surv_average <- exp(-exp(rowMeans(do.call("cbind", lapply(surv, function(x) {log(-log(x))})))))
  
  ### Return
  return(surv_average)
  
}


###
### Create dummy variables for smoking (at baseline, and during follow-up)
### This function is utilised in programs p4/p5.1 and p4/p6.1 (both prototype 3)
### The 'data' argument must contain "smoking" and "med_status_adj_smoking" variables
###
### These variables come from merging a baseline cohort (from one of the imputed datasets), 
### with cohort.split.times, which is the interval censored outcome data with time-varying variables
### for ah/statins/smoking.
###
### smoking is the smoking status at baseline variable (takes values non (never), ex, light, moderate, heavy)
### med_status_adj_smoking is the time-varying smoking status variable, takes values 1 (never), 2 (ex) or 3 (current)
###
### Note that med_status_adj_smoking is missing at baseline, because when deriving changes in smoking status during
### follow-up, we did not know value at baseline for some individuals. This comes from the imputed data.
###
### This program fills in the missing values, and creates dummies for smoking at baseline, and smoking during follow-up.
###
### At baseline: dummy1 = 1 if ex or current smoker, 0 if never smoker. 
### At baseline: dummy2 = 1 if ex-smoker, 0 otherwise.
### During follow-up, dummy1 = 0 if no change from baseline, = 1 if individual initates smoking
### During follow-up, dummy2 = -1 if an ex-smoker (at baseline) starts smoking, = 1 if a smoker (at baseline) stops smoking
###
### More details on the definition of the dummies is given inline
create_smoking_dummies <- function(data){
  
  ###
  ### Creation of baseline smoking variables
  ###
  
  ### The first dummy variable = 0 if never smoker, = 1 if current or past smoker.
  ### NB: This dummy effect (HR1) is the effect of current smoker vs never smoker
  ###
  ### The second dummy variable = 0 if never or current smoker, = 1 if past smoker
  ### NB: This dummy effect (HR2) is the effect of smoking cessation (stopping smoking, past smoker vs current smoker)
  ###
  ### These dummies are only appropriate when applied in tandem.
  ### A never smoker will have: LP*1*1
  ### A current smoker will have: LP*HR1*1
  ### A past smoker will have: LP*HR1*HR2, i.e., the effect of initiating smoking, then the effect of stopping
  ###
  data <-
    dplyr::mutate(data,
                  smoking_baseline = dplyr::case_when(smoking == "Non-smoker" ~ 0,
                                                      smoking == "Ex-smoker" ~ 1,
                                                      !(smoking %in% c("Non-smoker", "Ex-smoker")) ~ 2),
                  smoking_baseline_dummy1 = dplyr::case_when(smoking_baseline == 0 ~ 0,
                                                             TRUE ~ 1),
                  smoking_baseline_dummy2 = dplyr::case_when(smoking_baseline == 1 ~ 1,
                                                             TRUE ~ 0))
  
  
  ###
  ### Creation of change in smoking variables during follow-up (TIME VARYING)
  ###
  
  ### med_status_adj_smoking is a time-varying variable for smoking status during follow-up, however some are missing
  ### at baseline. These have been imputed, and the imputed values are in the "smoking" variable. They are missing from
  ### med_status_adj_smoking because we calculated the time-varying variable for all time points during follow-up (and not baseline), for 
  ### computational reasons.
  ###
  ### Add in imputed (or observed) smoking values at baseline to the time varying smoking variable
  ### If med_status_adj_smoking = NA, this is the first value of smoking status, for which we didn't know when
  ### create the cohort split times (as it was imputed for a large number of individual, and varies between imputed datasets)
  ### Therefore take the value of smoking, when med_status_adj_smoking is NA
  data <- 
    dplyr::mutate(data,
                  med_status_adj_smoking = 
                    dplyr::case_when(is.na(med_status_adj_smoking) & smoking == "Non-smoker" ~ 0,
                                     is.na(med_status_adj_smoking) & smoking == "Ex-smoker" ~ 1,
                                     is.na(med_status_adj_smoking) & !(smoking %in% c("Non-smoker", "Ex-smoker")) ~ 2,
                                     TRUE ~ med_status_adj_smoking))
  
  
  ### Now create dummies
  
  ### Similarly to the dummies at baseline
  ### The first dummy is whether someone has ever smoked
  ### The second dummy is whether someone has stopped smoking
  
  ### Note because dummaies are defined by med_status_adj_smoking, as opposed to smoking_baseline, this will change during follow-up
  data <- dplyr::mutate(data, 
                        med_status_adj_smoking_dummy1 = dplyr::case_when(med_status_adj_smoking == 0 ~ 0,
                                                                         TRUE ~ 1),
                        med_status_adj_smoking_dummy2 = dplyr::case_when(med_status_adj_smoking == 1 ~ 1,
                                                                         TRUE ~ 0))
  
  ### However, we want to adjust these dependent on individuals status at baseline
  ### We only want to change these variables if they change from the value at baseline
  
  ### i.e. med_status_adj_smoking_dummy1 and med_status_adj_smoking_dummy2 will both start at zero, and only change when
  ### the dummy status changes
  
  ### Dummy1
  ### The first dummy = 1 will only ever go from 0 to 1. We therefore want to set to zero, if it started at 1
  
  ### If smoking_baseline == 1 or 2 (smoking_baseline_dummy1 = 1), 
  ###   we want to make this a zero (as someone will never experience initiating smoking during follow-up).
  ### If smoking_baseline == 0 (smoking_baseline_dummy1 = 0), 
  ###   we leave it as it is.
  data <- dplyr::mutate(data, 
                        med_status_adj_smoking_dummy1 = dplyr::case_when(smoking_baseline_dummy1 == 1 ~ 0,
                                                                         TRUE ~ med_status_adj_smoking_dummy1))
  ### Dummy2
  
  ### The second dummy = 1 if a current smoker (baseline) stops smoking during follow-up, 
  ### The second dummy = -1 if an ex-smoker (at baseline) starts smoking during follow-up
  ### The second dummy = 1 if an never-smoker (at baseline) starts and then stops smoking during follow-up
  
  ### if an individuals starts smoking during follow-up.
  ### The same individual will never go from -1 to 1, it's only change from whatever their status is at baseline
  ###
  ### e.g. 
  ### If an individual is an ex-smoker at baseline, they have
  ### smoking_baseline == 1 (smoking_baseline_dummy2 = 1).
  ### If they start smoking again, we set med_status_adj_smoking_dummy2 = -1. If they stop, we set med_status_adj_smoking_dummy2 = 0. 
  ###
  ### If an individual is an smoker at baseline, they have
  ### smoking_baseline == 1 (smoking_baseline_dummy2 = 0).
  ### If they stop smoking again, we set med_status_adj_smoking_dummy2 = 1. If they start again, we set med_status_adj_smoking_dummy2 = 0.
  ### 
  ### Note that his means the sum of med_status_adj_smoking_dummy2 and smoking_baseline_dummy2 will always be either 0 or 1.
  ### This means we are never double applying the effects of starting or stopping smoking.
  ###
  data <- dplyr::mutate(data, 
                        med_status_adj_smoking_dummy2 = 
                          dplyr::case_when(smoking_baseline_dummy2 == 1 ~ med_status_adj_smoking_dummy2 - 1,
                                           TRUE ~ med_status_adj_smoking_dummy2))
  
  return(data)
}

###
### Function to prepare the data for prototype3, model B
###
### More details on the definition of the dummies is given inline
prepare_data_prototype3_modelB <- function(data){
  
  ###
  ### Creation of baseline smoking variables
  ###
  
  ### The first dummy variable = 0 if never smoker, = 1 if current or past smoker.
  ### NB: This dummy effect (HR1) is the effect of current smoker vs never smoker
  ###
  ### The second dummy variable = 0 if never or current smoker, = 1 if past smoker
  ### NB: This dummy effect (HR2) is the effect of smoking cessation (stopping smoking, past smoker vs current smoker)
  ###
  ### These dummies are only appropriate when applied in tandem.
  ### A never smoker will have: LP*1*1
  ### A current smoker will have: LP*HR1*1
  ### A past smoker will have: LP*HR1*HR2, i.e., the effect of initiating smoking, then the effect of stopping
  ###
  data <-
    dplyr::mutate(data,
                  smoking_baseline = dplyr::case_when(smoking == "Non-smoker" ~ 0,
                                                      smoking == "Ex-smoker" ~ 1,
                                                      !(smoking %in% c("Non-smoker", "Ex-smoker")) ~ 2),
                  smoking_baseline_dummy1 = dplyr::case_when(smoking_baseline == 0 ~ 0,
                                                             TRUE ~ 1),
                  smoking_baseline_dummy2 = dplyr::case_when(smoking_baseline == 1 ~ 1,
                                                             TRUE ~ 0))
  
  
  ###
  ### Creation of change in smoking variables during follow-up (TIME VARYING)
  ###
  
  ### med_status_adj_smoking is a time-varying variable for smoking status during follow-up, however some are missing
  ### at baseline. These have been imputed, and the imputed values are in the "smoking" variable. They are missing from
  ### med_status_adj_smoking because we calculated the time-varying variable for all time points during follow-up (and not baseline), for 
  ### computational reasons.
  ###
  ### Add in imputed (or observed) smoking values at baseline to the time varying smoking variable
  ### If med_status_adj_smoking = NA, this is the first value of smoking status, for which we didn't know when
  ### create the cohort split times (as it was imputed for a large number of individual, and varies between imputed datasets)
  ### Therefore take the value of smoking, when med_status_adj_smoking is NA
  data <- 
    dplyr::mutate(data,
                  med_status_adj_smoking = 
                    dplyr::case_when(is.na(med_status_adj_smoking) & smoking == "Non-smoker" ~ 0,
                                     is.na(med_status_adj_smoking) & smoking == "Ex-smoker" ~ 1,
                                     is.na(med_status_adj_smoking) & !(smoking %in% c("Non-smoker", "Ex-smoker")) ~ 2,
                                     TRUE ~ med_status_adj_smoking))
  
  
  ### Now create dummies
  
  ### Similarly to the dummies at baseline
  ### The first dummy is whether someone has ever smoked
  ### The second dummy is whether someone has stopped smoking
  
  ### Note because dummaies are defined by med_status_adj_smoking, as opposed to smoking_baseline, this will change during follow-up
  data <- dplyr::mutate(data, 
                        med_status_adj_smoking_dummy1 = dplyr::case_when(med_status_adj_smoking == 0 ~ 0,
                                                                         TRUE ~ 1),
                        med_status_adj_smoking_dummy2 = dplyr::case_when(med_status_adj_smoking == 1 ~ 1,
                                                                         TRUE ~ 0))
  
  ### However, we want to adjust these dependent on individuals status at baseline
  ### We only want to change these variables if they change from the value at baseline
  
  ### i.e. med_status_adj_smoking_dummy1 and med_status_adj_smoking_dummy2 will both start at zero, and only change when
  ### the dummy status changes
  
  ### Dummy1
  ### The first dummy = 1 will only ever go from 0 to 1. We therefore want to set to zero, if it started at 1
  
  ### If smoking_baseline == 1 or 2 (smoking_baseline_dummy1 = 1), 
  ###   we want to make this a zero (as someone will never experience initiating smoking during follow-up).
  ### If smoking_baseline == 0 (smoking_baseline_dummy1 = 0), 
  ###   we leave it as it is.
  data <- dplyr::mutate(data, 
                        med_status_adj_smoking_dummy1 = dplyr::case_when(smoking_baseline_dummy1 == 1 ~ 0,
                                                                         TRUE ~ med_status_adj_smoking_dummy1))
  ### Dummy2
  
  ### The second dummy = 1 if a current smoker (baseline) stops smoking during follow-up, 
  ### The second dummy = -1 if an ex-smoker (at baseline) starts smoking during follow-up
  ### The second dummy = 1 if an never-smoker (at baseline) starts and then stops smoking during follow-up
  
  ### if an individuals starts smoking during follow-up.
  ### The same individual will never go from -1 to 1, it's only change from whatever their status is at baseline
  ###
  ### e.g. 
  ### If an individual is an ex-smoker at baseline, they have
  ### smoking_baseline == 1 (smoking_baseline_dummy2 = 1).
  ### If they start smoking again, we set med_status_adj_smoking_dummy2 = -1. If they stop, we set med_status_adj_smoking_dummy2 = 0. 
  ###
  ### If an individual is an smoker at baseline, they have
  ### smoking_baseline == 1 (smoking_baseline_dummy2 = 0).
  ### If they stop smoking again, we set med_status_adj_smoking_dummy2 = 1. If they start again, we set med_status_adj_smoking_dummy2 = 0.
  ### 
  ### Note that his means the sum of med_status_adj_smoking_dummy2 and smoking_baseline_dummy2 will always be either 0 or 1.
  ### This means we are never double applying the effects of starting or stopping smoking.
  ###
  data <- dplyr::mutate(data, 
                        med_status_adj_smoking_dummy2 = 
                          dplyr::case_when(smoking_baseline_dummy2 == 1 ~ med_status_adj_smoking_dummy2 - 1,
                                           TRUE ~ med_status_adj_smoking_dummy2))
  
  
  ### SBP relative to 120
  ### BMI relative to 19.5
  ### nonhdl relative to 4
  data <- dplyr::mutate(data, 
                        sbp_adj = sbp - 120,
                        bmi_adj = bmi - 19.5,
                        nonhdl_adj = nonhdl - 4
  )
  
  ### Sort by patid and tstart, just for when looking at the dataset
  data <- dplyr::arrange(data, patid, tstart)
  
  ### Create appropriate offsets based on treatment effects
  data$offset_statins_timevar_lnHR <- lnHR_statins*data$med_status_adj_statins
  data$offset_ah_timevar_lnHR <- lnHR_ah*data$med_status_adj_ah
  data$offset_smoking_dummy1_lnHR <- lnHR_smoking_dummy1_direct*data$smoking_baseline_dummy1
  data$offset_smoking_dummy2_lnHR <- lnHR_smoking_dummy2_direct*data$smoking_baseline_dummy2
  data$offset_smoking_timevar_dummy1_lnHR <- lnHR_smoking_dummy1_total*data$med_status_adj_smoking_dummy1
  data$offset_smoking_timevar_dummy2_lnHR <- lnHR_smoking_dummy2_total*data$med_status_adj_smoking_dummy2
  data$offset_sbp_lnHR <- lnHR_sbp*data$sbp_adj
  data$offset_bmi_lnHR <- lnHR_bmi*data$bmi_adj
  data$offset_nonhdl_lnHR <- lnHR_nonhdl*data$nonhdl_adj
  
  return(data)
}
