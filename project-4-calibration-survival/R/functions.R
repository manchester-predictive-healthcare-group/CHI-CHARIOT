
#####################
### DGM FUNCTIONS ###
#####################

###
### General DGM
###
simulate_DGM <- function(m, 
                         cov,
                         lambda.y = 1, gamma.y = 1, lambda.cens = 1, gamma.cens = 1,
                         b.x1 = 0, b.x2 = 0, b.x3 = 0, b.x2sq = 0, b.x3log = 0, b.x4 = 0,
                         b.cens.x1 = 0, b.cens.x2 = 0, b.cens.x3 = 0, b.cens.x2sq = 0, b.cens.x3log = 0, b.cens.x4 = 0,
                         seed = NULL,
                         cens = TRUE){
  
  ### set seed
  if (!is.null(set.seed)){
    set.seed(seed)
  }
  
  ### Add some non-linear terms
  cov <- dplyr::mutate(cov, 
                       x2sq = x2^2,
                       x3log = x3^3)
  
  ### Simulate event times
  dat <- simsurv::simsurv(lambdas = lambda.y, 
                          gammas = gamma.y,
                          betas = c(x1 = b.x1, x2 = b.x2, x3 = b.x3, x2sq = b.x2sq, x3log = b.x3log, x4 = b.x4),
                          x = cov,
                          maxt = 100) |>
    dplyr::rename(time = eventtime)
  
  ### Add a censoring time
  ### Simulate event times
  if (cens == TRUE){
    cens.times <- simsurv::simsurv(lambdas = lambda.cens, 
                                   gammas = gamma.cens,
                                   betas = c(x1 = b.cens.x1, x2 = b.cens.x2, x3 = b.cens.x3,
                                             x2sq = b.cens.x2sq, x3log = b.cens.x3log,
                                             x4 = b.cens.x4),
                                   x = cov)
    
    cens.times <- dplyr::rename(cens.times, cens_time = eventtime) |>
      dplyr::select(-c(id,status))
    
    ### Combine with covariates and cens.times
    dat <- cbind(dat, cov, cens.times)
    
    ### Apply censoring time to the event time data
    dat <- dplyr::mutate(dat, 
                         status = dplyr::case_when(cens_time < time ~ 0,
                                                   TRUE ~ 1),
                         time = dplyr::case_when(cens_time < time ~ cens_time,
                                                 TRUE ~ time))
    
    ### Create a variable for censored or not censored 
    ### (if an individual doesn't have an event, they are censored at end of follow up)
    ### This will be used in the IPCW method
    dat$cens_time <- dat$time
    dat$cens_indicator <- 1 - dat$status
  } else {
    ### Add cens_time to be same as event time, and always uncensored
    dat$cens_time <- dat$time
    dat$cens_indicator <- 0
    ### Combine with covariates
    dat <- cbind(dat, cov)
  }
  
  ### Calculate true linear predictor and add to dat (this will be used to calculating true risks in simulation)
  dat$true.lp <- as.numeric(as.matrix(cov[,1:6]) %*% c(b.x1, b.x2, b.x3, b.x4, b.x2sq, b.x3log))
  dat$true.lp.cens <- as.numeric(as.matrix(cov[,1:6]) %*% c(b.cens.x1, b.cens.x2, b.cens.x3, b.cens.x4, b.cens.x2sq, b.cens.x3log))
  
  return(dat)
  
}


###
### DGM function for supplementary analyses, exploring how bias induced by unmeasured confounding differs
### depending on the amount of covariance in the linear predictor.
###
simulate_DGM_supplementary <- function(m, 
                                       cov,
                                       lambda.y = 1, gamma.y = 1, lambda.cens = 1, gamma.cens = 1,
                                       b.x1 = 0, b.x2 = 0, b.x3 = 0, b.x4 = 0, b.x5 = 0, b.x.cens = 0,
                                       seed = NULL,
                                       cens = TRUE){
  
  ### set seed
  if (!is.null(set.seed)){
    set.seed(seed)
  }
  
  ### Simulate event times
  dat <- simsurv::simsurv(lambdas = lambda.y, 
                          gammas = gamma.y,
                          betas = c(x1 = b.x1, x2 = b.x2, x3 = b.x3, x4 = b.x4, x5 = b.x5, x.cens = b.x.cens),
                          x = cov,
                          maxt = 100) |>
    dplyr::rename(time = eventtime)
  
  ### Add a censoring time
  ### Simulate event times
  if (cens == TRUE){
    cens.times <- simsurv::simsurv(lambdas = lambda.cens, 
                                   gammas = gamma.cens,
                                   betas = c(x.cens = b.x.cens),
                                   x = cov)
    
    cens.times <- dplyr::rename(cens.times, cens_time = eventtime) |>
      dplyr::select(-c(id,status))
    
    ### Combine with covariates and cens.times
    dat <- cbind(dat, cov, cens.times)
    
    ### Apply censoring time to the event time data
    dat <- dplyr::mutate(dat, 
                         status = dplyr::case_when(cens_time < time ~ 0,
                                                   TRUE ~ 1),
                         time = dplyr::case_when(cens_time < time ~ cens_time,
                                                 TRUE ~ time))
    
    ### Create a variable for censored or not censored 
    ### (if an individual doesn't have an event, they are censored at end of follow up)
    ### This will be used in the IPCW method
    dat$cens_time <- dat$time
    dat$cens_indicator <- 1 - dat$status
  } else {
    ### Add cens_time to be same as event time, and always uncensored
    dat$cens_time <- dat$time
    dat$cens_indicator <- 0
    ### Combine with covariates
    dat <- cbind(dat, cov)
  }
  
  ### Calculate true linear predictor and add to dat (this will be used to calculating true risks in simulation)
  dat$true.lp <- as.numeric(as.matrix(cov[,1:6]) %*% c(b.x1, b.x2, b.x3, b.x4, b.x5, b.x.cens))
  dat$true.lp.cens <- as.numeric(as.matrix(cov[,6]) %*% c(b.x.cens))
  
  return(dat)
  
}

###########################################
###########################################
### FUNCTIONS FOR ASSESSING CALIBRATION ###
###########################################
###########################################

###
### Write a function to estimate survival probabilities based on a fitted model (fit) and baseline hazard (bhaz).
### Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE).
###
est_surv <- function(newdata, fit, bhaz = NULL, t){
  
  ### Get bhaz if not specified
  if (is.null(bhaz)){bhaz <- survival::basehaz(fit, centered = TRUE)}
  
  ### Get the lp
  lp <- predict(fit, newdata = newdata, reference = "sample")
  
  ### Get the linear predictor for ne wdata
  surv <- as.numeric(exp(-exp(lp)*bhaz$hazard[max(which(bhaz$time <= t))]))
  
  return(surv)
  
}

###
### Assess calibration by estimating Kaplan-Meier (observed) and mean predicted risk (predicted) within subgroups, defined by predicted risk
###
est_calib_plot_group <- function(data, fit, bhaz, t, n.groups, surv = NULL){
  
  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Split data by surv
  data.split <- split(data, cut(data$surv, c(-Inf, quantile(data$surv, probs = 1:n.groups/n.groups))))
  
  ### Create vectors to store estimates
  obs <- vector("numeric", n.groups)
  pred <- vector("numeric", n.groups)
  
  ### Run through
  for (i in 1:n.groups){
    
    ### Get Kaplan-Meier estimate
    survobj <- survival::survfit(survival::Surv(time, status) ~ 1, data = data.split[[i]])
    obs[i] <- as.numeric(survobj$surv[max(which(survobj$time <= t))])
    
    ### Get predicted
    pred[i] <- 1 - (mean(data.split[[i]]$surv))
    
  }
  
  ### Create plot data
  plot.data <- data.frame("obs" = 1 - obs, "pred" = pred)
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed")
  
  ### Create output.object
  output.object <- list("plot" = plot)
  
  return(output.object)
  
}


###
### Assess mean calibration by estimating Kaplan-Meier (observed) and mean predicted risk (predicted)
###
est_calib_mean_km <- function(data, fit, bhaz, t, surv = NULL){
  
  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Create vectors to store estimates
  survobj <- survival::survfit(survival::Surv(time, status) ~ 1, data = data)
  obs <- 1 - as.numeric(survobj$surv[max(which(survobj$time <= t))])
  
  ### Get predicted
  pred <- 1 - (mean(data$surv))
  
  ### Calculate ratio
  ratio <- obs/pred
  
  return(ratio)
  
}


###
### Assessing calibration using proportional hazards regression approach (graphical calibration curves, Austin et al)
###
est_calib_ph <- function(data, fit, bhaz, t, surv = NULL, nk = 4, pred.plot.range = NULL, plot = TRUE){
  
  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Add complementary log-log of predicted survival probabilities to data
  data$cloglog <- log(-log(data$surv))
  
  ### Fit calibration model
  fit.calib <- survival::coxph(survival::Surv(time, status) ~ rms::rcs(cloglog, nk), data = data)
  bhaz.calib <- survival::basehaz(fit.calib, centered = TRUE)
  
  ###
  ### Generate predicted observed values,
  ###
  
  ### Create a dataframe with the calibration data
  pred.obs <- est_surv(newdata = data, fit = fit.calib, bhaz = bhaz.calib, t = t)
  output.data <- data.frame("id" = data$id, "pred.obs" = 1 - pred.obs, "pred" = 1 - data$surv)
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(output.data$pred.obs - output.data$pred))
  E50 <- median(abs(output.data$pred.obs - output.data$pred))
  E90 <- as.numeric(quantile(abs(output.data$pred.obs - output.data$pred), probs = .9, na.rm = TRUE))
  
  ### Create dataframe for plotting
  ### Do so over the range of values pred.plot.range, if specified
  if (!is.null(pred.plot.range)){
    ### Create temporary data frame
    tmp.data <- data.frame("pred" = pred.plot.range, cloglog = log(-log(1 - pred.plot.range)))
    
    ### Calculate predicted observed values over pred.plot.range
    pred.obs <- 1 - est_surv(newdata = tmp.data, fit = fit.calib, bhaz = bhaz.calib, t = t)
    
    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = tmp.data$pred)
    
  } else {
    
    ### Create plot data
    plot.data <- output.data
    
    ### Now set output.data to NA, 
    ### this is so we don't save it twice (it will be an element in the plot object, which will be outputted)
    output.data <- NA
  }
  
  ### Create plot
  if (plot == TRUE){
    plot <- ggplot2::ggplot(data = plot.data) +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") + 
      ggplot2::ggtitle("Proportional hazards") + 
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  } else {
    plot <- NULL
  }
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "calib_data" = output.data)
  
  return(output.object)
  
}

###
### Function to estimate calibration using IPCW
###

### First Write a function to estimate survival probability for a given model and fitted baseline hazard,
### but survival probability is probability at the time they had their event.
### Variable eventtime is a character string, of the name of the variable denoting the event time.
### Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE).
### This is used to calculate weights for individuals who have an event, and therefore we need to calculate
### probability of being uncensored at the time they had their event.
est_surv_eventtime <- function(newdata, bhaz, fit, eventtime, type = "coxph"){
  
  if (type == "coxph"){
    
    ### Get lp
    lp <- predict(fit, newdata = newdata, type = "lp", reference = "sample")
    
    ### Get bhaz specific for each patient at the time they had their event
    bhaz.pat <- unlist(lapply(1:nrow(newdata), function(x) {bhaz$hazard[max(which(bhaz$time <= newdata[x, eventtime]))]}))
    
    ### Get expected
    surv <- as.numeric(exp(-exp(lp)*bhaz.pat))
    
  } else if (type == "flexsurv"){
    
    # Calculate survival probabilities for these individuals
    surv <- unlist(lapply(1:nrow(newdata), 
                          function(x) {dplyr::pull(predict(fit, 
                                                           newdata = newdata[x, ], 
                                                           times = as.numeric(newdata[x, eventtime]), 
                                                           type = "survival"), 
                                                   .pred_survival)}
    ))
    
  }
  
  return(surv)
  
}

###
### Function to estimate survival time, at minimum of a variable (eventtime) or fixed numeric time (t)
### Note, this is same as est_surv_eventtime, but will estiamte survival probabilities at the minimum
### of the variable eventtime, or the numeric time t. This is used is est_calib_pv_eventglm
###
est_surv_eventtime_or_t <- function(newdata, bhaz, fit, eventtime, t, type = "coxph"){
  
  if (type == "coxph"){
    
    ### Get lp
    lp <- predict(fit, newdata = newdata, type = "lp", reference = "sample")
    
    ### Get bhaz specific for each patient at the time they had their event
    bhaz.pat <- unlist(lapply(1:nrow(newdata), function(x) {bhaz$hazard[max(which(bhaz$time <= min(newdata[x, eventtime], t)))]}))
    
    ### Get expected
    surv <- as.numeric(exp(-exp(lp)*bhaz.pat))
    
  } else if (type == "flexsurv"){
    
    # Calculate survival probabilities for these individuals
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
### Function to estimate calibration using IPCW
###
est_calib_ipcw <- function(data, 
                           fit, 
                           bhaz, 
                           t, 
                           surv = NULL, 
                           cens.max.follow = NULL, 
                           nk = 4, 
                           ipcw.formula, 
                           type = "coxph", 
                           pred.plot.range = NULL,
                           plot = TRUE){
  
  # data = dat.valid
  # fit = fit1
  # bhaz = bhaz1
  # t = t.eval
  # surv = NULL
  # nk = 3
  # cens.max.follow = NULL
  # ipcw.formula <- as.formula("Surv(cens_time, cens_indicator) ~ age + sex + ph.ecog + ph.karno")
  # 
  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Assign a variable for whether individuals have had an event by time of interest
  data <- dplyr::mutate(data,
                        status_time = dplyr::case_when(time <= t ~ status,
                                                       time > t ~ 0))
  
  ### Estimate weights
  weights <- est_ipcw(data = data, t = t, cens.max.follow = cens.max.follow, cens.formula = ipcw.formula, type = type)
  
  ### Merge with data
  data <- merge(data, weights, by.x = "id", by.y = "id")
  
  ### Convert predict risk onto logit scale
  data$pred <- 1 - data$surv
  data$pred.logit <- log(data$pred/(1 - data$pred))
  
  ### Reduce data to individuals uncensored at time t for fitting the model
  data.model <- subset(data, !is.na(ipcw))
  
  ### Fit weighted calibration model and generate predicted observed
  ## Fit model
  rcs.model.stab <- suppressWarnings(glm(status_time ~ rms::rcs(pred.logit, nk),
                                         family = binomial(link = "logit"),
                                         data = data.model,
                                         weights = data.model[, "ipcw.stab"]))
  ## Supress warnings due to using weights in binomial glm results in 
  ## "Warning message: In eval(family$initialize) : non-integer #successes in a binomial glm!"
  
  ###
  ### Generate predicted observed values,
  ###
  
  ### Create a dataframe with the calibration data
  pred.obs <- predict(rcs.model.stab, newdata = data, type = "response")
  output.data <- data.frame("id" = data$id, "pred.obs" = pred.obs, "pred" = data$pred)
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(output.data$pred.obs - output.data$pred))
  E50 <- median(abs(output.data$pred.obs - output.data$pred))
  E90 <- as.numeric(quantile(abs(output.data$pred.obs - output.data$pred), probs = .9, na.rm = TRUE))
  
  ### Create dataframe for plotting
  ### Do so over the range of values pred.plot.range, if specified
  if (!is.null(pred.plot.range)){
    ### Create temporary data frame
    tmp.data <- data.frame("pred" = pred.plot.range, "pred.logit" = log(pred.plot.range/(1-pred.plot.range)))
    
    ### Calculate predicted observed values over pred.plot.range
    pred.obs <- predict(rcs.model.stab, newdata = tmp.data, type = "response")
    
    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = tmp.data$pred)
    
  } else {
    
    ### Create plot data
    plot.data <- output.data
    
    ### Now set output.data to NA, 
    ### this is so we don't save it twice (it will be an element in the plot object, which will be outputted)
    output.data <- NA
  }
  
  ### Create plot
  if (plot == TRUE){
    plot <- ggplot2::ggplot(data = plot.data) +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") + 
      ggplot2::ggtitle("IPCW") + 
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  } else {
    plot <- NULL
  }
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "calib_data" = output.data)
  
  ### Return
  return(output.object)
  
}

###
### Write a function to estimate IPCW's that can be used in multiple other functions
###
est_ipcw <- function(data, t, cens.max.follow = NULL, cens.formula, type){
  
  ### If there is informative censoring (i.e. all individuals censored after 11 years)
  ### Don't want to model the effect of everybody becoming censored at the final follow up time,
  ### so set these individuals to be uncensored at time point of interest.
  ### This will not be used in the simulation, but can occur in real data so leaving in the functionality.
  if (!is.null(cens.max.follow)){
    data <- dplyr::mutate(data,
                          cens_indicator = dplyr::case_when(cens_time < cens.max.follow + 2 ~ cens_indicator,
                                                            cens_time >= cens.max.follow + 2 ~ 0),
                          cens_time = dplyr::case_when(cens_time < cens.max.follow + 2 ~ cens_time,
                                                       cens_time >= cens.max.follow + 2 ~ cens.max.follow + 2))
  }
  
  ### Create censoring model and create survival predictions from this model, probability of being uncensored at time t
  if (type == "coxph"){
    
    ### Create censoring model
    cens.model <- survival::coxph(cens.formula, data = data, model = TRUE)
    cens.model.int <- survival::coxph(survival::Surv(cens_time, cens_indicator) ~ 1, data = data)
    
    ### Create survival predictions from this model, probability of being uncensored at time t
    cens.bhaz <- survival::basehaz(cens.model, centered = TRUE)
    data$surv.cens <- est_surv(newdata = data, fit = cens.model, bhaz = cens.bhaz, t = t)
    
  } else if (type == "flexsurv"){
    
    ### Create censoring model
    cens.model <- flexsurv::flexsurvreg(cens.formula, data = data, dist = "weibull")
    cens.model.int <- flexsurv::flexsurvreg(survival::Surv(cens_time, cens_indicator) ~ 1, data = data, dist = "weibull")
    
    ### Create survival predictions from this model, probability of being uncensored at time t
    data$surv.cens <- dplyr::pull(predict(cens.model, newdata = data, times = t, type = "survival"), .pred_survival)
    
  }
  
  ### For individual who have an event prior to or equal to time t (and are therefore uncensored at time t),
  ### we want probability of being uncensored at t_event. For people censored before or on time t,
  ### they will not be included in the analysis, as they are censored.
  obs.event.prior <- which(data$cens_time <= t & data$cens_indicator == 0)
  obs.censored.prior <- which(data$cens_time <= t & data$cens_indicator == 1)
  
  ### Calculate weights for individuals who have events prior to time t
  surv.event.prior <- est_surv_eventtime(newdata = data[obs.event.prior, ], 
                                         bhaz = cens.bhaz, 
                                         fit = cens.model, 
                                         eventtime = "cens_time",
                                         type = type)
  
  ### Replace these values in dataset
  data$surv.cens[obs.event.prior] <- surv.event.prior
  
  ### Create weights
  data$ipcw <- 1/data$surv.cens
  
  ### Stabilise
  ### Create censoring model and create survival predictions from this model, probability of being uncensored at time t
  if (type == "coxph"){
    # Get bhaz
    cens.bhaz.int <- survival::basehaz(cens.model.int, centered = TRUE)
    # Get survival probabilities from intercept only model
    data$ipcw.numer <- est_surv(newdata = data, fit = cens.model.int, bhaz = cens.bhaz.int, t = t)
    # Apply stabilisation
    data$ipcw.stab <- data$ipcw.numer*data$ipcw
  } else if (type == "flexsurv"){
    # Get survival probabilities from intercept only model
    data$ipcw.numer <- dplyr::pull(predict(cens.model.int, newdata = data, times = t, type = "survival"), .pred_survival)
    # Apply stabilisation
    data$ipcw.stab <- data$ipcw.numer*data$ipcw
  }
  
  ### Cap the weights at p1 and p99
  data$ipcw <- pmax(data$ipcw, as.numeric(quantile(data$ipcw, p = .01)))
  data$ipcw <- pmin(data$ipcw, as.numeric(quantile(data$ipcw, p = .99)))
  
  data$ipcw.stab <- pmax(data$ipcw.stab, as.numeric(quantile(data$ipcw.stab, p = .01)))
  data$ipcw.stab <- pmin(data$ipcw.stab, as.numeric(quantile(data$ipcw.stab, p = .99)))
  
  ### Cap at 100
  data$ipcw <- pmin(data$ipcw, 100)
  data$ipcw.stab <- pmin(data$ipcw.stab, 100)
  
  ### Assign NA's to individuals censored prior to time t
  data$ipcw[obs.censored.prior] <- NA
  data$ipcw.stab[obs.censored.prior] <- NA
  
  ### Filter
  data <- dplyr::select(data, id, ipcw, ipcw.stab)
  
  return(data)
  
}


###
### Write a function to estimate IPCW's that can be used to group individuals in the pseudo-value approach
### Note here: we do not want to estimate the weights at min(t, t_abs)
### We are just interested in the probability of being censored assuming individuals are followed up for the same amount of time
###
est_ipcw_pv <- function(data, t, cens.max.follow = NULL, cens.formula, type){
  
  ### If there is informative censoring (i.e. all individuals censored after 11 years)
  ### Don't want to model the effect of everybody becoming censored at the final follow up time,
  ### so set these individuals to be uncensored at time point of interest.
  ### This will not be used in the simulation, but can occur in real data so leaving in the functionality.
  if (!is.null(cens.max.follow)){
    data <- dplyr::mutate(data,
                          cens_indicator = dplyr::case_when(cens_time < cens.max.follow + 2 ~ cens_indicator,
                                                            cens_time >= cens.max.follow + 2 ~ 0),
                          cens_time = dplyr::case_when(cens_time < cens.max.follow + 2 ~ cens_time,
                                                       cens_time >= cens.max.follow + 2 ~ cens.max.follow + 2))
  }
  
  ### Create censoring model and create survival predictions from this model, probability of being uncensored at time t
  if (type == "coxph"){
    
    ### Create censoring model
    cens.model <- survival::coxph(cens.formula, data = data, model = TRUE)
    cens.model.int <- survival::coxph(survival::Surv(cens_time, cens_indicator) ~ 1, data = data)
    
    ### Create survival predictions from this model, probability of being uncensored at time t
    cens.bhaz <- survival::basehaz(cens.model, centered = TRUE)
    data$surv.cens <- est_surv(newdata = data, fit = cens.model, bhaz = cens.bhaz, t = t)
    
  } else if (type == "flexsurv"){
    
    ### Create censoring model
    cens.model <- flexsurv::flexsurvreg(cens.formula, data = data, dist = "weibull")
    cens.model.int <- flexsurv::flexsurvreg(survival::Surv(cens_time, cens_indicator) ~ 1, data = data, dist = "weibull")
    
    ### Create survival predictions from this model, probability of being uncensored at time t
    data$surv.cens <- dplyr::pull(predict(cens.model, newdata = data, times = t, type = "survival"), .pred_survival)
    
  }
  
  ### Create weights
  data$ipcw <- 1/data$surv.cens
  
  ### Filter
  data <- dplyr::select(data, id, ipcw)
  
  return(data)
  
}



###
### Function to estimate calibration using pseudo-value appraoch, where pseudo-values estimated using
### the kaplan-meier estimate of survival
###

### First write a function which will estimate pseudo-values for a given dataset at a time t
est_pv <- function(df, t){
  
  ### Create prodlim object
  prodlim.obj <- prodlim::prodlim(prodlim::Hist(time, status) ~ 1, data = df)
  
  ### Calculate psuedo-values
  pv <- prodlim::jackknife(prodlim.obj, times = t)
  
  return(as.numeric(pv))
}

### Write a function to calculate a pseudo-value for a given individual, calculated within a window defined by the ipcws
est_pv_ind <- function(data, row, window_size, window_min_n){
  
  ### Get the weight of individual of interest
  ind_ipcw <- data$ipcw[row]
  ind_id <- data$id[row]
  
  ### Reduce data to individuals with a weight within window distance
  ind_data <- subset(data, ipcw > ind_ipcw - window_size & ipcw < ind_ipcw + window_size)
  ### Get this reduced cohort, excluding the indiviudal of interest
  ind_data_minus1 <- subset(data[-row,], ipcw > ind_ipcw - window_size & ipcw < ind_ipcw + window_size)
  
  ### If this cohort has less than 20 people, just pick the 20 that are closest?
  if (nrow(ind_data) < window_min_n){
    ind_data <- dplyr::mutate(data, distance = abs(ipcw - ind_ipcw)) |>
      dplyr::arrange(distance) |>
      dplyr::slice(1:window_min_n)
    ind_data_minus1 <- dplyr::filter(ind_data, id != ind_id)
  }
  
  ### Calculate and return the pseudo-value
  return(nrow(ind_data)*predict(prodlim::prodlim(prodlim::Hist(time, status) ~ 1, data = ind_data), times = t) - 
           nrow(ind_data_minus1)*predict(prodlim::prodlim(prodlim::Hist(time, status) ~ 1, data = ind_data_minus1), times = t))
}


### Now write the main function
est_calib_pv <- function(data, fit, bhaz, t, surv = NULL, nk = 4, split.n.groups = 1, 
                         group.by = c("predrisk", "ipcw"), 
                         ipcw.cens.max.follow = NULL, ipcw.formula = NULL, ipcw.type = "coxph",
                         rolling = FALSE, rolling_window_size = 0.05, rolling_window_min_n = 20,
                         group.var = NULL, pred.plot.range = NULL, plot = TRUE){
  
  # data = dat.valid
  # fit  = fit1
  # bhaz = bhaz1
  # t = t.eval
  # nk = 5
  # split.n.groups = 25
  # surv <- NULL
  # group.var = NULL
  # group.by = "ipcw"
  # ipcw.cens.max.follow = 0.5
  # ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ rms::rcs(x1,5) + rms::rcs(x2,5) + rms::rcs(x3,5)")
  # ipcw.type = "coxph"
  
  ### Create internal id
  ### This will be required when bootstrapping, as using patient ID's from the data will create issues when
  ### merging pseudo-values with the data, when some individual have been sampled more than once
  data$id_internal <- 1:nrow(data)
  
  ### Match.arg
  group.by <- match.arg(group.by)
  
  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Add predicted risk
  data$pred <- 1 - data$surv
  
  ### Split continuously if group.var is null
  if (rolling == FALSE){
    if (is.null(group.var)){
      ### Split data by survival probabilities
      if (group.by == "predrisk"){
        ### Split
        data.pv.split <- split(data, 
                               cut(data$surv, 
                                   c(-Inf, quantile(data$surv, probs = 1:split.n.groups/split.n.groups))))
        ### Split data by ipcws
      } else if (group.by == "ipcw"){
        ### Estimate weights
        weights <- est_ipcw_pv(data = data, t = t, cens.max.follow = ipcw.cens.max.follow, cens.formula = ipcw.formula, type = ipcw.type)
        
        ### Merge with data
        data <- merge(data, weights, by.x = "id_internal", by.y = "id_internal")
        
        ### Split
        data.pv.split <- split(data, 
                               cut(data$ipcw, 
                                   c(-Inf, quantile(data$ipcw, probs = 1:split.n.groups/split.n.groups))))
        ### Not uncommon for the final few groups to have no events
        ### Want to combine these, until the final group has an event
        ## Set counting variables
        split.n.groups.ticker <- split.n.groups
        n.events.final.group <- sum(data.pv.split[[split.n.groups.ticker]]$status)
        ## Create loop
        while(n.events.final.group <= 0){
          
          ## Reduce ticker by 1
          split.n.groups.ticker <- split.n.groups.ticker - 1
          ## combine final two groups
          data.pv.split[[split.n.groups.ticker]] <- rbind(data.pv.split[[split.n.groups.ticker]], data.pv.split[[split.n.groups.ticker+1]])
          ## Remove final group
          data.pv.split <- data.pv.split[-(split.n.groups.ticker+1)]
          ## Calculate n events
          n.events.final.group <- sum(data.pv.split[[split.n.groups.ticker]]$status)
        }
      }
    }
    
    ### Calculate pseudo-values for each data split (note est_pv calculates pseudo-value for survival prob, so want to take 1 - pv)
    pv <- lapply(data.pv.split, est_pv, t = t)
    pvs <- 1 - do.call("c", pv)
    # print(str(pv))
    
    ### Also get the corresponding ids
    ids <- lapply(data.pv.split, function(x) {x$id_internal})
    ids <- do.call("c", ids)
    
    ### Get in correct order
    pvs <- data.frame("pv" = pvs, "id_internal" = ids)
    pvs <- dplyr::arrange(pvs, id_internal) |> dplyr::pull(pv)
    
  } else if (rolling == TRUE){
    
    ### Estimate weights
    weights <- est_ipcw_pv(data = data, t = t, cens.max.follow = ipcw.cens.max.follow, cens.formula = ipcw.formula, type = ipcw.type)
    
    ### Merge with data
    data <- merge(data, weights, by.x = "id_internal", by.y = "id_internal")
    
    ### Apply this function to every row
    pvs <- lapply(1:nrow(data), function(x) {est_pv_ind(data = data, row = x, 
                                                        window_size = rolling_window_size, 
                                                        window_min_n = rolling_window_min_n)})
    
    ### Add the pseudo-values to validation dataset
    pvs <- unlist(pvs)
    
  }
  
  ### Add the pseudo-values to validation dataset
  data$pv <- pvs
  
  ### If pseudo-value is NA, this means all individual had had an event prior to time point t.eval
  ### By definition pseudo-values for these individuals are equal to 1, so assign these
  data$pv[is.na(data$pv)] <- 1
  
  ### Transform est.surv onto logit scale
  data$pred.logit <- log(data$pred/(1-data$pred))
  
  ### Fit the model using logit link function
  calib.model.pv <- stats::glm(pv ~ rms::rcs(pred.logit, nk), 
                               data = data, 
                               family = stats::gaussian(link = "logit"),
                               start = rep(0, nk))
  
  ###
  ### Generate predicted observed values,
  ###
  
  ### Create a dataframe with the calibration data
  pred.obs <- predict(calib.model.pv, newdata = data, type = "response")
  output.data <- data.frame("pred.obs" = pred.obs, "pred" = data$pred)
  if ("id" %in% colnames(data)){
    output.data <- dplyr::mutate(output.data, id = data$id) |>
      dplyr::relocate(id)
  }
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(output.data$pred.obs - output.data$pred))
  E50 <- median(abs(output.data$pred.obs - output.data$pred))
  E90 <- as.numeric(quantile(abs(output.data$pred.obs - output.data$pred), probs = .9, na.rm = TRUE))
  
  ### Create dataframe for plotting
  ### Do so over the range of values pred.plot.range, if specified
  if (!is.null(pred.plot.range)){
    ### Create temporary data frame
    tmp.data <- data.frame("pred" = pred.plot.range, "pred.logit" = log(pred.plot.range/(1-pred.plot.range)))
    
    ### Calculate predicted observed values over pred.plot.range
    pred.obs <- predict(calib.model.pv, newdata = tmp.data, type = "response")
    
    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = tmp.data$pred)
    
  } else {
    
    ### Create plot data
    plot.data <- output.data
    
    ### Now set output.data to NA, 
    ### this is so we don't save it twice (it will be an element in the plot object, which will be outputted)
    output.data <- NA
    
  }
  
  ### Create plot
  if (plot == TRUE){
    plot <- ggplot2::ggplot(data = plot.data) +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
      ggplot2::ggtitle(("Pseudo-value")) + 
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  } else {
    plot <- NULL
  }
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "calib_data" = output.data)
  
  return(output.object)
  
}


###
### Function to estimate calibration curve using pseudo-values with inverse probability of censoring weights,
### following the binder approach. 
###
### See:
### Binder; Pseudo-observations for competing risks with covariate dependent censoring. DOI: 10.1007/s10985-013-9257-7.
### Overgaard; Pseudo-observations under covariate-dependent censoring. DOI: 10.1016/j.jspi.2019.02.003
###
### These approaches are implemented using the eventglm package, which estimates pseudo-values using this approach.
### We estimate the weights using a cox proportional hazards model. We would like to use eventglm::pseudo_coxph directly,
### however this has computational issues in large datasets. See github.com/sachsmc/eventglm/issues/5.
### I believe this is because survival::survfit is used to estimate survival probabilities of being censored (the weights).
### 
### I have therefore taken the eventglm::pseudo_coxph function and edited it, so that the weights are estimated
### using manual functions, and it now runs in large datasets. Large amounts of the code below is taken directly
### from the eventglm GitHub page on 17/06/2025.
###
### There is an alternate function est_calib_pv_ipcw_eventglm, which uses eventglm::pseudo_coxph directly, which
### is fine in small datasets. I have written a test to showcase they lead to the same results.
###
est_calib_pv_ipcw <- function(data, fit, bhaz, t, surv = NULL, nk = 4,
                              ipcw.formula = NULL,
                              eventglm.pv.method = "binder",
                              pred.plot.range = NULL, plot = TRUE,
                              use.logit = FALSE){
  
  # data = dat.valid
  # fit  = fit1
  # bhaz = bhaz1
  # t = t.eval
  # nk = 5
  # split.n.groups = 25
  # surv <- NULL
  # group.var = NULL
  # group.by = "ipcw"
  # ipcw.cens.max.follow = 0.5
  # ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ rms::rcs(x1,5) + rms::rcs(x2,5) + rms::rcs(x3,5)")
  # ipcw.type = "coxph"
  # 
  
  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Add predicted risk
  data$pred <- 1 - data$surv
  data$pred.logit <- log(data$pred/(1 - data$pred))
  
  ###
  ### Estimate pseudo-values using binder approach from eventglm
  ### Do this manually, so we can alter the function for estimating the weights,
  ### so its more computationally efficient.
  ###
  ### The following code, is therefore predominately taken from the eventglm::pseudo_coxph function with a few edits
  ### see: github.com/sachsmc/eventglm
  ###
  
  ###
  ### Assign variable names to match those for the input parameters for eventglm::psudo_coxph
  ###
  formula = "Surv(time, status) ~ 1" #NB the formula argument is irrelevant as we are just estimating the pseudo-values
  time = t.eval
  cause = 1
  type = "cuminc"
  formula.censoring = paste("~", strsplit(as.character(ipcw.formula), "~")[[3]], sep = " ")
  ipcw.method = eventglm.pv.method
  
  ###
  ### Run code from eventglm::pseudo_coxph manually, until reach point of survfit, 
  ### where computational issues are encountered
  ###
  margformula <- update.formula(formula, . ~ 1)
  mr <- model.response(model.frame(margformula, data = data))
  stopifnot(attr(mr, "type") %in% c("right", "mright"))
  
  matcau <- match_cause(mr, cause)
  causen <- matcau$causen
  causec <- matcau$causec
  
  ###
  ### EDITS TO SOURCE CODE
  ###
  
  ### We run into issues with our data generating mechanism and eventglm.
  ### The censoring mechanism is estimated based off the outcome.
  ### However, we have applied a "stopped cox" to censor all individuals after the time point of interest.
  ### All these individuals will be viewed as "censored" when estimating the weights, and will massively impact
  ### the model. However, when estimating the probability of being censored, we don't want to include these. Its
  ### effectively administrative censoring that we have applied to the data. 
  ###
  ### NB: When using pseudo-values, you would never apply a stopped cox as its pointless, which is why simulating
  ### data this way isn't directly compatible with eventglm
  ###
  ### I am therefore going to replace .Ci and .Tci with our cens_time and cens_indicator variables
  ###
  
  # .Ci <- as.numeric(mr[, "status"] == 0)
  # .Tci <- mr[, "time"]
  
  .Ci <- data$cens_indicator
  .Tci <- data$cens_time
  
  ###
  ### EDITS OVER
  ###
  
  oldnames <- names(data)
  newnames <- make.unique(c(oldnames, c(".Ci", ".Tci")))
  
  add.nme <- newnames[length(newnames) - 1:0]
  data[[add.nme[1]]] <- c(.Ci)
  data[[add.nme[2]]] <- c(.Tci)
  
  if(is.null(formula.censoring)) {
    cens.formula <- update.formula(formula,
                                   as.formula(sprintf("survival::Surv(%s, %s) ~ .", add.nme[2], add.nme[1])))
    formula.censoring <- formula[-2]
  } else {
    cens.formula <- update.formula(formula.censoring,
                                   as.formula(sprintf("survival::Surv(%s, %s) ~ .", add.nme[2], add.nme[1])))
  }
  
  predmat <- model.matrix(cens.formula, data = data)
  
  fitcens <- survival::coxph(cens.formula, data = data, x = TRUE)
  if(!is.null(fitcens$na.action)) {
    stop("Missing data not allowed for covariates in the censoring model")
  }
  
  ###
  ### STOP running code from eventglm
  ### Estimate the weights (Gi) avoiding the use of survfit
  ###
  
  ### Extract baseline hazard in advance
  mybhaz <- basehaz(fitcens, centered = TRUE)
  
  ### Estimate Gi using function est_surv_eventtime_or_t, as opposed to survfit
  Gi.manual <- est_surv_eventtime_or_t(newdata = data, bhaz = mybhaz, fit = fitcens, eventtime = "time", t = time, type = "coxph")
  
  ### Calculate pseudo-values using unexported function from eventglm, calc_ipcw_pos
  data$pv <- calc_ipcw_pos(mr, time, causen, type, ipcw.method, Gi.manual)
  
  ### If pseudo-value is NA, this means all individual had had an event prior to time point t.eval
  ### By definition pseudo-values for these individuals are equal to 1, so assign these
  data$pv[is.na(data$pv)] <- 1
  
  ### Transform est.surv onto logit scale
  data$pred.logit <- log(data$pred/(1-data$pred))
  
  ### Fit the model using logit link function
  calib.model.pv <- stats::glm(pv ~ rms::rcs(pred.logit, nk), 
                               data = data, 
                               family = stats::gaussian(link = "logit"),
                               start = rep(0, nk))
  
  ###
  ### Generate predicted observed values,
  ###
  
  ### Create a dataframe with the calibration data
  pred.obs <- predict(calib.model.pv, newdata = data, type = "response")
  output.data <- data.frame("id" = data$id, "pred.obs" = pred.obs, "pred" = data$pred)
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(output.data$pred.obs - output.data$pred))
  E50 <- median(abs(output.data$pred.obs - output.data$pred))
  E90 <- as.numeric(quantile(abs(output.data$pred.obs - output.data$pred), probs = .9, na.rm = TRUE))
  
  ### Create dataframe for plotting
  ### Do so over the range of values pred.plot.range, if specified
  if (!is.null(pred.plot.range)){
    ### Create temporary data frame
    tmp.data <- data.frame("pred" = pred.plot.range, "pred.logit" = log(pred.plot.range/(1-pred.plot.range)))
    
    ### Calculate predicted observed values over pred.plot.range
    pred.obs <- predict(calib.model.pv, newdata = tmp.data, type = "response")
    
    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = tmp.data$pred)
    
  } else {
    
    ### Create plot data
    plot.data <- output.data
    
    ### Now set output.data to NA, 
    ### this is so we don't save it twice (it will be an element in the plot object, which will be outputted)
    output.data <- NA
    
  }
  
  ### Create plot
  if (plot == TRUE){
    plot <- ggplot2::ggplot(data = plot.data) +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
      ggplot2::ggtitle(("Pseudo-value")) + 
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  } else {
    plot <- NULL
  }
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "calib_data" = output.data)
  return(output.object)
  
}

###
### Function to estimate calibration curve using pseudo-values with inverse probability of censoring weights,
### following the binder approach. 
###
### See:
### Binder; Pseudo-observations for competing risks with covariate dependent censoring. DOI: 10.1007/s10985-013-9257-7.
### Overgaard; Pseudo-observations under covariate-dependent censoring. DOI: 10.1016/j.jspi.2019.02.003
###
est_calib_pv_ipcw_eventglm <- function(data, fit, bhaz, t, surv = NULL, nk = 4,
                                       eventglm.ipcw.type = "coxph",
                                       eventglm.ipcw.formula = NULL,
                                       eventglm.pv.method = "binder",
                                       pred.plot.range = NULL,
                                       use.logit = FALSE){
  
  # data = dat.valid
  # fit  = fit1
  # bhaz = bhaz1
  # t = t.eval
  # nk = 5
  # split.n.groups = 25
  # surv <- NULL
  # group.var = NULL
  # group.by = "ipcw"
  # ipcw.cens.max.follow = 0.5
  # ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ rms::rcs(x1,5) + rms::rcs(x2,5) + rms::rcs(x3,5)")
  # ipcw.type = "coxph"
  
  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Add predicted risk
  data$pred <- 1 - data$surv
  data$pred.logit <- log(data$pred/(1 - data$pred))
  
  ### Estimate pseudo-values
  ### NB: this uses the inverse of the outcome indicator to model censoring when estimating the weights.
  ### This can cause problems if censoring the data (i.e. stopped cox)
  pvs <- eventglm::pseudo_coxph(Surv(time, status) ~ 1, time = t, data = data, 
                                type = "cuminc",
                                formula.censoring = eventglm.ipcw.formula, 
                                ipcw.method = eventglm.pv.method)
  
  ### Add the pseudo-values to validation dataset
  data$pv <- pvs
  
  ### If pseudo-value is NA, this means all individual had had an event prior to time point t.eval
  ### By definition pseudo-values for these individuals are equal to 1, so assign these
  data$pv[is.na(data$pv)] <- 1
  
  ### Transform est.surv onto logit scale
  data$pred.logit <- log(data$pred/(1-data$pred))
  
  ### Fit the model using logit link function
  calib.model.pv <- stats::glm(pv ~ rms::rcs(pred.logit, nk), 
                               data = data, 
                               family = stats::gaussian(link = "logit"),
                               start = rep(0, nk))
  
  ###
  ### Generate predicted observed values,
  ###
  
  ### Create a dataframe with the calibration data
  pred.obs <- predict(calib.model.pv, newdata = data, type = "response")
  output.data <- data.frame("id" = data$id, "pred.obs" = pred.obs, "pred" = data$pred)
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(output.data$pred.obs - output.data$pred))
  E50 <- median(abs(output.data$pred.obs - output.data$pred))
  E90 <- as.numeric(quantile(abs(output.data$pred.obs - output.data$pred), probs = .9, na.rm = TRUE))
  
  ### Create dataframe for plotting
  ### Do so over the range of values pred.plot.range, if specified
  if (!is.null(pred.plot.range)){
    ### Create temporary data frame
    tmp.data <- data.frame("pred" = pred.plot.range, "pred.logit" = log(pred.plot.range/(1-pred.plot.range)))
    
    ### Calculate predicted observed values over pred.plot.range
    pred.obs <- predict(calib.model.pv, newdata = tmp.data, type = "response")
    
    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = tmp.data$pred)
    
  } else {
    
    ### Create plot data
    plot.data <- output.data
    
    ### Now set output.data to NA, 
    ### this is so we don't save it twice (it will be an element in the plot object, which will be outputted)
    output.data <- NA
    
  }
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle(("Pseudo-value")) + 
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "calib_data" = output.data)
  return(output.object)
  
}

#######################
#######################
### OTHER FUNCTIONS ###
#######################
#######################

###
### Function to calculate the true calibration curve for a model
### Note we must also be able to specify the true linear predictor for each individual, 
### and the cumulative bsaeline hazard at time t.eval, which are specific to each DGM
###
est_calib_true <- function(data, fit, bhaz, t, surv = NULL, true.lp, cumhaz.t, nk = 5){
  
  # data = dat.valid
  # fit = fit1
  # bhaz = bhaz1
  # t = t.eval
  # true.lp = dat.valid$true.lp
  # cumhaz.t = 0.5
  # surv = NULL
  # nk = 4
  
  ### Calculate survival probabilities from the model
  if (is.null(surv)){
    surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    surv <- as.numeric(surv)
  }
  
  ### Calculate the true survival probabilities
  true.surv <- exp(-cumhaz.t*exp(as.numeric(true.lp)))
  
  ### Combine with estimated survival probabilities
  true.dat <- data.frame("id" = data$id, 
                         "true.surv" = true.surv, "est.surv" = surv,
                         "true.pred" = 1 - true.surv, "est.pred" = 1 - surv)
  
  ### Transform est.surv onto logit scale
  true.dat$est.pred.logit <- log(true.dat$est.pred/(1-true.dat$est.pred))
  
  ### Fit the model using logit link function
  calib.model.true <- stats::glm(true.pred ~ rms::rcs(est.pred.logit, nk), 
                                 data = true.dat, 
                                 family = stats::gaussian(link = "logit"),
                                 start = rep(0,nk))
  
  ### Estimate predict observed values, add to true.dat, and select variables for output
  true.dat$pred.obs.true <- predict(calib.model.true, newdata = true.dat, type = "response")
  true.dat$Calibration <- "true"
  true.dat <- dplyr::select(true.dat, id, pred.obs.true, est.pred, Calibration) |>
    dplyr::rename(pred.obs = pred.obs.true, pred = est.pred)
  
  return(true.dat)
  
}

###
### Function to calculate the true mean calibration
###
est_calib_true_mean <- function(data, fit, bhaz, t, surv = NULL, true.lp, cumhaz.t, nk = 5){
  
  ### Calculate survival probabilities from the model
  if (is.null(surv)){
    surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    surv <- as.numeric(surv)
  }
  
  ### Calculate the true survival probabilities
  true.surv <- exp(-cumhaz.t*exp(as.numeric(true.lp)))
  
  ### Combine with estimated survival probabilities
  true.dat <- data.frame("id" = data$id, 
                         "true.surv" = true.surv, "est.surv" = surv,
                         "true.pred" = 1 - true.surv, "est.pred" = 1 - surv)
  
  ### Get ratio (mean observed risk over mean predicted risk)
  true.ratio <- mean(true.dat$true.pred)/mean(true.dat$est.pred)
  
  return(true.ratio)
  
}


###
### Write a function to add the true calibration curve to an existing calibration plot
###
add_calib_true <- function(calib.object, calib.true, legend.label = NULL){
  
  ### Extract the plot data
  plot.data <- calib.object[["plot"]]$data
  ### Add legend label
  if (is.null(legend.label)){
    plot.data$Calibration <- "Estimated"
  } else {
    plot.data$Calibration <- legend.label
  }
  
  ### Add true calibration curve
  plot.data <- rbind(plot.data, calib.true)
  
  ### Create plot
  calib.object[["plot"]] <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, color = Calibration)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle(calib.object[["plot"]]$labels$title) + 
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  
  ### Return calib.object
  return(calib.object)
  
}

###
### Write a function to add the confidence intervals to an existing calibration plot generated from the simulation
###
add_CI <- function(calib.object, calib.ci){
  
  ### Format the plot data
  plot.data <- calib.object[["plot"]]$data |>
    dplyr::mutate(plotnum = dplyr::case_when(Calibration == "Estimated" ~ 1,
                                             Calibration == "true" ~ 2),
                  colornum = dplyr::case_when(Calibration == "Estimated" ~ 1,
                                              Calibration == "true" ~ 2),
                  ltynum = 1)
  
  
  ### Format the calib.CI data
  calib.ci <- tidyr::pivot_longer(calib.ci,
                                  cols = c(pred.obs.lower, pred.obs.upper)) |>
    dplyr::mutate(name = dplyr::case_when(name == "pred.obs.lower" ~ 3,
                                          name == "pred.obs.upper" ~ 4),
                  colornum = 1,
                  ltynum = 2,
                  Calibration = "Estimated CI") |>
    dplyr::rename(plotnum = name, pred.obs = value)
  
  ### Add true calibration curve
  plot.data <- rbind(plot.data, calib.ci)
  
  ### Create plot
  calib.object[["plot"]] <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, group = plotnum, color = Calibration)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle(calib.object[["plot"]]$labels$title) + 
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  
  ### Return calib.object
  return(calib.object)
  
}


###
### Write a function to add the confidence intervals to an existing calibration plot, from one of the clinical examples
###
add_CI_ce <- function(calib.object, calib.ci){
  
  ### Format the plot data
  plot.data <- calib.object[["plot"]]$data |>
    dplyr::mutate(plotnum = 1,
                  ltynum = "Calibration")
  
  ### Format the calib.CI data
  calib.ci <- tidyr::pivot_longer(calib.ci,
                                  cols = c(pred.obs.lower, pred.obs.upper)) |>
    dplyr::mutate(name = dplyr::case_when(name == "pred.obs.lower" ~ 2,
                                          name == "pred.obs.upper" ~ 3),
                  ltynum = "CI") |>
    dplyr::rename(plotnum = name, pred.obs = value)
  
  ### Reduce to variables of interest, in case there were other varaibles in calib.ci (such as mean/median)
  plot.data <- dplyr::select(plot.data, pred, pred.obs, plotnum, ltynum)
  calib.ci <- dplyr::select(calib.ci, pred, pred.obs, plotnum, ltynum)
  
  ### Add true calibration curve
  plot.data <- rbind(plot.data, calib.ci)
  
  ### Create plot
  calib.object[["plot"]] <- ggplot2::ggplot(data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, group = plotnum, lty = ltynum), color = "red") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle(calib.object[["plot"]]$labels$title) + 
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
    ggplot2::theme(legend.position = "none")
  
  ### Return calib.object
  return(calib.object)
  
}
#######################################
#######################################
### Functions to run the simulation ###
#######################################
#######################################

###
### Write a function to run the simulation for model 1
### Model 1: Adjust for x1 only
### Note that the cumulative hazard at time t is an argument, which is used to estimate the true risks
### This might change from DGM to DGM, and hence is a separate argument. cumhaz.t should correspond to the value of t.eval,
### and the DGM used to generate the data.
### 
run_sim_large_sample <- function(dat.devel, dat.valid, t.eval, cumhaz.t, t.cens = TRUE, model, nk.fit = 5, 
                                 methods = c("ph", "pv.predrisk", "pv.binder", "ipcw.x1", "ipcw.flex")){
  
  # dat.devel = datC11
  # dat.valid = datC11
  # t.eval = 0.5
  # cumhaz.t = 0.5
  # R.boot = 20
  # plotname = "DGM1_C11"
  
  ### Assign the model formula
  if (model == 1){
    formula.rhs <- "x1"
  } else if (model == 2){
    formula.rhs <- "x1 + x2 + x3"
  } else if (model == 3){
    formula.rhs <- "rms::rcs(x1,nk.fit) + rms::rcs(x2,nk.fit) + rms::rcs(x3,nk.fit)"
  } else if (model == "supplementary"){
    formula.rhs <- "x1 + x2 + x3 + x4 + x5"
  }
  
  ### Censor validation data at time t if specified
  ### Some methods have had problems if not doing this
  if (t.cens == TRUE){
    dat.valid <- dplyr::mutate(dat.valid, 
                               status = dplyr::case_when(time > t.eval ~ 0,
                                                         TRUE ~ status),
                               time = dplyr::case_when(time > t.eval ~ t.eval,
                                                       TRUE ~ time))
  }
  
  ### Fit a model
  fit1 <- coxph(as.formula(paste("Surv(time, status) ~ ", formula.rhs, sep = "")), data = dat.devel)
  bhaz1 <- basehaz(fit1, centered = TRUE)
  
  ### Estimate true calibration curve
  calib.true <- est_calib_true(data = dat.valid, 
                               fit = fit1,
                               bhaz = bhaz1,
                               t = t.eval,
                               true.lp = dat.valid$true.lp, 
                               cumhaz.t = cumhaz.t)
  
  ### Evaluate calibration of this model using each of the functions that have been written
  ## Pseudo-value approach using eventglm
  print(paste("START", Sys.time()))
  
  ## Proportional hazards approach
  if ("ph" %in% methods){
    calib.ph <- est_calib_ph(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 5)
    print(paste("PH DONE", Sys.time()))
  }
  
  ## Pseudo-value approach, grouped by predicted risk
  if ("pv.predrisk" %in% methods){
    calib.pv <- est_calib_pv(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 5, split.n.groups = min(nrow(dat.valid)/500, 100),
                             group.by = "predrisk")
    print(paste("PV predrisk DONE", Sys.time()))
  }
  
  ## Pseudo-value approach using binder method (PV-IPCW from manuscript)
  if ("pv.binder" %in% methods){
    calib.pv.binder <- est_calib_pv_ipcw(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 5, 
                                         ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ rms::rcs(x1,5) + rms::rcs(x2,5) + rms::rcs(x3,5)"),
                                         eventglm.pv.method = "binder")
    print(paste("PV BINDER DONE", Sys.time()))
  }
  
  ## Pseudo-value approach with ipcw groupings
  ## NB: Do not confuse with PV-IPCW from the manuscript (which is what pv.binder is)
  ## This is a method we tried before discovering binder approach, which calculates pseudo-values within groups
  ## ranked by IPCWs
  # if ("pv.ipcw" %in% methods){
  #   calib.pv.ipcw <- est_calib_pv(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 5, split.n.groups = min(nrow(dat.valid)/500, 100),
  #                                 group.by = "ipcw",
  #                                 ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ rms::rcs(x1,5) + rms::rcs(x2,5) + rms::rcs(x3,5)"),
  #                                 ipcw.type = "coxph")
  #   print(paste("PV IPCW DONE", Sys.time()))
  # }
  
  ## IPCW approach
  if ("ipcw.x1" %in% methods){
    calib.ipcw.x1 <- est_calib_ipcw(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 5,
                                    ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ x1"))
    print(paste("IPCW-x1 DONE", Sys.time()))
  }
  
  ## IPCW-flexible approach
  if ("ipcw.flex" %in% methods){
    calib.ipcw.flex <- est_calib_ipcw(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 5,
                                      ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ rms::rcs(x1,5) + rms::rcs(x2,5) + rms::rcs(x3,5)"))
    print(paste("IPCW-flex DONE", Sys.time()))
  }
  
  ### Create empty list
  plot.data <- vector("list", 0)
  
  ### Add to list
  if (exists("calib.pv")){plot.data <- append(plot.data,  list(dplyr::mutate(calib.pv[["plot"]]$data, "Calibration" = "PV")))}
  if (exists("calib.pv.ipcw")){plot.data <- append(plot.data, list(dplyr::mutate(calib.pv.ipcw[["plot"]]$data, "Calibration" = "PV-ipcw")))}
  if (exists("calib.pv.binder")){plot.data <- append(plot.data, list(dplyr::mutate(calib.pv.binder[["plot"]]$data, "Calibration" = "PV-binder")))}
  if (exists("calib.ipcw.flex")){plot.data <- append(plot.data,  list(dplyr::mutate(calib.ipcw.flex[["plot"]]$data, "Calibration" = "IPCW-flex")))}
  if (exists("calib.ipcw.x1")){plot.data <- append(plot.data, list(dplyr::mutate(calib.ipcw.x1[["plot"]]$data, "Calibration" = "IPCW-x1")))}
  if (exists("calib.ph")){plot.data <- append(plot.data, list(dplyr::mutate(calib.ph[["plot"]]$data, "Calibration" = "PH")))}
  print(str(plot.data))
  ### rbind into a single dataset
  plot.data <- do.call("rbind", plot.data)
  print("PROGRESS")
  
  ### Add true calibration curve
  plot.data <- rbind(plot.data, calib.true)
  plot.data.reduced <- dplyr::group_by(plot.data, Calibration) |>
    dplyr::slice_head(n = 1000) |>
    as.data.frame()
  
  ### Create plot
  plot <- ggplot2::ggplot(data = plot.data.reduced) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, color = Calibration)) +
    ggplot2::geom_point(ggplot2::aes(x = pred, y = pred.obs), col = grDevices::rgb(0,0,0,alpha=0) , data = subset(plot.data.reduced, Calibration == "true")) +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
    ggplot2::theme(legend.position = "bottom") + 
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1))
  
  ### Add marginal density function
  plot <- ggExtra::ggMarginal(plot, type = "density", margins = "x", size = 6)
  
  ###
  ### Extract ICI, E50 and E90
  ###
  
  ### ICI
  ### Create empty list
  ICI <- vector("list", 0)
  
  ### Add to list
  if (exists("calib.pv")){ICI <- append(ICI,  list(calib.ph[["ICI"]]))}
  if (exists("calib.pv.ipcw")){ICI <- append(ICI, list(calib.pv.ipcw[["ICI"]]))}
  if (exists("calib.pv.binder")){ICI <- append(ICI, list(calib.pv.binder[["ICI"]]))}
  if (exists("calib.ipcw.flex")){ICI <- append(ICI,  list(calib.ipcw.flex[["ICI"]]))}
  if (exists("calib.ipcw.x1")){ICI <- append(ICI, list(calib.ipcw.x1[["ICI"]]))}
  if (exists("calib.ph")){ICI <- append(ICI, list(calib.ph[["ICI"]]))}
  
  ### E50
  ### Create empty list
  E50 <- vector("list", 0)
  
  ### Add to list
  if (exists("calib.pv")){E50 <- append(E50,  list(calib.ph[["E50"]]))}
  if (exists("calib.pv.ipcw")){E50 <- append(E50, list(calib.pv.ipcw[["E50"]]))}
  if (exists("calib.pv.binder")){E50 <- append(E50, list(calib.pv.binder[["E50"]]))}
  if (exists("calib.ipcw.flex")){E50 <- append(E50,  list(calib.ipcw.flex[["E50"]]))}
  if (exists("calib.ipcw.x1")){E50 <- append(E50, list(calib.ipcw.x1[["E50"]]))}
  if (exists("calib.ph")){E50 <- append(E50, list(calib.ph[["E50"]]))}
  
  ### E90
  ### Create empty list
  E90 <- vector("list", 0)
  
  ### Add to list
  if (exists("calib.pv")){E90 <- append(E90,  list(calib.ph[["E90"]]))}
  if (exists("calib.pv.ipcw")){E90 <- append(E90, list(calib.pv.ipcw[["E90"]]))}
  if (exists("calib.pv.binder")){E90 <- append(E90, list(calib.pv.binder[["E90"]]))}
  if (exists("calib.ipcw.flex")){E90 <- append(E90,  list(calib.ipcw.flex[["E90"]]))}
  if (exists("calib.ipcw.x1")){E90 <- append(E90, list(calib.ipcw.x1[["E90"]]))}
  if (exists("calib.ph")){E90 <- append(E90, list(calib.ph[["E90"]]))}
  
  ### Create output object
  output.object <- list("plot" = plot, "plotdata" = plot.data, "ICI" = ICI, "E50" = E50, "E90" = E90)
  
  return(output.object)
  
}


###
### General DGM
###
simulate_DGM_example <- function(m, 
                                 cov,
                                 lambda.y = 1, gamma.y = 1, lambda.cens = 1, gamma.cens = 1,
                                 b.x1 = 0, b.x2 = 0, b.x3 = 0, b.x2sq = 0, b.x3cu = 0,
                                 b.cens.x1 = 0, b.cens.x2 = 0, b.cens.x3 = 0, b.cens.x2sq = 0, b.cens.x3cu = 0,
                                 seed = NULL,
                                 cens = TRUE){
  
  ### set seed
  if (!is.null(set.seed)){
    set.seed(seed)
  }
  
  ### Add some non-linear terms
  cov <- dplyr::mutate(cov, 
                       x2sq = x2^2,
                       x3cu = x3^3)
  
  ### Simulate event times
  dat <- simsurv::simsurv(lambdas = lambda.y, 
                          gammas = gamma.y,
                          betas = c(x1 = b.x1, x2 = b.x2, x3 = b.x3, x2sq = b.x2sq, x3cu = b.x3cu),
                          x = cov,
                          maxt = 100) |>
    dplyr::rename(time = eventtime)
  
  ### Add a censoring time
  ### Simulate event times
  if (cens == TRUE){
    cens.times <- simsurv::simsurv(lambdas = lambda.cens, 
                                   gammas = gamma.cens,
                                   betas = c(x1 = b.cens.x1, x2 = b.cens.x2, x3 = b.cens.x3,
                                             x2sq = b.cens.x2sq, x3cu = b.cens.x3cu),
                                   x = cov)
    
    cens.times <- dplyr::rename(cens.times, cens_time = eventtime) |>
      dplyr::select(-c(id,status))
    
    ### Combine with covariates and cens.times
    dat <- cbind(dat, cov, cens.times)
    
    ### Apply censoring time to the event time data
    dat <- dplyr::mutate(dat, 
                         status = dplyr::case_when(cens_time < time ~ 0,
                                                   TRUE ~ 1),
                         time = dplyr::case_when(cens_time < time ~ cens_time,
                                                 TRUE ~ time))
    
    ### Create a variable for censored or not censored 
    ### (if an individual doesn't have an event, they are censored at end of follow up)
    ### This will be used in the IPCW method
    dat$cens_time <- dat$time
    dat$cens_indicator <- 1 - dat$status
  } else {
    ### Add cens_time to be same as event time, and always uncensored
    dat$cens_time <- dat$time
    dat$cens_indicator <- 0
    ### Combine with covariates
    dat <- cbind(dat, cov)
  }
  
  ### Calculate true linear predictor and add to dat (this will be used to calculating true risks in simulation)
  dat$true.lp <- as.numeric(as.matrix(cov[,1:5]) %*% c(b.x1, b.x2, b.x3, b.x2sq, b.x3cu))
  
  return(dat)
  
}



###          
### Function to calculate ICBI for a given method
###
get_ICBI <- function(calib.method, calib.data.in, ipcw.subset = FALSE){
  
  if (ipcw.subset == FALSE){
    
    ### Get the predicted observed for this method
    pred.obs.method <- subset(calib.data.in, Calibration == calib.method) |> dplyr::pull(pred.obs)
    
    ### Get the true pred.obs
    pred.obs.true <- subset(calib.data.in, Calibration == "true") |> dplyr::pull(pred.obs)
    
    ### Calculate mean absolute difference
    ICBI <- mean(abs(pred.obs.method - pred.obs.true))
    
    ### Calculate CB50 and CB90
    CB50 <- quantile(abs(pred.obs.method - pred.obs.true), p = 0.5)
    CB90 <- quantile(abs(pred.obs.method - pred.obs.true), p = 0.9)
    
  } else {
    
    ### Subset the data to those with id's that were included in the IPCW weighted analysis
    calib.data.in <- subset(calib.data.in, id %in% dplyr::pull(subset(calib.data.in, Calibration == "IPCW-flex"), id))
    
    ### Get the predicted observed for this method
    pred.obs.method <- subset(calib.data.in, Calibration == calib.method) |> dplyr::pull(pred.obs)
    
    ### Get the true pred.obs
    pred.obs.true <- subset(calib.data.in, Calibration == "true") |> dplyr::pull(pred.obs)
    
    ### Calculate mean absolute difference
    ICBI <- mean(abs(pred.obs.method - pred.obs.true))
    
    ### Calculate CB50 and CB90
    CB50 <- quantile(abs(pred.obs.method - pred.obs.true), p = 0.5)
    CB90 <- quantile(abs(pred.obs.method - pred.obs.true), p = 0.9)
    
  }
  
  ### Create output list
  output.object <- list("ICBI" = ICBI, "CB50" = as.numeric(CB50), "CB90" = as.numeric(CB90))
  return(output.object) 
  
}


#####################################################
### Unexported function from R package eventglm   ###
### SEE: github.com/sachsmc/eventglm              ###
#####################################################
match_cause <- function(mr, cause) {
  if (!is.null(attr(mr, "states"))) {
    states <- attr(mr, "states")
    if (is.numeric(cause)) {
      stopifnot(cause <= length(states))
      causec <- states[cause]
      causen <- cause
    } else {
      stopifnot(length(match(cause, states)) > 0)
      causen <- match(cause, states)[1]
      causec <- cause
    }
  } else {
    causen <- 1
    causec <- "dead"
  }
  list(causen = causen, causec = causec)
  
}

#####################################################
### Unexported function from R package eventglm   ###
### SEE: github.com/sachsmc/eventglm              ###
#####################################################
calc_ipcw_pos <- function(mr, time, causen, type, ipcw.method, Gi) {
  stopifnot(length(Gi) == nrow(mr))
  
  if (type == "cuminc") {
    Vi <- as.numeric(mr[, "time"] < time & mr[, "status"] == causen)
    
  } else if (type == "survival") {
    if (attr(mr, "type") != "right") {
      stop(
        "Survival estimand not available for outcome with censoring type",
        attr(mr, "type")
      )
    }
    
    Vi <-
      1 - as.numeric(mr[, "time"] < time & mr[, "status"] == causen)
    
  } else if (type == "rmean") {
    if (attr(mr, "type") == "mright") {
      Vi <-
        (time - pmin(mr[, "time"], time)) * as.numeric(mr[, "status"] == causen)
      
    } else {
      Vi <- pmin(mr[, "time"], time)
      
    }
  }
  
  Ii <- as.numeric(mr[, "time"] >= time | mr[, "status"] != 0)
  
  nn <- length(Vi)
  #theta.n <- sum(Ii * Vi / Gi) / sum(Ii / Gi)
  
  XXi <- Vi * Ii / Gi
  if (ipcw.method == "binder") {
    theta.n <- mean(Ii * Vi / Gi)
    POi <-
      theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i)
        mean(XXi[-i])))
    
    
  } else if (ipcw.method == "hajek") {
    theta.n <- sum(Ii * Vi / Gi) / sum(Ii / Gi)
    POi <- nn * theta.n - (nn - 1) *
      (sapply(1:length(XXi), function(i)
        sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))
    
    
  } else {
    stop(
      "Weighting method ",
      ipcw.method,
      " not available, options are 'binder' or 'hajek'"
    )
  }
  
  POi
}