###
### This file contains functions used to estimate calibration using IP weighting.
### These are primarily used as part of the sensitivity analyses to implement the K&VG
### artificial censoring approach to model validation.
###
### These functions have been taken and amended from the GitHub repo for the paper: 
### "Calibration curves for survival models, the role of covariate dependent censoring", currently under review.
### After publication, the most up to date version of these functions will be made available there:
### https://github.com/manchester-predictive-healthcare-group/CHI-CHARIOT/tree/main/project-4-calibration-survival
###

###
### Function to estimate calibration plots using BLR-IPCW approach
###
est_calib_ipcw_new <- function(data, 
                               fit, 
                               bhaz, 
                               t, 
                               surv = NULL, 
                               cens.max.follow = NULL, 
                               nk = 4, 
                               ipcw.formula, 
                               type = "coxph", 
                               pred.plot.range = NULL,
                               pred.cap.quantile = .99,
                               plot = TRUE,
                               weights_in = NULL){
  
  ### 1. Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, time = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### 2. Assign a variable for whether individuals have had an event by time of interest
  data <- dplyr::mutate(data,
                        status_time = dplyr::case_when(time <= t ~ status,
                                                       time > t ~ 0))
  
  ### 3. Create risks
  data$pred <- 1 - data$surv
  # Take predictions away from 0/1
  data$pred[data$pred == 1] <- 0.9999
  data$pred[data$pred == 0] <- 0.0001
  # Put onto logit scale
  data$pred.logit <- log(data$pred/(1 - data$pred))
  
  ### 4. Add weights
  if (is.null(weights_in)){
    ### Estimate weights
    weights <- est_ipcw(data = data, t = t, cens.max.follow = cens.max.follow, cens.formula = ipcw.formula, type = type)
    
    ### Merge with data
    data.model <- merge(data, weights, by.x = "id", by.y = "id")
    
  } else {
    
    ### Add weights
    ### Note, they are currently not stabilised, just rushing through this to see if it works
    data.model <- dplyr::mutate(data, ipcw.stab = weights_in)
  }
  
  ### 5. Fit weighted calibration model and generate predicted observed
  ## Fit model
  rcs.model.stab <- suppressWarnings(glm(status_time ~ rms::rcs(pred.logit, nk),
                                         family = binomial(link = "logit"),
                                         data = data.model,
                                         weights = data.model[, "ipcw.stab"]))
  ## Suppress warnings due to using weights in binomial glm results in 
  ## "Warning message: In eval(family$initialize) : non-integer #successes in a binomial glm!"
  
  ### 6. Generate predicted observed values and create output data
  data$pred.obs <- predict(rcs.model.stab, newdata = data, type = "response")
  if ("id" %in% colnames(data)){
    output.data <- data.frame("id" = data$id, "pred.obs" = data$pred.obs, "pred" = data$pred)
  } else {
    output.data <- data.frame("pred.obs" = data$pred.obs, "pred" = data$pred)
  }
  
  ### 7. Calculate ICI, E50 and E90
  ICI <- mean(abs(data$pred.obs - data$pred))
  E50 <- median(abs(data$pred.obs - data$pred))
  E90 <- as.numeric(quantile(abs(data$pred.obs - data$pred), probs = .9, na.rm = TRUE))
  
  ### 8. Create Plotting Data
  ### Cap plot range at xth percentile of predicted risks
  pred_cap <- as.numeric(stats::quantile(data$pred, probs = pred.cap.quantile, na.rm = TRUE))
  
  if (!is.null(pred.plot.range)) {
    plot_df <- data.frame(
      pred       = pred.plot.range,
      pred.logit = log(
        pred.plot.range /
          (1 - pred.plot.range)
      )
    )
    plot_df$pred.obs <- as.numeric(
      predict(rcs.model.stab, newdata = plot_df, type = "response")
    )
  } else {
    plot_df <- data.frame(
      pred = seq(min(data$pred), pred_cap, length.out = 100)
    ) |>
      dplyr::mutate(pred.logit = log(pred / (1 - pred)))
    plot_df$pred.obs <- as.numeric(
      predict(rcs.model.stab, newdata = plot_df, type = "response")
    )
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
      ggplot2::ggtitle("BLR in Weighted Uncensored Cohort") +
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
  
  return(output.object)
  
}


###
### Function to estimate calibration plots using smoothed-HT approach
###
est_calib_smoothed_ht <- function(data, 
                                  fit, 
                                  bhaz, 
                                  surv = NULL,
                                  t, 
                                  weights_in, 
                                  nk = 5, 
                                  pred.plot.range = NULL,
                                  pred.cap.quantile = .99,
                                  plot = TRUE) {
  
  ### 1. Calculate predicted risks at time t
  if (is.null(surv)){
    data$pred <- 1 - est_surv(newdata = data, fit = fit, bhaz = bhaz, time = t)
  } else {
    data$pred <- 1 - surv
  }
  
  # Take predictions away from 0/1
  data$pred[data$pred == 1] <- 0.9999
  data$pred[data$pred == 0] <- 0.0001
  
  ### 2. Create the HT-transformed outcome
  # I(T <= t & status == 1) is the numerator indicator
  # We multiply by the weights (IPCW * IPACW)
  # This is the 'pseudo-binary' outcome for the regression
  data$status_time <- ifelse(data$time <= t & data$status == 1, 1, 0)
  data$ht_outcome  <- data$status_time * weights_in
  
  ### 3. Fit the Smooth Calibration Model
  # We use a GAM with a spline (s). 
  # Note: family = gaussian is used because ht_outcome can be > 1.
  # The smoothing ensures we get the local expectation E[HT_outcome | Pred]
  calib_model <- mgcv::gam(ht_outcome ~ s(pred, k = nk), data = data)
  
  ### 4. Calculate predicted-observed values for metrics
  data$pred.obs <- as.numeric(predict(calib_model, newdata = data))
  
  ### 5. Calculate Metrics
  ICI <- mean(abs(data$pred.obs - data$pred))
  E50 <- median(abs(data$pred.obs - data$pred))
  E90 <- as.numeric(quantile(abs(data$pred.obs - data$pred), probs = .9, na.rm = TRUE))
  
  ### 6. Create output data
  if ("id" %in% colnames(data)){
    output.data <- data.frame("id" = data$id, "pred.obs" = data$pred.obs, "pred" = data$pred)
  } else {
    output.data <- data.frame("pred.obs" = data$pred.obs, "pred" = data$pred)
  }
  
  ### 7. Create Plotting Data
  ### Cap plot range at xth percentile of predicted risks
  pred_cap <- as.numeric(stats::quantile(data$pred, probs = pred.cap.quantile, na.rm = TRUE))
  
  if (!is.null(pred.plot.range)) {
    plot_df <- data.frame(pred = pred.plot.range)
    plot_df$pred.obs <- as.numeric(predict(calib_model, newdata = plot_df))
  } else {
    plot_df <- data.frame(pred = seq(min(data$pred), pred_cap, length.out = 100))
    plot_df$pred.obs <- as.numeric(predict(calib_model, newdata = plot_df))
  }
  
  ### 8. Generate Plot
  if (plot) {
    
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
      ggplot2::labs(
        title = paste0("HT-Weighted Regression"),
        x     = "Predicted Risk",
        y     = "Predicted-Observed Risk"
      ) +
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
  
  ### 9. Create and return output object
  return(list(
    "plot"       = plot_obj,
    "ICI"        = ICI,
    "E50"        = E50,
    "E90"        = E90,
    "calib_data" = output.data
  ))
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
                              weights_in,
                              ipcw.formula = NULL,
                              eventglm.pv.method = "binder",
                              pred.plot.range = NULL, 
                              pred.cap.quantile = .99,
                              plot = TRUE,
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
  
  ### 1. Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, time = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### 2. Add predicted risk
  data$pred <- 1 - data$surv
  # Take predictions away from 0/1
  data$pred[data$pred == 1] <- 0.9999
  data$pred[data$pred == 0] <- 0.0001
  # Put onto logit scale
  data$pred.logit <- log(data$pred/(1 - data$pred))
  
  ### 3. Estimate pseudo-values using binder approach from eventglm
  ### Do this manually, so we can alter the function for estimating the weights,
  ### so its more computationally efficient.
  ###
  ### The following code, is therefore predominately taken from the eventglm::pseudo_coxph function with a few edits
  ### see: github.com/sachsmc/eventglm
  ###
  
  ###
  ### Assign variable names to match those for the input parameters for eventglm::psudo_coxph
  ###
  formula = "survival::Surv(time, status) ~ 1" #NB the formula argument is irrelevant as we are just estimating the pseudo-values
  time = t
  cause = 1
  type = "cuminc"
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
  
  ### Weights put in manually
  ### (this should be survival probabilites, so going to use 1/weights)
  Gi.manual <- 1/weights_in
  
  ### Calculate pseudo-values using unexported function from eventglm, calc_ipcw_pos
  data$pv <- calc_ipcw_pos(mr, time, causen, type, ipcw.method, Gi.manual)
  
  ### If pseudo-value is NA, this means all individual had had an event prior to time point t.eval
  ### By definition pseudo-values for these individuals are equal to 1, so assign these
  data$pv[is.na(data$pv)] <- 1
  
  ### 4. Fit the calibration model
  calib.model.pv <- stats::glm(pv ~ rms::rcs(pred, nk), 
                               data = data, 
                               family = stats::gaussian(link = "identity"),
                               start = rep(0, nk))
  
  ### 5. Generate predicted observed values
  data$pred.obs <- predict(calib.model.pv, newdata = data, type = "response")
  
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
      pred.logit = log(
        pred.plot.range /
          (1 - pred.plot.range)
      )
    )
    plot_df$pred.obs <- as.numeric(
      predict(calib.model.pv, newdata = plot_df, type = "response")
    )
  } else {
    plot_df <- data.frame(
      pred = seq(min(data$pred), pred_cap, length.out = 100)
    ) |>
      dplyr::mutate(pred.logit = log(pred / (1 - pred)))
    plot_df$pred.obs <- as.numeric(
      predict(calib.model.pv, newdata = plot_df, type = "response")
    )
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
      ggplot2::ggtitle("Pseudo-value") +
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
  # Indicators
  Vi <- as.numeric(mr[, "time"] < time & mr[, "status"] == causen)
  Ii <- as.numeric(mr[, "time"] >= time | mr[, "status"] != 0)
  nn <- length(Vi)
  
  # Weighted values (Numerator terms)
  XXi <- (Vi * Ii) / Gi
  # Denominator terms
  Di <- Ii / Gi
  
  if (ipcw.method == "binder") {
    sum_XX <- sum(XXi)
    # Optimized O(n) Jackknife for Mean
    # Mean(-i) = (Total_Sum - Value_i) / (n - 1)
    theta_minus_i <- (sum_XX - XXi) / (nn - 1)
    theta_n <- mean(XXi)
    POi <- theta_n + (nn - 1) * (theta_n - theta_minus_i)
    
  } else if (ipcw.method == "hajek") {
    sum_XX <- sum(XXi)
    sum_D  <- sum(Di)
    # Optimized O(n) Jackknife for Hajek (Ratio of means)
    # Ratio(-i) = (Sum_Num - Num_i) / (Sum_Denom - Denom_i)
    theta_minus_i <- (sum_XX - XXi) / (sum_D - Di)
    theta_n <- sum_XX / sum_D
    POi <- nn * theta_n - (nn - 1) * theta_minus_i
  }
  
  return(POi)
}