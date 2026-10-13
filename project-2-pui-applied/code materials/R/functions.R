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
### Function to estimate survival probability for a given set of models and fitted baseline hazards,
### where each model was developed using a different imputed dataset.
### Survival probabilities are estimated using each model seperately, then combined on the cloglog scale
### NB: Baseline hazards must have been fitted using basehaz(surv.obj, centered = TRUE)
### Only used for non-PUI models
est_surv_mi <- function(newdata, fit.list, bhaz.list, time){
  
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
### Function to estimate calibration using proportional hazards approach 
### (Austin et al., graphical calibration curves. DOI: 10.1002/sim.8570)
###
est_calib_ph <- function(data, 
                         fit, 
                         bhaz, 
                         time, 
                         surv = NULL, 
                         select_knots = TRUE,  
                         knot_candidates = c(3, 4, 5),  
                         alpha_lrt = 0.05,
                         weights = NULL,
                         pred_plot_range = NULL,
                         pred_cap_quantile_lower = .01,
                         pred_cap_quantile_upper = .99,
                         plot = TRUE,
                         remove_low_risk = FALSE){
  
  ### Calculate survival probability
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv_offset_prototype3(newdata = data, fit = fit, bhaz = bhaz, time = time))
  } else {
    data$surv <- as.numeric(surv)
  }
  # Take predictions away from 0/1
  data$surv[data$surv == 1] <- 0.9999
  data$surv[data$surv == 0] <- 0.0001
  
  ### 2. Add predicted risk
  data$pred <- 1 - data$surv
  
  # Put surv onto cloglog scale
  data$cloglog <- log(-log(data$surv))
  
  ### ------------------------------------------------------------------
  ### NEW: Data-driven knot selection for the RCS calibration model
  ###
  ### Follows the procedure described by van Calster's group (e.g. the
  ### clustered/meta-analytic calibration paper, and the CalibrationCurves
  ### package docs): fit RCS calibration models with nk = 3, 4, 5, run a
  ### likelihood ratio test between every pair of candidates, 
  ### and pick the smallest nk that is not
  ### significantly outperformed by any larger candidate.
  ###
  ### Note: RCS bases with different nk use different quantile-based knot locations, so the
  ### candidate models are not strictly likelihood ratio tests for nested designs.
  ### ------------------------------------------------------------------
  if (select_knots) {
    
    knot_candidates <- sort(unique(knot_candidates))
    
    # Fit one Cox calibration model per candidate nk.
    # Wrapped in tryCatch so that a convergence failure at one nk (e.g.
    # from separation or a poorly populated region of cloglog) doesn't
    # kill the whole selection step -- we just drop that candidate.
    fits_by_nk <- list()
    # cycle through knots
    for (k in knot_candidates) {
      fit_k <- tryCatch(
        survival::coxph(
          survival::Surv(cvd_time, cvd_indicator) ~ rms::rcs(cloglog, k),
          data = data, weights = weights),
        error = function(e) NULL
      )
      if (!is.null(fit_k) && !anyNA(coef(fit_k))) {
        fits_by_nk[[as.character(k)]] <- fit_k
      } else {
        warning(sprintf("nk = %d failed to fit or has an NA coefficient; dropping from candidate set.", k))
      }
    }
    
    # List nk's that converged
    available_nk <- as.integer(names(fits_by_nk))
    
    if (length(available_nk) == 0) {
      stop("No candidate nk value produced a valid calibration model fit.")
    } else if (length(available_nk) == 1) {
      # Only one candidate survived -- nothing to compare, just use it
      nk <- available_nk
      message(sprintf("Only nk = %d converged; using it without an LRT comparison.", nk))
    } else {
      
      # All pairwise comparisons among the surviving candidates.
      pairs <- combn(available_nk, 2, simplify = FALSE)
      
      # Create table of comparisons
      lrt_table <- do.call(rbind, lapply(pairs, function(p) {
        k_small <- min(p)
        k_large <- max(p)
        fit_small <- fits_by_nk[[as.character(k_small)]]
        fit_large <- fits_by_nk[[as.character(k_large)]]
        
        # loglik[2] = model (partial) log-likelihood at convergence;
        lr_stat <- 2 * (fit_large$loglik[2] - fit_small$loglik[2])
        df_diff <- length(coef(fit_large)) - length(coef(fit_small))
        
        # Guard against a negative LR stat (can happen with
        # non-nested bases if the "larger" model happens to fit worse).
        # A negative or zero stat just means the small model is at least
        # as good, so treat it as trivially non-significant.
        p_val <- if (lr_stat <= 0) 1 else pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
        
        # Create table
        data.frame(k_small = k_small, k_large = k_large,
                   lr_stat = lr_stat, df = df_diff, p_value = p_val)
      }))
      
      # Walk candidates from fewest to most knots. Accept the first one
      # that is not significantly beaten by any larger candidate.
      chosen_nk <- max(available_nk)  # fallback if every comparison is significant
      for (k in sort(available_nk)) {
        rows_where_k_is_smaller <- lrt_table$k_small == k
        significantly_beaten <- any(lrt_table$p_value[rows_where_k_is_smaller] < alpha_lrt)
        if (!significantly_beaten) {
          chosen_nk <- k
          break
        }
      }
      
      # Define nk
      nk <- chosen_nk
      message(sprintf("Selected nk = %d via all-pairs LRT (alpha = %.2f).", nk, alpha_lrt))
      print(lrt_table)
    }
  } else {
    ### Use 5 by default
    nk <- 5
  }
  
  ### Fit calibration model (uses the fixed `nk` argument, or the value
  ### chosen above if select_knots = TRUE)
  if (select_knots){
    fit_calib <- fits_by_nk[[as.character(nk)]]
  } else {
    fit_calib <- survival::coxph(survival::Surv(cvd_time, cvd_indicator) ~ rms::rcs(cloglog, nk), data = data, weights = weights)
  }
  
  bhaz_calib <- survival::basehaz(fit_calib, centered = TRUE)
  
  ###
  ### Generate predicted observed values
  
  ### First we do this over range of predicted probabilities in the validation dataset
  data$pred.obs <- 1 - est_surv(newdata = data, fit = fit_calib, bhaz = bhaz_calib, time = time)
  
  ### 6. Calculate ICI, E50 and E90
  ICI <- mean(abs(data$pred.obs - data$pred), na.rm = TRUE)
  E50 <- median(abs(data$pred.obs - data$pred), na.rm = TRUE)
  E90 <- as.numeric(quantile(abs(data$pred.obs - data$pred), probs = .9, na.rm = TRUE))
  
  ### 7. Create output data
  if ("id" %in% colnames(data)){
    output_data <- data.frame("id" = data$id, "pred.obs" = data$pred.obs, "pred" = data$pred)
  } else {
    output_data <- data.frame("pred.obs" = data$pred.obs, "pred" = data$pred)
  }
  
  ### 8. Create Plotting Data
  
  if (!is.null(pred_plot_range)) {
    ### Create df
    plot_df <- data.frame(
      pred       = pred_plot_range,
      cloglog = log(-log(1-pred_plot_range))
    )
    plot_df$pred.obs <- 
      1 - as.numeric(est_surv(newdata = plot_df, fit = fit_calib, bhaz = bhaz_calib, time = time))
  } else {
    ### Cap plot range at xth percentile of predicted risks
    pred_cap_lower <- as.numeric(stats::quantile(data$pred, probs = pred_cap_quantile_lower, na.rm = TRUE))
    pred_cap_upper <- as.numeric(stats::quantile(data$pred, probs = pred_cap_quantile_upper, na.rm = TRUE))
    ### create df
    plot_df <- data.frame(
      pred = seq(pred_cap_lower, pred_cap_upper, length.out = 100)
    ) |>
      dplyr::mutate(cloglog = log(-log(1 - pred)))
    plot_df$pred.obs <- 
      1 - as.numeric(est_surv(newdata = plot_df, fit = fit_calib, bhaz = bhaz_calib, time = time))
  }
  
  ### 9. Generate Plot
  if (plot == TRUE) {
    
    ### Define axis max
    if (!is.null(pred_plot_range)){
      axis_max <- max(pred_plot_range) * 1.02
    } else {
      axis_max <- pred_cap_upper * 1.02
    }
    
    ### Subsample raw data for the marginal histogram (max 10,000 individuals)
    data_hist <- data[sample(nrow(data), min(10000L, nrow(data))), ]
    
    ### Create plot object
    plot_obj_base <- ggplot2::ggplot() +
      ### Invisible geom_point from subsampled raw data ? required by ggMarginal;
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
  output_object <- list("plot" = plot_obj,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "nk" = nk,  
                        "calib_data" = output_data,
                        "calib_data_plot" = plot_df)
  
  return(output_object)
  
}

###
### Write a function to estimate calibration plots for a dataset using Kaplan-Meier vs mean risk within subgroups.
### Calibration is assessed in one validation dataset, but the survival probabilities are estimated by averaging the predicted
### risks from each developed model (multiple imputation) using Rubins rules.
### data is the data in which calibration will be assessed
### fit.list is list of fitted Cox-models for which calibration will be assessed (each from a different imputed dataset)
### bhaz.list is corrresponding list of baseline hazards (centered)
###
### fit.type = "standard" is for a non-cox model
### fit.type = "prototype3" is for the final model used in CHARIOT (initial risk estimation layer)
est_calib_plot_group <- function(data, 
                                 fit, 
                                 bhaz, 
                                 time, 
                                 n.groups, 
                                 type = c("imp", "single"), 
                                 fit.type = c("standard","prototype3"),
                                 surv = NULL,
                                 single.CI = FALSE){
  
  ### Set args
  type <- match.arg(type)
  fit.type <- match.arg(fit.type)
  
  ### Get the survival probabilities
  if (is.null(surv)){
    if (fit.type == "standard"){
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
  
  ### Create output_object
  output_object <- list("plot.data" = plot.data)
  
  return(output_object)
  
}


#################################
### FUNCTIONS FOR PROTOTYPE 3 ###
#################################

###
### Prototype3 is the final pui model.
###

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
### Function to calculate survival probability for a single redacted model and baseline hazard.
###
### - fit is a single model
### - bhaz is a single cumulative baseline hazard at time t
### - offset_means is a vector of the mean offset values used for model development (on log-hazard scale)
est_surv_offset_prototype3_redacted <- function(newdata, fit, bhaz, offset_means){
  
  ### Calculate lp
  lp <- predict(fit, newdata = newdata, reference = "sample") - sum(offset_means)
  
  ### Calculate survival prob
  surv <- as.numeric(exp(-exp(lp) * bhaz))
  
  ### Return
  return(surv)
  
}

###
### Function to calculate survival probability for a list of models and baseline hazards, estimated from
### different multiple imputation chains, which have been redacted. The risk scores are calculated, and then combined on
### log-cumulative hazard scale.
###
### - fit is a list of models
### - bhaz is a list of cumulative baseline hazard at time t
### - offset_means is is a list of the mean offset values used for model development (on log-hazard scale)
est_surv_offset_prototype3_redacted_mi <- function(newdata, fit_list, bhaz_list, offset_means_list){
  
  if (length(fit_list) != length(bhaz_list) || length(fit_list) != length(offset_means_list)){
    stop("inputs should be same length")
  }
  
  ### Calculate surv probs
  surv_list <- lapply(seq_along(fit_list), function(x){
    est_surv_offset_prototype3_redacted(newdata = newdata, 
                                        fit = fit_list[[x]], 
                                        bhaz = bhaz_list[[x]], 
                                        offset_means = offset_means_list[[x]])
  })
  
  ### Put onto log cumulative hazard scale 
  surv_loglog_list <- lapply(surv_list, function(x){log(-log(x))})
  
  ### cbind and take average
  surv_loglog <- do.call("cbind", surv_loglog_list) 
  surv_loglog_av <- rowMeans(surv_loglog)
  
  ### Turn back to surv prob
  surv_av <- exp(-exp(surv_loglog_av))
  
  return(surv_av)
  
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


###########################################################
### Functions to estimate counterfactual survival times ###
###########################################################

### Function to extract cumulative baseline hazard at specified times
get_bhaz <- function(times, bhaz) {
  
  # For each requested time, find the index of the largest value
  # of bhaz$time that is less than or equal to that time
  idx <- findInterval(times, bhaz$time)
  
  # Initialise all cumulative hazard values as zero
  out <- numeric(length(times))
  
  # Identify times with a corresponding value in bhaz
  # idx = 0 occurs when a time precedes the first value of bhaz$time
  has_hazard <- idx > 0L
  
  # Extract cumulative hazard for valid indices
  out[has_hazard] <- bhaz$hazard[idx[has_hazard]]
  
  # Return cumulative hazard at each requested time
  return(out)
}

### Function to find the first time at which the cumulative baseline
### hazard is greater than or equal to a specified cumulative hazard
get_time_from_hazard <- function(cumhaz, bhaz) {
  
  # Find the first value of bhaz$hazard that is greater than
  # or equal to each supplied cumulative hazard
  idx <- findInterval(
    cumhaz,
    bhaz$hazard,
    left.open = TRUE
  ) + 1L
  
  # Initialise the output as missing
  out <- rep(NA_real_, length(cumhaz))
  
  # Identify cumulative hazards within the range of bhaz$hazard
  within_range <- idx <= nrow(bhaz)
  
  # Return the corresponding survival times for valid indices
  out[within_range] <- bhaz$time[idx[within_range]]
  
  # Return survival times
  return(out)
}

### Function to estimate counterfactual survival times.
estimate_cf_surv_times <- function(data_valid, cohort_split_times, bhaz, type) {
  
  ### Merge validation dataset with cohort.split.times
  if (type == "boot") {
    
    data_valid <- dplyr::left_join(
      dplyr::select(
        data_valid,
        c("patid", "patid_boot", "smoking")
      ),
      cohort_split_times,
      by = dplyr::join_by("patid"),
      relationship = "many-to-many"
    )
    
    # Use the bootstrap-specific identifier
    id_var <- "patid_boot"
    
  } else if (type == "orig") {
    
    data_valid <- dplyr::left_join(
      dplyr::select(
        data_valid,
        c("patid", "smoking")
      ),
      cohort_split_times,
      by = dplyr::join_by("patid")
    )
    
    # Use the original patient identifier
    id_var <- "patid"
    
  } else {
    
    stop("type must be either 'boot' or 'orig'.")
    
  }
  
  ###
  ### Create dummy variables for smoking at baseline and during follow-up
  ###
  data_valid <- create_smoking_dummies(data_valid)
  
  ### Read treatment-effect log hazard ratios
  lnHR_statins <- readRDS("data/p5/offsets_total_lnHR_statins.rds")
  lnHR_ah <- readRDS("data/p5/offsets_total_lnHR_ah.rds")
  lnHR_smoking_dummy1_total <- readRDS("data/p5/offsets_total_lnHR_smok_init.rds")
  lnHR_smoking_dummy2_total <- readRDS("data/p5/offsets_total_lnHR_smok_cess.rds")
  
  ### Check that variables used by findInterval() are ordered
  if (is.unsorted(bhaz$time)) {
    stop("bhaz$time must be sorted in non-decreasing order.")
  }
  
  if (is.unsorted(bhaz$hazard)) {
    stop("bhaz$hazard must be sorted in non-decreasing order.")
  }
  
  ### Create treatment-effect offsets
  data_valid <- dplyr::mutate(
    data_valid,
    
    offset_statins_timevar_lnHR =
      lnHR_statins * med_status_adj_statins,
    
    offset_ah_timevar_lnHR =
      lnHR_ah * med_status_adj_ah,
    
    offset_smoking_timevar_dummy1_lnHR =
      lnHR_smoking_dummy1_total *
      med_status_adj_smoking_dummy1,
    
    offset_smoking_timevar_dummy2_lnHR =
      lnHR_smoking_dummy2_total *
      med_status_adj_smoking_dummy2
  )
  
  ### Combine the four log hazard-ratio offsets
  ### exp(a + b + c + d) =
  ### exp(a) * exp(b) * exp(c) * exp(d)
  data_valid <- dplyr::mutate(
    data_valid,
    
    offset_timevar_lnHR =
      offset_statins_timevar_lnHR +
      offset_ah_timevar_lnHR +
      offset_smoking_timevar_dummy1_lnHR +
      offset_smoking_timevar_dummy2_lnHR
  )
  
  print("START")
  
  ### Start timing the counterfactual survival-time calculation
  time.in <- Sys.time()
  
  ### Calculate quantities for all follow-up intervals simultaneously
  data_valid <- dplyr::mutate(
    data_valid,
    
    # Obtain cumulative baseline hazard at the start of each interval
    cumhaz_tstart = get_bhaz(tstart, bhaz),
    
    # Obtain cumulative baseline hazard at the end of each interval
    cumhaz_tstop = get_bhaz(cvd_time, bhaz),
    
    # Convert the combined log hazard-ratio offset to a hazard ratio
    hazard_ratio = exp(offset_timevar_lnHR),
    
    # Calculate the adjusted cumulative-hazard contribution
    # from each follow-up interval
    interval_cumhaz_adj =
      (cumhaz_tstop - cumhaz_tstart) * hazard_ratio
  )
  
  ### Sum interval-level adjusted cumulative hazards within each individual
  cf_surv_times <- data_valid |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(id_var))
    ) |>
    dplyr::summarise(
      
      # Calculate total adjusted cumulative hazard
      cumhaz_adj = sum(interval_cumhaz_adj),
      
      # Retain the observed CVD indicator
      cvd_indicator_cf = max(cvd_indicator),
      
      .groups = "drop"
    )
  
  ### Convert adjusted cumulative hazards to counterfactual survival times
  ### for all individuals simultaneously
  cf_surv_times$cvd_time_cf <- get_time_from_hazard(
    cf_surv_times$cumhaz_adj,
    bhaz
  )
  
  ### Identify cumulative hazards exceeding the available hazard range
  beyond_range <- is.na(cf_surv_times$cvd_time_cf)
  
  ### Set their counterfactual survival time to maximum follow-up
  cf_surv_times$cvd_time_cf[beyond_range] <- max(bhaz$time)
  
  ### Treat these observations as censored
  cf_surv_times$cvd_indicator_cf[beyond_range] <- 0
  
  ### Retain the same output variables as the original function
  cf_surv_times <- cf_surv_times |>
    dplyr::select(
      dplyr::all_of(id_var),
      cvd_time_cf,
      cvd_indicator_cf
    ) |>
    dplyr::arrange(
      dplyr::across(dplyr::all_of(id_var))
    ) |>
    as.data.frame()
  
  ### Finish timing
  time.out <- Sys.time()
  print(time.out - time.in)
  
  warnings()
  
  return(cf_surv_times)
}
