#################
### Functions ###
#################

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
### Function to estimate survival probability for a given set of models and fitted baseline hazards,
### where each model was developed using a different imputed dataset.
### Survival probabilities are estimated using each model seperately, then combined on the cloglog scale
### NB: Baseline hazards must have been fitted using basehaz(surv.obj, centered = TRUE)
###
### Standard cox model
est_surv_mi <- function(newdata, fit_list, bhaz_list, time){
  
  #   newdata <- imps_valid[[1]][1:100, ]
  #   bhaz_list[[1]]$time
  #   time <- 10*365.25
  
  ### Get the lp
  lp <- lapply(fit_list, function(x) {predict(x, newdata = newdata, reference = "sample")})
  
  ### Get the linear predictor for ne wdata
  surv <- lapply(1:length(fit_list), function(x) {
    bhaz <- bhaz_list[[x]]
    as.numeric(exp(-exp(lp[[x]])*bhaz))
  })
  
  ### Now get average on cloglog scale
  surv_average <- exp(-exp(rowMeans(do.call("cbind", lapply(surv, function(x) {log(-log(x))})))))
  
  ### Return
  return(surv_average)
  
}

### Function to apply hazard ratio to a risk score
convert_risk_HR <- function(p, HR){
  
  ### Get survival probability
  surv <- 1 - p
  
  ### Get survival probability onto the hazard scale by taking logarithm
  log_surv <- log(surv)
  
  ### Apply hazard ratios
  log_surv_adjusted <- log_surv*HR
  
  ### Calculate survival probability (note, we do not take the negative, 
  ### as we did not do this when converted back onto the hazard scale)
  new_surv <- exp(log_surv_adjusted)
  
  ### Convert back onto risk scale
  new_risk <- 1 - new_surv
  
  return(new_risk)
}


###
### Function to create table of risks for a given patient at various ages
###
create_table <- function(pat,
                         fit_pui_list, bhaz_pui_list, means_pui_list,
                         HR_change_total,
                         input_type){
  
  # pat
  # fit_tn = fit_tn
  # bhaz_tn = bhaz_tn
  # means_statins = means_tn[["statins"]]
  # means_ah = means_tn[["ah"]]
  # fit_standard = fit_standard
  # bhaz_standard = bhaz_standard
  
  ### Define decimal places
  dp <- 2
  
  ### Create output data.frame
  # pat_risks <- data.frame(matrix(NA, nrow = 1, ncol = 5))
  pat_risks <- vector("character", 2)
  
  ###
  ### Estimate a risk using our offset model
  ###
  
  ### Assign treatment statuses to pat
  pat_pui_baseline <- dplyr::mutate(pat,
                                    offset_ah_timevar_lnHR = 0,
                                    offset_statins_timevar_lnHR = 0,
                                    offset_smoking_timevar_dummy1_lnHR = 0,
                                    offset_smoking_timevar_dummy2_lnHR = 0)
  
  ### Restrict baseline hazards to chosen time point
  bhaz_pui_list_y <- lapply(bhaz_pui_list, function(x){x[[10]]})
  
  ### Estimate risk from baseline layer
  baseline_risk  <- 1 - est_surv_offset_prototype3_redacted_mi(newdata = pat_pui_baseline,
                                                      fit_list = fit_pui_list,
                                                      bhaz_list = bhaz_pui_list_y,
                                                      offset_means_list = means_pui_list)
  pat_risks[1] <- paste(round(100*baseline_risk,dp), "%", sep = "")
  
  ### Adjust from interventional layer
  intervention_risk <- convert_risk_HR(baseline_risk, HR_change_total)
  
  ### Save it
  pat_risks[2] <- paste(round(100*intervention_risk,dp), "%", sep = "")
  
  ### Put into a table
  # pat_risks <- apply(round(100*pat_risks, 2), 2, function(x) paste(x, "%", sep = ""))
  pat_risks <- matrix(pat_risks, nrow = 1)
  pat_risks <- datatable(pat_risks, colnames = c("10-year risk with\nno change to treatment",
                                                 "10-year risk\nunder intervention"))
  
  return(pat_risks)
  
}
