#################
### Functions ###
#################

###
### Function to estimate survival probability for a given set of models and fitted baseline hazards,
### where each model was developed using a different imputed dataset.
### Survival probabilities are estimated using each model seperately, then combined on the cloglog scale
### NB: Baseline hazards must have been fitted using basehaz(surv.obj, centered = TRUE)
###
### Standard cox model
est_surv_mi <- function(newdata, fit.list, bhaz.list, time){

  #   newdata <- imps.valid[[1]][1:100, ]
  #   bhaz.list[[1]]$time
  #   time <- 10*365.25

  ### Get the lp
  lp <- lapply(fit.list, function(x) {predict(x, newdata = newdata, reference = "sample")})

  ### Get the linear predictor for ne wdata
  surv <- lapply(1:length(fit.list), function(x) {
    bhaz <- bhaz.list[[x]]
    as.numeric(exp(-exp(lp[[x]])*bhaz))
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
###
### Prototype 3 model
est_surv_offset_prototype3_mi <- function(newdata, fit.list, bhaz.list,
                                          means.statins, means.ah, means.smoking1, means.smoking2,
                                          time){

  #   newdata <- imps.valid[[1]][1:100, ]
  #   bhaz.list[[1]]$time
  #   time <- 10*365.25

  # newdata = pat.prototype3.baseline
  # fit.list = fit.prototype3
  # bhaz.list = bhaz.prototype3
  # means.statins = means.statins
  # means.ah = means.ah
  # means.smoking1 = means.smoking1
  # means.smoking2 = means.smoking2

  ### Get the lp
  lp <- lapply(1:length(fit.list), function(x) {

    ## Get lp
    lp.out <- predict(fit.list[[x]], newdata = newdata, reference = "sample")

    ## Adjust for offsets
    lp.out <- lp.out - means.statins[[x]] - means.ah[[x]] - means.smoking1[[x]] - means.smoking2[[x]]

    return(lp.out)
  })

  ### Get the survival probabilities
  surv <- lapply(1:length(fit.list), function(x) {
    bhaz <- bhaz.list[[x]]
    as.numeric(exp(-exp(lp[[x]])*bhaz))
  })

  ### Now get average on cloglog scale
  surv.average <- exp(-exp(rowMeans(do.call("cbind", lapply(surv, function(x) {log(-log(x))})))))

  ### Return
  return(surv.average)

}

###
### Function to adjust a risk score based on an odds ratio
###
convert_risk_odds <- function(p, OR){

  ### Get odds of risk score
  odds <- p/(1-p)

  ### Calculate new odds
  new_odds <- odds*OR

  ### Covert back onto risk scale
  new_risk <- new_odds/(1+new_odds)

  return(new_risk)
}

###
### Function to create table of risks for a given patient at various ages
###
create_table <- function(pat,
                         fit.prototype3, bhaz.prototype3,
                         means.statins, means.ah, means.smoking1, means.smoking2,
                         fit.standard, bhaz.standard,
                         or.change.total,
                         or.change.smoking,
                         or.change.sbp,
                         or.change.bmi,
                         or.change.nonhdl,
                         input_type){

  # pat
  # fit.tn = fit.tn
  # bhaz.tn = bhaz.tn
  # means.statins = means.tn[["statins"]]
  # means.ah = means.tn[["ah"]]
  # fit.standard = fit.standard
  # bhaz.standard = bhaz.standard

  ### Define decimal places
  dp <- 2

  ### Create output data.frame
  # pat.risks <- data.frame(matrix(NA, nrow = 1, ncol = 5))
  pat.risks <- vector("character", 4)

  ###
  ### Estimate a risk for this individual according to standard cox model
  ###

  ### Create new dataset for the old model
  ### TO DELETE ONCE OLD MODEL REFITTED WITH APPROPRIATE SMOKING TERMS
  # pat_original <- dplyr::mutate(pat, smoking = dplyr::case_when(smoking == "Current" ~ "Light",
  #                                                      TRUE ~ smoking))
  # risk <- 1 - est_surv_mi(newdata = pat_original,
  #                         fit.list = fit.standard,
  #                         bhaz.list = bhaz.standard,
  #                         time = 10*365.25)
  # pat.risks[1] <- paste(round(100*risk,dp), "%", sep = "")

  ###
  ### Estimate a risk using our offset model
  ###

  ### Assign treatment statuses to pat
  pat.prototype3.baseline <- dplyr::mutate(pat,
                                           offset_ah_timevar_lnHR = 0,
                                           offset_statins_timevar_lnHR = 0,
                                           offset_smoking_timevar_dummy1_lnHR = 0,
                                           offset_smoking_timevar_dummy2_lnHR = 0)
  # pat.on.statins <- dplyr::mutate(pat.off, med_status_statins = lnHR.statins)
  # pat.on.ah <- dplyr::mutate(pat.off, med_status_ah = lnHR.ah)
  # pat.on <- dplyr::mutate(pat.off, med_status_statins = lnHR.statins, med_status_ah = lnHR.ah)

  ### Estimate risk from baseline layer
  baseline.risk  <- 1 - est_surv_offset_prototype3_mi(newdata = pat.prototype3.baseline,
                                                      fit.list = fit.prototype3,
                                                      bhaz.list = bhaz.prototype3,
                                                      means.statins = means.statins,
                                                      means.ah = means.ah,
                                                      means.smoking1 = means.smoking1,
                                                      means.smoking2 = means.smoking2,
                                                      time = 10*365.25)
  pat.risks[1] <- paste(round(100*baseline.risk,dp), "%", sep = "")

  ### Adjust from interventional layer
  intervention.risk <- convert_risk_odds(baseline.risk, or.change.total)

  ### Save it
  pat.risks[2] <- paste(round(100*intervention.risk,dp), "%", sep = "")

  # ### Create predictors being on statins
  # if (pat$cholhdl_ratio > 3){
  #   risk <- 1 - est_surv_offset_mi(newdata = pat.on.statins,
  #                                  fit.list = fit.prototype3,
  #                                  bhaz.list = bhaz.prototype3,
  #                                  means.statins = means.statins,
  #                                  means.ah = means.ah,
  #                                  time = 10*365.25)
  #   pat.risks[3] <- paste(round(100*risk,dp), "%", sep = "")
  # } else {
  #   pat.risks[3] <- "Not eligible"
  # }
  #
  # ### Create predictors being on antihypertensives
  # if (pat$sbp > 135){
  #   risk <- 1 - est_surv_offset_mi(newdata = pat.on.ah,
  #                                  fit.list = fit.prototype3,
  #                                  bhaz.list = bhaz.prototype3,
  #                                  means.statins = means.statins,
  #                                  means.ah = means.ah,
  #                                  time = 10*365.25)
  #   pat.risks[4] <- paste(round(100*risk,dp), "%", sep = "")
  # } else {
  #   pat.risks[4] <- "Not eligible"
  # }
  #
  # ### Create predictors being on both treatment
  # if (pat$cholhdl_ratio > 3 & pat$sbp > 135){
  #   risk <- 1 - est_surv_offset_mi(newdata = pat.on,
  #                                  fit.list = fit.prototype3,
  #                                  bhaz.list = bhaz.prototype3,
  #                                  means.statins = means.statins,
  #                                  means.ah = means.ah,
  #                                  time = 10*365.25)
  #   pat.risks[5] <- paste(round(100*risk,dp), "%", sep = "")
  #
  # } else {
  #   pat.risks[5] <- "Not eligible"
  # }

  pat.risks[3] <- round(or.change.total, dp)
  pat.risks[4] <- round(or.change.sbp, dp)
  pat.risks[5] <- round(or.change.bmi, dp)
  pat.risks[6] <- round(or.change.nonhdl, dp)
  pat.risks[7] <- round(or.change.smoking, dp)

  ### Put into a table
  # pat.risks <- apply(round(100*pat.risks, 2), 2, function(x) paste(x, "%", sep = ""))
  pat.risks <- matrix(pat.risks, nrow = 1)
  pat.risks <- datatable(pat.risks, colnames = c("10-year risk with\nno change to treatment",
                                                 "10-year risk\nunder intervention",
                                                 "Relative risk reduction under intervention",
                                                 "Relative risk reduction due to change in SBP",
                                                 "Relative risk reduction due to change in BMI",
                                                 "Relative risk reduction due to change in Non-HDL cholesterol",
                                                 "Relative risk reduction due to change in smoking status"))

  return(pat.risks)

}
