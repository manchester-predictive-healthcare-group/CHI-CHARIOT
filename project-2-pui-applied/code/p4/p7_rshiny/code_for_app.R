###
### This will be the code that gets put into the rshiny app
###
library(shiny)
library(survival)
library(DT)
library(rms)
### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/code/p4/p7_rshiny")
getwd()

### Source functions
source("data/functions_shiny.R")

library(survival)

### Read in pat file
pat <- readRDS("data/fake.pat.rds")

    fit.prototype3 <- readRDS("data/fit.prototype3.female.rds")
    means.prototype3 <- readRDS("data/means.prototype3.female.rds")
    fit.standard <- readRDS("data/fit.standard.female.rds")
    bhaz.prototype3 <- readRDS("data/bhaz.prototype3.female.rds")
    bhaz.standard <- readRDS("data/bhaz.standard.female.rds")
  
  # ### Create pat file
  # pat$age <- input$age
  # pat$ethnicity[1] <- input$ethnicity
  # pat$bmi <- input$bmi
  # pat$sbp_unexposed <- input$sbp_unexposed
  # pat$cholhdl_ratio_unexposed <- input$cholhdl_ratio_unexposed
  # pat$smoking[1] <- input$smoking
  # pat$cortico[1] <- input$cortico
  # pat$antipsy[1] <- input$antipsy
  # pat$diabetes[1] <- input$diabetes
  # pat$hypertension[1] <- input$hypertension
  # pat$af[1] <- input$af
  # pat$ckd[1] <- input$ckd
  # pat$copd[1] <- input$copd
  # pat$fhcvd[1] <- input$fhcvd
  # pat$smi[1] <- input$smi
  # pat$migraine[1] <- input$migraine
  # pat$oral_cancer[1] <- input$oral_cancer
  # pat$brain_cancer[1] <- input$brain_cancer
  # pat$lung_cancer[1] <- input$lung_cancer
  # pat$blood_cancer[1] <- input$blood_cancer
  # pat$int_dis[1] <- input$int_dis
  # pat$downs[1] <- input$downs
  # pat$pre_eclampsia[1] <- input$pre_eclampsia
  # pat$postnatal_depression[1] <- input$postnatal_depression
  # pat$impotence[1] <- input$impotence
  # pat$IMD[1] <- input$IMD
  # 
    
    ###
    ### Estimate the relative risk change
    ###
    pat$smoking[1] <- "Current"
    
    ### First get values at baseline
    baseline.smoking <- pat$smoking
    baseline.sbp <- pat$sbp
    baseline.bmi <- pat$bmi
    baseline.nonhdl <- pat$nonhdl
    # baseline.smoking <- "Smoker"
    # baseline.sbp <- 160
    # baseline.bmi <- 25
    # baseline.nonhdl <- 4
    
    ### Assign values from intervention layer
    # intervention.smoking <- input$intervention.smoking
    # intervention.sbp <- input$intervention.sbp
    # intervention.bmi <- input$intervention.bmi
    # intervention.nonhdl <- input$intervention.nonhdl
    intervention.smoking <- "Smoker"
    intervention.sbp <- 130
    intervention.bmi <- 25
    intervention.nonhdl <- 4
    
    ### Apply limits
    if (baseline.sbp < 120){
      baseline.sbp <- 120
    }
    if (intervention.sbp < 120){
      intervention.sbp <- 120
    }
    if (baseline.bmi < 25){
      baseline.bmi <- 25
    }
    if (intervention.bmi < 25){
      intervention.bmi <- 25
    }
    if (baseline.nonhdl < 4){
      baseline.nonhdl <- 4
    }
    if (intervention.nonhdl < 4){
      intervention.nonhdl <- 4
    }
    
    ### Get the difference
    change.sbp <- intervention.sbp - baseline.sbp
    change.bmi <- intervention.bmi - baseline.bmi
    change.nonhdl <- intervention.nonhdl - baseline.nonhdl
    
    ### Read in the relative risks
    direct_RR_smoking_initiation <- readRDS("data/direct_RR_smoking_initiation")
    direct_RR_smoking_cessation <- readRDS("data/direct_RR_smoking_cessation")
    direct_RR_sbp <- readRDS("data/direct_RR_sbp")
    direct_RR_bmi <- readRDS("data/direct_RR_bmi")
    direct_RR_nonhdl <- readRDS("data/direct_RR_nonhdl")
    
    ### Apply the relative risk change
    rr.change.sbp <- direct_RR_sbp^change.sbp
    rr.change.bmi <- direct_RR_bmi^change.bmi
    rr.change.nonhdl <- direct_RR_nonhdl^change.nonhdl
    
    ### Also need to do the same for smoking, which is different because it's non-continuous
    if (baseline.smoking == intervention.smoking){
      rr.change.smoking <- 1
    } else if (baseline.smoking == "Ex-smoker" & intervention.smoking == "Smoker"){
      rr.change.smoking <- 1/direct_RR_smoking_cessation
    } else if (baseline.smoking == "Smoker" & intervention.smoking == "Ex-smoker"){
      rr.change.smoking <- direct_RR_smoking_cessation
    } else if (baseline.smoking == "Non-smoker" & intervention.smoking == "Smoker"){
      rr.change.smoking <- direct_RR_smoking_initiation
    }  else if (baseline.smoking == "Non-smoker" & intervention.smoking == "Ex-smoker"){
      rr.change.smoking <- direct_RR_smoking_initiation*direct_RR_smoking_cessation
    }
    
    ### Apply these RRs
    rr.change.total <-
      rr.change.smoking*
      rr.change.sbp*
      rr.change.bmi*
      rr.change.nonhdl
    
    ### Create output
    risk.table <- create_table(pat, 
                               fit.prototype3 = fit.prototype3, 
                               bhaz.prototype3 = bhaz.prototype3, 
                               means.statins = means.prototype3[["statins"]],
                               means.ah = means.prototype3[["ah"]],
                               means.smoking1 = means.prototype3[["smoking1"]],
                               means.smoking2 = means.prototype3[["smoking2"]],
                               fit.standard = fit.standard, 
                               bhaz.standard = bhaz.standard,
                               rr.change.total = rr.change.total)
  # pat
  # fit.prototype3 = fit.prototype3
  # bhaz.prototype3 = bhaz.prototype3
  # means.statins = means.prototype3[["statins"]]
  # means.ah = means.prototype3[["ah"]]
  # means.smoking1 = means.prototype3[["smoking1"]]
  # means.smoking2 = means.prototype3[["smoking2"]]
  # fit.standard = fit.standard
  # bhaz.standard = bhaz.standard
  
  ### Define decimal places
  dp <- 2
  
  ### Define HR for being on/off treatment
  direct_RR_smoking_initiation <- readRDS("data/direct_RR_smoking_initiation")
  direct_RR_smoking_cessation <- readRDS("data/direct_RR_smoking_cessation")
  direct_RR_nonhdl <- readRDS("data/direct_RR_nonhdl")
  direct_RR_bmi <- readRDS("data/direct_RR_bmi")
  direct_RR_sbp <- readRDS("data/direct_RR_sbp")
  
  ### Create output data.frame
  # pat.risks <- data.frame(matrix(NA, nrow = 1, ncol = 5))
  pat.risks <- vector("character", 3)
  
  ###
  ### Estimate a risk for this individual according to standard cox model
  ###
  risk <- 1 - est_surv_mi(newdata = pat, 
                          fit.list = fit.standard, 
                          bhaz.list = bhaz.standard, 
                          time = 10*365.25)
  pat.risks[1] <- paste(round(100*risk,dp), "%", sep = "")
  
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
  pat.risks[2] <- paste(round(100*baseline.risk,dp), "%", sep = "")
  
  ### Adjust from interventional layer
  
  ### First get values at baseline
  baseline.smoking <- pat.prototype3.baseline$smoking
  baseline.sbp <- pat.prototype3.baseline$sbp
  baseline.bmi <- pat.prototype3.baseline$bmi
  baseline.nonhdl <- pat.prototype3.baseline$nonhdl
  
  ### Assign values from intervention layer
  intervention.smoking <- input$intervention.smoking
  intervention.sbp <- input$intervention.sbp
  intervention.bmi <- input$intervention.bmi
  intervention.nonhdl <- input$intervention.nonhdl
  # intervention.smoking <- "Non-smoker"
  # intervention.sbp <- 130
  # intervention.bmi <- 25
  # intervention.nonhdl <- 4
  
  ### Apply limits
  if (baseline.sbp < 120){
    baseline.sbp <- 120
  }
  if (intervention.sbp < 120){
    intervention.sbp <- 120
  }
  if (baseline.bmi < 25){
    baseline.bmi <- 25
  }
  if (intervention.bmi < 25){
    intervention.bmi <- 25
  }
  if (baseline.nonhdl < 4){
    baseline.nonhdl <- 4
  }
  if (intervention.nonhdl < 4){
    intervention.nonhdl <- 4
  }
  
  ### Get the difference
  change.sbp <- intervention.sbp - baseline.sbp
  change.bmi <- intervention.bmi - baseline.bmi
  change.nonhdl <- intervention.nonhdl - baseline.nonhdl
  
  ### Read in the relative risks
  direct_RR_smoking_initiation <- readRDS("data/direct_RR_smoking_initiation")
  direct_RR_smoking_cessation <- readRDS("data/direct_RR_smoking_cessation")
  direct_RR_sbp <- readRDS("data/direct_RR_sbp")
  direct_RR_bmi <- readRDS("data/direct_RR_bmi")
  direct_RR_nonhdl <- readRDS("data/direct_RR_nonhdl")
  
  ### Apply the relative risk change
  rr.change.sbp <- direct_RR_sbp^change.sbp
  rr.change.bmi <- direct_RR_bmi^change.bmi
  rr.change.nonhdl <- direct_RR_nonhdl^change.nonhdl
  
  ### Also need to do the same for smoking, which is different because it's non-continuous
  if (baseline.smoking == intervention.smoking){
    rr.change.intervention <- 1
  } else if (baseline.smoking == "Ex-smoker" & intervention.smoking == "Smoker"){
    rr.change.intervention <- 1/direct_RR_smoking_cessation
  } else if (baseline.smoking == "Smoker" & intervention.smoking == "Ex-smoker"){
    rr.change.intervention <- direct_RR_smoking_cessation
  } else if (baseline.smoking == "Non-smoker" & intervention.smoking == "Smoker"){
    rr.change.intervention <- direct_RR_smoking_initiation
  }  else if (baseline.smoking == "Non-smoker" & intervention.smoking == "Ex-smoker"){
    rr.change.intervention <- direct_RR_smoking_initiation*direct_RR_smoking_cessation
  }
  
  ### Apply these RRs
  intervention.risk <- 
    baseline.risk*
    rr.change.intervention*
    rr.change.sbp*
    rr.change.bmi*
    rr.change.nonhdl
  
  ### Save it
  pat.risks[3] <- paste(round(100*intervention.risk,dp), "%", sep = "")
  
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
  
  ### Put into a table
  # pat.risks <- apply(round(100*pat.risks, 2), 2, function(x) paste(x, "%", sep = ""))
  pat.risks <- matrix(pat.risks, nrow = 1)
  pat.risks <- datatable(pat.risks, colnames = c("Non-causal\nmodel", "CHARIOT\n no change to treatment", "CHARIOT\under intervention"))
  
  output$riskTable <- renderDT({risk.table})
  
  
