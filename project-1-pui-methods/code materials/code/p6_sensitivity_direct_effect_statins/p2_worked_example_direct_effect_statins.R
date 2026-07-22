###
### Code for worked example (see relevant .Rmd file in p9_sensitivity_summaries)
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Define gender
gender <- 2

### Define model
model <- 4

### Load survival package
library(survival)

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define HR for offsets
lnHR_ah <- readRDS("data/offsets_total_lnHR_ah.rds")
lnHR_sbp <- readRDS("data/offsets_direct_lnHR_sbp.rds")
lnHR_bmi <- readRDS("data/offsets_direct_lnHR_bmi.rds")
lnHR_nonhdl <- readRDS("data/offsets_direct_lnHR_nonhdl.rds")

### Read in statin effects
lnHR_statins_direct <- readRDS("data/sens_pleiotropic_offsets_direct_lnHR_statins.rds")
lnHR_statins_total <- readRDS("data/sens_pleiotropic_offsets_total_lnHR_statins.rds")

### Read in model and baseline hazards from the model fitted assuming a pleiotropic effect of statins
fit_pleiotropic <- readRDS(paste("data/sens_pleiotropic_fit_", gender, "_model", model, ".rds", sep = ""))
bhaz_pleiotropic <- readRDS(paste("data/sens_pleiotropic_bhaz_", gender, "_model", model, ".rds", sep = ""))

### Read in model and baseline hazards from the model used in main analyses
fit <- readRDS(paste("data/fit_", gender, "_model", model, ".rds", sep = ""))
bhaz <- readRDS(paste("data/bhaz_", gender, "_model", model, ".rds", sep = ""))

### Create a data.frame with the inputs for the indiviudal we are intereseted in
df_valid <- readRDS(paste("data/df_imp_valid_", gender, sep = ""))[1,]
colnames(df_valid)

### Set predictors to be predictors of interest
df_valid <- dplyr::mutate(df_valid, 
                          age = 65,
                          sbp = 140, 
                          bmi = 30, 
                          nonhdl = 5, 
                          IMD = 10,
                          diabetes = factor("Type2", levels = levels(diabetes)),
                          ethnicity = factor("white", levels = levels(ethnicity)))

###
### Get predicted risk from the initial risk estimation layer using model from main analyses
###
pred_visit0 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0) |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit,
    bhaz = bhaz,
    time = 10*365.25)

###
### Get predicted risk from the initial risk estimation layer using model with pleiotropic statin effect
###
pred_visit0_pmodel <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0)  |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit_pleiotropic,
    bhaz = bhaz_pleiotropic,
    time = 10*365.25)

#############################################################
### Convert HR to OR to direct effect of statins    
### Needed to the estimate risks in the intervention layer 
#############################################################

###
### Define functions to convert HR to OR (see Van Der Weele, optimal minimax conversion)
###

### HR > 1, w < p0 < p1 < u
convert_HR_to_RR_HRgt1 <- function(HR, w_in, u_in){
  top <- 1 - (1-w_in)^HR
  bot <- 1 - (1-u_in)^(1/HR)
  out <- ((top/bot)*(u_in/w_in))^(1/2)
  return(out)
}

### HR < 1, w < p1 < p0 < u
convert_HR_to_RR_HRlt1 <- function(HR, w_in, u_in){
  top <- 1 - (1-u_in)^HR
  bot <- 1 - (1-w_in)^(1/HR)
  out <- ((top/bot)*(w_in/u_in))^(1/2)
  return(out)
}

### OR > 1, w < p0 < p1 < u
convert_OR_to_RR_ORgt1 <- function(OR, w_in, u_in){
  top <- OR*(u_in+OR-OR*u_in)
  bot <- 1 - w_in + OR*w_in
  out <- (top/bot)^(1/2)
  return(out)
}

### OR < 1, w < p1 < p0 < u
convert_OR_to_RR_ORlt1 <- function(OR, w_in, u_in){
  top <- OR*(w_in+OR-OR*w_in)
  bot <- 1 - u_in + OR*u_in
  out <- (top/bot)^(1/2)
  return(out)
}

###
### The process to convert RR to OR (inverse of the above formula)
###

### OR > 1
convert_RR_to_OR_gt1 <- function(RR, w_in, u_in) {
  a <- 1 - u_in
  b <- u_in - RR^2 * w_in
  c <- -RR^2 * (1 - w_in)
  
  D <- b^2 - 4 * a * c
  if (D < 0) {
    return(NA)  # No real solution
  } else {
    x1 <- (-b + sqrt(D)) / (2 * a)
    x2 <- (-b - sqrt(D)) / (2 * a)
    return(c(x1, x2))  # return both roots
  }
}

### OR < 1
convert_RR_to_OR_lt1 <- function(RR, w_in, u_in) {
  a <- 1 - w_in
  b <- w_in - RR^2 * u_in
  c <- -RR^2 * (1 - u_in)
  
  D <- b^2 - 4 * a * c
  if (D < 0) {
    return(NA)  # No real solution
  } else {
    x1 <- (-b + sqrt(D)) / (2 * a)
    x2 <- (-b - sqrt(D)) / (2 * a)
    return(c(x1, x2))  # return both roots
  }
}

###
### Write a function that will convert an HR to RR, and then RR to OR. 
### This is the two step process needed to convert HR to OR.

# When ratio is less than 1
convert_HR_to_OR_HRlt1 <- function(HR, w, u){
  # Convert to HR to RR
  rr <- convert_HR_to_RR_HRlt1(HR, w, u)
  # Convert RR to OR
  or <- convert_RR_to_OR_lt1(rr, w, u)
  # Return the positive element
  return(or[which(or > 0)])
}

# When ratio is bigger than 1
convert_HR_to_OR_HRgt1 <- function(HR, w, u){
  # Convert to HR to RR
  rr <- convert_HR_to_RR_HRgt1(HR, w, u)
  # Convert RR to OR
  or <- convert_RR_to_OR_gt1(rr, w, u)
  # Return
  return(or[which(or > 0)])
}

### Read in the hazard ratios (for now just assume same as the hazard ratio I have)
direct_HR_statin_pleiotropic <- exp(lnHR_statins_direct)

###
### Run this function to calculate the OR for direct effect of statins
direct_OR_statin_pleiotropic <- convert_HR_to_OR_HRlt1(direct_HR_statin_pleiotropic, w = 0.00001, u = 0.3)

### Read in other direct effects, calculated in file p1_data_prep/p8_calculate_treatment_effects.R
direct_OR_nonhdl <- readRDS("data/direct_OR_nonhdl.rds")
direct_OR_bmi <-  readRDS("data/direct_OR_bmi.rds")
direct_OR_sbp <- readRDS("data/direct_OR_sbp.rds")

###
### Function to apply odds ratio to a risk score, based on change in SBP, BMI, Non-HDL, and direct effect of statins
###
convert_risk_odds <- function(p, change_sbp = 0, change_bmi = 0, change_nonhdl = 0, change_statin_pleiotropic = 0){
  
  ### Apply the relative risk change
  or.change_sbp <- direct_OR_sbp^change_sbp
  or.change_bmi <- direct_OR_bmi^change_bmi
  or.change_nonhdl <- direct_OR_nonhdl^change_nonhdl
  or.change_statin_pleiotropic <- direct_OR_statin_pleiotropic^change_statin_pleiotropic
  
  ### Apply these RRs
  OR <-
    or.change_sbp*
    or.change_bmi*
    or.change_nonhdl*
    or.change_statin_pleiotropic
  
  ### Get odds of risk score
  odds <- p/(1-p)
  
  ### Calculate new odds
  new_odds <- odds*OR
  
  ### Convert back onto risk scale
  new_risk <- new_odds/(1+new_odds)
  
  return(new_risk)
}

###
### The effect of statins on BMI, SBP and non-HDL cholesterol
###
### Based on the derivation of the DAG document (see supplementary material)
###
### We expect an increase in BMI of 0.33
### We expect a reduction in nonhdl of 1.3268
### We expect a reduction in SBP of 2.62
###
### Plug these changes into convert_risk_odds
###

### Apply odds ratio to risk score to get hypothetical risk under intervention

# We expect an increase in BMI of 0.33
# We expect a reduction in nonhdl of 1.3268
# We expect a reduction in SBP of 2.62

pred_visit0_statins <- 
  convert_risk_odds(p = pred_visit0, 
                    change_sbp = -2.62,
                    change_bmi = 0.33, 
                    change_nonhdl = -1.3268, 
                    change_statin_pleiotropic = 0)

pred_visit0_pmodel_statins <- 
  convert_risk_odds(p = pred_visit0_pmodel,
                    change_sbp = -2.62,
                    change_bmi = 0.33, 
                    change_nonhdl = -1.3268, 
                    change_statin_pleiotropic = 1)

###
### Get risks at return visit
###

###
### Get predicted risk from the initial risk estimation layer using model from main analyses
### This is the risk if we had not intervened
###
pred_visit1 <- 1 - df_valid |>
  dplyr::mutate(
    age = age + 1,
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0)  |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit,
    bhaz = bhaz,
    time = 10*365.25)

###
### Get predicted risk from the initial risk estimation layer using model with pleiotropic statin effect
### This is the risk if we had not intervened
###
pred_visit1_pmodel <- 1 - df_valid |>
  dplyr::mutate(
    age = age + 1,
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0)  |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit_pleiotropic,
    bhaz = bhaz_pleiotropic,
    time = 10*365.25)


### Apply odds ratio to risk score to get factual risk based on achieved changes in MFRs
## statins uneffective
pred_visit1_statins_uneffective <- 
  convert_risk_odds(p = pred_visit1, 
                    change_sbp = 0,
                    change_bmi = 0, 
                    change_nonhdl = 0, 
                    change_statin_pleiotropic = 0)

pred_visit1_pmodel_statins_uneffective <- 
  convert_risk_odds(p = pred_visit1_pmodel,
                    change_sbp = 0,
                    change_bmi = 0, 
                    change_nonhdl = 0, 
                    change_statin_pleiotropic = 1)

## statins effective
# We expect an increase in BMI of 0.33
# We expect a reduction in nonhdl of 1.3268
# We expect a reduction in SBP of 2.62
pred_visit1_statins_effective <- 
  convert_risk_odds(p = pred_visit1, 
                    change_sbp = -2.62,
                    change_bmi = 0.33, 
                    change_nonhdl = -1.3268, 
                    change_statin_pleiotropic = 0)

pred_visit1_pmodel_statins_effective <- 
  convert_risk_odds(p = pred_visit1_pmodel,
                    change_sbp = -2.62,
                    change_bmi = 0.33, 
                    change_nonhdl = -1.3268, 
                    change_statin_pleiotropic = 1)

#########################
### Summarise results ###
#########################

### Risk at visit0 from each model
round(100*pred_visit0, 2)
round(100*pred_visit0_pmodel, 2)

### Risk at visit1 from each model, if we had not intervened
round(100*pred_visit1, 2)
round(100*pred_visit1_pmodel, 2)

### Risk at visit1 from each model, if we had not intervened, statins uneffective
round(100*pred_visit1_statins_uneffective, 2)
round(100*pred_visit1_pmodel_statins_uneffective, 2)

### Risk at visit1 from each model, if we had not intervened, statins effective
round(100*pred_visit1_statins_effective, 2)
round(100*pred_visit1_pmodel_statins_effective, 2)

### Save outputs for use in Rmd
saveRDS(pred_visit0, "data/sens_pleiotropic_pred_visit0.rds")
saveRDS(pred_visit0_pmodel, "data/sens_pleiotropic_pred_visit0_pmodel.rds")
saveRDS(pred_visit0_statins, "data/sens_pleiotropic_pred_visit0_statins.rds")
saveRDS(pred_visit0_pmodel_statins, "data/sens_pleiotropic_pred_visit0_pmodel_statins.rds")
saveRDS(pred_visit1, "data/sens_pleiotropic_pred_visit1.rds")
saveRDS(pred_visit1_pmodel, "data/sens_pleiotropic_pred_visit1_pmodel.rds")
saveRDS(pred_visit1_statins_uneffective, "data/sens_pleiotropic_pred_visit1_statins_uneffective.rds")
saveRDS(pred_visit1_pmodel_statins_uneffective, "data/sens_pleiotropic_pred_visit1_pmodel_statins_uneffective.rds")
saveRDS(pred_visit1_statins_effective, "data/sens_pleiotropic_pred_visit1_statins_effective.rds")
saveRDS(pred_visit1_pmodel_statins_effective, "data/sens_pleiotropic_pred_visit1_pmodel_statins_effective.rds")