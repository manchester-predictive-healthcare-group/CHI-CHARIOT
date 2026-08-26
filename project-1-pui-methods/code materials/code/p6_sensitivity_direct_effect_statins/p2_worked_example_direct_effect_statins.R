###
### Code for worked example
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

### Define HRs (not log) used for adjusting risks in intervention layer of two-component model
direct_HR_sbp <- exp(lnHR_sbp)
direct_HR_bmi <- exp(lnHR_bmi)
direct_HR_nonhdl <- exp(lnHR_nonhdl)
direct_HR_statin_pleiotropic <- exp(lnHR_statins_direct)

### Read in model and baseline hazards from the model fitted assuming a pleiotropic effect of statins
fit_pleiotropic <- readRDS(paste("data/sens_pleiotropic_fit_", gender, "_model", model, ".rds", sep = ""))
bhaz_pleiotropic <- readRDS(paste("data/sens_pleiotropic_bhaz_", gender, "_model", model, ".rds", sep = ""))

### Read in model and baseline hazards from the model used in main analyses
fit <- readRDS(paste("data/fit_", gender, "_model", model, ".rds", sep = ""))
bhaz <- readRDS(paste("data/bhaz_", gender, "_model", model, ".rds", sep = ""))

### Create a data.frame with the inputs for the indiviudal we are intereseted in
df_valid <- readRDS(paste("data/df_imp_valid_", gender, sep = ""))[0,]
colnames(df_valid)

### Set predictors to be predictors of interest (all medical history variable except diabetes are healthy, no previous history)
df_valid <- tibble::add_row(df_valid, 
                            age = 65,
                            sbp = 140, 
                            bmi = 30, 
                            nonhdl = 5, 
                            IMD = 10,
                            diabetes = factor("Type2", levels = levels(df_valid$diabetes)),
                            ethnicity = factor("white", levels = levels(df_valid$ethnicity)),
                            hypertension      = factor("Absent", levels = c("Absent", "Present")),
                            ra                = factor("Absent", levels = c("Absent", "Present")),
                            af                = factor("Absent", levels = c("Absent", "Present")),
                            ckd               = factor("Absent", levels = c("Absent", "Present")),
                            smi               = factor("Absent", levels = c("Absent", "Present")),
                            fhcvd             = factor("Absent", levels = c("Absent", "Present")),
                            migraine          = factor("Absent", levels = c("Absent", "Present")),
                            sle               = factor("Absent", levels = c("Absent", "Present")),
                            cortico           = factor("Absent", levels = c("Absent", "Present")),
                            antipsy           = factor("Absent", levels = c("Absent", "Present")),
                            copd              = factor("Absent", levels = c("Absent", "Present")),
                            int_dis           = factor("Absent", levels = c("Absent", "Present")),
                            downs             = factor("Absent", levels = c("Absent", "Present")),
                            oral_cancer       = factor("Absent", levels = c("Absent", "Present")),
                            brain_cancer      = factor("Absent", levels = c("Absent", "Present")),
                            lung_cancer       = factor("Absent", levels = c("Absent", "Present")),
                            blood_cancer      = factor("Absent", levels = c("Absent", "Present")),
                            pre_eclampsia     = factor("Absent", levels = c("Absent", "Present")),
                            postnatal_depression = factor("Absent", levels = c("Absent", "Present"))) |>
  as.data.frame()

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

### Function to apply hazard ratio to a risk score, based on change in SBP, BMI and Non-HDL
convert_risk_HR <- function(p, change_sbp = 0, change_bmi = 0, change_nonhdl = 0, change_statin_pleiotropic = 0){
  
  ### Get survival probability
  surv <- 1 - p
  
  ### Get survival probability onto the hazard scale by taking logarithm
  log_surv <- log(surv)
  
  ### Apply the relative risk change
  HR_change_sbp <- direct_HR_sbp^change_sbp
  HR_change_bmi <- direct_HR_bmi^change_bmi
  HR_change_nonhdl <- direct_HR_nonhdl^change_nonhdl
  HR_change_statin_pleiotropic <- direct_HR_statin_pleiotropic^change_statin_pleiotropic
  
  ### Apply these HRs
  HR <-
    HR_change_sbp*
    HR_change_bmi*
    HR_change_nonhdl*
    HR_change_statin_pleiotropic
  
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
### The effect of statins on BMI, SBP and non-HDL cholesterol
###
### Based on the derivation of the DAG document (see supplementary material)
###
### We expect an increase in BMI of 0.33
### We expect a reduction in nonhdl of 1.3268
### We expect a reduction in SBP of 2.62
###
### Plug these changes into convert_risk_HR
###

### Apply odds ratio to risk score to get hypothetical risk under intervention

# We expect an increase in BMI of 0.33
# We expect a reduction in nonhdl of 1.3268
# We expect a reduction in SBP of 2.62

pred_visit0_statins <- 
  convert_risk_HR(p = pred_visit0, 
                    change_sbp = -2.62,
                    change_bmi = 0.33, 
                    change_nonhdl = -1.3268, 
                    change_statin_pleiotropic = 0)

pred_visit0_pmodel_statins <- 
  convert_risk_HR(p = pred_visit0_pmodel,
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
  convert_risk_HR(p = pred_visit1, 
                    change_sbp = 0,
                    change_bmi = 0, 
                    change_nonhdl = 0, 
                    change_statin_pleiotropic = 0)

pred_visit1_pmodel_statins_uneffective <- 
  convert_risk_HR(p = pred_visit1_pmodel,
                    change_sbp = 0,
                    change_bmi = 0, 
                    change_nonhdl = 0, 
                    change_statin_pleiotropic = 1)

## statins effective
# We expect an increase in BMI of 0.33
# We expect a reduction in nonhdl of 1.3268
# We expect a reduction in SBP of 2.62
pred_visit1_statins_effective <- 
  convert_risk_HR(p = pred_visit1, 
                    change_sbp = -2.62,
                    change_bmi = 0.33, 
                    change_nonhdl = -1.3268, 
                    change_statin_pleiotropic = 0)

pred_visit1_pmodel_statins_effective <- 
  convert_risk_HR(p = pred_visit1_pmodel,
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