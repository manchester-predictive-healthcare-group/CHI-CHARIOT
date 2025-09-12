###
### Code for worked example
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Define gender
gender <- 2

### Load survival package
library(survival)

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define HR for offsets
lnHR_statins <- readRDS("data/offsets_total_lnHR_statins.rds")
lnHR_ah <- readRDS("data/offsets_total_lnHR_ah.rds")
lnHR_sbp <- readRDS("data/offsets_direct_lnHR_sbp.rds")
lnHR_bmi <- readRDS("data/offsets_direct_lnHR_bmi.rds")
lnHR_nonhdl <- readRDS("data/offsets_direct_lnHR_nonhdl.rds")

### Read in model and baseline hazards
fit_list <- lapply(1:4, function(x) {readRDS(paste("data/fit_", gender, "_model", x, ".rds", sep = ""))})
bhaz_list <- lapply(1:4, function(x) {readRDS(paste("data/bhaz_", gender, "_model", x, ".rds", sep = ""))})

### Add the non-causal model to the end of the list, to avoid confusion, wanted the model numbers 1 to 4 to correspond to the list element
fit_list[[length(fit_list) + 1]] <- readRDS(paste("data/fit_", gender, "_model", 0, ".rds", sep = ""))
bhaz_list[[length(bhaz_list) + 1]] <- readRDS(paste("data/bhaz_", gender, "_model", 0, ".rds", sep = ""))

########################################################
### Start by comparing risks for model 0 and model 1 ###
########################################################

### Create a data.frame with the inputs for the indiviudal we are intereseted in
df_valid <- readRDS(paste("data/df_imp_valid_", gender, sep = ""))[1,]
colnames(df_valid)

### Set predictors to be predictors of interest
df_valid <- dplyr::mutate(df_valid, 
                          age = 65,
                          sbp = 140, 
                          bmi = 30, 
                          nonhdl = 4, 
                          IMD = 10,
                          diabetes = factor("Type2", levels = levels(diabetes)),
                          ethnicity = factor("white", levels = levels(ethnicity)))


### Create datasets for on and off antihypertensives
df_valid_off_ah <- dplyr::mutate(df_valid, 
                                 offset_ah_timevar_lnHR = 0)

df_valid_on_ah <- dplyr::mutate(df_valid, 
                                offset_ah_timevar_lnHR = lnHR_ah)

###
### Get predicted risks
###

### Non-causal model
pred0 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define sbp
    sbp = 140) |>
  ### Estimate survival probability
  est_surv(fit = fit_list[[5]],
           bhaz = bhaz_list[[5]],
           time = 10*365.25)

### Visit 1, off antihyertensives
pred_model1_off_ah_visit1 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment off antihypertensives
    offset_ah_timevar_lnHR = 0, 
    ### Defin sbp
    sbp = 140) |>
  ### Estimate survival probability
  est_surv_offset(fit = fit_list[[1]],
                  bhaz = bhaz_list[[1]],
                  time = 10*365.25)

### Visit 1, on antihyertensives
pred_model1_on_ah_visit1 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment off antihypertensives
    offset_ah_timevar_lnHR = lnHR_ah, 
    ### Defin sbp
    sbp = 140) |>
  ### Estimate survival probability
  est_surv_offset(fit = fit_list[[1]],
                  bhaz = bhaz_list[[1]],
                  time = 10*365.25)

### Present risks
round(100*c(pred0, pred_model1_off_ah_visit1, pred_model1_on_ah_visit1), 2)


######################################################################################
### Individual initiates hypertensive therapy and achieves reduction of SBP to 120 ### 
######################################################################################

### Visit 2, off antihyertensives
pred_model1_off_ah_visit2 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment off antihypertensives
    offset_ah_timevar_lnHR = 0, 
    ### Defin sbp
    sbp = 120) |>
  ### Estimate survival probability
  est_surv_offset(fit = fit_list[[1]],
                  bhaz = bhaz_list[[1]],
                  time = 10*365.25)

### Visit 2, on antihyertensives
pred_model1_on_ah_visit2 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment off antihypertensives
    offset_ah_timevar_lnHR = lnHR_ah, 
    ### Defin sbp
    sbp = 120) |>
  ### Estimate survival probability
  est_surv_offset(fit = fit_list[[1]],
                  bhaz = bhaz_list[[1]],
                  time = 10*365.25)

### Compare risks to previous
round(100*c(pred_model1_off_ah_visit1, pred_model1_on_ah_visit1, pred_model1_off_ah_visit2, pred_model1_on_ah_visit2), 2)


##########################################################
### Now move onto model 2, which deals with this issue ###
##########################################################

### Visit 1, off antihypertensives
pred_model2_off_ah_visit1 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment off antihypertensives
    offset_ah_timevar_lnHR = 0, 
    ### Define sbp
    sbp_unexposed = 140) |>
  ### Estimate survival probability
  est_surv_offset(fit = fit_list[[2]],
                  bhaz = bhaz_list[[2]],
                  time = 10*365.25)

### Visit 1, on antihypertensives
pred_model2_on_ah_visit1 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment on antihypertensives
    offset_ah_timevar_lnHR = lnHR_ah, 
    ### Define sbp
    sbp_unexposed = 140) |>
  ### Estimate survival probability
  est_surv_offset(fit = fit_list[[2]],
                  bhaz = bhaz_list[[2]],
                  time = 10*365.25)

### There is a small increase in risk compared to model 1
round(100*c(pred_model1_off_ah_visit1, pred_model1_on_ah_visit1, pred_model2_off_ah_visit1, pred_model2_on_ah_visit1), 2)

### At visit 2, the unexposed SBP is the same, so no need to estimate any more risks

#####################################################################
### Now we move onto model 3, which deals with many interventions ###
#####################################################################

###
### Reductions of SBP from 140 to 130 and 120
###

### Get predicted risk for SBP = 140
pred_model3_sbp140 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0,
    ### Define sbp, bmi, and non-HDL with cutoffs
    sbp = 140,
    sbp_adj = dplyr::case_when(sbp < 120 ~ 0,
                               TRUE ~ sbp - 120),
    bmi_adj = dplyr::case_when(bmi < 25 ~ 0,
                               TRUE ~ bmi - 25),
    nonhdl_adj = dplyr::case_when(nonhdl < 2.6 ~ 0,
                                  TRUE ~ nonhdl - 2.6),
    ### Create offset terms for prediction
    offset_sbp_lnHR = lnHR_sbp*sbp_adj,
    offset_bmi_lnHR = lnHR_bmi*bmi_adj,
    offset_nonhdl_lnHR = lnHR_nonhdl*nonhdl_adj)  |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit_list[[3]],
    bhaz = bhaz_list[[3]],
    time = 10*365.25)

### Get predicted risk for SBP = 130
pred_model3_sbp130 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0,
    ### Define sbp, bmi, and non-HDL with cutoffs
    sbp = 130,
    sbp_adj = dplyr::case_when(sbp < 120 ~ 0,
                               TRUE ~ sbp - 120),
    bmi_adj = dplyr::case_when(bmi < 25 ~ 0,
                               TRUE ~ bmi - 25),
    nonhdl_adj = dplyr::case_when(nonhdl < 2.6 ~ 0,
                                  TRUE ~ nonhdl - 2.6),
    ### Create offset terms for prediction
    offset_sbp_lnHR = lnHR_sbp*sbp_adj,
    offset_bmi_lnHR = lnHR_bmi*bmi_adj,
    offset_nonhdl_lnHR = lnHR_nonhdl*nonhdl_adj)  |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit_list[[3]],
    bhaz = bhaz_list[[3]],
    time = 10*365.25)

### Get predicted risk for SBP = 120
pred_model3_sbp120 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0,
    ### Define sbp, bmi, and non-HDL with cutoffs
    sbp = 120,
    sbp_adj = dplyr::case_when(sbp < 120 ~ 0,
                               TRUE ~ sbp - 120),
    bmi_adj = dplyr::case_when(bmi < 25 ~ 0,
                               TRUE ~ bmi - 25),
    nonhdl_adj = dplyr::case_when(nonhdl < 2.6 ~ 0,
                                  TRUE ~ nonhdl - 2.6),
    ### Create offset terms for prediction
    offset_sbp_lnHR = lnHR_sbp*sbp_adj,
    offset_bmi_lnHR = lnHR_bmi*bmi_adj,
    offset_nonhdl_lnHR = lnHR_nonhdl*nonhdl_adj) |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit_list[[3]],
    bhaz = bhaz_list[[3]],
    time = 10*365.25)

### Compare risks
round(100*c(pred_model3_sbp140, pred_model3_sbp130, pred_model3_sbp120), 2)

###
### Lose 5kg intervention 
###

### Calculate the expected reduction in BMI
### We assume average UK height for women, and calculate the weight corresponding to BMI of 30
weight <- 30*1.6^2

### A reduction of 5 will result in a BMI of
bmi_new <- (weight-5)/1.6^2
bmi_reduction = 30 - bmi_new
bmi_reduction

### BMI -> non-HDL is 0.1043451 decrease
nonhdl_reduction <- bmi_reduction*0.1043451
nonhdl_reduction

### BMI -> SBP is 0.7 decrease
sbp_reduction <- bmi_reduction*0.7
sbp_reduction

###
### Get predicted risk for SBP = 140
pred_model3_bmi25 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0,
    ### Define sbp, bmi, and non-HDL with cutoffs
    sbp = 140 - sbp_reduction,
    bmi = 30 - bmi_reduction,
    nonhdl = 4 - nonhdl_reduction,
    sbp_adj = dplyr::case_when(sbp < 120 ~ 0,
                               TRUE ~ sbp - 120),
    bmi_adj = dplyr::case_when(bmi < 25 ~ 0,
                               TRUE ~ bmi - 25),
    nonhdl_adj = dplyr::case_when(nonhdl < 2.6 ~ 0,
                                  TRUE ~ nonhdl - 2.6),
    ### Create offset terms for prediction
    offset_sbp_lnHR = lnHR_sbp*sbp_adj,
    offset_bmi_lnHR = lnHR_bmi*bmi_adj,
    offset_nonhdl_lnHR = lnHR_nonhdl*nonhdl_adj)  |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit_list[[3]],
    bhaz = bhaz_list[[3]],
    time = 10*365.25)

### Compare risks
round(100*pred_model3_bmi25, 2)

###
### Reduction of all risk factors to healthy levels
###

### Get predicted risk for SBP = 140
pred_model3_healthy <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0,
    ### Define sbp, bmi, and non-HDL with cutoffs
    sbp = 120,
    bmi = 25,
    nonhdl = 2.6,
    sbp_adj = dplyr::case_when(sbp < 120 ~ 0,
                               TRUE ~ sbp - 120),
    bmi_adj = dplyr::case_when(bmi < 25 ~ 0,
                               TRUE ~ bmi - 25),
    nonhdl_adj = dplyr::case_when(nonhdl < 2.6 ~ 0,
                                  TRUE ~ nonhdl - 2.6),
    ### Create offset terms for prediction
    offset_sbp_lnHR = lnHR_sbp*sbp_adj,
    offset_bmi_lnHR = lnHR_bmi*bmi_adj,
    offset_nonhdl_lnHR = lnHR_nonhdl*nonhdl_adj)  |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit_list[[3]],
    bhaz = bhaz_list[[3]],
    time = 10*365.25)

### Compare risks
round(100*pred_model3_healthy, 2)

##############################################################################################
### Now we move onto model 4, which deals with many interventions, but retains performance ###
### We will consider the same interventions as in the previous section 
##############################################################################################

### Get predicted risk for SBP = 140 from the initial risk estimation layer
pred_model4_sbp140 <- 1 - df_valid |>
  dplyr::mutate(
    ### Define treatment to be 'continue as current'
    offset_ah_timevar_lnHR = 0, 
    offset_statins_timevar_lnHR = 0,
    ### Define sbp, bmi, and non-HDL with cutoffs
    sbp = 140,
    bmi = 30,
    nonhdl = 4,
    sbp_adj = dplyr::case_when(sbp < 120 ~ 0,
                               TRUE ~ sbp - 120),
    bmi_adj = dplyr::case_when(bmi < 25 ~ 0,
                               TRUE ~ bmi - 25),
    nonhdl_adj = dplyr::case_when(nonhdl < 2.6 ~ 0,
                                  TRUE ~ nonhdl - 2.6),
    ### Create offset terms for prediction
    offset_sbp_lnHR = lnHR_sbp*sbp_adj,
    offset_bmi_lnHR = lnHR_bmi*bmi_adj,
    offset_nonhdl_lnHR = lnHR_nonhdl*nonhdl_adj)  |>
  ### Estimate survival probability
  est_surv_offset(
    fit = fit_list[[4]],
    bhaz = bhaz_list[[4]],
    time = 10*365.25)

### Now we apply the relative risk reduction person A, using the intervention layer

### Read in the odds ratio (for now just assume same as the hazard ratio I have)
direct_OR_sbp <- exp(readRDS("data/offsets_direct_lnHR_sbp.rds"))
direct_OR_bmi <- exp(readRDS("data/offsets_direct_lnHR_bmi.rds"))
direct_OR_nonhdl <- exp(readRDS("data/offsets_direct_lnHR_nonhdl.rds"))

### Function to apply odds ratio to a risk score, based on change in SBP, BMI and Non-HDL
convert_risk_odds <- function(p, change_sbp = 0, change_bmi = 0, change_nonhdl = 0){
  
  ### Apply the relative risk change
  or.change_sbp <- direct_OR_sbp^change_sbp
  or.change_bmi <- direct_OR_bmi^change_bmi
  or.change_nonhdl <- direct_OR_nonhdl^change_nonhdl
  
  ### Apply these RRs
  OR <-
    or.change_sbp*
    or.change_bmi*
    or.change_nonhdl
  
  ### Get odds of risk score
  odds <- p/(1-p)
  
  ### Calculate new odds
  new_odds <- odds*OR
  
  ### Convert back onto risk scale
  new_risk <- new_odds/(1+new_odds)
  
  return(new_risk)
}

###
### Reductions of SBP from 140 to 130 and 120
###

### Apply odds ratio to risk score
pred_model4_sbp130 <- convert_risk_odds(p = pred_model4_sbp140, change_sbp = -10)
pred_model4_sbp120 <- convert_risk_odds(p = pred_model4_sbp140, change_sbp = -20)

### Compare risks
round(100*c(pred_model4_sbp140, pred_model4_sbp130, pred_model4_sbp120), 2)

###
### Reduction of BMI by 5
###
pred_model4_bmi25 <- convert_risk_odds(p = pred_model4_sbp140, change_sbp = -sbp_reduction, change_bmi = -bmi_reduction, change_nonhdl = -nonhdl_reduction)
pred_model4_healthy <- convert_risk_odds(p = pred_model4_sbp140, change_sbp = 20, change_bmi = -5, change_nonhdl = -1.4)

### Compare risks
round(100*pred_model4_bmi25, 2)

###
### Reduction of all values to healthy level
###

### Compare risks
round(100*pred_model4_healthy, 2)