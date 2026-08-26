### Code written by Bowen Jiang
# Goal: Estsimate the effects from the DAG presented in supplementary material, and used to drive the intervention layer of the CHARIOT model
# R script organization:
# 1. Define functions for effect size conversions (OR ↔ RR ↔ HR)
# 2. Estimate the direct effects of Modifiable Risk Factors (MRFs) on CVD, as odds ratios
# 3. Convert all MRF direct effect estimates to hazard ratios (HRs) for the modifiable risk factor model
# 4. Derivation of the total effects of statins and antihypertensives as hazards ratios,
#    used for adjusting for changes in treatment during follow-up #####
# 5. Perform external consistency check for statins effect

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

#################################################
### 1. Functions for converting RR, HR and OR ###
#################################################

###
### the process of converting HR and OR to RR, follow the formula in the Vander J paper (Optimal approximate conversions of odds ratios and hazard ratios to risk ratios)
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
### the process to convert RR to OR (inverse of the above formula)
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
### Function to do a numerical search to convert RR to HR (could not find formulaeic solution to inverse of above equation)
###

### Greater than 1
convert_RR_to_HR_gt1_numerical_search <- function(RR_in, w_in, u_in){

  ### Start by creating a vector of possible HRs
  HRs_possible_gt1 <- seq(1,2,0.00001)

  ### Convert this to a vector of possible RRs
  RRs_possible_gt1 <- sapply(HRs_possible_gt1, convert_HR_to_RR_HRgt1, w = w_in, u = u_in)

  ### Now find the first value that matches the target
  HR_out <- HRs_possible_gt1[min(which(RRs_possible_gt1 > RR_in))]

  return(HR_out)

}

### Smaller than 1
convert_RR_to_HR_lt1_numerical_search <- function(RR_in, w_in, u_in){

  ### Start by creating a vector of possible HRs
  HRs_possible_lt1 <- seq(0,1,0.000001)

  ### Convert this to a vector of possible RRs
  RRs_possible_lt1 <- sapply(HRs_possible_lt1, convert_HR_to_RR_HRlt1, w = w_in, u = u_in)

  ### Now find the first value that matches the target
  HR_out <- HRs_possible_lt1[min(which(RRs_possible_lt1 > RR_in))]

  return(HR_out)

}

###
### To convert HR to OR, first need to convert to RR, then convert RR to OR as shown in above.
###


######################################################
### 2. Estimation of the direct effects in the DAG ###
######################################################

### Define w and u for the conversions
u <- 0.3
w <- 0.00001

###
### Estimate the direct effect of SBP on CVD (section 2.1)
###

# Clinical trials reported RR = 0.8 for a 10 mmHg reduction in SBP (Ettehad et al., 2016).
# Here we scale this effect to 1 mmHg reduction: (0.8)^(1/10).
RR_SBP <- (0.8)^(1/10)
RR_SBP   # Risk ratio per 1 mmHg reduction in SBP

# Convert the RR to HR using the RR < 1 conversion function
HR_SBP <- convert_RR_to_HR_lt1_numerical_search(RR_SBP, w, u)[1]
HR_SBP #  Total effect of Non-HDL on CVD (odds ratio)

# total effect equals direct effect, and record per 1mmHg increase
HR_SBP_direct <- 1/HR_SBP
HR_SBP_direct #  #1.02463

###
### Estimate the direct effect of NonHDL on CVD (section 2.3)
###

# For a 1 mmol/L reduction in Non-HDL cholesterol,
# the reported total effect is HR = 0.8184262,
# calculated as (0.78)^(1/1.24) based on the JBS3 model.
HR_NonHDL_total <- (0.78)^(1/1.24)
HR_NonHDL_total # 0.8184262 total effect on CVD

# Express as the effect per 1 mmol/L increase in Non-HDL
HR_NonHDL_total <- 1/HR_NonHDL_total

# Estimate the indirect effect of Non-HDL on CVD via SBP:
# A 1 mmol/L reduction in Non-HDL corresponds to a SBP (derived as (0.4258/1.0244)*0.6615) reduction in SBP (section 2.2)
NonHDL_to_SBP <- (0.4258/1.0244)*0.6615
NonHDL_to_SBP # 0.2749577 mmHg decrease in SBP

# Therefore, the indirect effect of Non-HDL via SBP is:
HR_NonHDL_indirect <- (HR_SBP_direct)^(NonHDL_to_SBP)
HR_NonHDL_indirect # 1.006713

# The direct effect of Non-HDL on CVD is then:
# (Total effect) / (Indirect effect via SBP)
HR_NonHDL_direct <- HR_NonHDL_total / HR_NonHDL_indirect
HR_NonHDL_direct # 1.21371

# Express as the effect per 1 mmol/L reduction in Non-HDL
1 / HR_NonHDL_direct # 0.82392

###
### Estimate the direct effect of BMI on CVD (section 2.6)
###

# BMI differences between groups
# Overweight vs. Normal weight = 5 BMI units
# Obesity vs. Normal weight = 10 BMI units
HR_BMI_5 <- log(1.12)/5  # = 0.02266574. Convert HR (1.12, Overweight vs. Normal) into log-risk increase per 1 BMI unit
HR_BMI_10 <- log(1.22)/10 # = 0.01988509. Convert HR (1.22, Obesity vs. Normal) into log-risk increase per 1 BMI unit

HR_BMI_average <- (HR_BMI_5 + HR_BMI_10)/2 # Average log-risk increase per 1 BMI unit across both comparisons
HR_BMI_direct <- exp(HR_BMI_average)              # Convert back from log scale
HR_BMI_direct # HR = 1.021503 (per 1 BMI unit increase)
1/HR_BMI_direct # Inverse: HR for per 1 BMI unit decrease = 0.9789493

#########################################################################################################
### 3. Derivation of the total effects of statins and antihypertensives as hazards ratios, ###
### used for adjusting for changes in treatment during follow-up #####
#########################################################################################################

###
### To convert RR to either HR or OR, we will find the inverse solution to the above equation using a numerical search
### To do this, we create a vector of possible HRs or ORs, that we will put into the function, in order to find a matching
### RR to the RR we are trying to convert.
###

### Create a vector of possible HRs
HRs_possible_gt1 <- seq(1,2,0.00001)
HRs_possible_lt1 <- seq(0,1,0.00001)

### Convert these to RRs
RRs_possible_HRs_gt1 <- sapply(HRs_possible_gt1, convert_HR_to_RR_HRgt1, w = 0.00001, u = 0.3)
RRs_possible_HRs_lt1 <- sapply(HRs_possible_lt1, convert_HR_to_RR_HRlt1, w = 0.00001, u = 0.3)

### Create a vector of possible ORs
ORs_possible_gt1 <- seq(1,2,0.00001)
ORs_possible_lt1 <- seq(0,1,0.00001)

### Convert these to RRs
RRs_possible_ORs_gt1 <- sapply(ORs_possible_gt1, convert_OR_to_RR_HRgt1, w = 0.00001, u = 0.3)
RRs_possible_ORs_lt1 <- sapply(ORs_possible_lt1, convert_OR_to_RR_HRlt1, w = 0.00001, u = 0.3)

###
### Antihypertensives
### Want to convert to HR for adjusting for treatment drop in during model fitting
###

### We have RR from literature
RR_ah_total <- 0.74

### Now find the first value that matches the target
HR_ah_total <- HRs_possible_lt1[min(which(RRs_possible_HRs_lt1 > RR_ah_total))]
HR_ah_total # 0.72285

### Double check these conversions make sense, by converting the HR back to RR
convert_HR_to_RR_HRlt1(HR_ah_total, w = 0.00001, u = 0.3)
RR_ah_total

###
### Statins
### Want to convert to HR for adjusting for treatment drop in during model fitting
###

### We have RR from literature
RR_statins_total <- 0.75

### Now find the first value that matches the target
HR_statins_total <- HRs_possible_lt1[min(which(RRs_possible_HRs_lt1 > RR_statins_total))]
HR_statins_total # 0.73327

### Double check these conversions make sense, by converting the HR back to RR
convert_HR_to_RR_HRlt1(HR_statins_total, w = 0.00001, u = 0.3)
RR_statins_total


#############################
### 4. Consistency checks ###
#############################

###
### Statins (section 4.1)
###

# Consistency check for statin total effect using DAG-based decomposition

# --- Indirect effects through mediators ---
HR_statins_indirect_BMI <- (HR_BMI_direct)^(0.33)    
HR_statins_indirect_BMI # Effect through BMI = 1.007046
HR_statins_indirect_NonHDL <- (HR_NonHDL_direct)^(-1.3268)    
HR_statins_indirect_NonHDL # Effect through Non-HDL = 0.7733859
HR_statins_indirect_SBP <- (HR_SBP_direct)^(-2.62)      
HR_statins_indirect_SBP # Effect through SBP = 0.9382406

# --- Combine indirect effects to obtain total effect ---
HR_statins_total_DAG <- HR_statins_indirect_BMI * HR_statins_indirect_NonHDL * HR_statins_indirect_SBP
HR_statins_total_DAG # Total HR = 0.7307345

# --- Convert OR to RR for external comparison (JAMA meta-analysis RR = 0.75) ---
RR_statins_total_DAG <- convert_HR_to_RR_HRlt1(HR_statins_total_DAG, w, u)
RR_statins_total_DAG # 0.7475701

###########################################
### Validation of the inversion formula ###
###########################################

### We derived the inverse of the formula in Van der Weele and put it into a function, we should validate this

### lt1
# convert_RR_to_OR_lt1 is the inversion
# convert_OR_to_RR_ORlt1 is the original
new_RR <- convert_OR_to_RR_ORlt1(0.9, w, u)
convert_RR_to_OR_lt1(new_RR, w, u)[1]

# Check conversion back to RR gives correct answer
testthat::expect_equal(0.9, convert_RR_to_OR_lt1(new_RR, w, u)[1])

### gt1
# convert_RR_to_OR_lt1 is the inversion
# convert_OR_to_RR_ORlt1 is the original
new_RR <- convert_OR_to_RR_ORgt1(1.1, w, u)
convert_RR_to_OR_gt1(new_RR, w, u)[1]

# Check conversion back to RR gives correct answer
testthat::expect_equal(1.1, convert_RR_to_OR_gt1(new_RR, w, u)[1])

###############################################################################################################
### 5. Showcase arbitrary intervention (see section 3 of supplementary material file 3 - derivation of DAG) ###
###############################################################################################################

## sbp
arbitrary_change_sbp <- -7.5
## bmi 
arbitrary_change_bmi <- -5
## non-HDL cholesterol
arbitrary_change_nonhdl <- -1

### Get combined change
HR_comb <- (HR_SBP_direct)^arbitrary_change_sbp*(HR_BMI_direct)^arbitrary_change_bmi*(HR_NonHDL_direct)^arbitrary_change_nonhdl

### Get survival probability
surv <- 1 - 0.1688

### Get survival probability onto the hazard scale by taking logarithm
log_surv <- log(surv)

### Apply hazard ratios
log_surv_adjusted <- log_surv*HR_comb

### Calculate survival probability (note, we do not take the negative, 
### as we did not do this when converted back onto the hazard scale)
new_surv <- exp(log_surv_adjusted)

### Convert back onto risk scale
new_risk <- 1 - new_surv
new_risk # 0.1078426

###################
###################
### SAVE VALUES ###
###################
###################

###
### Save log-hazard ratios
###
saveRDS(log(HR_SBP_direct), "data/offsets_direct_lnHR_sbp.rds")
saveRDS(log(HR_BMI_direct), "data/offsets_direct_lnHR_bmi.rds")
saveRDS(log(HR_NonHDL_direct), "data/offsets_direct_lnHR_nonhdl.rds")
saveRDS(log(HR_statins_total), "data/offsets_total_lnHR_statins.rds")
saveRDS(log(HR_ah_total), "data/offsets_total_lnHR_ah.rds")