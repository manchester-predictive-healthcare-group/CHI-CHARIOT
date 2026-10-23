### Code written by Bowen Jiang, edited by Alex Pate
# Goal: Estsimate the effects from the DAG presented in supplementary material, and used to drive the intervention layer of the CHARIOT model
# R script organization:
# 1. Define functions for effect size conversions (OR ↔ RR ↔ HR)
# 2. Estimate the direct effects of Modifiable Risk Factors (MRFs) on CVD, as hazard ratios
# 3. Convert all MRF direct effect estimates to hazard ratios (HRs) for the modifiable risk factor model
# 4. Derivation of the total effects of statins and antihypertensives as hazards ratios,
#    used for adjusting for changes in treatment during follow-up #####
# 5. Perform external consistency check for statins effect

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

#################################################
### 1. Functions for converting RR, HR and OR ###
#################################################

###
### The process of converting HR and OR to RR, follow the formula in the Vander J paper (Optimal approximate conversions of odds ratios and hazard ratios to risk ratios)
### We optimise the conversion for risks between 0 and 0.3, assuming individuals with risks higher than this will need
### specialised treatment, and are less common
###

###
### HR to RR functions
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

###
### OR to RR functions
###

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
### Function to do a numerical search to convert RR to HR (could not find formulaeic solution to inverse of above equation)
###

### Greater than 1
convert_RR_to_HR_gt1_numerical_search <- function(RR_in, w_in, u_in){

  ### Start by creating a vector of possible HRs
  HRs_possible_gt1 <- seq(1,2,0.000001)

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

#######################################################################################################
### 2. Derivation of the total effects of statins, antihypertensives and smoking as hazards ratios, ###
### used for adjusting for changes in treatment during follow-up ########
#######################################################################################################

###
### Antihypertensives
### Want to convert to HR for adjusting for treatment drop in during model fitting
###

### We have RR from literature
RR_ah_total <- 0.74

### Now find the first value that matches the target
HR_ah_total <- convert_RR_to_HR_lt1_numerical_search(RR_ah_total, w_in = 0.000001, u_in = 0.3)

### Double check these conversions make sense, by converting the HR back to RR
convert_HR_to_RR_HRlt1(HR_ah_total, w = 0.000001, u = 0.3)
RR_ah_total

###
### Statins
### Want to convert to HR for adjusting for treatment drop in during model fitting
###

### We have RR from literature
RR_statins_total <- 0.75

### Now find the first value that matches the target
HR_statins_total <- convert_RR_to_HR_lt1_numerical_search(RR_statins_total, w_in = 0.000001, u_in = 0.3)

### Double check these conversions make sense, by converting the HR back to RR
convert_HR_to_RR_HRlt1(HR_statins_total, w = 0.000001, u = 0.3)
RR_statins_total

###
### Smoking
###

### We have the following effects from literature
OR_smok_init_total <- 1.44
OR_smok_cess_total <- 1/1.25

### Want to convert OR to HR

### Need to first convert OR to RR
RR_smok_init_total <- convert_OR_to_RR_ORgt1(OR = OR_smok_init_total, w = 0.000001, u = 0.3)
RR_smok_init_total

RR_smok_cess_total <- convert_OR_to_RR_ORlt1(OR = OR_smok_cess_total, w = 0.000001, u = 0.3)
RR_smok_cess_total

### Then convert RR to HR

### Now find the first value that matches the target
HR_smok_init_total <- convert_RR_to_HR_gt1_numerical_search(RR_smok_init_total, w_in = 0.000001, u_in = 0.3)
HR_smok_cess_total <- convert_RR_to_HR_lt1_numerical_search(RR_smok_cess_total, w_in = 0.000001, u_in = 0.3)
HR_smok_init_total # 1.4064
HR_smok_cess_total # 0.81217

### Double check these conversions make sense, by converting the HR back to RR
convert_HR_to_RR_HRgt1(HR_smok_init_total, w = 0.000001, u = 0.3)
RR_smok_init_total
convert_HR_to_RR_HRlt1(HR_smok_cess_total, w = 0.000001, u = 0.3)
RR_smok_cess_total


######################################################
### 3. Estimation of the direct effects in the DAG ###
######################################################

### Define w and u for the conversions
u <- 0.3
w <- 0.000001

###
### Estimate the direct effect of sbp on CVD (section 2.1)
###

# Clinical trials reported RR = 0.8 for a 10 mmHg reduction in sbp (Ettehad et al., 2016).
# Here we scale this effect to 1 mmHg reduction: (0.8)^(1/10).
RR_sbp <- (0.8)^(1/10)
RR_sbp   # Risk ratio per 1 mmHg reduction in sbp

# Convert the RR to HR using the RR < 1 conversion function
HR_sbp <- convert_RR_to_HR_lt1_numerical_search(RR_sbp, w, u)[1]
HR_sbp #  Total effect of Non-HDL on CVD (hazard ratio)

# total effect equals direct effect, and record per 1mmHg increase
HR_sbp_direct <- 1/HR_sbp
HR_sbp_direct #  #1.02463

###
### Estimate the direct effect of nonhdl on CVD (section 2.3)
###

# For a 1 mmol/L reduction in Non-HDL cholesterol,
# the reported total effect is HR = 0.8184262,
# calculated as (0.78)^(1/1.24) based on the JBS3 model.
HR_nonhdl_total <- (0.78)^(1/1.24)
HR_nonhdl_total # 0.8184262 total effect on CVD

# Express as the effect per 1 mmol/L increase in Non-HDL
HR_nonhdl_total <- 1/HR_nonhdl_total

# Estimate the indirect effect of Non-HDL on CVD via sbp:
# A 1 mmol/L reduction in Non-HDL corresponds to a sbp (derived as (0.4258/1.0244)*0.6615) reduction in sbp (section 2.2)
nonhdl_to_sbp <- (0.4258/1.0244)*0.6615
nonhdl_to_sbp # 0.2749577 mmHg decrease in sbp

# Therefore, the indirect effect of Non-HDL via sbp is:
HR_nonhdl_indirect <- (HR_sbp_direct)^(nonhdl_to_sbp)
HR_nonhdl_indirect # 1.006713

# The direct effect of Non-HDL on CVD is then:
# (Total effect) / (Indirect effect via sbp)
HR_nonhdl_direct <- HR_nonhdl_total / HR_nonhdl_indirect
HR_nonhdl_direct # 1.21371

# Express as the effect per 1 mmol/L reduction in Non-HDL
1 / HR_nonhdl_direct # 0.82392

###
### Estimate the direct effect of bmi on CVD (section 2.6)
###

# bmi differences between groups
# Overweight vs. Normal weight = 5 bmi units
# Obesity vs. Normal weight = 10 bmi units
HR_bmi_5 <- log(1.12)/5  # = 0.02266574. Convert HR (1.12, Overweight vs. Normal) into log-risk increase per 1 bmi unit
HR_bmi_10 <- log(1.22)/10 # = 0.01988509. Convert HR (1.22, Obesity vs. Normal) into log-risk increase per 1 bmi unit

HR_bmi_average <- (HR_bmi_5 + HR_bmi_10)/2 # Average log-risk increase per 1 bmi unit across both comparisons
HR_bmi_direct <- exp(HR_bmi_average)              # Convert back from log scale
HR_bmi_direct # HR = 1.021503 (per 1 bmi unit increase)
1/HR_bmi_direct # Inverse: HR for per 1 bmi unit decrease = 0.9789493


###
### Estimating the direct effect of smoking cessation on CVD (section 2.10)
###

# Step 1: Estimate changes in risk factors due to smoking cessation

# bmi
smok_cess_to_bmi <- 0.66 # average taken from literature (section 2.7)

# nonhdl
smok_cess_to_nonhdl <- 0.2*0.66/4.5*2.348
smok_cess_to_nonhdl # = 0.06887467 mmol/L increase in Non-HDL (section 2.8)

# sbp
smok_cess_to_sbp <- -3.5 # reported directly in literature (section 2.9)

# Step 2: Indirect effects of smoking cessation on CVD

# (a) via bmi
HR_smok_cess_indirect_bmi <- (HR_bmi_direct)^smok_cess_to_bmi
HR_smok_cess_indirect_bmi # = 1.014141 (HR mediated through bmi)

# (b) via Non-HDL
HR_smok_cess_indirect_nonhdl <- (HR_nonhdl_direct)^smok_cess_to_nonhdl
HR_smok_cess_indirect_nonhdl # = 1.013429 (HR mediated through Non-HDL)

# (c) via sbp
HR_smok_cess_indirect_sbp <- (HR_sbp_direct)^(smok_cess_to_sbp)
HR_smok_cess_indirect_sbp # = 0.9183647 (protective effect via sbp reduction)

# Step 3: Combine all indirect pathways
HR_smok_cess_indirect <- HR_smok_cess_indirect_bmi * HR_smok_cess_indirect_nonhdl * HR_smok_cess_indirect_sbp
HR_smok_cess_indirect # = 0.9438584 (total indirect effect)

# Step 4: Estimate direct effect
# Total effect from trials: OR = 0.8
# Direct effect = Total effect / Indirect effect
HR_smok_cess_direct <- HR_smok_cess_total / HR_smok_cess_indirect
HR_smok_cess_direct # = 0.8604713 (direct effect of smoking cessation on CVD)

###
### Direct effect of smoking initiation on CVD estimation (section 2.14)
###

# --- Step 1: Calculate indirect effects of smoking initiation ---

# effects of initiating smoking on bmi, nonhdl and sbp
# these all acts through bmi
smok_init_on_bmi <- -0.61 # decrease in bmi (units) reported in literature (section 2.11)

# the effect of bmi on nonhdl
bmi_on_nonhdl <- 0.2 * (1 / 4.5) * 2.348 # (section 2.5)
# the effect of smok init on nonhdl
smok_init_on_nonhdl <- smok_init_on_bmi * bmi_on_nonhdl # (section 2.12)
smok_init_on_nonhdl # estimated effect on nonhdl = 0.06365689 mmol/L decrease

# the effect of bmi on sbp
bmi_on_sbp <- 0.7 # (section 2.4)
# the effect of smok init on sbp
smok_init_on_sbp <- smok_init_on_bmi * bmi_on_sbp # (Section 2.13)
smok_init_on_sbp # estimated effect on sbp = 0.427 mmHg decrease

# Convert each pathway effect to hazard ratios (HRs):
HR_smok_init_indirect_bmi <- (HR_bmi_direct)^(smok_init_on_bmi)           # indirect effect through bmi → HR = 0.9871059
HR_smok_init_indirect_bmi
HR_smok_init_indirect_nonhdl <- (HR_nonhdl_direct)^(smok_init_on_nonhdl)     # indirect effect through Non-HDL → HR = 0.9877465
HR_smok_init_indirect_nonhdl
HR_smok_init_indirect_sbp <- (HR_sbp_direct)^(smok_init_on_sbp)          # indirect effect through sbp → HR = 0.9896642
HR_smok_init_indirect_sbp

# Combine all indirect effects:
HR_smok_init_indirect <- HR_smok_init_indirect_sbp * HR_smok_init_indirect_nonhdl * HR_smok_init_indirect_bmi  # total indirect effect HR = 0.9649328
HR_smok_init_indirect
# --- Step 2: Calculate direct effect ---
# Using total effect derived above
HR_smok_init_direct <- HR_smok_init_total / HR_smok_init_indirect
HR_smok_init_direct # direct effect HR = 1.457506


#############################
### 4. Consistency checks ###
#############################

###
### Statins (section 4.1)
###

# Consistency check for statin total effect using DAG-based decomposition

# --- Indirect effects through mediators ---
HR_statins_indirect_bmi <- (HR_bmi_direct)^(0.33)    
HR_statins_indirect_bmi # Effect through bmi = 1.007046
HR_statins_indirect_nonhdl <- (HR_nonhdl_direct)^(-1.3268)    
HR_statins_indirect_nonhdl # Effect through Non-HDL = 0.7733859
HR_statins_indirect_sbp <- (HR_sbp_direct)^(-2.62)      
HR_statins_indirect_sbp # Effect through sbp = 0.9382406

# --- Combine indirect effects to obtain total effect ---
HR_statins_total_DAG <- HR_statins_indirect_bmi * HR_statins_indirect_nonhdl * HR_statins_indirect_sbp
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

##########################################
### 5. Showcase arbitrary intervention ###
##########################################

## sbp
arbitrary_change_sbp <- -7.5
## bmi 
arbitrary_change_bmi <- -5
## non-HDL cholesterol
arbitrary_change_nonhdl <- -1

### Get combined change
HR_comb <- (HR_sbp_direct)^arbitrary_change_sbp*(HR_bmi_direct)^arbitrary_change_bmi*(HR_nonhdl_direct)^arbitrary_change_nonhdl

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

### Convert risk from exemplar
convert_risk_HR(1-surv, HR_comb) # 0.1078426


###################
###################
### SAVE VALUES ###
###################
###################

###
### Save log-hazard ratios
###
saveRDS(log(HR_sbp_direct), "data/p5/offsets_direct_lnHR_sbp.rds")
saveRDS(log(HR_bmi_direct), "data/p5/offsets_direct_lnHR_bmi.rds")
saveRDS(log(HR_nonhdl_direct), "data/p5/offsets_direct_lnHR_nonhdl.rds")
saveRDS(log(HR_smok_init_direct), "data/p5/offsets_direct_lnHR_smok_init.rds")
saveRDS(log(HR_smok_cess_direct), "data/p5/offsets_direct_lnHR_smok_cess.rds")
saveRDS(log(HR_statins_total), "data/p5/offsets_total_lnHR_statins.rds")
saveRDS(log(HR_ah_total), "data/p5/offsets_total_lnHR_ah.rds")
saveRDS(log(HR_smok_init_total), "data/p5/offsets_total_lnHR_smok_init.rds")
saveRDS(log(HR_smok_cess_total), "data/p5/offsets_total_lnHR_smok_cess.rds")

###
### Effects of modifiable risk factors, on the other modifiable risk factors (used in intervention layer and rshiny)
###

### Most have already been derived in order to estimate the direct effects as the difference between
### total and indirect effects.

### Some have not been derived yet
bmi_to_nonhdl <- 0.2*(1/4.5)*2.348 # 0.1043556 (section 2.5)
bmi_to_sbp <- 0.7 # (section 2.4)

saveRDS(bmi_to_nonhdl, "data/p5/slope_bmi_to_nonhdl.rds")
saveRDS(bmi_to_sbp, "data/p5/slope_bmi_to_sbp.rds")
saveRDS(nonhdl_to_sbp, "data/p5/slope_nonhdl_to_sbp.rds")
saveRDS(smok_cess_to_bmi, "data/p5/slope_smok_cess_to_bmi.rds")
saveRDS(smok_cess_to_nonhdl, "data/p5/slope_smok_cess_to_nonhdl.rds")
saveRDS(smok_cess_to_sbp, "data/p5/slope_smok_cess_to_sbp.rds")
saveRDS(smok_init_on_bmi, "data/p5/slope_smok_init_to_bmi.rds")
saveRDS(smok_init_on_nonhdl, "data/p5/slope_smok_init_to_nonhdl.rds")
saveRDS(smok_init_on_sbp, "data/p5/slope_smok_init_to_sbp.rds")


###
### Manual intervention (re-do calculation from exemplar)
###

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

### Get HR
HR_combined <- HR_sbp_direct^-7.5*HR_bmi_direct^-5*HR_nonhdl_direct^-1
HR_combined

### Convert risk from exemplar
convert_risk_HR(.1675, HR_combined)
