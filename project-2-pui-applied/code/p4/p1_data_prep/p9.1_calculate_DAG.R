### Code written by Bowen Jiang
# Goal: Estsimate the effects from the DAG presented in supplementary material, and used to drive the intervention layer of the CHARIOT model
# R script organization:
# 1. Define functions for effect size conversions (OR ↔ RR ↔ HR)
# 2. Estimate the direct effects of Modifiable Risk Factors (MRFs) on CVD,
#    including an internal consistency check for each MRF
# 3. Perform external consistency check for statin effects
# 4. Convert all MRF effect estimates to hazard ratios (HRs)

# the process to convert RR as OR for OR >1 , Two quadratic equations were derived based on the theorem presented in VanderWeele’s paper.
# Given values
u <- 0.3
w <- 0.00001

convert_RR_to_OR_gt1 <- function(RR, u, w) {
  a <- 1 - u
  b <- u - RR^2 * w
  c <- -RR^2 * (1 - w)
  
  D <- b^2 - 4 * a * c
  if (D < 0) {
    return(NA)  # No real solution
  } else {
    x1 <- (-b + sqrt(D)) / (2 * a)
    x2 <- (-b - sqrt(D)) / (2 * a)
    return(c(x1, x2))  # return both roots
  }
}

# the process to convert RR as OR for OR < 1 
convert_RR_to_OR_lt1 <- function(RR, u, w) {
  a <- 1 - w
  b <- w - RR^2 * u
  c <- -RR^2 * (1 - u)
  
  D <- b^2 - 4 * a * c
  if (D < 0) {
    return(NA)  # No real solution
  } else {
    x1 <- (-b + sqrt(D)) / (2 * a)
    x2 <- (-b - sqrt(D)) / (2 * a)
    return(c(x1, x2))  # return both roots
  }
}

# the process to convert HR to OR, first need to convert to RR, then convert RR to OR as shown in above. 
# the process of converting HR to RR, follow the formula in the Vander J paper

## HR > 1, w < p0 < p1 < u
convert_HR_to_RR_HRgt1 <- function(HR, w, u){
  top <- 1 - (1-w)^HR
  bot <- 1 - (1-u)^(1/HR) 
  out <- ((top/bot)*(u/w))^(1/2)
  return(out)
}

## HR < 1, w < p1 < p0 < u
convert_HR_to_RR_HRlt1 <- function(HR, w, u){
  top <- 1 - (1-u)^HR
  bot <- 1 - (1-w)^(1/HR) 
  out <- ((top/bot)*(w/u))^(1/2)
  return(out)
}


### OR to RR
###

## OR > 1, w < p0 < p1 < u
convert_OR_to_RR_ORgt1 <- function(OR, w, u){
  top <- OR*(u+OR-OR*u)
  bot <- 1 - w + OR*w
  out <- (top/bot)^(1/2)
  return(out)
}

## OR < 1, w < p1 < p0 < u
convert_OR_to_RR_ORlt1 <- function(OR, w, u){
  top <- OR*(w+OR-OR*w)
  bot <- 1 - u + OR*u
  out <- (top/bot)^(1/2)
  return(out)
}

###
### Function to do a numerical search to convert RR to HR
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


# Estimate the direct effect of SBP on CVD
# Clinical trials reported RR = 0.8 for a 10 mmHg reduction in SBP (Ettehad et al., 2016).
# Here we scale this effect to 1 mmHg reduction: (0.8)^(1/10).
RR_SBP <- (0.8)^(1/10)
RR_SBP   # Risk ratio per 1 mmHg reduction in SBP
# Convert the RR to OR using the RR < 1 conversion function
convert_RR_to_OR_lt1(RR_SBP, u, w)
# Express the effect per 1 mmHg increase in SBP (inverse of OR)
1 / 0.9741304 #1.026557



# For a 1 mmol/L reduction in Non-HDL cholesterol, 
# the reported total effect is HR = 0.8184262, 
# calculated as (0.78)^(1/1.24) based on the JBS3 model.
HR_NonHDL=0.8184262 #NonHDL total effect on CVD

# Convert HR to RR using the HR < 1 conversion function
convert_HR_to_RR_HRlt1(HR_NonHDL,w,u)
# (when further converted to OR, the value remains < 1).
RR_NonHDL= 0.8310591 # Total effect of Non-HDL on CVD (risk ratio)

# Convert RR to OR using the RR < 1 conversion function
convert_RR_to_OR_lt1(RR_NonHDL, u,w)
OR_NonHDL=0.8065872 #  Total effect of Non-HDL on CVD (odds ratio)

# Express as the effect per 1 mmol/L increase in Non-HDL 
1/0.8065872  # Equivalent OR = 1.239792

# Estimate the indirect effect of Non-HDL on CVD via SBP:
# A 1 mmol/L reduction in Non-HDL corresponds to a 
# 0.2749577 mmHg decrease in SBP (derived as (0.4258/1.0244)*0.6615).
NonHDL_SBP = (0.4258 / 1.0244) * 0.6615
NonHDL_SBP   # 0.2749577 mmHg reduction in SBP

# OR for CVD per 1 mmHg increase in SBP = 1.026577
# Therefore, the indirect effect of Non-HDL via SBP is:
OR_NonHDL_indirect = (1.026557)^(0.2749577)
OR_NonHDL_indirect   # 1.007233

# The direct effect of Non-HDL on CVD is then:
# (Total effect) / (Indirect effect via SBP)
OR_NonHDL_direct = 1.239792 / 1.007233
OR_NonHDL_direct

# Express as the effect per 1 mmol/L reduction in Non-HDL
1 / OR_NonHDL_direct # 0.812421

# Convert the direct effect of Non-HDL from OR to HR in two steps:
# Step 1: OR -> RR (apply conversion formula for OR < 1)
convert_OR_to_RR_ORlt1(0.812421, w, u)
# Result: RR_NonHDL_direct = 0.8362969
RR_NonHDL_direct <- 0.8362969

# Step 2: RR -> HR (numerical search, since RR-to-HR is not closed-form)
convert_RR_to_HR_lt1_numerical_search(RR_NonHDL_direct, w, u)
# Result: HR_NonHDL_direct = 0.823972
HR_NonHDL_direct <- 0.823972


# Internal consistency check:
# Total effect of Non-HDL on CVD, estimated through our DAG framework
# Direct effect of Non-HDL on CVD
1.230889   # OR_NonHDL_direct
# Indirect effect of Non-HDL via SBP:
#   Non-HDL increases SBP by 0.2749577 mmHg,
#   and the per-1 mmHg effect of SBP on CVD is OR = 1.026557
(1.026557)^(0.2749577)   # = 1.007233 (OR_NonHDL_indirect)
# Combine direct and indirect effects to obtain the total effect:
1.230889 * 1.007233      # = 1.239792 (OR_NonHDL_total)
# Express the total effect as the inverse (per 1 mmol/L reduction in Non-HDL):
1 / 0.8065872            # = 1.239792




# BMI differences between groups
# Overweight vs. Normal weight = 5 BMI units
# Obesity vs. Normal weight = 10 BMI units

log(1.12)/5  # = 0.02266574. Convert HR (1.12, Overweight vs. Normal) into log-risk increase per 1 BMI unit
log(1.22)/10 # = 0.01988509. Convert HR (1.22, Obesity vs. Normal) into log-risk increase per 1 BMI unit

(0.02266574 + 0.01988509)/2 # Average log-risk increase per 1 BMI unit across both comparisons
exp(0.02127541)              # Convert back from log scale to HR = 1.021503 (per 1 BMI unit increase)
1/1.021503                   # Inverse: HR for per 1 BMI unit decrease = 0.9789496

# Assign BMI direct effect HR
HR_BMI_direct = 0.9789496    

# Convert HR → RR → OR
convert_HR_to_RR_HRlt1(HR_BMI_direct, w, u)
RR_BMI_direct = 0.9806808
convert_RR_to_OR_lt1(RR_BMI_direct, u, w)
OR_BMI_direct = 0.9773420

1/OR_BMI_direct  # Express as OR for per 1 BMI unit increase = 1.023183


# Internal consistency :

# Compare BMI total effect estimated via DAGs with the total effect estimated from Causal Mediation analysis.
# Step 1: Estimate total effect of BMI on CVD using group differences
# BMI difference between groups:
#   Overweight vs. Normal weight = 5 BMI units
#   Obesity vs. Normal weight    = 10 BMI units
log(1.23)/5   # = 0.04140283. Convert HR=1.23 (Overweight vs. Normal) into log-risk increase per 1 BMI unit
log(1.43)/10  # = 0.03576744. Convert HR=1.43 (Obesity vs. Normal) into log-risk increase per 1 BMI unit
(0.04140283 + 0.03576744)/2  # Average log-risk increase per 1 BMI unit across both groups
exp(0.03858513)               # HR = 1.039339 per 1 BMI unit increase
1/1.039339                    # HR = 0.96215 per 1 BMI unit decrease
# Step 2: Assign BMI total effect HR
HR_BMI_total = 0.96215
HR_BMI_total
# Step 3: Convert HR → RR → OR
convert_HR_to_RR_HRlt1(HR_BMI_total, w, u)
RR_BMI_total = 0.96215  # (note: placeholder, same value carried forward here)
convert_RR_to_OR_lt1(RR_BMI_total, u, w)
OR_BMI_total = 0.9557408
# Step 4: Express effect as per 1 BMI unit increase
1/0.9557408  # = 1.0462 (OR for per 1 BMI unit increase)


# Total effect of BMI on CVD estimated through our "DAGs"

# Step 1: Direct effect of BMI on CVD (per 1 unit BMI increase)
1.023183  
# Step 2: Indirect effect via Non-HDL cholesterol
(1.230889)^0.1043451  # = 1.021913 (OR for 1 unit BMI increase mediated through Non-HDL)
# Step 3: Indirect effect via SBP
(1.026557)^(0.7)      # = 1.018517 (OR for 1 unit BMI increase mediated through SBP)
# Step 4: Combine direct and indirect effects to obtain total effect
1.023183 * 1.021913 * 1.018517  # = 1.064965 (OR for total effect of BMI on CVD)


# Estimating the direct effect of smoking cessation on CVD

# Step 1: Estimate changes in risk factors due to smoking cessation
0.2*0.66/4.5*2.348  # = 0.06887467 mmol/L increase in Non-HDL
# Step 2: Indirect effects of smoking cessation on CVD
# (a) via BMI
(1.023183)^0.66      # = 1.015241 (OR mediated through BMI)
# (b) via Non-HDL
(1.230889)^0.06887467  # = 1.014411 (OR mediated through Non-HDL)
# (c) via SBP
# Direct effect of SBP on CVD: OR = 1.026557 per 1 mmHg increase
(1.026557)^(-3.5)    # = 0.9123453 (protective effect via SBP reduction)
# Step 3: Combine all indirect pathways
1.015241 * 1.014411 * 0.9123453  # = 0.9395985 (total indirect effect)
# Step 4: Estimate direct effect
# Total effect from trials: OR = 0.8
# Direct effect = Total effect / Indirect effect
0.8 / 0.9395985  # = 0.8514275 (direct effect of smoking cessation on CVD)


# Total effect of smoking cessation on CVD estimated through our "DAGs"

# Step 1: Direct effect (from decomposition)
0.8514275  # direct effect of smoking cessation on CVD (OR)
# Step 2: Indirect effects via mediators
(1.023183)^(0.66)       # = 1.015241 (indirect effect via BMI)
(1.230889)^(0.06887467) # = 1.014411 (indirect effect via Non-HDL)
(1.026557)^(-3.5)       # = 0.9123453 (indirect effect via SBP)
# Step 3: Combine all pathways
0.8514275 * 1.015241 * 1.014411 * 0.9123453  
# = 0.8 (total effect of smoking cessation on CVD)


# Smoking cessation total effect estimation for two effect modifiers: CKD and T2D

# --- Step 1: Effect of SBP reduction on CVD for each subgroup ---

# CKD patients:
RR_SBP_CKD <- (0.84)^(1 / 10)   
# Trial evidence: RR = 0.84 for 10 mmHg SBP reduction
# Here: effect per 1 mmHg reduction in SBP
RR_SBP_CKD

# Convert RR -> OR for OR < 1
convert_RR_to_OR_lt1(RR_SBP_CKD, u, w)
OR_SBP_CKD = 0.979722  # OR per 1 mmHg reduction for CKD patients
1/OR_SBP_CKD           # OR per 1 mmHg increase = 1.020698

# Diabetes patients:
RR_SBP_diabetes <- (0.88)^(1 / 10)
# Trial evidence: RR = 0.88 for 10 mmHg SBP reduction
convert_RR_to_OR_lt1(RR_SBP_diabetes, u, w)
OR_SBP_diabetes = 0.9850869  # OR per 1 mmHg reduction for diabetes patients
1/OR_SBP_diabetes           # OR per 1 mmHg increase = 1.015139

# --- Step 2: Indirect effect of smoking cessation through SBP ---
(1.015139)^(-3.5)  # = 0.9487695 (diabetes)
(1.020698)^(-3.5)  # = 0.9308069 (CKD)

# --- Step 3: Total effect of smoking cessation on CVD ---
# Combine: direct effect (0.8514275) * indirect effects (BMI, Non-HDL, SBP)
0.8514275 * 1.015241 * 1.01441 * 0.9487695  # = 0.8319382 (diabetes)
0.8514275 * 1.015241 * 1.01441 * 0.9308069  # = 0.8161875 (CKD)


# Direct effect of smoking initiation on CVD estimation

# --- Step 1: Calculate indirect effects of smoking initiation ---
0.61                         # decrease in BMI (units)
0.2 * 0.61 / 4.5 * 2.348     # estimated effect on Non-HDL (≈ 0.06365689 mmol/L decrease)
0.61 * 0.7                   # estimated effect on SBP (≈ 0.427 mmHg decrease)
# Convert each pathway effect to odds ratios (ORs):
(1.023183)^(-0.61)           # indirect effect through BMI → OR = 0.9861171
(1.230889)^(-0.06365689)     # indirect effect through Non-HDL → OR = 0.9868632
(1.026557)^(-0.427)          # indirect effect through SBP → OR = 0.9888705
# Combine all indirect effects:
0.9861171 * 0.9868632 * 0.9888705  # total indirect effect OR = 0.9623319
# --- Step 2: Calculate direct effect ---
# Using total effect from literature (OR = 1.44), divide by total indirect effect:
1.44 / 0.9623319  # direct effect OR = 1.496365


# Total effect of smoking initiation estimated through our "DAGs"

1.496365  # Direct effect of smoking initiation on CVD (OR)
# --- Indirect effects via mediators ---
(1.023183)^(-0.61)        # BMI pathway → OR = 0.9861171
(1.230889)^(-0.06365689)  # Non-HDL pathway → OR = 0.9868632
(1.026557)^(-0.427)       # SBP pathway → OR = 0.9888705
# --- Combine direct and indirect effects ---
1.496365 * 0.9861171 * 0.9868635 * 0.9888705  
# Total effect ≈ 1.439566 (close to reported total effect of 1.44)


# Consistency check for statin total effect using DAG-based decomposition

# --- Indirect effects through mediators ---
(1.023183)^(0.33)       # Effect through BMI → OR = 1.007592
(1.230889)^(-1.3268)    # Effect through Non-HDL → OR = 0.7590975
(1.026557)^(-2.62)      # Effect through SBP → OR = 0.9336343

# --- Combine indirect effects to obtain total effect ---
1.007592 * 0.7590975 * 0.9336343  
# Total OR = 0.7141001

# --- Convert OR to RR for external comparison (JAMA meta-analysis RR = 0.75) ---
OR_statin = 0.7141001
convert_OR_to_RR_ORlt1(OR_statin, w, u)
RR_statin = 0.746848


# Section for calculating the HRs of Modifiable Risk Factors (MRFs)

# Direct effect of BMI on CVD
# Checked that the HR for BMI is consistent with the converted OR -> RR -> HR process
HR_BMI_direct # Original HR: 0.9789496
# Convert the direct effect OR back to RR
convert_OR_to_RR_ORlt1(OR_BMI_direct, w, u)
RR_BMI_direct = 0.9806808
# Convert the RR back to HR to verify consistency
convert_RR_to_HR_lt1_numerical_search(RR_BMI_direct, w, u)
HR_BMI_direct_Convert = 0.97895


# Converting the OR of the direct effect of NonHDL on CVD to HR
# Step 1: Start with the OR
OR_NonHDL_direct
# Step 2: Convert OR to RR
convert_OR_to_RR_ORlt1(OR_NonHDL_direct, w, u)
RR_NonHDL_direct = 0.8362969  # Resulting RR
# Step 3: Convert RR to HR
convert_RR_to_HR_lt1_numerical_search(RR_NonHDL_direct, w, u)
HR_NonHDL_direct = 0.823972  # Final HR


# Converting the OR of the direct effect of SBP on CVD to HR
# Step 1: Start with the OR
OR_SBP
# Step 2: Convert OR to RR
convert_OR_to_RR_ORlt1(OR_SBP, w, u)
RR_SBP_direct = 0.9779327  # Resulting RR
# Step 3: Convert RR to HR
convert_RR_to_HR_lt1_numerical_search(RR_SBP_direct, w, u)
HR_SBP_Direct = 0.975962  # Final HR

# Smoking cessation: Convert OR to HR
# Step 1: Direct effect OR
OR_Smoking_Cessation_Direct = 0.8514275
OR_Smoking_Cessation_total  = 0.8
# Step 2: Convert OR to RR
convert_OR_to_RR_ORlt1(OR_Smoking_Cessation_Direct, w, u)
RR_Smoking_Cessation_Direct = 0.8710619  # Resulting RR
# Step 3: Convert RR to HR
convert_RR_to_HR_lt1_numerical_search(RR_Smoking_Cessation_Direct, w, u)
HR_Smoking_Cessation_Direct = 0.860915  # Final HR


# Smoking initiation: Convert OR to HR
# Step 1: Direct effect OR
OR_Smoking_Initiation_Direct = 1.496365
# Step 2: Convert OR > 1 to RR
convert_OR_to_RR_ORgt1(OR_Smoking_Initiation_Direct, w, u)
RR_Smoking_Initiation_Direct = 1.419956  # Resulting RR
# Step 3: Convert RR to HR
convert_RR_to_HR_gt1_numerical_search(RR_Smoking_Initiation_Direct, w, u)
HR_Smoking_Initiation_Direct = 1.45827  # Final HR



