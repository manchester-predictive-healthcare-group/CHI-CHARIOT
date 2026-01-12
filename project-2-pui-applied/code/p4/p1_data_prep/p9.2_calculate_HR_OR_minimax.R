###
### Define treatment effects and save as reseearch object in R for later use
### Total effects of interventions (adjust for treatment drop-in)
### Direct effects of each MFR on CVD (used in intervention layer)
### Total effect of each MFR on the other MFRs, labelled as 'slopes' (used in intervention layer) 
### This program will apply the minimax conversion for between OR, RR and HR where relevant
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("")
getwd()

###
### First write functions to do the conversion
###

###
### Convert OR to RR

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
### Convert HR to RR

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
RRs_possible_ORs_gt1 <- sapply(ORs_possible_gt1, convert_HR_to_RR_HRgt1, w = 0.00001, u = 0.3)
RRs_possible_ORs_lt1 <- sapply(ORs_possible_lt1, convert_HR_to_RR_HRlt1, w = 0.00001, u = 0.3)

##############################################################
### Total effects used for adjusting for treatment drop-in ###
##############################################################

###
### Antihypertensives
### Want to convert to HR for adjusting for treatment drop in during model fitting
###

### We have RR from literature
RR_ah <- 0.74

### Now find the first value that matches the target
HR_ah <- HRs_possible_lt1[min(which(RRs_possible_HRs_lt1 > RR_ah))]
HR_ah

### Double check these conversions make sense, by converting the HR back to RR
convert_HR_to_RR_HRlt1(HR_ah, w = 0.00001, u = 0.3)
RR_ah

### Save this
saveRDS(log(HR_ah), "data/p4/offsets_lnHR_ah.rds")

###
### Statins
### Want to convert to HR for adjusting for treatment drop in during model fitting
###

### We have RR from literature
RR_statins <- 0.75

### Now find the first value that matches the target
HR_statins <- HRs_possible_lt1[min(which(RRs_possible_HRs_lt1 > RR_statins))]
HR_statins

### Double check these conversions make sense, by converting the HR back to RR
convert_HR_to_RR_HRlt1(HR_statins, w = 0.00001, u = 0.3)
RR_statins

### Save this
saveRDS(log(HR_statins), "data/p4/offsets_lnHR_statins.rds")


###
### Smoking
###

### We have the following effects from literature
OR_smoking_dummy1_total <- 1.44
OR_smoking_dummy2_total <- 1/1.25

### Want to convert OR to HR

### Need to first convert OR to RR
RR_smoking_dummy1_total <- convert_OR_to_RR_ORgt1(OR = OR_smoking_dummy1_total, w = 0.00001, u = 0.3)
RR_smoking_dummy1_total

RR_smoking_dummy2_total <- convert_OR_to_RR_ORlt1(OR = OR_smoking_dummy2_total, w = 0.00001, u = 0.3)
RR_smoking_dummy2_total

### Then convert RR to HR

### Now find the first value that matches the target
HR_smoking_dummy1_total <- HRs_possible_gt1[min(which(RRs_possible_HRs_gt1 > RR_smoking_dummy1_total))]
HR_smoking_dummy2_total <- HRs_possible_lt1[min(which(RRs_possible_HRs_lt1 > RR_smoking_dummy2_total))]
HR_smoking_dummy1_total
HR_smoking_dummy2_total

### Double check these conversions make sense, by converting the HR back to RR
convert_HR_to_RR_HRgt1(HR_smoking_dummy1_total, w = 0.00001, u = 0.3)
RR_smoking_dummy1_total
convert_HR_to_RR_HRlt1(HR_smoking_dummy2_total, w = 0.00001, u = 0.3)
RR_smoking_dummy2_total

saveRDS(log(HR_smoking_dummy1_total), "data/p4/offsets_lnHR_smoking_dummy1_total.rds")
saveRDS(log(HR_smoking_dummy2_total), "data/p4/offsets_lnHR_smoking_dummy2_total.rds")

######################
### Direct effects ###
######################

###
### Direct effects of each modifiable risk factor on CVD
###

### These have all been calculated by Bowen in his DAG and associated code, so I can read straight in
direct_OR_smoking_initiation <- 1.496817
direct_OR_smoking_cessation <- 0.8511499
direct_OR_nonhdl <- 1/0.812421
direct_OR_bmi <- 1/0.9768592
direct_OR_sbp <- 1/0.9741304

saveRDS(direct_OR_smoking_initiation, "data/p4/direct_OR_smoking_initiation.rds")
saveRDS(direct_OR_smoking_cessation, "data/p4/direct_OR_smoking_cessation.rds")
saveRDS(direct_OR_nonhdl, "data/p4/direct_OR_nonhdl.rds")
saveRDS(direct_OR_bmi, "data/p4/direct_OR_bmi.rds")
saveRDS(direct_OR_sbp, "data/p4/direct_OR_sbp.rds")

################################################################################
### Effects of modifiable risk factors, on the other modifiable risk factors ###
################################################################################

### These have been estimated in Bowens DAG, so can read straight in
bmi_to_nonhdl <- 0.1043451
bmi_to_sbp <- 0.7
nonhdl_to_sbp <- 0.2749577
smoking_cessation_to_bmi <- 0.66
smoking_cessation_to_nonhdl <- 0.068
smoking_cessation_to_sbp <- -3.5
smoking_initiation_to_bmi <- -0.61
smoking_initiation_to_nonhdl <- -0.063
smoking_initiation_to_sbp <- -0.427

saveRDS(bmi_to_nonhdl, "data/p4/slope_bmi_to_nondl.rds")
saveRDS(bmi_to_sbp, "data/p4/slope_bmi_to_sbp.rds")
saveRDS(nonhdl_to_sbp, "data/p4/slope_nonhdl_to_sbp.rds")
saveRDS(smoking_cessation_to_bmi, "data/p4/slope_smoking_cessation_to_bmi.rds")
saveRDS(smoking_cessation_to_nonhdl, "data/p4/slope_smoking_cessation_to_nonhdl.rds")
saveRDS(smoking_cessation_to_sbp, "data/p4/slope_smoking_cessation_to_sbp.rds")
saveRDS(smoking_initiation_to_bmi, "data/p4/slope_smoking_initiation_to_bmi.rds")
saveRDS(smoking_initiation_to_nonhdl, "data/p4/slope_smoking_initiation_to_nonhdl.rds")
saveRDS(smoking_initiation_to_sbp, "data/p4/slope_smoking_initiation_to_sbp.rds")