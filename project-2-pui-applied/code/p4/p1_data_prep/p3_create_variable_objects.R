###
### Create variable objects, such as:
### Lists of variables used in imputation model, final model, continuous/binary/polytomous, 
### variables that will be interacted with age
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

library(mice)
library(survival)
library(rms)

###
### Get vector of variables for male and female models, that will be used in the imputation
###
for (gender in 1:2){
  
  print(paste("GENDER = ", gender))
  
  ### Read in cohort
  cohort <- readRDS("data/p4/cohort_prototype3.rds")
  cohort <- data.table::as.data.table(cohort)
  
  # ### Reduce to gender of interest
  # cohort <- subset(cohort, gender == get("gender", pos = 1))
  # # cohort <- cohort[1:100000, ]
  
  ###
  ### Get variable names
  ###
  var.all <- colnames(cohort)
  
  ### Save the vector of all variables in the dataset
  saveRDS(var.all, "data/p4/var_all.rds")
  
  ###
  ### Define variables that will be used for imputation (for both genders)
  var.imp <- c("age",
               "ethnicity",
               "IMD",
               "hypertension",        
               "ra",
               "af", 
               "ckd", 
               "smi", 
               "fhcvd", 
               "migraine",
               "sle", 
               "diabetes", 
               "cortico",  
               "antipsy", 
               "bmi",     
               "sbp",  
               "nonhdl",
               "smoking",  
               "copd",
               "int_dis",
               "downs",
               "oral_cancer",
               "brain_cancer",
               "lung_cancer",
               "blood_cancer",
               "statins",  
               "antihypertensives") # Note cvd_indicator and cvd_time not added here, as they will be dealt with in the program that does the imputation
  
  ### For female cohort, add postnatal depression and pre_eclampsia
  ### For male cohort, add impotence
  if (gender == 1){
    var.imp <- c(var.imp,  "impotence")
  } else if (gender == 2){
    var.imp <- c(var.imp,  "pre_eclampsia", "postnatal_depression")
  }
  
  ### Save gender specific predictor variable lists
  saveRDS(var.imp, paste("data/p4/var_imp", gender, ".rds"))
  
  ### Remove variables that will be interacted with the spline of age (and remove age iteslf)
  inter.age.rcs <- var.imp[!grepl("age|bmi|sbp|chol|hdl|ldl|triglycerides|nonhdl|IMD", var.imp)]
  print(inter.age.rcs)
  saveRDS(inter.age.rcs, paste("data/p4/var_imp_inter_age_rcs", gender, ".rds"))
  
}

###
### Get vector of variables for male and female models, that will be used in the prediction model
###
for (gender in 1:2){
  
  print(paste("GENDER = ", gender))
  
  ### Read in cohort
  cohort <- readRDS("data/p4/cohort_prototype3.rds")
  cohort <- data.table::as.data.table(cohort)
  
  # ### Reduce to gender of interest
  # cohort <- subset(cohort, gender == get("gender", pos = 1))
  # # cohort <- cohort[1:100000, ]
  
  ###
  ### Get variable names
  ###
  var.all <- colnames(cohort)
  
  ###
  ### Define variables that will be used for model (for both genders)
  var.model <- c("age",
                 "ethnicity",
                 "IMD",
                 "bmi",
                 "sbp",
                 "nonhdl",
                 "smoking",
                 "hypertension",        
                 "ra",
                 "af", 
                 "ckd", 
                 "smi", 
                 "fhcvd", 
                 "migraine",
                 "sle", 
                 "diabetes", 
                 "cortico", 
                 "antipsy", 
                 "copd",
                 "int_dis",
                 "downs",
                 "oral_cancer",
                 "brain_cancer",
                 "lung_cancer",
                 "blood_cancer")
  
  ### For female cohort, add postnatal depression and pre_eclampsia
  ### For male cohort, add impotence
  if (gender == 1){
    var.model <- c(var.model,  "impotence")
  } else if (gender == 2){
    var.model <- c(var.model,  "pre_eclampsia", "postnatal_depression")
  }
  
  ### Save gender specific predictor variable lists
  saveRDS(var.model, paste("data/p4/var_model", gender, ".rds"))
  
  ### Remove variables that will be interacted with the spline of age (and remove age itself)
  inter.age.rcs <- var.model[!grepl("age|IMD|bmi|sbp|nonhdl|smoking", var.model)]
  print(inter.age.rcs)
  saveRDS(inter.age.rcs, paste("data/p4/var_model_inter_age_rcs", gender, ".rds"))
  
}

### Create a list of all variables
var.all <- c("age",
             "ethnicity",
             "IMD",
             "hypertension",        
             "ra",
             "af", 
             "ckd", 
             "smi", 
             "fhcvd", 
             "migraine",
             "sle", 
             "diabetes", 
             "cortico",  
             "antipsy", 
             "bmi",     
             "sbp",  
             "nonhdl",
             "smoking",  
             "copd",
             "int_dis",
             "downs",
             "oral_cancer",
             "brain_cancer",
             "lung_cancer",
             "blood_cancer",
             "pre_eclampsia", 
             "postnatal_depression",
             "impotence")

### Define continuous ones
var.cont <- c("IMD", "sbp", "bmi", "nonhdl")
saveRDS(var.cont, paste("data/p4/var_cont.rds", sep = ""))

### Define polytomous ones
var.poly <- c("ethnicity", "smoking", "diabetes")
saveRDS(var.poly, paste("data/p4/var_poly.rds", sep = ""))

### The rest are binary
var.binary <- var.all[!grepl("age|bmi|sbp|nonhdl|IMD|ethnicity|smoking|diabetes", var.all)]
saveRDS(var.binary, paste("data/p4/var_binary.rds", sep = ""))

###
### Define the hazard ratios for offsets in models
###

### THESE ARE NOW CALCULATED IN PROGRAM P9
# lnHR_statins <- log(0.75)
# saveRDS(lnHR_statins, "data/p4/offsets_lnHR_statins.rds")
# lnHR_ah <- log(0.73)
# saveRDS(lnHR_ah, "data/p4/offsets_lnHR_ah.rds")
# 
# lnHR_smoking_dummy1_total <- log(1.44)
# saveRDS(lnHR_smoking_dummy1_total, "data/p4/offsets_lnHR_smoking_dummy1_total.rds")
# lnHR_smoking_dummy2_total <- log(1/1.25)
# saveRDS(lnHR_smoking_dummy2_total, "data/p4/offsets_lnHR_smoking_dummy2_total.rds")
# 
# lnHR_smoking_dummy1_direct <- log(1.37)
# saveRDS(lnHR_smoking_dummy1_direct, "data/p4/offsets_lnHR_smoking_dummy1_direct.rds")
# lnHR_smoking_dummy2_direct <- log(1/1.32)
# saveRDS(lnHR_smoking_dummy2_direct, "data/p4/offsets_lnHR_smoking_dummy2_direct.rds")
# 
# lnHR_sbp <- log(1.028)
# saveRDS(lnHR_sbp, "data/p4/offsets_lnHR_sbp.rds")
# lnHR_bmi <- log(1.0174)
# saveRDS(lnHR_bmi, "data/p4/offsets_lnHR_bmi.rds")
# lnHR_nonhdl <- log(1/0.72)
# saveRDS(lnHR_nonhdl, "data/p4/offsets_lnHR_nonhdl.rds")