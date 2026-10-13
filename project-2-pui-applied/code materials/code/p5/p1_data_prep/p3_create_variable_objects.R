###
### Create variable objects, such as:
### Lists of variables used in imputation model, final model, continuous/binary/polytomous, 
### variables that will be interacted with age
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")

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
  cohort <- readRDS("data/p5/cohort_prototype3.rds")
  cohort <- data.table::as.data.table(cohort)
  
  # ### Reduce to gender of interest
  # cohort <- subset(cohort, gender == get("gender", pos = 1))
  # # cohort <- cohort[1:100000, ]
  
  ###
  ### Get variable names
  ###
  var.all <- colnames(cohort)
  
  ### Save the vector of all variables in the dataset
  saveRDS(var.all, "data/p5/var_all.rds")
  
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
  saveRDS(var.imp, paste0("data/p5/var_imp", gender, ".rds"))
  
  ### Remove variables that will be interacted with the spline of age (and remove age iteslf)
  inter.age.rcs <- var.imp[!grepl("age|bmi|sbp|chol|hdl|ldl|triglycerides|nonhdl|IMD", var.imp)]
  print(inter.age.rcs)
  saveRDS(inter.age.rcs, paste0("data/p5/var_imp_inter_age_rcs", gender, ".rds"))
  
}

###
### Get vector of variables for male and female models, that will be used in the prediction model
###
for (gender in 1:2){
  
  print(paste("GENDER = ", gender))
  
  ### Read in cohort
  cohort <- readRDS("data/p5/cohort_prototype3.rds")
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
  saveRDS(var.model, paste0("data/p5/var_model", gender, ".rds"))
  
  ### Remove variables that will be interacted with the spline of age (and remove age itself)
  inter.age.rcs <- var.model[!grepl("age|IMD|bmi|sbp|nonhdl|smoking", var.model)]
  print(inter.age.rcs)
  saveRDS(inter.age.rcs, paste0("data/p5/var_model_inter_age_rcs", gender, ".rds"))
  
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
saveRDS(var.cont, paste0("data/p5/var_cont.rds"))

### Define polytomous ones
var.poly <- c("ethnicity", "smoking", "diabetes")
saveRDS(var.poly, paste0("data/p5/var_poly.rds"))

### The rest are binary
var.binary <- var.all[!grepl("age|bmi|sbp|nonhdl|IMD|ethnicity|smoking|diabetes", var.all)]
saveRDS(var.binary, paste0("data/p5/var_binary.rds"))