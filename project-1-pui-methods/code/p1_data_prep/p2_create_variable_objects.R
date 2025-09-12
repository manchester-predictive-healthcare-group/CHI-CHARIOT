###
### Create variable objects, such as:
### Lists of variables used in imputation model, final model, continuous/binary/polytomous, 
### variables that will be interacted with age
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("")

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

library(mice)
library(survival)
library(rms)

### All variables that may be used in the model
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
             "blood_cancer")

### Define continuous ones
var.cont <- c("IMD", "sbp", "bmi", "nonhdl")
saveRDS(var.cont, paste("data/var_cont.rds", sep = ""))

### Define polytomous ones
var.poly <- c("ethnicity", "smoking", "diabetes")
saveRDS(var.poly, paste("data/var_poly.rds", sep = ""))

### The rest are binary
var.binary <- var.all[!grepl("age|bmi|sbp|nonhdl|IMD|ethnicity|smoking|diabetes", var.all)]
saveRDS(var.binary, paste("data/var_binary.rds", sep = ""))




###
### Get vector of variables for male and female models, that will be used in the imputation
###
for (gender in 1:2){
  
  print(paste("GENDER = ", gender))
  
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
  saveRDS(var.imp, paste("data/var_imp", gender, ".rds"))
  
  ### Remove variables that will be interacted with the spline of age (and remove age iteslf)
  inter.age.rcs <- var.imp[!grepl("age|bmi|sbp|chol|hdl|ldl|triglycerides|nonhdl|IMD", var.imp)]
  print(inter.age.rcs)
  saveRDS(inter.age.rcs, paste("data/var_imp_inter_age_rcs", gender, ".rds", sep = ""))
  
}

###
### Get vector of variables for male and female models, that will be used in the prediction model
### Note, this may change between the four models, but will act as a starting point that can be augmented.
###
for (gender in 1:2){
  
  print(paste("GENDER = ", gender))
  
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
  saveRDS(var.model, paste("data/var_model", gender, ".rds"))
  
  ### Remove variables that will be interacted with the spline of age (and remove age itself)
  inter.age.rcs <- var.model[!grepl("age|IMD|bmi|sbp|nonhdl|smoking", var.model)]
  print(inter.age.rcs)
  saveRDS(inter.age.rcs, paste("data/var_model_inter_age_rcs", gender, ".rds", sep = ""))
  
}


