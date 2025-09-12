###
### This program will get the vector of variables that is interacted with age in the imputation model
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

library(mice)
library(survival)
library(rms)

### Get vector of variables for male and female models
for (gender in 1:2){
  print(paste("GENDER = ", gender))
  
  ### Extract cohort
  cohort <- readRDS("data/extraction/cohort_baseline/cohort_var.rds")
  cohort <- data.table::as.data.table(cohort)
  
  ### Reduce to gender of interest
  cohort <- subset(cohort, gender == get("gender", pos = 1))
  # cohort <- cohort[1:100000, ]
  
  ###
  ### Get variable names
  ###
  var.all <- colnames(cohort)
  
  ### Save the vector of all variables in the dataset
  saveRDS(var.all, "data/extraction/cohort_baseline/var_all.rds")
  
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
               "sbp_var",  
               "cholhdl_ratio",
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
  saveRDS(var.imp, paste("data/extraction/cohort_baseline/var_imp", gender, ".rds", sep = ""))
  
  ### Remove variables that will be interacted with the spline of age (and remove age iteslf)
  inter.age.rcs <- var.imp[!grepl("age|bmi|sbp|cholhdl|IMD", var.imp)]
  print(inter.age.rcs)
  saveRDS(inter.age.rcs, paste("data/extraction/cohort_baseline/var_imp_inter_age_rcs", gender, ".rds", sep = ""))
  
}
