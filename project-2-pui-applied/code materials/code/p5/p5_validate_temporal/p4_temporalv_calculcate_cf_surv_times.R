###
### This program calculates the counterfactual survival times for the cohort when the index
### date is defined as 1/2/3/4/5 post baseline/start of follow-up.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- 180

### Write a function to estimate the counterfactual survival times for an imputed dataset/model
est_cf_surv_times_wrapper <- function(gender, t_fup_int){
  
  ### Define gender_char
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender, "t_fup_int =", t_fup_int, Sys.time()))
  
  ### Read in first dataset
  df_imp <- readRDS(paste0("data/p5/df_valid_temporalv_", gender, "_", t_fup_int, ".rds"))
  
  ### Read in the outcome data, which contains the split survival times
  cohort_split_times <- readRDS(paste0("data/p5/temporalv_cohort_split_times_smoking_statins_antihypertensives_burnout", 
                                       burnout, "_t", 
                                       t_fup_int, ".rds"))
  
  ### Read in baseline hazard from the first prediction model
  bhaz <- readRDS(paste0("data/p5/model_cox_bhaz_", gender, "_imp", 1, ".rds"))
  
  ### Estimate counterfactual survival times
  cf_surv_times <- estimate_cf_surv_times(df_imp, cohort_split_times, bhaz, type = "orig")
  
  ### Save the counterfactual survival times
  saveRDS(cf_surv_times, paste0("data/p5/temporalv_cf_surv_times_", gender, "_t", t_fup_int, ".rds"))
  
}

### Run function
lapply(1:2, function(x){
  lapply(1:5, function(y){
    est_cf_surv_times_wrapper(gender = x, t_fup_int = y)
  })
})

print(paste("FINISHED", Sys.time()))