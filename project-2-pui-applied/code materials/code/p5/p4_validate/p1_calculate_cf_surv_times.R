###
### Program to estimate counterfactual survival times.
### These are survival times if individuals had of maintained treatment at baseline over the prediction interval
### This is done seperately for every imputed dataset.
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

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- 180
print(paste("burnout = ", burnout))

### Write a function to estimate the counterfactual survival times for an imputed dataset/model
est_cf_surv_times_wrapper <- function(gender, imp_num){
  
  ### Define gender_char
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char, imp_num, Sys.time()))
  
  ### Read in imputed dataset
  df_imp <- readRDS(paste0("data/p5/dfs_all_", gender_char, imp_num, ".rds"))
  
  ### Read in the outcome data, which contains the split survival times
  cohort_split_times <- readRDS(paste0("data/p5/cohort_split_times_smoking_statins_antihypertensives_burnout", burnout, ".rds"))
  
  ### Read in baseline hazard from the prediction model
  bhaz <- readRDS(paste0("data/p5/model_cox_bhaz_", gender, "_imp", imp_num, ".rds"))
  
  ### Estimate counterfactual survival times
  cf_surv_times <- estimate_cf_surv_times(df_imp, cohort_split_times, bhaz, type = "orig")
  
  ### Save the counterfactual survival times
  saveRDS(cf_surv_times, paste0("data/p5/cf_surv_times_", gender, "_imp", imp_num, ".rds"))
  
}

### Calculate counterfactual survival times
lapply(1:10, function(x) {est_cf_surv_times_wrapper(1, x)})
lapply(1:10, function(x) {est_cf_surv_times_wrapper(2, x)})
print(paste("FINISHED", Sys.time()))