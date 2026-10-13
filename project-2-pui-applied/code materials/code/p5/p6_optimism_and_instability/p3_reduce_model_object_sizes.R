###
### This program will take the models fitted in optimism loop, and save them as smaller objects, making 
### future code run more efficiently
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Define filepath to file directory system containing extracted data.
common.data.dir <- file.path("..", "..")

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define burnout and t_eval
t_eval <- 10*365.25
burnout <- as.numeric(180)

### Load packages
library(survival)

### Write a function to read in data, calculate offsets, remove the development dataset from the object, and save it
reduce_model_size <- function(gender, iter){
  
  ### Read in fit and bhaz
  fit <- readRDS(paste("data/p5/optimism_boot_cox_", gender, "_iter", iter, ".rds", sep = ""))
  bhaz <- readRDS(paste("data/p5/optimism_boot_cox_bhaz_", gender, "_iter", iter, ".rds", sep = ""))
  
  ### Calculate offset means (required for prediction)
  offset_means <- c("statins" = mean(fit$model[,"offset(offset_statins_timevar_lnHR)"]),
                    "ah" = mean(fit$model[,"offset(offset_ah_timevar_lnHR)"]),
                    "smoking1" = mean(fit$model[,"offset(offset_smoking_timevar_dummy1_lnHR)"]),
                    "smoking2" = mean(fit$model[,"offset(offset_smoking_timevar_dummy2_lnHR)"]))
  
  ### Get rid of residuals and linear predictors to reduce file size
  fit[names(fit) == "model"] <- "DELETED"
  fit[names(fit) == "y"] <- "DELETED"
  fit[names(fit) == "linear.predictors"] <- "DELETED"
  fit[names(fit) == "residuals"] <- "DELETED"
  
  ### Reduce the baseline hazard to one observation
  bhaz <- bhaz$hazard[max(which(bhaz$time <= t_eval))]
  
  ### Save
  saveRDS(fit, paste0("data/p5/optimism_boot_reduced_cox_", gender, "_iter", iter, ".rds"))
  saveRDS(bhaz, paste0("data/p5/optimism_boot_reduced_cox_bhaz_", gender, "_iter", iter, ".rds"))
  saveRDS(offset_means, paste0("data/p5/optimism_boot_reduced_cox_offset_means_", gender, "_iter", iter, ".rds"))
  
}

### Apply this function
lapply(c(2,1), function(y){
  lapply(1:500, function(x){
    print(paste("gender = ", y, "iter =", x))
    reduce_model_size(gender = y, iter = x)
  })
})