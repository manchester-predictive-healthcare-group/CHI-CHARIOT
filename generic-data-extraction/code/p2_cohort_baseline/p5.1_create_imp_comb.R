###
### Fit a standard cox model, do a sample size calculations, assess non-linear curves and interactions
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

### Extract chain seed and gender from command line
# args <- commandArgs(trailingOnly = T)
# gender <- as.numeric(args[1])

for (gender in c(2,1)){
  
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  
  ### Read in the first imputed dataset
  imp.comb <- readRDS(paste("data/extraction/cohort_baseline/mice_mids_", gender, "_", 1, ".rds", sep = ""))
  
  ### Create one large mids object
  for (i in c(2:10)){
    print(paste("COMBINE", i, Sys.time()))
    imp.new <- readRDS(paste("data/extraction/cohort_baseline/mice_mids_", gender, "_", i, ".rds", sep = ""))
    imp.comb <- ibind(imp.comb, imp.new)
  }
  
  ### Save the combined imputation object
  saveRDS(imp.comb, paste("data/extraction/cohort_baseline/imp_comb", gender, ".rds", sep = ""))
  print(paste("FINISHED", gender, Sys.time()))
        
}

print(paste("FINISHED", Sys.time()))