###
### Sample size calculation for model
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
R.func.sources = list.files("Aurum_Jun2021_extract/R", full.names = TRUE)
sapply(R.func.sources, source)

###
### Write a function to get required number of variables for female or male models
###
get_nparms <- function(gender){
  
  ### Read in model 
  fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))

  return(length(fit$coefficients))
  
}

### Calculate number of parameters to be estimated
nparms_female <- get_nparms(2)
nparms_male <- get_nparms(1)
print("The number of coefficients for female model is:")
nparms_female
print("The number of coefficients for male model is:")
nparms_male



### Write a function to do the sample size calculation for female model or male model
est_sample_size <- function(gender.var, nparms, cvd.rate){
  
  ### Read in 1 of the imputed development datasets
  gender_char <- c("male", "female")[gender.var]
  cohort <- readRDS(paste("data/p4/dfs_devel_", gender_char, ".rds", sep = ""))[[1]]
  
  ### Calculate required inputs for sample size calculation
  E <- sum(cohort$cvd_indicator)
  n <- nrow(cohort)
  lnLnull <- E*log(E/n) + (n-E)*log(1 - E/n)
  maxR2cs <- 1 - exp(2*lnLnull/n)
  R2cs <- 0.15*maxR2cs
  mean_fup <- mean(cohort$cvd_time)/365.25
  
  ### Calculate pmsampsize according to an R2 NAGELKERKE of 0.15
  sampsize <- pmsampsize::pmsampsize(type = "s", 
                                     nagrsquared = 0.15, 
                                     parameters = nparms, 
                                     shrinkage = 0.99, 
                                     rate = cvd.rate, 
                                     meanfup = mean_fup,
                                     timepoint = 10)
  
  return(sampsize)
  
}

### Estimate required sample sizes
sampsize_female <- est_sample_size(gender.var = 2, nparms = nparms_female, cvd.rate = 4.91/1000)
sampsize_male <- est_sample_size(gender.var = 1, nparms = nparms_male, cvd.rate = 6.54/1000)
print("The number of observations required for female model is:")
sampsize_female
print("The number of observations required for male model is:")
sampsize_male


