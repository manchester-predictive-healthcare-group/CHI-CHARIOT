###
### This program will check that the split survival times, extracted for index dates
### 1/2/3/4/5 years post baseline, contains the same set of individuals in the cohorts
### that we will be testing validation in. For a given number of years post follow-up, 
### we have to exclude individuals that have had an event or been censored prior to this
### time.
###
### This check is sensible, because when we dervied the split survival times, this was
### done in a seperate program compared to where we defined the validation cohort, 
### and the exclusion criteria (for individuals censored prior to
### the new index date), was applied on the interval censored data, as opposed to the
### non-interval censored data.
###
### The split survival times were derived in :
### /project2/code/p4/p5_temporal_validation/p2_layer_cohort_split_times_fup.R
###
### The predictors at the relevant index dates were derived in:
### /project2/code/p4/p5_temporal_validation/p1.2_temporalv_create_individual_cohorts.R
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Source functions
R.func.sources = list.files(file.path(common.data.dir, "Aurum_Jun2021_extract/R"), full.names = TRUE)
sapply(R.func.sources, source)

### Function to run tests
runtest <- function(gender_in, fup){
  gender_in <- 2
  fup <- 1
  ### Define follow-up time
  t_fup <- round(365.25*fup)
  
  ### Split survival times
  cohort_split_times <- readRDS(paste("data/p4/validation_cohort_split_times_smoking_statins_antihypertensives_burnout", 
                                           180, "_t", 
                                      t_fup, ".rds", sep = ""))
  
  ### Read in cohort at baseline to get gender
  cohort_baseline <- readRDS("data/p4/cohort_prototype3.rds") |>
    dplyr::select(patid, gender)
  ### Merge and reduceto gender of interest
  cohort_split_times <- dplyr::left_join(cohort_split_times, cohort_baseline, by = dplyr::join_by(patid)) |>
    dplyr::filter(gender == gender_in)
  
  ### Get n
  npat_interval <- length(unique(cohort_split_times$patid))

  ### Reduce to these patients only
  df_valid <- readRDS(paste("data/p4/df_valid_temporalv_", gender_in, "_", t_fup, ".rds", sep = ""))
  npat_non_interval <- nrow(df_valid)
  
  testthat::expect_equal(npat_interval, npat_non_interval)
  
}

### ADD SOME TESTS THAT CVD_TIME MATCHES THAT FROM DF_VALID?
# runtest(2,1)
runtest(2,2)
runtest(2,3)
runtest(2,4)
runtest(2,5)
# runtest(1,1)
runtest(1,2)
runtest(1,3)
runtest(1,4)
runtest(1,5)
print("FINISHED")