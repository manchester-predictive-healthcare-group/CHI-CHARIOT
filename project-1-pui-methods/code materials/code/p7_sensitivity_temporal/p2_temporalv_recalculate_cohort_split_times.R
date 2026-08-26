###
### We need to estimate counterfactual survival times for the cohort when the index
### date is defined as 1/2/3/4/5 post baseline/start of follow-up.
###
### When estimating counterfactual survival times, we do so under the treatment strategy 
### of "do nothing", adjuting for changes in antihypertensive and statin comparative to baseline.
###
### In order to do this, we need to get the intervention status based on the
### new index dates. If the intervention status has changed
### at the new index date, so do all the subsequent values when they change, which is what we
### are adjusting for.
###
### This program adjusts the interval censored data, relative to the new index dates
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- 180

### Define model = 4 (two component model) for these analyses
model <- 4

### Extract gender from command line
args <- commandArgs(trailingOnly = T)
gender <- as.numeric(args[1])
gender_char <- c("male", "female")[gender]
print(paste("gender = ", gender_char))

### Extract followup time from command line
t_fup <- round(365.25*as.numeric(args[2]))
print(paste("t_fup = ", t_fup))

### Read in the validation dataset
df_valid <- readRDS(paste("data/sens_temporalv_df_valid_", gender, "_", t_fup, ".rds", sep = ""))

### Read in split survival times
cohort_split_times <- readRDS(paste("data/cohort_split_times_augmented_burnout", burnout, ".rds", sep = ""))

### Reduce to patients in df_valid
print(paste("before reducing to validation, number of individuals = ", length(unique(cohort_split_times$patid))))
cohort_split_times <- cohort_split_times[!is.na(fastmatch::fmatch(cohort_split_times$patid, df_valid$patid)), ]
print(paste("after reducing to validation, number of individuals = ", length(unique(cohort_split_times$patid))))

### Reduce tstart and cvd_time by the number of years
cohort_split_times <- dplyr::mutate(cohort_split_times, 
                                    tstart = tstart - t_fup, 
                                    cvd_time = cvd_time - t_fup)

### Remove rows happening entirely prior to new index date
cohort_split_times <- subset(cohort_split_times, cvd_time > 0)

### Set tstart to zero if negative
cohort_split_times <- dplyr::mutate(cohort_split_times, tstart = 
                                      dplyr::case_when(tstart < 0 ~ 0,
                                                       TRUE ~ tstart))


### Write function.
### For each patient, if first value is a 0, do nothing. 
### If first value is 1, they are now on treatment at baseline, so take 1 off both "on" and "off"
### (individual should move between 0 (on treatment) and -1 (off treatment))
### If first value is -1,they are now off treatment at baseline, so add 1 on both "on" and "off"
### (individual should move between 0 (off treatment) and 1 (on treatment))
my_mutate <- function(df, id){
  
  if (df[1,"med_status_statins"] == 1){
    df <- dplyr::mutate(df, med_status_statins = med_status_statins - 1)
  } else if (df[1,"med_status_statins"] == -1){
    df <- dplyr::mutate(df, med_status_statins = med_status_statins + 1)
  }
  
  if (df[1,"med_status_ah"] == 1){
    df <- dplyr::mutate(df, med_status_ah = med_status_ah - 1)
  } else if (df[1,"med_status_ah"] == -1){
    df <- dplyr::mutate(df, med_status_ah = med_status_ah + 1)
  }
  
  return(df)
  
}

### Group data
cohort_split_times <- dplyr::arrange(cohort_split_times, patid, tstart) |>
  dplyr::group_by(patid)

### Apply function
print(paste("Start group modify", Sys.time()))
cohort_split_times <- dplyr::group_modify(.data = cohort_split_times, 
                                          .f = my_mutate) |>
  as.data.frame()
print(str(cohort_split_times))

### Save output
saveRDS(cohort_split_times, paste("data/sens_temporalv_cohort_split_times_augmented_", gender, 
                                  "_burnout", burnout, 
                                  "_t", t_fup, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))

