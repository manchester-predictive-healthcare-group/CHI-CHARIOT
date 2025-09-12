###
### We need to estimate counterfactual survival times for the cohort when the index
### date is defined as 1/2/3/4/5 post baseline/start of follow-up.
###
### When estimating counterfactual survival times, we do so under the treatment strategy 
### of "suppose the individual had not initiated any interventions", of which we consider
### antihypertensives, statins, or change in smoking status.
###
### In order to do this, we need to get the intervention status based on the
### new index dates. In prototype3, we only consider changes in the intervention status
### post index date. Therefore, if the intervention status has changed
### at the new index date, so do all the subsequent values when they change, which is what we
### are adjusting for.
###

### Because smoking is dependent on the value at baseline, which changes depending on the
### imputation, this gets derived relative to the baseline value after merging with the imputed
### data. Therefore the only variables that need to be adjusted are the statin and antihypertensive use

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

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
t_fup <- round(365.25*as.numeric(args[1]))
burnout <- 180

### Read in split survival times
cohort_split_times <- readRDS(paste("data/p4/cohort_split_times_smoking_statins_antihypertensives_burnout", burnout, ".rds", sep = ""))

### We only need to do this for individuals in the validation cohorts
### Read in patids
patids.valid.female <- readRDS(paste("data/p4/patids_valid_", 2, sep = ""))
patids.valid.male <- readRDS(paste("data/p4/patids_valid_", 1, sep = ""))
patids <- c(patids.valid.female, patids.valid.male)

### Reduce to these patients only
print(paste("before reducing to validation, number of individuals = ", length(unique(cohort_split_times$patid))))
cohort_split_times <- cohort_split_times[!is.na(fastmatch::fmatch(cohort_split_times$patid, patids)), ]
print(paste("after reducing to validation, number of individuals = ", length(unique(cohort_split_times$patid))))


### Reduce tstart and cvd_time by the number of years
cohort_split_times <- dplyr::mutate(cohort_split_times, 
                                    tstart = tstart - t_fup, 
                                    cvd_time = cvd_time - t_fup)

### Remove rows happening entirely prior to new index date
cohort_split_times <- subset(cohort_split_times, cvd_time > 0)
print(paste("after removing individuals censored before fup time, number of individuals = ", length(unique(cohort_split_times$patid))))

### set tstart to zero if negative
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
  
  if (df[1,"med_status_adj_statins"] == 1){
    df <- dplyr::mutate(df, med_status_adj_statins = med_status_adj_statins - 1)
  } else if (df[1,"med_status_adj_statins"] == -1){
    df <- dplyr::mutate(df, med_status_adj_statins = med_status_adj_statins + 1)
  }
  
  if (df[1,"med_status_adj_ah"] == 1){
    df <- dplyr::mutate(df, med_status_adj_ah = med_status_adj_ah - 1)
  } else if (df[1,"med_status_adj_ah"] == -1){
    df <- dplyr::mutate(df, med_status_adj_ah = med_status_adj_ah + 1)
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
saveRDS(cohort_split_times, paste("data/p4/validation_cohort_split_times_smoking_statins_antihypertensives_burnout", 
                                  burnout, "_t", 
                                  t_fup, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))

print(paste("final number of individuals = ", length(unique(cohort_split_times$patid))))
