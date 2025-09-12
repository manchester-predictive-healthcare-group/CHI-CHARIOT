###
### This program will create variables for statin and antihypertensive use during follow-up, required for models 3 and 4.
### We will therefore assign a 0 at baseline, and a -1 or 1 depending on when an individual stops/starts treatment.
### i.e., the derived variable will represent medication status relative to status at baseline.
###
### In a previous program (p6.2), we created variables for statin and antihypertensive use during follow-up.
### These variables were just 0 (off treatment) and 1 (on treatment). Given all the times at which individuals start/stop are the same,
### it will be computationally quicker to augment these previous variables, rather than extracting these afresh.
###
### The program in which the original variables were extracted in is: p6.2_extract_medication_status.R. This got medication status for
### statins and antihypertensives separately. This program will read in the data created in programs p6.2, and augment this to be in 
### the desired format outlined in lines 2 - 4.
###
### NB: These survival times will be used for models 3 and 4.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("")
getwd()

### Extract chain seed and gender from command line
args <- commandArgs(trailingOnly = T)
med <- args[1]
burnout <- 180

print(paste("prescription = ", med))
print(paste("burnout = ", burnout))

### Read in split survival times based on medication
cohort_split_times <- readRDS(paste("data/cohort_split_times_", med, "_burnout", burnout, ".rds", sep = ""))

### Arrange and group by patid
cohort_split_times <- dplyr::arrange(cohort_split_times, patid, tstart) |>
  dplyr::group_by(patid)

### Write function.
### For each patient, if first value is a 0, do nothing. If first value is 1, take 1 off both "on" and "off"
my_mutate <- function(df, id){
  
  if (df[1,"med_status"] == 1){
    df <- dplyr::mutate(df, med_status_adj = med_status - 1)
  } else {
    df <- dplyr::mutate(df, med_status_adj = med_status)
  }
  
  return(df)
  
}

### Apply function
print(paste("Start group modify", Sys.time()))
cohort_split_times <- dplyr::group_modify(.data = cohort_split_times, 
                           .f = my_mutate) |>
  as.data.frame() |>
  dplyr::arrange(patid, tstart)
print(str(cohort_split_times))

###
### Save
###
saveRDS(cohort_split_times, paste("data/cohort_split_times_augmented_", med, "_burnout", burnout, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))