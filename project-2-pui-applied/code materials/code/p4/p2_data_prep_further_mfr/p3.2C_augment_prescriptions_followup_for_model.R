###
### This program will create variables for statin and antihypertensive use during follow-up.
### We will therefore assign a 0 at baseline, and a -1 or 1 depending on when an individual stops/starts treatment.
###
### In a previous program (p3.2A), we created variables for statin and antihypertensive use during follow-up.
### These variables were just 0 (off treatment) and 1 (on treatment). Given all the times at which individuals start/stop are the same,
### it will be computationally quicker to augment these previous variables, rather than extracting these afresh.
###
### The program in which the original variables were extracted in is: p3.2A_extract_prescriptions_followup_for_model.R. That program
### was written to be parallelised and run in batches, as the process was slow. These were then combined in:
### p3.2B_combine_prescriptions_followup_for_model.R.
###
### This program will read in the file created in p3.2A and p3.2B, and then create the variables needed for this analysis, as described above,
### by substracting or adding 1, depending on if the individual is on or off treatment at baseline.
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

### Extract chain seed and gender from command line
args <- commandArgs(trailingOnly = T)
med <- args[1]
burnout <- 180

print(paste("prescription = ", med))
print(paste("burnout = ", burnout))

### Read in cohort
cohort <- readRDS("data/p4/cohort_prototype3.rds") |>
  dplyr::select(patid, cvd_time, cvd_indicator, fup_start)

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
  as.data.frame()
print(str(cohort_split_times))

###
### Save
###
saveRDS(cohort_split_times, paste("data/p4/cohort_split_times_", med, "_burnout", burnout, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))