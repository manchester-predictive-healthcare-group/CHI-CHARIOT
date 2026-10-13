###
### Program to combine the event times created in p3.1A
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Extract chain seed and gender from command line
args <- commandArgs(trailingOnly = T)
num.groups <- 100
burnout <- 180
print(paste("num.groups = ", num.groups))

### Read in cohort
cohort <- readRDS("data/p5/cohort_prototype3.rds") |>
  dplyr::select(patid, cvd_time, cvd_indicator, fup_start)

###
### Read in smoking data for individuals during follow-up
###
cohort.smoking.split.times <- lapply(1:num.groups, function(group.id) {
  readRDS(paste("data/p5/cohort_split_times_smoking_id", group.id, ".rds", sep = ""))}
)
cohort.smoking.split.times <- do.call("rbind", cohort.smoking.split.times)
print("SMOKING DONE")
str(cohort.smoking.split.times)

###
### Save
saveRDS(cohort.smoking.split.times, paste("data/p5/cohort_split_times_smoking.rds", sep = ""))
print(paste("FINISHED", Sys.time()))