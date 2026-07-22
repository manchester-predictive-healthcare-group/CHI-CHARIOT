###
### Program to combine the event times created in p3.2A_extract_prescriptions_followup_for_model.R, along with data formatted the same
### way for all the individuals who didn't have a prescription.
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
med <- args[1]
burnout <- as.numeric(args[2])
num.groups <- as.numeric(args[3])
# med <- "statins"
# burnout <- 180
# num.groups <- 10
print(paste("medication = ", med))
print(paste("burnout = ", burnout))
print(paste("num.groups = ", num.groups))

### Extract cohort
cohort <- readRDS(file.path(common.data.dir, "Aurum_Jun2021_extract/data/extraction/cohort_baseline/cohort_var.rds"))

###
### Read in data for individuals who had a prescription and have had their event times split, and combine into a single dataset
###
cohort.pres.split.times <- lapply(1:num.groups, function(group.id) {
  readRDS(paste("data/cohort_split_times_", med, "_burnout", burnout, "_id", group.id, ".rds", sep = ""))}
  )

cohort.pres.split.times <- do.call("rbind", cohort.pres.split.times)
str(cohort.pres.split.times)

###
### For individuals with no prescriptions, just need to change variables names and set med_status = 0
###

### Read in the db.qry of the prescriptions of interest
db.qry <- readRDS(paste("data/db_qry_", med, ".rds", sep = ""))

### Get patids of patients that have prescriptions
patids.pres <- unique(db.qry$patid)
rm(db.qry)

### Get cohort of individuals with no prescriptions
cohort.nopres.split.times <- cohort[is.na(fastmatch::fmatch(cohort$patid, patids.pres)), ]

### Create new required variables so structure is same as cohort.pres.split.times
cohort.nopres.split.times <- 
  dplyr::select(cohort.nopres.split.times, patid, cvd_time, cvd_indicator) |> 
  dplyr::mutate(tstart = 0, 
                med_status = 0)
str(cohort.nopres.split.times)

###
### Combine the two datasets
###
cohort.split.times <- rbind(cohort.pres.split.times, cohort.nopres.split.times)
cohort.split.times <- data.frame(cohort.split.times)
str(cohort.split.times)

###
### Save
saveRDS(cohort.split.times, paste("data/cohort_split_times_", med, "_burnout", burnout, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))

print(paste("number of individuals with prescriptions = ", length(unique(cohort.pres.split.times$patid))))
print(paste("number of individuals without prescriptions = ", length(unique(cohort.nopres.split.times$patid))))