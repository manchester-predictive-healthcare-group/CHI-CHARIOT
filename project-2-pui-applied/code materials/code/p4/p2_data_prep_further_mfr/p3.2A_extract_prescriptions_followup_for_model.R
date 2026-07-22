###
### Function to create a numeric variable which is when an individual is "on" and "off" treatment
### An individual will assumed to be "off", 180 days after their last prescription.

### Steps will be to
### 1) Read in cohort and all prescriptions (prescriptions have been extracted and saved in extract_presciptions.R)
### 2) For individuals with a prescription, identify times where a "change" happens
### 3) Split survival dataset, and create a predictor for on/off statins

### The query of prescriptions will be split into a number of smaller datasets before applying the operations from 2) and 3).
### This is for comptational reasons.
### These will then be combined in a seperate program, extract_medication_status2.R, along with data formatted the same
### way for all the individuals who didn't have a prescription.
### To do this we set number of groups (num.groups) to be 100, and cycle the task ID (group.id) from 1 - 100 
### They will need to be recombined after (done in program p3.2B).
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
burnout <- as.numeric(args[2])
num.groups <- as.numeric(args[3])
group.id <- as.numeric(args[4])

print(paste("medication = ", med))
print(paste("burnout = ", burnout))
print(paste("num.groups = ", num.groups))
print(paste("group.id = ", group.id))

### Extract cohorts
## This contains cvd_time and cvd_indicator
cohort <- readRDS(file.path(common.data.dir, "Aurum_Jun2021_extract/data/extraction/cohort_baseline/cohort_var.rds"))
cohort <- dplyr::select(cohort, patid, cvd_time, cvd_indicator, fup_start)
print(paste("cohort", Sys.time()))
str(cohort)

### Read in the db.qry of the prescriptions of interest (contains all prescription dates of interest)
db.qry <- readRDS(paste("data/db_qry_", med, ".rds", sep = ""))

###
### NB: I think this splitting/grouping of patients has been done in a roundabout way, as I was trying to avoid
### searching %in%, which has been solved by fastmatch::fmatch. It works, but see how this is done in p4 (prototype3 code)
### for a simpler way of how to do this.
###

### Split db.qry into different data frames, ensuring no patient is split over two dataframes
## Get patids of patients that have prescriptions
patids.pres <- unique(db.qry$patid)
## Assign a group variable for the patids, by which we will split the dataset
groups <- data.frame("patid" = patids.pres, "group" = cut(seq_along(patids.pres), breaks = num.groups, labels = FALSE))
## Merge db.qry with the group variable
db.qry <- merge(db.qry, groups, by.x = "patid", by.y = "patid")
## Reduce to group of interest
db.qry <- db.qry[!is.na(fastmatch::fmatch(db.qry$group, group.id)), ]
str(db.qry)

##############################################################
### Function to identify change of medication status times ###
##############################################################

###
### Write a function when an individuals changes from "on" to "off" or "off" to "on"
### Writing a function that can be utilised with dplyr::group_modify and dplyr::group_map
### burnout is number of days of no prescription to be considered off treatment
###
get_change_status_times <- function(data, id, burnout = 180){
  
  ### Sort by issuedate
  data <- dplyr::arrange(data, issuedate)
  
  ### Create empty vector
  change_status_times <- vector("numeric", 0)
  change_status_times <- append(change_status_times, data$issuedate[1])
  
  ### If more than 1 prescription
  if (nrow(data) > 1){
    
    for (i in 2:nrow(data)){
      ### Check if next issuedate is within burnout days of last issue date
      if (data$issuedate[i] - data$issuedate[i-1] <= burnout){
        ### Do nothing
        ### Else
      } else {
        ### Append last issuedate + burnout days to the change_status_times, for when last episode finished
        change_status_times <- append(change_status_times, data$issuedate[i-1] + burnout)
        ### Append the new issuedate to the change_status_times, for when new episode starts
        change_status_times <- append(change_status_times, data$issuedate[i])
      }
      
      ### If i is the last prescription, append prescription + burnout days, to denote patient finishing
      if (i == nrow(data)){
        change_status_times <- append(change_status_times, data$issuedate[i] + burnout)
      }
    }
    
  } else {
    ### If just 1 prescription
    change_status_times <- c(data$issuedate, data$issuedate + burnout)
  }
  
  ### Create output data.frame
  out <- data.frame("times" = change_status_times)
  
  return(tibble::tibble(out))
  
}


####################################################################################
### Function to split the survival times based on medication status change times ###
####################################################################################

### Function will
### 1) Split survival times based on med.status.times
### 2) Define approriate medication status for each interval
### Writing a function that can be utilised with dplyr::group_modify and dplyr::group_map
get_cut_surv_times <- function(data, id, mst){
  
  ### Reduce mst to just patient of interest
  mstpat <- mst[!is.na(fastmatch::fmatch(mst$patid, id)), c("times_adj", "med_status")]
  
  ### Get splits
  splits <- mstpat$times_adj
  
  ### Apply survsplit
  data.split <- survival::survSplit(data = data, cut = splits, end = "cvd_time", event = "cvd_indicator")
  data.split <- data.split[,c("tstart", "cvd_time", "cvd_indicator")]
  
  ### Merge by tstat and times_adj
  data.split <- merge(data.split, mstpat[, c("times_adj", "med_status")], 
                      by.x = c("tstart"), 
                      by.y = c("times_adj"), 
                      all.x = TRUE)
  
  ### If the first row is NA (i.e. didn't have a "change on day 0", which is the case for most individuals)
  ### Set the appropriate value
  if (is.na(data.split$med_status[1])){
    ## If all data points happen prior to time zero, individual is off medication at baseline 
    ## (note we have a data point for individual transitioning off at the end)
    ### If individualss have all data points after time zero, they are also off medication at baseline
    ### (first data point is always first prescription)
    if (0 > max(mstpat$times_adj) | 0 < min(mstpat$times_adj)){
      data.split$med_status[1] <- 0
    } else {
      ## Else, individual has a data points before and after time zero, and could be on or off at baseline
      data.split$med_status[1] <- mstpat$med_status[max(which(mstpat$times_adj < 0))]
    } 
  }
  
  return(tibble::tibble(data.split))
  
}


#################################################################################################
### Function to covert the data into split survival times, calling on the above two functions ###
#################################################################################################

### Write a function covert the data
convert_data <- function(df, cohort){

  ### Group data
  df <- dplyr::arrange(df, patid) |>
    dplyr::group_by(patid)
  
  ## Run function and arrange
  med.status.times <- dplyr::group_modify(.data = df, .f = get_change_status_times, burnout = burnout)
  med.status.times <- dplyr::arrange(med.status.times, patid, times)
  
  ### Merge med.status.times with cohort, so we can amend status change times to be in line with cvd_time, relative to the index date
  med.status.times <- merge(med.status.times, cohort[, c("patid", "fup_start")], by.x = "patid", by.y = "patid", all.x = TRUE)
  med.status.times$times_adj <- med.status.times$times - as.numeric(med.status.times$fup_start)
  
  ### Add medication status indicator
  med.status.times <- 
    dplyr::group_by(med.status.times, patid) |>
    dplyr::mutate(med_status = dplyr::case_when((dplyr::row_number() %% 2) == 1 ~ 1, 
                                                (dplyr::row_number() %% 2) == 0 ~ 0))
  
  ### Clear space
  print(paste("med.status.times", Sys.time()))
  str(med.status.times)
  
  ### Get patids of patients from the db.qry
  patids.df <- unique(med.status.times$patid)
  print(paste("patids.pres", Sys.time()))
  str(patids.df)
  
  ### Group cohort by patid
  cohort <- dplyr::group_by(cohort, patid) |>
    dplyr::select(-fup_start)
  
  ### Apply function to these individuals
  cohort.pres.split.times <- dplyr::group_modify(.data = cohort[!is.na(fastmatch::fmatch(cohort$patid, patids.df)), ], 
                                                 .f = get_cut_surv_times, 
                                                 mst = med.status.times)
  print(paste("cohort.pres.split.times", Sys.time()))
  
  return(cohort.pres.split.times)
  
}

### Run function
print(paste("start conversion", Sys.time()))
df.converted <- convert_data(df = db.qry, cohort)
print(paste("end conversion", Sys.time()))

### Save data
saveRDS(df.converted, paste("data/cohort_split_times_", med, "_burnout", burnout, "_id", group.id, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))