###
### Function to create a numeric variable which is when an individual is "on" and "off" treatment
### An individual will assumed to be "off", 180 days after their last prescription.
###
### Steps will be to
### 1) Read in cohort and all prescriptions (prescriptions have been extracted and saved in extract_presciptions.R)
### 2) For individuals with a prescription, identify times where a "change" happens
### 3) Split survival dataset, and create a predictor for on/off the medication
###
### NB: These survival times will be used for models 1 and 2.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Load paralleisation packages
library(foreach)
library(doParallel)
library(doFuture)

### Extract chain seed and gender from command line
args <- commandArgs(trailingOnly = T)
med <- args[1]
burnout <- as.numeric(180)
print(paste("medication = ", med))
print(paste("burnout = ", burnout))

### Extract cohorts
## This contains cvd_time and cvd_indicator
cohort <- readRDS("data/cohort_pui.rds")
cohort <- dplyr::select(cohort, patid, cvd_time, cvd_indicator, fup_start)
print(paste("cohort", Sys.time()))
str(cohort)

### Read in the db.qry of the prescriptions of interest (contains all prescription dates of interest)
db_qry <- readRDS(paste("data/db_qry_", med, ".rds", sep = ""))

### Split cohort into smaller chunks
cohort_split <- split(cohort, cut(seq_len(nrow(cohort)), 10, labels = FALSE))

### Split db_qry
db_qry_split <- lapply(1:10, function(x) {db_qry[!is.na(fastmatch::fmatch(db_qry$patid, cohort_split[[x]]$patid)), ]})
str(db_qry_split)

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
### 1) Split survival times based on med_status_times
### 2) Define approriate medication status for each interval
### Writing a function that can be utilised with dplyr::group_modify and dplyr::group_map
get_cut_surv_times <- function(data, id, mst){
  
  #   id <- 1028382820230
  #   data <- subset(cohort[,c("patid", "cvd_time", "cvd_indicator")], patid == id)
  #   mst <- med_status_times
  
#   id <- "1028382820230"
#   data = subset(cohort, patid == id)
#   .f = get_cut_surv_times, 
#   mst = med_status_times
  
  ### Reduce mst to just patient of interest
  mstpat <- mst[!is.na(fastmatch::fmatch(mst$patid, id)), c("times_adj", "med_status")]
  
  ### Get splits
  splits <- mstpat$times_adj
  
  ### Apply survsplit
  data_split <- survival::survSplit(data = data, cut = splits, end = "cvd_time", event = "cvd_indicator")
  data_split <- data_split[,c("tstart", "cvd_time", "cvd_indicator")]
  
  ### Merge by tstat and times_adj
  data_split <- merge(data_split, mstpat[, c("times_adj", "med_status")], 
                      by.x = c("tstart"), 
                      by.y = c("times_adj"), 
                      all.x = TRUE)
  
  ### If the first row is NA (i.e. didn't have a "change on day 0", which is the case for most individuals)
  ### Set the appropriate value
  if (is.na(data_split$med_status[1])){
    ## If all data points happen prior to time zero, individual is off medication at baseline 
    ## (note we have a data point for individual transitioning off at the end)
    ### If individualss have all data points after time zero, they are also off medication at baseline
    ### (first data point is always first prescription)
    if (0 > max(mstpat$times_adj) | 0 < min(mstpat$times_adj)){
      data_split$med_status[1] <- 0
    } else {
      ## Else, individual has a data points before and after time zero, and could be on or off at baseline
      data_split$med_status[1] <- mstpat$med_status[max(which(mstpat$times_adj < 0))]
    } 
  }
  
  return(tibble::tibble(data_split))
  
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
  med_status_times <- dplyr::group_modify(.data = df, .f = get_change_status_times, burnout = burnout)
  med_status_times <- dplyr::arrange(med_status_times, patid, times)
  
  ### Merge med_status_times with cohort, so we can amend status change times to be in line with cvd_time, relative to the index date
  med_status_times <- merge(med_status_times, cohort[, c("patid", "fup_start")], by.x = "patid", by.y = "patid", all.x = TRUE)
  med_status_times$times_adj <- med_status_times$times - as.numeric(med_status_times$fup_start)
  
  ### Add medication status indicator
  med_status_times <- 
    dplyr::group_by(med_status_times, patid) |>
    dplyr::mutate(med_status = dplyr::case_when((dplyr::row_number() %% 2) == 1 ~ 1, 
                                                (dplyr::row_number() %% 2) == 0 ~ 0))
  
  ### Clear space
  print(paste("med_status_times", Sys.time()))
  str(med_status_times)
  
  ### Get patids of patients from the db_qry
  patids_df <- unique(med_status_times$patid)
  print(paste("patids_pres", Sys.time()))
  str(patids_df)
  
  ### Group cohort by patid
  cohort <- dplyr::group_by(cohort, patid) |>
    dplyr::select(-fup_start)
  
  ### Reduce to relevant individuals
  cohort <- cohort[!is.na(fastmatch::fmatch(cohort$patid, patids_df)), ]
  
  ### Apply function to these individuals
  cohort_split_times <- dplyr::group_modify(.data = cohort, 
                                            .f = get_cut_surv_times, 
                                            mst = med_status_times)
  print(paste("cohort_split_times", Sys.time()))
  
  return(cohort_split_times)
  
}

### Run function to get split times for people with statin/ah prescriptions during followup
print(paste("start conversion", Sys.time()))

### Apply this functions over the 10 splits cohort datasets and db queries
cl <- parallel::makeCluster(11)
doParallel::registerDoParallel(cl)
foreach::getDoParWorkers()
cohort_pres_split_times <- (foreach(x = 1:10, .combine = list, .multicombine = TRUE, 
                                    .packages = c("dplyr", "fastmatch"), 
                                    .export = c("get_cut_surv_times", "get_change_status_times")) %dopar% {
                                      convert_data(df = db_qry_split[[x]], cohort = cohort_split[[x]])
                                    })
stopCluster(cl)
print(paste("end conversion", Sys.time()))
str(cohort_pres_split_times)

### Recombine
print(paste("Recombine", Sys.time()))
cohort_pres_split_times <- do.call("rbind", cohort_pres_split_times) |>
  dplyr::arrange(patid, tstart)
str(cohort_pres_split_times)

###
### For individuals with no prescriptions, just need to change variables names and set med_status = 0
###

### Get patids of patients that have prescriptions
patids_pres <- unique(cohort_pres_split_times$patid)

### Get cohort of individuals with no prescriptions
cohort_nopres_split_times <- cohort[is.na(fastmatch::fmatch(cohort$patid, patids_pres)), ]

### Create new required variables so structure is same as cohort_pres_split_times
cohort_nopres_split_times <- 
  dplyr::select(cohort_nopres_split_times, patid, cvd_time, cvd_indicator) |> 
  dplyr::mutate(tstart = 0, 
                med_status = 0)
str(cohort_nopres_split_times)

###
### Combine the two datasets
###
cohort_split_times <- rbind(cohort_pres_split_times, cohort_nopres_split_times)
cohort_split_times <- data.frame(cohort_split_times) |>
  dplyr::arrange(patid, tstart)


### Save data
saveRDS(cohort_split_times, paste("data/cohort_split_times_", med, "_burnout", burnout, ".rds", sep = ""))

### Check this has same number of people as cohort
testthat::expect_equal(length(unique(cohort_split_times$patid)), length(unique(cohort$patid)))
testthat::expect_equal(sum(cohort_split_times$cvd_indicator), sum(cohort$cvd_indicator))
print(paste("FINISHED", Sys.time()))