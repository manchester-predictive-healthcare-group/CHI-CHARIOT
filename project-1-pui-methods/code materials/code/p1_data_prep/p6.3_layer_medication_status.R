###
### Program to combine the interval censored survival times for statin and antihypertensive use during followup, 
### created in p6.2_extract_medication_status.R. The files created in p6.2 split at intervals when an individual changes
### treatment. We combine these into a single dataset, for when an individual changes either treatment (statins or antihypertensives).
### Note, an individual may change treatment status with respect to both on the same day.
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
burnout <- as.numeric(180)
print(paste("burnout = ", burnout))

### Read in the outcome data
cohort_split_times_statins <- readRDS(paste("data/cohort_split_times_", "statins", "_burnout", burnout, ".rds", sep = ""))
cohort_split_times_ah <- readRDS(paste("data/cohort_split_times_", "antihypertensives", "_burnout", burnout, ".rds", sep = ""))

### Create treatment indicator so we can combine datasets
## Statins
cohort_split_times_statins <- dplyr::mutate(cohort_split_times_statins, treatment = "statins")

## Anti-hypertensives
cohort_split_times_ah <- dplyr::mutate(cohort_split_times_ah, treatment = "ah")

###
### Create treatment indicator so we can combine datasets, and reduce dataset to subset for running on CSF

### Extract cohort (used to define list of variables for adjustment)
cohort <- readRDS("data/cohort_pui.rds")
cohort <- data.table::as.data.table(cohort)

### Merge with the reduced cohort
## Statins
cohort_split_times_statins <- merge(dplyr::select(cohort, c("patid")),
                                    cohort_split_times_statins,
                                    by.x = "patid",
                                    by.y = "patid") |>
  dplyr::mutate(treatment = "statins")

## Anti-hypertensives
cohort_split_times_ah <- merge(dplyr::select(cohort, c("patid")),
                               cohort_split_times_ah,
                               by.x = "patid",
                               by.y = "patid") |>
  dplyr::mutate(treatment = "ah")

###
### Write a function to get split survival times, with intervals for change in either medication status
### Fuction written to be applied in conjunction with dplyr::group_modify
###
get_cut_surv_times_stat_ah <- function(data, id){
  #     id <-1011847520626
  #     data <- subset(data.valid.split, patid == id)
  ### Break the input data dependent on whether its statins or ah
  pat_statins <- subset(data, treatment == "statins")
  pat_ah <- subset(data, treatment == "ah")
  
  ### If neither prescription, create output data frame with just one row
  ### With times over the follow up interval for that individual
  
  ### Get status change times for ah and statins
  times_ah <- c(0, pat_ah$cvd_time)
  times_statins <- c(0, pat_statins$cvd_time)
  
  ### Create a combined times vector, which is now what we will split cvd_event_times on
  tstart_comb <- c(times_ah, times_statins) |>
    unique() |> 
    sort() |> 
    head(-1)
  
  cvd_time_comb <- c(times_ah, times_statins) |>
    unique() |> 
    sort()
  cvd_time_comb <- cvd_time_comb[-1]
  
  ### Create vectors for on and off treatment
  
  ## Get initial status
  status_statins <- pat_statins$med_status[1]
  status_ah <- pat_ah$med_status[1]
  
  ## Get med_status for all combined change times
  if (length(tstart_comb) > 1){
    
    ## Statins
    for (i in 2:length(tstart_comb)){
      if (tstart_comb[i] %in% pat_statins$tstart){
        status_statins <- append(status_statins, 1-status_statins[i-1])
      } else {
        status_statins <- append(status_statins, status_statins[i-1])
      }
    }
    
    ### AH
    for (i in 2:length(tstart_comb)){
      if (tstart_comb[i] %in% pat_ah$tstart){
        status_ah <- append(status_ah, 1-status_ah[i-1])
      } else {
        status_ah <- append(status_ah, status_ah[i-1])
      }
    }
    
  }
  
  ### Create data frame with split times
  out <- data.frame("tstart" = tstart_comb,
                    "cvd_time" = cvd_time_comb,
                    "med_status_statins" = status_statins,
                    "med_status_ah" = status_ah)
  
  ### Return as a tibble
  return(tibble::tibble(out))
  
}

### Combine the split survival time datasets, with an indicator for whether treatment is statins or antihypertensives
### This step is purely so the input dataset will conform with dplyr:group_modify
cvd_split_times <- rbind(cohort_split_times_statins, cohort_split_times_ah) |>
  dplyr::arrange(patid, tstart, cvd_time) |>
  dplyr::group_by(patid)

### Reduce cohort to just cvd_time, cvd_indicator and patid
cohort_reduced <- cohort[,c("patid", "cvd_time", "cvd_indicator")]

###
### Write a function to layer the medication status times using above functions, then merge with cohort to add back in cvd_indicator
###
convert_data2 <- function(df){
  
  ### Layer the event times
  cohort_split_times <- dplyr::group_modify(.data = df, 
                                            .f = get_cut_surv_times_stat_ah) |>
    as.data.frame() |>
    dplyr::arrange(patid, tstart)
  print(paste("end conversion", Sys.time()))
  
  ### Merge with cohort_reduced by patid and cvd_time to add cvd_indicator back in
  ### Note cvd_time in cohort will always match the last interval in cohort_split_times
  cohort_split_times <- merge(cohort_split_times, cohort_reduced, by.x = c("patid", "cvd_time"), by.y = c("patid", "cvd_time"), all.x = TRUE)
  
  ### Observations with missing cvd_indicator, set to zero
  cohort_split_times <- dplyr::mutate(cohort_split_times, 
                                      cvd_indicator = dplyr::case_when(is.na(cvd_indicator) ~ 0,
                                                                       TRUE ~ cvd_indicator)
  )
  
  ### Return 
  return(cohort_split_times)
  
}

### Cut the split survival times into smaller chunks
### We can't split the cvd_split_times times data frame directly given multiple observations per patient,
### which we don't want to be split up

### First cut the cohort
cohort_cut <- split(cohort, cut(seq_len(nrow(cohort)), 10, labels = FALSE))

### Split cvd_split_times based on the cohort groupings
cvd_split_times_cut <- lapply(1:10, function(x) {cvd_split_times[!is.na(fastmatch::fmatch(cvd_split_times$patid, cohort_cut[[x]]$patid)), ]})
str(cvd_split_times_cut)

### Run function using parallelisation
cl <- parallel::makeCluster(11)
doParallel::registerDoParallel(cl)
foreach::getDoParWorkers()
cohort_split_times <- (foreach(x = 1:10, .combine = list, .multicombine = TRUE, 
                               .packages = c("dplyr", "fastmatch"), 
                               .export = c("get_cut_surv_times_stat_ah")) %dopar% {
                                 convert_data2(df = cvd_split_times_cut[[x]])
                               })
stopCluster(cl)
print(paste("end conversion", Sys.time()))
str(cohort_split_times)

### Combine into single dataset and sort
cohort_split_times <- do.call("rbind", cohort_split_times) |>
  dplyr::arrange(patid, tstart, cvd_time)
str(cohort_split_times)

### Save adjsuted survival times, as it may take a while to run
saveRDS(cohort_split_times, paste("data/cohort_split_times_burnout", burnout, ".rds", sep = ""))

### Test number of events matches that in cohort
testthat::expect_equal(sum(cohort_split_times$cvd_indicator), sum(cohort$cvd_indicator))

print("FINISHED")
