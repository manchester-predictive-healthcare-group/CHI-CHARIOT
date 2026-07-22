###
### Program to combine the interval censored survival times for statin and antihypertensive use during followup, 
### created in p7.1_augment_medication_status. The files created in p7.1 split at interval when an individual changes
### treatment. We combine these into a single dataset, for when an individual changes either treatment (statins or antihypertensives).
### Note, an individual may change treatment status with respect to both on the same day.
###
### NB: These survival times will be used for models 3 and 4.
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
cvd_split_times_statins <- readRDS(paste("data/cohort_split_times_augmented_", "statins", "_burnout", burnout, ".rds", sep = ""))
cvd_split_times_ah <- readRDS(paste("data/cohort_split_times_augmented_", "antihypertensives", "_burnout", burnout, ".rds", sep = ""))

### Extract cohort (used to define list of variables for adjustment)
cohort <- readRDS("data/cohort_pui.rds")
cohort <- data.table::as.data.table(cohort)

###
### Start by combining statins and antihypertensives
###

###
### Create treatment indicator so we can combine datasets
## Statins
cvd_split_times_statins <- dplyr::mutate(cvd_split_times_statins, treatment = "statins") |>
  dplyr::select(-med_status)

## Anti-hypertensives
cvd_split_times_ah <- dplyr::mutate(cvd_split_times_ah, treatment = "ah") |>
  dplyr::select(-med_status)

###
### Write a function to get split survival times, with intervals for change in either medication status
### Function written to be applied in conjunction with dplyr::group_modify
###
get_cut_surv_times <- function(data, id){
    # id <-1000065520274
    # id <- 1000003520274
  # id <- 1000002920274
  #   data <- subset(temp, patid == id)
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
  
  ### Note, med_status_adj = 0 at start, and will take values of -1 or 1 depending on if an individual
  ### stops (-1) or initiates (1) treatment from that point forward. med_status just = 0 if individual off treatment, 
  ### and = 1 if individual is on treatment
  
  ## Get initial status, which is always zero
  status_adj_statins <- pat_statins$med_status_adj[1]
  status_adj_ah <- pat_ah$med_status_adj[1]
  
  ### Get med_status_adj for every time med status changes.
  
  ### If med_status has changed from 0 to 1, or from 1 to 0 (i.e., the time in tstart_comb is a time when medication
  ### status changed, we use the new value of med_status_adj
  ### If med_status has not changed, we make no change
  ### Note, we know status only changes if tstart_comb[i] is in pat_statins$tstart or pat_ah$tstart respectively
  if (length(tstart_comb) > 1){
    
    ## Statins
    j <- 2
    for (i in 2:length(tstart_comb)){
      if (tstart_comb[i] %in% pat_statins$tstart){
        status_adj_statins <- 
          append(status_adj_statins, 
                 pat_statins$med_status_adj[j])
        j <- j + 1
      } else {
        status_adj_statins <- append(status_adj_statins, status_adj_statins[i-1])
      }
    }
    
    ### AH
    j <- 2
    for (i in 2:length(tstart_comb)){
      if (tstart_comb[i] %in% pat_ah$tstart){
        status_adj_ah <- 
          append(status_adj_ah, 
                 pat_ah$med_status_adj[j])
        j <- j + 1
      } else {
        status_adj_ah <- append(status_adj_ah, status_adj_ah[i-1])
      }
    }
    
  }

  ### Create data frame with split times
  out <- data.frame("tstart" = tstart_comb,
                    "cvd_time" = cvd_time_comb,
                    "med_status_statins" = status_adj_statins,
                    "med_status_ah" = status_adj_ah)
  
  ### Return as a tibble
  return(tibble::tibble(out))
  
}

### Combine the split survival time datasets, with an indicator for whether treatment is statins or antihypertensives
### This step is purely so the input dataset will conform with dplyr:group_modify
cvd_split_times <- rbind(cvd_split_times_statins, cvd_split_times_ah) |>
  dplyr::arrange(patid, treatment, tstart, cvd_time) |>
  dplyr::group_by(patid)

### Reduce cohort to just cvd_time, cvd_indicator and patid
cohort_reduced <- cohort[,c("patid", "cvd_time", "cvd_indicator")]

###
### Write a function to layer the medication status times using above functions, then merge with cohort to add back in cvd_indicator
###
convert_data3 <- function(df){
  
  ### Layer the event times
  cohort_split_times <- dplyr::group_modify(.data = df, 
                                            .f = get_cut_surv_times) |>
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
print(paste("start conversion", Sys.time()))

cl <- parallel::makeCluster(11)
doParallel::registerDoParallel(cl)
foreach::getDoParWorkers()
cohort_split_times <- (foreach(x = 1:10, .combine = list, .multicombine = TRUE, 
                               .packages = c("dplyr", "fastmatch"), 
                               .export = c("get_cut_surv_times")) %dopar% {
                                 convert_data3(df = cvd_split_times_cut[[x]])
                               })
stopCluster(cl)
print(paste("end conversion", Sys.time()))
str(cohort_split_times)

### Combine into single dataset and sort
cohort_split_times <- do.call("rbind", cohort_split_times) |>
  dplyr::arrange(patid, tstart, cvd_time)
str(cohort_split_times)

### Save adjsuted survival times, as it may take a while to run
saveRDS(cohort_split_times, paste("data/cohort_split_times_augmented_burnout", burnout, ".rds", sep = ""))
print("FINISHED")