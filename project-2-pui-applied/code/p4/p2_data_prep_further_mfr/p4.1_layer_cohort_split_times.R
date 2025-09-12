###
### Program to layer the split survival times from smoking, statins and antihypertensives,
### i.e. create interval censored data with intervals whenever an individual changes status with
### respect to smoking, statins or antihypertensives.
###
### This program is written to be parallelised, i.e. will do the layer for a subgroup of patients
### They will need to be recombined after (done in program p4.2).
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
group.id <- as.numeric(args[1])
burnout <- as.numeric(180)
num.groups <- 100
print(paste("burnout = ", burnout))
print(paste("group.id = ", group.id))
print(paste("num.groups = ", num.groups))

### Read in the outcome data
cvd.split.times.smoking <- readRDS(paste("data/p4/cohort_split_times_smoking.rds", sep = ""))
cvd.split.times.statins <- readRDS(paste("data/p4/cohort_split_times_statins_burnout", burnout, ".rds", sep = ""))
cvd.split.times.ah <- readRDS(paste("data/p4/cohort_split_times_antihypertensives_burnout", burnout, ".rds", sep = ""))

### Read in cohort
cohort <- readRDS("data/p4/cohort_prototype3.rds") |>
  dplyr::select(patid, cvd_time, cvd_indicator, fup_start) |>
  dplyr::arrange(patid)
# cohort <- readRDS(file.path(common.data.dir, "Aurum_Jun2021_extract/data/extraction/cohort_baseline/cohort_var.rds"))
# cohort <- dplyr::select(cohort, patid, cvd_time, cvd_indicator, fup_start)

### Split up patids
chunk.size <- nrow(cohort)/num.groups
patids.split <- split(cohort$patid, ceiling(seq_along(cohort$patid) / chunk.size))
# sum(unlist(lapply(patids.split, length))) 

### Reduce cohort split times
cvd.split.times.smoking <- cvd.split.times.smoking[!is.na(fastmatch::fmatch(cvd.split.times.smoking$patid, patids.split[[group.id]])), ]
cvd.split.times.statins <- cvd.split.times.statins[!is.na(fastmatch::fmatch(cvd.split.times.statins$patid, patids.split[[group.id]])), ]
cvd.split.times.ah <- cvd.split.times.ah[!is.na(fastmatch::fmatch(cvd.split.times.ah$patid, patids.split[[group.id]])), ]

### Just check same number of individuals in each (i.e. no-one has been missed)
testthat::expect_equal(length(unique(cvd.split.times.smoking$patid)), length(unique(cvd.split.times.statins$patid)))
testthat::expect_equal(length(unique(cvd.split.times.smoking$patid)), length(unique(cvd.split.times.ah$patid)))

###
### Start by combining statins and antihypertensives
###

###
### Create treatment indicator so we can combine datasets
## Statins
cvd.split.times.statins <- dplyr::mutate(cvd.split.times.statins, treatment = "statins") |>
  dplyr::select(-med_status)

## Anti-hypertensives
cvd.split.times.ah <- dplyr::mutate(cvd.split.times.ah, treatment = "ah") |>
  dplyr::select(-med_status)

## smoking
cvd.split.times.smoking <- dplyr::mutate(cvd.split.times.smoking, treatment = "smoking") |>
  dplyr::rename(med_status_adj = smoking) |>
  dplyr::select(-fup_start)

###
### Write a function to get split survival times, with intervals for change in either medication status
### Function written to be applied in conjunction with dplyr::group_modify
###
get_cut_surv_times_stat_ah_smoking <- function(data, id){

  ### Break the input data dependent on whether its statins or ah
  pat.statins <- subset(data, treatment == "statins")
  pat.ah <- subset(data, treatment == "ah")
  pat.smoking <- subset(data, treatment == "smoking")
  
  ### If neither prescription, create output data frame with just one row
  ### With times over the follow up interval for that individual
  
  ### Get status change times for ah and statins
  times.ah <- c(0, pat.ah$cvd_time)
  times.statins <- c(0, pat.statins$cvd_time)
  times.smoking <- c(0, pat.smoking$cvd_time)
  
  ### Create a combined times vector, which is now what we will split cvd_event_times on
  tstart.comb <- c(times.ah, times.statins, times.smoking) |>
    unique() |> 
    sort() |> 
    head(-1)
  
  cvd_time.comb <- c(times.ah, times.statins, times.smoking) |>
    unique() |> 
    sort()
  cvd_time.comb <- cvd_time.comb[-1]
  
  ### Create vectors for on and off treatment
  
  ### Note, med_status_adj = 0 at start, and will take values of -1 or 1 depending on if an individual
  ### stops (-1) or initiates (1) treatment from that point forward. med_status just = 0 if individual off treatment, 
  ### and = 1 if individual is on treatment
  
  ## Get initial status, which is always zero
  status_adj_statins <- pat.statins$med_status_adj[1]
  status_adj_ah <- pat.ah$med_status_adj[1]
  status_adj_smoking <- pat.smoking$med_status_adj[1]
  
  ### Get med_status_adj for every time med status changes.
  
  ### If med_status has changed from 0 to 1, or from 1 to 0 (i.e., the time in tstart.comb is a time when medication
  ### status changed, we use the new value of med_status_adj
  ### If med_status has not changed, we make no change
  ### Note, we know status only changes if tstart.comb[i] is in pat.statins$tstart or pat.ah$tstart respectively
  if (length(tstart.comb) > 1){
    
    ## Statins
    j <- 2
    for (i in 2:length(tstart.comb)){
      if (tstart.comb[i] %in% pat.statins$tstart){
        status_adj_statins <- 
          append(status_adj_statins, 
                 pat.statins$med_status_adj[j])
        j <- j + 1
      } else {
        status_adj_statins <- append(status_adj_statins, status_adj_statins[i-1])
      }
    }
    
    ### AH
    j <- 2
    for (i in 2:length(tstart.comb)){
      print(paste("i=",i))
      if (tstart.comb[i] %in% pat.ah$tstart){
        status_adj_ah <- 
          append(status_adj_ah, 
                 pat.ah$med_status_adj[j])
        j <- j + 1
      } else {
        status_adj_ah <- append(status_adj_ah, status_adj_ah[i-1])
      }
    }
    
    ## Smoking
    j <- 2
    for (i in 2:length(tstart.comb)){
      if (tstart.comb[i] %in% pat.smoking$tstart){
        status_adj_smoking <- 
          append(status_adj_smoking, 
                 pat.smoking$med_status_adj[j])
        j <- j + 1
      } else {
        status_adj_smoking <- append(status_adj_smoking, status_adj_smoking[i-1])
      }
    }
    
  }

  ### Create data frame with split times
  out <- data.frame("tstart" = tstart.comb,
                    "cvd_time" = cvd_time.comb,
                    "med_status_adj_statins" = status_adj_statins,
                    "med_status_adj_ah" = status_adj_ah,
                    "med_status_adj_smoking" = status_adj_smoking)
  
  ### Return as a tibble
  return(tibble::tibble(out))
  
}

### Combine data.valid.split datasets, with an indicator for whether treatment is statins, antihypertensives, or smoking
### This step is purely so the input dataset will conform with dplyr:group_modify
cvd.split.times <- rbind(cvd.split.times.statins, cvd.split.times.ah, cvd.split.times.smoking) |>
  dplyr::arrange(patid, treatment, tstart, cvd_time) |>
  dplyr::group_by(patid)

### Run function
print(paste("start conversion", Sys.time()))
cohort.split.times <- dplyr::group_modify(.data = cvd.split.times, 
                                          .f = get_cut_surv_times_stat_ah_smoking) |>
  as.data.frame()
print(paste("end conversion", Sys.time()))

### Save adjsuted survival times, as it may take a while to run
saveRDS(cohort.split.times, paste("data/p4/cohort_split_times_smoking_statins_antihypertensives_burnout", burnout, "_id", group.id, ".rds", sep = ""))
