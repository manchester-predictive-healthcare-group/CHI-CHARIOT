###
### We need to estimate counterfactual survival times for the cohort when the index
### date is defined as 1/2/3/4/5 years post baseline/start of follow-up.
###
### This version uses data.table to adjust statin and antihypertensive status
### relative to treatment status at the new index date.
###

### Clear workspace
rm(list = ls())
print(Sys.time())

### Set working directory
setwd()
print(getwd())

### Define filepath to file directory system containing extracted data and functions
common.data.dir <- file.path("..", "..")

### Source functions
R.func.sources = list.files(file.path(common.data.dir, "Aurum_Jun2021_extract/R"), full.names = TRUE)
sapply(R.func.sources, source)

### Define burnout
burnout <- 180

### Read in split survival times
cohort_split_times <- readRDS(paste0("data/p5/cohort_split_times_smoking_statins_antihypertensives_burnout",burnout,".rds"))

### Convert to data.table and modify by reference to limit unnecessary copies
cohort_split_times <- data.table::as.data.table(cohort_split_times)

### Function to re-estimate split survival times based on new index date
calc_new_split_times <- function(cohort_split_times, t_fup_int){
  
  ### Define t_fup
  t_fup <- round(365.25*t_fup_int)
  
  ### Create new dataset
  cohort_split_times_new <- data.table::copy(cohort_split_times)
  
  ### Move the time origin to the new index date
  cohort_split_times_new[
    ,
    `:=`(
      tstart = tstart - t_fup,
      cvd_time = cvd_time - t_fup
    )
  ]
  
  ### Remove rows happening entirely before the new index date
  cohort_split_times_new <- cohort_split_times_new[cvd_time > 0]
  
  print(paste(
    "after removing individuals censored before fup time, number of individuals =",
    data.table::uniqueN(cohort_split_times_new$patid)
  ))
  
  ### Set tstart to zero where the interval began before the new index date
  cohort_split_times_new[tstart < 0, tstart := 0]
  
  ### Sort before identifying treatment status at the new index date
  ### data.table ordering is stable, so existing order is retained when tstart is tied
  ### after negative values have been set to zero.
  data.table::setorder(cohort_split_times_new, patid, tstart)
  
  ### Confirm that the first post-index status is valid for every patient.
  ### The original function assumes these values are limited to -1, 0 and 1 and
  ### would not handle a missing first value.
  invalid_baseline_status <- cohort_split_times_new[
    ,
    .(
      first_statin = med_status_adj_statins[1L],
      first_ah = med_status_adj_ah[1L]
    ),
    by = patid
  ][
    is.na(first_statin) |
      is.na(first_ah) |
      !first_statin %in% c(-1, 0, 1) |
      !first_ah %in% c(-1, 0, 1)
  ]
  
  if (nrow(invalid_baseline_status) > 0L) {
    stop(
      "At least one patient has a missing or invalid first post-index treatment status."
    )
  }
  
  ### Adjust each treatment-status history relative to its value at the new index date
  ### This is equivalent to the original rules:
  ###   first value  0: no change
  ###   first value  1: subtract 1
  ###   first value -1: add 1
  print(paste("Start data.table adjustment", Sys.time()))
  
  cohort_split_times_new[
    ,
    `:=`(
      med_status_adj_statins =
        med_status_adj_statins - med_status_adj_statins[1L],
      med_status_adj_ah =
        med_status_adj_ah - med_status_adj_ah[1L]
    ),
    by = patid
  ]
  
  ### Convert back to a standard data.frame as required
  cohort_split_times_new <- as.data.frame(cohort_split_times_new)
  
  ### Save output
  saveRDS(cohort_split_times_new, paste0("data/p5/temporalv_cohort_split_times_smoking_statins_antihypertensives_burnout", 
                                        burnout, "_t", 
                                        t_fup_int, ".rds"))
  print(paste("SAVED", t_fup, Sys.time()))
  
}

### Run the function
lapply(1:5, function(x) {calc_new_split_times(cohort_split_times, x)})
