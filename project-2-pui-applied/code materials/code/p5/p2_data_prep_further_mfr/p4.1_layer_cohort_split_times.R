###
### Program to layer the split survival times from smoking, statins and antihypertensives,
### i.e. create split survival data with intervals whenever an individual changes status with
### respect to smoking, statins or antihypertensives.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("PATH/TO/ROOT")
getwd()

### Define burnout period
burnout <- as.numeric(180)
print(paste("burnout =", burnout))

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Read in the outcome data
cvd.split.times.smoking <- readRDS(paste("data/p5/cohort_split_times_smoking.rds", sep = ""))
cvd.split.times.statins <- readRDS(paste("data/p5/cohort_split_times_statins_burnout", burnout, ".rds", sep = ""))
cvd.split.times.ah <- readRDS(paste("data/p5/cohort_split_times_antihypertensives_burnout", burnout, ".rds", sep = ""))

### Read in cohort
cohort <- readRDS("data/p5/cohort_prototype3.rds") |>
  dplyr::select(patid, cvd_time, cvd_indicator, fup_start) |>
  dplyr::arrange(patid)

###
### Prepare the three source datasets
###

### Convert inputs to data.table copies so the original objects are not modified
cohort_dt <- data.table::as.data.table(data.table::copy(cohort))
smoking_dt <- data.table::as.data.table(data.table::copy(cvd.split.times.smoking))
statins_dt <- data.table::as.data.table(data.table::copy(cvd.split.times.statins))
ah_dt <- data.table::as.data.table(data.table::copy(cvd.split.times.ah))

### Retain the interval starts and status from each source
smoking_dt <- smoking_dt[, .(patid, tstart, med_status_adj_smoking = smoking)]
statins_dt <- statins_dt[, .(patid, tstart, med_status_adj_statins = med_status_adj)]
ah_dt <- ah_dt[, .(patid, tstart, med_status_adj_ah = med_status_adj)]

### Check that the cohort contains one valid row per patient
stopifnot(
  !anyNA(cohort_dt$patid),
  !anyDuplicated(cohort_dt$patid),
  all(is.finite(cohort_dt$cvd_time)),
  all(cohort_dt$cvd_time > 0),
  all(cohort_dt$cvd_indicator %in% c(0, 1))
)

### Check that each source has unique patient and interval-start combinations
stopifnot(
  !anyDuplicated(smoking_dt, by = c("patid", "tstart")),
  !anyDuplicated(statins_dt, by = c("patid", "tstart")),
  !anyDuplicated(ah_dt, by = c("patid", "tstart")),
  all(is.finite(smoking_dt$tstart)),
  all(is.finite(statins_dt$tstart)),
  all(is.finite(ah_dt$tstart)),
  all(smoking_dt$tstart >= 0),
  all(statins_dt$tstart >= 0),
  all(ah_dt$tstart >= 0)
)

### Check that all three sources contain the same cohort patients
stopifnot(
  setequal(unique(smoking_dt$patid), cohort_dt$patid),
  setequal(unique(statins_dt$patid), cohort_dt$patid),
  setequal(unique(ah_dt$patid), cohort_dt$patid)
)

### Check the expected status coding
stopifnot(
  all(is.na(smoking_dt$med_status_adj_smoking) | smoking_dt$med_status_adj_smoking %in% c(0, 1, 2)),
  all(is.na(statins_dt$med_status_adj_statins) | statins_dt$med_status_adj_statins %in% c(-1, 0, 1)),
  all(is.na(ah_dt$med_status_adj_ah) | ah_dt$med_status_adj_ah %in% c(-1, 0, 1))
)

###
### Create the combined interval-start grid
###

### Take the union of all interval starts across the three time-varying variables
cohort.split.times <- unique(
  data.table::rbindlist(
    list(
      smoking_dt[, .(patid, tstart)],
      statins_dt[, .(patid, tstart)],
      ah_dt[, .(patid, tstart)]
    ),
    use.names = TRUE
  )
)

### Add the status recorded at each exact interval start
cohort.split.times[smoking_dt, med_status_adj_smoking := i.med_status_adj_smoking, on = .(patid, tstart)]
cohort.split.times[statins_dt, med_status_adj_statins := i.med_status_adj_statins, on = .(patid, tstart)]
cohort.split.times[ah_dt, med_status_adj_ah := i.med_status_adj_ah, on = .(patid, tstart)]

### Order interval starts chronologically within patients
data.table::setorder(cohort.split.times, patid, tstart)

### Carry each status forward until that status changes again
### A leading missing smoking status remains missing until the first observed
### follow-up smoking record; downstream code then uses baseline smoking status
cohort.split.times[, med_status_adj_smoking := data.table::nafill(med_status_adj_smoking, type = "locf"), by = patid]
cohort.split.times[, med_status_adj_statins := data.table::nafill(med_status_adj_statins, type = "locf"), by = patid]
cohort.split.times[, med_status_adj_ah := data.table::nafill(med_status_adj_ah, type = "locf"), by = patid]

###
### Define the combined survival intervals
###

### Add each patient's original outcome or censoring time and event indicator
cohort.split.times <- cohort_dt[
  cohort.split.times,
  on = "patid",
  nomatch = 0L
]

### Retain interval starts strictly before the end of follow-up
cohort.split.times <- cohort.split.times[tstart < cvd_time]

### Confirm that every patient has an initial interval at time zero
stopifnot(
  cohort.split.times[, sum(tstart == 0), by = patid][, all(V1 == 1L)]
)

### Use the next interval start as the end of the current interval
data.table::setorder(cohort.split.times, patid, tstart)
cohort.split.times[, interval_end := data.table::shift(tstart, type = "lead"), by = patid]

### The final interval ends at the patient's original event or censoring time
cohort.split.times[is.na(interval_end), interval_end := cvd_time]

### Only the final interval can contain the CVD event
cohort.split.times[interval_end < cvd_time, cvd_indicator := 0L]

### Replace the patient-level follow-up time with the interval-specific end time
cohort.split.times[, cvd_time := interval_end]
cohort.split.times[, interval_end := NULL]

###
### Validate the final layered dataset
###

### Check interval validity, patient coverage and unique interval starts
stopifnot(
  !anyNA(cohort.split.times$patid),
  all(is.finite(cohort.split.times$tstart)),
  all(is.finite(cohort.split.times$cvd_time)),
  all(cohort.split.times$cvd_time > cohort.split.times$tstart),
  !anyDuplicated(cohort.split.times, by = c("patid", "tstart")),
  data.table::uniqueN(cohort.split.times$patid) == nrow(cohort_dt),
  sum(cohort.split.times$cvd_indicator) == sum(cohort_dt$cvd_indicator)
)

### Check that consecutive intervals are contiguous with no gaps or overlaps
cohort.split.times[, previous_end := data.table::shift(cvd_time), by = patid]
stopifnot(cohort.split.times[tstart > 0, all(tstart == previous_end)])
cohort.split.times[, previous_end := NULL]

### Return the column structure expected by model development
cohort.split.times <- cohort.split.times[,.(
  patid,
  tstart,
  cvd_time,
  cvd_indicator,
  med_status_adj_statins,
  med_status_adj_ah,
  med_status_adj_smoking)]

cohort.split.times <- as.data.frame(cohort.split.times) |>
  dplyr::arrange(patid, tstart)

###
### Save the complete layered split-survival dataset
###
saveRDS(cohort.split.times,paste0("data/p5/cohort_split_times_smoking_statins_antihypertensives_burnout",burnout,".rds"))

print(paste("FINISHED", Sys.time()))
print(paste("number of individuals =", length(unique(cohort.split.times$patid))))