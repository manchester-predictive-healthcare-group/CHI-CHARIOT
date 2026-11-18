###
### This program will extract smoking status during follow-up and create split-time survival data based on changes
### in smoking status during follow-up.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("./")
getwd()

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Source functions
R.func.sources = list.files(file.path(common.data.dir, "Aurum_Jun2021_extract/R"), full.names = TRUE)
sapply(R.func.sources, source)

###
### First extract all smoking records and define in same way as smoking was defined at baseline
###

### Read in cohort
cohort <- readRDS("data/p5/cohort_prototype3.rds") |>
  dplyr::select(patid, cvd_time, cvd_indicator, fup_start)

### Read in the db_qry of the events of interest (contains all smoking records of interest)
db_qry_smoking_non <- readRDS(paste("data/p5/db_qry_", "smoking_non", ".rds", sep = ""))
db_qry_smoking_ex <- readRDS(paste("data/p5/db_qry_", "smoking_ex", ".rds", sep = ""))
db_qry_smoking_light <- readRDS(paste("data/p5/db_qry_", "smoking_light", ".rds", sep = ""))
db_qry_smoking_mod <- readRDS(paste("data/p5/db_qry_", "smoking_mod", ".rds", sep = ""))
db_qry_smoking_heavy <- readRDS(paste("data/p5/db_qry_", "smoking_heavy", ".rds", sep = ""))

### Reduce to patients of interest
db_qry_smoking_non <- db_qry_smoking_non[!is.na(fastmatch::fmatch(db_qry_smoking_non$patid, cohort$patid)), ]
db_qry_smoking_ex <- db_qry_smoking_ex[!is.na(fastmatch::fmatch(db_qry_smoking_ex$patid, cohort$patid)), ]
db_qry_smoking_light <- db_qry_smoking_light[!is.na(fastmatch::fmatch(db_qry_smoking_light$patid, cohort$patid)), ]
db_qry_smoking_mod <- db_qry_smoking_mod[!is.na(fastmatch::fmatch(db_qry_smoking_mod$patid, cohort$patid)), ]
db_qry_smoking_heavy <- db_qry_smoking_heavy[!is.na(fastmatch::fmatch(db_qry_smoking_heavy$patid, cohort$patid)), ]

### Define smoking status as 0 = never, 1 = ex, 2 = current
db_qry_smoking_non$value <- 0
db_qry_smoking_ex$value <- 1

### Fucntion to use recorded cigarette count for current-smoking codelists, if recorded
recode_current_smoking <- function(x) {
  dplyr::case_when(
    is.na(x) ~ 2,
    x == 0 ~ 0,
    x > 0 & x <= 100 ~ 2,
    x > 100 ~ NA_real_
  )
}

### Apply function
db_qry_smoking_light$value <- recode_current_smoking(db_qry_smoking_light$value)
db_qry_smoking_mod$value <- recode_current_smoking(db_qry_smoking_mod$value)
db_qry_smoking_heavy$value <- recode_current_smoking(db_qry_smoking_heavy$value)

### Remove implausible values that were set to missing
db_qry_smoking_light <- db_qry_smoking_light[!is.na(db_qry_smoking_light$value),]
db_qry_smoking_mod <- db_qry_smoking_mod[!is.na(db_qry_smoking_mod$value),]
db_qry_smoking_heavy <- db_qry_smoking_heavy[!is.na(db_qry_smoking_heavy$value),]

### Combine, and take highest smoking status if multiple on same day
db_qry <- rbind(db_qry_smoking_non, db_qry_smoking_ex, db_qry_smoking_light, db_qry_smoking_mod, db_qry_smoking_heavy) |>
  dplyr::filter(!is.na(obsdate)) |>
  dplyr::arrange(patid, obsdate, dplyr::desc(value)) |>
  dplyr::group_by(patid, obsdate) |>
  dplyr::filter(dplyr::row_number() == 1) |>
  dplyr::ungroup() 


### Ensure records are ordered chronologically
db_qry <- data.table::as.data.table(db_qry)
data.table::setorder(db_qry, patid, obsdate)

### Identify whether Ex or Current has occurred previously
db_qry[
  ,
  previous_smoking := data.table::shift(
    cummax(value),
    fill = 0
  ),
  by = patid
]

### A later Never status becomes Ex if previous smoking exists
db_qry[
  value == 0 & previous_smoking %in% c(1, 2),
  value := 1
]

### Remove temporary variable
db_qry[, previous_smoking := NULL]

### Arrange and return to data frame
db_qry <- dplyr::arrange(db_qry, patid, obsdate) |>
  as.data.frame()

### Merge with cohort to get relative time to index date
db_qry <- dplyr::left_join(db_qry, cohort[,c("patid", "fup_start")], by = dplyr::join_by(patid)) |>
  dplyr::mutate(times_adj = obsdate - as.numeric(fup_start))

###
### Test to ensure a smoker/ex-smoker never goes back to being a non(never)-smoker

### Copy only the columns needed for the check
smoking_history_check <- data.table::as.data.table(
  db_qry
)[, .(patid, obsdate, value)]

### Order observations chronologically within each patient
data.table::setorder(smoking_history_check, patid, obsdate)

### Count never-smoker records occurring after previous smoking evidence
smoking_history_failures <- smoking_history_check[
  ,
  .(failures = sum(
    value == 0 &
      data.table::shift(
        cummax(as.integer(value %in% c(1, 2))),
        fill = 0L
      ) == 1L,
    na.rm = TRUE
  )),
  by = patid
][, sum(failures)]

### Remove the temporary data before checking the result
rm(smoking_history_check)

cat(
  "\nNever-smoker records after previous smoking evidence:",
  smoking_history_failures,
  "\n"
)

stopifnot(smoking_history_failures == 0L)

### Test over
### 

### Remove observations that happened prior to index date
db_qry_after <- dplyr::filter(db_qry, times_adj > 0)

########################################################################
### Split the survival times based on medication status change times ###
########################################################################

### Create split-survival dataset that contains one or more observation intervals.
### All patients receive an initial starting intervention, additional intervals begin when
### smoking status changes during follow-up.

### The first post-baseline smoking observation always starts a new interval, irrespective
### of whether it is different from baseline (we do not look at baseline value). This is because
### there is missing data on smoking at baseline, and therefore we cannot know if it has changed
### from baseline or not, until combining with the imputed data. Subsequent intervals only start
### when there is a change from the first post-baseline record.

### Convert inputs to data.table without modifying the original objects
cohort_for_split <- data.table::as.data.table(data.table::copy(cohort))
smoking_changes <- data.table::as.data.table(data.table::copy(db_qry_after))

### Retain only the columns required for interval construction
cohort_for_split <- cohort_for_split[,.(patid, cvd_time, cvd_indicator, fup_start)]
smoking_changes <- smoking_changes[,.(patid, times_adj, value)]

### Check the input structure
stopifnot(
  !anyNA(cohort_for_split$patid),
  !anyDuplicated(cohort_for_split$patid),
  all(is.finite(cohort_for_split$cvd_time)),
  all(cohort_for_split$cvd_time > 0),
  all(cohort_for_split$cvd_indicator %in% c(0, 1)),
  all(
    is.na(smoking_changes$value) |
      smoking_changes$value %in% c(0, 1, 2)
  )
)

### Ensure smoking observations are chronological
data.table::setorder(smoking_changes,patid,times_adj)

### Retain the first follow-up smoking observation and subsequent
### observations where smoking status changes
smoking_changes[,keep_change := seq_len(.N) == 1L | value != data.table::shift(value), by = patid]
smoking_changes <- smoking_changes[keep_change == TRUE]
smoking_changes[,keep_change := NULL]

### Add each patient's end of follow-up
smoking_changes <- cohort_for_split[smoking_changes, on = "patid", nomatch = 0L]

### Retain valid cut points occurring strictly within follow-up
smoking_changes <- smoking_changes[times_adj > 0 & times_adj < cvd_time]

### Create the initial interval start for every patient
initial_intervals <- cohort_for_split[,
  .(
    patid,
    tstart = 0,
    cvd_time,
    cvd_indicator,
    smoking = NA_real_,
    fup_start
  )
]

### Create an interval start at every smoking-status change
smoking_intervals <- smoking_changes[,
  .(
    patid,
    tstart = times_adj,
    cvd_time,
    cvd_indicator,
    smoking = value,
    fup_start
  )
]

### Combine initial intervals and smoking-change intervals
cohort_split_times <- data.table::rbindlist(
  list(
    initial_intervals,
    smoking_intervals
  ),
  use.names = TRUE,
  fill = FALSE
)

### Order interval starts within each patient
data.table::setorder(cohort_split_times, patid, tstart)

### Use the next interval start as the current interval end
cohort_split_times[, interval_end := data.table::shift(tstart,type = "lead"), by = patid]

### The final interval ends at the patient's observed end of follow-up
cohort_split_times[is.na(interval_end),interval_end := cvd_time]

### Only the final interval can contain the CVD event
cohort_split_times[interval_end < cvd_time, cvd_indicator := 0]

### Replace the original follow-up end with the interval-specific end
cohort_split_times[,cvd_time := interval_end]
cohort_split_times[,interval_end := NULL]

### Confirm that all resulting intervals are valid
stopifnot(
  !anyNA(cohort_split_times$patid),
  all(is.finite(cohort_split_times$tstart)),
  all(is.finite(cohort_split_times$cvd_time)),
  all(cohort_split_times$cvd_time > cohort_split_times$tstart),
  !anyDuplicated(
    cohort_split_times,
    by = c("patid", "tstart")
  ),
  data.table::uniqueN(cohort_split_times$patid) ==
    nrow(cohort_for_split)
)

### Return keep neccesary columns
cohort_split_times <- cohort_split_times[,.(patid, tstart, cvd_time, cvd_indicator, smoking, fup_start)]

### Convert to data frame
cohort_split_times <- as.data.frame(cohort_split_times) |> 
  dplyr::arrange(patid, tstart)
  

### Save data
saveRDS(cohort_split_times, "data/p5/cohort_split_times_smoking.rds")
print(paste("FINISHED", Sys.time()))