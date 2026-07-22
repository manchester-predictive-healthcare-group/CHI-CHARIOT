###
### Program to combine the interval censored event times created in p4.1_layer_cohort_split_times.R,
### into a single dataset.
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
burnout <- as.numeric(args[1])
num.groups <- as.numeric(args[2])
burnout <- 180
num.groups <- 100
print(paste("burnout = ", burnout))
print(paste("num.groups = ", num.groups))

### Read in data and combine
cohort.split.times <- lapply(1:num.groups, function(group.id) {
  readRDS(paste("data/p4/cohort_split_times_smoking_statins_antihypertensives_burnout", burnout, "_id", group.id, ".rds", sep = ""))}
)

cohort.split.times <- do.call("rbind", cohort.split.times)
print(paste("number of individuals = ", length(unique(cohort.split.times$patid))))
str(cohort.split.times)

### Read in cohort, as need to add cvd_indicator
cohort <- readRDS("data/p4/cohort_prototype3.rds") |>
  dplyr::select(patid, cvd_time, cvd_indicator) |>
  dplyr::arrange(patid)

### Merge the two by patid and cvd_time
### Note cvd_time in cohort will always match the last interval in cohort.split.times
cohort.split.times <- merge(cohort.split.times, cohort, by.x = c("patid", "cvd_time"), by.y = c("patid", "cvd_time"), all.x = TRUE)
testthat::expect_equal(sum(!is.na(cohort.split.times$cvd_indicator)), nrow(cohort))

### Observations with missing cvd_indicator, set to zero
cohort.split.times <- dplyr::mutate(cohort.split.times, 
                                    cvd_indicator = dplyr::case_when(is.na(cvd_indicator) ~ 0,
                                                                     TRUE ~ cvd_indicator)
                                    )

### Re-order the variables
cohort.split.times <- dplyr::relocate(cohort.split.times, patid, tstart, cvd_time, cvd_indicator)

### Test number of events matches that in cohort
testthat::expect_equal(sum(cohort.split.times$cvd_indicator), sum(cohort$cvd_indicator))
str(cohort.split.times)

###
### Save
saveRDS(cohort.split.times, paste("data/p4/cohort_split_times_smoking_statins_antihypertensives_burnout", burnout, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))

print(paste("number of individuals = ", length(unique(cohort.split.times$patid))))
