###
### This program calculates the counterfactual survival times for the cohort when the index
### date is defined as 1/2/3/4/5 post baseline/start of follow-up.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- 180

### Define model = 4 (two component model) for these analyses
model <- 4

### Extract gender from command line
args <- commandArgs(trailingOnly = T)
gender <- as.numeric(args[1])
gender_char <- c("male", "female")[gender]
print(paste("gender = ", gender_char))

### Extract followup time from command line
t_fup <- round(365.25*as.numeric(args[2]))
print(paste("t_fup = ", t_fup))

### Read in the outcome data, which contains the split survival times
cohort_split_times <- readRDS(paste("data/sens_temporalv_cohort_split_times_augmented_", gender, 
                                    "_burnout", burnout, 
                                    "_t", t_fup, ".rds", sep = ""))

### Read in the dataset with index date at t_fup years post start of follow-up.
df_valid <- readRDS(paste("data/sens_temporalv_df_valid_", gender, "_", t_fup, ".rds", sep = ""))

### Test
testthat::expect_equal(length(unique(df_valid$patid)), length(unique(cohort_split_times$patid)))

# ### Merge df_valid with cohort_split_times (I think this can be removed? cohort_split_times already reduced to relevant individuals, 
# ### and previously wanted to include smoking status at baseline)
# df_valid <- dplyr::left_join(cohort_split_times,
#                              dplyr::select(df_valid, c("patid")),
#                              by = dplyr::join_by("patid"))

### Define HR for offsets
lnHR_statins <- readRDS("data/offsets_total_lnHR_statins.rds")
lnHR_ah <- readRDS("data/offsets_total_lnHR_ah.rds")

### Create appropriate offsets based on treatment effects
cohort_split_times$offset_statins_timevar_lnHR <- lnHR_statins*cohort_split_times$med_status_statins
cohort_split_times$offset_ah_timevar_lnHR <- lnHR_ah*cohort_split_times$med_status_ah

### Read in the baseline hazard
bhaz <- readRDS(paste("data/bhaz_", gender, "_model", model, ".rds", sep = ""))

### Write a function to extract bhaz for a specific time
get.bhaz <- function(time, bhaz){
  if (time == 0){
    out <- 0
  } else {
    out <- as.numeric(bhaz$hazard[max(which(bhaz$time <= time))])
  } 
  return(out)
}

### Function to get counterfactual adjsuted survival times
### Function written so it can be applied in conjunction with group_modify
### B is the log(hazard) of the effect of statins or ah
get_survtimes_adj <- function(data, id){
  
  ### Calculate cumulative hazard at transition times
  data$cumhaz_tstart <- unlist(lapply(data$tstart, get.bhaz, bhaz = bhaz))
  data$cumhaz_tstop <- unlist(lapply(data$cvd_time, get.bhaz, bhaz = bhaz))
  
  ### Intervals which are "on treatment" get multiplied by exp(B)
  data <- dplyr::mutate(data,
                        cumhaz_tstart = cumhaz_tstart*
                          exp(offset_statins_timevar_lnHR)*
                          exp(offset_ah_timevar_lnHR),
                        cumhaz_tstop = cumhaz_tstop*
                          exp(offset_statins_timevar_lnHR)*
                          exp(offset_ah_timevar_lnHR))
  
  ### Then calculate total hazard
  cumhaz_adj <- sum(data$cumhaz_tstop - data$cumhaz_tstart)
  
  ### Now see what value of t corresponds to a hazard of this
  ### Want minimum time, for which this amount of hazard is reached
  t.out <- as.numeric(bhaz$time[min(which(bhaz$hazard >= cumhaz_adj))])
  
  ### Get cvd_indicator
  cvd_indicator <- max(data$cvd_indicator)
  
  ### GET NA values for when we increase hazard over the max (this happens to people with long survival times,
  ### and increased risk factors. We can set these to the max followup time)
  if (is.na(t.out)){
    t.out <- max(bhaz$time)
    cvd_indicator <- 0}
  
  ### Create output tibble
  out <- data.frame("cvd_time_cf" = t.out, "cvd_indicator_cf" = cvd_indicator)
  
  return(tibble::tibble(out))
  
}

### Group data
cohort_split_times <- dplyr::arrange(cohort_split_times, patid, tstart) |>
  dplyr::group_by(patid)

## Run function and arrange
print("START GROUP MODIFY")
time.in <- Sys.time()
cf_surv_times <- dplyr::group_modify(.data = cohort_split_times, .f = get_survtimes_adj)
time.out <- Sys.time()
time.out - time.in

warnings()

### Save adjusted survival times, as it may take a while to run
saveRDS(cf_surv_times, paste("data/sens_temporalv__cf_surv_times_", gender, 
                             "_t", t_fup, 
                             ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))