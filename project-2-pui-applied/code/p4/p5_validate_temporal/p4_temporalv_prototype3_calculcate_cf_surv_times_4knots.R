###
### This program calculates the counterfactual survival times for the cohort when the index
### date is defined as 1/2/3/4/5 post baseline/start of follow-up.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- 180

### Extract gender from command line
args <- commandArgs(trailingOnly = T)
gender <- as.numeric(args[1])
gender_char <- c("male", "female")[gender]
print(paste("gender = ", gender_char))

### Extract followup time from command line
t_fup <- round(365.25*as.numeric(args[2]))
print(paste("t_fup = ", t_fup))

### Define validation and development datasets (we are only doing temporal validation for 1)
valid.num <- 1
devel.num <- 1

### Read in the outcome data, which contains the split survival times
cohort.split.times <- readRDS(paste("data/p4/validation_cohort_split_times_smoking_statins_antihypertensives_burnout", 
                                    burnout, "_t", 
                                    t_fup, ".rds", sep = ""))

### Read in the dataset with index date at t_fup years post start of follow-up.
### NB: This contains a smoking variable with is the most recently recorded smoking status value.
### If this is missing, it takes the imputed value from the first imputed cohort with index date at start of follow-up.
### (see p1.2_temporalv_create_individual_cohorts.R)
data.valid <- readRDS(paste("data/p4/df_valid_temporalv_", gender, "_", t_fup, ".rds", sep = ""))

### Merge imputed dataset with cohort.split.times
data.valid <- dplyr::left_join(cohort.split.times,
                               dplyr::select(data.valid, c("patid", "smoking")),
                               by = dplyr::join_by("patid"))

### Note, cohort.split.times only contains individuals that are not censored prior
### to t_fup. The way this merge is done, will remove the individuals from data.valid
### that we no longer need.

###
### Create dummy variables for smoking (at baseline, and during follow-up)
### 
### There are seperate dummies at baseline and the time-varying variable during follow-up
### See R/functions.R for a deeper explanation of these
###
data.valid <- create_smoking_dummies(data.valid)

### Define HR for offsets
lnHR_statins <- readRDS("data/p4/offsets_lnHR_statins.rds")
lnHR_ah <- readRDS("data/p4/offsets_lnHR_ah.rds")
lnHR_smoking_dummy1_total <- readRDS("data/p4/offsets_lnHR_smoking_dummy1_total.rds")
lnHR_smoking_dummy2_total <- readRDS("data/p4/offsets_lnHR_smoking_dummy2_total.rds")
lnHR_smoking_dummy1_direct <- readRDS("data/p4/offsets_lnHR_smoking_dummy1_direct.rds")
lnHR_smoking_dummy2_direct <- readRDS("data/p4/offsets_lnHR_smoking_dummy2_direct.rds")
lnHR_sbp <- readRDS("data/p4/offsets_lnHR_sbp.rds")
lnHR_bmi <- readRDS("data/p4/offsets_lnHR_bmi.rds")
lnHR_nonhdl <- readRDS("data/p4/offsets_lnHR_nonhdl.rds")

### Create appropriate offsets based on treatment effects
data.valid$offset_statins_timevar_lnHR <- lnHR_statins*data.valid$med_status_adj_statins
data.valid$offset_ah_timevar_lnHR <- lnHR_ah*data.valid$med_status_adj_ah
data.valid$offset_smoking_timevar_dummy1_lnHR <- lnHR_smoking_dummy1_total*data.valid$med_status_adj_smoking_dummy1
data.valid$offset_smoking_timevar_dummy2_lnHR <- lnHR_smoking_dummy2_total*data.valid$med_status_adj_smoking_dummy2

### Read in the baseline hazard
bhaz <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", devel.num, ".rds", sep = ""))

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
  data$cumhaz.tstart <- unlist(lapply(data$tstart, get.bhaz, bhaz = bhaz))
  data$cumhaz.tstop <- unlist(lapply(data$cvd_time, get.bhaz, bhaz = bhaz))
  
  ### Intervals which are "on treatment" get multiplied by exp(B)
  data <- dplyr::mutate(data,
                        cumhaz.tstart = cumhaz.tstart*
                          exp(offset_statins_timevar_lnHR)*
                          exp(offset_ah_timevar_lnHR)*
                          exp(offset_smoking_timevar_dummy1_lnHR)*
                          exp(offset_smoking_timevar_dummy2_lnHR),
                        cumhaz.tstop = cumhaz.tstop*
                          exp(offset_statins_timevar_lnHR)*
                          exp(offset_ah_timevar_lnHR)*
                          exp(offset_smoking_timevar_dummy1_lnHR)*
                          exp(offset_smoking_timevar_dummy2_lnHR))
  
  ### Then calculate total hazard
  cumhaz_adj <- sum(data$cumhaz.tstop - data$cumhaz.tstart)
  
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
data.valid <- dplyr::arrange(data.valid, patid, tstart) |>
  dplyr::group_by(patid)

## Run function and arrange
print("START GROUP MODIFY")
time.in <- Sys.time()
cf.surv.times <- dplyr::group_modify(.data = data.valid, .f = get_survtimes_adj)
time.out <- Sys.time()
time.out - time.in

warnings()

### Save adjusted survival times, as it may take a while to run
saveRDS(cf.surv.times, paste("data/p4/prototype3_4knots_cf_surv_times_", gender, 
                             "_devel", devel.num, 
                             "_valid", valid.num, 
                             "_t", t_fup, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))