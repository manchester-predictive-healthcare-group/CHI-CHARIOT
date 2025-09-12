###
### Program to estimate counterfactual treatment times for individuals on the current treatment strategy
### Done in groups for parallelisation
###

###
### We adjust survival times for changes in treatment (statins/ah/smoking) use during follow-up
### We do not adjust for different BMI/nonhdl/sbp/smoking at baseline, as these are just what will be made for the prediction
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

### Define list elements to cycle through
crossing <- tidyr::expand_grid(
  "devel.num" = 1:10,
  "valid.num" = 1:10,
  "gender" = 1:2)

### Extract relevant scenario from command line
row <- as.numeric(commandArgs(trailingOnly = T)[1])
# row <- 1
print(crossing[row,])

### Define variables that change from scenario to scenario

### Extract imputation number for development dataset
gender <- as.numeric(crossing[row, "gender"])
gender_char <- c("male", "female")[gender]
print(paste("gender = ", gender_char))

### Extract imputation number for development dataset
devel.num <- as.numeric(crossing[row, "devel.num"])
print(paste("devel.num = ", devel.num))

### Extract imputation number for validation dataset
valid.num <- as.numeric(crossing[row, "valid.num"])
print(paste("valid.num = ", valid.num))

### Define burnout
burnout <- 180
print(paste("burnout = ", burnout))

### Read in the outcome data, which contains the split survival times
cohort.split.times <- readRDS(paste("data/p4/cohort_split_times_smoking_statins_antihypertensives_burnout", burnout, ".rds", sep = ""))

### Read in imputed validation datasets
imp.list <- readRDS(paste("data/p4/dfs_valid_", gender_char, ".rds", sep = ""))

### Pick one imputed dataset, as the adjusted times will be the same for every validation dataset
### NB: Counterfactual survival tims will be different dependent on the development dataset, at this alters the baseline hazard
data.valid <- imp.list[[valid.num]]
nrows_test <- nrow(data.valid)

### Merge imputed dataset with cohort.split.times
data.valid <- dplyr::left_join(dplyr::select(data.valid, c("patid", "smoking")),
                               cohort.split.times,
                               by = dplyr::join_by("patid"))

testthat::expect_equal(length(unique(data.valid$patid)), nrows_test)
# str(data.valid)
# str(cohort.split.times)

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

### Create appropriate offsets based on treatment effects
data.valid$offset_statins_timevar_lnHR <- lnHR_statins*data.valid$med_status_adj_statins
data.valid$offset_ah_timevar_lnHR <- lnHR_ah*data.valid$med_status_adj_ah
data.valid$offset_smoking_timevar_dummy1_lnHR <- lnHR_smoking_dummy1_total*data.valid$med_status_adj_smoking_dummy1
data.valid$offset_smoking_timevar_dummy2_lnHR <- lnHR_smoking_dummy2_total*data.valid$med_status_adj_smoking_dummy2

### Read in the baseline hazard
bhaz <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", devel.num, ".rds", sep = ""))

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
time.in <- Sys.time()
cf.surv.times <- dplyr::group_modify(.data = data.valid, .f = get_survtimes_adj)
time.out <- Sys.time()
time.out - time.in

warnings()

### Save adjusted survival times, as it may take a while to run
saveRDS(cf.surv.times, paste("data/p4/prototype3_cf_surv_times_", gender, "_devel", devel.num, "_valid", valid.num, ".rds", sep = ""))
