###
### Program to estimate counterfactual treatment times.
### Note they are calculated differently for each model
### For model 1...
### For model 2...
### For model 3...
### For model 4...
###

###
### We adjust survival times for changes in treatment (statins/ah/smoking) use during follow-up
### We do not adjust for different BMI/nonhdl/sbp/smoking at baseline, as these are just what will be made for the prediction
### 

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Extract imputation number for development dataset
args <- commandArgs(trailingOnly = T)
gender <- as.numeric(args[1])
model <- as.numeric(args[2])
gender_char <- c("male", "female")[gender]
print(paste("gender = ", gender_char))
print(paste("model = ", model))

### Define burnout
burnout <- 180
print(paste("burnout = ", burnout))

### Define HR for offsets
lnHR_statins <- readRDS("data/offsets_total_lnHR_statins.rds")
lnHR_ah <- readRDS("data/offsets_total_lnHR_ah.rds")
lnHR_sbp <- readRDS("data/offsets_direct_lnHR_sbp.rds")
lnHR_bmi <- readRDS("data/offsets_direct_lnHR_bmi.rds")
lnHR_nonhdl <- readRDS("data/offsets_direct_lnHR_nonhdl.rds")

### Read in the interval censored outcome times NB CHANGE WHEN FINISHED
### Note the time varying medication status variables are defined differently for models 1 and 2 compared to models 3 and 4.
if (model %in% c(0,1,2)){
  cohort_split_times <- readRDS(paste("data/cohort_split_times_burnout", burnout, ".rds", sep = ""))
} else if (model %in% c(3,4,5,6,7)){
  cohort_split_times <- readRDS(paste("data/cohort_split_times_augmented_burnout", burnout, ".rds", sep = ""))
}

### Read in validation data
df_valid <- readRDS(paste("data/df_imp_valid_", gender, sep = ""))

### Merge imputed dataset with cohort_split_times
df_valid <- dplyr::left_join(dplyr::select(df_valid, c("patid")),
                               cohort_split_times,
                               by = dplyr::join_by("patid"))

### Create appropriate offsets based on treatment effects
df_valid$offset_statins_timevar_lnHR <- lnHR_statins*df_valid$med_status_statins
df_valid$offset_ah_timevar_lnHR <- lnHR_ah*df_valid$med_status_ah

### Read in the baseline hazard
bhaz <- readRDS(paste("data/bhaz_", gender, "_model", model, ".rds", sep = ""))

### Group data
df_valid <- dplyr::arrange(df_valid, patid, tstart) |>
  dplyr::group_by(patid)

## Run function and arrange
time.in <- Sys.time()
if (model %in% c(1,2)){
  cf_surv_times <- dplyr::group_modify(.data = df_valid, .f = get_survtimes_adj_model12)
} else if (model %in% c(3,4,5,6,7)){
  cf_surv_times <- dplyr::group_modify(.data = df_valid, .f = get_survtimes_adj_model34567)
}
time.out <- Sys.time()
time.out - time.in

warnings()

### Save adjusted survival times, as it may take a while to run
saveRDS(cf_surv_times, paste("data/cf_surv_times_", gender, "_model", model, ".rds", sep = ""))
