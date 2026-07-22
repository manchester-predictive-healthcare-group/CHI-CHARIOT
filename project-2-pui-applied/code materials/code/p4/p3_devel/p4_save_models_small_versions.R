### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

###
### Read in female and male models
###

### Read in prototye3 female
fit_prototype3_female <- lapply(1:10, function(x) {readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", 2, "_imp", x, ".rds", sep = ""))})
bhaz_prototype3_female <- lapply(1:10, function(x) {readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", 2, "_imp", x, ".rds", sep = ""))})
# bhaz_uncent_prototype3_female <- lapply(1:10, function(x) {readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_uncent_", 2, "_imp", x, ".rds", sep = ""))})

### Read in prototype3
fit_prototype3_male <- lapply(1:10, function(x) {readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", 1, "_imp", x, ".rds", sep = ""))})
bhaz_prototype3_male <- lapply(1:10, function(x) {readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", 1, "_imp", x, ".rds", sep = ""))})
# bhaz_uncent_prototype3_male <- lapply(1:10, function(x) {readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_uncent_", 1, "_imp", x, ".rds", sep = ""))})

### Extract development means for offset models
means_prototype3_female <- vector("list", 4)
means_prototype3_female[[1]] <- lapply(fit_prototype3_female, function(x) {mean(x$model[,"offset(offset_statins_timevar_lnHR)"])})
means_prototype3_female[[2]] <- lapply(fit_prototype3_female, function(x) {mean(x$model[,"offset(offset_ah_timevar_lnHR)"])})
means_prototype3_female[[3]] <- lapply(fit_prototype3_female, function(x) {mean(x$model[,"offset(offset_smoking_timevar_dummy1_lnHR)"])})
means_prototype3_female[[4]] <- lapply(fit_prototype3_female, function(x) {mean(x$model[,"offset(offset_smoking_timevar_dummy2_lnHR)"])})
names(means_prototype3_female) <- c("statins", "ah", "smoking1", "smoking2")

means_prototype3_male <- vector("list", 4)
means_prototype3_male[[1]] <- lapply(fit_prototype3_male, function(x) {mean(x$model[,"offset(offset_statins_timevar_lnHR)"])})
means_prototype3_male[[2]] <- lapply(fit_prototype3_male, function(x) {mean(x$model[,"offset(offset_ah_timevar_lnHR)"])})
means_prototype3_male[[3]] <- lapply(fit_prototype3_male, function(x) {mean(x$model[,"offset(offset_smoking_timevar_dummy1_lnHR)"])})
means_prototype3_male[[4]] <- lapply(fit_prototype3_male, function(x) {mean(x$model[,"offset(offset_smoking_timevar_dummy2_lnHR)"])})
names(means_prototype3_male) <- c("statins", "ah", "smoking1", "smoking2")

### Get rid of patient level data from model objects
fit_prototype3_male <- lapply(fit_prototype3_male, function(x) {
  x[names(x) == "model"] <- "DELETED"
  x[names(x) == "y"] <- "DELETED"
  x[names(x) == "linear.predictors"] <- "DELETED"
  x[names(x) == "residuals"] <- "DELETED"
  return(x)})

fit_prototype3_female <- lapply(fit_prototype3_female, function(x) {
  x[names(x) == "model"] <- "DELETED"
  x[names(x) == "y"] <- "DELETED"
  x[names(x) == "linear.predictors"] <- "DELETED"
  x[names(x) == "residuals"] <- "DELETED"
  return(x)})

### Reduce size of bhaz object to just baseline hazard at 1 year intervals
## Female
bhaz_prototype3_female <- lapply(bhaz_prototype3_female, function(x) {
  lapply(1:10, function(y){
    return(x$hazard[max(which(x$time <= y*365.25))])
  })
})

# bhaz_uncent_prototype3_female <- lapply(bhaz_uncent_prototype3_female, function(x) {
#   lapply(1:10, function(y){
#     return(x$hazard[max(which(x$time <= y*365.25))])
#   })
# })

## Male
bhaz_prototype3_male <- lapply(bhaz_prototype3_male, function(x) {
  lapply(1:10, function(y){
    return(x$hazard[max(which(x$time <= y*365.25))])
  })
})

# bhaz_uncent_prototype3_male <- lapply(bhaz_uncent_prototype3_male, function(x) {
#   lapply(1:10, function(y){
#     return(x$hazard[max(which(x$time <= y*365.25))])
#   })
# })

### Save these
saveRDS(fit_prototype3_female, "data/p4/small_fit_prototype3_female.rds")
saveRDS(means_prototype3_female, "data/p4/small_means_prototype3_female.rds")
saveRDS(bhaz_prototype3_female, "data/p4/small_bhaz_prototype3_female.rds")
# saveRDS(bhaz_uncent_prototype3_female, "data/p4/small_bhaz_uncent_prototype3_female.rds")

saveRDS(fit_prototype3_male, "data/p4/small_fit_prototype3_male.rds")
saveRDS(means_prototype3_male, "data/p4/small_means_prototype3_male.rds")
saveRDS(bhaz_prototype3_male, "data/p4/small_bhaz_prototype3_male.rds")
# saveRDS(bhaz_uncent_prototype3_male, "data/p4/small_bhaz_uncent_prototype3_male.rds")

print(paste("FINISHED", Sys.time()))

### Read in validation datasets
imps_valid <- readRDS(paste("data/p4/dfs_valid_", "female", ".rds", sep = ""))

### Pick one
data_valid <- imps_valid[[1]]

### Obtain a patient
pat <- data_valid[1, ]

### Assign all variables and de-identfiy...
pat <- dplyr::select(pat, -c(cvd_time, cvd_indicator, cvd_time_prim, cvd_indicator_prim,
                             cvd_time_hes, cvd_indicator_hes, cvd_time_death, cvd_indicator_death,
                             patid, pracid, gender, cvd_ev_prim_aj, statins, antihypertensives))
pat$age <- 70
pat$bmi <- 30
pat$sbp <- 160
pat$nonhdl <- 10
pat$smoking[1] <- "Current"
pat$ethnicity[1] <- "black caribbean"
pat$IMD <- 10
data_valid$caltime <- 5537

saveRDS(pat, "data/p4/worked_example_patient.rds")


# 
# ### The prototype3 dataset didn't have sbp variability
# ### I am therefore going to put in a mean value for when estimating the risks of the standard model, meaning
# ### the two aren't directly comparable
# imps.valid.standard <- readRDS(paste("data/dfs_valid_", "female", ".rds", sep = ""))
# pat$sbp_var <- mean(imps.valid.standard[[1]]$sbp_var)
# pat$cholhdl_ratio <- mean(imps.valid.standard[[1]]$cholhdl_ratio)
# 
# saveRDS(pat, "code/p4/p7_rshiny/data/fake.pat.rds")
# 
# 
# ### Save the HRs/ORs
# 
# ### Total effects
# total_hr_ah <- readRDS("data/p4/offsets_lnHR_ah.rds")
# total_hr_statins <- readRDS("data/p4/offsets_lnHR_statins.rds")
# total_hr_smoking1 <- readRDS("data/p4/offsets_lnHR_smoking_dummy1_total.rds")
# total_hr_smoking2 <- readRDS("data/p4/offsets_lnHR_smoking_dummy2_total.rds")
# 
# saveRDS(total_hr_ah, "code/p4/p7_rshiny/data/offsets_lnHR_ah.rds")
# saveRDS(total_hr_statins, "code/p4/p7_rshiny/data/offsets_lnHR_statins.rds")
# saveRDS(total_hr_smoking1, "code/p4/p7_rshiny/data/offsets_lnHR_smoking_dummy1_total.rds")
# saveRDS(total_hr_smoking2, "code/p4/p7_rshiny/data/offsets_lnHR_smoking_dummy2_total.rds")
# 
# ### Direct effects
# direct_RR_smoking_initiation <- readRDS("data/p4/direct_RR_smoking_initiation")
# direct_RR_smoking_cessation <- readRDS("data/p4/direct_RR_smoking_cessation")
# direct_RR_nonhdl <- readRDS("data/p4/direct_RR_nonhdl")
# direct_RR_bmi <- readRDS("data/p4/direct_RR_bmi")
# direct_RR_sbp <- readRDS("data/p4/direct_RR_sbp")
# 
# saveRDS(direct_RR_smoking_initiation, "code/p4/p7_rshiny/data/direct_RR_smoking_initiation.rds")
# saveRDS(direct_RR_smoking_cessation, "code/p4/p7_rshiny/data/direct_RR_smoking_cessation.rds")
# saveRDS(direct_RR_nonhdl, "code/p4/p7_rshiny/data/direct_RR_nonhdl.ds")
# saveRDS(direct_RR_bmi, "code/p4/p7_rshiny/data/direct_RR_bmi.rds")
# saveRDS(direct_RR_sbp, "code/p4/p7_rshiny/data/direct_RR_sbp.rds")
# print("FINISHED")
