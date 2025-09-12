### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

###
### Read in female and male models
###

### Read in prototye3
fit.prototype3.female <- lapply(1:1, function(x) {readRDS(paste("data/p4/prototype3_4knots_cox_", 2, "_imp", x, ".rds", sep = ""))})
bhaz.prototype3.female <- lapply(1:1, function(x) {readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", 2, "_imp", x, ".rds", sep = ""))})

### Standard model
fit.standard.female <- lapply(1:1, function(x) {readRDS(paste("data/standard_cox_", 2, "_imp", x, ".rds", sep = ""))})
bhaz.standard.female <- lapply(1:1, function(x) {readRDS(paste("data/standard_cox_bhaz_", 2, "_imp", x, ".rds", sep = ""))})

### Read in prototype3
fit.prototype3.male <- lapply(1:1, function(x) {readRDS(paste("data/p4/prototype3_4knots_cox_", 1, "_imp", x, ".rds", sep = ""))})
bhaz.prototype3.male <- lapply(1:1, function(x) {readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", 1, "_imp", x, ".rds", sep = ""))})

### Standard model
fit.standard.male <- lapply(1:1, function(x) {readRDS(paste("data/standard_cox_", 1, "_imp", x, ".rds", sep = ""))})
bhaz.standard.male <- lapply(1:1, function(x) {readRDS(paste("data/standard_cox_bhaz_", 1, "_imp", x, ".rds", sep = ""))})
  

### Extract development means for offset models
means.prototype3.female <- vector("list", 4)
means.prototype3.female[[1]] <- lapply(fit.prototype3.female, function(x) {mean(x$model[,"offset(offset_statins_timevar_lnHR)"])})
means.prototype3.female[[2]] <- lapply(fit.prototype3.female, function(x) {mean(x$model[,"offset(offset_ah_timevar_lnHR)"])})
means.prototype3.female[[3]] <- lapply(fit.prototype3.female, function(x) {mean(x$model[,"offset(offset_smoking_timevar_dummy1_lnHR)"])})
means.prototype3.female[[4]] <- lapply(fit.prototype3.female, function(x) {mean(x$model[,"offset(offset_smoking_timevar_dummy2_lnHR)"])})
names(means.prototype3.female) <- c("statins", "ah", "smoking1", "smoking2")

means.prototype3.male <- vector("list", 4)
means.prototype3.male[[1]] <- lapply(fit.prototype3.male, function(x) {mean(x$model[,"offset(offset_statins_timevar_lnHR)"])})
means.prototype3.male[[2]] <- lapply(fit.prototype3.male, function(x) {mean(x$model[,"offset(offset_ah_timevar_lnHR)"])})
means.prototype3.male[[3]] <- lapply(fit.prototype3.male, function(x) {mean(x$model[,"offset(offset_smoking_timevar_dummy1_lnHR)"])})
means.prototype3.male[[4]] <- lapply(fit.prototype3.male, function(x) {mean(x$model[,"offset(offset_smoking_timevar_dummy2_lnHR)"])})
names(means.prototype3.male) <- c("statins", "ah", "smoking1", "smoking2")

### Get rid of patient level data from model objects
fit.prototype3.male <- lapply(fit.prototype3.male, function(x) {
  x[names(x) == "model"] <- "DELETED"
  x[names(x) == "y"] <- "DELETED"
  x[names(x) == "linear.predictors"] <- "DELETED"
  x[names(x) == "residuals"] <- "DELETED"
  return(x)})

fit.standard.male <- lapply(fit.standard.male, function(x) {
  x[names(x) == "model"] <- "DELETED"
  x[names(x) == "y"] <- "DELETED"
  x[names(x) == "linear.predictors"] <- "DELETED"
  x[names(x) == "residuals"] <- "DELETED"
  return(x)})

fit.prototype3.female <- lapply(fit.prototype3.female, function(x) {
  x[names(x) == "model"] <- "DELETED"
  x[names(x) == "y"] <- "DELETED"
  x[names(x) == "linear.predictors"] <- "DELETED"
  x[names(x) == "residuals"] <- "DELETED"
  return(x)})

fit.standard.female <- lapply(fit.standard.female, function(x) {
  x[names(x) == "model"] <- "DELETED"
  x[names(x) == "y"] <- "DELETED"
  x[names(x) == "linear.predictors"] <- "DELETED"
  x[names(x) == "residuals"] <- "DELETED"
  return(x)})


### Reduce size of bhaz object
bhaz.prototype3.male <- lapply(bhaz.prototype3.male, function(x) {
  return(x$hazard[max(which(x$time <= 10*365.25))])}
  )

bhaz.standard.male <- lapply(bhaz.standard.male, function(x) {
  return(x$hazard[max(which(x$time <= 10*365.25))])})

bhaz.prototype3.female <- lapply(bhaz.prototype3.female, function(x) {
  return(x$hazard[max(which(x$time <= 10*365.25))])})

bhaz.standard.female <- lapply(bhaz.standard.female, function(x) {
  return(x$hazard[max(which(x$time <= 10*365.25))])})


### Save these to Rshiny folder
saveRDS(fit.prototype3.female, "code/p4/p7_rshiny/data/fit.prototype3.female.rds")
saveRDS(means.prototype3.female, "code/p4/p7_rshiny/data/means.prototype3.female.rds")
saveRDS(fit.standard.female, "code/p4/p7_rshiny/data/fit.standard.female.rds")
saveRDS(bhaz.prototype3.female, "code/p4/p7_rshiny/data/bhaz.prototype3.female.rds")
saveRDS(bhaz.standard.female, "code/p4/p7_rshiny/data/bhaz.standard.female.rds")

saveRDS(fit.prototype3.male, "code/p4/p7_rshiny/data/fit.prototype3.male.rds")
saveRDS(means.prototype3.male, "code/p4/p7_rshiny/data/means.prototype3.male.rds")
saveRDS(fit.standard.male, "code/p4/p7_rshiny/data/fit.standard.male.rds")
saveRDS(bhaz.prototype3.male, "code/p4/p7_rshiny/data/bhaz.prototype3.male.rds")
saveRDS(bhaz.standard.male, "code/p4/p7_rshiny/data/bhaz.standard.male.rds")

### Read in validation datasets
imps.valid <- readRDS(paste("data/p4/dfs_valid_", "female", ".rds", sep = ""))

### Pick one
data.valid <- imps.valid[[1]]

### Obtain a patient
pat <- data.valid[1, ]

### Assign all variables and de-identfiy...
pat <- dplyr::select(pat, -c(cvd_time, cvd_indicator, cvd_time_prim, cvd_indicator_prim,
                        cvd_time_hes, cvd_indicator_hes, cvd_time_death, cvd_indicator_death, 
                        patid, pracid, gender, cvd_ev_prim_aj, statins, antihypertensives))
pat$age <- 50
pat$bmi <- 25
pat$sbp_unexposed <- 130
pat$nonhdl <- 3
pat$smoking[1] <- "Non-smoker"
pat$ethnicity[1] <- "white"
pat$IMD <- 10

### The prototype3 dataset didn't have sbp variability
### I am therefore going to put in a mean value for when estimating the risks of the standard model, meaning
### the two aren't directly comparable
imps.valid.standard <- readRDS(paste("data/dfs_valid_", "female", ".rds", sep = ""))
pat$sbp_var <- mean(imps.valid.standard[[1]]$sbp_var)
pat$cholhdl_ratio <- mean(imps.valid.standard[[1]]$cholhdl_ratio)

saveRDS(pat, "code/p4/p7_rshiny/data/fake.pat.rds")


### Save the HRs/ORs

### Total effects
total_hr_ah <- readRDS("data/p4/offsets_lnHR_ah.rds")
total_hr_statins <- readRDS("data/p4/offsets_lnHR_statins.rds")
total_hr_smoking1 <- readRDS("data/p4/offsets_lnHR_smoking_dummy1_total.rds")
total_hr_smoking2 <- readRDS("data/p4/offsets_lnHR_smoking_dummy2_total.rds")

saveRDS(total_hr_ah, "code/p4/p7_rshiny/data/offsets_lnHR_ah.rds")
saveRDS(total_hr_statins, "code/p4/p7_rshiny/data/offsets_lnHR_statins.rds")
saveRDS(total_hr_smoking1, "code/p4/p7_rshiny/data/offsets_lnHR_smoking_dummy1_total.rds")
saveRDS(total_hr_smoking2, "code/p4/p7_rshiny/data/offsets_lnHR_smoking_dummy2_total.rds")

### Direct effects
direct_RR_smoking_initiation <- readRDS("data/p4/direct_RR_smoking_initiation")
direct_RR_smoking_cessation <- readRDS("data/p4/direct_RR_smoking_cessation")
direct_RR_nonhdl <- readRDS("data/p4/direct_RR_nonhdl")
direct_RR_bmi <- readRDS("data/p4/direct_RR_bmi")
direct_RR_sbp <- readRDS("data/p4/direct_RR_sbp")

saveRDS(direct_RR_smoking_initiation, "code/p4/p7_rshiny/data/direct_RR_smoking_initiation")
saveRDS(direct_RR_smoking_cessation, "code/p4/p7_rshiny/data/direct_RR_smoking_cessation")
saveRDS(direct_RR_nonhdl, "code/p4/p7_rshiny/data/direct_RR_nonhdl")
saveRDS(direct_RR_bmi, "code/p4/p7_rshiny/data/direct_RR_bmi")
saveRDS(direct_RR_sbp, "code/p4/p7_rshiny/data/direct_RR_sbp")
print("FINISHED")
