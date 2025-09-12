###
### This program will verify the functions used to estimate risks in the rshiny app and worked example
### This is neccesary, becuase in the rshiny app, we cannot export the full model fit.
### In the normal est_surv function (est_surv_offset_prototype3), the offsets (which are deducted from the linear predictor) 
### are calculated from the fitted object which contains the  development data frame. 
### This has to be removed when sharing the model fit with other users and for the rshiny.
### We also want the worked example to be reproducible.
### In this alternate est_surv function (est_surv_offset_prototype3_mi_manual_offset), 
### the offsets are inputted manually, rather than estimating them
### from the fitted object. I want to verify this leads to the same results.
###
### I will do two tests.
### 1) I will compare risk using normal est_surv functions with risk from the manual offset est_surv functions.
### 2) I will check the function used in rshiny, which averages results over 10 multiply imputed models, is similar to the 
### function operating on one model

###
### Test 1
###

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Read in the rshiny models
fit_prototype3 <- readRDS("code/p4/p7_rshiny/data/fit.prototype3.female.rds")
means_prototype3 <- readRDS("code/p4/p7_rshiny/data/means.prototype3.female.rds")
bhaz_prototype3 <- readRDS("code/p4/p7_rshiny/data/bhaz.prototype3.female.rds")

### Read in pat file
pat <- readRDS("code/p4/p7_rshiny/data/fake.pat.rds") |> dplyr::mutate(offset_ah_timevar_lnHR = 0,
                                                     offset_statins_timevar_lnHR = 0,
                                                     offset_smoking_timevar_dummy1_lnHR = 0,
                                                     offset_smoking_timevar_dummy2_lnHR = 0)

### Give them some elevated risk factors
pat$age <- 70
pat$bmi <- 30
pat$sbp <- 160
pat$nonhdl <- 10
pat$smoking[1] <- "Current"
pat$ethnicity[1] <- "black caribbean"
pat$IMD <- 10

### Read in model fit and bhaz that has been used for model evaluation
fit <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", 2, "_imp", 1, ".rds", sep = ""))
bhaz <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", 2, "_imp", 1, ".rds", sep = ""))

### Estimate predicted risk using function that has been used to do all the model evaluation
surv1 <- est_surv_offset_prototype3(newdata = pat, fit = fit, bhaz = bhaz, time = 10*365.25)

### Estimate risk using the function that will be used in rshiny and worked example, but just for the first imputed dataset
surv2 <- est_surv_offset_prototype3_mi_manual_offset(newdata = pat,
                                                     fit_list = fit_prototype3[1],
                                                     bhaz_list = bhaz_prototype3[1],
                                                     means_statins = means_prototype3[["statins"]][1],
                                                     means_ah = means_prototype3[["ah"]][1],
                                                     means_smoking1 = means_prototype3[["smoking1"]][1],
                                                     means_smoking2 = means_prototype3[["smoking2"]][1],
                                                     time = 10*265.25)

testthat::expect_equal(surv1, surv2)


###
### Now compare mi vs one dataset using only the manual offset versions
###
est_surv_offset_prototype3_mi_manual_offset(newdata = pat,
                                                fit_list = fit_prototype3,
                                                bhaz_list = bhaz_prototype3,
                                                means_statins = means_prototype3[["statins"]],
                                                means_ah = means_prototype3[["ah"]],
                                                means_smoking1 = means_prototype3[["smoking1"]],
                                                means_smoking2 = means_prototype3[["smoking2"]],
                                                time = 10*265.25)

lapply(1:10, function(x) {est_surv_offset_prototype3_mi_manual_offset(newdata = pat,
                                                                          fit_list = fit_prototype3[x],
                                                                          bhaz_list = bhaz_prototype3[x],
                                                                          means_statins = means_prototype3[["statins"]][x],
                                                                          means_ah = means_prototype3[["ah"]][x],
                                                                          means_smoking1 = means_prototype3[["smoking1"]][x],
                                                                          means_smoking2 = means_prototype3[["smoking2"]][x],
                                                                          time = 10*265.25)})

### There is the expected amount of variation in the individual models, and average equals approx the result from the multiple imputation function 