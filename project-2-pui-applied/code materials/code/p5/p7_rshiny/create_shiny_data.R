### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

###
### Read in female and male models
###

### Read in 
### - pui model (lists of 10)
### - bhaz (list of 10 vectors, each containing baseline hazard at years 1:10)
### - mean offsets (used to estimate risk scores)
### These are the small versions of the model objects, which are lightweight and redacted

# Female
fit_pui_female <- readRDS("data/p5/small_fits_model_female.rds")
bhaz_pui_female <- readRDS("data/p5/small_bhaz_model_female.rds")
means_pui_female <- readRDS("data/p5/small_means_model_female.rds")

# Male
fit_pui_male <- readRDS("data/p5/small_fits_model_male.rds")
bhaz_pui_male <- readRDS("data/p5/small_bhaz_model_male.rds")
means_pui_male <- readRDS("data/p5/small_means_model_male.rds")

### Save these to Rshiny folder
saveRDS(fit_pui_female, "code/p5/p7_rshiny/data/fit_pui_female.rds")
saveRDS(bhaz_pui_female, "code/p5/p7_rshiny/data/bhaz_pui_female.rds")
saveRDS(means_pui_female, "code/p5/p7_rshiny/data/means_pui_female.rds")

saveRDS(fit_pui_male, "code/p5/p7_rshiny/data/fit_pui_male.rds")
saveRDS(bhaz_pui_male, "code/p5/p7_rshiny/data/bhaz_pui_male.rds")
saveRDS(means_pui_male, "code/p5/p7_rshiny/data/means_pui_male.rds")

###
### Create a default patient
###

### Read in dataset to use for formatting (select no rows/data)
pat <- readRDS("data/p5/cohort_prototype3.rds")[0,]

### Create fake patient
pat <- tibble::add_row(pat, 
                       age = 60,
                       sbp = 130, 
                       bmi = 25, 
                       nonhdl = 3, 
                       IMD = 10,
                       diabetes = factor("Absent", levels = levels(pat$diabetes)),
                       smoking = factor("Non-smoker", levels = levels(pat$smoking)),
                       ethnicity = factor("white", levels = levels(pat$ethnicity)),
                       hypertension      = factor("Absent", levels = c("Absent", "Present")),
                       ra                = factor("Absent", levels = c("Absent", "Present")),
                       af                = factor("Absent", levels = c("Absent", "Present")),
                       ckd               = factor("Absent", levels = c("Absent", "Present")),
                       smi               = factor("Absent", levels = c("Absent", "Present")),
                       fhcvd             = factor("Absent", levels = c("Absent", "Present")),
                       migraine          = factor("Absent", levels = c("Absent", "Present")),
                       sle               = factor("Absent", levels = c("Absent", "Present")),
                       cortico           = factor("Absent", levels = c("Absent", "Present")),
                       antipsy           = factor("Absent", levels = c("Absent", "Present")),
                       copd              = factor("Absent", levels = c("Absent", "Present")),
                       int_dis           = factor("Absent", levels = c("Absent", "Present")),
                       downs             = factor("Absent", levels = c("Absent", "Present")),
                       oral_cancer       = factor("Absent", levels = c("Absent", "Present")),
                       brain_cancer      = factor("Absent", levels = c("Absent", "Present")),
                       lung_cancer       = factor("Absent", levels = c("Absent", "Present")),
                       blood_cancer      = factor("Absent", levels = c("Absent", "Present")),
                       pre_eclampsia     = factor("Absent", levels = c("Absent", "Present")),
                       postnatal_depression = factor("Absent", levels = c("Absent", "Present"))) |>
  dplyr::mutate(caltime = 5537) |>
  as.data.frame()

## Remove unnecessary variable (these will be missing anyway)
pat <- dplyr::select(pat, -c(cvd_time, cvd_indicator, cvd_time_prim, cvd_indicator_prim,
                             cvd_time_hes, cvd_indicator_hes, cvd_time_death, cvd_indicator_death,
                             patid, pracid, gender, statins, antihypertensives))

### Set the offsets for predicting continue on current treatment
pat <- dplyr::mutate(pat, 
                     offset_ah_timevar_lnHR = 0,
                     offset_statins_timevar_lnHR = 0,
                     offset_smoking_timevar_dummy1_lnHR = 0,
                     offset_smoking_timevar_dummy2_lnHR = 0)

saveRDS(pat, "code/p5/p7_rshiny/data/fake_pat.rds")


###
### Read in HRs
###

### Total effects
HR_ah_total <- exp(readRDS("data/p5/offsets_total_lnHR_ah.rds"))
HR_statins_total <- exp(readRDS("data/p5/offsets_total_lnHR_statins.rds"))
HR_smok_init_total <- exp(readRDS("data/p5/offsets_total_lnHR_smok_init.rds"))
HR_smok_cess_total <- exp(readRDS("data/p5/offsets_total_lnHR_smok_cess.rds"))

### Direct effects
HR_sbp_direct <- exp(readRDS("data/p5/offsets_direct_lnHR_sbp.rds"))
HR_nonhdl_direct <- exp(readRDS("data/p5/offsets_direct_lnHR_nonhdl.rds"))
HR_bmi_direct <- exp(readRDS("data/p5/offsets_direct_lnHR_bmi.rds"))
HR_smok_init_direct <- exp(readRDS("data/p5/offsets_direct_lnHR_smok_init.rds"))
HR_smok_cess_direct <- exp(readRDS("data/p5/offsets_direct_lnHR_smok_cess.rds"))

###
### Save HRs to shiny data folder (not log HRs)
### 

### Total effects
saveRDS(HR_statins_total, "code/p5/p7_rshiny/data/offsets_total_HR_statins.rds")
saveRDS(HR_ah_total, "code/p5/p7_rshiny/data/offsets_total_HR_ah.rds")
saveRDS(HR_smok_init_total, "code/p5/p7_rshiny/data/offsets_total_HR_smok_init.rds")
saveRDS(HR_smok_cess_total, "code/p5/p7_rshiny/data/offsets_total_HR_smok_cess.rds")

### Direct effects
saveRDS(HR_sbp_direct, "code/p5/p7_rshiny/data/offsets_direct_HR_sbp.rds")
saveRDS(HR_nonhdl_direct, "code/p5/p7_rshiny/data/offsets_direct_HR_nonhdl.rds")
saveRDS(HR_bmi_direct, "code/p5/p7_rshiny/data/offsets_direct_HR_bmi.rds")
saveRDS(HR_smok_init_direct, "code/p5/p7_rshiny/data/offsets_direct_HR_smok_init.rds")
saveRDS(HR_smok_cess_direct, "code/p5/p7_rshiny/data/offsets_direct_HR_smok_cess.rds")
print("FINISHED")
