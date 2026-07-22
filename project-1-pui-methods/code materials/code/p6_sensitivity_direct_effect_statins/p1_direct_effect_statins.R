###
### This file will repeat the entire analysis of the two component model with
### the inclusion of a direct effect of statins.
### It is broken up into different sections determined by the files within the code/p2_devel folder.
###

# preliminaries -----------------------------------------------------------

#####################
### preliminaries ###
#####################

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Load packages
library(survival)

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Extract arguments from command line
args <- commandArgs(trailingOnly = T)
gender <- as.numeric(args[1])
print(paste("gender = ", gender))

### Define model number
model <- 4

### Define burnout
burnout <- 180

### Define pleiotropic effect as a percentage increase the total effect
pleio_effect_multiplier <- 1.2

### Read in HR for antihypertensives
lnHR_ah <- readRDS("data/offsets_total_lnHR_ah.rds")

### For statins, add pleiotropic effect
## Read in total effect identified from literature and converted to HR,
### and remove from log scale to simplify calculations
lnHR_statins_lit <- readRDS("data/offsets_total_lnHR_statins.rds")
HR_statins_lit <- exp(lnHR_statins_lit)

## Make this pleio_effect_multiplier times bigger 
## (i.e., the reduction is pleio_effect_multiplier times bigger, rather than the relative risk value itself)
HR_statins_total <- 1 - pleio_effect_multiplier*(1 - HR_statins_lit)
testthat::expect_equal(pleio_effect_multiplier, (1 - HR_statins_total)/(1 - HR_statins_lit))

## Get the direct effect of statins as the difference between the two
HR_statins_direct <- HR_statins_total/HR_statins_lit

## Return to log scale for application in models
lnHR_statins_direct <- log(HR_statins_direct)
lnHR_statins_total <- log(HR_statins_total)
saveRDS(lnHR_statins_direct, "data/sens_pleiotropic_offsets_direct_lnHR_statins.rds")
saveRDS(lnHR_statins_total, "data/sens_pleiotropic_offsets_total_lnHR_statins.rds")

### Define lnHR_statins, used through code as 'total effect'
lnHR_statins <- lnHR_statins_total

# p1_fit_model ------------------------------------------------------------

print(paste("Start p1_fit_model", Sys.time()))

######################
### p1_fit_model.R ###
######################

### Read in development data
df_devel <- readRDS(paste("data/df_imp_devel_", gender, sep = ""))

### Read in the interval censored outcome times
cohort_split_times <- readRDS(paste("data/cohort_split_times_augmented_burnout", burnout, ".rds", sep = ""))

### Merge imputed data with interval censored data
df_devel <- merge(dplyr::select(df_devel, -c("cvd_time", "cvd_indicator",
                                             "cvd_time_prim", "cvd_time_hes", "cvd_time_death",
                                             "cvd_indicator_prim", "cvd_indicator_hes", "cvd_indicator_death")),
                  cohort_split_times,
                  by.x = "patid",
                  by.y = "patid")

###
### Create offset terms
###

###
### For statin and antihypertensive use, we multiply the indicator by the log-hazard ratio
### These are total effects
### The total effect of statins includes the pleiotropic effect
df_devel$offset_statins_timevar_lnHR <- lnHR_statins_total*df_devel$med_status_statins
df_devel$offset_ah_timevar_lnHR <- lnHR_ah*df_devel$med_status_ah

### Sort by patid and tstart
df_devel <- dplyr::arrange(df_devel, patid, tstart)

### Read in variables that will be interacted with the spline for age
inter_age_rcs <- readRDS(paste("data/var_model_inter_age_rcs", gender, ".rds", sep = ""))

### Define formula (create vector with all terms for the formula, which will be turned into a formula object)
if (model == 4){
  full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))", 
                        paste("age*", inter_age_rcs, sep = ""),
                        "age*rms::rcs(IMD, c(1,10,20))", 
                        "age*rms::rcs(sbp, 4)",
                        "age*rms::rcs(bmi, 4)",
                        "age*rms::rcs(nonhdl, 4)",
                        "offset(offset_statins_timevar_lnHR)",
                        "offset(offset_ah_timevar_lnHR)")
} 

### Create formula
model_formula_full <- as.formula(
  paste("survival::Surv(tstart, cvd_time, cvd_indicator) ~ ", paste(full_formula_vec, sep = "", collapse = "+"), sep = "", collapse = "")
)

### Fit a cox model and get baseline hazard
fit <- survival::coxph(model_formula_full, data = df_devel, model = TRUE)
bhaz <- survival::basehaz(fit, centered = TRUE)
bhaz_uncent <- survival::basehaz(fit, centered = FALSE)

### Save fit
saveRDS(fit, paste("data/sens_pleiotropic_fit_", gender, "_model", model, ".rds", sep = ""))
saveRDS(bhaz, paste("data/sens_pleiotropic_bhaz_", gender, "_model", model, ".rds", sep = ""))
saveRDS(bhaz_uncent, paste("data/sens_pleiotropic_bhaz_uncent_", gender, "_model", model, ".rds", sep = ""))


# p2_calculate_cf_surv_times ------------------------------------------------------------

print(paste("Start p2_calculate_cf_surv_times", Sys.time()))

####################################
### p2_calculate_cf_surv_times.R ###
####################################

### Program to estimate counterfactual survival times if individuals had followed treatment strategy we are
### predicting under.

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
saveRDS(cf_surv_times, paste("data/sens_pleiotropic_cf_surv_times_", gender, "_model", model, ".rds", sep = ""))


# p3_calibration_ph -------------------------------------------------------

print(paste("Start p3_calibration_ph", Sys.time()))

###########################
### p3_calibration_ph.R ###
###########################

###
### Estimate calibration
###

### Read in validation data fresh
df_valid <- readRDS(paste("data/df_imp_valid_", gender, sep = ""))

### Replace cvd_time with counterfactual survival times
df_valid <- dplyr::select(df_valid, -cvd_time)
df_valid <- merge(df_valid, cf_surv_times, by.x = "patid", by.y = "patid") |>
  dplyr::rename(cvd_time = cvd_time_cf)

### Assign time and status variables for compatibility with est_calib_ph (used in boot_func_for_calib_metrics_CI)
df_valid <- dplyr::mutate(df_valid, time = cvd_time, status = cvd_indicator)

### Create offset variable with appropriate names
df_valid_pred <- 
  dplyr::mutate(df_valid,
                offset_statins_timevar_lnHR = 0,
                offset_ah_timevar_lnHR = 0)

##########################
### Assess calibration ###
##########################

##############################
### PH regression approach ###
##############################
print(paste("PH", Sys.time()))
calib_ph <- est_calib_ph(data = df_valid_pred, 
                         fit = fit, 
                         bhaz = bhaz, 
                         time = round(10*365.25))

#####################################
### Get CI's for ICI, E50 and E90 ###
#####################################
print(paste("PH boot", Sys.time()))
calib_ph_boot <- boot::boot(data = df_valid_pred, 
                            statistic = boot_func_for_calib_metrics_CI, 
                            R = 1000, fit = fit, bhaz = bhaz, time = round(10*365.25), nk = 4)

###########################
### KM grouped approach ###
###########################
print(paste("KM group", Sys.time()))
calib_km_group <- est_calib_plot_group(data = df_valid_pred, 
                                       fit = fit, 
                                       bhaz = bhaz,  
                                       time = round(10*365.25),
                                       n.groups = 50,
                                       CI = TRUE)

### Extract plot data and save
df_calib_smooth <- calib_ph[["calib_data"]]
df_calib_grouped <- calib_km_group[["plot"]]$data
saveRDS(df_calib_smooth, paste("data/sens_pleiotropic_calib_ph_df_smooth_", gender, "_model", model, ".rds", sep = ""))
saveRDS(df_calib_grouped, paste("data/sens_pleiotropic_calib_ph_df_grouped_", gender, "_model", model, ".rds", sep = ""))

### Print and save ICI, E50, E90 and the CI
print(paste("ICI = ", calib_ph[["ICI"]]))
print(paste("E50 = ", calib_ph[["E50"]]))
print(paste("E90 = ", calib_ph[["E90"]]))
saveRDS(calib_ph[["ICI"]], paste("data/sens_pleiotropic_calib_ph_ICI_", gender, "_model", model, ".rds", sep = ""))
saveRDS(calib_ph[["E50"]], paste("data/sens_pleiotropic_calib_ph_E50_", gender, "_model", model, ".rds", sep = ""))
saveRDS(calib_ph[["E90"]], paste("data/sens_pleiotropic_calib_ph_E90_", gender, "_model", model, ".rds", sep = ""))
saveRDS(calib_ph_boot$t, paste("data/sens_pleiotropic_calib_ph_CI_ICI_E50_E90_", gender, "_model", model, ".rds", sep = ""))

# p4_discrimination -------------------------------------------------------

print(paste("Start p4_discrimination", Sys.time()))

###########################
### p4_discrimination.R ###
###########################

### Estimate risks
### Note that changing the time, which alter the risks, but not the order of the risks, and therefore
### will not effect the estimation of Harrels C
surv <- as.numeric(est_surv_offset(newdata = df_valid_pred, fit = fit, bhaz = bhaz, time = 10*365.25))

### Estimate C-statistics
Cstat_object <- intsurv::cIndex(time = df_valid_pred$cvd_time, 
                                event = df_valid_pred$cvd_indicator, 
                                risk_score = 1 - surv)

### Calculate confidence interval and add to object
Cstat <- Cstat_object[["index"]]

## Get logit
logit_Cstat <- log(Cstat/(1-Cstat))
## Use delta method to get standard error (sd) of the logit of the C-statistic (which is a proportion)
logit_se <- sqrt(1/(Cstat*(1-Cstat)*Cstat_object[["comparable"]]))
## Get confidence interval by transformation logit_Cstat and upper and lower bounds back onto proper scale
C_lower <- 1/(1 + exp(-(logit_Cstat - qnorm(0.975)*logit_se)))
C_upper <- 1/(1 + exp(-(logit_Cstat + qnorm(0.975)*logit_se)))

### Add to object
Cstat_object[["index_lower"]] <- C_lower
Cstat_object[["index_upper"]] <- C_upper
print(Cstat_object)

### Save output
saveRDS(Cstat_object, paste("data/sens_pleiotropic_discrim_", gender, "_model", model, ".rds", sep = ""))

# clinical exemplar -------------------------------------------------------

print(paste("Start clinical exemplar", Sys.time()))

### I now need to exemplify the differences from the other analyses