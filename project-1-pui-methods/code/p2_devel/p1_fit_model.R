###
### Program to fit models 0 - 4
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Extract model number from command line
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])
print(paste("model = ", model))

### Define burnout
burnout <- 180

### Define HR for offsets
lnHR_statins <- readRDS("data/offsets_total_lnHR_statins.rds")
lnHR_ah <- readRDS("data/offsets_total_lnHR_ah.rds")
lnHR_sbp <- readRDS("data/offsets_direct_lnHR_sbp.rds")
lnHR_bmi <- readRDS("data/offsets_direct_lnHR_bmi.rds")
lnHR_nonhdl <- readRDS("data/offsets_direct_lnHR_nonhdl.rds")

###
### Write a function to fit the model for female/male cohorts
###
fit_model <- function(gender){
  
  ### Read in development data
  df_devel <- readRDS(paste("data/df_imp_devel_", gender, sep = ""))
  
  ### Read in the interval censored outcome times NB CHANGE WHEN FINISHED
  ### Note the time varying medication status variables are defined differently for models 1 and 2 compared to models 3 and 4.
  if (model %in% c(0,1,2)){
    cohort_split_times <- readRDS(paste("data/cohort_split_times_burnout", burnout, ".rds", sep = ""))
  } else if (model %in% c(3,4,5,6,7)){
    cohort_split_times <- readRDS(paste("data/cohort_split_times_augmented_burnout", burnout, ".rds", sep = ""))
  }

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
  ### This is used in all four models (noting that the indicators are different for models 1 and 2, compared to 3 and 4)
  df_devel$offset_statins_timevar_lnHR <- lnHR_statins*df_devel$med_status_statins
  df_devel$offset_ah_timevar_lnHR <- lnHR_ah*df_devel$med_status_ah
  
  ###
  ### For model model 3, apply an offset for differences in sbp, nonhdl and BMI.
  ### We create adjusted versions of SBP, nonhdl and BMI to be relative to some baseline, to apply the offsets
  ### These are also capped:
  ### Assume no benefit for SBP below 120, so it is capped there
  ### Assume no benefit for BMI below 25, so it is capped there
  ###
  
  ### SBP relative to 120 (no benefit for being lower than 120)
  ### BMI relative to 25 (no benefit for being lower than 25)
  ### nonhdl relative to 4
  df_devel <- dplyr::mutate(df_devel,
                            sbp_adj = dplyr::case_when(sbp < 120 ~ 0,
                                                       TRUE ~ sbp - 120),
                            bmi_adj = dplyr::case_when(bmi < 25 ~ 0,
                                                       TRUE ~ bmi - 25),
                            nonhdl_adj = dplyr::case_when(nonhdl < 2.6 ~ 0,
                                                          TRUE ~ nonhdl - 2.6)
  )
  
  ### Create appropriate offsets for baseline variables (i.e. for dropping in effect of variables at baseline)
  df_devel$offset_sbp_lnHR <- lnHR_sbp*df_devel$sbp_adj
  df_devel$offset_bmi_lnHR <- lnHR_bmi*df_devel$bmi_adj
  df_devel$offset_nonhdl_lnHR <- lnHR_nonhdl*df_devel$nonhdl_adj
  
  ### Sort by patid and tstart
  df_devel <- dplyr::arrange(df_devel, patid, tstart)

  ### Read in variables that will be interacted with the spline for age
  inter_age_rcs <- readRDS(paste("data/var_model_inter_age_rcs", gender, ".rds", sep = ""))
  
  ### Define formula (create vector with all terms for the formula, which will be turned into a formula object)
  if (model == 1){
    full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))", 
                          paste("age*", inter_age_rcs, sep = ""), 
                          "age*rms::rcs(IMD, c(1,10,20))",
                          "age*rms::rcs(sbp, 4)",
                          "age*rms::rcs(bmi, 4)",
                          "age*rms::rcs(nonhdl, 4)", 
                          "offset(offset_ah_timevar_lnHR)")
  } else if (model == 2){
    full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))",
                          paste("age*", inter_age_rcs, sep = ""), 
                          "age*rms::rcs(IMD, c(1,10,20))", 
                          "age*rms::rcs(sbp_unexposed, 4)",
                          "age*rms::rcs(bmi, 4)",
                          "age*rms::rcs(nonhdl, 4)", 
                          "offset(offset_ah_timevar_lnHR)")
  } else if (model == 3){
    full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))", 
                          paste("age*", inter_age_rcs, sep = ""),
                          "age*rms::rcs(IMD, c(1,10,20))", 
                          "offset(offset_statins_timevar_lnHR)",
                          "offset(offset_ah_timevar_lnHR)",
                          "offset(offset_bmi_lnHR)", 
                          "offset(offset_nonhdl_lnHR)",
                          "offset(offset_sbp_lnHR)")
  } else if (model == 4){
    full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))", 
                          paste("age*", inter_age_rcs, sep = ""),
                          "age*rms::rcs(IMD, c(1,10,20))", 
                          "age*rms::rcs(sbp, 4)",
                          "age*rms::rcs(bmi, 4)",
                          "age*rms::rcs(nonhdl, 4)",
                          "offset(offset_statins_timevar_lnHR)",
                          "offset(offset_ah_timevar_lnHR)")
  } else if (model == 5){
    full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))", 
                          paste("age*", inter_age_rcs, sep = ""),
                          "age*rms::rcs(IMD, c(1,10,20))", 
                          "age*rms::rcs(bmi, 4)",
                          "age*rms::rcs(nonhdl, 4)", 
                          "offset(offset_statins_timevar_lnHR)",
                          "offset(offset_ah_timevar_lnHR)",
                          "offset(offset_sbp_lnHR)")
  } else if (model == 6){
    full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))", 
                          paste("age*", inter_age_rcs, sep = ""),
                          "age*rms::rcs(IMD, c(1,10,20))", 
                          "age*rms::rcs(sbp, 4)",
                          "age*rms::rcs(nonhdl, 4)", 
                          "offset(offset_statins_timevar_lnHR)",
                          "offset(offset_ah_timevar_lnHR)",
                          "offset(offset_bmi_lnHR)")
  } else if (model == 7){
    full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))", 
                          paste("age*", inter_age_rcs, sep = ""),
                          "age*rms::rcs(IMD, c(1,10,20))", 
                          "age*rms::rcs(sbp, 4)",
                          "age*rms::rcs(bmi, 4)", 
                          "offset(offset_statins_timevar_lnHR)",
                          "offset(offset_ah_timevar_lnHR)",
                          "offset(offset_nonhdl_lnHR)")
  } else if (model == 0){
    full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))", 
                          paste("age*", inter_age_rcs, sep = ""),
                          "age*rms::rcs(IMD, c(1,10,20))",
                          "age*rms::rcs(sbp, 4)",
                          "age*rms::rcs(bmi, 4)",
                          "age*rms::rcs(nonhdl, 4)")
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
  saveRDS(fit, paste("data/fit_", gender, "_model", model, ".rds", sep = ""))
  saveRDS(bhaz, paste("data/bhaz_", gender, "_model", model, ".rds", sep = ""))
  saveRDS(bhaz_uncent, paste("data/bhaz_uncent_", gender, "_model", model, ".rds", sep = ""))
  
}

### Apply function
for (gender in 1:2){
  print(paste("Gender = ", gender, Sys.time()))
  fit_model(gender)
}

### Finished
print(paste("FINISHED", Sys.time()))


