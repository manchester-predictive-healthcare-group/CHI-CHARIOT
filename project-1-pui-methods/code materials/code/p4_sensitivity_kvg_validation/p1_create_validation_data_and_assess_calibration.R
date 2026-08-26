###
### Estimate the IPACW's for K&VG validation approach
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

### Load survival
library(survival)

### Assign burnout
burnout <- 180

### And follow up time
t_fup <- 10*365.25

### Extract arguments from command line
args <- commandArgs(trailingOnly = T)
gender_in <- as.numeric(args[1])
print(paste("gender = ", gender_in))

### Read in validation dataset
df_valid <- readRDS(paste("data/df_imp_valid_", gender_in, sep = ""))

###########################
### Estimate the IPCW's ###
###########################

### Create indicators for estimating probability of being censored (ipcw's)
df_valid <- dplyr::mutate(df_valid, 
                          ipcw_indicator = dplyr::case_when(cvd_time <= t_fup ~ 1 - cvd_indicator,
                                                            cvd_time > t_fup ~ 0),
                          ipcw_time = dplyr::case_when(cvd_time <= t_fup ~ cvd_time,
                                                       cvd_time > t_fup ~ t_fup)
)

###
### Create formula for estimating ipcw's 

### Read in variables that were interacted with the spline for age in prediction model
inter_age_rcs <- readRDS(paste("data/var_model_inter_age_rcs", gender_in, ".rds", sep = ""))

### Create a vector of terms
formula_vec_ipcw <- 
  c("rms::rcs(age, c(25, 40, 57.5, 75))", 
    paste("age*", inter_age_rcs, sep = ""),
    "age*rms::rcs(IMD, c(1,10,20))", 
    "age*rms::rcs(sbp, 4)",
    "age*rms::rcs(bmi, 4)",
    "age*rms::rcs(nonhdl, 4)")

### Define formula for estimating ipcw's
formula_ipcw <- stats::as.formula(
  paste("survival::Surv(ipcw_time, ipcw_indicator) ~ ", paste(formula_vec_ipcw, collapse = " + "))
)

### Fit model to estimate ipcws
fit_est_ipcw <- survival::coxph(formula_ipcw, data = df_valid, model = TRUE)
bhaz_est_ipcw <- survival::basehaz(fit_est_ipcw, centered = TRUE)
bhaz_uncent_est_ipcw <- survival::basehaz(fit_est_ipcw, centered = FALSE)

### Estimate the IPCWs
df_valid$ipcw_weight <- 1/est_surv_eventtime_or_t(newdata = df_valid,
                                                  fit = fit_est_ipcw,
                                                  bhaz = bhaz_est_ipcw,
                                                  eventtime = "cvd_time", # note cvd_time == ipcw_time (prior to merge)
                                                  t = t_fup)


###########################################
### Prepare data for estimating weights ###
###########################################

###
### Read in times at which individuals change medication status, and the time-varying variable for medication status
### NB: We don't use the augmented medication times (used for fitting the prediction model), where
### individuals medication always takes value 0 at zero, and changes when individual deviates from medication
### status at baseline. Instead, we use the original medication status times, as we want to stratify dependent
### on treatment status at baseline when estimating the ipacw weights.
medication_times <- readRDS(paste("data/cohort_split_times_burnout", burnout, ".rds", sep = ""))

### Create indicators for original cvd_time and cvd_indcator in validation dataset
df_valid <- dplyr::mutate(df_valid, time_original = cvd_time, status_original = cvd_indicator)

### Merge with validation data frames and add gender, then recombine
df_ipacw <- 
  merge(dplyr::select(df_valid, -c(cvd_time, cvd_indicator, cvd_ev_prim_aj)),
        medication_times,
        by.x = "patid",
        by.y = "patid") |>
  dplyr::arrange(gender, patid, tstart)
colnames(medication_times)
### Convert to data.table for efficiency
data.table::setDT(df_ipacw)

### Create ipacw indicator and time, and reduce dataset to first row for each individual, 
### which is the time they deviate from baseline treatment strategy
df_ipacw <- df_ipacw[
  # Create ipacw event and time
  , `:=`(
    ipacw_indicator = as.integer(.N > 1),
    ipacw_time = cvd_time[1]
  ),
  by = patid][tstart == 0] |>
  as.data.frame()

### View these
# View(df_ipacw[1:20,])

### Note that in this dataset, 'cvd_time' and 'cvd_indicator' variables are already correctly defined for
### assessing calibration in the artificially censored cohort. cvd_indicator = 1 if the individual had
### an event before being censored by either mechanism, and = 0 otherwise. cvd_time = time until either
### event, or censored by either mechanism. This is because they were taken from the 'medication_times'
### dataset when merging above, not from the df_valid dataset.

### Just need to rename these to 'time' and 'status' to align with functions used later.
df_ipacw <- dplyr::mutate(df_ipacw, time = cvd_time, status = cvd_indicator)

###
### NB:
### ipacw_time = time until min of artificial censoring, censoring, event
### ipacw_indicator = indicator = 1 when artificial censoring occurs first
### ipcw_time = time until min of censoring, event
### ipcw_indicator = indicator = 1 when censoring occurs first
### time = time until min of artificial censoring, censoring, event
### status = 1 when event occurs first
### time_original = time until min of censoring, event
### status_original = 1 when event occurs first
###

# ### Estimate the IPCWs at 'ipacw time'
# ### This doesn't actually effect the results. The only individuals for which ipacw_time and ipcw_time
# ### differ, are those who are aritificially censored, don't have an event, and contribute zero to the outcome
# ### (look specifically at the functions used to estimate calibration to understand this)
# df_ipacw$ipcw_weight <- 1/est_surv_eventtime_or_t(newdata = df_ipacw,
#                                                   fit = fit_est_ipcw,
#                                                   bhaz = bhaz_est_ipcw,
#                                                   eventtime = "ipacw_time", # note cvd_time == ipcw_time (prior to merge)
#                                                   t = t_fup)

#########################
### Estimate IPACW's  ###
#########################

### We are going to build stratified models for estimating the IPACW's, given the probability of remaining
### off treatment (if off treatment at baseline), will be distributed differently than the probability of
### remaining on treatment (if on treatment at baseline). We will consider being on statins, antihypertensives, or
### both (at baseline), separately.

### Given the higher number of individuals off treatment at baseline, we will use a different set of predictors for
### this model. The other three models for on treatment at baseline, will share the same, reduced set of predictors.

###
### Sample size calculations for these models (to help decide how many predictors to include in the models for estimating
### the weights for those on treatment at baseline)
###

### Formula for off treatment at baseline (same predictors at the prediction model, 
### which are different for female and male cohorts)
formula_vec_off_treatment <-
  c("rms::rcs(age, c(25, 40, 57.5, 75))", 
    paste("age*", inter_age_rcs, sep = ""),
    "age*rms::rcs(IMD, c(1,10,20))", 
    "age*rms::rcs(sbp, 4)",
    "age*rms::rcs(bmi, 4)",
    "age*rms::rcs(nonhdl, 4)")

### Formula for on treatment at baseline (reduced number of predictors for sample size reasons)
formula_vec_on_treatment <- c("rms::rcs(age, c(25, 40, 57.5, 75))",
                              "hypertension",
                              "ra",
                              "af",
                              "ckd",
                              "fhcvd",
                              "diabetes",
                              "IMD", 
                              "rms::rcs(sbp, 4)",
                              "rms::rcs(bmi, 4)",
                              "rms::rcs(nonhdl, 4)")

### Get number of coefficients to be estimated for sample size calculation

## Combine formula vector into a full formula
# Off treament
formula_off_treatment <- stats::as.formula(
  paste("survival::Surv(ipacw_time, ipacw_indicator) ~ ", 
        paste(formula_vec_off_treatment, collapse = " + ")))

# On treatment
formula_on_treatment <- stats::as.formula(
  paste("survival::Surv(ipacw_time, ipacw_indicator) ~ ", 
        paste(formula_vec_on_treatment, collapse = " + "))
)

## Count coefficients (subtract 1 for intercept if not in model)
n_coeff_off <- ncol(model.matrix(formula_off_treatment, data = df_ipacw[1:50000, ]))
n_coeff_on <- ncol(model.matrix(formula_on_treatment, data = df_ipacw[1:50000, ]))

## Put into vectors, to be used in lapply, for sample size calculation for the four models
n_coeff <- c(n_coeff_off, rep(n_coeff_on, 3))

### Get mean_fup and rate for sample size calculation

### This needs to be done separately for gender, but also for each stratified model
ss_mean_fup <- dplyr::group_by(df_ipacw, med_status_statins, med_status_ah) |>
  dplyr::summarise(mean_fup = mean(ipacw_time)/365.25)

ss_rate <- dplyr::group_by(df_ipacw, med_status_statins, med_status_ah) |>
  dplyr::summarise(
    events = sum(ipacw_indicator),
    person_time = sum(ipacw_time)/365.25,
    inc_rate = events/person_time
  )


### Estimate minimum require sample size for each group (defined by treatment status at baseline).
### based off number of predictors, mean follow-up time, and event (artificial censoring) rate
min_ss_obj <- lapply(c(1,2,3,4), function(x){
  pmsampsize::pmsampsize(type = "s", 
                         nagrsquared = 0.15, 
                         shrinkage = 0.9, 
                         parameters = n_coeff[x], 
                         timepoint = 10, 
                         rate = as.numeric(ss_rate[x, "inc_rate"]), 
                         meanfup = as.numeric(ss_mean_fup[x, "mean_fup"]))
})

### Extract the number from the object
min_ss <- lapply(min_ss_obj, function(x){x[["sample_size"]]})
min_ss

### Minimum sample size is well below what we have. So proceed.

###
### Development of models used to estimate IPACWs
### Need to develop 4 models (one for each medication group, separately for female/male).
###

### Create cohorts
df_ipacw_stratified_list <- 
  lapply(c(1,2,3,4), function(baseline_group){
    df_out <- df_ipacw
    if (baseline_group == 1){
      df_out <- subset(df_out, med_status_statins == 0 & med_status_ah == 0)
    } else if (baseline_group == 2){
      df_out <- subset(df_out, med_status_statins == 0 & med_status_ah == 1)
    } else if (baseline_group == 3){
      df_out <- subset(df_out, med_status_statins == 1 & med_status_ah == 0)
    } else if (baseline_group == 4){
      df_out <- subset(df_out, med_status_statins == 1 & med_status_ah == 1)
    }
  })

### Check they sum to 1,000,000
testthat::expect_equal(1000000, nrow(do.call("rbind", df_ipacw_stratified_list)))
### And print sizes of each cohort
lapply(1:4, function(baseline_group){
  nrow(df_ipacw_stratified_list[[baseline_group]])
})
### These are bigger than the minimum sample sizes

###
### Fit models for estimating the IPACWs
###
print("FIT MODELS")
for (baseline_group in 1:4){
  
  print(paste("gender = ", gender_in, ", baseline_group = ", baseline_group, ",", Sys.time(), sep = ""))
  
  ### Extract correct formula
  if (baseline_group == 1){
    ipacw_formula <- formula_off_treatment
  } else {
    ipacw_formula <- formula_on_treatment
  }
  
  ### Extract relevant development dataset
  df_ipacw_stratified <- df_ipacw_stratified_list[[baseline_group]]
  
  ### Fit model
  fit_est_ipacw <- survival::coxph(ipacw_formula, data = df_ipacw_stratified, model = TRUE)
  bhaz_est_ipacw <- survival::basehaz(fit_est_ipacw, centered = TRUE)
  bhaz_uncent_est_ipacw <- survival::basehaz(fit_est_ipacw, centered = FALSE)
  
  ### Save fit
  saveRDS(fit_est_ipacw, paste("data/fit_est_ipacw_", gender_in, "_bgroup", baseline_group, ".rds", sep = ""))
  saveRDS(bhaz_est_ipacw, paste("data/bhaz_est_ipacw_", gender_in, "_bgroup", baseline_group, ".rds", sep = ""))
  saveRDS(bhaz_uncent_est_ipacw, paste("data/bhaz_uncent_est_ipacw_", gender_in, "_bgroup", baseline_group, ".rds", sep = ""))
  
  rm(df_ipacw_stratified)
  
}

###
### Calculate the IPACWs
###

### The weight, is the probability of having not been artificially censored by 10-years
### Run through the seperate stratified datasets and calculate the weights by estimating survival 
### probabilities from the models fitted above.
print("ESTIMATE IPACWs")
for (baseline_group in 1:4){
  
  print(paste("gender = ", gender_in, ", baseline_group = ", baseline_group, ",", Sys.time(), sep = ""))
  
  ### Extract data 
  df_temp <- df_ipacw_stratified_list[[baseline_group]]
  
  ### Read in fit and basehaz for estimating the weights
  fit_temp <- readRDS(paste("data/fit_est_ipacw_", gender_in, "_bgroup", baseline_group, ".rds", sep = ""))
  bhaz_temp <- readRDS(paste("data/bhaz_est_ipacw_", gender_in, "_bgroup", baseline_group, ".rds", sep = ""))
  
  ### The output from est_surv is a survival probability, so the probability of having not been artificially censored by time t
  ipacw_temp <- est_surv_eventtime_or_t(newdata = df_temp,
                                        fit = fit_temp,
                                        bhaz = bhaz_temp,
                                        eventtime = "ipacw_time",
                                        t = t_fup)
  
  ### The weights are 1/surv
  ### Assign these to the dataset
  df_ipacw_stratified_list[[baseline_group]]$ipacw_weight <- 1/ipacw_temp
  
}

### Recombine into a single dataset and save
df_valid_ipacw <- do.call("rbind", df_ipacw_stratified_list) |> as.data.frame()

### Get combined weights
df_valid_ipacw$comb_weight <- df_valid_ipacw$ipacw_weight * df_valid_ipacw$ipcw_weight

##########################################################
### Assess calibration in artificially censored cohort ###
##########################################################

###
### First need to generate predictions
###

### Read in prediction model (note, we are focused on model 4 here, the two component model)
fit <- readRDS(paste("data/fit_", gender_in, "_model", 4, ".rds", sep = ""))
bhaz <- readRDS(paste("data/bhaz_", gender_in, "_model", 4, ".rds", sep = ""))

### Create offset variables with appropriate names and values, for predictions under treatment
### strategy of 'do nothing'.
df_valid_ipacw <-
  dplyr::mutate(df_valid_ipacw,
                offset_statins_timevar_lnHR = 0,
                offset_ah_timevar_lnHR = 0)

### Get the survival probabilities and risks
df_valid_ipacw$surv <- as.numeric(est_surv_offset(newdata = df_valid_ipacw, fit = fit, bhaz = bhaz, time = t_fup))
df_valid_ipacw$pred <- 1 - df_valid_ipacw$surv

### Cap the weights
weight_cap <- 10
df_valid_ipacw <- dplyr::mutate(df_valid_ipacw, 
                                comb_weight_cap = dplyr::case_when(comb_weight > weight_cap ~ weight_cap,
                                                                   TRUE ~ comb_weight))

### Create an uncensored at t_fup cohort
df_valid_ipacw_uncensored <- subset(df_valid_ipacw, 
                                    time >= t_fup | status == 1 & time < t_fup)


### Save cohort
saveRDS(df_valid_ipacw, paste("data/df_valid_ipacw", gender_in, ".rds"))

###
### Mean calibration using Horvitz-Thompson estimator
###
ht_estimate <- mean(
  (df_valid_ipacw$status == 1 & df_valid_ipacw$time <= t_fup) * 
    (df_valid_ipacw$comb_weight)
)
ht_estimate
mean(df_valid_ipacw$pred)
mean(df_valid_ipacw$pred)/ht_estimate
# This is better, but still doesn't work!!

### Calculate calibration plot using IPCW approach in uncensored cohort
calib_ipcw_uncensored <- est_calib_ipcw_new(data = df_valid_ipacw_uncensored, 
                                            surv = df_valid_ipacw_uncensored$surv, 
                                            t = t_fup,
                                            nk = 7,
                                            weights_in = df_valid_ipacw_uncensored$comb_weight)

### Calculate calibration plot using IPCW approach using the smoothed Horvitz-Thompson approach
calib_smoothed_ht <- est_calib_smoothed_ht(data = df_valid_ipacw, 
                                           surv = df_valid_ipacw$surv, 
                                           t = t_fup,
                                           nk = 7,
                                           weights_in = df_valid_ipacw$comb_weight)


### Calculate calibration plot using weighted pseudo-value approach
calib_pv_ipcw <- est_calib_pv_ipcw(data = df_valid_ipacw, 
                                   surv = df_valid_ipacw$surv, 
                                   t = t_fup,
                                   nk = 7,
                                   weights_in = df_valid_ipacw$comb_weight)

calib_ipcw_uncensored[["ICI"]]
calib_smoothed_ht[["ICI"]]
calib_pv_ipcw[["ICI"]]
calib_ipcw_uncensored[["E50"]]
calib_smoothed_ht[["E50"]]
calib_pv_ipcw[["E50"]]
calib_ipcw_uncensored[["E90"]]
calib_smoothed_ht[["E90"]]
calib_pv_ipcw[["E90"]]
calib_ipcw_uncensored[["plot"]]
calib_smoothed_ht[["plot"]]
calib_pv_ipcw[["plot"]]

### Write a function to save these plots
save_calibration_plot <- function(calib_obj,
                                  filename,
                                  width = 7,
                                  height = 7,
                                  res = 300) {
  
  ragg::agg_png(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = res
  )
  
  print(calib_obj$plot)
  
  dev.off()
}

### Save the plots
save_calibration_plot(calib_ipcw_uncensored, paste0("figures/calibration_moderate_kvg_ipcw_uncesnored", gender_in, ".png"))
save_calibration_plot(calib_smoothed_ht, paste0("figures/calibration_moderate_kvg_ht", gender_in, ".png"))
save_calibration_plot(calib_pv_ipcw, paste0("figures/calibration_moderate_kvg_pv", gender_in, ".png"))

###
### Function to create a summary table of calibration metrics
###
create_calibration_summary <- function(...,
                                       model_names = NULL,
                                       digits = 3) {
  
  ### Store all calibration objects supplied via ...
  calib_list <- list(...)
  
  ### If model names are not supplied,
  ### create default names
  if (is.null(model_names)) {
    model_names <- paste0(
      "Model_",
      seq_along(calib_list)
    )
  }
  
  ### Check the number of names matches the number of calibration objects
  if (length(model_names) != length(calib_list)) {
    stop(
      "Length of model_names must equal the number of calibration objects."
    )
  }
  
  ### Extract ICI, E50, and E90 from each calibration object and combine into a single data frame
  summary_df <- do.call(
    rbind,
    lapply(
      seq_along(calib_list),
      function(i) {
        
        ### Create one row per calibration object
        data.frame(
          Method = model_names[i],
          
          ### Round metrics to requested number of digits
          ICI = round(
            calib_list[[i]]$ICI,
            digits = digits
          ),
          
          E50 = round(
            calib_list[[i]]$E50,
            digits = digits
          ),
          
          E90 = round(
            calib_list[[i]]$E90,
            digits = digits
          )
        )
      }
    )
  )
  
  ### Remove row names created during rbind
  rownames(summary_df) <- NULL
  
  ### Return summary table
  return(summary_df)
}


###
### Example usage
###

###
### Create summary table
###
results_table <- create_calibration_summary(
  calib_ipcw_uncensored,
  calib_smoothed_ht,
  calib_pv_ipcw,
  model_names = c(
    "BLR",
    "Smoothed HT",
    "Pseudo-value"
  ),
  digits = 3
)

###
### View results
###
print(results_table)
saveRDS(results_table, paste0("data/calibration_moderate_kvg_metric_table", gender_in, ".rds"))
write.csv(results_table, paste0("data/calibration_moderate_kvg_metric_table", gender_in, ".csv"))
print(paste("FINISHED", Sys.time()))
